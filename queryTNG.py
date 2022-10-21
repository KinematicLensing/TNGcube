import requests
import os
import numpy as np
import h5py
import pickle
from TNGcube import Subhalo

baseURL = 'http://www.tng-project.org/api/'
headers = {"api-key": "b703779c1f099efed6f47b91607b1bb1"}


def available_simNames(baseURL=baseURL, headers=headers):
    rInfo = get(baseURL)
    simNames = [sim['name'] for sim in rInfo['simulations']]
    return simNames


def get_snapNum(redshift, simName='TNG50-1', baseURL=baseURL):
    '''Find the Snapshat number that are closest to the input redshift, given simName.'''
    simSnapsInfo = get(baseURL+simName+'/snapshots/')
    snaps_z_list = [simSnapsInfo[i]['redshift'] for i in range(len(simSnapsInfo))]
    snapNum = min(range(len(snaps_z_list)), key=lambda i: abs(snaps_z_list[i]-redshift))
    return snapNum


class QueryTNG():
    '''The manager object to handle the entire query process.'''
    def __init__(self, simName, snapNum, baseURL=baseURL, headers=headers):
        '''
            Args:
                simName : [str] 
                    The simulation name tag associated with the TNG database.
                    available simName : 'TNG50-1', 'TNG100-1' 
        '''

        self.simURL = baseURL+simName+'/'
        self.snapURL = self.simURL+f'snapshots/{int(snapNum)}/'
        self.subhalosURL = self.snapURL + 'subhalos/'

        self.simInfo = get(self.simURL)
        self.snapInfo = get(self.snapURL)

        self.h = self.simInfo['hubble'] # hubble constant of TNG simulations
        self.redshift = self.snapInfo['redshift']
        self.a = 1./(1.+self.redshift)

        self.simName = simName
        self.snapNum = snapNum
    
    def query_subhaloCat(self, mass_range=[1e8, 1e15], vmax_range=[50., 1050.], limit=10, pageID=0):
        '''Query a catalog of available subhalo IDs whose masses and rotation curves are above the input thresholds. 
            Args:
                mass_range : [unit: Msun/h]
                vmax_range : [unit: km/s]
                    query subhalos within the given vmax range 
                    vmax is the maximum value of the spherically-averaged rotation curve
                limit : int
                    number of subhalos to query for each pageID
                pageID : int
                    e.g. if set limie=10:
                        pageID = 0, return the first 10 subhalo IDs indexed from 0~9
                        pageID = 1, return the next 10 subhalos IDs indexed from 10~19
            Returns:
                suhaloIDs : list
                    a list of suhaloIDs satisfying the query condition.
        '''
        mass_range = np.array(mass_range)/1e10*self.h  # change unit to [unit: 1e10 Msun/h]
  
        search_query = f'?limit={limit}&offset={int(pageID*limit)}&mass__gt={mass_range[0]}&mass__lt={mass_range[1]}&vmax__gt={vmax_range[0]}&vmax__lt={vmax_range[1]}'

        searchurl = self.subhalosURL + search_query
        print('search url:', searchurl)
        self.subhaloCat = get(searchurl)
        subhaloIDs = [self.subhaloCat['results'][i]['id'] for i in range(limit)]
        return subhaloIDs

    def query_subhaloInfo(self, subhaloID):
        '''query the overall properties of the subhalo, given subhaloID'''

        subhaloURL = self.subhalosURL + f'{int(subhaloID)}/'

        subhaloInfo_all = get(subhaloURL)

        keys_takeout = ['snap', 'id', 'mass', 'stellarphotometrics_r', 'vmax', 'vmaxrad', 'mass_log_msun']

        subhaloInfo = {key: subhaloInfo_all[key] for key in keys_takeout}

        key_vec = ['cm', 'pos', 'spin', 'vel']
        _xyz = ['_x', '_y', '_z']

        for key in key_vec:
            subhaloInfo[key] = np.array([subhaloInfo_all[key+i] for i in _xyz])

        return subhaloInfo
    
    def query_subhaloCutout(self, subhaloID):
        '''Download the snapshot particle information in the given subhalo field.
            Output file: 
                e.g. if subhaloID=46, output a file at './cutout_46.hdf5'
        '''

        cutout_request = {'gas': 'Coordinates,Masses,Velocities,StarFormationRate,InternalEnergy',
                          'stars': 'Coordinates,Masses,Velocities,GFM_StellarPhotometrics' }
                        # 'dm'   :'Coordinates,Velocities'

        subhaloURL = self.subhalosURL + f'{int(subhaloID)}/'
        print(subhaloURL+'cutout.hdf5')
        fname_cutout = get(subhaloURL+'cutout.hdf5', cutout_request)
        return fname_cutout
    
    def _preprocess_snap_arrs(self, f_hdf5, ptl_type, subhaloInfo):
        ''' 
            Args:
                ptl_type : 'PartType0' (for gas) or 'PartType4' (for stars)
            Returns:
                snap['pos']  : ndarray Nptlx3 [unit: ckpc/h]
                    particle position w.r.t. the subhalo central of mass coordinate
                snap['mass'] : 1d array       [unit: 10^10 Msun/h]
                    particle mass 
                snap['vel']  : ndarray Nptlx3 [unit: km/s]
                    particle velocities w.r.t. the overall systematic velocity of the subhalo
                    (Note: the sqrt(a) factor in the original unit is taken away after processing.)
                snap['SFR']  : 1d array     [unit: Msun/yr]
                    Instantaneous star formation rate of this gas cell.
                snap['InternalEnergy'] : 1d array  [unit: (km/s)^2] (only available for gas ptl.)
                    Internal (thermal) energy per unit mass for this gas cell.
                snap['GFM_StellarPhotometrics'] : ndarray Nptlx8 [unit: mag]
                    Stellar magnitudes in eight bands: U, B, V, K, g, r, i, z.
            
            Reference: 
                https://www.tng-project.org/data/docs/specifications/
        '''

        try :
            snap = {    'pos': f_hdf5[ptl_type]['Coordinates'][:, :], 
                        'vel': f_hdf5[ptl_type]['Velocities'][:, :]*np.sqrt(self.a), 
                            # change the unit of velocity from [sqrt(a) km/s] to [km/s]
                        'mass': f_hdf5[ptl_type]['Masses'][:]   }
        except KeyError:
            print(f'SubhaloID %s does not have any {ptl_type} particles.' % subhaloInfo['id'])
            return None

        # modify the ptl. position and velocity vectors such that they are w.r.t. the subhalo c.m. and system velocity
        for j in range(3):
            snap['pos'][:, j] -= subhaloInfo['cm'][j]
            snap['vel'][:, j] -= subhaloInfo['vel'][j]
        
        if ptl_type == 'PartType0':
            snap['InternalEnergy'] = f_hdf5[ptl_type]['InternalEnergy'][:]
            snap['SFR'] = f_hdf5[ptl_type]['StarFormationRate'][:]
        
        if ptl_type == 'PartType4':
            snap['GFM_StellarPhotometrics'] = f_hdf5[ptl_type]['GFM_StellarPhotometrics'][:,:]
        
        # coordinate re-adjustment on subhaloInfo (w.r.t. the subhalo c.m. and system velocity)
        # note that this is not implemented yet.
        
        return snap
    
    def download_subhalos(self, subhaloIDs):
        '''Download and pre-process the queried subhalo information, and save it as a dictionary that can work with TNGcube.
            Output file:
                e.g. subhaloIDs = [46, 47]
                save the subhalo info. as './cutout_46.pkl' and './cutout_47.pkl'
                to load subhalo info. back: 
                    subhaloInfo = pickle.load(open('cutout_46.pkl', 'rb'))
        '''

        catalog = {}
        for ID in subhaloIDs:
            subhaloInfo = self.query_subhaloInfo(ID)
            fhdf5_cutout = self.query_subhaloCutout(ID)
            
            with h5py.File(fhdf5_cutout, 'r') as f_hdf5:
                snapInfo = {    'gas' :  self._preprocess_snap_arrs(f_hdf5, 'PartType0', subhaloInfo),
                                'stars' : self._preprocess_snap_arrs(f_hdf5, 'PartType4', subhaloInfo)  }

            subhalo = Subhalo(info=subhaloInfo, snap=snapInfo)

            os.remove(fhdf5_cutout)
            fpkl_cutout = f'subhalo_{ID}.pkl'
            pickle.dump(subhalo, open(fpkl_cutout, 'wb'))
            catalog[ID] = subhalo

        return catalog
    

def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json()  # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename  # return the filename string

    return r


if __name__ == '__main__':
    Q_TNG = QueryTNG(simName='TNG50-1', snapNum=75)
