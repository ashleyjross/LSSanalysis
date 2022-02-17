import os
import argparse
import logging
import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt

# To install pycorr: 
#pip uninstall Corrfunc pycorr
#module load gsl # required by Corrfunc
#python -m pip install git+https://github.com/adematti/pycorr@DESI-1.0.0#egg=pycorr[corrfunc]

from pycorr import TwoPointCorrelationFunction, TwoPointEstimator, utils, project_to_multipoles, project_to_wp, setup_logging
from Cosmo import distance

ds = distance(.31,.69) #eboss fiducial cosmo

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected",default='QSO')
parser.add_argument("--basedir", help="where to find catalogs",default='/global/cscratch1/sd/ajross/ebosscat/')
parser.add_argument("--outdir", help="where to output results",default='/global/cscratch1/sd/ajross/ebossxi/')
parser.add_argument("--bintype",help="log or lin",default='lin')
parser.add_argument("--nthreads",help="number of threads for parallel comp",default=32,type=int)
parser.add_argument("--vis",help="set to y to plot each xi ",default='n')

#only relevant for reconstruction
parser.add_argument("--rectype",help="IFT or MG supported so far",default='IFT')
parser.add_argument("--convention",help="recsym or reciso supported so far",default='reciso')

setup_logging()
args = parser.parse_args()

ttype = args.type
lssdir = args.basedir

if args.bintype == 'log':
    bine = np.logspace(-1.5, 2.2, 80)
if args.bintype == 'lin':
    bine = np.linspace(1e-4, 200, 201)

dirxi = args.outdir


dirname = lssdir 

def compute_correlation_function(mode, tracer='QSO', region='NGC', zlim=(0., np.inf), nthreads=8, dtype='f8', wang=None,fnroot=''):
    data_fn = os.path.join(dirname, 'eBOSS_{}_clustering_data-{}-vDR16.fits'.format(tracer, region))
    data = Table.read(data_fn)

    shifted = None

    randoms_fn = os.path.join(dirname, 'eBOSS_{}_clustering_random-{}-vDR16.fits'.format(tracer, region)) 
    randoms = Table.read(fn) 
  
    corrmode = mode
    if mode == 'wp':
        corrmode = 'rppi'
    if mode == 'multi':
        corrmode = 'smu'
    if corrmode == 'smu':
        edges = (bine, np.linspace(-1., 1., 201)) #s is input edges and mu evenly spaced between 0 and 1
    if corrmode == 'rppi':
        edges = (bine, bine) #transverse and radial separations are  coded to be the same here
        if mode == 'wp':
            edges = (bins,np.linspace(0,40.,41)) #if you want wp, only go out to pi = 40; consider setting pi_max as argument
    def get_positions_weights(catalog):
        mask = (catalog['Z'] >= zlim[0]) & (catalog['Z'] < zlim[1])
        print('Using {:d} rows for {}'.format(mask.sum(),name))
        positions = [catalog['RA'][mask], catalog['DEC'][mask], ds.dc(catalog['Z'][mask])]
        #if weight_type is None:
        #    weights = None
        #else:
        weights = catalog['WEIGHT_FKP'][mask]*catalog['WEIGHT_SYSTOT'][mask]*catalog['WEIGHT_CP'][mask]*catalog['WEIGHT_NOZ'][mask] 
                
        return positions, weights
    
    data_positions, data_weights = get_positions_weights(data)
    randoms_positions, randoms_weights = get_positions_weights(randoms)
    shifted_positions, shifted_weights = None, None
    if shifted is not None:
        shifted_positions, shifted_weights = get_positions_weights(shifted, name='shifted')

    kwargs = {}

    result = TwoPointCorrelationFunction(corrmode, edges, data_positions1=data_positions, data_weights1=data_weights,
                                         randoms_positions1=randoms_positions, randoms_weights1=randoms_weights,
                                         shifted_positions1=shifted_positions, shifted_weights1=shifted_weights,
                                         engine='corrfunc', position_type='rdd', nthreads=nthreads, dtype=dtype, **kwargs)
    #save paircounts
    fn = dirxi+'paircounts_'+fnroot+'.npy'
    result.save(fn)
    return result

ranwt1=False

regl = ['NGC','SGC']

tcorr = ttype
tw = ttype

bsl = [1,4,5,10]

if ttype == 'QSO':
    zmin = 0.8
    zmax = 2.2

for reg in regl:
	print(reg)
	#(sep, xiell), wang = compute_correlation_function('multi', bs, tracer=tcorr, region=reg, nrandoms=args.nran, zlim=(zmin,zmax), weight_type=weight_type,nthreads=args.nthreads)
	fnroot = tw+'eboss'+reg+'_'+str(zmin)+str(zmax)+'DR16'+'_'+args.bintype
	pfn = dirxi+'paircounts_'+fnroot+'.npy'
	result = compute_correlation_function('multi', tracer=tcorr, region=reg, zlim=(zmin,zmax),nthreads=args.nthreads,fnroot=fnroot)
	for bs in bsl:
		result = TwoPointEstimator.load(pfn)
		result.rebin((bs, 1))
		sep,xiell = project_to_multipoles(result)#, wang
		fo = open(dirxi+'xi024'+fnroot+str(bs)+'.dat','w')
		for i in range(0,len(sep)):
			fo.write(str(sep[i])+' '+str(xiell[0][i])+' '+str(xiell[1][i])+' '+str(xiell[2][i])+'\n')
		fo.close()
		if args.vis == 'y':
			if args.bintype == 'log':
				plt.loglog(sep,xiell[0])
			if args.bintype == 'lin':
				plt.plot(sep,sep**2.*xiell[0])
			plt.title(ttype+' '+str(zmin)+'<z<'+str(zmax)+' in '+reg)
			plt.show()    
