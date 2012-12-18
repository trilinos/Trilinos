import os
from distutils.core import setup, Extension

exopy_PreinstallSeacas = \
  Extension('exopy',
            ['exopy.c'],                 
            libraries = ['exoIIv2c', 'netcdf', 'hdf5_hl', 'hdf5', 'z', 'm'],
            include_dirs = ['/projects/seacas/linux_rhel5/current/include'],
            library_dirs = ['/projects/seacas/linux_rhel5/current/lib'],
           )


relativeLocLen = -len('/TPLs_src/Trilinos/packages/seacas/applications/exopy')
codePath = os.getcwd()[:relativeLocLen]

exopy_Mac = \
  Extension('exopy',
            ['exopy.c'],
            libraries = ['exodus', 'netcdf', 'hl', 'hdf5', 'mpich', 'tmpe', 'pmpich', 'mpl'],
            include_dirs = [
              ( codePath + '/seacas/libraries/exodus/cbind/include' ),
              ( codePath + '/TPLs_src/netcdf/src/include' ),
              ( codePath + '/TPLs_src/hdf5/hl/src' ),
              ( codePath + '/TPLs_src/hdf5/hdf5-1.8.7/src' ),
              '/opt/local/include/mpich2'
            ],
            library_dirs = [
              ( codePath + '/objs/apps/seacas/votd/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared' ),
              ( codePath + '/objs/tpls/netcdf/4.1.3/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared' ),
              ( codePath + '/objs/tpls/hdf5/1.8.7/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared' ),
              ( codePath + '/objs/apps/MPIH/votd/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared' ),
              '/opt/local/lib'
            ]
           )

sierra_machine = os.getenv('SIERRA_MACHINE')
exopy_module = exopy_PreinstallSeacas
if sierra_machine == 'darwin':
  exopy_module = exopy_Mac

setup (
       name = 'ExodusPython',
       version = '1.0',
       description = 'This a python interface to the ExodusII library.',
       author='Mike Veilleux, Timothy Shelton, and Kendall Pierson',
       author_email='sierra-help@sandia.gov',
       ext_modules = [ exopy_module ]
      )

