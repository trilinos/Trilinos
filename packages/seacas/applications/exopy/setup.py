from distutils.core import setup, Extension

exopy_PreinstallSeacas = \
  Extension('exopy',
            ['exopy.c'],                 
            libraries = ['exoIIv2c', 'netcdf', 'hdf5_hl', 'hdf5', 'z', 'm'],
            include_dirs = ['/projects/seacas/linux_rhel5/current/include'],
            library_dirs = ['/projects/seacas/linux_rhel5/current/lib'],
           )

exopy_CodeSeacas = \
  Extension('exopy',
            ['exopy.c'],
            libraries = ['exodus', 'netcdf', 'hl', 'hdf5', 'mpich', 'mpl', 'tmpe', 'pmpich'], # 'mpih'],
            include_dirs = [
              '/Users/trshelt/sierra/codeMac/seacas/libraries/exodus/cbind/include',
              '/Users/trshelt/sierra/codeMac/TPLs_src/netcdf/src/include',
              '/Users/trshelt/sierra/codeMac/TPLs_src/hdf5/hl/src',
              '/Users/trshelt/sierra/codeMac/TPLs_src/hdf5/hdf5-1.8.7/src',
              '/opt/local/include/mpich2'
            ],
            library_dirs = [
              '/Users/trshelt/sierra/codeMac/objs/apps/seacas/votd/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared',
              '/Users/trshelt/sierra/codeMac/objs/tpls/netcdf/4.1.3/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared',
              '/Users/trshelt/sierra/codeMac/objs/tpls/hdf5/1.8.7/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared',
              '/Users/trshelt/sierra/codeMac/objs/apps/MPIH/votd/darwin-4.4.macports/release/address-model-64/mpi-mpich/runtime-link-shared',
              '/opt/local/lib'
            ]
           )

setup (
       name = 'ExodusPython',
       version = '1.0',
       description = 'This a python interface to the ExodusII library.',
       author='Mike Veilleux, Timothy Shelton, and Kendall Pierson',
       author_email='sierra-help@sandia.gov',
       ext_modules = [ exopy_PreinstallSeacas ]
       #ext_modules = [ exopy_CodeSeacas ]
      )

