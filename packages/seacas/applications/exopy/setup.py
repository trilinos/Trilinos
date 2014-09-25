C Copyright (c) 2014, Sandia Corporation.
C Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 
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

