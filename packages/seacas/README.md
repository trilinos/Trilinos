# SEACAS  [[Documentation](http://gsjaardema.github.io/seacas/)]

* [Build instructions](#build-instructions)
* [Configure, Build, and Install SEACAS](#configure-build-and-install-seacas)
* [Testing](#testing)
* [Contact information](#contact-information)

## Build instructions

### Download and build dependencies

There are a few externally developed third-party libraries (TPL) that are required to
build SEACAS

 * [Zoltan](#zoltan) -- required, part of Trilinos
 * [HDF5](#hdf5) -- optional, but highly recommended
 * [Parallel-NetCDF](#parallel-netcdf) -- optional for parallel
 * [NetCDF](#netcdf) -- required with modifications
 * [MatIO](#matio) -- optional

#### Zoltan
Zoltan is a package in Trilinos and it must be enabled for a SEACAS build.

#### HDF5 

If you are using the netcdf-4 capability in the netcdf library or are
using the MatIO library for conversion of exodus to/from matlab
format, then you will need the hdf5 library.

The hdf5 library is used for the netcdf4 capability in netcdf which in
turn is used by exodus.  The netcdf4 capability is typically used for
large models (>150 million elements); if you are not planning to
create or read models of this size and do not want compression
support, you do not have to build hdf5.

   * Download HDF5 from http://www.hdfgroup.org/HDF5/release/obtain5.html
   * untar it, creating a directory will will refer to as `hdf5-X.X.X`
   * `cd` to that directory and enter the command:
     * Serial:
    ```
    ./configure --prefix=${WHERE_TO_INSTALL} --enable-shared --enable-production --enable-debug=no --enable-static-exec
    ```
     * Parallel: 
    ```
    CC=mpicc ./configure --prefix=${WHERE_TO_INSTALL} --enable-shared --enable-production --enable-debug=no --enable-static-exec --enable-parallel
    ```
   * `make && make install`

#### Parallel-NetCDF
  For a parallel build of Trilinos, especially to use the
  auto-decomposition support of the Ioss library, you will need the
  parallel-netcdf library, also known as pnetcdf.

  * Download http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.6.1.tar.gz
  * `tar zxvf parallel-netcdf-1.6.1.tar.gz`
  * `cd` to the `parallel-netcdf-1.6.1` directory and enter the command:
  ```
  CC=mpicc ./configure --disable-fortran --prefix ${WHERE_TO_INSTALL}
  ```
  * `make && make install`
  
#### NetCDF
The most recent released version is recommended. For use with Exodus, some local modifications to the netcdf.h include file are required.  See [NetCDF-Mapping.md](NetCDF-Mapping.md) for an explanation of why these modifications are required (or highly recommended)

 * Download the latest netcdf-c release from http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
 * `tar zxvf netcdf-4.4.0.tar.gz`  (or whatever the latest version is)
 * Modify the following defines in netcdf-4.4.0/include/netcdf.h.

    ```
    #define NC_MAX_DIMS     65536    /* max dimensions per file */
    #define NC_MAX_VARS     524288   /* max variables per file */
    ```

 * `cd netcdf-4.4.0` and enter the command:
    * serial
    ```
    CFLAGS="-I${WHERE_TO_INSTALL}/include" \
    CPPFLAGS="-DNDEBUG" LDFLAGS="-L${WHERE_TO_INSTALL}/lib" \
    ./configure --enable-netcdf-4  \
      --disable-fsync --prefix ${WHERE_TO_INSTALL} \
      --disable-dap --disable-v2
    ```

    * parallel
    ```
    CC='mpicc' CFLAGS="-I${WHERE_TO_INSTALL}/include" \
    CPPFLAGS="-DNDEBUG" LDFLAGS="-L${WHERE_TO_INSTALL}/lib" \
    ./configure --enable-netcdf-4  --enable-pnetcdf \
      --disable-fsync --prefix ${WHERE_TO_INSTALL} \
      --disable-dap --disable-v2
    ```

 * Check the results of the configure and make sure that the listings
   under features are similar to:

   ```
   # Features
   --------
   NetCDF-2 API:	no
   NetCDF-4 API:	yes
   CDF-5 Support:	yes
   HDF5 Support:	yes
   PNetCDF Support:	yes
   NC-4 Parallel Support:	yes
   ```
   There will be other features, but their settings are not important
   for seacas. For a serial build, `PNetCDF` and `NC-4 Parallel Support`
   should be `no`

 * `make && make install`

#### MatIO
The MatIO library is used in the `exo2mat` and `mat2exo` programs which convert an exodus file to and from a MATLAB binary file.  To use this do:

 * Download matio via git:
 * `git clone https://github.com/tbeu/matio.git`
 * `cd matio` and enter the command:
    ```
    ./autogen.sh
    # The -L is to find the hdf5 library...
    export LDFLAGS="-L${WHERE_TO_INSTALL}/lib"
    ./configure --with-hdf5=${WHERE_TO_INSTALL} --enable-mat73 --enable-shared --prefix=${WHERE_TO_INSTALL}
    ```
    
 * `make && make install`

## Configure, Build, and Install Trilinos
At this time, you should have all external TPL libraries built and
installed into `${WHERE_TO_INSTALL}/lib` and `${WHERE_TO_INSTALL}/include`. You are now ready
to configure the SEACAS portion of the Trilinos cmake build.

The relavant defines for SEACAS are listed below:

```
 -D Trilinos_ENABLE_SEACAS:BOOL=ON
 -D TPL_ENABLE_Netcdf:BOOL=ON 
 -D Netcdf_LIBRARY_DIRS:PATH=${WHERE_TO_INSTALL}/lib 
 -D TPL_Netcdf_INCLUDE_DIRS:PATH=${WHERE_TO_INSTALL}/include
 -D TPL_Netcdf_Enables_Netcdf4:BOOL=ON  (if built with hdf5 libraries which give netcdf-4 capability)
 -D TPL_Netcdf_Enables_PNetcdf:BOOL=ON  (if built with parallel-netcdf which gives parallel I/O capability)
 -D TPL_ENABLE_Matio:BOOL=ON  (set to OFF if not available)
 -D Matio_LIBRARY_DIRS:PATH=${WHERE_TO_INSTALL}/lib
 -D TPL_Matio_INCLUDE_DIRS:PATH=${WHERE_TO_INSTALL}/include
 -D TPL_X11_INCLUDE_DIRS:PATH=/usr/X11R6/include  (SVDI, blot, fastq require X11 includes and libs)
 -D SEACAS_ENABLE_TESTS=ON
 -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L${WHERE_TO_INSTALL}/lib -lpnetcdf -lhdf5_hl -lhdf5 -lz" 
```

If you are using a NetCDF version prior to 4.3 and NetCDF was compiled
with parallel enabled, then you also need to add:

 -D TPL_Netcdf_PARALLEL:BOOL=ON

This is automatically determined for NetCDF-4.3 and later.

## Testing
There are a few unit tests for exodus, and aprepro that can be run via `make test` if you configured with `-D SEACAS_ENABLE_TESTS=ON`.

## Contact information

 Greg Sjaardema  (gsjaardema@gmail.com, gdsjaar@sandia.gov)
      
