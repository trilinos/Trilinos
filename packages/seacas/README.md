# Copyright(C) 1999-2020 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details
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
 * [DataWarehouse](#data_warehouse) -- optional

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

   * Download HDF5 from <http://www.hdfgroup.org/HDF5/release/obtain5.html>

   * untar it, creating a directory will will refer to as `hdf5-X.X.X`

   * `cd` to that directory and enter the command:
     * Serial:
       ```bash
       ./configure --prefix=${WHERE_TO_INSTALL} --enable-shared --enable-production --enable-debug=no --enable-static-exec
       ```

     * Parallel:
       ```bash
       CC=mpicc ./configure --prefix=${WHERE_TO_INSTALL} --enable-shared --enable-production --enable-debug=no --enable-static-exec --enable-parallel
       ```

   * `make && make install`

#### Parallel-NetCDF
  For a parallel build of Trilinos, especially to use the
  auto-decomposition support of the Ioss library, you will need the
  parallel-netcdf library, also known as pnetcdf.

  * Download <http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.6.1.tar.gz>

  * `tar zxvf parallel-netcdf-1.6.1.tar.gz`

  * `cd` to the `parallel-netcdf-1.6.1` directory and enter the command:
    ```bash
    CC=mpicc ./configure --disable-fortran --prefix ${WHERE_TO_INSTALL}
    ```

  * `make && make install`

#### NetCDF
The most recent released version is recommended. For use with Exodus, some local modifications to the netcdf.h include file are required.  See [NetCDF-Mapping.md](NetCDF-Mapping.md) for an explanation of why these modifications are required (or highly recommended)

 * Download the latest netcdf-c release from <http://www.unidata.ucar.edu/downloads/netcdf/index.jsp>

 * `tar zxvf netcdf-4.6.3.tar.gz`  (or whatever the latest version is)

 * If the version is *prior* to 4.5.1, then you need to modify the
   following defines in
   seacas/TPL/netcdf/netcdf-4.6.3/include/netcdf.h.  Versions *4.5.1 or
   later* do not check these limits and can be run unmodified.

    ```c
    #define NC_MAX_DIMS     65536    /* max dimensions per file */
    #define NC_MAX_VARS     524288   /* max variables per file */
    ```

 * `cd netcdf-4.6.3` and enter the command:

    * serial
      ```bash
      CFLAGS="-I${WHERE_TO_INSTALL}/include" \
      CPPFLAGS="-DNDEBUG" LDFLAGS="-L${WHERE_TO_INSTALL}/lib" \
      ./configure --enable-netcdf-4  \
        --disable-fsync --prefix ${WHERE_TO_INSTALL} \
        --disable-dap --disable-v2
      ```

    * parallel
      ```bash
      CC='mpicc' CFLAGS="-I${WHERE_TO_INSTALL}/include" \
      CPPFLAGS="-DNDEBUG" LDFLAGS="-L${WHERE_TO_INSTALL}/lib" \
      ./configure --enable-netcdf-4  --enable-pnetcdf \
        --disable-fsync --prefix ${WHERE_TO_INSTALL} \
        --disable-dap --disable-v2
      ```

 * Check the results of the configure and make sure that the listings
   under features are similar to:

   ```bash
   # Features
   --------
   NetCDF-2 API:        no
   NetCDF-4 API:        yes
   CDF-5 Support:       yes
   HDF5 Support:        yes
   PNetCDF Support:     yes
   NC-4 Parallel Support:       yes
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
   ```bash
   ./autogen.sh
   # The -L is to find the hdf5 library...
   export LDFLAGS="-L${WHERE_TO_INSTALL}/lib"
   ./configure --with-hdf5=${WHERE_TO_INSTALL} --enable-mat73 --enable-shared --prefix=${WHERE_TO_INSTALL}
   ```

 * `make && make install`

#### DataWarehouse
The Data Warehouse is a collection of data management tools that Sandia is currently developing to improve how datasets migrate between memory and storage resources in a distributed system. While the software is currently only available internally to Sandia, there will be an open source release of the tools during Fall of 2017. The Data Warehouse extensions to SEACAS are intended to serve as a placeholder for the upcoming release and are not intended for use by external developers at this time.

The repository [data-warehouse-release](https://gitlab.sandia.gov/nessie-dev/data-warehouse-release) is a superbuild for the Data Warehouse tools and was created to make the build process as easy as possible. It includes the files INSTALL.md and INSTALL_TPL.md which contain instructions for building Data Warehouse and it's TPLs: Boost, googletest, libfabric, and libhio (optional). These builds are straightforward so a "runconfigure.sh" script is left to the end-user. Note that it's possible to supply your own build of these tools. Following the SEACAS pattern for building TPLs:

 * `cd TPL`
 * `git clone git@gitlab.sandia.gov:nessie-dev/data-warehouse-release.git`
 * follow the instructions in INSTALL_TPL.md and then INSTALL.md
 * consider installing data-warehouse-release and it's TPLs to the directory set in  the $ACCESS env. var.

To build SEACAS with an installation of data-warehouse-release and it's TPLs, add the following lines to the list of cmake command arguments as found in the top-level `cmake-config` file. Where DataWarehouse_PATH in environment variable that contains the path to the top-level install directory for the Data Warehouse, and HAVE_DATA_WAREHOUSE={ON|OFF} is a variable defined in `cmake-config`.

```bash
    -DTPL_ENABLE_DATAWAREHOUSE:BOOL=${HAVE_DATAWAREHOUSE}           \
    -DDataWarehouse_LIBRARY_DIRS:PATH=${DataWarehouse_PATH}/lib     \
    -DDataWarehouse_INCLUDE_DIRS:PATH=${DataWarehouse_PATH}/include \
```

## Configure, Build, and Install Trilinos
At this time, you should have all external TPL libraries built and
installed into `${WHERE_TO_INSTALL}/lib` and `${WHERE_TO_INSTALL}/include`. You are now ready
to configure the SEACAS portion of the Trilinos cmake build.

The relevant defines for SEACAS are:
```bash
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

## Testing
There are a few unit tests for exodus, and aprepro that can be run via `make test` if you configured with `-D SEACAS_ENABLE_TESTS=ON`.

## Contact information

 Greg Sjaardema  (<mailto:gsjaardema@gmail.com>, <mailto:gdsjaar@sandia.gov>)
