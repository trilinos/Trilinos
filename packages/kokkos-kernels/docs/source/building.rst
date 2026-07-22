Building
########

Kokkos Kernels is a stand-alone library in the Kokkos Ecosystem, as well as a package within the `Trilinos Project <https://github.com/trilinos/Trilinos>`_.

Building Kokkos Kernels as a stand-alone library requires CMake.
Using Kokkos Kernels as a package within Trilinos additionally requires the TriBITS build system.

General Requirements
====================

  * Compatible versions of Kokkos and Kokkos Kernels cloned or downloaded from https://github.com/kokkos/kokkos-kernels.
  * Supported compiler and computing hardware - see Kokkos' `README <https://github.com/kokkos/kokkos/blob/master/README.md>`_ for what is currently tested.
  * CUDA builds require use of the ``nvcc_wrapper`` script provided by Kokkos, unless using Clang-CUDA.

Basic steps for building stand-alone Kokkos Kernels:
----------------------------------------------------

#. Create a build directory ``<BUILD_DIR>`` (different from source and install locations).

.. code-block:: shell

  > mkdir <BUILD_DIR>
  cd <BUILD_DIR>

#. Run ``cmake -S ${SOURCE_DIR} <...>`` where ``SOURCE_DIR`` is the location of the Kokkos Kernels source. ``<...>`` is a list of CMake options given as ``-D{OPTION}={VALUE}``
#. Build and install the library, depending on the generator (make is default):

.. code-block:: shell

  > make install
  or
  > ninja install


To use the Ninja build system add ``-G Ninja``.

A full listing of CMake options is given below. Another way to get a list of options and documentation is to the use ``ccmake`` utility.

.. code-block:: shell

  ccmake <SOURCE_DIR>

which brings up a user interface listing all the options, their default values, and associated documentation.


Sample CMake
------------

Below is a list of example CMake configurations.
Kokkos Kernels requires first building Kokkos (or including it as a subproject).
To see a full list of options for building Kokkos, see [BUILD](https://github.com/kokkos/kokkos/blob/master/BUILD.md)

OpenMP backend, g++ compiler, Intel Skylake architecture:
---------------------------------------------------------

First install Kokkos:

.. code-block:: shell
  
  > cmake \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_INSTALL_PREFIX=${HOME}/kokkos-install \
    -DKokkos_ENABLE_OPENMP=ON \
    -DKokkos_ARCH_SKX=ON \
    <KOKKOS_SOURCE>
  > make install

Then build Kokkos Kernels, pointing to Kokkos:

.. code-block:: shell

  > cmake \
    -S <KOKKOS_KERNELS_SOURCE> \
    -B <KOKKOS_KERNELS_BUILD_DIRECTORY> \
    -DKokkos_ROOT=${HOME}/kokkos-install \
    -DCMAKE_CXX_COMPILER=g++ 
  > cmake --build <KOKKOS_KERNELS_BUILD_DIRECTORY> --parallel

Cuda and Serial backends, nvcc_wrapper compiler, Power8 and Volta sm_70 architectures, various compilation flags

First install Kokkos:

.. code-block:: shell

  > cmake \
    -S <KOKKOS_SOURCE> \
    -B <KOKKOS_BUILD_DIR> \
    -DCMAKE_CXX_COMPILER=<KOKKOS_SOURCE>/bin/nvcc_wrapper \
    -DCMAKE_INSTALL_PREFIX=${HOME}/kokkos-install \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ARCH_VOLTA70=ON \
    -DKokkos_ARCH_POWER8=ON
  > cmake --build <KOKKOS_BUILD_DIR> --parallel
  > cmake --install <KOKKOS_BUILD_DIR>

Then build Kokkos Kernels, pointing to Kokkos:

.. code-block:: shell

  > cmake \
    -S <KOKKOS_KERNELS_SOURCE> \
    -B <KOKKOS_KERNELS_BUILD_DIR> \
    -DKokkos_ROOT=${HOME}/kokkos-install \
    -DCMAKE_CXX_COMPILER=${HOME}/kokkos-install/bin/nvcc_wrapper \
  > cmake --build <KOKKOS_KERNELS_BUILD_DIR>

If you wish to enable certain CUDA third-party libraries (TPLs), you can also configure with

.. code-block:: shell

  > cmake \
    -S <KOKKOS_KERNELS_SOURCE> \
    -B <KOKKOS_KERNELS_BUILD_DIR> \
    -DKokkos_ROOT=${HOME}/kokkos-install \
    -DCMAKE_CXX_COMPILER=${HOME}/kokkos-install/bin/nvcc_wrapper \
    -DKokkosKernels_ENABLE_TPL_CUSPARSE=ON

Required
--------

  * CMake >= 3.10
  * Compatible compiler and hardware

#### Trilinos

If building with Trilinos, the same set of CMake options apply. The only difference is you must enable KokkosKernels:

.. code-block:: shell

  > cmake \
    -D Trilinos_ENABLE_KokkosKernels:BOOL=ON \
    ...


## Running tests:

Note, no tests will be available unless ``-DKokkosKernels_ENABLE_TESTS=ON`` is in your cmake command.
To run the tests, simply execute a CMake build and then run:

.. code-block:: shell

  > make test

To limit the tests, one can ``cd`` into either ``unit_test`` or ``perf_test`` and also run ``make test``. To show full detail of all tests, you can run ``ctest --extra-verbose``.
You can filter exactly which tests are run based on regular expressions with

.. code-block:: shell

  > ctest -R <match_string>


Tests are grouped into individual executables.  You can run the executable for one of the enabled backends based on your configuration, for example if OpenMP is enabled:

.. code-block:: shell

  > ./KokkosKernels_UnitTest_OpenMP

To run a specific test in the executable use the ``--gtest_filter`` flag:

.. code-block:: shell

  > ./KokkosKernels_UnitTest_OpenMP --gtest_filter=openmp.dot_double`

