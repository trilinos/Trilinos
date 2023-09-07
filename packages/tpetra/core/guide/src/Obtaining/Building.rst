
.. _building_tpetra:

Building Tpetra
###############

Requirements
============

Tpetra requires a C++11 compatible compiler for compilation. The minimum required version of compilers are GCC (4.7.2 and later), Intel (13 and later), and clang (3.5 and later).

Dependencies
============

Dependencies for Tpetra are listed below.  Certain dependencies are optional, whereas others are required. Furthermore, Tpetra’s tests depend on certain libraries that are not required if you only want to link against the Tpetra library and do not want to compile its tests. Additionally, some functionality in Tpetra may depend on other Trilinos packages that may require additional dependencies. We refer to the documentation of those packages for a full list of dependencies.

.. tabularcolumns:: |l|c|c|c|c|

+----------------+----------------------+----------------------+
| **Dependency** |    **Library**       |       **Testing**    |
+                +-----------+----------+-----------+----------+
|                | Required  | Optional | Required  | Optional |
+----------------+-----------+----------+-----------+----------+
| Teuchos        |  x        |          |     x     |          |
+----------------+-----------+----------+-----------+----------+
| Kokkos         |  x        |          |     x     |          |
+----------------+-----------+----------+-----------+----------+

Configuration
=============

The preferred way to configure and build Tpetra is to do so outside of the source directory.  For example, if the Trilinos source code is located in ``/path/to/Trilinos``, Tpetra can be built in ``/path/to/TrilinosBuild``.  Here we provide a sample configure script that should be executed in ``/path/to/TrilinosBuild`` that will enable Tpetra and all of its optional dependencies:

.. code-block:: sh

   export TRILINOS_HOME=/path/to/Trilinos/
   cmake -D BUILD_SHARED_LIBS:BOOL=ON \
         -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
         -D CMAKE_CXX_FLAGS:STRING="-g" \
         -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
         -D Trilinos_ENABLE_TESTS:BOOL=OFF \
         -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
         -D Trilinos_ENABLE_Tpetra:BOOL=ON \
         -D Tpetra_ENABLE_TESTS:STRING=ON \
         -D Tpetra_ENABLE_EXAMPLES:STRING=ON \
         -D TPL_ENABLE_BLAS:BOOL=ON \
         -D TPL_ENABLE_MPI:BOOL=ON \
         ${TRILINOS_HOME}

More configure examples can be found in ``Trilinos/sampleScripts``. For more information on configuring, see the `TRILINOS Cmake Quickstart guide <https://trilinos.org/oldsite/TrilinosBuildQuickRef.html>`_.

Enabling Templated Types
========================

.. note::

   The following instructions are a summary of the Tpetra `CMake specificaton <https://github.com/trilinos/Trilinos/blob/develop/packages/tpetra/CMakeLists.txt>`_ which should be consulted for more details.

For

* ``Scalar`` types,
* (``LocalOrdinal``, ``GlobalOrdinal``) type pairs, and
* ``Node`` types,

Tpetra defines macros which say

1. whether Tpetra instantiates (if explicit template instantiations ``ETI`` is enabled) and/or tests its objects for that type (or type pair)

2. whether the type (or type pair) is available

If ``ETI`` is enabled, both are identical.  If not, the settings of the latter depend on both Kokkos settings and Trilinos' settings.

For all templated types, the CMake options have the pattern: ``Tpetra_INST_<TYPE>``, i.e., to enable a templated type add the following to the Tpetra configuration:

.. code-block:: sh

   -D Tpetra_INST_<TYPE>:BOOL=ON

+---------------+------------------------------+
| Template type | CMake option ``TYPE``        |
+===============+==============================+
|               |   ``DOUBLE``,                |
|               |   ``FLOAT``,                 |
| ``Scalar``    |   ``COMPLEX_DOUBLE``,        |
|               |   ``COMPLEX_FLOAT``,         |
+---------------+------------------------------+
|               |   ``INT_INT``,               |
|               |   ``INT_UNSIGNED``,          |
| ``Ordinal``   |   ``INT_LONG``,              |
|               |   ``INT_LONG_LONG``,         |
|               |   ``INT_UNSIGNED_LONG_LONG`` |
+---------------+------------------------------+
|               |   ``SERIAL``,                |
|               |   ``OPENMP``,                |
| ``Node``      |   ``PTHREAD``,               |
|               |   ``CUDA``                   |
+---------------+------------------------------+

.. note::

   Concerning (1), Tpetra also adds enabled ``GlobalOrdinal`` types internally to the ``Scalar`` types for instantiation of ``MultiVector``.

.. note::

   Concerning (2), some of the types can only be enabled if certain requirements are fulfilled.


Kokkos Execution Space Types
----------------------------

If ``ETI`` is ``ON``, "using an execution space ``ExecutionSpace``" means that Tpetra objects with ``Node = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<ExecutionSpace>`` will

1. get instantiated explicitly for that ``Node`` type, and

2. get tested for that ``Node`` type, if their test is also templated on ``Node`` type.

If ``ETI`` is ``OFF``, 1 no longer holds, but 2 still holds.

Tpetra uses exactly one Kokkos execution space by default, whether ``ETI`` is ``ON`` or ``OFF``.  This keeps build times for tests down in the non-``ETI`` case, and library sizes small in the ``ETI`` case.

The best option for building Tpetra is to enable the execution space type at the Trilinos level (OpenMP, CUDA, etc.) and let the defaults propagate into Tpetra. For example:

- Enabling CUDA (by using NVCC and ``nvcc_wrapper``) makes CUDA Tpetra's default execution space.
- Enabling OpenMP (``Trilinos_ENABLE_OpenMP:BOOL=ON``), but not enabling CUDA, makes OpenMP Tpetra's default execution space.
- The Pthreads (``Kokkos::Threads``) back-end is a special case; it does not get enabled by default. This avoids surprises, because Trilinos enables its Pthreads TPL by default as long as it can detect it. Users may set the CMake option ``Kokkos_ENABLE_THREADS:BOOL=ON`` to enable use of Pthreads in Tpetra, and to make it default.

