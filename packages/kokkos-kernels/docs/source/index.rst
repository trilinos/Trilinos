.. role:: raw-html-m2r(raw)
   :format: html

.. include:: mydefs.rst

Kokkos Kernels: Portable Math Kernels
=====================================

.. admonition:: :medium:`C++ Performance Portable Math Library`
    :class: important

    :medium:`Kokkos Kernels provides math kernels for dense, batched and sparse linear algebra as well as graph computations. Our interfaces and implementations are portable to major HPC architectures and based on the Kokkos programming model. The library supports vendor implementation when available and more performant than our native algorithms.`

The `Kokkos Ecosystem <https://github.com/kokkos>`_ includes:

.. list-table::
   :widths: 30 50 20
   :header-rows: 1
   :align: left

   * - Name
     - Info
     -

   * - ``kokkos``
     - Programming Model - Parallel Execution and Memory Abstraction
     - `GitHub link <https://github.com/kokkos/kokkos>`__

   * - ``kokkos-kernels``
     - (this library) Sparse, dense, batched math kernels
     - `GitHub link <https://github.com/kokkos/kokkos-kernels>`__

   * - ``kokkos-tools``
     - Profiling and debugging tools
     - `GitHub link <https://github.com/kokkos/kokkos-tools>`__

   * - ``pykokkos``
     - Provides Python bindings to the Kokkos performance portable parallel programming.
     - `GitHub link <https://github.com/kokkos/pykokkos>`__

   * - ``kokkos-remote-spaces``
     - Shared memory semantics across multiple processes
     - `GitHub link <https://github.com/kokkos/kokkos-remote-spaces>`__

   * - ``kokkos-resilience``
     - Resilience and Checkpointing Extensions for Kokkos
     - `GitHub link <https://github.com/kokkos/kokkos-resilience>`__

Questions?
----------

Find us on Slack: https://kokkosteam.slack.com #kokkos-kernels or
open an issue on `github <https://github.com/kokkos/kokkos-kernels/issues>`_.

Website Content
---------------

.. toctree::
   :maxdepth: 1

   quick_start
   requirements
   building
   cmake-keywords
   ./API/common-index
   ./API/blas-index
   ./API/batched-index
   ./API/lapack-index
   ./API/sparse-index
   ./API/graph-index
   ./known_issues
   deprecation_page
   Tutorials <https://github.com/kokkos/kokkos-tutorials>
   GitHub Repo <https://github.com/kokkos/kokkos-kernels>
   Contributing
   publications
   license
