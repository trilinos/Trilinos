API: Batched
============

.. toctree::
    :maxdepth: 2
    :hidden:

    batched/dense-index

.. toctree::
    :maxdepth: 2
    :hidden:

    batched/sparse-index

Overview
--------

The Kokkos batched interface provides multiple functor-level interfaces for dense linear algebra (DLA) or sparse linear algebra (SLA),
which correspond to Kokkos hierarchical parallelism. Unlike other batched BLAS and LAPACK interfaces, we do not provide a front-level (or subroutine) interface that launches a streaming parallel kernel.
Instead, we provide a functor-level interface that can be used in Kokkos parallel expressions (e.g., parallel for, reduce and scan).
The advantage of this approach is that a user can compose various batched linear algebra kernels, exploiting temporal locality via the functor-level interfaces.

- :doc:`batched dense linear algebra (DLA) <batched/dense-index>`
- :doc:`batched sparse linear algebra (SLA) <batched/sparse-index>`
