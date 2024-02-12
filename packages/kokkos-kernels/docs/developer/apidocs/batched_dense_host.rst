BATCHED -- KokkosKernels batched host-level interfaces
=========================================================

BatchedGemm
-----------
.. doxygenfunction:: KokkosBatched::BatchedGemm(BatchedGemmHandleType *const handle, const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C)
.. doxygenclass:: KokkosBatched::BatchedGemmHandle
    :members: