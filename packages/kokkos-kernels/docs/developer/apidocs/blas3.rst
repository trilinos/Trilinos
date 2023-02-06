BLAS3 -- KokkosKernels blas3 interfaces
=======================================

gemm
----
.. doxygenfunction:: KokkosBlas::gemm(const char transA[], const char transB[], typename AViewType::const_value_type &alpha, const AViewType &A, const BViewType &B, typename CViewType::const_value_type &beta, const CViewType &C)
.. doxygenfunction:: KokkosBlas::gemm(const typename CViewType::execution_space &space, const char transA[], const char transB[], typename AViewType::const_value_type &alpha, const AViewType &A, const BViewType &B, typename CViewType::const_value_type &beta, const CViewType &C)
