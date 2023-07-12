BLAS3 -- KokkosKernels blas3 interfaces
=======================================

gemm
----
.. doxygenfunction:: KokkosBlas::gemm(const execution_space &space, const char transA[], const char transB[], typename AViewType::const_value_type &alpha, const AViewType &A, const BViewType &B, typename CViewType::const_value_type &beta, const CViewType &C)
.. doxygenfunction:: KokkosBlas::gemm(const char transA[], const char transB[], typename AViewType::const_value_type &alpha, const AViewType &A, const BViewType &B, typename CViewType::const_value_type &beta, const CViewType &C)

trmm
----  
.. doxygenfunction:: KokkosBlas::trmm(const execution_space& space, const char side[], const char uplo[], const char trans[], const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B)
.. doxygenfunction:: KokkosBlas::trmm(const char side[], const char uplo[], const char trans[], const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B)

trtri
-----
.. doxygenfunction:: KokkosBlas::trtri

trsm
----
.. doxygenfunction:: KokkosBlas::trsm(const execution_space& space, const char side[], const char uplo[], const char trans[], const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B)
.. doxygenfunction:: KokkosBlas::trsm(const char side[], const char uplo[], const char trans[], const char diag[], typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B)
