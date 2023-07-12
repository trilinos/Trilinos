BLAS2 -- KokkosKernels blas2 interfaces
=======================================

gemv
----
.. doxygenfunction:: KokkosBlas::gemv(const char trans[], typename AViewType::const_value_type &alpha, const AViewType &A, const XViewType &x, typename YViewType::const_value_type &beta, const YViewType &y)
.. doxygenfunction:: KokkosBlas::gemv(const execution_space &space, const char trans[], typename AViewType::const_value_type &alpha, const AViewType &A, const XViewType &x, typename YViewType::const_value_type &beta, const YViewType &y)

ger
----
.. doxygenfunction:: KokkosBlas::ger(const ExecutionSpace& space, const char trans[], const typename AViewType::const_value_type& alpha, const XViewType& x, const YViewType& y, const AViewType& A)
.. doxygenfunction:: KokkosBlas::ger(const char trans[], const typename AViewType::const_value_type& alpha, const XViewType& x, const YViewType& y, const AViewType& A)
