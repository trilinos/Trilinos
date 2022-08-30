BLAS2 -- KokkosKernels blas2 interfaces
=======================================

gemv
----
.. doxygenfunction:: KokkosBlas::gemv(const char trans[], typename AViewType::const_value_type &alpha, const AViewType &A, const XViewType &x, typename YViewType::const_value_type &beta, const YViewType &y)
.. doxygenfunction:: KokkosBlas::gemv(const typename AViewType::execution_space &space, const char trans[], typename AViewType::const_value_type &alpha, const AViewType &A, const XViewType &x, typename YViewType::const_value_type &beta, const YViewType &y)
