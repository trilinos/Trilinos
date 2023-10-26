BLAS1 -- KokkosKernels blas1 interfaces
=======================================

abs
---
.. doxygenfunction:: KokkosBlas::abs(const execution_space& space, const RMV& R, const XMV& X)
.. doxygenfunction:: KokkosBlas::abs(const RMV& R, const XMV& X)

axpby
-----
.. doxygenfunction:: KokkosBlas::axpby(const execution_space& space, const AV& a, const XMV& X, const BV& b, const YMV& Y)
.. doxygenfunction:: KokkosBlas::axpby(const AV& a, const XMV& X, const BV& b, const YMV& Y)

dot
---
.. doxygenfunction:: KokkosBlas::dot(const RV &, const XMV &, const YMV &, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)
.. doxygenfunction:: KokkosBlas::dot(const XVector &, const YVector &)

fill
----
.. doxygenfunction:: KokkosBlas::fill(const execution_space& space, const XMV& X, const typename XMV::non_const_value_type& val)
.. doxygenfunction:: KokkosBlas::fill(const XMV& X, const typename XMV::non_const_value_type& val)

mult
----
.. doxygenfunction:: KokkosBlas::mult(const execution_space& space, typename YMV::const_value_type& gamma, const YMV& Y, typename AV::const_value_type& alpha, const AV& A, const XMV& X)
.. doxygenfunction:: KokkosBlas::mult(typename YMV::const_value_type& gamma, const YMV& Y, typename AV::const_value_type& alpha, const AV& A, const XMV& X)

nrm1
----
.. doxygenfunction:: KokkosBlas::nrm1(const RV &, const XMV &, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)
.. doxygenfunction:: KokkosBlas::nrm1(const XVector &)

nrm2
----
.. doxygenfunction:: KokkosBlas::nrm2(const RV &R, const XMV &X, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)
.. doxygenfunction:: KokkosBlas::nrm2(const XVector &x)

nrm2w
-----
.. doxygenfunction:: KokkosBlas::nrm2w(const RV &R, const XMV &X, const XMV &W, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)
.. doxygenfunction:: KokkosBlas::nrm2w(const XVector &x, const XVector &w)

nrminf
------
.. doxygenfunction:: KokkosBlas::nrminf(const RV &R, const XMV &X, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)
.. doxygenfunction:: KokkosBlas::nrminf(const XVector &x)

reciprocal
----------
.. doxygenfunction:: KokkosBlas::reciprocal(const execution_space& space, const RMV& R, const XMV& X)
.. doxygenfunction:: KokkosBlas::reciprocal(const RMV& R, const XMV& X)

scal
----
.. doxygenfunction:: KokkosBlas::scal(const execution_space& space, const RMV& R, const AV& a, const XMV& X)
.. doxygenfunction:: KokkosBlas::scal(const RMV& R, const AV& a, const XMV& X)

sum
---
.. doxygenfunction:: KokkosBlas::sum(const RV &R, const XMV &X, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)

swap
----
.. doxygenfunction:: KokkosBlas::swap(execution_space const&, XVector const&, YVector const&)
.. doxygenfunction:: KokkosBlas::swap(const XVector&, const YVector&)

update
------
.. doxygenfunction:: KokkosBlas::update(const execution_space& space, const typename XMV::non_const_value_type& alpha, const XMV& X, const typename YMV::non_const_value_type& beta, const YMV& Y, const typename ZMV::non_const_value_type& gamma, const ZMV& Z)
.. doxygenfunction:: KokkosBlas::update(const typename XMV::non_const_value_type& alpha, const XMV& X, const typename YMV::non_const_value_type& beta, const YMV& Y, const typename ZMV::non_const_value_type& gamma, const ZMV& Z)
