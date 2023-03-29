BLAS1 -- KokkosKernels blas1 interfaces
=======================================

axpby
-----
.. doxygenfunction:: KokkosBlas::axpby

dot
---
.. doxygenfunction:: KokkosBlas::dot(const RV &, const XMV &, const YMV &, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)
.. doxygenfunction:: KokkosBlas::dot(const XVector &, const YVector &)

fill
----
.. doxygenfunction:: KokkosBlas::fill

mult
----
.. doxygenfunction:: KokkosBlas::mult

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
.. doxygenfunction:: KokkosBlas::reciprocal

scal
----
.. doxygenfunction:: KokkosBlas::scal

sum
---
.. doxygenfunction:: KokkosBlas::sum(const RV &R, const XMV &X, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0)

swap
---
.. doxygenfunction:: KokkosBlas::swap(execution_space const& space, XVector const& X, YVector const& Y)
.. doxygenfunction:: KokkosBlas::swap(XVector const& X, YVector const& Y)

update
------
.. doxygenfunction:: KokkosBlas::update
