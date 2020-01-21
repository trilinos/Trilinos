#ifndef BLASWRAPPER_IAMAX_HPP_
#define BLASWRAPPER_IAMAX_HPP_

#include<BlasWrapper_iamax_spec.hpp>
#include<KokkosKernels_helpers.hpp>

namespace BlasWrapper {

/// \brief Return the (smallest) index of the element of the maximum magnitude of the vector x. 
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
///
/// \return The (smallest) index of the element of the maximum magnitude; a single value.
///         Note: Returned index is 1-based for compatibility with Fortran.    
template<class XVector>
typename XVector::size_type iamax (const XVector& x)
{
  static_assert (Kokkos::Impl::is_view<XVector>::value,
                 "BlasWrapper::iamax: XVector must be a Kokkos::View.");
  static_assert (XVector::rank == 1, "BlasWrapper::iamax: "
                 "Both Vector inputs must have rank 1.");

  typedef typename XVector::size_type index_type;

  typedef Kokkos::View<typename XVector::const_value_type*,
    typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
    typename XVector::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XVector_Internal;

  typedef Kokkos::View<index_type,
    Kokkos::LayoutLeft,
    Kokkos::HostSpace,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > RVector_Internal;

  index_type result;
  RVector_Internal R = RVector_Internal(&result);
  XVector_Internal X = x;

  Impl::Iamax<RVector_Internal,XVector_Internal>::iamax (R,X);
  Kokkos::fence();
  return result;
}

}

#endif // BLASWRAPPER_IAMAX_HPP_
