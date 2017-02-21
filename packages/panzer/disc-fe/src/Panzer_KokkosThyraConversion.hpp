#ifndef __Panzer_KokkosThyraConversion_hpp__
#define __Panzer_KokkosThyraConversion_hpp__

#include "Kokkos_View.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_DefaultSpmdVector.hpp"

#include "Epetra_Vector.h"

/** Convert a non-blocked thyra vector into a Kokkos view 
  */

template <typename V>
class Traits { };

template < >
class Traits<Epetra_Vector> {
public:
  typedef Kokkos::View<double*,Kokkos::HostSpace ,Kokkos::MemoryTraits<Kokkos::Unmanaged > > View;
  typedef Thyra::VectorBase<double> ThyraVector;
};

template < >
class Traits<const Epetra_Vector> {
public:
  typedef Kokkos::View<const double*,Kokkos::HostSpace ,Kokkos::MemoryTraits<Kokkos::Unmanaged > > View;
  typedef const Thyra::VectorBase<double> ThyraVector;
};

template <typename VectorType>
typename Traits<VectorType>::View
getKokkosView(typename Traits<VectorType>::ThyraVector & v);

 
template < >
typename Traits<Epetra_Vector>::View
getKokkosView<Epetra_Vector>(typename Traits<Epetra_Vector>::ThyraVector & v)
{
  auto values = Teuchos::ptr_dynamic_cast<Thyra::DefaultSpmdVector<double> >(Teuchos::ptrFromRef(v))->getRCPtr();

  Traits<Epetra_Vector>::View view(values.get(),values.size());
 
  return view;
}

template < >
typename Traits<const Epetra_Vector>::View
getKokkosView<const Epetra_Vector>(typename Traits<const Epetra_Vector>::ThyraVector & v)
{
  auto values = Teuchos::ptr_dynamic_cast<const Thyra::DefaultSpmdVector<double> >(Teuchos::ptrFromRef(v))->getRCPtr();

  Traits<const Epetra_Vector>::View view(values.get(),values.size());
 
  return view;
}

#endif
