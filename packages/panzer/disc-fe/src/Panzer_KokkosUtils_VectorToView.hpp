// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_KokkosThyraConversion_hpp__
#define __Panzer_KokkosThyraConversion_hpp__

#include "Kokkos_Core.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_DefaultSpmdVector.hpp"

#include "Epetra_Vector.h"

namespace panzer {
namespace kokkos_utils {

/** Convert a non-blocked thyra vector into a Kokkos view 
  */

template <typename V>
class VectorToViewTraits { };

template < >
class VectorToViewTraits<Epetra_Vector> {
public:
  typedef Kokkos::View<double*,Kokkos::HostSpace ,Kokkos::MemoryTraits<Kokkos::Unmanaged > > View;
  typedef Thyra::VectorBase<double> ThyraVector;
};

template < >
class VectorToViewTraits<const Epetra_Vector> {
public:
  typedef Kokkos::View<const double*,Kokkos::HostSpace ,Kokkos::MemoryTraits<Kokkos::Unmanaged > > View;
  typedef const Thyra::VectorBase<double> ThyraVector;
};

template <typename VectorType>
inline
typename VectorToViewTraits<VectorType>::View
getView(typename VectorToViewTraits<VectorType>::ThyraVector & v);
 
template < >
inline
typename VectorToViewTraits<Epetra_Vector>::View
getView<Epetra_Vector>(typename VectorToViewTraits<Epetra_Vector>::ThyraVector & v)
{
  auto values = Teuchos::ptr_dynamic_cast<Thyra::DefaultSpmdVector<double> >(Teuchos::ptrFromRef(v))->getRCPtr();

  VectorToViewTraits<Epetra_Vector>::View view(values.get(),values.size());
 
  return view;
}

template < >
inline
typename VectorToViewTraits<const Epetra_Vector>::View
getView<const Epetra_Vector>(typename VectorToViewTraits<const Epetra_Vector>::ThyraVector & v)
{
  auto values = Teuchos::ptr_dynamic_cast<const Thyra::DefaultSpmdVector<double> >(Teuchos::ptrFromRef(v))->getRCPtr();

  VectorToViewTraits<const Epetra_Vector>::View view(values.get(),values.size());
 
  return view;
}

}
}

#endif
