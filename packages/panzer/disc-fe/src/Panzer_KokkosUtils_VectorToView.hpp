// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
