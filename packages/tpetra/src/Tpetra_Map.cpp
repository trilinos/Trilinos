/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#include "Tpetra_Map.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "Tpetra_ETIHelperMacros.h"
#include "Tpetra_Map_def.hpp"

namespace Tpetra {

  // for default node
  template Teuchos::RCP< const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> >
    createNonContigMap<int,int>(const Teuchos::ArrayView<const int> &elementList, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
  template Teuchos::RCP< const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> >
    createContigMap<int,int>(size_t numElements, size_t numLocalElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
  template Teuchos::RCP< const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> >
    createLocalMap<int,int>(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
  template Teuchos::RCP< const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> >
    createUniformContigMap<int,int>(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
#ifdef HAVE_TPETRA_INST_INT_LONG
  template Teuchos::RCP< const Map<int,long,Kokkos::DefaultNode::DefaultNodeType> >
    createNonContigMap<int,long>(const Teuchos::ArrayView<const long> &elementList, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
  template Teuchos::RCP< const Map<int,long,Kokkos::DefaultNode::DefaultNodeType> >
    createContigMap<int,long>(size_t numElements, size_t numLocalElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
  template Teuchos::RCP< const Map<int,long,Kokkos::DefaultNode::DefaultNodeType> >
    createLocalMap<int,long>(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
  template Teuchos::RCP< const Map<int,long,Kokkos::DefaultNode::DefaultNodeType> >
    createUniformContigMap<int,long>(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);
#endif

  TPETRA_INSTANTIATE_LGN(TPETRA_MAP_INSTANT)

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
