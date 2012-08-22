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

#ifndef TPETRA_VIEWACCEPTER_HPP
#define TPETRA_VIEWACCEPTER_HPP

#include <Teuchos_ArrayRCP.hpp>
#include "Tpetra_ConfigDefs.hpp"

#include <Kokkos_SerialNode.hpp>
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include <Kokkos_TBBNode.hpp>
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include <Kokkos_TPINode.hpp>
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#include <Kokkos_OpenMPNode.hpp>
#endif

namespace Tpetra {

  namespace details {

  //! \brief ViewAccepter provides support for compile-time detection based on Node type, of support for accepting user views.
  template <class Node>
  class ViewAccepter {
    private:
      // no support for object instantiation
      ViewAccepter();
    public:
      template <class T>
      inline static ArrayRCP<T> acceptView(const ArrayRCP<T> &view) {
        Node::this_node_type_does_not_support_object_construction_from_user_views();
        return null;
      }
  };

  //! Implementation 
  class ViewAccepterSupportedNode {
    private:
      // no support for object instantiation
      ViewAccepterSupportedNode();
    public:
      template <class T>
      inline static ArrayRCP<T> acceptView(const ArrayRCP<T> &view) {
        return view;
      }
  };

  template <>
  class ViewAccepter<Kokkos::SerialNode> : public ViewAccepterSupportedNode {};
#ifdef HAVE_KOKKOSCLASSIC_TBB
  template <>
  class ViewAccepter<Kokkos::TBBNode> : public ViewAccepterSupportedNode {};
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  template <>
  class ViewAccepter<Kokkos::TPINode> : public ViewAccepterSupportedNode {};
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  template <>
  class ViewAccepter<Kokkos::OpenMPNode> : public ViewAccepterSupportedNode {};
#endif

  } // end of namespace Tpetra::details

} // end of namespace Tpetra

#endif
