//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
//@HEADER

#ifndef KOKKOS_NODE_HELPERS_HPP_
#define KOKKOS_NODE_HELPERS_HPP_

#include <Kokkos_NodeAPIConfigDefs.hpp>
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace Kokkos {

  /** \brief A class to assist in readying buffers via the Node::readyBuffer() method. 
      \ingroup kokkos_node_api
    */
  template <class Node>
  class ReadyBufferHelper {
    public:
      /*! The node via which buffers are being readied. */
      ReadyBufferHelper(RCP<Node> node);

      /*! Destructor. */
      virtual ~ReadyBufferHelper();

      /*! Begin the ready-buffer transaction. */
      void begin();

      /*! Add a const buffer. */
      template <class T>
      const T* addConstBuffer(ArrayRCP<const T> buff);

      /*! Add a non-const buffer. */
      template <class T>
      T* addNonConstBuffer(ArrayRCP<T> buff);

      /*! End the ready-buffer transaction. */
      void end();


    protected:
      RCP<Node> node_;
      Array< ArrayRCP<const char> >  cbufs_;
      Array< ArrayRCP<      char> > ncbufs_;
  };

  template <class Node>
  ReadyBufferHelper<Node>::ReadyBufferHelper(RCP<Node> node)
  : node_(node) {
  }

  template <class Node>
  ReadyBufferHelper<Node>::~ReadyBufferHelper() {
  }

  template <class Node>
  void ReadyBufferHelper<Node>::begin() {
    cbufs_.clear();
    ncbufs_.clear();
  }

  template <class Node>
  template <class T>
  const T* ReadyBufferHelper<Node>::addConstBuffer(ArrayRCP<const T> buff) {
    cbufs_.push_back( arcp_reinterpret_cast<const char>(buff) );
    return buff.get();
  }

  template <class Node>
  template <class T>
  T* ReadyBufferHelper<Node>::addNonConstBuffer(ArrayRCP<T> buff) {
    cbufs_.push_back( arcp_reinterpret_cast<char>(buff) );
    return buff.get();
  }

  template <class Node>
  void ReadyBufferHelper<Node>::end() {
    node_->readyBuffers(cbufs_(), ncbufs_());  
  }

  template <class Node>
  class ArrayOfViewsHelper {
    public:
      template <class T>
      static ArrayRCP<ArrayRCP<T> > getArrayOfNonConstViews(const RCP<Node> &node, ReadWriteOption rw, const ArrayRCP<ArrayRCP<T> > &arrayOfBuffers);
    private:
      /*! Cannot allocate object; all static */
      ArrayOfViewsHelper();
      ~ArrayOfViewsHelper();
  };

  template <class Node>
  template <class T>
  ArrayRCP<ArrayRCP<T> > 
  ArrayOfViewsHelper<Node>::getArrayOfNonConstViews(const RCP<Node> &node, ReadWriteOption rw, const ArrayRCP<ArrayRCP<T> > &arrayOfBuffers) {
    ArrayRCP< ArrayRCP<T> > arrayofviews;
    const size_t numBufs = arrayOfBuffers.size();
    if (numBufs > 0) {
      arrayofviews = arcp< ArrayRCP<T> >(numBufs);
      for (size_t i=0; i < numBufs; ++i) {
        if (arrayOfBuffers[i].size() > 0) {
          arrayofviews[i] = node->template viewBufferNonConst<T>(rw,arrayOfBuffers[i].size(),arrayOfBuffers[i]);
        }
      }
    }
    return arrayofviews;
  }

  /* A trivial implementation of ArrayOfViewsHelper, for CPU-only nodes. */
  template <class Node>
  class ArrayOfViewsHelperTrivialImpl {
    public:
      template <class T>
      static ArrayRCP<ArrayRCP<T> > getArrayOfNonConstViews(const RCP<Node> &node, ReadWriteOption rw, const ArrayRCP<ArrayRCP<T> > &arrayOfBuffers);
    private:
      /*! Cannot allocate object; all static */
      ArrayOfViewsHelperTrivialImpl();
      ~ArrayOfViewsHelperTrivialImpl();
  };

  template <class Node>
  template <class T>
  ArrayRCP<ArrayRCP<T> > 
  ArrayOfViewsHelperTrivialImpl<Node>::getArrayOfNonConstViews(const RCP<Node> &node, ReadWriteOption rw, const ArrayRCP<ArrayRCP<T> > &arrayOfBuffers) {
    (void)node;
    (void)rw;
    return arrayOfBuffers;
  }

} // end of namespace Kokkos


#endif
