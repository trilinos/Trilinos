//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

#ifndef KOKKOS_NODE_HELPERS_HPP_
#define KOKKOS_NODE_HELPERS_HPP_

#include <Kokkos_NodeAPIConfigDefs.hpp>
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace Kokkos {

  /// \class ReadyBufferHelper
  /// \brief A class to assist in readying buffers via the Node::readyBuffers() method.
  /// \ingroup kokkos_node_api
  ///
  /// \tparam Node The Kokkos Node type.
  ///
  /// Kokkos asks that you call the Node's readyBuffers() method on
  /// all compute buffers that you plan to pass into a kernel.  In a
  /// debug build, this checks to make sure that the ArrayRCP objects
  /// are really compute buffers.  ReadyBufferHelper provides a
  /// transaction-like interface to help you with this.  The
  /// ReadyBufferHelper also keeps the ArrayRCP objects for you during
  /// its lifetime.  This prevents them from falling out of scope and
  /// getting deleted, thus ensuring that when you extract the raw
  /// pointer to give to the Kokkos kernel, the pointer will be valid
  /// throughout the lifetime of the kernel.
  ///
  /// Here is an example of how to use ReadyBufferHelper with a simple
  /// Kokkos kernel.  Suppose you have a \c parallel_for kernel Op
  /// that takes a const array and a nonconst array:
  /// \code
  /// template<class T>
  /// class Op {
  /// public:
  ///   Op (const T* x, T* y, size_t n) :
  ///     x_ (x), y_ (y), n_ (n) {}
  ///
  ///   inline KERNEL_PREFIX void execute (size_t i) {
  ///     y_[i] += 2.0 * x_[i];
  ///   }
  /// private:
  ///   const T* x_;
  ///   T* y_;
  ///   size_t n_;
  /// };
  /// \endcode
  /// Here's how you would use ReadyBufferHelper when invoking the
  /// kernel in parallel:
  /// \code
  /// using Teuchos::ArrayRCP;
  /// using Teuchos::RCP;
  ///
  /// // Suppose node_type is the Kokkos Node type.
  /// RCP<node_type> node = ...;
  ///
  /// // x and y are compute buffers that you got elsewhere.
  /// ArrayRCP<const T> x = ...;
  /// ArrayRCP<T> y = ...;
  /// const size_t n = x.size ();
  ///
  /// ReadyBufferHelper<node_type> rbh (node);
  /// rbh.begin (); // Start adding compute buffers
  ///
  /// // Create the Kokkos kernel "work-data pair" (Kokkos' phrase for
  /// // what other programming languages call "closure").  When adding
  /// // a const buffer (ArrayRCP<const T>), omit the "const" when naming
  /// // the type in the template method call.
  /// //
  /// // If T is a template parameter in this scope, you need the "template"
  /// // keyword when invoking ReadyBufferHelper template methods.  You don't
  /// // need it if T is a concrete type in this scope.
  /// Op wdp (rbh.template addConstBuffer<T> (x),
  ///         rbh.template addNonConstBuffer<T> (y),
  ///         n);
  /// rbh.end (); // Done adding compute buffers
  ///
  /// // Invoke the kernel in parallel over the range 0, ..., n-1.
  /// // The ReadyBufferHelper will protect x and y from falling out
  /// // of scope.
  /// //
  /// // If T is a template parameter in this scope, you need the "template"
  /// // keyword.  You don't need it if T is a concrete type in this scope.
  /// node->template parallel_for<Op<T> > (0, x.size (), wdp);
  /// \endcode
  template <class Node>
  class ReadyBufferHelper {
  public:
    /*! The node via which buffers are being readied. */
    ReadyBufferHelper(RCP<Node> node);

    /*! Destructor. */
    virtual ~ReadyBufferHelper();

    //! Tell this object that you are ready to start adding buffers.
    void begin();

    //! Add a const buffer, and return its raw pointer.
    template <class T>
    const T* addConstBuffer(ArrayRCP<const T> buff);

    //! Add a non-const buffer, and return its raw pointer.
    template <class T>
    T* addNonConstBuffer(ArrayRCP<T> buff);

    //! Tell this object that you are done adding buffers.
    void end();

  protected:
    //! The Kokkos Node instance for which to add buffers.
    RCP<Node> node_;
    //! List of compute buffers that were added with addConstBuffer().
    Array<ArrayRCP<const char> >  cbufs_;
    //! List of compute buffers that were added with addNonConstBuffer().
    Array<ArrayRCP<      char> > ncbufs_;
  };

  template <class Node>
  ReadyBufferHelper<Node>::ReadyBufferHelper (RCP<Node> node)
    : node_(node)
  {}

  template <class Node>
  ReadyBufferHelper<Node>::~ReadyBufferHelper()
  {}

  template <class Node>
  void ReadyBufferHelper<Node>::begin() {
    cbufs_.clear ();
    ncbufs_.clear ();
  }

  template <class Node>
  template <class T>
  const T* ReadyBufferHelper<Node>::addConstBuffer (ArrayRCP<const T> buff) {
    cbufs_.push_back (arcp_reinterpret_cast<const char> (buff));
    return buff.get ();
  }

  template <class Node>
  template <class T>
  T* ReadyBufferHelper<Node>::addNonConstBuffer (ArrayRCP<T> buff) {
    cbufs_.push_back (arcp_reinterpret_cast<char> (buff));
    return buff.get ();
  }

  template <class Node>
  void ReadyBufferHelper<Node>::end () {
    node_->readyBuffers (cbufs_ (), ncbufs_ ());
  }

  /// \class ArrayOfViewsHelper
  /// \brief Helper class for getting an array of views.
  /// \ingroup kokkos_node_api
  ///
  /// \tparam Node The Kokkos Node type.
  ///
  /// The class method getArrayOfNonConstView takes an array of
  /// (device) buffers, and returns an array of (host) views of those
  /// buffers.  The host views may be either write-only (meaning that
  /// their contents on the host before being written are undefined),
  /// or read-and-write (meaning that they have valid contents on the
  /// host).
  ///
  /// All methods of ArrayOfViewsHelper are class (i.e., static)
  /// methods.  It doesn't make sense to make an instance of
  /// ArrayOfViewsHelper.  We forbid this syntactically by declaring
  /// the (unimplemented) constructor private.
  template <class Node>
  class ArrayOfViewsHelper {
  public:
    /// \brief Invoke the Node's viewBufferNonConst() method to get an array of views.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    ///
    /// \param rw [in] Whether to create read-and-write views, or
    ///   write-only views.  Read-and-write views may be safely read
    ///   before they are written.  Write-only views have undefined
    ///   contents before they are written.  In both cases, views are
    ///   to be read and written on the host.  If the device has a
    ///   separate memory space, then each view is copied back to the
    ///   device (namely, to their corresponding buffers in
    ///   arrayOfBuffers) once its reference count falls to zero.
    ///
    /// \param arrayOfBuffers [in/out] The array of buffers for which
    ///   to create views.
    template <class T>
    static ArrayRCP<ArrayRCP<T> >
    getArrayOfNonConstViews (const RCP<Node> &node,
                             ReadWriteOption rw,
                             const ArrayRCP<ArrayRCP<T> > &arrayOfBuffers);
  private:
    //! Constructor is private and undeclared, so you can't call it.
    ArrayOfViewsHelper();
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

  //! A trivial implementation of ArrayOfViewsHelper, for CPU-only nodes.
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
