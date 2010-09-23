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
