#ifndef KOKKOS_NODE_HELPERS_HPP_
#define KOKKOS_NODE_HELPERS_HPP_

#include <Kokkos_NodeAPIConfigDefs.hpp>
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace Kokkos {

  /*! A class to assist in readying buffers via the Node::readyBuffer() method. */
  template <class Node>
  class ReadyBufferHelper {
    public:
      /*! The node via which buffers are being readied. */
      ReadyBufferHelper(Teuchos::RCP<Node> node);

      /*! Destructor. */
      virtual ~ReadyBufferHelper();

      /*! Begin the ready-buffer transaction. */
      void begin();

      /*! Add a const buffer. */
      template <class T>
      const T* addConstBuffer(Teuchos::ArrayRCP<const T> buff);

      /*! Add a non-const buffer. */
      template <class T>
      T* addNonConstBuffer(Teuchos::ArrayRCP<T> buff);

      /*! End the ready-buffer transaction. */
      void end();


    protected:
      Teuchos::RCP<Node> node_;
      Teuchos::Array< Teuchos::ArrayRCP<const char> >  cbufs_;
      Teuchos::Array< Teuchos::ArrayRCP<      char> > ncbufs_;
  };

  template <class Node>
  ReadyBufferHelper<Node>::ReadyBufferHelper(Teuchos::RCP<Node> node)
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
  const T* ReadyBufferHelper<Node>::addConstBuffer(Teuchos::ArrayRCP<const T> buff) {
    cbufs_.push_back( Teuchos::arcp_reinterpret_cast<const char>(buff) );
    return buff.get();
  }

  template <class Node>
  template <class T>
  T* ReadyBufferHelper<Node>::addNonConstBuffer(Teuchos::ArrayRCP<T> buff) {
    cbufs_.push_back( Teuchos::arcp_reinterpret_cast<char>(buff) );
    return buff.get();
  }

  template <class Node>
  void ReadyBufferHelper<Node>::end() {
    node_->readyBuffers(cbufs_(), ncbufs_());  
  }

  /*! A class to efficient get an array of views */
  template <class Node>
  class ArrayOfViewsHelper {
    public:
      template <class T>
      static Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > getArrayOfNonConstViews(const Teuchos::RCP<Node> &node, ReadWriteOption rw, const Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > &arrayOfBuffers);
    private:
      /*! Cannot allocate object; all static */
      ArrayOfViewsHelper();
      ~ArrayOfViewsHelper();
  };

  template <class Node>
  template <class T>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > 
  ArrayOfViewsHelper<Node>::getArrayOfNonConstViews(const Teuchos::RCP<Node> &node, ReadWriteOption rw, const Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > &arrayOfBuffers) {
    Teuchos::ArrayRCP< Teuchos::ArrayRCP<T> > arrayofviews;
    const size_t numBufs = arrayOfBuffers.size();
    if (numBufs > 0) {
      arrayofviews = Teuchos::arcp< Teuchos::ArrayRCP<T> >(numBufs);
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
      static Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > getArrayOfNonConstViews(const Teuchos::RCP<Node> &node, ReadWriteOption rw, const Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > &arrayOfBuffers);
    private:
      /*! Cannot allocate object; all static */
      ArrayOfViewsHelperTrivialImpl();
      ~ArrayOfViewsHelperTrivialImpl();
  };

  template <class Node>
  template <class T>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > 
  ArrayOfViewsHelperTrivialImpl<Node>::getArrayOfNonConstViews(const Teuchos::RCP<Node> &node, ReadWriteOption rw, const Teuchos::ArrayRCP<Teuchos::ArrayRCP<T> > &arrayOfBuffers) {
    (void)node;
    (void)rw;
    return arrayOfBuffers;
  }

} // end of namespace Kokkos


#endif
