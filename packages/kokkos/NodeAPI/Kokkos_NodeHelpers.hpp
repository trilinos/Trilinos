#ifndef KOKKOS_NODE_HELPERS_HPP_
#define KOKKOS_NODE_HELPERS_HPP_

#include <Kokkos_ConfigDefs.hpp>
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace Kokkos {

  /*! A class to assist in readying buffers via the Node::readyBuffer() method. */
  template <class Node>
  class ReadyBufferHelper {
    public:
      /*! The node via which buffers are being readied. */
      ReadyBufferHelper(Node &node);

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
      Node &node_;
      Teuchos::Array< Teuchos::ArrayRCP<const char> >  cbufs_;
      Teuchos::Array< Teuchos::ArrayRCP<      char> > ncbufs_;
  };

  template <class Node>
  ReadyBufferHelper<Node>::ReadyBufferHelper(Node &node)
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
    node_.readyBuffers(cbufs_(), ncbufs_());  
  }

} // end of namespace Kokkos


#endif
