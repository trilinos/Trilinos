#ifndef XPETRA_SERIALPLATFORM_HPP
#define XPETRA_SERIALPLATFORM_HPP

#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"

namespace Xpetra {

  //! \brief A implementation of the Platform class for serial platforms.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template<class Node=Kokkos::DefaultNode::DefaultNodeType>
  class SerialPlatform : public Teuchos::Describable {
  public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
    typedef Node NodeType;
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    explicit SerialPlatform(const Teuchos::RCP<Node> &node);

    //! Destructor
    ~SerialPlatform();

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    const Teuchos::RCP< const Teuchos::SerialComm<int> > getComm() const;

    //! Get Get a node for parallel computation.
    const Teuchos::RCP<Node> getNode() const;

    //@}
  private:
    SerialPlatform(const SerialPlatform<Node> &platform);

  protected: 
    //! Teuchos::Comm object instantiated for the platform.
    Teuchos::RCP<const Teuchos::SerialComm<int> > comm_;
    //! Node object instantiated for the platform.
    Teuchos::RCP<Node> node_;
  };

  template<class Node>
  SerialPlatform<Node>::SerialPlatform(const Teuchos::RCP<Node> &node) 
    : node_(node) { 
    comm_ = Teuchos::rcp(new Teuchos::SerialComm<int>() );
  }

  template<class Node>
  SerialPlatform<Node>::~SerialPlatform() {  }

  template<class Node>
  const Teuchos::RCP< const Teuchos::SerialComm<int> >
  SerialPlatform<Node>::getComm() const { 
    return comm_;
  }

  template<class Node>
  const Teuchos::RCP< Node >
  SerialPlatform<Node>::getNode() const { 
    return node_; 
  }

} // namespace Xpetra

#endif // XPETRA_SERIALPLATFORM_HPP
