#ifndef CTHULHU_MPIPLATFORM_HPP
#define CTHULHU_MPIPLATFORM_HPP

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Debug.hpp"

namespace Cthulhu {

	//! \brief A implementation of the Platform class for MPI-based platforms.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class Node=Kokkos::DefaultNode::DefaultNodeType>
	class MpiPlatform : public Teuchos::Describable {
    public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
      typedef Node NodeType;
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor
      explicit MpiPlatform(Teuchos::RCP<Node> node);

      //! Constructor
      MpiPlatform(Teuchos::RCP<Node> node, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm);

      //! Destructor
      ~MpiPlatform();

      //@}

      //! @name Class Creation and Accessor Methods
      //@{ 

      //! Comm Instance
      Teuchos::RCP< const Teuchos::Comm<int> > getComm() const;

      //! Get Get a node for parallel computation.
      Teuchos::RCP<Node> getNode() const;

      //@}

    protected: 
      //! Node object instantiated for the platform.
      Teuchos::RCP<Node> node_;

    private:
      Teuchos::RCP<Teuchos::MpiComm<int> > comm_;
      MpiPlatform(const MpiPlatform<Node> &platform);
  };

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> node, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm)
  : node_(node) { CTHULHU_DEBUG_ME;
    comm_ = Teuchos::createMpiComm<int>(rawMpiComm);
  }

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> node)
  : node_(node) { CTHULHU_DEBUG_ME;
    comm_ = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  } 

  template <class Node>
  MpiPlatform<Node>::~MpiPlatform() { CTHULHU_DEBUG_ME; }

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(const MpiPlatform<Node> &platform) { CTHULHU_DEBUG_ME;
    comm_ = platform.comm_;
  }

  template <class Node>
  Teuchos::RCP< const Teuchos::Comm<int> > 
  MpiPlatform<Node>::getComm() const { CTHULHU_DEBUG_ME;
    return comm_;
  }

  template <class Node>
  Teuchos::RCP<Node> MpiPlatform<Node>::getNode() const 
  { CTHULHU_DEBUG_ME; return node_; }

} // namespace Cthulhu

#endif // CTHULHU_MPIPLATFORM_HPP

