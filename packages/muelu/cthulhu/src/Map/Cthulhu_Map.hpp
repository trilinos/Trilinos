#ifndef CTHULHU_MAP_HPP
#define CTHULHU_MAP_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include "Cthulhu_ConfigDefs.hpp"

namespace Cthulhu {

  enum UnderlyingLib {
    UseEpetra,
    UseTpetra
  };

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Map
    : public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    // !Destructor.
    virtual ~Map() { }

    //@}

    //! @name Map Attribute Methods
    //@{

    //! Returns the number of elements in this Map.
    virtual global_size_t getGlobalNumElements() const = 0;

    //! Returns the number of elements belonging to the calling node.
    virtual size_t getNodeNumElements() const = 0;

    //! Returns the index base for this Map.
    virtual GlobalOrdinal getIndexBase() const = 0;

    //! Returns minimum local index.
    virtual LocalOrdinal getMinLocalIndex() const = 0;

    //! Returns maximum local index.
    virtual LocalOrdinal getMaxLocalIndex() const = 0;

    //! Returns minimum global index owned by this node.
    virtual GlobalOrdinal getMinGlobalIndex() const = 0;

    //! Returns maximum global index owned by this node.
    virtual GlobalOrdinal getMaxGlobalIndex() const = 0;

    //! Return the minimum global index over all nodes.
    virtual GlobalOrdinal getMinAllGlobalIndex() const = 0;

    //! Return the maximum global index over all nodes.
    virtual GlobalOrdinal getMaxAllGlobalIndex() const = 0;

    //! Return the local index for a given global index.
    virtual LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const = 0;

    //! Return the global index for a given local index.
    virtual GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const = 0;

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const = 0;

    //! Returns the node IDs for a given list of global indices.
    virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const = 0;

    //! Return a list of the global indices owned by this node.
    virtual Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const = 0;

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    virtual bool isNodeLocalElement(LocalOrdinal localIndex) const = 0;

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    virtual bool isNodeGlobalElement(GlobalOrdinal globalIndex) const = 0;

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    virtual bool isContiguous() const = 0;

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    virtual bool isDistributed() const = 0;

    //@}

    //! @name Boolean Tests
    //@{

    //! Returns true if map is compatible with this Map.
    virtual bool isCompatible(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const = 0;

    //! Returns true if map is identical to this Map.
    virtual bool isSameAs(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const = 0;

    //@}

    //! @name Implements Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    virtual std::string description() const = 0;

    //! Print the object with some verbosity level to a FancyOStream object.
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

    //@}

    //! @name Misc.
    //@{

    //! Get the Comm object for this Map.
    virtual const Teuchos::RCP< const Teuchos::Comm< int > > getComm() const = 0;

    //! Get the Node object for this Map.
    virtual const Teuchos::RCP< Node > getNode() const = 0;

    //@}

    //! @name Cthulhu specific
    //@{

    //! Get the library used by this object (Tpetra or Epetra?)
    virtual UnderlyingLib lib() const = 0;

    //@}

  }; // Map class

} // Cthulhu namespace

#define CTHULHU_MAP_SHORT
#endif // CTHULHU_MAP_HPP
