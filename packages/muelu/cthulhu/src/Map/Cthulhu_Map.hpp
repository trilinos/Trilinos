#ifndef CTHULHU_MAP_HPP
#define CTHULHU_MAP_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

#include "Tpetra_ConfigDefs.hpp" //Tpetra::LookupStatus

// enums and defines
#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Debug.hpp"

/** \file Cthulhu_Map.hpp 

    The declarations for the class Cthulhu::Map and related non-member constructors.
 */

namespace Cthulhu {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward dec
  template <class LO, class GO, class N> class Directory;
#endif

  /** \brief A class for partitioning distributed objects.

   This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
   The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Map : public Teuchos::Describable {

  public:

//     //! @name Constructor/Destructor Methods
//     //@{ 

//     /** \brief Map constructor with Cthulhu-defined contiguous uniform distribution.
//      *   The elements are distributed among nodes so that the subsets of global elements
//      *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
//      *   possible.
//      */
//     Map(global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
//         LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

//     /** \brief Map constructor with a user-defined contiguous distribution.
//      *  The elements are distributed among the nodes so that the subsets of global elements
//      *  are non-overlapping and contiguous 
//      *  
//      *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
//      *  Otherwise, it must be equal to the sum of the local elements across all 
//      *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
//      *  If this verification fails, a std::invalid_argument exception will be thrown.
//      */
//     Map(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, 
//         const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());
        

//     /** \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
//      *  
//      *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
//      *  Otherwise, it must be equal to the sum of the local elements across all 
//      *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
//      *  If this verification fails, a std::invalid_argument exception will be thrown.
//      */
//     Map(global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, 
//         const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

//     //! Map destructor. 
    virtual ~Map() { CTHULHU_DEBUG_ME; };

//     //@}

    //! @name Map Attribute Methods
    //@{ 

    //! Returns the number of elements in this Map.
    virtual global_size_t getGlobalNumElements() const = 0;

    //! Returns the number of elements belonging to the calling node.
    virtual size_t getNodeNumElements() const = 0;

    //! Returns the index base for this Map.
    virtual GlobalOrdinal getIndexBase() const = 0;

    //! Returns minimum local index
    virtual LocalOrdinal getMinLocalIndex() const = 0;

    //! Returns maximum local index
    virtual LocalOrdinal getMaxLocalIndex() const = 0;

    //! Returns minimum global index owned by this node
    virtual GlobalOrdinal getMinGlobalIndex() const = 0;

    //! Returns maximum global index owned by this node
    virtual GlobalOrdinal getMaxGlobalIndex() const = 0;

    //! Return the minimum global index over all nodes
    virtual GlobalOrdinal getMinAllGlobalIndex() const = 0;

    //! Return the maximum global index over all nodes
    virtual GlobalOrdinal getMaxAllGlobalIndex() const = 0;

    //! \brief Return the local index for a given global index
    /** If the global index is not owned by this node, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid(). */
    virtual LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const = 0;

    //! Return the global index for a given local index
    /** If the local index is not valid for this node, returns Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(). */
    virtual GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const = 0;

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    //TODO: LookupStatus
    virtual Tpetra::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                                                    const Teuchos::ArrayView<                int> & nodeIDList, 
                                                                    const Teuchos::ArrayView<       LocalOrdinal> & LIDList) const = 0;

    //! Returns the node IDs for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    virtual Tpetra::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                                                    const Teuchos::ArrayView<                int> & nodeIDList) const = 0;

    //! Return a list of the global indices owned by this node.
    virtual Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList() const = 0;

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

    //! Returns true if \c map is compatible with this Map.
    virtual bool isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const = 0;

    //! Returns true if \c map is identical to this Map.
    virtual bool isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const = 0;

    //@}

    //@{ Misc. 

    //! Get the Comm object for this Map
    virtual const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const = 0;

    //! Get the Node object for this Map
    virtual const Teuchos::RCP<Node> & getNode() const = 0;

    //@}

    //@{ Implements Teuchos::Describable 

    //! \brief Return a simple one-line description of this object.
    virtual std::string description() const = 0;

    //! Print the object with some verbosity level to a \c FancyOStream object.
    virtual void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const = 0;

    //@}

  }; // Map class

} // Cthulhu namespace

// /** \brief  Returns true if \c map is identical to this map. Implemented in Cthulhu::Map::isSameAs().
//     \relates Cthulhu::Map */
// template <class LocalOrdinal, class GlobalOrdinal, class Node>
// bool operator== (const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> &map1, const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> &map2);

// /** \brief Returns true if \c map is not identical to this map. Implemented in Cthulhu::Map::isSameAs().
//     \relates Cthulhu::Map */
// template <class LocalOrdinal, class GlobalOrdinal, class Node>
// bool operator!= (const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> &map1, const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> &map2);

#endif // CTHULHU_MAP_DECL_HPP
