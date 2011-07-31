#ifndef CTHULHU_TPETRAMAP_HPP
#define CTHULHU_TPETRAMAP_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

// enums and defines
#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Map.hpp"

//#ifdef HAVE_CTHULHU_TPETRA
//#include "Tpetra_Map_decl.hpp"
#include "Tpetra_Map.hpp"
//#endif

#include "Cthulhu_Debug.hpp"

/** \file Cthulhu_TpetraMap.hpp 

    The declarations for the class Cthulhu::TpetraMap and related non-member constructors.
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
  class TpetraMap : public Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    /** \brief TpetraMap constructor with Cthulhu-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    //TODO: replace Tpetra::LocalGlobal by Cthulhu::LocalGlobal
    TpetraMap(global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
              Cthulhu::LocalGlobal lg=Cthulhu::GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) //TODO: convert Cthulhu::LocalGlobal to Tpetra::LocalGlobal
      : map_(Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, indexBase, comm))) { CTHULHU_DEBUG_ME; } //, lg, node TODO TODO TODO
      
    /** \brief TpetraMap constructor with a user-defined contiguous distribution.
     *  The elements are distributed among the nodes so that the subsets of global elements
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
     TpetraMap(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, 
               const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) 
       : map_(Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, numLocalElements, indexBase, comm, node))) { CTHULHU_DEBUG_ME; }
        
    /** \brief TpetraMap constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
     TpetraMap(global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, 
               const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) 
       : map_(Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, elementList, indexBase, comm, node))) { CTHULHU_DEBUG_ME; }

    /** \brief TpetraMap constructor to wrap a Tpetra::Map object
     */
    TpetraMap(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &map) : map_(map) { CTHULHU_DEBUG_ME; }

    //! TpetraMap destructor. 
    ~TpetraMap() { CTHULHU_DEBUG_ME; }

    //@}

    //! @name TpetraMap Attribute Methods
    //@{ 

    //! Returns the number of elements in this Map.
    inline global_size_t getGlobalNumElements() const { CTHULHU_DEBUG_ME; return map_->getGlobalNumElements(); }

    //! Returns the number of elements belonging to the calling node.
    inline size_t getNodeNumElements() const { CTHULHU_DEBUG_ME; return map_->getNodeNumElements(); }

    //! Returns the index base for this Map.
    inline GlobalOrdinal getIndexBase() const { CTHULHU_DEBUG_ME; return map_->getIndexBase(); }

    //! Returns minimum local index
    inline LocalOrdinal getMinLocalIndex() const { CTHULHU_DEBUG_ME; return map_->getMinLocalIndex(); }

    //! Returns maximum local index
    inline LocalOrdinal getMaxLocalIndex() const { CTHULHU_DEBUG_ME; return map_->getMaxLocalIndex(); }

    //! Returns minimum global index owned by this node
    inline GlobalOrdinal getMinGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->getMinGlobalIndex(); }

    //! Returns maximum global index owned by this node
    inline GlobalOrdinal getMaxGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->getMaxGlobalIndex(); }

    //! Return the minimum global index over all nodes
    inline GlobalOrdinal getMinAllGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->getMinAllGlobalIndex(); }

    //! Return the maximum global index over all nodes
    inline GlobalOrdinal getMaxAllGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->getMaxAllGlobalIndex(); }

    //! \brief Return the local index for a given global index
    /** If the global index is not owned by this node, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid(). */
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const { CTHULHU_DEBUG_ME; return map_->getLocalElement(globalIndex); };

    //! Return the global index for a given local index
    /** If the local index is not valid for this node, returns Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(). */
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const { CTHULHU_DEBUG_ME; return map_->getGlobalElement(localIndex); };

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    //TODO: LookupStatus convert
    Cthulhu::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                             const Teuchos::ArrayView<                int> & nodeIDList, 
                                             const Teuchos::ArrayView<       LocalOrdinal> & LIDList) const { CTHULHU_DEBUG_ME; map_->getRemoteIndexList(GIDList, nodeIDList, LIDList); return AllIDsPresent; }; //TODO: convert Cthulhu::LookupStatus

    //! Returns the node IDs for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    Cthulhu::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                             const Teuchos::ArrayView<                int> & nodeIDList) const { CTHULHU_DEBUG_ME; map_->getRemoteIndexList(GIDList, nodeIDList); return AllIDsPresent; }; //TODO: convert Cthulhu::LookupStatus

    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList() const { CTHULHU_DEBUG_ME; return map_->getNodeElementList(); };

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isNodeLocalElement(LocalOrdinal localIndex) const { CTHULHU_DEBUG_ME; return map_->isNodeLocalElement(localIndex); };

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const { CTHULHU_DEBUG_ME; return map_->isNodeGlobalElement(globalIndex); };

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const { CTHULHU_DEBUG_ME; return map_->isContiguous(); };

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const { CTHULHU_DEBUG_ME; return map_->isDistributed(); };

    //@}

    //! @name Boolean Tests
    //@{ 

    //! Returns true if \c map is compatible with this Map.
    /** Note that a TpetraMap is never compatible with an EpetraMap. **/
    bool isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const { CTHULHU_DEBUG_ME; 
      try
	{
          const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> & tpetraMap = dynamic_cast<const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &>(map);
          return map_->isCompatible(*tpetraMap.getTpetra_Map()); 
	}
      catch (const std::bad_cast& e)
	{
          // We consider that an TpetraMap is never compatible with a map stored in another format (ie: a EpetraMap).
          return false;
	}
    };

    //! Returns true if \c map is identical to this Map.
    /** Note that an TpetraMap is never the 'same as' a EpetraMap. **/
    bool isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const { CTHULHU_DEBUG_ME; 
      try
	{
          const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> & tpetraMap = dynamic_cast<const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &>(map);
          return map_->isSameAs(*tpetraMap.getTpetra_Map()); 
	}
      catch (const std::bad_cast& e)
	{
          // We consider that an TpetraMap is never the 'same as' a map stored in an other formats (ie: EpetraMap).
          return false;
	}
    }
    //@}

    //@{ Misc. 

    //! Get the Comm object for this Map
    const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const { CTHULHU_DEBUG_ME; return map_->getComm(); }; // removed &

    //! Get the Node object for this Map
    const Teuchos::RCP<Node> getNode() const { CTHULHU_DEBUG_ME; return map_->getNode(); };

    //@}

    //@{ Implements Teuchos::Describable 

    //! \brief Return a simple one-line description of this object.
    std::string description() const { CTHULHU_DEBUG_ME; return map_->description(); };

    //! Print the object with some verbosity level to a \c FancyOStream object.
    void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const { CTHULHU_DEBUG_ME; map_->describe(out, verbLevel); };

    //@}

    inline UnderlyingLib lib() const { return Cthulhu::UseTpetra; };

    const RCP< const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > & getTpetra_Map() const { CTHULHU_DEBUG_ME; return map_; }

  private:

    const RCP< const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map_;

  }; // TpetraMap class

  namespace useTpetra {

    /** \brief Non-member function to create a locally replicated Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    template <class LocalOrdinal, class GlobalOrdinal>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
    createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
    }

    /** \brief Non-member function to create a locally replicated Map with a specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >
    createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;
      return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
    }

    /** \brief Non-member function to create a uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >
    createUniformContigMapWithNode(global_size_t numElements,
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;
      return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
    }

    /** \brief Non-member function to create a uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    template <class LocalOrdinal, class GlobalOrdinal>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
    createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createUniformContigMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    template <class LocalOrdinal, class GlobalOrdinal>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
    createContigMap(global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createContigMap<LocalOrdinal,GlobalOrdinal>(numElements, localNumElements, comm)));
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >
    createContigMapWithNode(global_size_t numElements, size_t localNumElements, 
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;
      return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, localNumElements, comm, node)));
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    /** \brief Non-member function to create a contiguous Map with user-defined weights and a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates TpetraMap
    */
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >
    createWeightedContigMapWithNode(int thisNodeWeight, global_size_t numElements, 
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) { CTHULHU_DEBUG_ME;
      return rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createWeightedContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(thisNodeWeight, numElements, comm, node)));
    }
#endif

  } // useTpetra namespace

} // Cthulhu namespace

/** \brief  Returns true if \c map is identical to this map. Implemented in Cthulhu::TpetraMap::isSameAs().
    \relates Cthulhu::TpetraMap */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator== (const Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map1, const Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map2) { CTHULHU_DEBUG_ME;
  return map1.isSameAs(map2);
}

/** \brief Returns true if \c map is not identical to this map. Implemented in Cthulhu::TpetraMap::isSameAs().
    \relates Cthulhu::TpetraMap */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator!= (const Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map1, const Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map2) { CTHULHU_DEBUG_ME;
  return !map1.isSameAs(map2);
}

#define CTHULHU_TPETRAMAP_SHORT
#endif // CTHULHU_TPETRAMAP_HPP

// NOTE: not copy constructor for Tpetra::Map ?

//TODO: why there is createBlaBla and createBlaBleWithNode in Tpetra ? Can we factorize that ?
