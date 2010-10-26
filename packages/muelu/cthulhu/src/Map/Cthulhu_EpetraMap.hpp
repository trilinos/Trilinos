#ifndef CTHULHU_EPETRAMAP_HPP
#define CTHULHU_EPETRAMAP_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Debug.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_EpetraComm.hpp"

#include <Epetra_Map.h>

/** \file Cthulhu_EpetraMap.hpp 

    The declarations for the class Cthulhu::EpetraMap and related non-member constructors.
 */

namespace Cthulhu {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward dec
  template <class LO, class GO, class N> class Directory;
#endif

  /** \brief A class for partitioning distributed objects.

   This class is templated on \c int and \c int. 
   The \c int type, if omitted, defaults to the \c int type.
  */
  class EpetraMap : public Cthulhu::Map<int,int> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    // Implementation note for constructors: the Epetra_Comm is cloned in the constructor of Epetra_Map. We don't need to keep a reference on it.

    /** \brief EpetraMap constructor with Cthulhu-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    EpetraMap(global_size_t numGlobalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
              LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) 
      :  map_(rcp(new Epetra_Map(numGlobalElements, indexBase, *Teuchos_Comm2Epetra_Comm(comm)))) { CTHULHU_DEBUG_ME; }
     
    /** \brief EpetraMap constructor with a user-defined contiguous distribution.
     *  The elements are distributed among the nodes so that the subsets of global elements
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    EpetraMap(global_size_t numGlobalElements, size_t numLocalElements, int indexBase, 
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode())
      : map_(rcp(new Epetra_Map(numGlobalElements, numLocalElements, indexBase, *Teuchos_Comm2Epetra_Comm(comm)))) { CTHULHU_DEBUG_ME;}
        
    /** \brief EpetraMap constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    EpetraMap(global_size_t numGlobalElements, const Teuchos::ArrayView<const int> &elementList, int indexBase, 
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode())
      : map_(rcp(new Epetra_Map(numGlobalElements, elementList.size(), elementList.getRawPtr(), indexBase, *Teuchos_Comm2Epetra_Comm(comm)))) { CTHULHU_DEBUG_ME;}

    /** \brief EpetraMap constructor to wrap a Epetra_Map object.
     */
    EpetraMap(const Teuchos::RCP<const Epetra_Map > &map) : map_(map) { CTHULHU_DEBUG_ME;}

    //! EpetraMap destructor. 
    ~EpetraMap() { CTHULHU_DEBUG_ME;}

    //@}

    //! @name EpetraMap Attribute Methods
    //@{ 

    //! Returns the number of elements in this Map.
    inline global_size_t getGlobalNumElements() const { CTHULHU_DEBUG_ME; return map_->NumGlobalElements(); }

    //! Returns the number of elements belonging to the calling node.
    inline size_t getNodeNumElements() const { CTHULHU_DEBUG_ME; return map_->NumMyElements(); }

    //! Returns the index base for this Map.
    inline int getIndexBase() const { CTHULHU_DEBUG_ME; return map_->IndexBase(); }

    //! Returns minimum local index
    inline int getMinLocalIndex() const { CTHULHU_DEBUG_ME; return map_->MinLID(); }

    //! Returns maximum local index
    inline int getMaxLocalIndex() const { CTHULHU_DEBUG_ME; return map_->MaxLID(); }

    //! Returns minimum global index owned by this node
    inline int getMinGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->MinMyGID(); }

    //! Returns maximum global index owned by this node
    inline int getMaxGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->MaxMyGID(); }

    //! Return the minimum global index over all nodes
    inline int getMinAllGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->MinAllGID(); }

    //! Return the maximum global index over all nodes
    inline int getMaxAllGlobalIndex() const { CTHULHU_DEBUG_ME; return map_->MaxAllGID(); }

    //! \brief Return the local index for a given global index
    /** If the global index is not owned by this node, returns Teuchos::OrdinalTraits<int>::invalid(). */
    int getLocalElement(int globalIndex) const { CTHULHU_DEBUG_ME; return map_->LID(globalIndex); };

    //! Return the global index for a given local index
    /** If the local index is not valid for this node, returns Teuchos::OrdinalTraits<int>::invalid(). */
    int getGlobalElement(int localIndex) const { CTHULHU_DEBUG_ME; return map_->GID(localIndex); };

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    Tpetra::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const int> & GIDList, 
                                            const Teuchos::ArrayView<      int> & nodeIDList, 
                                            const Teuchos::ArrayView<      int> & LIDList) const { CTHULHU_DEBUG_ME; 

      map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), LIDList.getRawPtr()); 
      
      return Tpetra::AllIDsPresent; // JG TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
    };

    //! Returns the node IDs for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    Tpetra::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const int> & GIDList, 
                                            const Teuchos::ArrayView<      int> & nodeIDList) const { CTHULHU_DEBUG_ME; 

      // JG Note: It's not on the documentation of Epetra_Map but it is in fact safe to call
      // EpetraMap RemoteIDList with LIDList == 0.
      // (because RemoteIDList only call directory->GetDirectoryEntries and this method accept LocalEntries=0)

      map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), 0); 

      return Tpetra::AllIDsPresent; // JG TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
    };

    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView<const int> getNodeElementList() const { CTHULHU_DEBUG_ME;
      int* nodeElementList = map_->MyGlobalElements(); // Pointer to *internal* array containing list of global IDs assigned to the calling processor. 
      int  numMyElements   = map_->NumMyElements();    // Number of elements on the calling processor.
      
      // JG Note: this method return a const array, so it is safe to use directly the internal array.

      return ArrayView<const int>(nodeElementList, numMyElements);
    };

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isNodeLocalElement(int localIndex) const { CTHULHU_DEBUG_ME; 
      return map_->LID(localIndex); 
    };

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(int globalIndex) const { CTHULHU_DEBUG_ME; return map_->GID(globalIndex); };

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const { CTHULHU_DEBUG_ME; return map_->LinearMap(); };

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const { CTHULHU_DEBUG_ME; return map_->DistributedGlobal(); };

    //@}

    //! @name Boolean Tests
    //@{ 

    //! Returns true if \c map is compatible with this Map.
    /** Note that an EpetraMap is never compatible with an TpetraMap. **/
    bool isCompatible (const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> &map) const { CTHULHU_DEBUG_ME; 
      try
	{
          const EpetraMap & epetraMap = dynamic_cast<const EpetraMap &>(map);
          return map_->PointSameAs(epetraMap.getEpetra_Map()); 
	}
      catch (const std::bad_cast& e)
	{
          // We consider that an EpetraMap is never compatible with a map stored in another format (ie: a TpetraMap).
          return false;
        }
    }

    //! Returns true if \c map is identical to this Map.
    /** Note that an EpetraMap is never the 'same as' a TpetraMap. **/
    bool isSameAs (const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> &map) const { CTHULHU_DEBUG_ME; 
      try
	{
          const EpetraMap & epetraMap = dynamic_cast<const EpetraMap &>(map);
          return map_->SameAs(epetraMap.getEpetra_Map()); 
	}
      catch (const std::bad_cast& e)
	{
          // We consider that an EpetraMap is never the 'same as' a map stored in an other formats (ie: TpetraMap).
          return false;
	}
    }

    //@}

    //@{ Misc. 

    //! Get the Comm object for this Map
    const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const { CTHULHU_DEBUG_ME;  //removed &
      RCP<const Epetra_Comm> rcpComm = rcpFromRef(map_->Comm());
      const Teuchos::RCP<const Teuchos::Comm<int> > r = Epetra_Comm2Teuchos_Comm(rcpComm);
      return r;
    };

    //! Get the Node object for this Map
#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> getNode() const { CTHULHU_DEBUG_ME; 
      TEST_FOR_EXCEPTIONS(1, Cthulhu::Exceptions::NotImplemented, 'Cthulhu::EpetraMap->getNode()');
      const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> & r = Teuchos::null;
      return r;
    };
#endif
    //@}

    //@{ Implements Teuchos::Describable 

    //! \brief Return a simple one-line description of this object.
    std::string description() const { 
      CTHULHU_DEBUG_ME; 
    
      // This implementation come from Tpetra_Map_def.hpp (without modification)
      std::ostringstream oss;
      oss << Teuchos::Describable::description();
      oss << "{getGlobalNumElements() = " << getGlobalNumElements()
          << ", getNodeNumElements() = " << getNodeNumElements()
          << ", isContiguous() = " << isContiguous()
          << ", isDistributed() = " << isDistributed()
          << "}";
      return oss.str();
    };

    //! Print the object with some verbosity level to a \c FancyOStream object.
    void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const { 
      CTHULHU_DEBUG_ME; 

      const Teuchos::RCP<const Teuchos::Comm<int> > comm_ = getComm();

      // This implementation come from Tpetra_Map_def.hpp (without modification)
      using std::endl;
      using std::setw;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;
      
      const size_t nME = getNodeNumElements();
      Teuchos::ArrayView<const int> myEntries = getNodeElementList();
      int myImageID = comm_->getRank();
      int numImages = comm_->getSize();
      
      Teuchos::EVerbosityLevel vl = verbLevel;
      if (vl == VERB_DEFAULT) vl = VERB_LOW;
      
      size_t width = 1;
      for (size_t dec=10; dec<getGlobalNumElements(); dec *= 10) {
        ++width;
      }
      width = std::max<size_t>(width,12) + 2;
      
      Teuchos::OSTab tab(out);
      
      if (vl == VERB_NONE) {
        // do nothing
      }
      else if (vl == VERB_LOW) {
        out << this->description() << endl;
      }
      else {  // MEDIUM, HIGH or EXTREME
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            if (myImageID == 0) { // this is the root node (only output this info once)
              out << endl 
                  << "Number of Global Entries = " << getGlobalNumElements()  << endl
                  << "Maximum of all GIDs      = " << getMaxAllGlobalIndex() << endl
                  << "Minimum of all GIDs      = " << getMinAllGlobalIndex() << endl
                  << "Index Base               = " << getIndexBase()         << endl;
            }
            out << endl;
            if (vl == VERB_HIGH || vl == VERB_EXTREME) {
              out << "Number of Local Elements   = " << nME           << endl
                  << "Maximum of my GIDs         = " << getMaxGlobalIndex() << endl
                  << "Minimum of my GIDs         = " << getMinGlobalIndex() << endl;
              out << endl;
            }
            if (vl == VERB_EXTREME) {
              out << std::setw(width) << "Node ID"
                  << std::setw(width) << "Local Index"
                  << std::setw(width) << "Global Index"
                  << endl;
              for (size_t i=0; i < nME; i++) {
                out << std::setw(width) << myImageID 
                    << std::setw(width) << i
                    << std::setw(width) << myEntries[i]
                    << endl;
              }
              out << std::flush;
            }
          }
          // Do a few global ops to give I/O a chance to complete
          comm_->barrier();
          comm_->barrier();
          comm_->barrier();
        }
      }
    }

    //@}

    const Epetra_Map& getEpetra_Map() const { CTHULHU_DEBUG_ME; return *map_; }

  private:

    RCP<const Epetra_Map> map_;

  }; // EpetraMap class

  namespace useEpetra {
    /** \brief Non-member function to create a locally replicated Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    
    // JG TODO: mv decl/defs
    Teuchos::RCP< const EpetraMap > createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node);
    Teuchos::RCP< const EpetraMap > createContigMapWithNode(global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node);

    Teuchos::RCP< const EpetraMap >
    createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return createLocalMapWithNode(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
    }

    /** \brief Non-member function to create a locally replicated Map with a specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) { CTHULHU_DEBUG_ME;
      Teuchos::RCP< EpetraMap > map;
      map = Teuchos::rcp( new EpetraMap((Tpetra::global_size_t)numElements, // num elements, global and local
                                        0,                                  // index base is zero
                                        comm, LocallyReplicated, node));
      return map.getConst();
    }

    /** \brief Non-member function to create a uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createUniformContigMapWithNode(global_size_t numElements,
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) { CTHULHU_DEBUG_ME;
      Teuchos::RCP< EpetraMap > map;
      map = Teuchos::rcp( new EpetraMap(numElements,        // num elements, global and local
                                        0,                  //index base is zero
                                        comm, GloballyDistributed, node));
      return map.getConst();
    }

    /** \brief Non-member function to create a uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return createUniformContigMapWithNode(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createContigMap(global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return createContigMapWithNode(numElements, localNumElements, comm, Kokkos::DefaultNode::getDefaultNode() );
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createContigMapWithNode(global_size_t numElements, size_t localNumElements, 
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) { CTHULHU_DEBUG_ME;
      Teuchos::RCP< EpetraMap > map;
      map = Teuchos::rcp( new EpetraMap(numElements,localNumElements,
                                        0,  // index base is zero
                                        comm, node) );
      return map.getConst();
    }

    /** \brief Non-member function to create a contiguous Map with user-defined weights and a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createWeightedContigMapWithNode(int myWeight, global_size_t numElements, 
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) { CTHULHU_DEBUG_ME;
      Teuchos::RCP< EpetraMap > map;
      int sumOfWeights, elemsLeft, localNumElements;
      const int numImages = comm->getSize(), 
        myImageID = comm->getRank();
      Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,myWeight,Teuchos::outArg(sumOfWeights));
      const double myShare = ((double)myWeight) / ((double)sumOfWeights);
      localNumElements = (int)std::floor( myShare * ((double)numElements) );
      // std::cout << "numElements: " << numElements << "  myWeight: " << myWeight << "  sumOfWeights: " << sumOfWeights << "  myShare: " << myShare << std::endl;
      Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,localNumElements,Teuchos::outArg(elemsLeft));
      elemsLeft = numElements - elemsLeft;
      // std::cout << "(before) localNumElements: " << localNumElements << "  elemsLeft: " << elemsLeft << std::endl;
      // i think this is true. just test it for now.
      TEST_FOR_EXCEPT(elemsLeft < -numImages || numImages < elemsLeft);
      if (elemsLeft < 0) {
        // last elemsLeft nodes lose an element
        if (myImageID >= numImages-elemsLeft) --localNumElements;
      }
      else if (elemsLeft > 0) {
        // first elemsLeft nodes gain an element
        if (myImageID < elemsLeft) ++localNumElements;
      }
      // std::cout << "(after) localNumElements: " << localNumElements << std::endl;
      return createContigMapWithNode(numElements,localNumElements,comm,node);
    }

  } // useEpetra namespace

} // Cthulhu namespace

/** \brief  Returns true if \c map is identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
    \relates Cthulhu::EpetraMap */
bool operator== (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) { CTHULHU_DEBUG_ME;
  return map1.isSameAs(map2);
}

/** \brief Returns true if \c map is not identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
    \relates Cthulhu::EpetraMap */
bool operator!= (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) { CTHULHU_DEBUG_ME;
  return !map1.isSameAs(map2);
}

#endif // CTHULHU_EPETRAMAP_HPP
