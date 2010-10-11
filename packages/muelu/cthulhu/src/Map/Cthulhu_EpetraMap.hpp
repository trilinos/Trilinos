#ifndef CTHULHU_EPETRAMAP_HPP
#define CTHULHU_EPETRAMAP_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

// enums and defines
#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Map.hpp"

#include "Tpetra_Map.hpp" //tmp
#include "Epetra_Map.h"

#include <Cthulhu_EpetraComm.hpp>

#include <Epetra_SerialComm.h> //tmp
/** \file Cthulhu_EpetraMap.hpp 

    The declarations for the class Cthulhu::EpetraMap and related non-member constructors.
 */

namespace Cthulhu {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward dec
  template <class LO, class GO, class N> class Directory;
#endif

  /** \brief A class for partitioning distributed objects.

   This class is templated on \c int and \c GlobalOrdinal. 
   The \c GlobalOrdinal type, if omitted, defaults to the \c int type.
  */
  class EpetraMap : public Cthulhu::Map<int,int> {

  public:

    // JG Note: we don't need this bunch constructor and they complicate the code (conversion of Comm, Vector etc.). Let see if I keep them.

    //! @name Constructor/Destructor Methods
    //@{ 

    /** \brief EpetraMap constructor with Cthulhu-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    EpetraMap(global_size_t numGlobalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
              LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode())
      : map_(rcp(new Epetra_Map(numGlobalElements, indexBase, *Teuchos_Comm2Epetra_Comm(comm)))) {}
      // JG Note: epetraComm is cloned in the constructor of Epetra_Map. We don't need to keep a reference on epetraComm.

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
      : map_(rcp(new Epetra_Map(numGlobalElements, numLocalElements, indexBase, *Teuchos_Comm2Epetra_Comm(comm)))) {}
        
    /** \brief EpetraMap constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    EpetraMap(global_size_t numGlobalElements, const Teuchos::ArrayView<const int> &elementList, int indexBase, 
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode())
      : map_(rcp(new Epetra_Map(numGlobalElements, elementList.size(), elementList.getRawPtr(), indexBase, *Teuchos_Comm2Epetra_Comm(comm)))) {}

    /** \brief EpetraMap constructor to wrap a Epetra_Map object
     */

    EpetraMap(const Teuchos::RCP<const Epetra_Map > &map) : map_(map) {}

    //! EpetraMap destructor. 
    ~EpetraMap() {}

    //@}

    //! @name EpetraMap Attribute Methods
    //@{ 

    //! Returns the number of elements in this Map.
    inline global_size_t getGlobalNumElements() const { return map_->NumGlobalElements(); }

    //! Returns the number of elements belonging to the calling node.
    inline size_t getNodeNumElements() const { return map_->NumMyElements(); }

    //! Returns the index base for this Map.
    inline int getIndexBase() const { return map_->IndexBase(); }

    //! Returns minimum local index
    inline int getMinLocalIndex() const { return map_->MinLID(); }

    //! Returns maximum local index
    inline int getMaxLocalIndex() const { return map_->MaxLID(); }

    //! Returns minimum global index owned by this node
    inline int getMinGlobalIndex() const { return map_->MinMyGID(); }

    //! Returns maximum global index owned by this node
    inline int getMaxGlobalIndex() const { return map_->MaxMyGID(); }

    //! Return the minimum global index over all nodes
    inline int getMinAllGlobalIndex() const { return map_->MinAllGID(); }

    //! Return the maximum global index over all nodes
    inline int getMaxAllGlobalIndex() const { return map_->MaxAllGID(); }

    //! \brief Return the local index for a given global index
    /** If the global index is not owned by this node, returns Teuchos::OrdinalTraits<int>::invalid(). */
    int getLocalElement(int globalIndex) const { return map_->LID(globalIndex); };

    //! Return the global index for a given local index
    /** If the local index is not valid for this node, returns Teuchos::OrdinalTraits<int>::invalid(). */
    int getGlobalElement(int localIndex) const { return map_->GID(localIndex); };

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    Tpetra::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const int> & GIDList, 
                                            const Teuchos::ArrayView<      int> & nodeIDList, 
                                            const Teuchos::ArrayView<      int> & LIDList) const { 

      map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), LIDList.getRawPtr()); 
      
      return Tpetra::AllIDsPresent; // JG TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
    };

    //! Returns the node IDs for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    Tpetra::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const int> & GIDList, 
                                            const Teuchos::ArrayView<      int> & nodeIDList) const { 

      // JG Note: It's not on the documentation of Epetra_Map but it is in fact safe to call
      // EpetraMap RemoteIDList with LIDList == 0.
      // (because RemoteIDList only call directory->GetDirectoryEntries and this method accept LocalEntries=0)

      map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), 0); 

      return Tpetra::AllIDsPresent; // JG TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
    };

    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView<const int> getNodeElementList() const {
      int* nodeElementList = map_->MyGlobalElements(); // Pointer to *internal* array containing list of global IDs assigned to the calling processor. 
      int  numMyElements   = map_->NumMyElements();    // Number of elements on the calling processor.
      
      // JG Note: this method return a const array, so it is safe to use directly the internal array.

      return ArrayView<const int>(nodeElementList, numMyElements);
    };

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isNodeLocalElement(int localIndex) const { 
      return map_->LID(localIndex); 
    };

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(int globalIndex) const { return map_->GID(globalIndex); };

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const { return map_->LinearMap(); };

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const { return map_->DistributedGlobal(); };

    //@}

    //! @name Boolean Tests
    //@{ 

    //! Returns true if \c map is compatible with this Map.
    bool isCompatible (const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> &map) const { 
      try
	{
          const EpetraMap & epetraMap = dynamic_cast<const EpetraMap &>(map);
          return map_->PointSameAs(epetraMap.getEpetra_Map()); 
	}
      catch (const std::bad_cast& e)
	{
          // Let say that a EpetraMap is not compatible with other format of Map.
          return false;
        }
    }

    //! Returns true if \c map is identical to this Map.
    bool isSameAs (const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> &map) const { 
      try
	{
          const EpetraMap & epetraMap = dynamic_cast<const EpetraMap &>(map);
          return map_->SameAs(epetraMap.getEpetra_Map()); 
	}
      catch (const std::bad_cast& e)
	{
          // Let say that a EpetraMap is not compatible with other format of Map.
          return false;
	}
    }

    //@}

    //@{ Misc. 

    //! Get the Comm object for this Map
    const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const { 
      RCP<const Epetra_Comm> rcpComm = rcpFromRef(map_->Comm());
      const Teuchos::RCP<const Teuchos::Comm<int> > & r = Epetra_Comm2Teuchos_Comm(rcpComm);
      return r;
    };

    //! Get the Node object for this Map
    const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> & getNode() const { 
      const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> & r = Teuchos::null;
      return r;
    }; // JG TODO

    //@}

    //@{ Implements Teuchos::Describable 

    //! \brief Return a simple one-line description of this object.
    std::string description() const { return NULL; };

    //! Print the object with some verbosity level to a \c FancyOStream object.
    // TODO   void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const { map_->describe(); };

    //@}

    const Epetra_Map& getEpetra_Map() const { return *map_; } //JG TODO: wrap output in a RCP ?

  private:

    RCP< const Epetra_Map > map_; // TODO JG const ?

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
    createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
      return createLocalMapWithNode(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
    }

    /** \brief Non-member function to create a locally replicated Map with a specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) {
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
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) {
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
    createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
      return createUniformContigMapWithNode(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createContigMap(global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
      return createContigMapWithNode(numElements, localNumElements, comm, Kokkos::DefaultNode::getDefaultNode() );
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createContigMapWithNode(global_size_t numElements, size_t localNumElements, 
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) {
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
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) {
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
bool operator== (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) {
  return map1.isSameAs(map2);
}

/** \brief Returns true if \c map is not identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
    \relates Cthulhu::EpetraMap */
bool operator!= (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) {
  return !map1.isSameAs(map2);
}

#endif // CTHULHU_EPETRAMAP_HPP
