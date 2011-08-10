#ifndef CTHULHU_EPETRAMAP_HPP
#define CTHULHU_EPETRAMAP_HPP

#include "Cthulhu_EpetraConfigDefs.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Kokkos_DefaultNode.hpp>

#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>

#include "Cthulhu_Map.hpp"

#include "Cthulhu_EpetraExceptions.hpp"
#include "Cthulhu_LookupStatus.hpp"
#include "Cthulhu_Comm.hpp"

namespace Tpetra { //TODO to be removed
  typedef size_t global_size_t;
}

namespace Cthulhu {

  const Epetra_Map & toEpetra(const Cthulhu::Map<int, int> &map);

  class EpetraMap
    : public Cthulhu::Map<int,int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

  public:

    //! @name Constructor/Destructor Methods
    //@{

    // Implementation note for constructors: the Epetra_Comm is cloned in the constructor of Epetra_BlockMap. We don't need to keep a reference on it.

    /** \brief EpetraMap constructor with Cthulhu-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as
     *   possible.
     */
    EpetraMap(global_size_t numGlobalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode())
    {
            

      // This test come from Tpetra (Epetra doesn't check if numGlobalElements,indexBase are equivalent across images).
      // In addition, for the test TEST_THROW(M map((myImageID == 0 ? GSTI : 0),0,comm), std::invalid_argument), only one node throw an exception and there is a dead lock.
      std::string errPrefix;
      errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,indexBase,comm,lOrG): ";
     
      if (lg == GloballyDistributed) {
        const int myImageID = comm->getRank();
       
        // check that numGlobalElements,indexBase is equivalent across images
        global_size_t rootNGE = numGlobalElements;
        int rootIB  = indexBase;
        Teuchos::broadcast<int,global_size_t>(*comm,0,&rootNGE);
        Teuchos::broadcast<int,int>(*comm,0,&rootIB);
        int localChecks[2], globalChecks[2];
        localChecks[0] = -1;   // fail or pass
        localChecks[1] = 0;    // fail reason
        if (numGlobalElements != rootNGE) {
          localChecks[0] = myImageID;
          localChecks[1] = 1;
        }
        else if (indexBase != rootIB) {
          localChecks[0] = myImageID;
          localChecks[1] = 2;
        }
        // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass, as well as the reason
        // these will be -1 and 0 if all procs passed
        Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
        if (globalChecks[0] != -1) {
          if (globalChecks[1] == 1) {
            TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal must be the same on all nodes (examine node " << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 2) {
            TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "indexBase must be the same on all nodes (examine node " << globalChecks[0] << ").");
          }
          else {
            // logic error on our part
            TEST_FOR_EXCEPTION(true,std::logic_error,
                               errPrefix << "logic error. Please contact the Tpetra team.");
          }
        }
      }
     
      // Note: validity of numGlobalElements checked by Epetra.

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(numGlobalElements, 1, indexBase, *Teuchos2Epetra_Comm(comm))))));
    }

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
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode())
    {
      

      // This test come from Tpetra
      using Teuchos::outArg;

      const size_t  L0 = Teuchos::OrdinalTraits<size_t>::zero();
      const size_t  L1 = Teuchos::OrdinalTraits<size_t>::one();
      const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
      const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
      const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid();

      std::string errPrefix;
      errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,numLocal,indexBase,platform): ";
     
      // get a internodal communicator from the Platform
      const int myImageID = comm->getRank();
     
      global_size_t global_sum;
      { // begin scoping block
        // for communicating failures
        int localChecks[2], globalChecks[2];
        /* compute the global size
           we are computing the number of global elements because exactly ONE of the following is true:
           - the user didn't specify it, and we need it
           - the user did specify it, but we need to
           + validate it against the sum of the local sizes, and
           + ensure that it is the same on all nodes
        */
        Teuchos::reduceAll<int,global_size_t>(*comm,Teuchos::REDUCE_SUM,
                                              Teuchos::as<global_size_t>(numLocalElements),outArg(global_sum));
        /* there are three errors we should be detecting:
           - numGlobalElements != invalid() and it is incorrect/invalid
           - numLocalElements invalid (<0)
        */
        localChecks[0] = -1;
        localChecks[1] = 0;
        if (numLocalElements < L1 && numLocalElements != L0) {
          // invalid
          localChecks[0] = myImageID;
          localChecks[1] = 1;
        }
        else if (numGlobalElements < GST1 && numGlobalElements != GST0 && numGlobalElements != GSTI) {
          // invalid
          localChecks[0] = myImageID;
          localChecks[1] = 2;
        }
        else if (numGlobalElements != GSTI && numGlobalElements != global_sum) {
          // incorrect
          localChecks[0] = myImageID;
          localChecks[1] = 3;
        }
        // now check that indexBase is equivalent across images
        int rootIB = indexBase;
        Teuchos::broadcast<int,int>(*comm,0,&rootIB);   // broadcast one ordinal from node 0
        if (indexBase != rootIB) {
          localChecks[0] = myImageID;
          localChecks[1] = 4;
        }
        // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
        // this will be -1 if all procs passed
        Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
        if (globalChecks[0] != -1) {
          if (globalChecks[1] == 1) {
            TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numLocal is not valid on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 2) {
            TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal is not valid on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 3) {
            TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal doesn't match sum of numLocal (== "
                               << global_sum << ") on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 4) {
            TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "indexBase is not the same on all nodes (examine node "
                               << globalChecks[0] << ").");
          }
          else {
            // logic error on my part
            TEST_FOR_EXCEPTION(true,std::logic_error,
                               errPrefix << "logic error. Please contact the Tpetra team.");
          }
        }
       
      }

      // set numGlobalElements
      if (numGlobalElements == GSTI) {
        numGlobalElements = global_sum;
      }

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(numGlobalElements, numLocalElements, 1, indexBase, *Teuchos2Epetra_Comm(comm))))));
    }
       
    /** \brief EpetraMap constructor with user-defined non-contiguous (arbitrary) distribution.
     * 
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
// TODO: UnitTest FAILED
    EpetraMap(global_size_t numGlobalElements, const Teuchos::ArrayView<const int> &elementList, int indexBase,
	      const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode())
    {
      

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(numGlobalElements, elementList.size(), elementList.getRawPtr(), 1, indexBase, *Teuchos2Epetra_Comm(comm))))));
    }
   
    /** \brief EpetraMap constructor to wrap a Epetra_BlockMap object.
     */
    EpetraMap(const Teuchos::RCP<const Epetra_BlockMap > &map) : map_(map) { }

    //! EpetraMap destructor.
    ~EpetraMap() { }

    //@}

    //! @name Map Attribute Methods
    //@{

    //! Returns the number of elements in this Map.
    global_size_t getGlobalNumElements() const { return map_->NumGlobalElements(); }

    //! Returns the number of elements belonging to the calling node.
    size_t getNodeNumElements() const { return map_->NumMyElements(); }

    //! Returns the index base for this Map.
    int getIndexBase() const { return map_->IndexBase(); }

    //! Returns minimum local index.
    int getMinLocalIndex() const { return map_->MinLID(); }

    //! Returns maximum local index.
    int getMaxLocalIndex() const { return map_->MaxLID(); }

    //! Returns minimum global index owned by this node.
    int getMinGlobalIndex() const { return map_->MinMyGID(); }

    //! Returns maximum global index owned by this node.
    int getMaxGlobalIndex() const { return map_->MaxMyGID(); }

    //! Return the minimum global index over all nodes.
    int getMinAllGlobalIndex() const { return map_->MinAllGID(); }

    //! Return the maximum global index over all nodes.
    int getMaxAllGlobalIndex() const { return map_->MaxAllGID(); }

    //! Return the local index for a given global index.
    int getLocalElement(GlobalOrdinal globalIndex) const { return map_->LID(globalIndex); }

    //! Return the global index for a given local index.
    int getGlobalElement(LocalOrdinal localIndex) const { return map_->GID(localIndex); }

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isNodeLocalElement(LocalOrdinal localIndex) const { return map_->MyLID(localIndex); }

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const { return map_->MyGID(globalIndex); }

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const { return map_->LinearMap(); }

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const { return map_->DistributedGlobal(); }

    //@}

    //! @name Boolean Tests
    //@{

    //! Returns true if map is compatible with this Map.
    bool isCompatible(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { return map_->PointSameAs(toEpetra(map)); }

    //! Returns true if map is identical to this Map.
    bool isSameAs(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { return map_->SameAs(toEpetra(map)); }

    //@}

    //! @name Map Attribute Methods
    //@{

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< int > &LIDList) const { return toCthulhu(map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), LIDList.getRawPtr())); }
    
    //! Returns the node IDs for a given list of global indices.
    // Note: It's not on the documentation of Epetra_BlockMap but it is in fact safe to call EpetraMap RemoteIDList(...,LIDList) with LIDList == 0.
    // (because RemoteIDList only call directory->GetDirectoryEntries and this method accept LocalEntries=0)
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const int > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const { return toCthulhu(map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), 0)); }
    
    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView< const int > getNodeElementList() const {
      int* nodeElementList = map_->MyGlobalElements(); // Pointer to *internal* array containing list of global IDs assigned to the calling processor.
      int  numMyElements   = map_->NumMyElements();    // Number of elements on the calling processor.
     
      // Note: this method return a const array, so it is safe to use directly the internal array.

      return ArrayView< const int >(nodeElementList, numMyElements);
    }

    //@}

    //! @name Implements Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const {
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

    //! @name Misc.
    //@{

    //! Get the Comm object for this Map
    const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const { //removed &
      RCP<const Epetra_Comm> rcpComm = rcpFromRef(map_->Comm());
      const Teuchos::RCP<const Teuchos::Comm<int> > r = Epetra2Teuchos_Comm(rcpComm);
      return r;
    }
    
    //! Get the Node object for this Map
    const Teuchos::RCP<Node> getNode() const { //removed &
      return Kokkos::DefaultNode::getDefaultNode();
    }
    
    //@}

    //! @name Cthulhu specific
    //@{

    //! Get the library used by this object (Tpetra or Epetra?)
    UnderlyingLib lib() const { return Cthulhu::UseEpetra; }

    //! Get the underlying Epetra map
    const Epetra_BlockMap& getEpetra_BlockMap() const { return *map_; }
    const Epetra_Map& getEpetra_Map() const { return (Epetra_Map &)*map_; } // Ugly, but the same is done in Epetra_CrsMatrix.h to get the map.

    //@}

  private:

    RCP<const Epetra_BlockMap> map_;
    //const RCP< const Epetra::Map< LocalOrdinal, GlobalOrdinal, Node > > map_;

  }; // EpetraMap class

} // Cthulhu namespace

// /** \brief  Returns true if \c map is identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
//     \relates Cthulhu::EpetraMap */
// bool operator== (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2);

// /** \brief Returns true if \c map is not identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
//     \relates Cthulhu::EpetraMap */
// bool operator!= (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2);

#endif // CTHULHU_EPETRAMAP_HPP
