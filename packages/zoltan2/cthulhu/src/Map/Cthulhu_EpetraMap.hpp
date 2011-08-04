/*
Important note: 

While Tpetra object take in argument a Tpetra::Map, Epetra objects take in argument a BlockMap.
So this adapter wraps an Epetra_BlockMap. In the current Epetra implementation, an Epetra_Map is just a Epetra_BlockMap.

More details:
-------------
There is an important difference between Epetra and Tpetra. 
Tpetra uses Map everywhere but in VbrMatrix. In Epetra, a lot of
thinks use BlockMap instead of a Map (basically, all the constructors
of Vector, Graph, Import etc...). So I think that the Cthulhu::Map
should be a Epetra_BlockMap when we use Epetra in Cthulhu.

But I faced a problem with the CrsMatrix class. Most Epetra_CrsMatrix
constructors need an Epetra_Map (ie: not a BlockMap) as input
argument. It is not consistent with the rest of the code. For
instance, that prevent you to create a CrsMatrix from the Map of any
DistObject because by inheritence, a Map is a BlockMap but the inverse
statment is false. CrsGraph constructors use BlockMap and thus can
either take an Epetra_Map or an Epetra_BlockMap. It seems very strange
to me that the interface of CrsMatrix constructors are not the same as
the interface of CrsGraph.

I took a look on the implementation of Epetra_CrsMatrix. Here are some comments:
- Epetra_CrsMatrix don't use specifically the fact that it is a Map
  and not a BlockMap. The map is just given as input argument to the
  graph constructor.
- You can create a CrsMatrix from a CrsGraph using another
  constructor. In this case, you don't need a Map. Note that CrsGraph
  don't know anything about Map but stores a BlockMap instead. So a
  CrsMatrix don't really use a Map.
- Epetra_CrsMatrix::ColMap() return a Map just by casting the BlockMap
  of the underlying graph object:
    const Epetra_Map& ColMap() const {return((Epetra_Map &) Graph_.ColMap());}
- I also checked the implementation of Epetra_Map. Basically, there is
  no difference at all between a Map and a BlockMap. So casting a
  BlockMap to a Map should be safe.
Finally, I think that all of that is just a design problem of Epetra.

To solve this problem, I cast the Epetra_BlockMap to a Epetra_Map like
it is done in Epetra_CrsMatrix::ColMap(). I need to do that only in
two spot: the constructors of CrsMatrix and
CrsMatrix::fillComplete(map,map). I think it is OK but maybe I need
to do more tests.
*/

#ifndef CTHULHU_EPETRAMAP_HPP
#define CTHULHU_EPETRAMAP_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Debug.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_Comm.hpp"
#include "Cthulhu_EpetraExceptions.hpp"

#include <Epetra_BlockMap.h>
#include <Epetra_Map.h>

namespace Tpetra { //TODO
typedef size_t global_size_t;
}

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

    // Implementation note for constructors: the Epetra_Comm is cloned in the constructor of Epetra_BlockMap. We don't need to keep a reference on it.

    /** \brief EpetraMap constructor with Cthulhu-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    EpetraMap(global_size_t numGlobalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
              LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) 
    {
      CTHULHU_DEBUG_ME;       

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
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode())
    {
      CTHULHU_DEBUG_ME; 

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
	      const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode())
    { 
      CTHULHU_DEBUG_ME; 

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(numGlobalElements, elementList.size(), elementList.getRawPtr(), 1, indexBase, *Teuchos2Epetra_Comm(comm))))));
    }
    
    /** \brief EpetraMap constructor to wrap a Epetra_BlockMap object.
     */
    EpetraMap(const Teuchos::RCP<const Epetra_BlockMap > &map) : map_(map) { CTHULHU_DEBUG_ME;}

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
    Cthulhu::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const int> & GIDList, 
                                            const Teuchos::ArrayView<      int> & nodeIDList, 
                                            const Teuchos::ArrayView<      int> & LIDList) const { CTHULHU_DEBUG_ME; 

      map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), LIDList.getRawPtr()); 
      
      return Cthulhu::AllIDsPresent; // JG TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
    };

    //! Returns the node IDs for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    Cthulhu::LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const int> & GIDList, 
                                            const Teuchos::ArrayView<      int> & nodeIDList) const { CTHULHU_DEBUG_ME; 

      // JG Note: It's not on the documentation of Epetra_BlockMap but it is in fact safe to call
      // EpetraMap RemoteIDList with LIDList == 0.
      // (because RemoteIDList only call directory->GetDirectoryEntries and this method accept LocalEntries=0)

      map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), 0); 

      return Cthulhu::AllIDsPresent; // JG TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
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
      return map_->MyLID(localIndex); 
    };

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(int globalIndex) const { CTHULHU_DEBUG_ME; return map_->MyGID(globalIndex); };

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
          return map_->PointSameAs(epetraMap.getEpetra_BlockMap()); 
	}
      catch (const std::bad_cast& e)
	{
          // We consider that an EpetraMap is never compatible with a map stored in another format (ie: a TpetraMap).
          // TODO: or throw exception ?
          return false;
        }
    }

    //! Returns true if \c map is identical to this Map.
    /** Note that an EpetraMap is never the 'same as' a TpetraMap. **/
    bool isSameAs (const Map<int,int,Kokkos::DefaultNode::DefaultNodeType> &map) const { CTHULHU_DEBUG_ME; 
      try
	{
          const EpetraMap & epetraMap = dynamic_cast<const EpetraMap &>(map);
          return map_->SameAs(epetraMap.getEpetra_BlockMap()); 
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
      const Teuchos::RCP<const Teuchos::Comm<int> > r = Epetra2Teuchos_Comm(rcpComm);
      return r;
    };

    //! Get the Node object for this Map
    const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> getNode() const { CTHULHU_DEBUG_ME;  //removed &
      return Kokkos::DefaultNode::getDefaultNode();
    };

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

    /** \brief Get the underlying Epetra map.. */
    const Epetra_BlockMap& getEpetra_BlockMap() const { CTHULHU_DEBUG_ME; return *map_; }
    
    /** \brief Get the underlying Epetra map. This method cast the underlying EpetraBlockMap to an EpetraMap.
        Try to use getEpetra_BlockMap() instead of getEpetra_Map() as much as possible. See the note in the top of the file Cthulhu_EpetraMap.hpp for more information. 
    */
    const Epetra_Map& getEpetra_Map() const { CTHULHU_DEBUG_ME; return (Epetra_Map &)*map_; } //TODO: write a note about that. It's the same in Epetra_CrsMatrix.h to get the map.

    inline UnderlyingLib lib() const { return Cthulhu::UseEpetra; };

  private:

    RCP<const Epetra_BlockMap> map_;

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

    inline Teuchos::RCP< const EpetraMap >
    createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return createLocalMapWithNode(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
    }

    /** \brief Non-member function to create a locally replicated Map with a specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    inline Teuchos::RCP< const EpetraMap >
    createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) { CTHULHU_DEBUG_ME;
      Teuchos::RCP< EpetraMap > map;
      map = Teuchos::rcp( new EpetraMap((Cthulhu::global_size_t)numElements, // num elements, global and local
                                        0,                                   // index base is zero
                                        comm, LocallyReplicated, node));
      return map.getConst();
    }

    /** \brief Non-member function to create a uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    inline Teuchos::RCP< const EpetraMap >
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
    inline Teuchos::RCP< const EpetraMap >
    createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return createUniformContigMapWithNode(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.

    This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    inline Teuchos::RCP< const EpetraMap >
    createContigMap(global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) { CTHULHU_DEBUG_ME;
      return createContigMapWithNode(numElements, localNumElements, comm, Kokkos::DefaultNode::getDefaultNode() );
    }

    /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    inline Teuchos::RCP< const EpetraMap >
    createContigMapWithNode(global_size_t numElements, size_t localNumElements, 
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node) { CTHULHU_DEBUG_ME;
      Teuchos::RCP< EpetraMap > map;
      map = Teuchos::rcp( new EpetraMap(numElements,localNumElements,
                                        0,  // index base is zero
                                        comm, node) );
      return map.getConst();
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    /** \brief Non-member function to create a contiguous Map with user-defined weights and a user-specified node.

    The Map is configured to use zero-based indexing.

    \relates EpetraMap
    */
    Teuchos::RCP< const EpetraMap >
    createWeightedContigMapWithNode(int myWeight, global_size_t numElements, 
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Kokkos::DefaultNode::DefaultNodeType > &node);
#endif

  } // useEpetra namespace

} // Cthulhu namespace

/** \brief  Returns true if \c map is identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
    \relates Cthulhu::EpetraMap */
inline bool operator== (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) { CTHULHU_DEBUG_ME;
  return map1.isSameAs(map2);
}

/** \brief Returns true if \c map is not identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
    \relates Cthulhu::EpetraMap */
inline bool operator!= (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) { CTHULHU_DEBUG_ME;
  return !map1.isSameAs(map2);
}

#endif // CTHULHU_EPETRAMAP_HPP
