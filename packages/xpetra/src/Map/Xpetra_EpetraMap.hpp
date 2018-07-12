// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_EPETRAMAP_HPP
#define XPETRA_EPETRAMAP_HPP


#include "Xpetra_EpetraConfigDefs.hpp"

#include "Xpetra_Map.hpp"

#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>

#include "Xpetra_Utils.hpp"
#include "Xpetra_EpetraUtils.hpp"
#include "Xpetra_EpetraExceptions.hpp"

#include "Xpetra_ConfigDefs.hpp"

namespace Xpetra {

  // TODO: move that elsewhere
  template<class GlobalOrdinal, class Node>
  const Epetra_Map & toEpetra(const Map<int,GlobalOrdinal, Node> &);
  template<class GlobalOrdinal, class Node>
  const Epetra_Map & toEpetra(const RCP< const Map<int, GlobalOrdinal, Node> > &);
  //template<class GlobalOrdinal>
  //const RCP< const Map<int, GlobalOrdinal> > toXpetra(const RCP< const Epetra_Map > &);
  template<class GlobalOrdinal, class Node>
  const RCP< const Map<int, GlobalOrdinal, Node> > toXpetra(const Epetra_BlockMap &);

  // stub implementation for EpetraMapT
  template<class GlobalOrdinal, class Node>
  class EpetraMapT
    : public virtual Map<int, GlobalOrdinal, Node>
  {
    typedef int LocalOrdinal;

  public:
    typedef int local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

    //! @name Constructors and destructor
    //@{

    //! Constructor with Tpetra-defined contiguous uniform distribution.
    EpetraMapT(global_size_t numGlobalElements,
               GlobalOrdinal indexBase,
               const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
               LocalGlobal lg=GloballyDistributed,
               const Teuchos::RCP< Node > &node = Teuchos::rcp (new Node))
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
        "Xpetra::EpetraMap only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
    }

    //! Constructor with a user-defined contiguous distribution.
    EpetraMapT(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::rcp(new Node)) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
        "Xpetra::EpetraMap only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
    }

    //! Constructor with user-defined arbitrary (possibly noncontiguous) distribution.
    EpetraMapT(global_size_t numGlobalElements,
        const Teuchos::ArrayView< const GlobalOrdinal > &elementList,
        GlobalOrdinal indexBase,
        const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
        const Teuchos::RCP< Node > &node = Teuchos::rcp(new Node)) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
        "Xpetra::EpetraMap only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
    }

    //@}

    //! @name Attributes
    //@{

    //! The number of elements in this Map.
    global_size_t getGlobalNumElements() const { return 0; }

    //! The number of elements belonging to the calling process.
    size_t getNodeNumElements() const { return 0; }

    //! The index base for this Map.
    GlobalOrdinal getIndexBase() const { return 0; }

    //! The minimum local index.
    LocalOrdinal getMinLocalIndex() const { return 0; }

    //! The maximum local index on the calling process.
    LocalOrdinal getMaxLocalIndex() const { return 0; }

    //! The minimum global index owned by the calling process.
    GlobalOrdinal getMinGlobalIndex() const { return 0; }

    //! The maximum global index owned by the calling process.
    GlobalOrdinal getMaxGlobalIndex() const { return 0; }

    //! The minimum global index over all processes in the communicator.
    GlobalOrdinal getMinAllGlobalIndex() const { return 0; }

    //! The maximum global index over all processes in the communicator.
    GlobalOrdinal getMaxAllGlobalIndex() const { return 0; }

    //! The local index corresponding to the given global index.
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const { return 0; }

    //! Return the process ranks and corresponding local indices for the given global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const { return Xpetra::IDNotPresent; }

    //! Return the process ranks for the given global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const { return Xpetra::IDNotPresent; }

    //! Return a view of the global indices owned by this process.
    //Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const;

    Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const { return ArrayView< const GlobalOrdinal >(); }
    //@}

    //! @name Boolean tests
    //@{

    //! Whether the given local index is valid for this Map on this process.
    bool isNodeLocalElement(LocalOrdinal localIndex) const { return false; }

    //! Whether the given global index is valid for this Map on this process.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const { return false; }

    //! True if this Map is distributed contiguously, else false.
    bool isContiguous() const { return false; }

    //! Whether this Map is globally distributed or locally replicated.
    bool isDistributed() const { return false; }

    //! True if and only if map is compatible with this Map.
    bool isCompatible(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { return false; }

    //! True if and only if map is identical to this Map.
    bool isSameAs(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { return false; }

    //@}

    //! @name
    //@{

    //! Get this Map's Comm object.
    Teuchos::RCP< const Teuchos::Comm< int > > getComm() const { return Teuchos::null; }

    //! Get this Map's Node object.
    Teuchos::RCP< Node > getNode() const {
      XPETRA_MONITOR("EpetraMapT<GlobalOrdinal>::getNode");
      return Teuchos::rcp (new Node);
    }

    //@}

    //! @name
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const { return std::string(""); }

    //! Print this object with the given verbosity level to the given Teuchos::FancyOStream.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  }

    //@}

    //! @name
    //@{

    //! Return a new Map with processes with zero elements removed.
    RCP<const Map<int,GlobalOrdinal,Node> > removeEmptyProcesses() const { return Teuchos::null; }

    //! Replace this Map's communicator with a subset communicator.
    RCP<const Map<int,GlobalOrdinal,Node> > replaceCommWithSubset(const Teuchos::RCP< const Teuchos::Comm< int > > &newComm) const { return Teuchos::null; }

    //@}

    //! Return the global index for a given local index.  Note that this returns -1 if not found on this processor.  (This is different than Epetra's behavior!)
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const { return -1; }

    //! @name Xpetra specific
    //@{

    //! Destructor.
    virtual ~EpetraMapT() {}

    //! EpetraMapT constructor to wrap a Epetra_Map object
    EpetraMapT(const Teuchos::RCP<const Epetra_BlockMap> &map)
      : map_(map) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
        "Xpetra::EpetraMap only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
    }

    //! Get the library used by this object (Epetra or Epetra?)
    UnderlyingLib lib() const { return Xpetra::UseEpetra; }

    //! Get the underlying Epetra map
    //const RCP< const Epetra_Map > & getEpetra_Map() const { return map_; }
    const Epetra_BlockMap& getEpetra_BlockMap() const { return *map_; }
    const Epetra_Map& getEpetra_Map() const { return (Epetra_Map &)*map_; } // Ugly, but the same is done in Epetra_CrsMatrix.h to get the map.

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    using local_map_type = typename Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type;
    /// \brief Get the local Map for Kokkos kernels.
    local_map_type getLocalMap () const {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
        "Xpetra::EpetraMap only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
      TEUCHOS_UNREACHABLE_RETURN(local_map_type());
    }
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif
#endif

    //@}

  protected:

    RCP<const Epetra_BlockMap> map_;
  }; // EpetraMapT class

  // specialization on GO=int and EpetraNode
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
  template<>
  class EpetraMapT<int, EpetraNode>
    : public virtual Map<int, int, EpetraNode>
  {
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

    //! @name Constructors and destructor
    //@{

    //! Constructor with Tpetra-defined contiguous uniform distribution.
    EpetraMapT(global_size_t numGlobalElements,
               GlobalOrdinal indexBase,
               const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
               LocalGlobal lg=GloballyDistributed,
               const Teuchos::RCP< Node > &node = Teuchos::rcp(new Node))
    {
      // This test come from Tpetra (Epetra doesn't check if numGlobalElements,indexBase are equivalent across images).
      // In addition, for the test TEST_THROW(M map((myImageID == 0 ? GSTI : 0),0,comm), std::invalid_argument), only one node throw an exception and there is a dead lock.
      std::string errPrefix;
      errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,indexBase,comm,lOrG): ";

      if (lg == GloballyDistributed) {
        const int myImageID = comm->getRank();

        // check that numGlobalElements,indexBase is equivalent across images
        global_size_t rootNGE = numGlobalElements;
        GlobalOrdinal rootIB  = indexBase;
        Teuchos::broadcast<int,global_size_t>(*comm,0,&rootNGE);
        Teuchos::broadcast<int,GlobalOrdinal>(*comm,0,&rootIB);
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
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal must be the same on all nodes (examine node " << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "indexBase must be the same on all nodes (examine node " << globalChecks[0] << ").");
          }
          else {
            // logic error on our part
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               errPrefix << "logic error. Please contact the Tpetra team.");
          }
        }
      }

      // Note: validity of numGlobalElements checked by Epetra.

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(static_cast<GlobalOrdinal>(numGlobalElements), 1, indexBase, *toEpetra(comm))))));
    }

    //! Constructor with a user-defined contiguous distribution.
    EpetraMapT(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::rcp(new Node))
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
        GlobalOrdinal rootIB = indexBase;
        Teuchos::broadcast<int,GlobalOrdinal>(*comm,0,&rootIB);   // broadcast one ordinal from node 0
        if (indexBase != rootIB) {
          localChecks[0] = myImageID;
          localChecks[1] = 4;
        }
        // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
        // this will be -1 if all procs passed
        Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
        if (globalChecks[0] != -1) {
          if (globalChecks[1] == 1) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numLocal is not valid on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal is not valid on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 3) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal doesn't match sum of numLocal (== "
                               << global_sum << ") on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 4) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "indexBase is not the same on all nodes (examine node "
                               << globalChecks[0] << ").");
          }
          else {
            // logic error on my part
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               errPrefix << "logic error. Please contact the Tpetra team.");
          }
        }

      }

      // set numGlobalElements
      if (numGlobalElements == GSTI) {
        numGlobalElements = global_sum;}

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(static_cast<GlobalOrdinal>(numGlobalElements), static_cast<int>(numLocalElements), 1, indexBase, *toEpetra(comm))))));
    }

    //! Constructor with user-defined arbitrary (possibly noncontiguous) distribution.
    EpetraMapT(global_size_t numGlobalElements,
               const Teuchos::ArrayView< const GlobalOrdinal > &elementList,
               GlobalOrdinal indexBase,
               const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
               const Teuchos::RCP< Node > &node = Teuchos::rcp(new Node))
    {
      if (numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid()) {
        IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(-1, static_cast<int>(elementList.size()), elementList.getRawPtr(), 1, indexBase, *toEpetra(comm))))));
      } else {
        IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(static_cast<int>(numGlobalElements), static_cast<int>(elementList.size()), elementList.getRawPtr(), 1, indexBase, *toEpetra(comm))))));
      }
    }

    //@}

    //! @name Attributes
    //@{

    //! The number of elements in this Map.
    global_size_t getGlobalNumElements() const { XPETRA_MONITOR("EpetraMapT::getGlobalNumElements"); return map_->NumGlobalElements64(); }

    //! The number of elements belonging to the calling process.
    size_t getNodeNumElements() const { XPETRA_MONITOR("EpetraMapT::getNodeNumElements"); return map_->NumMyElements(); }

    //! The index base for this Map.
    GlobalOrdinal getIndexBase() const { XPETRA_MONITOR("EpetraMapT::getIndexBase"); return (GlobalOrdinal) map_->IndexBase64(); }

    //! The minimum local index.
    LocalOrdinal getMinLocalIndex() const { XPETRA_MONITOR("EpetraMapT::getMinLocalIndex"); return map_->MinLID(); }

    //! The maximum local index on the calling process.
    LocalOrdinal getMaxLocalIndex() const { XPETRA_MONITOR("EpetraMapT::getMaxLocalIndex"); return map_->MaxLID(); }

    //! The minimum global index owned by the calling process.
    GlobalOrdinal getMinGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMinGlobalIndex"); return (GlobalOrdinal) map_->MinMyGID64(); }

    //! The maximum global index owned by the calling process.
    GlobalOrdinal getMaxGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMaxGlobalIndex"); return (GlobalOrdinal) map_->MaxMyGID64(); }

    //! The minimum global index over all processes in the communicator.
    GlobalOrdinal getMinAllGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMinAllGlobalIndex"); return (GlobalOrdinal) map_->MinAllGID64(); }

    //! The maximum global index over all processes in the communicator.
    GlobalOrdinal getMaxAllGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMaxAllGlobalIndex"); return (GlobalOrdinal) map_->MaxAllGID64(); }

    //! The local index corresponding to the given global index.
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const { XPETRA_MONITOR("EpetraMapT::getLocalElement"); return map_->LID(globalIndex); }

    //! Return the process ranks and corresponding local indices for the given global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const { XPETRA_MONITOR("EpetraMapT::getRemoteIndexList"); return toXpetra(map_->RemoteIDList(static_cast<int>(GIDList.size()), GIDList.getRawPtr(), nodeIDList.getRawPtr(), LIDList.getRawPtr())); }

    //! Return the process ranks for the given global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const { XPETRA_MONITOR("EpetraMapT::getRemoteIndexList"); return toXpetra(map_->RemoteIDList(static_cast<int>(GIDList.size()), GIDList.getRawPtr(), nodeIDList.getRawPtr(), 0)); }

    //! Return a view of the global indices owned by this process.
    Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const { XPETRA_MONITOR("EpetraMapT::getNodeElementList"); return ArrayView< const int >(map_->MyGlobalElements(), map_->NumMyElements());  }
    //@}

    //! @name Boolean tests
    //@{

    //! Whether the given local index is valid for this Map on this process.
    bool isNodeLocalElement(LocalOrdinal localIndex) const { XPETRA_MONITOR("EpetraMapT::isNodeLocalElement"); return map_->MyLID(localIndex); }

    //! Whether the given global index is valid for this Map on this process.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const { XPETRA_MONITOR("EpetraMapT::isNodeGlobalElement"); return map_->MyGID(globalIndex); }

    //! True if this Map is distributed contiguously, else false.
    bool isContiguous() const { XPETRA_MONITOR("EpetraMapT::isContiguous"); return map_->LinearMap(); }

    //! Whether this Map is globally distributed or locally replicated.
    bool isDistributed() const { XPETRA_MONITOR("EpetraMapT::isDistributed"); return map_->DistributedGlobal(); }

    //! True if and only if map is compatible with this Map.
    bool isCompatible(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { XPETRA_MONITOR("EpetraMapT::isCompatible"); return map_->PointSameAs(toEpetra<GlobalOrdinal,Node>(map)); }

    //! True if and only if map is identical to this Map.
    bool isSameAs(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { XPETRA_MONITOR("EpetraMapT::isSameAs"); return map_->SameAs(toEpetra<GlobalOrdinal,Node>(map)); }

    //@}

    //! @name
    //@{

    //! Get this Map's Comm object.
    Teuchos::RCP< const Teuchos::Comm< int > > getComm() const { XPETRA_MONITOR("EpetraMapT::getComm"); return toXpetra(map_->Comm()); }

    //! Get this Map's Node object.
    Teuchos::RCP< Node > getNode() const {
      XPETRA_MONITOR("EpetraMapT<GlobalOrdinal>::getNode");
      return Teuchos::rcp (new Node);
    }

    //@}

    //! @name
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const {
      XPETRA_MONITOR("EpetraMapT::description");

      // This implementation come from Tpetra_Map_def.hpp (without modification)
      std::ostringstream oss;
      oss << Teuchos::Describable::description();
      oss << "{getGlobalNumElements() = " << getGlobalNumElements()
          << ", getNodeNumElements() = " << getNodeNumElements()
          << ", isContiguous() = " << isContiguous()
          << ", isDistributed() = " << isDistributed()
          << "}";
      return oss.str();
    }

    //! Print this object with the given verbosity level to the given Teuchos::FancyOStream.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {
      XPETRA_MONITOR("EpetraMapT::describe");

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
      Teuchos::ArrayView<const GlobalOrdinal> myEntries = getNodeElementList();
      int myImageID = comm_->getRank();
      int numImages = comm_->getSize();

      Teuchos::EVerbosityLevel vl = verbLevel;
      if (vl == VERB_DEFAULT) vl = VERB_LOW;

      size_t width = 1;
      for (size_t dec=10; dec<getGlobalNumElements(); dec *= 10) {
        ++width;
      }
      width = ::std::max<size_t>(width, (size_t) 12) + 2; // casting to size_t to avoid ambiguity error when compiling Sacado.

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

    //! @name
    //@{

    //! Return a new Map with processes with zero elements removed.
    RCP<const Map<int,GlobalOrdinal,Node> > removeEmptyProcesses() const {
      const Epetra_BlockMap * NewMap = map_->RemoveEmptyProcesses();
       if (!NewMap) {
         return Teuchos::null;
       } else {
         const RCP< const Map<int, GlobalOrdinal, Node> >  NewMapX = toXpetra<GlobalOrdinal, Node>(*NewMap);
         delete NewMap;   // NOTE: toXpetra *copys* the epetra map rather than wrapping it, so we have to delete NewMap to avoid a memory leak.
         return NewMapX;
       }
    }

    //! Replace this Map's communicator with a subset communicator.
    RCP<const Map<int,GlobalOrdinal,Node> > replaceCommWithSubset(const Teuchos::RCP< const Teuchos::Comm< int > > &newComm) const {
      throw std::runtime_error("Xpetra::EpetraMapT::replaceCommWithSubset has not yet been implemented.");
      TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
    }

    //@}

    //! Return the global index for a given local index.  Note that this returns -1 if not found on this processor.  (This is different than Epetra's behavior!)
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const {
      XPETRA_MONITOR("EpetraMapT::getGlobalElement");

      GlobalOrdinal gid = (GlobalOrdinal) map_->GID64(localIndex);
      if (gid == map_->IndexBase64()-1) return (-1);
      else                              return (gid);
     }

    //! @name Xpetra specific
    //@{

    //! Destructor.
    virtual ~EpetraMapT() {}

    //! EpetraMapT constructor to wrap a Epetra_Map object
    EpetraMapT(const Teuchos::RCP<const Epetra_BlockMap> &map)
      : map_(map) {
      TEUCHOS_TEST_FOR_EXCEPTION(!map->GlobalIndicesIsType<GlobalOrdinal>(), std::runtime_error, "Xpetra::EpetraMapT: GlobalOrdinal mismatch.");
    }

    //! Get the library used by this object (Epetra or Epetra?)
    UnderlyingLib lib() const { return Xpetra::UseEpetra; }

    //! Get the underlying Epetra map
    //const RCP< const Epetra_Map > & getEpetra_Map() const { return map_; }
    const Epetra_BlockMap& getEpetra_BlockMap() const { return *map_; }
    const Epetra_Map& getEpetra_Map() const { return (Epetra_Map &)*map_; } // Ugly, but the same is done in Epetra_CrsMatrix.h to get the map.

    //@}

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    using local_map_type = typename Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type;
    /// \brief Get the local Map for Kokkos kernels.
    local_map_type getLocalMap () const {
      if (isInitializedLocalMap_)
        return localMap_;

      typedef GlobalOrdinal GO;
      typedef LocalOrdinal  LO;

      typedef typename Node::device_type DeviceType;

      typedef Tpetra::Details::FixedHashTable<GO, LO, DeviceType>   glMapType;
      typedef ::Kokkos::View<GO*, ::Kokkos::LayoutLeft, DeviceType> lgMapType;

      GO   indexBase  = getIndexBase();
      GO   myMinGID   = getMinGlobalIndex();
      GO   myMaxGID   = getMaxGlobalIndex();
      bool contiguous = isContiguous();

      size_t numLocalElements = map_->NumMyElements();

      GO firstContiguousGID, lastContiguousGID;
      glMapType glMap;
      lgMapType lgMap;
      if (numLocalElements) {
        if (contiguous) {
          // Don't need to initialize glMap and lgMap
          firstContiguousGID = myMinGID;
          lastContiguousGID  = myMaxGID+1;

        } else {

          auto GIDs = getNodeElementList();

          lgMap = lgMapType("lgMap", numLocalElements);

          firstContiguousGID = GIDs[0];
          lastContiguousGID  = firstContiguousGID+1;

          size_t i = 1;
          lgMap(0) = firstContiguousGID;
          for ( ; i < numLocalElements; ++i) {
            const GO curGID = GIDs[i];
            const LO curLid = Teuchos::as<LO> (i);

            if (lastContiguousGID != curGID)
              break;

            // Add the entry to the LID->GID table only after we know that
            // the current GID is in the initial contiguous sequence, so
            // that we don't repeat adding it in the first iteration of
            // the loop below over the remaining noncontiguous GIDs.
            lgMap(curLid) = curGID;
            ++lastContiguousGID;
          }
          --lastContiguousGID;

          ::Kokkos::View<GO*, ::Kokkos::LayoutLeft, DeviceType>
            nonContigGIDs ("nonContigGIDs", numLocalElements-i);
          // FIXME_KOKKOS: relies on UVM
          for (size_t idx = 0; idx < numLocalElements - i; idx++)
            nonContigGIDs(idx) = GIDs[i+idx];

          glMap = glMapType (nonContigGIDs,
                             firstContiguousGID,
                             lastContiguousGID,
                             static_cast<LO> (i));

          for ( ; i < numLocalElements; ++i) {
            const GO curGID = GIDs[i];
            const LO curLid = Teuchos::as<LO> (i);

            lgMap(curLid) = curGID;
          }
        }

      } else {
        // No elements
        firstContiguousGID = indexBase+1;
        lastContiguousGID  = indexBase;
      }
      localMap_ = local_map_type(glMap, lgMap, indexBase, myMinGID, myMaxGID, firstContiguousGID, lastContiguousGID, numLocalElements, contiguous);

      isInitializedLocalMap_ = true;

      return localMap_;
    }

  private:
    mutable local_map_type localMap_;
    mutable bool isInitializedLocalMap_ = false; // It's OK to use C++11 when Tpetra is enabled
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif
#endif

  protected:

    RCP<const Epetra_BlockMap> map_;
}; // EpetraMapT class
#endif // #ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

// specialization on GO=long long and EpetraNode
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
  template<>
  class EpetraMapT<long long, EpetraNode>
    : public virtual Map<int, long long, EpetraNode>
  {
    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

    //! @name Constructors and destructor
    //@{

    //! Constructor with Tpetra-defined contiguous uniform distribution.
    EpetraMapT(global_size_t numGlobalElements,
               GlobalOrdinal indexBase,
               const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
               LocalGlobal lg=GloballyDistributed,
               const Teuchos::RCP< Node > &node = Teuchos::rcp(new Node)) {
      // This test come from Tpetra (Epetra doesn't check if numGlobalElements,indexBase are equivalent across images).
      // In addition, for the test TEST_THROW(M map((myImageID == 0 ? GSTI : 0),0,comm), std::invalid_argument), only one node throw an exception and there is a dead lock.
      std::string errPrefix;
      errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,indexBase,comm,lOrG): ";

      if (lg == GloballyDistributed) {
        const int myImageID = comm->getRank();

        // check that numGlobalElements,indexBase is equivalent across images
        global_size_t rootNGE = numGlobalElements;
        GlobalOrdinal rootIB  = indexBase;
        Teuchos::broadcast<int,global_size_t>(*comm,0,&rootNGE);
        Teuchos::broadcast<int,GlobalOrdinal>(*comm,0,&rootIB);
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
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal must be the same on all nodes (examine node " << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "indexBase must be the same on all nodes (examine node " << globalChecks[0] << ").");
          }
          else {
            // logic error on our part
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               errPrefix << "logic error. Please contact the Tpetra team.");
          }
        }
      }

      // Note: validity of numGlobalElements checked by Epetra.

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(static_cast<GlobalOrdinal>(numGlobalElements), 1, indexBase, *toEpetra(comm))))));
    }

    //! Constructor with a user-defined contiguous distribution.
    EpetraMapT(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::rcp(new Node)) {
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
        GlobalOrdinal rootIB = indexBase;
        Teuchos::broadcast<int,GlobalOrdinal>(*comm,0,&rootIB);   // broadcast one ordinal from node 0
        if (indexBase != rootIB) {
          localChecks[0] = myImageID;
          localChecks[1] = 4;
        }
        // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
        // this will be -1 if all procs passed
        Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
        if (globalChecks[0] != -1) {
          if (globalChecks[1] == 1) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numLocal is not valid on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal is not valid on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 3) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "numGlobal doesn't match sum of numLocal (== "
                               << global_sum << ") on at least one node (possibly node "
                               << globalChecks[0] << ").");
          }
          else if (globalChecks[1] == 4) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               errPrefix << "indexBase is not the same on all nodes (examine node "
                               << globalChecks[0] << ").");
          }
          else {
            // logic error on my part
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               errPrefix << "logic error. Please contact the Tpetra team.");
          }
        }

      }

      // set numGlobalElements
      if (numGlobalElements == GSTI) {
        numGlobalElements = global_sum;}

      IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(static_cast<GlobalOrdinal>(numGlobalElements), numLocalElements, 1, indexBase, *toEpetra(comm))))));
    }

    //! Constructor with user-defined arbitrary (possibly noncontiguous) distribution.
    EpetraMapT(global_size_t numGlobalElements,
               const Teuchos::ArrayView< const GlobalOrdinal > &elementList,
               GlobalOrdinal indexBase,
               const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
               const Teuchos::RCP< Node > &node = Teuchos::rcp(new Node)) {
      if (numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid()) {
        IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(-1, elementList.size(), elementList.getRawPtr(), 1, indexBase, *toEpetra(comm))))));
      } else {
        IF_EPETRA_EXCEPTION_THEN_THROW_GLOBAL_INVALID_ARG((map_ = (rcp(new Epetra_BlockMap(numGlobalElements, elementList.size(), elementList.getRawPtr(), 1, indexBase, *toEpetra(comm))))));
      }
    }

    //@}

    //! @name Attributes
    //@{

    //! The number of elements in this Map.
    global_size_t getGlobalNumElements() const { XPETRA_MONITOR("EpetraMapT::getGlobalNumElements"); return map_->NumGlobalElements64(); }

    //! The number of elements belonging to the calling process.
    size_t getNodeNumElements() const { XPETRA_MONITOR("EpetraMapT::getNodeNumElements"); return map_->NumMyElements(); }

    //! The index base for this Map.
    GlobalOrdinal getIndexBase() const { XPETRA_MONITOR("EpetraMapT::getIndexBase"); return (GlobalOrdinal) map_->IndexBase64(); }

    //! The minimum local index.
    LocalOrdinal getMinLocalIndex() const { XPETRA_MONITOR("EpetraMapT::getMinLocalIndex"); return map_->MinLID(); }

    //! The maximum local index on the calling process.
    LocalOrdinal getMaxLocalIndex() const { XPETRA_MONITOR("EpetraMapT::getMaxLocalIndex"); return map_->MaxLID(); }

    //! The minimum global index owned by the calling process.
    GlobalOrdinal getMinGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMinGlobalIndex"); return (GlobalOrdinal) map_->MinMyGID64(); }

    //! The maximum global index owned by the calling process.
    GlobalOrdinal getMaxGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMaxGlobalIndex"); return (GlobalOrdinal) map_->MaxMyGID64(); }

    //! The minimum global index over all processes in the communicator.
    GlobalOrdinal getMinAllGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMinAllGlobalIndex"); return (GlobalOrdinal) map_->MinAllGID64(); }

    //! The maximum global index over all processes in the communicator.
    GlobalOrdinal getMaxAllGlobalIndex() const { XPETRA_MONITOR("EpetraMapT::getMaxAllGlobalIndex"); return (GlobalOrdinal) map_->MaxAllGID64(); }

    //! The local index corresponding to the given global index.
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const { XPETRA_MONITOR("EpetraMapT::getLocalElement"); return map_->LID(globalIndex); }

    //! Return the process ranks and corresponding local indices for the given global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const { XPETRA_MONITOR("EpetraMapT::getRemoteIndexList"); return toXpetra(map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), LIDList.getRawPtr())); }

    //! Return the process ranks for the given global indices.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const { XPETRA_MONITOR("EpetraMapT::getRemoteIndexList"); return toXpetra(map_->RemoteIDList(GIDList.size(), GIDList.getRawPtr(), nodeIDList.getRawPtr(), 0)); }

    //! Return a view of the global indices owned by this process.
    Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const { XPETRA_MONITOR("EpetraMapT::getNodeElementList"); return ArrayView< const long long >(map_->MyGlobalElements64(), map_->NumMyElements()); }
    //@}

    //! @name Boolean tests
    //@{

    //! Whether the given local index is valid for this Map on this process.
    bool isNodeLocalElement(LocalOrdinal localIndex) const { XPETRA_MONITOR("EpetraMapT::isNodeLocalElement"); return map_->MyLID(localIndex); }

    //! Whether the given global index is valid for this Map on this process.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const { XPETRA_MONITOR("EpetraMapT::isNodeGlobalElement"); return map_->MyGID(globalIndex); }

    //! True if this Map is distributed contiguously, else false.
    bool isContiguous() const { XPETRA_MONITOR("EpetraMapT::isContiguous"); return map_->LinearMap(); }

    //! Whether this Map is globally distributed or locally replicated.
    bool isDistributed() const { XPETRA_MONITOR("EpetraMapT::isDistributed"); return map_->DistributedGlobal(); }

    //! True if and only if map is compatible with this Map.
    bool isCompatible(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { XPETRA_MONITOR("EpetraMapT::isCompatible"); return map_->PointSameAs(toEpetra<GlobalOrdinal,Node>(map)); }

    //! True if and only if map is identical to this Map.
    bool isSameAs(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { XPETRA_MONITOR("EpetraMapT::isSameAs"); return map_->SameAs(toEpetra<GlobalOrdinal,Node>(map)); }

    //@}

    //! @name
    //@{

    //! Get this Map's Comm object.
    Teuchos::RCP< const Teuchos::Comm< int > > getComm() const { XPETRA_MONITOR("EpetraMapT::getComm"); return toXpetra(map_->Comm()); }

    //! Get this Map's Node object.
    Teuchos::RCP< Node > getNode() const {
      XPETRA_MONITOR("EpetraMapT<GlobalOrdinal>::getNode");
      return Teuchos::rcp (new Node);
    }

    //@}

    //! @name
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const {
      XPETRA_MONITOR("EpetraMapT::description");

      // This implementation come from Tpetra_Map_def.hpp (without modification)
      std::ostringstream oss;
      oss << Teuchos::Describable::description();
      oss << "{getGlobalNumElements() = " << getGlobalNumElements()
          << ", getNodeNumElements() = " << getNodeNumElements()
          << ", isContiguous() = " << isContiguous()
          << ", isDistributed() = " << isDistributed()
          << "}";
      return oss.str();
    }

    //! Print this object with the given verbosity level to the given Teuchos::FancyOStream.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {
      XPETRA_MONITOR("EpetraMapT::describe");

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
      Teuchos::ArrayView<const GlobalOrdinal> myEntries = getNodeElementList();
      int myImageID = comm_->getRank();
      int numImages = comm_->getSize();

      Teuchos::EVerbosityLevel vl = verbLevel;
      if (vl == VERB_DEFAULT) vl = VERB_LOW;

      size_t width = 1;
      for (size_t dec=10; dec<getGlobalNumElements(); dec *= 10) {
        ++width;
      }
      width = ::std::max<size_t>(width, (size_t) 12) + 2; // casting to size_t to avoid ambiguity error when compiling Sacado.

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

    //! @name
    //@{

    //! Return a new Map with processes with zero elements removed.
    RCP<const Map<int,GlobalOrdinal,Node> > removeEmptyProcesses() const {
      const Epetra_BlockMap * NewMap = map_->RemoveEmptyProcesses();
       if (!NewMap) {
         return Teuchos::null;
       } else {
         const RCP< const Map<int, GlobalOrdinal, Node> >  NewMapX = toXpetra<GlobalOrdinal, Node>(*NewMap);
         delete NewMap;   // NOTE: toXpetra *copys* the epetra map rather than wrapping it, so we have to delete NewMap to avoid a memory leak.
         return NewMapX;
       }
    }

    //! Replace this Map's communicator with a subset communicator.
    RCP<const Map<int,GlobalOrdinal,Node> > replaceCommWithSubset(const Teuchos::RCP< const Teuchos::Comm< int > > &newComm) const {
      throw std::runtime_error("Xpetra::EpetraMapT::replaceCommWithSubset has not yet been implemented.");
      // return Teuchos::null; // unreachable
    }

    //@}

    //! Return the global index for a given local index.  Note that this returns -1 if not found on this processor.  (This is different than Epetra's behavior!)
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const {
      XPETRA_MONITOR("EpetraMapT::getGlobalElement");

      GlobalOrdinal gid = (GlobalOrdinal) map_->GID64(localIndex);
      if (gid == map_->IndexBase64()-1) return (-1);
      else                              return (gid);
     }

    //! @name Xpetra specific
    //@{

    //! Destructor.
    virtual ~EpetraMapT() {}

    //! EpetraMapT constructor to wrap a Epetra_Map object
    EpetraMapT(const Teuchos::RCP<const Epetra_BlockMap> &map)
      : map_(map) {
      TEUCHOS_TEST_FOR_EXCEPTION(!map->GlobalIndicesIsType<GlobalOrdinal>(), std::runtime_error, "Xpetra::EpetraMapT: GlobalOrdinal mismatch.");
    }

    //! Get the library used by this object (Epetra or Epetra?)
    UnderlyingLib lib() const { return Xpetra::UseEpetra; }

    //! Get the underlying Epetra map
    //const RCP< const Epetra_Map > & getEpetra_Map() const { return map_; }
    const Epetra_BlockMap& getEpetra_BlockMap() const { return *map_; }
    const Epetra_Map& getEpetra_Map() const { return (Epetra_Map &)*map_; } // Ugly, but the same is done in Epetra_CrsMatrix.h to get the map.

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    using local_map_type = typename Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type;
    /// \brief Get the local Map for Kokkos kernels.
    local_map_type getLocalMap () const {
      if (isInitializedLocalMap_)
        return localMap_;

      typedef GlobalOrdinal GO;
      typedef LocalOrdinal  LO;

      typedef typename Node::device_type DeviceType;

      typedef Tpetra::Details::FixedHashTable<GO, LO, DeviceType>   glMapType;
      typedef ::Kokkos::View<GO*, ::Kokkos::LayoutLeft, DeviceType> lgMapType;

      GO   indexBase  = getIndexBase();
      GO   myMinGID   = getMinGlobalIndex();
      GO   myMaxGID   = getMaxGlobalIndex();
      bool contiguous = isContiguous();

      size_t numLocalElements = map_->NumMyElements();

      GO firstContiguousGID, lastContiguousGID;
      glMapType glMap;
      lgMapType lgMap;
      if (numLocalElements) {
        if (contiguous) {
          // Don't need to initialize glMap and lgMap
          firstContiguousGID = myMinGID;
          lastContiguousGID  = myMaxGID+1;

        } else {

          auto GIDs = getNodeElementList();

          lgMap = lgMapType("lgMap", numLocalElements);

          firstContiguousGID = GIDs[0];
          lastContiguousGID  = firstContiguousGID+1;

          size_t i = 1;
          lgMap(0) = firstContiguousGID;
          for ( ; i < numLocalElements; ++i) {
            const GO curGID = GIDs[i];
            const LO curLid = Teuchos::as<LO> (i);

            if (lastContiguousGID != curGID)
              break;

            // Add the entry to the LID->GID table only after we know that
            // the current GID is in the initial contiguous sequence, so
            // that we don't repeat adding it in the first iteration of
            // the loop below over the remaining noncontiguous GIDs.
            lgMap(curLid) = curGID;
            ++lastContiguousGID;
          }
          --lastContiguousGID;

          ::Kokkos::View<GO*, ::Kokkos::LayoutLeft, DeviceType>
            nonContigGIDs ("nonContigGIDs", numLocalElements-i);
          // FIXME_KOKKOS: relies on UVM
          for (size_t idx = 0; idx < numLocalElements - i; idx++)
            nonContigGIDs(idx) = GIDs[i+idx];

          glMap = glMapType (nonContigGIDs,
                             firstContiguousGID,
                             lastContiguousGID,
                             static_cast<LO> (i));

          for ( ; i < numLocalElements; ++i) {
            const GO curGID = GIDs[i];
            const LO curLid = Teuchos::as<LO> (i);

            lgMap(curLid) = curGID;
          }
        }

      } else {
        // No elements
        firstContiguousGID = indexBase+1;
        lastContiguousGID  = indexBase;
      }
      localMap_ = local_map_type(glMap, lgMap, indexBase, myMinGID, myMaxGID, firstContiguousGID, lastContiguousGID, numLocalElements, contiguous);

      isInitializedLocalMap_ = true;

      return localMap_;
    }

  private:
    mutable local_map_type localMap_;
    mutable bool isInitializedLocalMap_ = false; // It's OK to use C++11 when Tpetra is enabled
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif
#endif

    //@}

  protected:

    RCP<const Epetra_BlockMap> map_;
}; // EpetraMapT class
#endif // #ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

} // Xpetra namespace

#endif // XPETRA_EPETRAMAP_HPP
