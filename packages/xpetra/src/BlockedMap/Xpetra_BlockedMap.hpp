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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_HPP_
#define PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_HPP_

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Map.hpp"
//#include "Xpetra_MapFactory.hpp"
#include "Xpetra_ImportFactory.hpp"
//#include "Xpetra_MapUtils.hpp"

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Vector, needed to prevent circular inclusions
//  template<class S, class LO, class GO, class N> class Vector;
  template<class LO, class GO, class N> class MapFactory;
#endif

  template <class LocalOrdinal = Map<>::local_ordinal_type,
            class GlobalOrdinal = typename Map<LocalOrdinal>::global_ordinal_type,
            class Node = typename Map<LocalOrdinal, GlobalOrdinal>::node_type>
  class BlockedMap
    : public Map< LocalOrdinal, GlobalOrdinal, Node >
  {
  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

  private:
#undef XPETRA_BLOCKEDMAP_SHORT
#include "Xpetra_UseShortNamesOrdinal.hpp"

  public:
    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor
    /*!
     */
    BlockedMap() {
      bThyraMode_ = false;
    }

    //! BlockedMap basic constructor
    //!
    //! @param[in] fullmap Full map containing all GIDs throughout the full vector. This parameter is only important if bThyraMode == false (see below)
    //! @param[in] maps    Vector containing submaps. The set of all GIDs stored in the submaps should be the same than stored in fullmap, if bThyraMode == false. In Thyra mode, the submaps should contain consecutive GIDs starting with 0 in each submap.
    //! @param[in] bThyraMode Flag which allows to switch between generating a MapExtractor in Thyra mode or Xpetra mode
    //!
    //! In Thyra mode, fullmap is not important as a fullmap with unique blocked GIDs is automatically generated which map the GIDs of the submaps
    //! to uniquely defined GIDs in the fullmap. The user has to provide a fullmap in Thyra mode to specify the underlying linear algebra library
    //! (Epetra or Tpetra) and some other map information (e.g. indexBase). This could be fixed.
    //!
    //! In Xpetra mode, the fullmap has to be the same as the union of the GIDs stored in the submaps in maps. The intersection of the GIDs of the sub-
    //! maps in maps must be empty.
    BlockedMap(const RCP<const Map>& fullmap, const std::vector<RCP<const Map> >& maps, bool bThyraMode = false) {
      bThyraMode_ = bThyraMode;

      if(bThyraMode == false) {
        // use Xpetra-style numbering for sub-block maps
        // That is, all sub-block maps have unique GIDs which may not be contiguous and start with GIDs different than zero.

        // plausibility check
        size_t numAllElements = 0;
        for(size_t v = 0; v < maps.size(); ++v) {
          numAllElements += maps[v]->getGlobalNumElements();
        }
        TEUCHOS_TEST_FOR_EXCEPTION(fullmap->getGlobalNumElements() != numAllElements, std::logic_error,
                                   "logic error. full map and sub maps have not same number of elements (" << fullmap->getGlobalNumElements() <<" versus " << numAllElements << "). We cannot build MapExtractor with Xpetra-style numbering. Please make sure that you want Xpetra-style numbering instead of Thyra-style numbering.");

        fullmap_ = fullmap;
        maps_ = maps;
      } else {
        //std::cout << "Create Map Extractor in Thyra Mode!!! " << std::endl;
        // use Thyra-style numbering for sub-block maps
        // That is, all sub-block maps start with zero as GID and are contiguous

        // plausibility check
        for(size_t v = 0; v < maps.size(); ++v) {
          TEUCHOS_TEST_FOR_EXCEPTION(maps[v]->getMinAllGlobalIndex() != 0, std::logic_error,
                                             "logic error. When using Thyra-style numbering all sub-block maps must start with zero as GID. Map block " << v << " starts with GID " << maps[v]->getMinAllGlobalIndex());
        }

        // store submaps in Thyra-style ordering
        thyraMaps_ = maps;

        // get offsets
        std::vector<GlobalOrdinal> gidOffsets(maps.size(),0);
        for(size_t v = 1; v < maps.size(); ++v) {
          gidOffsets[v] = maps[v-1]->getMaxAllGlobalIndex() + gidOffsets[v-1] + 1;
        }

        // build submaps
        maps_.resize(maps.size());
        std::vector<GlobalOrdinal> fullMapGids;
        const GO INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
        for(size_t v = 0; v < maps.size(); ++v) {
          size_t myNumElements = maps[v]->getNodeNumElements();
          std::vector<GlobalOrdinal> subMapGids(myNumElements,0);
          for (LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(myNumElements); ++l) {
            GlobalOrdinal myGid = maps[v]->getGlobalElement(l);
            subMapGids[l] = myGid + gidOffsets[v];
            fullMapGids.push_back(myGid + gidOffsets[v]);
          }
          //std::sort(subMapGids.begin(), subMapGids.end());
          //subMapGids.erase(std::unique(subMapGids.begin(), subMapGids.end()), subMapGids.end());

          Teuchos::ArrayView<GlobalOrdinal> subMapGidsView(&subMapGids[0], subMapGids.size());
          Teuchos::RCP<Map> mySubMap = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(maps[v]->lib(), INVALID, subMapGidsView, maps[v]->getIndexBase(), maps[v]->getComm());
          maps_[v] = mySubMap;
        }

        //const GO INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
        //std::sort(coarseMapGids.begin(), coarseMapGids.end());
        //coarseMapGids.erase(std::unique(coarseMapGids.begin(), coarseMapGids.end()), coarseMapGids.end());
        //Teuchos::ArrayView<GO> coarseMapGidsView(&coarseMapGids[0], coarseMapGids.size());
        //std::sort(fullMapGids.begin(), fullMapGids.end());
        //fullMapGids.erase(std::unique(fullMapGids.begin(), fullMapGids.end()), fullMapGids.end());

        Teuchos::ArrayView<GlobalOrdinal> fullMapGidsView(&fullMapGids[0], fullMapGids.size());
        fullmap_ = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(fullmap->lib(), INVALID, fullMapGidsView, fullmap->getIndexBase(), fullmap->getComm());

        // plausibility check
        size_t numAllElements = 0;
        for(size_t v = 0; v < maps_.size(); ++v) {
          numAllElements += maps_[v]->getGlobalNumElements();
        }
        TEUCHOS_TEST_FOR_EXCEPTION(fullmap_->getGlobalNumElements() != numAllElements, std::logic_error,
                                   "logic error. full map and sub maps have not same number of elements. This cannot be. Please report the bug to the Xpetra developers!");
      }

      // build importers for sub maps
      importers_.resize(maps_.size());
      for (unsigned i = 0; i < maps_.size(); ++i)
        if (maps[i] != null)
          importers_[i] = Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(fullmap_, maps_[i]);
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, std::logic_error,
                                 "logic error. full map and sub maps are inconsistently distributed over the processors.");

    }

    //! Expert constructor for Thyra maps
    BlockedMap(const std::vector<RCP<const Map> >& maps, const std::vector<RCP<const Map> >& thyramaps) {
      bThyraMode_ = true;

      // plausibility check
      TEUCHOS_TEST_FOR_EXCEPTION(thyramaps.size() != maps.size(), std::logic_error, "logic error. The number of submaps must be identical!");
      for(size_t v = 0; v < thyramaps.size(); ++v) {
        TEUCHOS_TEST_FOR_EXCEPTION(thyramaps[v]->getMinAllGlobalIndex() != 0, std::logic_error,
                                           "logic error. When using Thyra-style numbering all sub-block maps must start with zero as GID.");
        XPETRA_TEST_FOR_EXCEPTION(thyramaps[v]->getNodeNumElements() != maps[v]->getNodeNumElements(), std::logic_error,
                                           "logic error. The size of the submaps must be identical (same distribution, just different GIDs)");
      }

      // store user-provided maps and thyramaps
      thyraMaps_ = thyramaps;
      maps_      = maps;

      fullmap_ = this->concatenateMaps(maps);

      // plausibility check
      size_t numAllElements = 0;
      for(size_t v = 0; v < maps_.size(); ++v) {
        numAllElements += maps_[v]->getGlobalNumElements();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(fullmap_->getGlobalNumElements() != numAllElements, std::logic_error,
                                 "logic error. full map and sub maps have not same number of elements. This cannot be. Please report the bug to the Xpetra developers!");

      // build importers for sub maps
      importers_.resize(maps_.size());
      for (unsigned i = 0; i < maps_.size(); ++i)
        if (maps[i] != null)
          importers_[i] = Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(fullmap_, maps_[i]);
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, std::logic_error,
                                 "logic error. full map and sub maps are inconsistently distributed over the processors.");
    }

    //! copy constructor
    BlockedMap(const BlockedMap& input) {
      bThyraMode_ = input.getThyraMode();
      fullmap_ = Teuchos::null;
      maps_.resize(input.getNumMaps(), Teuchos::null);
      thyraMaps_.resize(input.getNumMaps(), Teuchos::null);
      this->assign(input);
    }

    //! Destructor.
    virtual ~BlockedMap() {

      // make sure all RCP's are freed
      for(size_t v = 0; v < maps_.size(); ++v) {
        maps_[v] = Teuchos::null;
        if(bThyraMode_ == true)
          thyraMaps_[v] = Teuchos::null;
        importers_[v] = Teuchos::null;
      }

      fullmap_ = Teuchos::null;
    }

    //! @name Attributes
    //@{

    //! The number of elements in this Map.
    virtual global_size_t getGlobalNumElements() const { return fullmap_->getGlobalNumElements(); };

    //! The number of elements belonging to the calling process.
    virtual size_t getNodeNumElements() const { return fullmap_->getNodeNumElements(); };

    //! The index base for this Map.
    virtual GlobalOrdinal getIndexBase() const  { return fullmap_->getIndexBase(); };

    //! The minimum local index.
    virtual LocalOrdinal getMinLocalIndex() const  { return fullmap_->getMinLocalIndex(); };

    //! The maximum local index on the calling process.
    virtual LocalOrdinal getMaxLocalIndex() const  { return fullmap_->getMaxLocalIndex(); };

    //! The minimum global index owned by the calling process.
    virtual GlobalOrdinal getMinGlobalIndex() const  { return fullmap_->getMinGlobalIndex(); };

    //! The maximum global index owned by the calling process.
    virtual GlobalOrdinal getMaxGlobalIndex() const  { return fullmap_->getMaxGlobalIndex(); };

    //! The minimum global index over all processes in the communicator.
    virtual GlobalOrdinal getMinAllGlobalIndex() const  { return fullmap_->getMinAllGlobalIndex(); };

    //! The maximum global index over all processes in the communicator.
    virtual GlobalOrdinal getMaxAllGlobalIndex() const  { return fullmap_->getMaxAllGlobalIndex(); };

    //! The local index corresponding to the given global index.
    virtual LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const  { return fullmap_->getLocalElement(globalIndex); };

    //! The global index corresponding to the given local index.
    virtual GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const { return fullmap_->getGlobalElement(localIndex); };

    //! Return the process ranks and corresponding local indices for the given global indices.
    virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &/* GIDList */, const Teuchos::ArrayView< int > &/* nodeIDList */, const Teuchos::ArrayView< LocalOrdinal > &/* LIDList */) const {
      throw Xpetra::Exceptions::RuntimeError("BlockedMap::getRemoteIndexList: routine not implemented.");
      TEUCHOS_UNREACHABLE_RETURN(IDNotPresent);
    };

    //! Return the process ranks for the given global indices.
    virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &/* GIDList */, const Teuchos::ArrayView< int > &/* nodeIDList */) const {
      throw Xpetra::Exceptions::RuntimeError("BlockedMap::getRemoteIndexList: routine not implemented.");
      TEUCHOS_UNREACHABLE_RETURN(IDNotPresent);
    };

    //! Return a view of the global indices owned by this process.
    virtual Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const {
      return fullmap_->getNodeElementList();
    };

    //@}

    //! @name Boolean tests
    //@{

    //! Whether the given local index is valid for this Map on this process.
    virtual bool isNodeLocalElement(LocalOrdinal localIndex) const {return fullmap_->isNodeLocalElement(localIndex);};

    //! Whether the given global index is valid for this Map on this process.
    virtual bool isNodeGlobalElement(GlobalOrdinal globalIndex) const {return fullmap_->isNodeGlobalElement(globalIndex);};

    //! True if this Map is distributed contiguously, else false.
    virtual bool isContiguous() const {
      throw Xpetra::Exceptions::RuntimeError("BlockedMap::isContiguous: routine not implemented.");
      TEUCHOS_UNREACHABLE_RETURN(false);
    };

    //! Whether this Map is globally distributed or locally replicated.
    virtual bool isDistributed() const {return fullmap_->isDistributed();};

    //! True if and only if map is compatible with this Map.
    virtual bool isCompatible(const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map) const {
      RCP<const Map> rcpMap = Teuchos::rcpFromRef(map);
      RCP<const BlockedMap> rcpBMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rcpMap);
      if(rcpBMap.is_null() == true) return false;

      for(size_t v = 0; v < maps_.size(); ++v) {
        bool bSame = getMap(v,false)->isCompatible(*(rcpBMap->getMap(v,false)));
        if (bSame == false) return false;
        if (bThyraMode_) {
          bSame = getMap(v,true)->isCompatible(*(rcpBMap->getMap(v,true)));
        }
      }
      return true;
    };

    //! True if and only if map is identical to this Map.
    virtual bool isSameAs(const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map) const {
      RCP<const Map> rcpMap = Teuchos::rcpFromRef(map);
      RCP<const BlockedMap> rcpBMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rcpMap);
      if(rcpBMap.is_null() == true) {
        // If this is a blocked map with > 1 blocks but "map" is a plain map they can't be the same
        if (this->getNumMaps() > 1)
          return false;
        // special case: this is a single blocked map and "map" is a plain map object
        bool bSame = getMap(0,bThyraMode_)->isSameAs(*rcpMap);
        return bSame;
      }

      for(size_t v = 0; v < maps_.size(); ++v) {
        bool bSame = getMap(v,false)->isSameAs(*(rcpBMap->getMap(v,false)));
        if (bSame == false) return false;
        if (bThyraMode_) {
          bSame = getMap(v,true)->isSameAs(*(rcpBMap->getMap(v,true)));
        if (bSame == false) return false;
        }
      }
      return true;
    };

    //@}

    //! @name
    //@{

    //! Get this Map's Comm object.
    virtual Teuchos::RCP< const Teuchos::Comm< int > > getComm() const { return fullmap_->getComm(); } ;

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    //! Get this Map's Node object.
    virtual Teuchos::RCP< Node > getNode() const { return fullmap_->getNode();};
#endif // TPETRA_ENABLE_DEPRECATED_CODE

    //@}

    /// \brief Assignment operator: Does a deep copy.
    ///
    /// The assignment operator does a deep copy, just like
    /// subclasses' copy constructors.
    ///
    /// \note This currently only works if both <tt>*this</tt> and the
    ///   input argument are instances of the same subclass.
    BlockedMap<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const BlockedMap& rhs) {
      assign (rhs); // dispatch to protected virtual method
      return *this;
    }

    //@}

    //! @name Attribute access functions
    //@{
    //! Local number of rows on the calling process.
    /*virtual size_t getLocalLength() const {
      throw Xpetra::Exceptions::RuntimeError("BlockedMap::getLocalLength: routine not implemented.");
      return 0;
    }*/

    //! Global number of rows in the multivector.
    /*virtual global_size_t getGlobalLength() const {
      throw Xpetra::Exceptions::RuntimeError("BlockedMap::getGlobalLength: routine not implemented.");
      return 0;
    }*/

    //! returns true if internally stored sub maps are in Thyra mode (i.e. start all with GIDs=0)
    virtual bool getThyraMode() const {
      return bThyraMode_;
    }
    //@}

    /** \name Maps */
    //@{

    //! Return a new Map with processes with zero elements removed.
    virtual RCP< const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > removeEmptyProcesses() const {
      throw Xpetra::Exceptions::RuntimeError("BlockedMap::removeEmptyProcesses: routine not implemented.");
    }


    //! Replace this Map's communicator with a subset communicator.
    virtual RCP< const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > replaceCommWithSubset(const Teuchos::RCP< const Teuchos::Comm< int > > &/* newComm */) const {
      throw Xpetra::Exceptions::RuntimeError("BlockedMap::replaceCommWithSubset: routine not implemented.");
    }

    //@}

    //! @name Xpetra specific
    //@{

    //! Get the library used by this object (Tpetra or Epetra?)
    virtual UnderlyingLib lib() const { return fullmap_->lib(); } ;

    // TODO: find a better solution for this hack
    // The problem is that EpetraMap, TpetraMap and StridedMap all inherit Map. To have proper toEpetra() we
    // need to understand the type of underlying matrix. But in src/Map we have no knowledge of StridedMaps, so
    // we cannot check for it by casting. This function allows us to avoid the restriction, as StridedMap redefines
    // it to return the base map.
    virtual RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getMap() const { return getFullMap(); }

    //@}


    /// number of partial maps
    size_t getNumMaps() const { return maps_.size(); }

    /// get the map
    /// returns the sub map i from list of sub maps
    /// depending on the parameter bThyraMode the sub map that is returned uses Thyra or Xpetra numbering
    /// Note: Thyra-numbering is only allowed if the BlockedMap is also constructed using Thyra numbering
    const RCP<const Map> getMap(size_t i, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION( i >= getNumMaps(), Xpetra::Exceptions::RuntimeError, "BlockedMap::getMap: tried to access block " << i << ", but BlockedMap has only " << getNumMaps() << " blocks! Block indices must be between 0 and " << getNumMaps() - 1 << ".");
      if(bThyraMode_ == true && bThyraMode == true)
        return thyraMaps_[i];
      XPETRA_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "BlockedMap::getMap: cannot return sub map in Thyra-style numbering if BlockedMap object is not created using Thyra-style numbered submaps.");
      return maps_[i];
    }

    /// get the importer between full map and partial map
    const RCP<Import> getImporter(size_t i) const {
      XPETRA_TEST_FOR_EXCEPTION( i >= getNumMaps(), Xpetra::Exceptions::RuntimeError, "BlockedMap::getImporter: tried to access block " << i << ", but BlockedMap has only " << getNumMaps() << " blocks! Block indices must be between 0 and " << getNumMaps() - 1 << ".");
      return  importers_[i];
    }

    /// the full map
    const RCP<const Map> getFullMap() const { return fullmap_; }

    /// returns map index in map extractor which contains GID
    size_t getMapIndexForGID(GlobalOrdinal gid) const {
      for (size_t i = 0; i < getNumMaps(); i++)
        if (getMap(i)->isNodeGlobalElement(gid) == true)
          return i;

      TEUCHOS_TEST_FOR_EXCEPTION(false, Xpetra::Exceptions::RuntimeError,
                                 "getMapIndexForGID: GID " << gid << " is not contained by a map in mapextractor." );
      return 0;
    }

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    using local_map_type = typename Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type;
    /// \brief Get the local Map for Kokkos kernels.
    local_map_type getLocalMap () const {
      return fullmap_->getLocalMap();
    }
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif
#endif

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    virtual std::string description() const {
      return std::string("BlockedMap");
    }

    //! Print the object with the given verbosity level to a FancyOStream.
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {
      out << "------------- Blocked Map -----------" << std::endl;
      out << description() << std::endl;
      out << "Thyra mode: " << getThyraMode() << std::endl;
      out << "No of submaps: " << getNumMaps() << std::endl;
      Teuchos::OSTab tab(out);
      for(size_t r = 0; r < getNumMaps(); r++) {
        std::cout << "MAP " << r << "/" << getNumMaps() - 1 << std::endl;
        getMap(r,false)->describe(out, verbLevel);
      }
      if(getThyraMode() == true) {
        for(size_t r = 0; r < getNumMaps(); r++) {
          std::cout << "Thyra MAP " << r << "/" << getNumMaps() - 1 << std::endl;
          getMap(r,true)->describe(out, verbLevel);
        }
      }
      out << "-------------------------------------" << std::endl;
    }


    //@}


  protected:
    /// \brief Implementation of the assignment operator (operator=);
    ///   does a deep copy.
    virtual void assign (const BlockedMap& input) {
      // TODO check implementation, simplify copy constructor
      bThyraMode_ = input.getThyraMode();

      fullmap_ = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(input.getFullMap(),1);

      maps_.resize(input.getNumMaps(), Teuchos::null);
      if(bThyraMode_ == true)
        thyraMaps_.resize(input.getNumMaps(), Teuchos::null);
      for(size_t i = 0; i < input.getNumMaps(); ++i) {
        maps_[i] = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(input.getMap(i,false),1);
        if(bThyraMode_ == true)
          thyraMaps_[i] = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(input.getMap(i,true),1);
      }

      // plausibility check
      size_t numAllElements = 0;
      for(size_t v = 0; v < maps_.size(); ++v) {
        numAllElements += maps_[v]->getGlobalNumElements();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(fullmap_->getGlobalNumElements() != numAllElements, std::logic_error,
                                 "logic error. full map and sub maps have not same number of elements. This cannot be. Please report the bug to the Xpetra developers!");

      // build importers for sub maps
      importers_.resize(maps_.size());
      for (unsigned i = 0; i < maps_.size(); ++i)
        if (maps_[i] != null)
          importers_[i] = Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(fullmap_, maps_[i]);
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, std::logic_error,
                                 "logic error. full map and sub maps are inconsistently distributed over the processors.");
    }

    /*! @brief Helper function to concatenate several maps

      @param  subMaps    vector of maps which are concatenated
      @return            concatenated map

      The routine builds a global map by concatenating all provided maps in the ordering defined by the vector.
      The GIDs are just appended in the same ordering as in the subMaps. No reordering or sorting is performed.
      This routine is supposed to generate the full map in an Xpetra::MapExtractor for a block operator. Note, it
      should not be used for strided maps since the GIDs are not reordered.

      Example: subMap[0] = { 0, 1, 3, 4 };
               subMap[1] = { 2, 5 };
               concatenated map = { 0, 1, 3, 4, 2 ,5 };
      */
    static Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > concatenateMaps(const std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > & subMaps) {

      // merge submaps to global map
      std::vector<GlobalOrdinal> gids;
      for(size_t tt = 0; tt<subMaps.size(); ++tt) {
        Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > subMap = subMaps[tt];
#if 1
        Teuchos::ArrayView< const GlobalOrdinal > subMapGids = subMap->getNodeElementList();
        gids.insert(gids.end(), subMapGids.begin(), subMapGids.end());
#else
        size_t myNumElements = subMap->getNodeNumElements();
        for(LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(myNumElements); ++l) {
          GlobalOrdinal gid = subMap->getGlobalElement(l);
          gids.push_back(gid);
        }
#endif
      }

      const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
      //std::sort(gids.begin(), gids.end());
      //gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
      Teuchos::ArrayView<GlobalOrdinal> gidsView(&gids[0], gids.size());
      Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > fullMap = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(subMaps[0]->lib(), INVALID, gidsView, subMaps[0]->getIndexBase(), subMaps[0]->getComm());
      return fullMap;
    }

  private:
    bool CheckConsistency() const {
      const RCP<const Map> fullMap = getFullMap();

      for (size_t i = 0; i < getNumMaps(); i++) {
        const RCP<const Map> map = getMap(i);

        ArrayView<const GlobalOrdinal> mapGids = map->getNodeElementList();
        for (typename ArrayView< const GlobalOrdinal >::const_iterator it = mapGids.begin(); it != mapGids.end(); it++)
          if (fullMap->isNodeGlobalElement(*it) == false)
            return false; // Global ID (*it) not found locally on this proc in fullMap -> error
      }
      return true;
    }

  private:
    RCP<const Map >               fullmap_;
    std::vector<RCP<const Map> >  maps_;
    std::vector<RCP<Import > >    importers_;
    bool                          bThyraMode_; //< boolean flag: use Thyra numbering for local sub-block maps. default = false (for Xpetra mode)
    std::vector<RCP<const Map> >  thyraMaps_;  //< store Thyra-style numbering maps here in Thyra mode. In Xpetra mode this vector is empty.
  }; // BlockedMap class

} // Xpetra namespace

#define XPETRA_BLOCKEDMAP_SHORT

#endif /* PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_HPP_ */
