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
#ifndef XPETRA_MAPEXTRACTOR_HPP_
#define XPETRA_MAPEXTRACTOR_HPP_

#include <map>

#include <iostream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapUtils.hpp>
#include <Xpetra_MultiVector.hpp>
//#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of BlockedMultiVector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class BlockedMultiVector;
#endif

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class MapExtractor : public Teuchos::Describable {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

  private:
#undef XPETRA_MAPEXTRACTOR_SHORT
#include "Xpetra_UseShortNames.hpp"

  public:

    //! MapExtractor basic constructor
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
    MapExtractor(const RCP<const Map>& fullmap, const std::vector<RCP<const Map> >& maps, bool bThyraMode = false) {
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
                                   "logic error. full map and sub maps have not same number of elements. We cannot build MapExtractor with Xpetra-style numbering. Please make sure that you want Xpetra-style numbering instead of Thyra-style numbering.");

        fullmap_ = fullmap;
        maps_ = maps;
      } else {
        //std::cout << "Create Map Extractor in Thyra Mode!!! " << std::endl;
        // use Thyra-style numbering for sub-block maps
        // That is, all sub-block maps start with zero as GID and are contiguous

        // plausibility check
        for(size_t v = 0; v < maps.size(); ++v) {
          TEUCHOS_TEST_FOR_EXCEPTION(maps[v]->getMinAllGlobalIndex() != 0, std::logic_error,
                                             "logic error. When using Thyra-style numbering all sub-block maps must start with zero as GID.");
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
          std::vector<GlobalOrdinal> subMapGids(maps[v]->getNodeNumElements(),0);
          for (LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(maps[v]->getNodeNumElements()); ++l) {
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
        fullmap_ = MapFactory::Build(fullmap->lib(), INVALID, fullMapGidsView, fullmap->getIndexBase(), fullmap->getComm());

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
    MapExtractor(const std::vector<RCP<const Map> >& maps, const std::vector<RCP<const Map> >& thyramaps) {
      bThyraMode_ = true;


      // plausibility check
      TEUCHOS_TEST_FOR_EXCEPTION(thyramaps.size() != maps.size(), std::logic_error, "logic error. The number of submaps must be identical!");
      for(size_t v = 0; v < thyramaps.size(); ++v) {
        TEUCHOS_TEST_FOR_EXCEPTION(thyramaps[v]->getMinAllGlobalIndex() != 0, std::logic_error,
                                           "logic error. When using Thyra-style numbering all sub-block maps must start with zero as GID.");
        TEUCHOS_TEST_FOR_EXCEPTION(thyramaps[v]->getNodeNumElements() != maps[v]->getNodeNumElements(), std::logic_error,
                                           "logic error. The size of the submaps must be identical (same distribution, just different GIDs)");
      }

      // store user-provided maps and thyramaps
      thyraMaps_ = thyramaps;
      maps_      = maps;

      fullmap_ = Xpetra::MapUtils<LocalOrdinal,GlobalOrdinal,Node>::concatenateMaps(maps);

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
    MapExtractor(const MapExtractor& input) {
      bThyraMode_ = input.getThyraMode();

      fullmap_ = MapFactory::Build(input.getFullMap(),1);

      maps_.resize(input.NumMaps(), Teuchos::null);
      if(bThyraMode_ == true)
        thyraMaps_.resize(input.NumMaps(), Teuchos::null);
      for(size_t i = 0; i < input.NumMaps(); ++i) {
        maps_[i] = MapFactory::Build(input.getMap(i,false),1);
        if(bThyraMode_ == true)
          thyraMaps_[i] = MapFactory::Build(input.getMap(i,true),1);
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

    /** \name Extract subblocks from full map */
    //@{
    void ExtractVector(const Vector& full, size_t block, Vector& partial) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");

      partial.doImport(full, *importers_[block], Xpetra::INSERT);
    }
    void ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      partial.doImport(full, *importers_[block], Xpetra::INSERT);
    }
    void ExtractVector(RCP<const      Vector>& full, size_t block, RCP<     Vector>& partial) const { ExtractVector(*full, block, *partial); }
    void ExtractVector(RCP<           Vector>& full, size_t block, RCP<     Vector>& partial) const { ExtractVector(*full, block, *partial); }
    void ExtractVector(RCP<const MultiVector>& full, size_t block, RCP<MultiVector>& partial) const { ExtractVector(*full, block, *partial); }
    void ExtractVector(RCP<      MultiVector>& full, size_t block, RCP<MultiVector>& partial) const { ExtractVector(*full, block, *partial); }

    RCP<     Vector> ExtractVector(RCP<const      Vector>& full, size_t block, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      // extract sub vector
      const RCP<Vector> xpetraVec = VectorFactory::Build(getMap(block,false), true);
      ExtractVector(*full, block, *xpetraVec);
      if(bThyraMode == false) return xpetraVec;
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      const RCP<Vector> thyraVec = VectorFactory::Build(getMap(block,true), true);
      // TODO introduce Kokkos version of this.
      Teuchos::ArrayRCP<Scalar> thyraVecData  = thyraVec->getDataNonConst(0);
      Teuchos::ArrayRCP<const Scalar> xpetraVecData = xpetraVec->getData(0);

      for(size_t i=0; i < xpetraVec->getLocalLength(); i++) {
        thyraVecData[i] = xpetraVecData[i];
      }
      return thyraVec;
    }
    RCP<     Vector> ExtractVector(RCP<           Vector>& full, size_t block, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      const RCP<Vector> xpetraVec = VectorFactory::Build(getMap(block,false), true);
      ExtractVector(*full, block, *xpetraVec);
      if(bThyraMode == false) return xpetraVec;
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      const RCP<Vector> thyraVec = VectorFactory::Build(getMap(block,true), true);
      // TODO introduce Kokkos version of this.
      Teuchos::ArrayRCP<Scalar> thyraVecData  = thyraVec->getDataNonConst(0);
      Teuchos::ArrayRCP<const Scalar> xpetraVecData = xpetraVec->getData(0);

      for(size_t i=0; i < xpetraVec->getLocalLength(); i++) {
        thyraVecData[i] = xpetraVecData[i];
      }
      return thyraVec;
    }
    RCP<MultiVector> ExtractVector(RCP<const MultiVector>& full, size_t block, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      const RCP<MultiVector> xpetraVec = MultiVectorFactory::Build(getMap(block,false), full->getNumVectors(), true);
      ExtractVector(*full, block, *xpetraVec);
      if(bThyraMode == false) return xpetraVec;
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      const RCP<MultiVector> thyraVec = MultiVectorFactory::Build(getMap(block,true), xpetraVec->getNumVectors(), true);
      // TODO introduce Kokkos version of this.
      for(size_t k=0; k < xpetraVec->getNumVectors(); k++) {
        Teuchos::ArrayRCP<Scalar> thyraVecData  = thyraVec->getDataNonConst(k);
        Teuchos::ArrayRCP<const Scalar> xpetraVecData = xpetraVec->getData(k);
        for(size_t i=0; i < xpetraVec->getLocalLength(); i++) {
          thyraVecData[i] = xpetraVecData[i];
        }
      }
      return thyraVec;
    }
    RCP<MultiVector> ExtractVector(RCP<      MultiVector>& full, size_t block, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      const RCP<MultiVector> xpetraVec = MultiVectorFactory::Build(getMap(block,false), full->getNumVectors(), true);
      ExtractVector(*full, block, *xpetraVec);
      if(bThyraMode == false) return xpetraVec;
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      const RCP<MultiVector> thyraVec = MultiVectorFactory::Build(getMap(block,true), xpetraVec->getNumVectors(), true);
      // TODO introduce Kokkos version of this.
      for(size_t k=0; k < xpetraVec->getNumVectors(); k++) {
        Teuchos::ArrayRCP<Scalar> thyraVecData  = thyraVec->getDataNonConst(k);
        Teuchos::ArrayRCP<const Scalar> xpetraVecData = xpetraVec->getData(k);
        for(size_t i=0; i < xpetraVec->getLocalLength(); i++) {
          thyraVecData[i] = xpetraVecData[i];
        }
      }
      return thyraVec;
    }
    RCP<MultiVector> ExtractVector(RCP<const BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& full, size_t block, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      TEUCHOS_TEST_FOR_EXCEPTION(NumMaps() != full->getMapExtractor()->NumMaps(), Xpetra::Exceptions::RuntimeError,
            "ExtractVector: Number of blocks in map extractor is " << NumMaps() << " but should be " << full->getMapExtractor()->NumMaps() << " (number of blocks in BlockedMultiVector)");
      Teuchos::RCP<MultiVector> vv = full->getMultiVector(block);
      TEUCHOS_TEST_FOR_EXCEPTION(vv->getMap()->isSameAs(*(getMap(block,bThyraMode_))) == false, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: partial map of BlockedMultiVector and MapExtractor are incompatible");
      return vv;
    }
    RCP<MultiVector> ExtractVector(RCP<      Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& full, size_t block, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      TEUCHOS_TEST_FOR_EXCEPTION(NumMaps() != full->getMapExtractor()->NumMaps(), Xpetra::Exceptions::RuntimeError,
            "ExtractVector: Number of blocks in map extractor is " << NumMaps() << " but should be " << full->getMapExtractor()->NumMaps() << " (number of blocks in BlockedMultiVector)");
      Teuchos::RCP<MultiVector> vv = full->getMultiVector(block);
      TEUCHOS_TEST_FOR_EXCEPTION(vv->getMap()->isSameAs(*(getMap(block,bThyraMode_))) == false, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: partial map of BlockedMultiVector and MapExtractor are incompatible");
      return vv;
    }
    //@}

    /** \name Insert subblocks into full map */
    //@{
    void InsertVector(const Vector& partial, size_t block, Vector& full, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: maps_[" << block << "] is null");
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      if(bThyraMode) {
        const RCP<Vector> xpetraVec = VectorFactory::Build(getMap(block,false), true); // get sub vector in xpetra-style numbering
        // TODO introduce Kokkos version of this.
        Teuchos::ArrayRCP<const Scalar> thyraVecData  = partial.getData(0);
        Teuchos::ArrayRCP<Scalar> xpetraVecData = xpetraVec->getDataNonConst(0);
        for(size_t i=0; i < xpetraVec->getLocalLength(); i++) {
          xpetraVecData[i] = thyraVecData[i];
        }
        full.doExport(*xpetraVec, *importers_[block], Xpetra::INSERT);
      } else {
        // Xpetra style numbering
        full.doExport(partial, *importers_[block], Xpetra::INSERT);
      }
    }
    void InsertVector(const MultiVector& partial, size_t block, MultiVector& full, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: maps_[" << block << "] is null");
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      if(bThyraMode) {
        const RCP<MultiVector> xpetraVec = MultiVectorFactory::Build(getMap(block,false), partial.getNumVectors(), true); // get sub vector in xpetra-style numbering
        // TODO introduce Kokkos version of this.
        for(size_t k=0; k < partial.getNumVectors(); k++) {
          Teuchos::ArrayRCP<const Scalar> thyraVecData  = partial.getData(k);
          Teuchos::ArrayRCP<Scalar> xpetraVecData = xpetraVec->getDataNonConst(k);
          for(size_t i=0; i < xpetraVec->getLocalLength(); i++) {
            xpetraVecData[i] = thyraVecData[i];
          }
        }
        full.doExport(*xpetraVec, *importers_[block], Xpetra::INSERT);
      } else {
        // Xpetra style numbering
        full.doExport(partial, *importers_[block], Xpetra::INSERT);
      }
    }

    void InsertVector(RCP<const      Vector> partial, size_t block, RCP<     Vector> full, bool bThyraMode = false) const { InsertVector(*partial, block, *full, bThyraMode); }
    void InsertVector(RCP<           Vector> partial, size_t block, RCP<     Vector> full, bool bThyraMode = false) const { InsertVector(*partial, block, *full, bThyraMode); }
    void InsertVector(RCP<const MultiVector> partial, size_t block, RCP<MultiVector> full, bool bThyraMode = false) const { InsertVector(*partial, block, *full, bThyraMode); }
    void InsertVector(RCP<      MultiVector> partial, size_t block, RCP<MultiVector> full, bool bThyraMode = false) const { InsertVector(*partial, block, *full, bThyraMode); }

    void InsertVector(RCP<const MultiVector> partial, size_t block, RCP<Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > full, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: maps_[" << block << "] is null");
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      full->setMultiVector(block, partial);
    }
    void InsertVector(RCP<      MultiVector> partial, size_t block, RCP<Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > full, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: maps_[" << block << "] is null");
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      full->setMultiVector(block, partial);
    }

    //@}

    RCP<     Vector> getVector(size_t i, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      return VectorFactory::Build(getMap(i, bThyraMode), true);
    }
    RCP<MultiVector> getVector(size_t i, size_t numvec, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      return MultiVectorFactory::Build(getMap(i, bThyraMode), numvec, true);
    }

    /// returns true, if sub maps are stored in Thyra-style numbering
    bool getThyraMode() const { return bThyraMode_; }

    /** \name Maps */
    //@{

    /// number of partial maps
    size_t NumMaps() const { return maps_.size(); }

    /// get the map
    /// returns the sub map i from list of sub maps
    /// depending on the parameter bThyraMode the sub map that is returned uses Thyra or Xpetra numbering
    /// Note: Thyra-numbering is only allowed if the MapExtractor is also constructed using Thyra numbering
    const RCP<const Map> getMap(size_t i, bool bThyraMode = false) const {
      TEUCHOS_TEST_FOR_EXCEPTION( i >= NumMaps(), Xpetra::Exceptions::RuntimeError, "MapExtractor::getMap: tried to access block " << i << ", but MapExtractor has only " << NumMaps() << " blocks! Block indices must be between 0 and " << NumMaps() - 1 << ".");
      TEUCHOS_TEST_FOR_EXCEPTION( i < 0, Xpetra::Exceptions::RuntimeError, "MapExtractor::getMap: A negative block index " << i << " is invalid. Block indices must be between 0 and " << NumMaps() - 1 << ".");
      if(bThyraMode_ == true && bThyraMode == true)
        return thyraMaps_[i];
      TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::getMap: cannot return sub map in Thyra-style numbering if MapExtractor object is not created using Thyra-style numbered submaps.");
      return maps_[i];
    }

    /// the full map
    const RCP<const Map> getFullMap() const { return fullmap_; }

    /// returns map index in map extractor which contains GID
    size_t getMapIndexForGID(GlobalOrdinal gid) const {
      for (size_t i = 0; i < NumMaps(); i++)
        if (getMap(i)->isNodeGlobalElement(gid) == true)
          return i;

      TEUCHOS_TEST_FOR_EXCEPTION(false, Xpetra::Exceptions::RuntimeError,
                                 "getMapIndexForGID: GID " << gid << " is not contained by a map in mapextractor." );
      return 0;
    }

    //@}

  private:
    bool CheckConsistency() const {
      const RCP<const Map> fullMap = getFullMap();

      for (size_t i = 0; i < NumMaps(); i++) {
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
  };
}

#define XPETRA_MAPEXTRACTOR_SHORT
#endif /* XPETRA_MAPEXTRACTOR_HPP_ */
