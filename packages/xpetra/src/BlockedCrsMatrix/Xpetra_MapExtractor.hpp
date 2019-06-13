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
      map_ = Teuchos::rcp(new BlockedMap(fullmap, maps, bThyraMode));
    }

    //! Expert constructor for Thyra maps
    MapExtractor(const std::vector<RCP<const Map> >& maps, const std::vector<RCP<const Map> >& thyramaps) {
      map_ = Teuchos::rcp(new BlockedMap(maps, thyramaps));
    }

    /*!
     * Constructor which accepts a const version
     * of a blocked map
     *
     * \param map BlockedMap defining the block structure of the multi vector

     */
    MapExtractor(const Teuchos::RCP< const BlockedMap > &map) :
          map_(map) { }

    //! copy constructor
    MapExtractor(const MapExtractor& input) {
      map_ = Teuchos::rcp(new BlockedMap(*(input.getBlockedMap())));
    }

    //! Destructor.
    virtual ~MapExtractor() {
      map_ = Teuchos::null;
    }

    /** \name Extract subblocks from full map */
    //@{
    void ExtractVector(const Vector& full, size_t block, Vector& partial) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getMap(" << block << ",false) is null");

      partial.doImport(full, *(map_->getImporter(block)), Xpetra::INSERT);
    }
    void ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getMap(" << block << ",false) is null");

      partial.doImport(full, *(map_->getImporter(block)), Xpetra::INSERT);
    }
    void ExtractVector(RCP<const      Vector>& full, size_t block, RCP<     Vector>& partial) const { ExtractVector(*full, block, *partial); }
    void ExtractVector(RCP<           Vector>& full, size_t block, RCP<     Vector>& partial) const { ExtractVector(*full, block, *partial); }
    void ExtractVector(RCP<const MultiVector>& full, size_t block, RCP<MultiVector>& partial) const { ExtractVector(*full, block, *partial); }
    void ExtractVector(RCP<      MultiVector>& full, size_t block, RCP<MultiVector>& partial) const { ExtractVector(*full, block, *partial); }

    RCP<     Vector> ExtractVector(RCP<const      Vector>& full, size_t block, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getMap(" << block << ",false) is null");
      // first extract partial vector from full vector (using xpetra style GIDs)
      const RCP<Vector> vv = VectorFactory::Build(getMap(block,false), false);
      ExtractVector(*full, block, *vv);
      if(bThyraMode == false) return vv;
      TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      vv->replaceMap(getMap(block,true)); // switch to Thyra-style map
      return vv;
    }

    RCP<     Vector> ExtractVector(RCP<           Vector>& full, size_t block, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getmap(" << block << ",false) is null");
      // first extract partial vector from full vector (using xpetra style GIDs)
      const RCP<Vector> vv = VectorFactory::Build(getMap(block,false), false);
      ExtractVector(*full, block, *vv);
      if(bThyraMode == false) return vv;
      TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      vv->replaceMap(getMap(block,true)); // switch to Thyra-style map
      return vv;
    }

    RCP<MultiVector> ExtractVector(RCP<const MultiVector>& full, size_t block, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getmap(" << block << ",false) is null");
      RCP<const BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(full);
      if(bfull.is_null() == true) {
        // standard case: full is not of type BlockedMultiVector
        // first extract partial vector from full vector (using xpetra style GIDs)
        const RCP<MultiVector> vv = MultiVectorFactory::Build(getMap(block,false), full->getNumVectors(), false);
        //if(bThyraMode == false) {
        //  ExtractVector(*full, block, *vv);
        //  return vv;
        //} else {
          RCP<const Map> oldThyMapFull = full->getMap(); // temporarely store map of full
          RCP<MultiVector> rcpNonConstFull = Teuchos::rcp_const_cast<MultiVector>(full);
          rcpNonConstFull->replaceMap(map_->getImporter(block)->getSourceMap());
          ExtractVector(*rcpNonConstFull, block, *vv);
          TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                     "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
          if(bThyraMode == true)
            vv->replaceMap(getMap(block,true)); // switch to Thyra-style map
          rcpNonConstFull->replaceMap(oldThyMapFull);
          return vv;
        //}
      } else {
        // special case: full is of type BlockedMultiVector
        XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(), Xpetra::Exceptions::RuntimeError,
              "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be " << bfull->getBlockedMap()->getNumMaps() << " (number of blocks in BlockedMultiVector)");
        return bfull->getMultiVector(block,bThyraMode);
      }
    }

    RCP<MultiVector> ExtractVector(RCP<      MultiVector>& full, size_t block, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getmap(" << block << ",false) is null");
      RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
      if(bfull.is_null() == true) {
        // standard case: full is not of type BlockedMultiVector
        // first extract partial vector from full vector (using xpetra style GIDs)
        const RCP<MultiVector> vv = MultiVectorFactory::Build(getMap(block,false), full->getNumVectors(), false);
        //if(bThyraMode == false) {
        //  ExtractVector(*full, block, *vv);
        //  return vv;
        //} else {
          RCP<const Map> oldThyMapFull = full->getMap(); // temporarely store map of full
          full->replaceMap(map_->getImporter(block)->getSourceMap());
          ExtractVector(*full, block, *vv);
          TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                     "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
          if(bThyraMode == true)
            vv->replaceMap(getMap(block,true)); // switch to Thyra-style map
          full->replaceMap(oldThyMapFull);
          return vv;
        //}
      } else {
        // special case: full is of type BlockedMultiVector
        XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(), Xpetra::Exceptions::RuntimeError,
              "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be " << bfull->getBlockedMap()->getNumMaps() << " (number of blocks in BlockedMultiVector)");
        return bfull->getMultiVector(block,bThyraMode);
      }
    }

    RCP<MultiVector> ExtractVector(RCP<const BlockedMultiVector>& full, size_t block, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getmap(" << block << ",false) is null");
      XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != full->getBlockedMap()->getNumMaps(), Xpetra::Exceptions::RuntimeError,
            "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be " << full->getBlockedMap()->getNumMaps() << " (number of blocks in BlockedMultiVector)");
      Teuchos::RCP<MultiVector> vv = full->getMultiVector(block,bThyraMode);
      return vv;
    }

    RCP<MultiVector> ExtractVector(RCP<      Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& full, size_t block, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getmap(" << block << ",false) is null");
      XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != full->getBlockedMap()->getNumMaps(), Xpetra::Exceptions::RuntimeError,
            "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be " << full->getBlockedMap()->getNumMaps() << " (number of blocks in BlockedMultiVector)");
      Teuchos::RCP<MultiVector> vv = full->getMultiVector(block,bThyraMode);
      return vv;
    }
    //@}

    /** \name Insert subblocks into full map */
    //@{
    void InsertVector(const Vector& partial, size_t block, Vector& full, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getmap(" << block << ",false) is null");
      XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      if(bThyraMode) {
        // NOTE: the importer objects in the BlockedMap are always using Xpetra GIDs (or Thyra style Xpetra GIDs)
        // The source map corresponds to the full map (in Xpetra GIDs) starting with GIDs from zero. The GIDs are consecutive in Thyra mode
        // The target map is the partial map (in the corresponding Xpetra GIDs)

        // TODO can we skip the Export call in special cases (i.e. Src = Target map, same length, etc...)

        // store original GIDs (could be Thyra GIDs)
        RCP<const MultiVector> rcpPartial = Teuchos::rcpFromRef(partial);
        RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
        RCP<const Map> oldThyMapPartial = rcpNonConstPartial->getMap(); // temporarely store map of partial
        RCP<const Map> oldThyMapFull = full.getMap(); // temporarely store map of full

        // check whether getMap(block,false) is identical to target map of importer
        XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false)->isSameAs(*(map_->getImporter(block)->getTargetMap()))==false, Xpetra::Exceptions::RuntimeError,
                   "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of partial vector are not identical to target Map of Importer. This should not be.");

        //XPETRA_TEST_FOR_EXCEPTION(full.getMap()->isSameAs(*(map_->getImporter(block)->getSourceMap()))==false, Xpetra::Exceptions::RuntimeError,
        //           "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of full vector are not identical to source Map of Importer. This should not be.");

        rcpNonConstPartial->replaceMap(getMap(block,false)); // temporarely switch to xpetra-style map
        full.replaceMap(map_->getImporter(block)->getSourceMap());     // temporarely switch to Xpetra GIDs

        // do the Export
        full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Xpetra::INSERT);

        // switch back to original maps
        full.replaceMap(oldThyMapFull); // reset original map (Thyra GIDs)
        rcpNonConstPartial->replaceMap(oldThyMapPartial); // change map back to original map
      } else {
        // Xpetra style numbering
        full.doExport(partial, *(map_->getImporter(block)), Xpetra::INSERT);
      }
    }
    void InsertVector(const MultiVector& partial, size_t block, MultiVector& full, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range, "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: map_->getmap(" << block << ",false) is null");
      XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      if(bThyraMode) {
        // NOTE: the importer objects in the BlockedMap are always using Xpetra GIDs (or Thyra style Xpetra GIDs)
        // The source map corresponds to the full map (in Xpetra GIDs) starting with GIDs from zero. The GIDs are consecutive in Thyra mode
        // The target map is the partial map (in the corresponding Xpetra GIDs)

        // TODO can we skip the Export call in special cases (i.e. Src = Target map, same length, etc...)

        // store original GIDs (could be Thyra GIDs)
        RCP<const MultiVector> rcpPartial = Teuchos::rcpFromRef(partial);
        RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
        RCP<const Map> oldThyMapPartial = rcpNonConstPartial->getMap(); // temporarely store map of partial
        RCP<const Map> oldThyMapFull = full.getMap(); // temporarely store map of full

        // check whether getMap(block,false) is identical to target map of importer
        XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false)->isSameAs(*(map_->getImporter(block)->getTargetMap()))==false, Xpetra::Exceptions::RuntimeError,
                   "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of partial vector are not identical to target Map of Importer. This should not be.");

        //XPETRA_TEST_FOR_EXCEPTION(full.getMap()->isSameAs(*(map_->getImporter(block)->getSourceMap()))==false, Xpetra::Exceptions::RuntimeError,
        //           "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of full vector are not identical to source Map of Importer. This should not be.");

        rcpNonConstPartial->replaceMap(getMap(block,false)); // temporarely switch to xpetra-style map
        full.replaceMap(map_->getImporter(block)->getSourceMap());     // temporarely switch to Xpetra GIDs

        // do the Export
        full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Xpetra::INSERT);

        // switch back to original maps
        full.replaceMap(oldThyMapFull); // reset original map (Thyra GIDs)
        rcpNonConstPartial->replaceMap(oldThyMapPartial); // change map back to original map
      } else {
        // Xpetra style numbering
        full.doExport(partial, *(map_->getImporter(block)), Xpetra::INSERT);
      }
    }

    void InsertVector(RCP<const      Vector> partial, size_t block, RCP<     Vector> full, bool bThyraMode = false) const { InsertVector(*partial, block, *full, bThyraMode); }
    void InsertVector(RCP<           Vector> partial, size_t block, RCP<     Vector> full, bool bThyraMode = false) const { InsertVector(*partial, block, *full, bThyraMode); }
    void InsertVector(RCP<const MultiVector> partial, size_t block, RCP<MultiVector> full, bool bThyraMode = false) const {
      RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
      if(bfull.is_null() == true)
        InsertVector(*partial, block, *full, bThyraMode);
      else {
        XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
              "InsertVector: map_->getmap(" << block << ",false) is null");
        full->setMultiVector(block, partial, bThyraMode);
      }
    }
    void InsertVector(RCP<      MultiVector> partial, size_t block, RCP<MultiVector> full, bool bThyraMode = false) const {
      RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
      if(bfull.is_null() == true)
        InsertVector(*partial, block, *full, bThyraMode);
      else {
        XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
              "InsertVector: map_->getmap(" << block << ",false) is null");
        bfull->setMultiVector(block, partial, bThyraMode);
      }
    }
    void InsertVector(RCP<const MultiVector> partial, size_t block, RCP<Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > full, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: map_->getmap(" << block << ",false) is null");
      full->setMultiVector(block, partial, bThyraMode);
    }
    void InsertVector(RCP<      MultiVector> partial, size_t block, RCP<Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > full, bool bThyraMode = false) const {
      XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false) == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: map_->getmap(" << block << ",false) is null");
      full->setMultiVector(block, partial, bThyraMode);
    }

    //@}

    RCP<     Vector> getVector(size_t i, bool bThyraMode = false, bool bZero = true) const {
      XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      // TODO check whether this can return a blocked multivector
      return VectorFactory::Build(getMap(i, bThyraMode), bZero);
    }
    RCP<MultiVector> getVector(size_t i, size_t numvec, bool bThyraMode = false, bool bZero = true) const {
      XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, Xpetra::Exceptions::RuntimeError,
                 "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using Thyra-style numbered submaps.");
      // TODO check whether this can return a blocked multivector
      return MultiVectorFactory::Build(getMap(i, bThyraMode), numvec, bZero);
    }

    /// returns true, if sub maps are stored in Thyra-style numbering
    bool getThyraMode() const { return map_->getThyraMode(); }

    /** \name Maps */
    //@{

    /// number of partial maps
    size_t NumMaps() const { return map_->getNumMaps(); }

    /// get the map
    /// returns the sub map i from list of sub maps
    /// depending on the parameter bThyraMode the sub map that is returned uses Thyra or Xpetra numbering
    /// Note: Thyra-numbering is only allowed if the MapExtractor is also constructed using Thyra numbering
    const RCP<const Map> getMap(size_t i, bool bThyraMode = false) const {
      return map_->getMap(i,bThyraMode);
    }

    /// get the underlying BlockedMap object (as Map)
    const RCP<const Map> getMap() const { return map_; }

    /// get the underlying BlockedMap object (as BlockedMap)
    const RCP<const BlockedMap> getBlockedMap() const { return map_; }

    /// the full map
    const RCP<const Map> getFullMap() const { return map_->getFullMap(); }

    /// returns map index in map extractor which contains GID
    size_t getMapIndexForGID(GlobalOrdinal gid) const {
      return map_->getMapIndexForGID(gid);
    }

    //@}

  private:
    Teuchos::RCP<const BlockedMap> map_;         ///< blocked map containing the sub block maps (either thyra or xpetra mode)
  };
}

#define XPETRA_MAPEXTRACTOR_SHORT
#endif /* XPETRA_MAPEXTRACTOR_HPP_ */
