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

#include "Xpetra_Import.hpp"
#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"


namespace Xpetra {

  template <class Scalar = MultiVector<>::scalar_type,
            class LocalOrdinal = Map<>::local_ordinal_type,
            class GlobalOrdinal = typename Map<LocalOrdinal>::global_ordinal_type,
            class Node = typename Map<LocalOrdinal, GlobalOrdinal>::node_type>
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
    MapExtractor(const RCP<const Map>& fullmap, const std::vector<RCP<const Map> >& maps) {
      fullmap_ = fullmap;
      maps_ = maps;

      importers_.resize(maps_.size());
      for (unsigned i = 0; i < maps_.size(); ++i)
        if (maps[i] != null)
          importers_[i] = ImportFactory::Build(fullmap_, maps[i]);

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

    RCP<     Vector> ExtractVector(RCP<const      Vector>& full, size_t block) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      const RCP<Vector> ret = VectorFactory::Build(getMap(block), true);
      ExtractVector(*full, block, *ret);
      return ret;
    }
    RCP<     Vector> ExtractVector(RCP<           Vector>& full, size_t block) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      const RCP<Vector> ret = VectorFactory::Build(getMap(block), true);
      ExtractVector(*full, block, *ret);
      return ret;
    }
    RCP<MultiVector> ExtractVector(RCP<const MultiVector>& full, size_t block) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      const RCP<Vector> ret = VectorFactory::Build(getMap(block), true);
      ExtractVector(*full, block, *ret);
      return ret;
    }
    RCP<MultiVector> ExtractVector(RCP<      MultiVector>& full, size_t block) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "] is null");
      const RCP<Vector> ret = VectorFactory::Build(getMap(block), true);
      ExtractVector(*full, block, *ret);
      return ret;
    }
    //@}

    /** \name Insert subblocks into full map */
    //@{
    void InsertVector(const Vector& partial, size_t block, Vector& full) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: maps_[" << block << "] is null");

      full.doExport(partial, *importers_[block], Xpetra::INSERT);
    }
    void InsertVector(const MultiVector& partial, size_t block, MultiVector& full) const {
      TEUCHOS_TEST_FOR_EXCEPTION(maps_[block] == null, Xpetra::Exceptions::RuntimeError,
            "InsertVector: maps_[" << block << "] is null");

      full.doExport(partial, *importers_[block], Xpetra::INSERT);
    }

    void InsertVector(RCP<const      Vector>& partial, size_t block, RCP<     Vector>& full) const { InsertVector(*partial, block, *full); }
    void InsertVector(RCP<           Vector>& partial, size_t block, RCP<     Vector>& full) const { InsertVector(*partial, block, *full); }
    void InsertVector(RCP<const MultiVector>& partial, size_t block, RCP<MultiVector>& full) const { InsertVector(*partial, block, *full); }
    void InsertVector(RCP<      MultiVector>& partial, size_t block, RCP<MultiVector>& full) const { InsertVector(*partial, block, *full); }


    //@}

    RCP<     Vector> getVector(size_t i) const                { return      VectorFactory::Build(getMap(i), true); }
    RCP<MultiVector> getVector(size_t i, size_t numvec) const { return MultiVectorFactory::Build(getMap(i), numvec, true); }

    /** \name Maps */
    //@{

    /// number of partial maps
    size_t NumMaps() const { return maps_.size(); }

    /// get the map
    const RCP<const Map> getMap(size_t i) const { return maps_[i]; }

    /// the full map
    const RCP<const Map> getFullMap() const { return fullmap_; }

    /// returns map index in map extractor which contains GID or -1 otherwise
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
    std::vector<RCP<const Map > > maps_;
    std::vector<RCP<Import > >    importers_;
  };
}

#define XPETRA_MAPEXTRACTOR_SHORT
#endif /* XPETRA_MAPEXTRACTOR_HPP_ */
