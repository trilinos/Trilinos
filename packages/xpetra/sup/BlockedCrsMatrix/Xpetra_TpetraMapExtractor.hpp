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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * Xpetra_TpetraMapExtractor.hpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_TPETRAMAPEXTRACTOR_HPP_
#define XPETRA_TPETRAMAPEXTRACTOR_HPP_

#include <Xpetra_MapExtractor.hpp>

namespace Xpetra
{
  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class TpetraMapExtractor : public Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  {
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> MapClass;
    typedef TpetraMap<LocalOrdinal,GlobalOrdinal,Node> TpetraMapClass;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> VectorClass;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVectorClass;
    typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraVectorClass;
    typedef TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraMultiVectorClass;
    typedef Import<LocalOrdinal,GlobalOrdinal,Node> ImportClass;
    typedef TpetraImport<LocalOrdinal,GlobalOrdinal,Node> TpetraImportClass;

  public:
    //! TpetraMapExtractor basic constructor
    TpetraMapExtractor(const Teuchos::RCP<const MapClass>& fullmap, const std::vector<Teuchos::RCP<const MapClass> >& maps)
    {
      XPETRA_RCP_DYNAMIC_CAST(const TpetraMapClass, fullmap, tfullmap, "Xpetra::TpetraMapextractor constructors only accept Xpetra::TpetraMap as input arguments.");

      fullmap_ = tfullmap;

      unsigned int numMaps = maps.size();
      maps_.empty();
      importer_.resize(numMaps);
      for (unsigned i=0; i<numMaps; ++i)
      {
        if (maps[i]!=Teuchos::null)
        {
          XPETRA_RCP_DYNAMIC_CAST(const TpetraMapClass, maps[i], tmapi, "Xpetra::TpetraMapextractor constructors only accept Xpetra::TpetraMap as input arguments.");
          maps_.push_back(tmapi);
          importer_[i] = Teuchos::rcp(new Xpetra::TpetraImport<LocalOrdinal,GlobalOrdinal,Node>(tfullmap,tmapi));
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency()==false,std::logic_error,"logic error. full map and sub maps are inconsistently distributed over the processors.");
    }

    //! TpetraMapExtractor constructor
    TpetraMapExtractor(const Teuchos::RCP<const TpetraMapClass>& fullmap, const std::vector<Teuchos::RCP<const TpetraMapClass> >& maps)
    {
      XPETRA_RCP_DYNAMIC_CAST(const TpetraMapClass, fullmap, tfullmap, "Xpetra::TpetraMapextractor constructors only accept Xpetra::TpetraMap as input arguments.");

      fullmap_ = tfullmap;
      maps_    = maps;

      importer_.resize(maps_.size());
      for (unsigned i=0; i< importer_.size(); ++i)
      {
        if (maps_[i]!=Teuchos::null)
        {
          Teuchos::RCP<TpetraImportClass> tImport = Teuchos::rcp_dynamic_cast<TpetraImportClass>(Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(fullmap_,maps_[i]));
          if (tImport == Teuchos::null) std::cout << "no TpetraImportFactory?" << std::endl;
          importer_[i] = tImport;
          std::cout << *importer_[i] << std::endl;
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency()==false,std::logic_error,"logic error. full map and sub maps are inconsistently distributed over the processors.");
    }

    //! Destructor
    inline virtual ~TpetraMapExtractor() {};

    void InsertVector(const TpetraVectorClass& partial, size_t block, TpetraVectorClass& full) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block" << block << std::endl;

      full.doExport(partial,*importer_[block],Xpetra::INSERT);
    }

    void InsertVector(Teuchos::RCP<const VectorClass>& partial, size_t block, Teuchos::RCP<VectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(const TpetraVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");

      InsertVector(*tpartial,block,*tfull);
    }

    void InsertVector(Teuchos::RCP<VectorClass>& partial, size_t block, Teuchos::RCP<VectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accept Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(TpetraVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accept Xpetra::TpetraMap as input arguments.");

      InsertVector(*tpartial,block,*tfull);
    }

    void InsertVector(const TpetraMultiVectorClass& partial, size_t block, TpetraMultiVectorClass& full) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block" << block << std::endl;

      full.doExport(partial,*importer_[block],Xpetra::INSERT);
    }

    void InsertVector(Teuchos::RCP<const MultiVectorClass>& partial, size_t block, Teuchos::RCP<MultiVectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraMultiVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(const TpetraMultiVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");

      InsertVector(*tpartial,block,*tfull);
    }

    void InsertVector(Teuchos::RCP<MultiVectorClass>& partial, size_t block, Teuchos::RCP<MultiVectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraMultiVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accept Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(TpetraMultiVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accept Xpetra::TpetraMap as input arguments.");

      InsertVector(*tpartial,block,*tfull);
    }

    void ExtractVector(const TpetraVectorClass& full, size_t block, TpetraVectorClass& partial) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block " << block << std::endl;
      partial.doImport(full, *importer_[block], Xpetra::INSERT);
    }

    void ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const TpetraVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(      TpetraVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");

      ExtractVector(*tfull,block,*tpartial);
    }

    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const TpetraVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<TpetraVectorClass> ret = Teuchos::rcp(new TpetraVectorClass(getMap(block),true));
      ExtractVector(*tfull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<VectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<TpetraVectorClass> ret = Teuchos::rcp(new TpetraVectorClass(getMap(block),true));
      ExtractVector(*tfull,block,*ret);
      return ret;
    }

    void ExtractVector(Teuchos::RCP<VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(TpetraVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");

      ExtractVector(*tfull,block,*tpartial);
    }

    void ExtractVector(const TpetraMultiVectorClass& full, size_t block, TpetraMultiVectorClass& partial) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block " << block << std::endl;
      partial.doImport(full, *importer_[block], Xpetra::INSERT);
    }

    virtual void ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const TpetraMultiVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(      TpetraMultiVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");

      ExtractVector(*tfull,block,*tpartial);
    }

    virtual void ExtractVector(Teuchos::RCP<      MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(      TpetraMultiVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(      TpetraMultiVectorClass, partial, tpartial, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");

      ExtractVector(*tfull,block,*tpartial);
    }

    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const TpetraMultiVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<TpetraMultiVectorClass> ret = Teuchos::rcp(new TpetraMultiVectorClass(getMap(block),full->getNumVectors(),true));
      ExtractVector(*tfull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<MultiVectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraMultiVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<TpetraMultiVectorClass> ret = Teuchos::rcp(new TpetraMultiVectorClass(getMap(block),full->getNumVectors(),true));
      ExtractVector(*tfull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<VectorClass> getVector(size_t i) const
    {
      return Teuchos::rcp(new TpetraVectorClass(getMap(i),true));
    }

    virtual Teuchos::RCP<MultiVectorClass> getVector(size_t i, size_t numvec) const
    {
      return Teuchos::rcp(new TpetraMultiVectorClass(getMap(i),numvec,true));
    }


    /** \name Maps */
    //@{

    /// number of partial maps
    size_t NumMaps() const { return maps_.size(); }

    /// get the map
    const Teuchos::RCP<const MapClass> getMap(size_t i) const { return maps_[i]; }

    /// the full map
    const Teuchos::RCP<const MapClass> getFullMap() const { return fullmap_; }

    size_t getMapIndexForGID(GlobalOrdinal gid) const {

      for(size_t i = 0; i < NumMaps(); i++) {
        if(getMap(i)->isNodeGlobalElement(gid) == true)
          return i;
      }
      TEUCHOS_TEST_FOR_EXCEPTION( false, Xpetra::Exceptions::RuntimeError,
                  "getMapIndexForGID: GID " << gid << " is not contained by a map in mapextractor." );
      return 0;
    }

    //@}
  private:

    virtual bool CheckConsistency() const {
      const Teuchos::RCP<const MapClass> fullMap = getFullMap();

      for(size_t i = 0; i<NumMaps(); i++) {
        const Teuchos::RCP<const MapClass> map = getMap(i);

        Teuchos::ArrayView< const GlobalOrdinal > mapGids = map->getNodeElementList();
        typename Teuchos::ArrayView< const GlobalOrdinal >::const_iterator it;
        for(it = mapGids.begin(); it!=mapGids.end(); it++) {
          if(fullMap->isNodeGlobalElement(*it)==false)
            return false; // Global ID (*it) not found locally on this proc in fullMap -> error
        }
      }
      return true;
    }
  protected:
    std::vector<Teuchos::RCP<const TpetraMapClass > > maps_;
    Teuchos::RCP<const TpetraMapClass > fullmap_;
    std::vector<Teuchos::RCP<TpetraImportClass > > importer_;
  };
}


#endif /* XPETRA_TPETRAMAPEXTRACTOR_HPP_ */
