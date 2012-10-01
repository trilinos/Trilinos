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
 * Xpetra_EpetraMapExtractor.hpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_EPETRAMAPEXTRACTOR_HPP_
#define XPETRA_EPETRAMAPEXTRACTOR_HPP_

#include <Xpetra_MapExtractor.hpp>

namespace Xpetra
{
  class EpetraMapExtractor : public Xpetra::MapExtractor<double,int,int>
  {
    typedef Xpetra::Map<int,int> MapClass;
    typedef Xpetra::Vector<double,int,int> VectorClass;
    typedef Xpetra::MultiVector<double,int,int> MultiVectorClass;
    typedef Xpetra::Import<int,int> ImportClass;

  public:
    //! EpetraMapExtractor basic constructor
    EpetraMapExtractor(const Teuchos::RCP<const MapClass>& fullmap, const std::vector<Teuchos::RCP<const MapClass> >& maps)
    {
      XPETRA_RCP_DYNAMIC_CAST(const EpetraMap, fullmap, efullmap, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      fullmap_ = efullmap;

      unsigned int numMaps = maps.size();
      maps_.empty();
      importer_.resize(numMaps);
      for (unsigned i=0; i<numMaps; ++i)
      {
        if (maps[i]!=Teuchos::null)
        {
          XPETRA_RCP_DYNAMIC_CAST(const EpetraMap, maps[i], emapi, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
          maps_.push_back(emapi);
          importer_[i] = Xpetra::ImportFactory<int,int>::Build(efullmap,emapi);
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency()==false,std::logic_error,"logic error. full map and sub maps are inconsistently distributed over the processors.");
    }

    //! EpetraMapExtractor constructor
    EpetraMapExtractor(const Teuchos::RCP<const Xpetra::EpetraMap>& fullmap, const std::vector<Teuchos::RCP<const Xpetra::EpetraMap> >& maps)
    {
      XPETRA_RCP_DYNAMIC_CAST(const EpetraMap, fullmap, efullmap, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      fullmap_ = efullmap;
      maps_    = maps;

      importer_.resize(maps_.size());
      for (unsigned i=0; i< importer_.size(); ++i)
      {
        if (maps_[i]!=Teuchos::null)
        {
          importer_[i] = Xpetra::ImportFactory<int,int>::Build(fullmap,maps_[i]);
          std::cout << *importer_[i] << std::endl;
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency()==false,std::logic_error,"logic error. full map and sub maps are inconsistently distributed over the processors.");
    }

    //! Destructor
    virtual inline ~EpetraMapExtractor() {};

    void InsertVector(const Xpetra::EpetraVector& partial, size_t block, Xpetra::EpetraVector& full) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block" << block << std::endl;

      full.doExport(partial,*importer_[block],Xpetra::INSERT);
    }

    void InsertVector(Teuchos::RCP<const VectorClass>& partial, size_t block, Teuchos::RCP<VectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(EpetraVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(const EpetraVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      InsertVector(*epartial,block,*efull);
    }

    void InsertVector(Teuchos::RCP<VectorClass>& partial,size_t block, Teuchos::RCP<VectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(EpetraVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(EpetraVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      InsertVector(*epartial,block,*efull);
    }

    void InsertVector(const Xpetra::EpetraMultiVector& partial,size_t block, Xpetra::EpetraMultiVector& full) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block" << block << std::endl;

      full.doExport(partial,*importer_[block],Xpetra::INSERT);
    }

    void InsertVector(Teuchos::RCP<const MultiVectorClass>& partial,size_t block, Teuchos::RCP<MultiVectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(EpetraMultiVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(const EpetraMultiVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      InsertVector(*epartial,block,*efull);
    }

    void InsertVector(Teuchos::RCP<MultiVectorClass>& partial,size_t block, Teuchos::RCP<MultiVectorClass>& full) const
    {
      XPETRA_RCP_DYNAMIC_CAST(EpetraMultiVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(EpetraMultiVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      InsertVector(*epartial,block,*efull);
    }

    void ExtractVector(const Xpetra::EpetraVector& full, size_t block, Xpetra::EpetraVector& partial) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block " << block << std::endl;
      partial.doImport(full, *importer_[block], Xpetra::INSERT);
    }

    void ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const EpetraVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(      EpetraVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      ExtractVector(*efull,block,*epartial);
    }

    void ExtractVector(Teuchos::RCP<VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(EpetraVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(EpetraVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      ExtractVector(*efull,block,*epartial);
    }

    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const EpetraVector, full, efull, "Xpetra::EpetraMapextractor constructors only accepts Xpetra::EpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<EpetraVector> ret = Teuchos::rcp(new Xpetra::EpetraVector(getMap(block),true));
      ExtractVector(*efull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<VectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(EpetraVector, full, efull, "Xpetra::EpetraMapextractor constructors only accepts Xpetra::EpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<EpetraVector> ret = Teuchos::rcp(new Xpetra::EpetraVector(getMap(block),true));
      ExtractVector(*efull,block,*ret);
      return ret;
    }

    void ExtractVector(const Xpetra::EpetraMultiVector& full, size_t block, Xpetra::EpetraMultiVector& partial) const
    {
      if (maps_[block] == Teuchos::null)
        std::cout << "null map at block " << block << std::endl;
      partial.doImport(full, *importer_[block], Xpetra::INSERT);
    }

    virtual void ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const EpetraMultiVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(      EpetraMultiVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      ExtractVector(*efull,block,*epartial);
    }

    virtual void ExtractVector(Teuchos::RCP<      MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const
    {
      XPETRA_RCP_DYNAMIC_CAST(      EpetraMultiVector, full, efull, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");
      XPETRA_RCP_DYNAMIC_CAST(      EpetraMultiVector, partial, epartial, "Xpetra::EpetraMapextractor constructors only accept Xpetra::EpetraMap as input arguments.");

      ExtractVector(*efull,block,*epartial);
    }

    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(const EpetraMultiVector, full, efull, "Xpetra::EpetraMapextractor constructors only accepts Xpetra::EpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<EpetraMultiVector> ret = Teuchos::rcp(new Xpetra::EpetraMultiVector(getMap(block),full->getNumVectors(),true));
      ExtractVector(*efull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<MultiVectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(EpetraMultiVector, full, efull, "Xpetra::EpetraMapextractor constructors only accepts Xpetra::EpetraMap as input arguments.");
      TEUCHOS_TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<EpetraMultiVector> ret = Teuchos::rcp(new Xpetra::EpetraMultiVector(getMap(block),full->getNumVectors(),true));
      ExtractVector(*efull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<VectorClass> getVector(size_t i) const
    {
      return Teuchos::rcp(new Xpetra::EpetraVector(getMap(i),true));
    }

    virtual Teuchos::RCP<MultiVectorClass> getVector(size_t i, size_t numvec) const
    {
      return Teuchos::rcp(new Xpetra::EpetraMultiVector(getMap(i),numvec,true));
    }

    /** \name Maps */
    //@{

    /// number of partial maps
    size_t NumMaps() const { return maps_.size(); }

    /// get the map
    const Teuchos::RCP<const MapClass> getMap(size_t i) const { return maps_[i]; }

    /// the full map
    const Teuchos::RCP<const MapClass> getFullMap() const { return fullmap_; }

    size_t getMapIndexForGID(int gid) const {

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

        Teuchos::ArrayView< const int > mapGids = map->getNodeElementList();
        Teuchos::ArrayView< const int >::const_iterator it;
        for(it = mapGids.begin(); it!=mapGids.end(); it++) {
          if(fullMap->isNodeGlobalElement(*it)==false)
            return false; // Global ID (*it) not found locally on this proc in fullMap -> error
        }
      }
      return true;
    }

  protected:
    std::vector<Teuchos::RCP<const Xpetra::EpetraMap > > maps_;
    Teuchos::RCP<const EpetraMap > fullmap_;
    std::vector<Teuchos::RCP<ImportClass > > importer_;
  };
}


#endif /* XPETRA_EPETRAMAPEXTRACTOR_HPP_ */
