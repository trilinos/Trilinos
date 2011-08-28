/*
 * Xpetra_TpetraMapExtractor.hpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

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
          if (tImport == Teuchos::null) cout << "no TpetraImportFactory?" << endl;
          importer_[i] = tImport;
          std::cout << *importer_[i] << std::endl;
        }
      }
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
      TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<TpetraVectorClass> ret = Teuchos::rcp(new TpetraVectorClass(getMap(block),true));
      ExtractVector(*tfull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<VectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
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
      TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
            "ExtractVector: maps_[" << block << "]=Teuchos::null" );
      const Teuchos::RCP<TpetraMultiVectorClass> ret = Teuchos::rcp(new TpetraMultiVectorClass(getMap(block),full->getNumVectors(),true));
      ExtractVector(*tfull,block,*ret);
      return ret;
    }

    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<MultiVectorClass>& full, size_t block) const
    {
      XPETRA_RCP_DYNAMIC_CAST(TpetraMultiVectorClass, full, tfull, "Xpetra::TpetraMapextractor constructors only accepts Xpetra::TpetraMap as input arguments.");
      TEST_FOR_EXCEPTION( maps_[block] == Teuchos::null, Xpetra::Exceptions::RuntimeError,
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
    const Teuchos::RCP<const MapClass> FullMap() const { return fullmap_; }

    //@}

  protected:
    std::vector<Teuchos::RCP<const TpetraMapClass > > maps_;
    Teuchos::RCP<const TpetraMapClass > fullmap_;
    std::vector<Teuchos::RCP<TpetraImportClass > > importer_;
  };
}


#endif /* XPETRA_TPETRAMAPEXTRACTOR_HPP_ */
