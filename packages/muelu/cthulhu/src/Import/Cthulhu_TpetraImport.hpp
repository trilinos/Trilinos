#ifndef CTHULHU_TPETRA_IMPORT_HPP
#define CTHULHU_TPETRA_IMPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include <iterator>
#include "Cthulhu_Import.hpp"
#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class TpetraImport: public Import<LocalOrdinal, GlobalOrdinal, Node> {

    // The following typedef are used by the CTHULHU_DYNAMIC_CAST() macro.
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraImport<LocalOrdinal, GlobalOrdinal, Node> TTpetraImport;

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Import object from the source and target Maps.
    TpetraImport(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & source, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & target)       
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, source, tSource, "Cthulhu::TpetraImport constructors only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, target, tTarget, "Cthulhu::TpetraImport constructors only accept Cthulhu::TpetraMap as input arguments.");
      import_ = rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(tSource->getTpetra_Map(), tTarget->getTpetra_Map()));
    }
 
    //! copy constructor. 
    TpetraImport(const Import<LocalOrdinal,GlobalOrdinal,Node> & import) {
      CTHULHU_RCP_DYNAMIC_CAST(const TTpetraImport, import, tImport, "Cthulhu::TpetraImport copy constructors only accept Cthulhu::TpetraImport as input arguments.");
      import_ = rcp(new TpetraImport<LocalOrdinal,GlobalOrdinal,Node>(tImport->getTpetra_Map()));
    }

    TpetraImport(const Teuchos::RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import) : import_(import) {  }

    //! destructor.
    virtual ~TpetraImport() {  }

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    inline size_t getNumSameIDs() const { return import_->getNumSameIDs(); }

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    inline size_t getNumPermuteIDs() const { return import_->getNumPermuteIDs(); }

    //! List of entries in the source Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getPermuteFromLIDs() const { return import_->getPermuteFromLIDs(); }

    //! List of entries in the target Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getPermuteToLIDs() const { return import_->getPermuteToLIDs(); }

    //! Returns the number of entries that are not on the calling image.
    inline size_t getNumRemoteIDs() const { return import_->getNumRemoteIDs(); }

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getRemoteLIDs() const { return import_->getRemoteLIDs(); }

    //! Returns the number of entries that must be sent by the calling image to other images.
    inline size_t getNumExportIDs() const { return import_->getNumExportIDs(); }

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getExportLIDs() const { return import_->getExportLIDs(); }

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportImageIDs() const { return import_->getExportImageIDs(); }

    //! Returns the Source Map used to construct this importer.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getSourceMap() const { return rcp( new TpetraMapClass(import_->getSourceMap())); };

    //! Returns the Target Map used to construct this importer.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getTargetMap() const { return rcp( new TpetraMapClass(import_->getTargetMap())); };

    RCP< const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_Import() const { return import_; }
    
  private:
    
    RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > import_;

  };

} // namespace Cthulhu

#define CTHULHU_TPETRAIMPORT_SHORT
#endif // CTHULHU_TPETRAIMPORT_HPP
