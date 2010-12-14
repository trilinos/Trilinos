#ifndef CTHULHU_EPETRA_IMPORT_HPP
#define CTHULHU_EPETRA_IMPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include <iterator>

#include "Cthulhu_Import.hpp"

namespace Cthulhu {

  //! \brief This class builds an object containing information necesary for efficiently importing off-processor entries.
  /*! Import is used to construct a communication plan that can be called repeatedly by computational
      classes to efficiently import entries from other nodes.
      For example, an exporter is used when we start out with a multiple-ownership distribution,
      and we want to merge that into a uniquely-owned distribution.

      This class currently has one constructor, taking two Map objects
      specifying the distributions of the distributed objects on which the Export class will operate.

      This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
      The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class EpetraImport<int,int>: public Import<int,int> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Import object from the source and target Maps.
    EpetraImport(const Teuchos::RCP<const Map<int,int> > & source, 
            const Teuchos::RCP<const Map<int,int> > & target)       
    { CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, source, tSource, "Cthulhu::EpetraImport constructors only accept Cthulhu::EpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, target, tTarget, "Cthulhu::EpetraImport constructors only accept Cthulhu::EpetraMap as input arguments.");
      import_ = rcp(new Epetra::Import<int,int>(tSource->getEpetra_Map(), tTarget->getEpetra_Map()));
    }
 
    //! copy constructor. 
    EpetraImport(const Import<int,int> & import) {
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraImport, import, timport, "Cthulhu::EpetraImport copy constructors only accept Cthulhu::EpetraImport as input arguments.");
      import_ = rcp(new Epetra::CrsGraph<int,int>(tImport->getEpetra_Map()));
    }

    EpetraImport(const Teuchos::RCP<Epetra::Import<int, int> > &import) : import_(import) { CTHULHU_DEBUG_ME; }

    //! destructor.
    ~EpetraImport() {}

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    inline size_t getNumSameIDs() const { CTHULHU_DEBUG_ME; return import_->getNumSameIDs(); }

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    inline size_t getNumPermuteIDs() const { CTHULHU_DEBUG_ME; return import_->getNumPermuteIDs(); }

    //! List of entries in the source Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const int> getPermuteFromLIDs() const { CTHULHU_DEBUG_ME; return import_->getPermuteFromLIDs(); }

    //! List of entries in the target Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const int> getPermuteToLIDs() const { CTHULHU_DEBUG_ME; return import_->getPermuteToLIDs(); }

    //! Returns the number of entries that are not on the calling image.
    inline size_t getNumRemoteIDs() const { CTHULHU_DEBUG_ME; return import_->getNumRemoteIDs(); }

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    inline Teuchos::ArrayView<const int> getRemoteLIDs() const { CTHULHU_DEBUG_ME; return import_->getRemoteLIDs(); }

    //! Returns the number of entries that must be sent by the calling image to other images.
    inline size_t getNumExportIDs() const { CTHULHU_DEBUG_ME; return import_->getNumExportIDs(); }

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportLIDs() const { CTHULHU_DEBUG_ME; return import_->getExportLIDs(); }

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportImageIDs() const { CTHULHU_DEBUG_ME; return import_->getExportImageIDs(); }

    //! Returns the Source Map used to construct this importer.
    inline const Teuchos::RCP<const Map<int,int> > getSourceMap() const { CTHULHU_DEBUG_ME; return import_->getSourceMap(); }

    //! Returns the Target Map used to construct this importer.
    inline const Teuchos::RCP<const Map<int,int> > getTargetMap() const { CTHULHU_DEBUG_ME; return import_->getTargetMap(); }

    inline Distributor & getDistributor() const { CTHULHU_DEBUG_ME; return import_->getDistributor(); }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Assignment operator
    inline Import<int,int>& operator = (const Import<int,int> & Source) { CTHULHU_DEBUG_ME; return }

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method 
    inline void print(std::ostream& os) const { CTHULHU_DEBUG_ME; } 

    //@}
#endif

    RCP< const Epetra::Import<int, int> > getEpetra_Import() const { CTHULHU_DEBUG_ME; return import_; }
    
  private:
    
    RCP<Epetra::Import> import_;

  };

} // namespace Cthulhu

#endif // CTHULHU_EPETRAIMPORT_HPP
