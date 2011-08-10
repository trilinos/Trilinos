#ifndef CTHULHU_EPETRA_IMPORT_HPP
#define CTHULHU_EPETRA_IMPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include <iterator>

#include "Cthulhu_Import.hpp"
#include "Epetra_Import.h"

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
  class EpetraImport: public Import<int,int> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Import object from the source and target Maps.
    EpetraImport(const Teuchos::RCP<const Map<int,int> > & source, 
		 const Teuchos::RCP<const Map<int,int> > & target)       
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, source, tSource, "Cthulhu::EpetraImport constructors only accept Cthulhu::EpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, target, tTarget, "Cthulhu::EpetraImport constructors only accept Cthulhu::EpetraMap as input arguments.");
      import_ = rcp(new Epetra_Import(tTarget->getEpetra_BlockMap(), tSource->getEpetra_BlockMap())); // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)
    }
 
    //! copy constructor. 
    EpetraImport(const Import<int,int> & import) 
      : import_(EpetraImportConstructorHelper(import))
    {
      
    }

    static RCP<Epetra_Import> EpetraImportConstructorHelper(const Import<int,int> & import) {
      
      CTHULHU_DYNAMIC_CAST(const EpetraImport, import, tImport, "Cthulhu::EpetraImport copy constructors only accept Cthulhu::EpetraImport as input arguments.");
      return rcp(new Epetra_Import(*tImport.getEpetra_Import()));
    }

    EpetraImport(const Teuchos::RCP<const Epetra_Import> &import) : import_(import) {  }

    //! destructor.
    ~EpetraImport() {}

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    inline size_t getNumSameIDs() const {  return import_->NumSameIDs(); }

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    inline size_t getNumPermuteIDs() const {  return import_->NumPermuteIDs(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! List of entries in the source Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const int> getPermuteFromLIDs() const { 
       
      return import_->PermuteFromLIDs(); 
    }

    //! List of entries in the target Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const int> getPermuteToLIDs() const { 
       
      return import_->PermuteToLIDs(); 
    }

    //! Returns the number of entries that are not on the calling image.
    inline size_t getNumRemoteIDs() const {  return import_->NumRemoteIDs(); }

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    inline Teuchos::ArrayView<const int> getRemoteLIDs() const { 
       
      return import_->RemoteLIDs(); 
    }

    //! Returns the number of entries that must be sent by the calling image to other images.
    inline size_t getNumExportIDs() const {  return import_->getNumExportIDs(); }

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportLIDs() const { 
       
      return import_->ExportLIDs(); 
    }

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportImageIDs() const { 
       
      return import_->ExportImageIDs(); 
    }
#endif

    //! Returns the Source Map used to construct this importer.
    inline const Teuchos::RCP<const Map<int,int> > getSourceMap() const { 
       

      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(import_->SourceMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the Target Map used to construct this importer.
    inline const Teuchos::RCP<const Map<int,int> > getTargetMap() const { 
       

      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(import_->TargetMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    inline Distributor & getDistributor() const {  return import_->Distributor(); }
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Assignment operator
    inline Import<int,int>& operator = (const Import<int,int> & Source) {  return }

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method 
    inline void print(std::ostream& os) const {  } 

    //@}
#endif

    RCP< const Epetra_Import > getEpetra_Import() const {  return import_; }
    
  private:
    
    RCP<const Epetra_Import> import_;

  };

} // namespace Cthulhu

#endif // CTHULHU_EPETRAIMPORT_HPP
