#ifndef CTHULHU_EPETRA_EXPORT_HPP
#define CTHULHU_EPETRA_EXPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include <iterator>

#include "Cthulhu_Export.hpp"
#include "Epetra_Export.h"

namespace Cthulhu {

  //! \brief This class builds an object containing information necesary for efficiently exporting off-processor entries.
  /*! Export is used to construct a communication plan that can be called repeatedly by computational
      classes to efficiently export entries from other nodes.
      For example, an exporter is used when we start out with a multiple-ownership distribution,
      and we want to merge that into a uniquely-owned distribution.

      This class currently has one constructor, taking two Map objects
      specifying the distributions of the distributed objects on which the Export class will operate.

      This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
      The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
  */
  class EpetraExport: public Export<int,int> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Export object from the source and target Maps.
    EpetraExport(const Teuchos::RCP<const Map<int,int> > & source, 
		 const Teuchos::RCP<const Map<int,int> > & target)       
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, source, tSource, "Cthulhu::EpetraExport constructors only accept Cthulhu::EpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, target, tTarget, "Cthulhu::EpetraExport constructors only accept Cthulhu::EpetraMap as input arguments.");
      export_ = rcp(new Epetra_Export(tTarget->getEpetra_BlockMap(), tSource->getEpetra_BlockMap())); // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)
    }
 
    //! copy constructor. 
    EpetraExport(const Export<int,int> & export2) 
      : export_(EpetraExportConstructorHelper(export2))
    {
      
    }

    static RCP<Epetra_Export> EpetraExportConstructorHelper(const Export<int,int> & export2) {
      
      CTHULHU_DYNAMIC_CAST(const EpetraExport, export2, tExport, "Cthulhu::EpetraExport copy constructors only accept Cthulhu::EpetraExport as input arguments.");
      return rcp(new Epetra_Export(*tExport.getEpetra_Export()));
    }

    EpetraExport(const Teuchos::RCP<const Epetra_Export> &export2) : export_(export2) {  }

    //! destructor.
    ~EpetraExport() {}

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    inline size_t getNumSameIDs() const {  return export_->NumSameIDs(); }

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    inline size_t getNumPermuteIDs() const {  return export_->NumPermuteIDs(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! List of entries in the source Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const int> getPermuteFromLIDs() const { 
       
      return export_->PermuteFromLIDs(); 
    }

    //! List of entries in the target Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const int> getPermuteToLIDs() const { 
       
      return export_->PermuteToLIDs(); 
    }

    //! Returns the number of entries that are not on the calling image.
    inline size_t getNumRemoteIDs() const {  return export_->NumRemoteIDs(); }

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    inline Teuchos::ArrayView<const int> getRemoteLIDs() const { 
       
      return export_->RemoteLIDs(); 
    }

    //! Returns the number of entries that must be sent by the calling image to other images.
    inline size_t getNumExportIDs() const {  return export_->getNumExportIDs(); }

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportLIDs() const { 
       
      return export_->ExportLIDs(); 
    }

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportImageIDs() const { 
       
      return export_->ExportImageIDs(); 
    }
#endif

    //! Returns the Source Map used to construct this exporter.
    inline const Teuchos::RCP<const Map<int,int> > getSourceMap() const { 
       

      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(export_->SourceMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the Target Map used to construct this exporter.
    inline const Teuchos::RCP<const Map<int,int> > getTargetMap() const { 
       

      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(export_->TargetMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    inline Distributor & getDistributor() const {  return export_->Distributor(); }
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Assignment operator
    inline Export<int,int>& operator = (const Export<int,int> & Source) {  return }

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method 
    inline void print(std::ostream& os) const {  } 

    //@}
#endif

    RCP< const Epetra_Export > getEpetra_Export() const {  return export_; }
    
  private:
    
    RCP<const Epetra_Export> export_;

  };

} // namespace Cthulhu

#endif // CTHULHU_EPETRAEXPORT_HPP
