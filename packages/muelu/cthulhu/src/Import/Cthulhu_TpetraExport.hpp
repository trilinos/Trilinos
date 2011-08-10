#ifndef CTHULHU_TPETRA_EXPORT_HPP
#define CTHULHU_TPETRA_EXPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include <iterator>

#include "Cthulhu_Export.hpp"
#include "Cthulhu_Exceptions.hpp"

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
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class TpetraExport: public Export<LocalOrdinal, GlobalOrdinal, Node> {

    // The following typedef are used by the CTHULHU_DYNAMIC_CAST() macro.
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraExport<LocalOrdinal, GlobalOrdinal, Node> TTpetraExport;

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Export object from the source and target Maps.
    TpetraExport(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & source, 
            const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & target)       
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, source, tSource, "Cthulhu::TpetraExport constructors only accept Cthulhu::TpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMapClass, target, tTarget, "Cthulhu::TpetraExport constructors only accept Cthulhu::TpetraMap as input arguments.");
      export_ = rcp(new Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node>(tSource->getTpetra_Map(), tTarget->getTpetra_Map()));
    }
 
    //! copy constructor. 
    TpetraExport(const Export<LocalOrdinal,GlobalOrdinal,Node> & export2) {
      
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraExport, export2, tExport, "Cthulhu::TpetraExport copy constructors only accept Cthulhu::TpetraExport as input arguments.");
      export_ = rcp(new TpetraExport<LocalOrdinal,GlobalOrdinal,Node>(tExport->getTpetra_Map()));
    }

    TpetraExport(const Teuchos::RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node> > &export2) : export_(export2) {  }

    //! destructor.
    virtual ~TpetraExport() {  }

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    inline size_t getNumSameIDs() const {  return export_->getNumSameIDs(); }

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    inline size_t getNumPermuteIDs() const {  return export_->getNumPermuteIDs(); }

    //! List of entries in the source Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getPermuteFromLIDs() const {  return export_->getPermuteFromLIDs(); }

    //! List of entries in the target Map that are permuted. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getPermuteToLIDs() const {  return export_->getPermuteToLIDs(); }

    //! Returns the number of entries that are not on the calling image.
    inline size_t getNumRemoteIDs() const {  return export_->getNumRemoteIDs(); }

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getRemoteLIDs() const {  return export_->getRemoteLIDs(); }

    //! Returns the number of entries that must be sent by the calling image to other images.
    inline size_t getNumExportIDs() const {  return export_->getNumExportIDs(); }

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    inline Teuchos::ArrayView<const LocalOrdinal> getExportLIDs() const {  return export_->getExportLIDs(); }

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    inline Teuchos::ArrayView<const int> getExportImageIDs() const {  return export_->getExportImageIDs(); }

    //! Returns the Source Map used to construct this exporter.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getSourceMap() const {  return rcp( new TpetraMapClass(export_->getSourceMap())); };

    //! Returns the Target Map used to construct this exporter.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getTargetMap() const {  return rcp( new TpetraMapClass(export_->getTargetMap())); };

#ifdef CTHULHU_NOT_IMPLEMENTED
    inline Distributor & getDistributor() const {  return export_->getDistributor(); }
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Assignment operator
    inline Export<LocalOrdinal,GlobalOrdinal,Node>& operator = (const Export<LocalOrdinal,GlobalOrdinal,Node> & Source) {  return }

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method 
    inline void print(std::ostream& os) const {  } 

    //@}
#endif

    RCP< const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_Export() const {  return export_; }
    
  private:
    
    RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node> > export_;

  };

} // namespace Cthulhu

#define CTHULHU_TPETRAEXPORT_SHORT
#endif // CTHULHU_TPETRAEXPORT_HPP
