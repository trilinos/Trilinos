#ifndef CTHULHU_IMPORT_HPP
#define CTHULHU_IMPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include <iterator>

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
  class Import: public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! destructor.
    virtual ~Import() { CTHULHU_DEBUG_ME; }

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    virtual size_t getNumSameIDs() const = 0;

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    virtual size_t getNumPermuteIDs() const = 0;

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! List of entries in the source Map that are permuted. (non-persisting view)
    virtual Teuchos::ArrayView<const LocalOrdinal> getPermuteFromLIDs() const = 0;

    //! List of entries in the target Map that are permuted. (non-persisting view)
    virtual Teuchos::ArrayView<const LocalOrdinal> getPermuteToLIDs() const = 0;

    //! Returns the number of entries that are not on the calling image.
    virtual size_t getNumRemoteIDs() const = 0;

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    virtual Teuchos::ArrayView<const LocalOrdinal> getRemoteLIDs() const = 0;

    //! Returns the number of entries that must be sent by the calling image to other images.
    virtual size_t getNumExportIDs() const = 0;

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    virtual Teuchos::ArrayView<const LocalOrdinal> getExportLIDs() const = 0;

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    virtual Teuchos::ArrayView<const int> getExportImageIDs() const = 0;
#endif

    //! Returns the Source Map used to construct this importer.
    virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getSourceMap() const = 0;

    //! Returns the Target Map used to construct this importer.
    virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getTargetMap() const = 0;

#ifdef CTHULHU_NOT_IMPLEMENTED
    virtual Distributor & getDistributor() const = 0;
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Assignment operator
    virtual Import<LocalOrdinal,GlobalOrdinal,Node>& operator = (const Import<LocalOrdinal,GlobalOrdinal,Node> & Source) = 0;

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method 
    virtual void print(std::ostream& os) const = 0;

    //@}
#endif

  };

} // namespace Cthulhu

#define CTHULHU_IMPORT_SHORT
#endif // CTHULHU_IMPORT_HPP
