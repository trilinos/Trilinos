#ifndef CTHULHU_IMPORT_HPP
#define CTHULHU_IMPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include <iterator>

namespace Cthulhu {

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Import
    : public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    // !Destructor.
    virtual ~Import() { }

    //@}

    //! @name Export Attribute Methods
    //@{

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    virtual size_t getNumSameIDs() const = 0;

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    virtual size_t getNumPermuteIDs() const = 0;

    //! Returns the Source Map used to construct this importer.
    virtual const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getSourceMap() const = 0;

    //! Returns the Target Map used to construct this importer.
    virtual const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getTargetMap() const = 0;

    //@}

  }; // CrsMatrix class

} // Cthulhu namespace

#define CTHULHU_IMPORT_SHORT
#endif // CTHULHU_IMPORT_HPP
