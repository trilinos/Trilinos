#ifndef CTHULHU_DISTOBJECT_HPP
#define CTHULHU_DISTOBJECT_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_Import.hpp"
#include "Cthulhu_Export.hpp"
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

namespace Cthulhu {

  template <class Packet, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class DistObject
    : virtual public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Destructor.
    virtual ~DistObject() { }

    //@}

    //! @name Import/Export Methods
    //@{

    //! Import.
    virtual void doImport(const DistObject< Packet, LocalOrdinal, GlobalOrdinal, Node > &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)= 0;

    //! Export.
    virtual void doExport(const DistObject< Packet, LocalOrdinal, GlobalOrdinal, Node > &dest, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)= 0;

    //! Import (using an Exporter).
    virtual void doImport(const DistObject< Packet, LocalOrdinal, GlobalOrdinal, Node > &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)= 0;

    //! Export (using an Importer).
    virtual void doExport(const DistObject< Packet, LocalOrdinal, GlobalOrdinal, Node > &dest, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)= 0;

    //@}

    //! @name Attribute Accessor Methods
    //@{

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    virtual const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getMap() const = 0;

    //@}

  }; // DistObject class

} // Cthulhu namespace

#define CTHULHU_DISTOBJECT_SHORT
#endif // CTHULHU_DISTOBJECT_HPP
