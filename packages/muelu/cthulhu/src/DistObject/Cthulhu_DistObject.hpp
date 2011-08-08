#ifndef CTHULHU_DISTOBJECT_HPP
#define CTHULHU_DISTOBJECT_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_Import.hpp"
#include "Cthulhu_Export.hpp"
//#include "Cthulhu_Distributor.hpp"

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

namespace Cthulhu {

  //! A base class for distributed objects that support import and export operations.

  /*! The DistObject is a base class for all Cthulhu distributed global objects.  It provides the basic
      mechanisms and interface specifications for importing and exporting operations using Cthulhu::Import and
      Cthulhu::Export objects.
    
    <b> Distributed Global vs. Replicated Local.</b>
    
    <ul>
    <li> Distributed Global objects - In most instances, a distributed object will be partitioned
    across multiple memory images associated with multiple processors.  In this case, there is 
    a unique copy of each element and elements are spread across all images specified by 
    the Teuchos::Comm object.
    <li> Replicated Local Objects - Some algorithms use objects that are too small to
    be distributed across all processors, the Hessenberg matrix in a GMRES
    computation.  In other cases, such as with block iterative methods,  block dot product 
    functions produce small dense matrices that are required by all images.  
    Replicated local objects handle these types of situation.
    </ul>
  */

  template <class Packet, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class DistObject : virtual public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //     //! constructor
    //     explicit DistObject(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map);
    
    //     //! copy constructor
    //     DistObject(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source);

    //! destructor
    virtual ~DistObject() { CTHULHU_DEBUG_ME; }

    //@}

    //! @name Import/Export Methods
    //@{ 

    //! Import
    virtual void doImport(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source, 
                          const Import<LocalOrdinal,GlobalOrdinal,Node> &importer, CombineMode CM) = 0;

    //! Export
    virtual void doExport(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &dest, 
                          const Export<LocalOrdinal,GlobalOrdinal,Node> &exporter, CombineMode CM) = 0;

    //! Import (using an Exporter)
    virtual void doImport(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source,
                          const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, CombineMode CM) = 0;

    //! Export (using an Importer)
    virtual void doExport(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &dest,
                          const Import<LocalOrdinal,GlobalOrdinal,Node>& importer, CombineMode CM) = 0;

    //@}

    //! @name Attribute Accessor Methods
    //@{ 

    //! Accessor for whether or not this is a global object
//     virtual bool isDistributed() const = 0;

    //! Access function for the Cthulhu::Map this DistObject was constructed with.
    virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > /*&*/ getMap() const = 0;

    //@}

    //! @name I/O methods
    //@{ 

    //! Print method.

//     virtual void print(std::ostream &os) const = 0;

    //@} 

  }; // class DistObject

} // namespace Cthulhu

#endif /* CTHULHU_DISTOBJECT_HPP */
