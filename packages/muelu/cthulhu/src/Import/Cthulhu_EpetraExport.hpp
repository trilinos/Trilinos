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
  class EpetraExport
    : public Export<int,int>
  {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

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

    //! destructor.
    ~EpetraExport() {}

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    size_t getNumSameIDs() const {  return export_->NumSameIDs(); }

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    size_t getNumPermuteIDs() const {  return export_->NumPermuteIDs(); }

    //! Returns the Source Map used to construct this exporter.
    const Teuchos::RCP<const Map<int,int> > getSourceMap() const { 
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(export_->SourceMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the Target Map used to construct this exporter.
    const Teuchos::RCP<const Map<int,int> > getTargetMap() const { 
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(export_->TargetMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! @name Cthulhu specific
    //@{

    //! EpetraExport constructor to wrap a Epetra_Export object
    EpetraExport(const Teuchos::RCP<const Epetra_Export> &export2) : export_(export2) {  }

    //! Get the underlying Epetra export
    RCP< const Epetra_Export > getEpetra_Export() const {  return export_; }

    //@}

    
  private:
    
    RCP<const Epetra_Export> export_;

  };

} // namespace Cthulhu

#endif // CTHULHU_EPETRAEXPORT_HPP
