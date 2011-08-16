#ifndef CTHULHU_EPETRA_IMPORT_HPP
#define CTHULHU_EPETRA_IMPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Cthulhu_Map.hpp"
#include "Cthulhu_EpetraMap.hpp"
#include <iterator>

#include "Cthulhu_Import.hpp"
#include "Epetra_Import.h"

namespace Cthulhu {

  // TODO: move that elsewhere
  //   const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & toTpetra(const Import<LocalOrdinal,GlobalOrdinal,Node> &import);

  //  RCP< const Import<int, int > > toCthulhu(const Epetra_Import &import);
  RCP< const Import<int, int > > toCthulhu(const Epetra_Import *import);
  //

  //! \brief This class builds an object containing information necesary for efficiently importing off-processor entries.
  class EpetraImport
    : public Import<int,int>
  {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Import object from the source and target Maps.
    EpetraImport(const Teuchos::RCP<const Map<int,int> > & source, const Teuchos::RCP<const Map<int,int> > & target)       
    { 
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, source, tSource, "Cthulhu::EpetraImport constructors only accept Cthulhu::EpetraMap as input arguments.");
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, target, tTarget, "Cthulhu::EpetraImport constructors only accept Cthulhu::EpetraMap as input arguments.");
      import_ = rcp(new Epetra_Import(tTarget->getEpetra_BlockMap(), tSource->getEpetra_BlockMap())); // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)
    }
 
    //! copy constructor. 
    EpetraImport(const Import<int,int> & import) : import_(EpetraImportConstructorHelper(import)) { }
    static RCP<Epetra_Import> EpetraImportConstructorHelper(const Import<int,int> & import) {
      CTHULHU_DYNAMIC_CAST(const EpetraImport, import, tImport, "Cthulhu::EpetraImport copy constructors only accept Cthulhu::EpetraImport as input arguments.");
      return rcp(new Epetra_Import(*tImport.getEpetra_Import()));
    }

    //! destructor.
    ~EpetraImport() {}

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    size_t getNumSameIDs() const {  return import_->NumSameIDs(); }

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    size_t getNumPermuteIDs() const {  return import_->NumPermuteIDs(); }

    //! Returns the Source Map used to construct this importer.
    const Teuchos::RCP<const Map<int,int> > getSourceMap() const { 
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(import_->SourceMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //! Returns the Target Map used to construct this importer.
    const Teuchos::RCP<const Map<int,int> > getTargetMap() const { 
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(import_->TargetMap()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    //@}

    //! @name Cthulhu specific
    //@{

    //! EpetraImport constructor to wrap a Epetra_Import object
    EpetraImport(const Teuchos::RCP<const Epetra_Import> &import) : import_(import) { }

    //! Get the underlying Epetra import
    RCP< const Epetra_Import > getEpetra_Import() const {  return import_; }

    //@}
    
  private:
    
    RCP<const Epetra_Import> import_;

  };

} // namespace Cthulhu

#endif // CTHULHU_EPETRAIMPORT_HPP
