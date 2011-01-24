#ifndef __Panzer_EpetraLinearObjFactory_hpp__
#define __Panzer_EpetraLinearObjFactory_hpp__

#include <map>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "Panzer_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"

#include "Teuchos_RCP.hpp"

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace panzer {

template <typename LocalOrdinalT>
class EpetraLinearObjFactory : public LinearObjFactory {
public:

   EpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider);

   virtual ~EpetraLinearObjFactory();

   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > getGhostedVector() const;
   virtual Teuchos::RCP<Thyra::LinearOpBase<double> > getGhostedMatrix() const;

   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > getVector() const;
   virtual Teuchos::RCP<Thyra::LinearOpBase<double> > getMatrix() const;

   virtual void ghostToGlobalMatrix(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & ghostA, 
                                      const Teuchos::RCP<Thyra::LinearOpBase<double> > & A) const;

   virtual void ghostToGlobalVector(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & ghostA, 
                                    const Teuchos::RCP<Thyra::MultiVectorBase<double> > & A) const;

   virtual void globalToGhostVector(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & A, 
                                    const Teuchos::RCP<Thyra::MultiVectorBase<double> > & ghostA) const;

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getMap() const;

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedMap() const;

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGraph() const;

   //! get the ghosted graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGhostedGraph() const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getGhostedImport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getGhostedExport() const;

protected:
   // get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> buildMap() const;
   virtual const Teuchos::RCP<Epetra_Map> buildGhostedMap() const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGraph() const;
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGhostedGraph() const;

   // storage for Epetra graphs and maps
   Teuchos::RCP<const Epetra_Comm> comm_;
   mutable Teuchos::RCP<Epetra_Map> map_;
   mutable Teuchos::RCP<Epetra_Map> ghostedMap_;
   mutable Teuchos::RCP<Epetra_CrsGraph> graph_;
   mutable Teuchos::RCP<Epetra_CrsGraph> ghostedGraph_;

   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > gidProvider_;
};

}

#include "Panzer_EpetraLinearObjFactoryT.hpp"

#endif
