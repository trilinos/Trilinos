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
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_ScatterResidual_Epetra.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra.hpp"
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

template <typename Traits,typename LocalOrdinalT>
class EpetraLinearObjFactory : public LinearObjFactory<Traits> {
public:

   EpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider);

   virtual ~EpetraLinearObjFactory();

/*************** Linear object factory methods *******************/

   virtual Teuchos::RCP<LinearObjContainer> buildLinearObjContainer() const;

   virtual Teuchos::RCP<LinearObjContainer> buildGhostedLinearObjContainer() const;

   virtual void globalToGhostContainer(const LinearObjContainer & container,
                                       LinearObjContainer & ghostContainer) const;
   virtual void ghostToGlobalContainer(const LinearObjContainer & ghostContainer,
                                       LinearObjContainer & container) const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   { return Teuchos::rcp(new ScatterResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return Teuchos::rcp(new ScatterDirichletResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

/*************** Epetra based methods *******************/

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
   Teuchos::RCP<Epetra_Vector> getGhostedEpetraVector() const;
   Teuchos::RCP<Epetra_Vector> getEpetraVector() const;
   Teuchos::RCP<Epetra_CrsMatrix> getEpetraMatrix() const;
   Teuchos::RCP<Epetra_CrsMatrix> getGhostedEpetraMatrix() const;

   void ghostToGlobalEpetraVector(const Epetra_Vector in,Epetra_Vector & out) const;
   void ghostToGlobalEpetraMatrix(const Epetra_CrsMatrix in,Epetra_CrsMatrix & out) const;
   void globalToGhostEpetraVector(const Epetra_Vector in,Epetra_Vector & out) const;

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
