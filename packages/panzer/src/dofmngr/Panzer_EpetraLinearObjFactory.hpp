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
#include "Panzer_ScatterResidual_Epetra.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra.hpp"
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Teuchos_RCP.hpp"

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace panzer {

template <typename Traits,typename LocalOrdinalT>
class EpetraLinearObjFactory : public LinearObjFactory<Traits> {
public:

   EpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider);

   virtual ~EpetraLinearObjFactory();

   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > getGhostedVector() const;
   virtual Teuchos::RCP<Thyra::LinearOpBase<double> > getGhostedMatrix() const;

   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > getVector() const;
   virtual Teuchos::RCP<Thyra::LinearOpBase<double> > getMatrix() const;

   virtual void ghostToGlobalMatrix(const Thyra::LinearOpBase<double> & ghostA, 
                                    Thyra::LinearOpBase<double> & A) const;

   virtual void ghostToGlobalVector(const Thyra::MultiVectorBase<double> & ghostA, 
                                    Thyra::MultiVectorBase<double> & A) const;

   virtual void globalToGhostVector(const Thyra::MultiVectorBase<double> & A, 
                                    Thyra::MultiVectorBase<double> & ghostA) const;

   /** Do a simple assignment to a linear operator. The intention is that this
     * sets up the operator to be filled.
     */
   virtual void assignToMatrix(Thyra::LinearOpBase<double> & oper,double value);

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

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   { return Teuchos::rcp(new ScatterResidual_Epetra<EvalT,Traits>); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits>); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return Teuchos::rcp(new ScatterDirichletResidual_Epetra<EvalT,Traits>); }

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
