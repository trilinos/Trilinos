#ifndef PANZER_MODEL_EVALUATOR_EPETRA_HPP
#define PANZER_MODEL_EVALUATOR_EPETRA_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"

namespace panzer {

  template<typename, typename>  class FieldManagerBuilder;
  template<typename, typename>  class EpetraLinearObjFactory;
  class EpetraLinearObjContainer;

  class ModelEvaluator_Epetra : public EpetraExt::ModelEvaluator {
  public:
    
    ModelEvaluator_Epetra(Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb,
			  Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > lof);
    
    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{
    
    Teuchos::RCP<const Epetra_Map> get_x_map() const;
    Teuchos::RCP<const Epetra_Map> get_f_map() const;
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;
    Teuchos::RCP<Epetra_Operator> create_W() const;
    InArgs createInArgs() const;
    OutArgs createOutArgs() const;
    void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;
    
    //@}

  private:

    // /////////////////////////////////////
    // Private member data
    
    double    d_;
    
    bool      isInitialized_;
    
    Teuchos::RCP<const Epetra_Map>   map_x_;
    
    Teuchos::RCP<Epetra_Vector> x0_;
    Teuchos::RCP<Epetra_Vector> p_;
  
    Teuchos::RCP<Epetra_CrsGraph>  W_graph_;
    
    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb_;
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > lof_;

    mutable panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm_;
    mutable Teuchos::RCP<panzer::AssemblyEngineInArgs> ae_inargs_;
    mutable Teuchos::RCP<EpetraLinearObjContainer> ghostedContainer_, container_;
  };
  
}

#endif
