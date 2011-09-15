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
#include <vector>
#include <string>

namespace panzer {

  template<typename, typename>  class FieldManagerBuilder;
  template<typename, typename>  class EpetraLinearObjFactory;
  template<typename, typename>  class SGEpetraLinearObjFactory;
  template<typename>  class ResponseLibrary;
  class EpetraLinearObjContainer;
  class SGEpetraLinearObjContainer;

  class ModelEvaluator_Epetra : public EpetraExt::ModelEvaluator {
  public:
    
    ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder<int,int> >& fmb,
                          const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
			  const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof,
			  const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
			  bool build_transient_support);

    ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder<int,int> >& fmb,
                          const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
			  const Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> >& sg_lof,
			  const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
			  bool build_transient_support);
    
    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{
    
    Teuchos::RCP<const Epetra_Map> get_x_map() const;
    Teuchos::RCP<const Epetra_Map> get_f_map() const;
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;
    Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;
    Teuchos::RCP<Epetra_Operator> create_W() const;
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
    Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
    Teuchos::RCP<const Epetra_Map> get_g_map(int l) const;
    InArgs createInArgs() const;
    OutArgs createOutArgs() const;
    void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

    
    //@}

  private:

    // /////////////////////////////////////
    // Private methods

    /** Initialize Epetra linear objects.
      *
      * \note Requires lof_ to be set.
      */
    void initializeEpetraObjs();

    // /////////////////////////////////////
    // Private evaluation methods

    //! for evaluation and handling of normal quantities, x,f,W, etc
    void evalModel_basic( const InArgs& inArgs, const OutArgs& outArgs ) const; 

    /** handles evaluation of responses g, dgdx
      *
      * \note This method should (basically) be a no-op if <code>required_basic_g(outArgs)==false</code>.
      *       However, for efficiency this is not checked.
      */
    void evalModel_basic_g(AssemblyEngineInArgs ae_inargs,const InArgs & inArgs,const OutArgs & outArgs) const;

    //! Are their required responses in the out args? g (and soon DgDx) 
    bool required_basic_g(const OutArgs & outArgs) const;

    #ifdef HAVE_STOKHOS
       /** for evaluation and handling of Stochastic Galerkin quantities, x_sg, f_sg, W_sg, etc
         *
         * \note A precondition for this is that <code>sg_lof_</code> has been initialized
         *       with a call to the appropriate constructor.
         */
       void evalModel_sg( const InArgs& inArgs, const OutArgs& outArgs ) const;

       //! Are their required responses in the out args? g (and soon DgDx) 
       bool required_sg_g(const OutArgs & outArgs) const;

       /** handles evaluation of responses g, dgdx
         *
         * \note This method should (basically) be a no-op if <code>required_basic_g(outArgs)==false</code>.
         *       However, for efficiency this is not checked.
         */
       void evalModel_sg_g(AssemblyEngineInArgs ae_inargs,const InArgs & inArgs,const OutArgs & outArgs) const;
    #endif

    // /////////////////////////////////////
    // Private member data
    
    /** \defgroup EpetraObjs Underlying epetra types
      * @{ 
      */
    Teuchos::RCP<const Epetra_Map>   map_x_;
    Teuchos::RCP<Epetra_Vector> x0_;
    Teuchos::RCP<Epetra_Vector> x_dot_init_;
    mutable Teuchos::RCP<Epetra_Vector> dummy_f_;    
    
    // parameters
    std::vector<Teuchos::RCP<Epetra_Map> > p_map_;
    std::vector<Teuchos::RCP<Epetra_Vector> > p_init_;

    // responses
    std::vector<Teuchos::RCP<Epetra_Map> > g_map_;
    
    Teuchos::RCP<Epetra_CrsGraph>  W_graph_;

    /** @} */
    
    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb_;
    mutable Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > responseLibrary_; // These objects are basically the same
    mutable panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm_;   // they control and provide access to evaluate
    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names_;
    bool build_transient_support_;

    // basic specific linear object objects
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > lof_;
    mutable Teuchos::RCP<EpetraLinearObjContainer> ghostedContainer_;

    // sg specific linear object objects
    Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> > sg_lof_;
    mutable Teuchos::RCP<SGEpetraLinearObjContainer> sg_ghostedContainer_;
  };
  
}

#endif
