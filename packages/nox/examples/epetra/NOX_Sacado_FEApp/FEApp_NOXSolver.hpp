#ifndef FEAPP_NOXSOLVER_H
#define FEAPP_NOXSOLVER_H

#include "EpetraExt_ModelEvaluator.h"
#include "LOCA.H"
#include "LOCA_Epetra.H"

namespace FEApp {

  //! Version of Albany's ENAT_NOXSolver modified for SG
  /*
   * This class converts a ModelEvaluator for f(x,p) = 0, g = g(x,p) to
   * g = g(p) by elminating x in f(x,p) = 0 with the nonlinear solver NOX.
   */
  class NOXSolver : public EpetraExt::ModelEvaluator {
  public:

    /** \name Constructors/initializers */
    //@{
    
    /** \brief Takes the number of elements in the discretization . */
    NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
	      Teuchos::RCP<EpetraExt::ModelEvaluator> model,
	      Teuchos::RCP<Epetra_Operator> M = Teuchos::null);
    
    
    //@}
    
    ~NOXSolver();
    
    
    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{
    
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
    //! Return array of parameter names
    Teuchos::RCP<const Teuchos::Array<std::string> > 
    get_p_names(int l) const;
    /** \brief . */
    //  Teuchos::RCP<Epetra_Operator> create_W() const;
    /** \brief . */
    EpetraExt::ModelEvaluator::InArgs createInArgs() const;
    /** \brief . */
    EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
    /** \brief . */
    void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

    //! Get solution from solver
    Teuchos::RCP<const Epetra_Vector> getFinalSolution() const;

    //! Get solver status
    NOX::StatusTest::StatusType getSolverStatus() const;
    
  private:
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_x_map() const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_f_map() const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;
    
    
    //@}
    
  private:
    
    //These are set in the constructor and used in evalModel
    mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
    Teuchos::RCP<EpetraExt::ModelEvaluator> model;
    
    Teuchos::RCP<NOX::Solver::Generic> solver;
    Teuchos::RCP<const Epetra_Vector> p_init;
    Teuchos::RCP<const Epetra_Map> g_map;
    Teuchos::RCP<NOX::Epetra::Vector> currentSolution;
    Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface;
    Teuchos::RCP<LOCA::ParameterVector> pVector;
    Teuchos::RCP<NOX::Epetra::Group> grp;
    mutable Teuchos::RCP<const Epetra_Vector> finalSolution;
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
    mutable NOX::StatusTest::StatusType status;
  };

}
#endif
