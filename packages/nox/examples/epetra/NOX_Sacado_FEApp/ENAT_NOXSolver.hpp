#ifndef ENAT_NOXSOLVER_H
#define ENAT_NOXSOLVER_H

#include "NOX.H"
#include "NOX_Epetra.H"
#include "EpetraExt_ModelEvaluator.h"

/** \brief Epetra-based Model Evaluator subclass for Charon!
 *
 * This class will support a wide number of different types of abstract
 * problem types that will allow NOX, LOCA, Rythmos, Aristos, and MOOCHO to
 * solve different types of problems with Charon.
 * 
 * ToDo: Finish documentation!
 */

namespace ENAT {
class NOXSolver
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
            Teuchos::RCP<EpetraExt::ModelEvaluator> model,
            Teuchos::RCP<Epetra_Operator> M = Teuchos::null
            );


  //@}

  ~NOXSolver();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
//  Teuchos::RCP<Epetra_Operator> create_W() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::InArgs createInArgs() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  private:
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \brief . */
  void setProblemParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  void setSolverParamDefaults(Teuchos::ParameterList* appParams_);

  //@}

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;

   Teuchos::RCP<NOX::Solver::Generic> solver;
   Teuchos::RCP<NOX::Epetra::Vector> currentSolution;
   Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface;
   int num_p;
   int num_g;

   Teuchos::RCP<NOX::Epetra::Group> grp;
};

}
#endif
