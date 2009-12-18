#ifndef PIRO_EPETRA_RYTHMOSSOLVER_H
#define PIRO_EPETRA_RYTHMOSSOLVER_H

#include <iostream>

#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_ModelEvaluator.h"

#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"


/** \brief Epetra-based Model Evaluator subclass for Charon!
 *
 * This class will support a wide number of different types of abstract
 * problem types that will allow NOX, LOCA, Rythmos, Aristos, and MOOCHO to
 * solve different types of problems with Charon.
 * 
 * ToDo: Finish documentation!
 */

namespace Piro {
namespace Epetra {

class RythmosSolver
    : public EpetraExt::ModelEvaluator
{

  typedef double Scalar;

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  RythmosSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
                Teuchos::RCP<EpetraExt::ModelEvaluator> model,
                Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer = Teuchos::null
                );

  //@}

  ~RythmosSolver();


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
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidRythmosParameters() const;

  //@}

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;

   int num_p;
   int num_g;

   Teuchos::RCP<Rythmos::StepperBase<Scalar> > fwdStateStepper;
   Teuchos::RCP<Teuchos::FancyOStream> out;
   Teuchos::EVerbosityLevel solnVerbLevel;
   Teuchos::RCP<Rythmos::DefaultIntegrator<Scalar> > fwdStateIntegrator;
   Teuchos::RCP<Thyra::ModelEvaluator<double> > fwdStateModel;
   Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<double> > fwdTimeStepSolver;
   Scalar t_final;

};

}
}
#endif
