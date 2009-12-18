#ifndef PIRO_EPETRA_LOCASOLVER_H
#define PIRO_EPETRA_LOCASOLVER_H

#include <iostream>

#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"
#include "LOCA_Epetra_ModelEvaluatorInterface.H"
#include <NOX_Epetra_MultiVector.H>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_ModelEvaluator.h"
#include "Piro_Epetra_NOXObserver.hpp"

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

class LOCASolver
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  LOCASolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
            Teuchos::RCP<EpetraExt::ModelEvaluator> model,
            Teuchos::RCP<NOXObserver> observer = Teuchos::null
            );


  //@}

  ~LOCASolver();


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
   Teuchos::RCP<NOXObserver> observer;
   NOX::Utils utils;

   Teuchos::RCP<LOCA::Stepper> stepper;
   Teuchos::RCP<NOX::Epetra::Vector> currentSolution;
   Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> interface;
   Teuchos::RCP<LOCA::ParameterVector> pVector;
   int num_p;
   int num_g;

   Teuchos::RCP<LOCA::Epetra::Group> grp;
};

}
}
#endif
