#ifndef PIRO_EPETRA_STOKHOSNOXSOLVER_H
#define PIRO_EPETRA_STOKHOSNOXSOLVER_H

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
#include "Piro_Epetra_StokhosNOXObserver.hpp"

#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Quadrature.hpp"

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
  class StokhosNOXSolver : public EpetraExt::ModelEvaluator {
  public:

    /** \name Constructors/initializers */
    //@{

    /** \brief Takes the number of elements in the discretization . */
    StokhosNOXSolver(const Teuchos::RCP<Teuchos::ParameterList>& appParams,
		const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
 //		const Teuchos::RCP<FEApp::Application>& app,
		const Teuchos::RCP<const Epetra_Comm>& comm,
                Teuchos::RCP<NOXObserver> noxObserver);


    //@}

    ~StokhosNOXSolver();


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

    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> >
    getBasis() const { return basis; }
    Teuchos::RCP<const Stokhos::Quadrature<int,double> >
    getQuad() const { return quad; }
    
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

    Teuchos::RCP<const Teuchos::ParameterList>
     getValidSGParameters() const;

    //@}
    
  private:
    
    enum SG_METHOD {
      SG_AD,
      SG_GLOBAL,
      SG_NI
    };

    //These are set in the constructor and used in evalModel
    Teuchos::RCP<EpetraExt::ModelEvaluator> sg_solver;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
  };
  
}
}
#endif
