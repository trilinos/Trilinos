#ifndef EPETRA_MODEL_EVAL_2D_SIM_HPP
#define EPETRA_MODEL_EVAL_2D_SIM_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

/** \brief Simple example ModelEvaluator subclass for a 2x2 set of
 * parameterized nonlinear equations.
 *
 * The equations modeled are:
 \verbatim

    f[0] =       x[0]      + x[1]*x[1] - p[0];
    f[1] = d * ( x[0]*x[0] - x[1]      - p[1] );

 \endverbatim
 */
class EpetraModelEval2DSim : public EpetraExt::ModelEvaluator {
public:

  /** \brief . */
  EpetraModelEval2DSim(
		 const Teuchos::RCP<const Epetra_Comm>& comm
		,const double         d                   = 10.0
		,const double        p0                  = 2.0
		,const double        p1                  = 0.0
		,const double        x00                 = 1.0
		,const double        x01                 = 1.0
    );

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private member data

  double    d_;

  bool      isInitialized_;
  
  Teuchos::RCP<const Epetra_Comm>  epetra_comm_;
  Teuchos::RCP<const Epetra_Map>   map_x_;
  
  Teuchos::RCP<Epetra_Vector> x0_;
  Teuchos::RCP<Epetra_Vector> p_;
  
  Teuchos::RCP<Epetra_CrsGraph>  W_graph_;
  
};

#endif // EPETRA_MODEL_EVAL_2D_SIM_HPP
