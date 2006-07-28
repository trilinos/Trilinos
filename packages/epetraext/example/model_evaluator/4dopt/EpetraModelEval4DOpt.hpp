#ifndef EPETRA_MODEL_EVAL_4D_OPT_HPP
#define EPETRA_MODEL_EVAL_4D_OPT_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

/** \brief An upgrade of <tt>EpetraModelEval2DSim</tt> that add the ability to
 * manipulate the parameters, adds a response function and includes first
 * derivatives.
 *
 * ToDo: Finish Documentation!
 */
class EpetraModelEval4DOpt : public EpetraExt::ModelEvaluator {
public:

  // Constructor
  EpetraModelEval4DOpt(
		const double         xt0         = 1.0
		,const double        xt1         = 1.0
		,const double        pt0         = 2.0
		,const double        pt1         = 0.0
		,const double        d           = 10.0
		,const double        x00         = 1.0
		,const double        x01         = 1.0
		,const double        p00         = 2.0
		,const double        p01         = 0.0
    );

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_lower_bounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_upper_bounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_lower_bounds(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_upper_bounds(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
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

	double    xt0_;
	double    xt1_;
	double    pt0_;
	double    pt1_;
  double    d_;

	bool      isInitialized_;

  Teuchos::RefCountPtr<const Epetra_Comm>  epetra_comm_;
	Teuchos::RefCountPtr<const Epetra_Map>   map_x_;
	Teuchos::RefCountPtr<const Epetra_Map>   map_p_;
	Teuchos::RefCountPtr<const Epetra_Map>   map_g_;

	Teuchos::RefCountPtr<Epetra_Vector> xL_;
	Teuchos::RefCountPtr<Epetra_Vector> xU_;
	Teuchos::RefCountPtr<Epetra_Vector> pL_;
	Teuchos::RefCountPtr<Epetra_Vector> pU_;
	Teuchos::RefCountPtr<Epetra_Vector> gL_;
	Teuchos::RefCountPtr<Epetra_Vector> gU_;
	Teuchos::RefCountPtr<Epetra_Vector> x0_;
	Teuchos::RefCountPtr<Epetra_Vector> p0_;

  Teuchos::RefCountPtr<Epetra_CrsGraph>  W_graph_;

};

#endif // EPETRA_MODEL_EVAL_4D_OPT_HPP
