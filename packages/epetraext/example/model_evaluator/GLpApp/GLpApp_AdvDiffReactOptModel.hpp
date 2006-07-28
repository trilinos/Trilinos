#ifndef GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP
#define GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "GLpApp_GLpYUEpetraDataPool.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Array.hpp"

namespace GLpApp {

/** \brief 
 *
 * ToDo: Finish Documentation!
 */
class AdvDiffReactOptModel
  : public EpetraExt::ModelEvaluator
  , public Teuchos::VerboseObject<AdvDiffReactOptModel>
{
public:

  /** \brief Constructor. */
  AdvDiffReactOptModel(
    const Teuchos::RefCountPtr<const Epetra_Comm>  &comm
    ,const double                                  beta
    ,const double                                  len_x     // Ignored if meshFile is *not* empty
    ,const double                                  len_y     // Ignored if meshFile is *not* empty
    ,const int                                     local_nx  // Ignored if meshFile is *not* empty
    ,const int                                     local_ny  // Ignored if meshFile is *not* empty
    ,const char                                    meshFile[]
    ,const int                                     np
    ,const double                                  x0
    ,const double                                  p0
    ,const double                                  reactionRate
    ,const bool                                    normalizeBasis
    );

  /** \brief . */
  void set_q( Teuchos::RefCountPtr<const Epetra_Vector> const& q );

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
  Teuchos::RefCountPtr<Epetra_Operator> create_DfDp_op(int l) const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private types

  typedef Teuchos::Array<Teuchos::RefCountPtr<const Epetra_Map> >  RCP_Eptra_Map_Array_t;
  typedef Teuchos::Array<Teuchos::RefCountPtr<Epetra_Vector> >     RCP_Eptra_Vector_Array_t;

  // /////////////////////////////////////
  // Private member data

  static const int Np_         = 2; // Number of axiliary parameters
  static const int p_bndy_idx  = 0; // index for boundary flux parameters
  static const int p_rx_idx    = 1; // index for reaction rate parameter

  bool      isInitialized_;

  Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool>   dat_;
  int                                                 np_;
  Teuchos::RefCountPtr<const Epetra_Vector>           q_;

  Teuchos::RefCountPtr<const Epetra_Map>              map_p_bar_;
  Teuchos::RefCountPtr<Epetra_MultiVector>            B_bar_;

  Teuchos::RefCountPtr<const Epetra_Comm>  epetra_comm_;
  Teuchos::RefCountPtr<const Epetra_Map>   map_x_;
  RCP_Eptra_Map_Array_t                    map_p_;
  Teuchos::RefCountPtr<const Epetra_Map>   map_f_;
  Teuchos::RefCountPtr<const Epetra_Map>   map_g_;

  Teuchos::RefCountPtr<Epetra_Vector> x0_;
  Teuchos::RefCountPtr<Epetra_Vector> xL_;
  Teuchos::RefCountPtr<Epetra_Vector> xU_;
  RCP_Eptra_Vector_Array_t            p0_;
  RCP_Eptra_Vector_Array_t            pL_;
  RCP_Eptra_Vector_Array_t            pU_;
  Teuchos::RefCountPtr<Epetra_Vector> gL_;
  Teuchos::RefCountPtr<Epetra_Vector> gU_;

  Teuchos::RefCountPtr<Epetra_CrsGraph>  W_graph_;

};

} // namespace GLpApp

#endif // GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP
