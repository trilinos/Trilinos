// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MASSSPRINGDAMPERMODEL_H
#define MASSSPRINGDAMPERMODEL_H

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "MatrixBased_LOWS.hpp"


using LO = Tpetra::Map<>::local_ordinal_type;
using GO = Tpetra::Map<>::global_ordinal_type;
typedef Tpetra::Map<LO,GO>  Tpetra_Map;
typedef Tpetra::Vector<double,LO,GO>  Tpetra_Vector;
typedef Tpetra::MultiVector<double,LO,GO>  Tpetra_MultiVector;
typedef Tpetra::Operator<double,LO,GO>  Tpetra_Operator;
typedef Tpetra::CrsGraph<LO,GO>  Tpetra_CrsGraph;
typedef Tpetra::CrsMatrix<double,LO,GO>  Tpetra_CrsMatrix;
typedef Thyra::TpetraOperatorVectorExtraction<
    double, LO, GO> ConverterT;

/** \brief Mass-spring-damper model problem for transient PDE constrained optimization.
 * 
 * This problem is subject to the constraint of a mass-spring-damper differential equation
 *   \f[
 *   m \ddot{x} + 2 \sqrt{m\,k} \dot{x} + k x - F = 0
 *   \f]
 * where \f$k\f$ and \f$m\f$ are respectivelly the first and second parameter.
 * 
 * Those parameters have initial values set to \f$1.\f$ for both of them and there 
 * are bound constraints that prevent them from leaving the \f$\left[0.5, 1.5\right]\f$ range 
 * for both of them.
 * 
 * The objective function that we try to minimize is:
 *   \f[
 *   s_{g_x} * (( x - t_x )^2 + s ( \dot{x} - t_{\dot{x}} )^2) + s_{g_p} * (( k - t_k )^2 + ( m - t_m )^2)
 *   \f]
 * where \f$s_{g_x}\f$, \f$s\f$, and \f$s_{g_p}\f$ are scaling factors and 
 * where \f$t_x\f$, \f$t_{\dot{x}}\f$, \f$t_k\f$, and \f$t_m\f$ are targeted
 * values for \f$x\f$, \f$\dot{x}\f$, \f$k\f$, and \f$m\f$ respecivelly.
 * 
 * The second order in time equation is solved by writing the residual 
 * using \f$\boldsymbol{u}=\left[ x, \dot{x}\right]\f$.
 * 
 * Therefore, the objective function can either uses \f$u_1\f$ or \f$\dot{u}_0\f$ to 
 * represent \f$\dot{x}\f$. Both of these options can be used and can be selected by choosing
 * use_x_dot_in_g = false or true respectively during the construction of the problem class.
*/

class MassSpringDamperModel
    : public Thyra::ModelEvaluatorDefaultBase<double>
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MassSpringDamperModel(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, bool adjoint = false, const Teuchos::RCP<Teuchos::ParameterList>& problemList = Teuchos::null, bool hessianSupport = false, bool use_x_dot_in_g = false);

  //@}

  ~MassSpringDamperModel();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double> getLowerBounds() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double> getUpperBounds() const;

  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<double>>
  create_W_op() const;

  /** \brief . */
  Teuchos::RCP<Thyra::PreconditionerBase<double>>
  create_W_prec() const;

  /** \brief . */
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
  get_W_factory() const;

  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<double>>
  create_hess_g_pp( int j, int l1, int l2 ) const;

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double>
  createInArgs() const;

  /** \brief . */
  void
  reportFinalPoint(
      const Thyra::ModelEvaluatorBase::InArgs<double>& finalPoint,
      const bool wasSolved);

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>  get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>  get_f_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> get_p_space(int l) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> get_g_space(int j) const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  Teuchos::ArrayView<const std::string> get_g_names(int j) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "not implemented");
  }
  //@}

  protected:

  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<double>
  createOutArgsImpl() const;

  /** \brief . */
  void
  evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const;
  //@}


  private:

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double>
  createInArgsImpl() const;

   //These are set in the constructor and used in evalModel
  Teuchos::RCP<const Tpetra_Map> x_map;
  Teuchos::RCP<const Tpetra_Map> p_map;
  Teuchos::RCP<const Tpetra_Map> g_map;
  Teuchos::RCP<Tpetra_CrsGraph> crs_graph;
  Teuchos::RCP<Tpetra_CrsGraph> hess_crs_graph;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  Teuchos::RCP<Tpetra_Vector> p_vec;
  Teuchos::RCP<Tpetra_Vector> x_vec;
  Teuchos::RCP<Tpetra_Vector> x_dot_vec;

  //! Cached nominal values and lower/upper bounds
  Thyra::ModelEvaluatorBase::InArgs<double> nominalValues;
  Thyra::ModelEvaluatorBase::InArgs<double> lowerBounds;
  Thyra::ModelEvaluatorBase::InArgs<double> upperBounds;

  //whether hessian is supported 
  bool hessSupport, adjoint_, use_x_dot_in_g_;

  //Problem parameter list
  Teuchos::RCP<Teuchos::ParameterList> probList_;

  double target_x_, target_x_dot_, target_k_, target_m_;
  double scaling_, scaling_g_x_, scaling_g_p_;
};

#endif // MASSSPRINGDAMPERMODEL_H
