// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MOCKMODELEVAL_C_TPETRA_H
#define MOCKMODELEVAL_C_TPETRA_H

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Piro_TransientDecorator.hpp"
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

/** \brief Concrete Tpetra-based Model Evaluator
 *
 * Concrete model evaluator for the solution of the following PDE-Constrained problem:
 *
 * solve
 * u_tt + p = 0
 * p = 1
 * 
 * g=u
 *
 * solution is u(t) = 0.5*t*(2-t).
 */

class MockModelEval_C_Tpetra
    : public Piro::TransientDecorator<double>
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MockModelEval_C_Tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, bool adjoint = false, const Teuchos::RCP<Teuchos::ParameterList>& problemList = Teuchos::null, bool hessianSupport = false);

  //@}

  ~MockModelEval_C_Tpetra();


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
  Teuchos::RCP<Tpetra_Vector> x_dotdot_vec;

   //! Cached nominal values and lower/upper bounds
   Thyra::ModelEvaluatorBase::InArgs<double> nominalValues;
   Thyra::ModelEvaluatorBase::InArgs<double> lowerBounds;
   Thyra::ModelEvaluatorBase::InArgs<double> upperBounds;

   //whether hessian is supported 
   bool hessSupport;

   //Problem parameter list
   Teuchos::RCP<Teuchos::ParameterList> probList_;

};

#endif // SIMPLE_MODELEVAL_H
