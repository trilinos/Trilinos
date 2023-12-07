// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef MOCKMODELEVAL_H_TPETRA_H
#define MOCKMODELEVAL_H_TPETRA_H

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

/** \brief Concrete Tpetra-based Model Evaluator
 *
 * Concrete model evaluator for the solution of the following PDE-Constrained problem:
 *
 * find p that minimizes
 * g0 = 0.5 (x - T)' (x-T),  T = (1+X)^3
 * subject to:
 * f = x - p^3
 * 
 * p_opt = (1+X)
 * x_opt = (1+X)^3
 * 
 * f_true  = x - p^3 + 0.2 p^2
 *
 * p_sample_0 = p_opt
 * p_sample_1 = X + X^2
 */


  class DfDpOp : public Thyra::LinearOpBase<double> {
      public:

    // Constructor
    DfDpOp(
      const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space,
      const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space);

    //! Destructor
    virtual ~DfDpOp() = default;

    inline void set(const Teuchos::RCP<const Tpetra_Vector> p_vec) {p_vec_ = p_vec;}

    //! Overrides Thyra::LinearOpBase purely virtual method
    inline Teuchos::RCP<const Thyra::VectorSpaceBase<double>> domain() const {
      return p_space_;
    }

    //! Overrides Thyra::LinearOpBase purely virtual method
    inline Teuchos::RCP<const Thyra::VectorSpaceBase<double>> range() const {
      return x_space_;
    }

    //@}

  protected:
    //! Overrides Thyra::LinearOpBase purely virtual method
    inline bool opSupportedImpl(Thyra::EOpTransp /*M_trans*/) const {
      // The underlying scalar type is not complex, and we support transpose, so we support everything.
      return true;
    }

    //! Overrides Thyra::LinearOpBase purely virtual method
    void applyImpl (const Thyra::EOpTransp M_trans,
                    const Thyra::MultiVectorBase<double>& X,
                    const Teuchos::Ptr<Thyra::MultiVectorBase<double>>& Y,
                    const double /* alpha */,
                    const double /* beta */) const; 

    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space_;
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space_;
    Teuchos::RCP<const Tpetra_Vector> p_vec_;

  };


class MockModelEval_H_Tpetra
    : public Thyra::ModelEvaluatorDefaultBase<double>
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MockModelEval_H_Tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >  appComm, bool adjoint = false, const Teuchos::RCP<Teuchos::ParameterList>& problemList = Teuchos::null, bool hessianSupport = false);

  //@}

  ~MockModelEval_H_Tpetra();




  /** \brief . */
  Teuchos::RCP<Thyra::VectorBase<double>> get_true_p_opt(int j) const;


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double> getLowerBounds() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<double> getUpperBounds() const;

  /** \brief . */
  Teuchos::RCP<Thyra::VectorBase<double>> get_p_opt(int j) const;
  /** \brief . */
  Teuchos::RCP<Thyra::VectorBase<double>> get_param_samples(int k) const;
  /** \brief . */
  Teuchos::RCP<Thyra::VectorBase<double>> get_solution_diff_at_samples(int k) const;

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
  Teuchos::RCP<Thyra::LinearOpBase<double>>
  create_DfDp_op_impl(int j) const;

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
  Teuchos::RCP<Tpetra_CrsGraph> hess_crs_graph_p;
  Teuchos::RCP<Tpetra_CrsGraph> hess_crs_graph_x;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  Teuchos::RCP<Tpetra_Vector> p_vec_0;
  Teuchos::RCP<Tpetra_Vector> p_vec_1;
  Teuchos::RCP<Tpetra_Vector> x_vec;
  Teuchos::RCP<Tpetra_Vector> x_dot_vec;
  Teuchos::RCP<Tpetra_Vector> coords_vec;

  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> p_space;
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> x_space;

   //! Cached nominal values and lower/upper bounds
   Thyra::ModelEvaluatorBase::InArgs<double> nominalValues;
   Thyra::ModelEvaluatorBase::InArgs<double> lowerBounds;
   Thyra::ModelEvaluatorBase::InArgs<double> upperBounds;

   //whether hessian is supported 
   bool hessSupport;

   //discretization step
   double h;

   //Problem parameter list
   Teuchos::RCP<Teuchos::ParameterList> probList_;

};

#endif // SIMPLE_MODELEVAL_H
