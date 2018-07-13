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

#ifndef MOCKMODELEVAL_A_TPETRA_H
#define MOCKMODELEVAL_A_TPETRA_H

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"


typedef Tpetra::Map<int,int>  Tpetra_Map;
typedef Tpetra::Vector<double,int,int>  Tpetra_Vector;
typedef Tpetra::MultiVector<double,int,int>  Tpetra_MultiVector;
typedef Tpetra::Operator<double,int,int>  Tpetra_Operator;
typedef Tpetra::CrsGraph<int,int>  Tpetra_CrsGraph;
typedef Tpetra::CrsMatrix<double,int,int>  Tpetra_CrsMatrix;
typedef Thyra::TpetraOperatorVectorExtraction<
    double, int, int> ConverterT;

/** \brief Concrete Tpetra-based Model Evaluator
 *
 * Concrete model evaluator for the solution of the following PDE-Constrained problem:
 *
 * find (p_0,p_1) that minimizes
 * g = 0.5*(Sum(x)-Sum(p)-12)^2 + 0.5*(p0-1)^2
 * subject to:
 * f_0 = (x_0)^2 - p_0 = 0
 * f_i = x_i^2 - (i+p_1)^2 (for i != 0), for i = 1,2,3,4
 *
 * solution is p = (1,3).
 */


class MockModelEval_A_Tpetra
    : public Thyra::ModelEvaluatorDefaultBase<double>
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MockModelEval_A_Tpetra(const MPI_Comm appComm);

  //@}

  ~MockModelEval_A_Tpetra();


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
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "not impl'ed");
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
  Teuchos::RCP<Teuchos::Comm<int> > comm;

  Teuchos::RCP<Tpetra_Vector> p_vec;
  Teuchos::RCP<Tpetra_Vector> x_vec;
  Teuchos::RCP<Tpetra_Vector> x_dot_vec;

   //! Cached nominal values and lower/upper bounds
   Thyra::ModelEvaluatorBase::InArgs<double> nominalValues;
   Thyra::ModelEvaluatorBase::InArgs<double> lowerBounds;
   Thyra::ModelEvaluatorBase::InArgs<double> upperBounds;

};

#endif // SIMPLE_MODELEVAL_H
