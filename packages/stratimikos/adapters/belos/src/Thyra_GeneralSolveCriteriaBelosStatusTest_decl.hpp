/*
// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#ifndef THYRA_GENERAL_SOLVE_CRITERIA_BELOS_STATUS_TEST_DECL_HPP
#define THYRA_GENERAL_SOLVE_CRITERIA_BELOS_STATUS_TEST_DECL_HPP

#include "Thyra_SolveSupportTypes.hpp"
#include "BelosStatusTest.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Thyra {


/** \brief Subclass of Belos::StatusTest that implements every possible form
 * of SolveCriteria that exists by forcing the computation of <tt>X</tt> and
 * <tt>R</tt>.
 */
template<class Scalar>
class GeneralSolveCriteriaBelosStatusTest
  : public Belos::StatusTest<Scalar, MultiVectorBase<Scalar>, LinearOpBase<Scalar> >,
    public Teuchos::VerboseObject<GeneralSolveCriteriaBelosStatusTest<Scalar> >
{
public:

  /** \name Public typdefs. */
  //@{
  /** \brief . */
  typedef MultiVectorBase<Scalar> MV;
  /** \brief . */
  typedef LinearOpBase<Scalar> OP;
  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  //@}

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  GeneralSolveCriteriaBelosStatusTest();

  /** \brief . */
  void setSolveCriteria(const SolveCriteria<Scalar> &solveCriteria,
    const int convergenceTestFrequency);

  /** \brief . */
  ArrayView<const ScalarMag> achievedTol() const;

  //@}
  
  /** \name Overridden public functions from Belos::StatusTest. */
  //@{ 
  /** \brief . */
  virtual Belos::StatusType checkStatus(Belos::Iteration<Scalar,MV,OP> *iSolver);
  /** \brief . */
  virtual Belos::StatusType getStatus() const;
  /** \brief . */
  virtual void reset();
  /** \brief . */
  virtual void print(std::ostream& os, int indent) const;
  //@}

private:

  SolveCriteria<Scalar> solveCriteria_;
  int convergenceTestFrequency_;

  bool compute_x_;
  bool compute_r_;

  Array<ScalarMag> r0_nrm_;
  Array<ScalarMag> b_nrm_;
  Array<ScalarMag> lastNumerator_;
  Array<ScalarMag> lastDenominator_;
  Array<ScalarMag> lastAchievedTol_;
  int lastCurrIter_;
  Belos::StatusType lastRtnStatus_;

  // Private member functions

  ScalarMag computeReductionFunctional(ESolveMeasureNormType measureType,
    const Ptr<const ReductionFunctional<Scalar> > &reductFunc,
    const Ptr<const VectorBase<Scalar> > &x,
    const Ptr<const VectorBase<Scalar> > &r
    ) const;

  void printRhsStatus(const int currIter, const int j, std::ostream &out,
    int indent = 0) const;

};


/** \brief Nonmember constructor.
 *
 * \relates GeneralSolveCriteriaBelosStatusTest
 */
template<class Scalar>
RCP<GeneralSolveCriteriaBelosStatusTest<Scalar> >
createGeneralSolveCriteriaBelosStatusTest(
  const SolveCriteria<Scalar> &solveCriteria,
  const int convergenceTestFrequency
  )
{
  RCP<GeneralSolveCriteriaBelosStatusTest<Scalar> >
    gscbst = Teuchos::rcp(new GeneralSolveCriteriaBelosStatusTest<Scalar>);
  gscbst->setSolveCriteria(solveCriteria, convergenceTestFrequency);
  return gscbst;
}


} // namespace Thyra


#endif	// THYRA_GENERAL_SOLVE_CRITERIA_BELOS_STATUS_TEST_DECL_HPP
