// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef BELOS_STATUS_TEST_OUTPUTTER_HPP
#define BELOS_STATUS_TEST_OUTPUTTER_HPP

/*!
  \file BelosStatusTestOutputter.hpp
  \brief Belos::StatusTest for affecting printing of Belos::StatusTest objects.
*/

#include "BelosStatusTestResNorm.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace Belos {

/** \brief StatusTest implementation that does nothing but outputting.
 *
 * ToDo: Finish documentation.
 */
template <class ScalarType, class MV, class OP>
class StatusTestOutputter: public StatusTest<ScalarType,MV,OP> {
public:

  /** \name Public types */
  //@{
  
  /** \brief . */
  typedef Teuchos::ScalarTraits<ScalarType> ST;
  /** \brief . */
  typedef typename ST::magnitudeType MagnitudeType;
  /** \brief . */
  typedef StatusTestResNorm<ScalarType,MV,OP> StatusTestResNorm_t;
  
  /** \name Constructors/initializers/accessors */
  //@{
  
  /** \brief Set the number of iterations skipped between outputting.
   * A value <tt><=0</tt> means no outputting */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, outputFrequency );
    
  /** \brief Determine if only max residual is printed or not for active RHS. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, outputMaxResOnly );
    
  /** \brief Determine if only max residual is printed or not for active RHS. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, resString );
    
  /** \brief Set the residual norm status test object [Required] */
  STANDARD_COMPOSITION_MEMBERS( StatusTestResNorm_t, resNormStatusTest );
  // Note: Above the typedef is needed for this macro to work
    
  /** \brief Set the output manager [Required] */
  STANDARD_COMPOSITION_MEMBERS( OutputManager<ScalarType>, outputManager );
    
  /** \brief . */
  StatusTestOutputter(
    const int                        outputFrequency    = -1
    ,const bool                      outputMaxResOnly   = true
    ,const std::string               &resString         = "||A*x-b||/||b||"
    );
  
  //@}

  /** \name Overridden from StatusTests */
  //@{

  /** \brief. */
  StatusType checkStatus(Iteration<ScalarType,MV,OP>* iSolver);
  /** \brief. */
  StatusType getStatus() const;
  /** \brief. */
  void reset();
  /** \brief. */
  bool residualVectorRequired() const;
  /** \brief. */
  ostream& print(ostream& os, int indent) const;

  //@}

private:

  int  callsSinceLastOutput_;
  bool iterZeroWasOutput_;
  int  lastRestart_;
  
};

// /////////////////////////
// Implementation

// Constructors/initializers/accessors

template <class ScalarType, class MV, class OP>
StatusTestOutputter<ScalarType,MV,OP>::StatusTestOutputter(
  const int                        outputFrequency
  ,const bool                      outputMaxResOnly
  ,const std::string               &resString
  )
  :outputFrequency_(outputFrequency)
  ,outputMaxResOnly_(outputMaxResOnly)
  ,resString_(resString)
  ,callsSinceLastOutput_(10000)
  ,iterZeroWasOutput_(false)
{}

// Overridden from StatusTests

template <class ScalarType, class MV, class OP>
StatusType StatusTestOutputter<ScalarType,MV,OP>::checkStatus(Iteration<ScalarType,MV,OP>* iSolver)
{
  typedef MultiVecTraits<ScalarType,MV>  MVT;
  StatusType status = resNormStatusTest_->checkStatus(iSolver);
  RefCountPtr<LinearProblem<ScalarType,MV,OP> > lp = iSolver->getProblem();
  const int currIter = iSolver->GetNumIters();
  const int currRestart = iSolver->GetNumRestarts();
  if(currIter==0 && !iterZeroWasOutput_)
    lastRestart_ = -1;
  ++callsSinceLastOutput_;
  if(
    ( outputFrequency() > 0 )
    &&
    (
      status==Passed
      ||
      callsSinceLastOutput_>=outputFrequency()
      ||
      currRestart > lastRestart_
      ||
      ( currIter==0 && !iterZeroWasOutput_ )
      )
    )
  {
    const int currRhsOffset = lp->GetRHSIndex();
    const int currNumRhs = lp->GetNumToSolve();
    TEST_FOR_EXCEPT(resNormStatusTest_->GetTestValue()==NULL);
    const std::vector<MagnitudeType> &resTestValuesVector = *resNormStatusTest_->GetTestValue();
    std::ostream &out = *outputManager_->GetOStream();
    if(status==Passed)
      out << "[Converged]";
    out << "iter="<<currIter<<", restart="<<currRestart;
    const MagnitudeType maxRelRes = *std::max_element(
      resTestValuesVector.begin()+currRhsOffset,resTestValuesVector.begin()+currRhsOffset+currNumRhs
      );
    if(currNumRhs==1) {
      out << ", "<<resString()<<"="<<maxRelRes<<"\n";
    }
    else {
      if(outputMaxResOnly()) {
        out << ", max{"<<resString()<<",i="<<currRhsOffset<<"..."<<currRhsOffset+currNumRhs-1<<"}=" << maxRelRes << "\n";
      }
      else {
        out << ", "<<resString()<<":\n";
        for(int i=0; i<currNumRhs; ++i) {
          out <<"  "<<(i+currRhsOffset)<<": "<<resTestValuesVector[i+currRhsOffset]<<"\n"; 
        }
      }
    }
    if(currIter==0) {
      iterZeroWasOutput_ = true;
    }
    else {
      iterZeroWasOutput_ = false;
    }
    lastRestart_ = currRestart;
    callsSinceLastOutput_ = 0;
  }
  return status;
}

template <class ScalarType, class MV, class OP>
StatusType StatusTestOutputter<ScalarType,MV,OP>::getStatus() const
{
  return resNormStatusTest_->getStatus();
}

template <class ScalarType, class MV, class OP>
void StatusTestOutputter<ScalarType,MV,OP>::reset()
{
  return resNormStatusTest_->reset();
  callsSinceLastOutput_ = 10000;
  iterZeroWasOutput_ = false;
}

template <class ScalarType, class MV, class OP>
bool StatusTestOutputter<ScalarType,MV,OP>::residualVectorRequired() const
{
  return resNormStatusTest_->residualVectorRequired();
}

template <class ScalarType, class MV, class OP>
ostream& StatusTestOutputter<ScalarType,MV,OP>::print(ostream& os, int indent) const
{
  return resNormStatusTest_->print(os,indent);
}

} // end namespace Belos

#endif // BELOS_STATUS_TEST_OUTPUTTER_HPP
