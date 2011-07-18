// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_SolverCore_def.hpp
  \class  Amesos2::SolverCore
  \author Eric T Bavier <etbavie@sandia.gov>
  \date   Thu May 27 14:02:35 CDT 2010

  \brief  Templated core-functionality class for Amesos2 solvers.  Definitions.
*/

#ifndef AMESOS2_SOLVERCORE_DEF_HPP
#define AMESOS2_SOLVERCORE_DEF_HPP


namespace Amesos2 {


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverCore<ConcreteSolver,Matrix,Vector>::SolverCore(
  Teuchos::RCP<Matrix> A,
  Teuchos::RCP<Vector> X,
  Teuchos::RCP<Vector> B )
  : matrixA_(createMatrixAdapter<Matrix>(A))
  , multiVecX_(X) // may be null
  , multiVecB_(B) // may be null
  , globalNumRows_(matrixA_->getGlobalNumRows())
  , globalNumCols_(matrixA_->getGlobalNumCols())
  , globalNumNonZeros_(matrixA_->getGlobalNNZ())
  , status_(matrixA_->getComm())
{
    TEST_FOR_EXCEPTION(
    !matrixShapeOK(),
    std::invalid_argument,
    "Matrix shape inappropriate for this solver");
}


/// Destructor
template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverCore<ConcreteSolver,Matrix,Vector>::~SolverCore( )
{
  // TODO: The below code is still dependent on there being some kind
  // of status and control code available to Amesos2::SolverCore,
  // either in the form of a shared library or private inheritance of
  // classes that provide that functionality.

  // print out some information if required by the user
  // if ((control_.verbose_ && control_.printTiming_) || control_.verbose_ == 2) printTiming();
  // if ((control_.verbose_ && control_.printStatus_) || control_.verbose_ == 2) printStatus();

}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::preOrdering()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  static_cast<solver_type*>(this)->preOrdering_impl();
  ++status_.numPreOrder_;
  status_.preOrderingDone_ = true;
  status_.symbolicFactorizationDone_ = false;
  status_.numericFactorizationDone_  = false;

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::symbolicFactorization()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  if( !preOrderingDone() ){
    preOrdering();
  }

  static_cast<solver_type*>(this)->symbolicFactorization_impl();
  ++status_.numSymbolicFact_;
  status_.symbolicFactorizationDone_ = true;
  status_.numericFactorizationDone_  = false;

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::numericFactorization()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  if( !symbolicFactorizationDone() ){
    symbolicFactorization();
  }

  static_cast<solver_type*>(this)->numericFactorization_impl();
  ++status_.numNumericFact_;
  status_.numericFactorizationDone_ = true;

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::solve()
{
  solve(multiVecX_.ptr(), multiVecB_.ptr());
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::solve(const Teuchos::Ptr<Vector> X,
                                                const Teuchos::Ptr<const Vector> B) const
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  X.assert_not_null();
  B.assert_not_null();

  const Teuchos::RCP<MultiVecAdapter<Vector> > x =
    createMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(X));
  const Teuchos::RCP<const MultiVecAdapter<Vector> > b =
    createConstMultiVecAdapter<Vector>(Teuchos::rcpFromPtr(B));
    
#ifdef HAVE_AMESOS2_DEBUG
  // Check some required properties of X and B
  TEST_FOR_EXCEPTION(x->getGlobalLength() != matrixA_->getGlobalNumCols(),
                     std::invalid_argument,
                     "MultiVector X must have length equal to the number of "
                     "global columns in A");

  TEST_FOR_EXCEPTION(b->getGlobalLength() != matrixA_->getGlobalNumRows(),
                     std::invalid_argument,
                     "MultiVector B must have length equal to the number of "
                     "global rows in A");
  
  TEST_FOR_EXCEPTION(x->getGlobalNumVectors() != b->getGlobalNumVectors(),
                     std::invalid_argument,
                     "X and B MultiVectors must have the same number of vectors");
#endif  // HAVE_AMESOS2_DEBUG
  
  if( !numericFactorizationDone() ){
    // This casting-away of constness is probably OK because this
    // function is meant to be "logically const"
    const_cast<type*>(this)->numericFactorization();
  }
  
  static_cast<const solver_type*>(this)->solve_impl(Teuchos::outArg(*x), Teuchos::ptrInArg(*b));
  ++status_.numSolve_;
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::solve(Vector* X, const Vector* B) const
{
  solve(Teuchos::ptr(X), Teuchos::ptr(B));
}

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
bool
SolverCore<ConcreteSolver,Matrix,Vector>::matrixShapeOK()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  return( static_cast<solver_type*>(this)->matrixShapeOK_impl() );
}

  // The RCP should probably be to a const Matrix, but that throws a
  // wrench in some of the traits mechanisms that aren't expecting
  // const types
template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::setA( const Teuchos::RCP<Matrix> a )
{
  matrixA_ = createMatrixAdapter(a);
  if( !matrixA_.is_null() ){
    // The user passed a non-null rcp, reset state
    //
    // Resetting the done-ness of each phase will force each to be
    // repeated from scratch.  The solver interface is responsible
    // for assuring that no dependencies linger between one matrix
    // and another.
    
    status_.preOrderingDone_ = false;
    status_.symbolicFactorizationDone_ = false;
    status_.numericFactorizationDone_  = false;

    globalNumNonZeros_ = matrixA_->getGlobalNNZ();
    globalNumCols_     = matrixA_->getGlobalNumCols();
    globalNumRows_     = matrixA_->getGlobalNumRows();
  }
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<Matrix,Vector>&
SolverCore<ConcreteSolver,Matrix,Vector>::setParameters(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1(timers_.totalTime_);
#endif

  if( parameterList->name() == "Amesos2" ){
    // Do everything here that is for generic status and control parameters
    control_.setControlParameters(parameterList);
    status_.setStatusParameters(parameterList);
    
    // Finally, hook to the implementation's parameter list parser
    // First check if there is a dedicated sublist for this solver and use that if there is
    if( parameterList->isSublist(name()) ){
      static_cast<solver_type*>(this)->setParameters_impl(Teuchos::sublist(parameterList, name()));
    } else {
      // See if there is anything in the parent sublist that the solver recognizes
      static_cast<solver_type*>(this)->setParameters_impl(parameterList);
    }
  }

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Teuchos::RCP<const Teuchos::ParameterList>
SolverCore<ConcreteSolver,Matrix,Vector>::getValidParameters() const
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor LocalTimer1( timers_.totalTime_ );
#endif

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<ParameterList> control_params = rcp(new ParameterList("Amesos2 Control"));
  control_params->set("Transpose",false);
  control_params->set("AddToDiag","");
  control_params->set("AddZeroToDiag",false);
  control_params->set("MatrixProperty","general");
  control_params->set("ScaleMethod",0);

  Teuchos::RCP<ParameterList> status_params = rcp(new ParameterList("Amesos2 Status"));
  status_params->set("PrintTiming",false);
  status_params->set("PrintStatus",false);
  status_params->set("ComputeVectorNorms",false);
  status_params->set("ComputeTrueResidual",false);

  RCP<const ParameterList>
    solver_params = static_cast<const solver_type*>(this)->getValidParameters_impl();

  RCP<ParameterList> amesos2_params = rcp(new ParameterList("Amesos2"));
  amesos2_params->setParameters(*control_params);
  amesos2_params->setParameters(*status_params);
  amesos2_params->set(name(), *solver_params);

  return amesos2_params;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
std::string
SolverCore<ConcreteSolver,Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << name() << " solver interface";
  return oss.str();
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
{
  if( matrixA_.is_null() || (status_.myPID_ != 0) ){ return; }
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getComm();
  size_t width = 1;
  for( size_t dec = 10; dec < globalNumRows_; dec *= 10 ) {
    ++width;
  }
  width = std::max<size_t>(width,11) + 2;
  Teuchos::OSTab tab(out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium: print O(P) info, num entries per node
  //    high: print O(N) info, num entries per row
  // extreme: print O(NNZ) info: print indices and values
  //
  // for medium and higher, print constituent objects at specified verbLevel
  if( vl != VERB_NONE ) {
    std::string p = name();
    Util::printLine(out);
    out << this->description() << std::endl << std::endl;

    out << p << "Matrix has " << globalNumRows_ << "rows"
        << " and " << globalNumNonZeros_ << "nonzeros"
        << std::endl;
    if( vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME ){
      out << p << "Nonzero elements per row = "
          << globalNumNonZeros_ / globalNumRows_
          << std::endl;
      out << p << "Percentage of nonzero elements = "
          << 100.0 * globalNumNonZeros_ / (globalNumRows_ * globalNumCols_)
          << std::endl;
    }
    if( vl == VERB_HIGH || vl == VERB_EXTREME ){
      out << p << "Use transpose = " << control_.useTranspose_
          << std::endl;
    }
    if ( vl == VERB_EXTREME ){
      printTiming(out,vl);
    }
    Util::printLine(out);
  }
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::printTiming(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  if( matrixA_.is_null() || (status_.myPID_ != 0) ){ return; }

  double preTime  = timers_.preOrderTime_.totalElapsedTime();
  double symTime  = timers_.symFactTime_.totalElapsedTime();
  double numTime  = timers_.numFactTime_.totalElapsedTime();
  double solTime  = timers_.solveTime_.totalElapsedTime();
  double totTime  = timers_.totalTime_.totalElapsedTime();
  double overhead = totTime - (preTime + symTime + numTime + solTime);

  std::string p = name() + " : ";
  Util::printLine(out);

  out << p << "Time to convert matrix to implementation format = "
      << timers_.mtxConvTime_.totalElapsedTime() << " (s)"
      << std::endl;
  out << p << "Time to redistribute matrix = "
      << timers_.mtxRedistTime_.totalElapsedTime() << " (s)"
      << std::endl;

  out << p << "Time to convert vectors to implementation format = "
      << timers_.vecConvTime_.totalElapsedTime() << " (s)"
      << std::endl;
  out << p << "Time to redistribute vectors = "
      << timers_.vecRedistTime_.totalElapsedTime() << " (s)"
      << std::endl;

  out << p << "Number of pre-orderings = "
      << getNumPreOrder()
      << std::endl;
  out << p << "Time for pre-ordering = "
      << preTime << " (s), avg = "
      << preTime / getNumPreOrder() << " (s)"
      << std::endl;

  out << p << "Number of symbolic factorizations = "
      << getNumSymbolicFact()
      << std::endl;
  out << p << "Time for sym fact = "
      << symTime << " (s), avg = "
      << symTime / getNumSymbolicFact() << " (s)"
      << std::endl;

  out << p << "Number of numeric factorizations = "
      << getNumNumericFact()
      << std::endl;
  out << p << "Time for num fact = "
      << numTime << " (s), avg = "
      << numTime / getNumNumericFact() << " (s)"
      << std::endl;

  out << p << "Number of solve phases = "
      << getNumSolve()
      << std::endl;
  out << p << "Time for solve = "
      << solTime << " (s), avg = "
      << solTime / getNumSolve() << " (s)"
      << std::endl;

  out << p << "Total time spent in Amesos2 = "
      << totTime << " (s)"
      << std::endl;
  out << p << "Total time spent in the Amesos2 interface = "
      << overhead << " (s)"
      << std::endl;
  out << p << "  (the above time does not include solver time)"
      << std::endl;
  out << p << "Amesos2 interface time / total time = "
      << overhead / totTime
      << std::endl;
  Util::printLine(out);
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
SolverCore<ConcreteSolver,Matrix,Vector>::getTiming(
  Teuchos::ParameterList& timingParameterList) const
{}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
std::string
SolverCore<ConcreteSolver,Matrix,Vector>::name() const {
  std::string solverName = solver_type::name;
  return solverName;
}


} // end namespace Amesos2

#endif  // AMESOS2_SOLVERCORE_DEF_HPP
