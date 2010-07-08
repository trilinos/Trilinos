// @HEADER
// ***********************************************************************
//
//                Amesos: Direct Sparse Solver Package
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

/**
  \file   Amesos2_Solver_def.hpp
  \class  Amesos::Solver
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu May 27 14:02:35 CDT 2010
  
  \brief  Templated class for Amesos2 solvers.  Definition.
*/

#ifndef AMESOS2_SOLVER_DEF_HPP
#define AMESOS2_SOLVER_DEF_HPP

#include <Teuchos_TimeMonitor.hpp>

#include "Amesos2_MultiVecAdapter.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_SolverBase_decl.hpp"
#include "Amesos2_Util.hpp"

//#include "Amesos2_Superlu_decl.hpp"

namespace Amesos {


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<ConcreteSolver,Matrix,Vector>::Solver(
  Teuchos::RCP<Matrix> A,
  Teuchos::RCP<Vector> X,
  Teuchos::RCP<Vector> B )
  :	matrixA_(new MatrixAdapter<Matrix>(A))
  , multiVecX_(new MultiVecAdapter<Vector>(X))
  , multiVecB_(new MultiVecAdapter<Vector>(B))
  , globalNumNonZeros_(matrixA_->getGlobalNNZ())
{
  globalNumRows_     = matrixA_->getGlobalNumRows();
  globalNumCols_     = matrixA_->getGlobalNumCols();
//  globalNumNonZeros_ = matrixA_->getGlobalNNZ();
}


/// Destructor
template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
Solver<ConcreteSolver,Matrix,Vector>::~Solver( )
{
  // TODO: The below code is still dependent on there being some
  // kind of status and control code available to Amesos2::Solver,
  // either in the form of a shared library or private inheritance
  // of classes that provide that functionality.
      
  // print out some information if required by the user
  // if ((control_.verbose_ && control_.printTiming_) || control_.verbose_ == 2) printTiming();
  // if ((control_.verbose_ && control_.printStatus_) || control_.verbose_ == 2) printStatus();

}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverBase&
Solver<ConcreteSolver,Matrix,Vector>::preOrdering()
{
  // Teuchos::TimeMonitor LocalTimer1(timers_.overheadTime_);
  // static_cast<solver_type*>(this)->preOrdering_impl();

  return *this;
}
    

template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverBase&
Solver<ConcreteSolver,Matrix,Vector>::symbolicFactorization()
{
  ++status_.numSymbolicFact_;
  Teuchos::TimeMonitor LocalTimer1(timers_.overheadTime_);
  Teuchos::TimeMonitor LocalTimer2(timers_.symFactTime_);
  static_cast<solver_type*>(this)->symbolicFactorization_impl();

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverBase&
Solver<ConcreteSolver,Matrix,Vector>::numericFactorization()
{
  ++status_.numNumericFact_;
  Teuchos::TimeMonitor LocalTimer1(timers_.overheadTime_);
  Teuchos::TimeMonitor LocalTimer2(timers_.numFactTime_);
  static_cast<solver_type*>(this)->numericFactorization_impl();

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
Solver<ConcreteSolver,Matrix,Vector>::solve()
{
  ++status_.numSolve_;
  Teuchos::TimeMonitor LocalTimer1(timers_.overheadTime_);
  Teuchos::TimeMonitor LocalTimer2(timers_.solveTime_);

  static_cast<solver_type*>(this)->solve_impl();
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
bool 
Solver<ConcreteSolver,Matrix,Vector>::matrixShapeOK()
{
  return( static_cast<solver_type*>(this)->matrixShapeOK_impl() );
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
SolverBase&
Solver<ConcreteSolver,Matrix,Vector>::setParameters(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  // The setParameters method should be consitent over all concrete solvers.
  // It will accept general status and control parameters, as well as
  // parameters specific to a solver.  If the solver does not recognize the
  // parameter, then it will simply be ignored
 
  // Do everything here that is for generic status and control parameters
  control_.setControlParameters(parameterList);
  status_.setStatusParameters(parameterList);

  // Finally, hook to the implementation's parameter list parser
  static_cast<solver_type*>(this)->setParameters_impl(parameterList);

  return *this;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
const Teuchos::RCP<const Teuchos::Comm<int> >&
Solver<ConcreteSolver,Matrix,Vector>::getComm() const
{
  return matrixA_->getComm();
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
std::string
Solver<ConcreteSolver,Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << name() << "solver interface";
  return oss.str();
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
Solver<ConcreteSolver,Matrix,Vector>::describe(
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
    Util::printLine();
    out << this->description() << std::endl << std::endl;
    
    out << p << "Matrix has " << globalNumRows_ << "rows"
        << " and " << globalNumNonZeros_ << "nonzeros"
        << std::endl;
    if( vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME ){
      out << p << "Nonzero elements per row = "
          << double(globalNumNonZeros_)/globalNumRows_
          << std::endl;
      out << p << "Percentage of nonzero elements = "
          << 100.0 * globalNumNonZeros_/(pow(double(globalNumRows_),double(2.0)))
          << std::endl;
    }
    if( vl == VERB_HIGH || vl == VERB_EXTREME ){
      out << p << "Use transpose = " << control_.useTranspose_
          << std::endl;
    }
    if ( vl == VERB_EXTREME ){
      printTiming(out,vl);
    }
    Util::printLine();
  }
  return;
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
Solver<ConcreteSolver,Matrix,Vector>::printTiming(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  if( matrixA_.is_null() || (status_.myPID_ != 0) ){ return; }

  double numTime   = timers_.numFactTime_.totalElapsedTime();
  double solTime   = timers_.solveTime_.totalElapsedTime();
  double overhead  = timers_.overheadTime_.totalElapsedTime();
  double totalTime = numTime + solTime;
	
  std::string p = name() + " : ";
  Util::printLine();

  out << p << "Time to convert matrix to implementation format = "
      << timers_.mtxConvTime_.totalElapsedTime() << " (s)"
      << std::endl;
  out << p << "Time to redistribute matrix = "
      << timers_.mtxRedistTime_.totalElapsedTime() << " (s)"
      << std::endl;
  out << p << "Time to redistribute vectors = "
      << timers_.vecRedistTime_.totalElapsedTime() << " (s)"
      << std::endl;
  out << p << "Number of numeric factorizations = "
      << status_.numNumericFact_
      << std::endl;
  out << p << "Time for num fact = "
      << numTime << " (s), avg = "
      << numTime/status_.numNumericFact_ << " (s)"
      << std::endl;
  out << p << "Number of solve phases = "
      << status_.numSolve_
      << std::endl;
  out << p << "Time for solve = "
      << solTime << " (s), avg = "
      << solTime/status_.numSolve_ << " (s)"
      << std::endl;
  out << p << "Total time spent in Amesos2 = "
      << totalTime << " (s)"
      << std::endl;
  out << p << "Total time spent in the Amesos2 interface = "
      << overhead << " (s)"
      << std::endl;
  out << p << "(the above time does not include solver time)"
      << std::endl;
  out << p << "Amesos2 interface time / total time = "
      << overhead/totalTime
      << std::endl;
  Util::printLine();
}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
void
Solver<ConcreteSolver,Matrix,Vector>::getTiming(
  Teuchos::ParameterList& timingParameterList) const
{}


template <template <class,class> class ConcreteSolver, class Matrix, class Vector >
std::string
Solver<ConcreteSolver,Matrix,Vector>::name() const {
  std::string solverName = solver_type::name;
  return solverName;
}


} // end namespace Amesos

#endif	// AMESOS2_SOLVER_DEF_HPP
