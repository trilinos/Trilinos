// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
   
#ifndef BELOS_STATUS_TEST_LOGRESNORM_HPP
#define BELOS_STATUS_TEST_LOGRESNORM_HPP

/*!
  \file BelosStatusTestLogResNorm.hpp
  \brief Belos::StatusTest debugging class for storing the absolute residual norms generated during a solve.
*/

#include "BelosStatusTest.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RCP.hpp"

/*! \class Belos::StatusTestLogResNorm: 
    \brief A Belos::StatusTest debugging class for storing the absolute residual norms generated during a solve.

    This implementation of the Belos::StatusTest base class stores the absolute residual norm provided by the iteration
    throughout an entire solve.
  
    \note This debugging status test will only work for single-rhs, single block solves.
*/

namespace Belos {

template <class ScalarType, class MV, class OP>
class StatusTestLogResNorm: public StatusTest<ScalarType,MV,OP> {

public:
  //! The type of the magnitude (absolute value) of a ScalarType.
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

private:
  //! @name Abbreviations for method implementations
  //@{
  typedef MultiVecTraits<ScalarType,MV> MVT;
  //@}

 public:

   //! @name Constructor/Destructor.
  //@{ 

  //! Constructor
  StatusTestLogResNorm(int maxIters);

  //! Destructor
  virtual ~StatusTestLogResNorm() {};
  //@}

  //! @name Status methods
  //@{ 

  //! Check convergence status of the iterative solver: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met using the current information from the 
    iterative solver.
  */
  StatusType checkStatus(Iteration<ScalarType,MV,OP> *iSolver );

  //! Return the result of the most recent CheckStatus call.
  StatusType getStatus() const {return(Undefined);}

  //@}

  //! @name Reset methods
  //@{ 

  //! Resets the status test to the initial internal state.
  void reset();

  //! Sets the maximum number of iterations allowed so internal storage can be resized.
  void setMaxIters(int maxIters) { maxIters_ = maxIters; logResNorm_.reserve( maxIters_ ); }

  //@}

  //! @name Accessor methods
  //@{ 

  //! Returns the maximum number of iterations set in the constructor.
  int getMaxIters() const { return(maxIters_); }

  //! Returns the current number of iterations from the most recent StatusTest call.
  int getNumIters() const { return(nIters_); }

  //! Returns the log of the absolute residual norm from the iteration.
  const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>&
      getLogResNorm() const { return(logResNorm_); }

  //@}

  //! @name Print methods
  //@{ 

  //! Output formatted description of stopping test to output stream.
  void print(std::ostream& os, int indent = 0) const;

  //! Print message for each status specific to this stopping test.
  void printStatus(std::ostream& os, StatusType type) const;

  //@}
 
  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief Method to return description of the debugging status test  */
  std::string description() const 
  {  
    std::ostringstream oss; 
    oss << "Belos::StatusTestLogResNorm<>: [ " << getNumIters() << " / " << getMaxIters() << " ]"; 
    return oss.str();
  }
  //@} 

private:

  //! @name Private data members.
  //@{ 
  //! Maximum number of iterations allowed
  int maxIters_;

  //! Current number of iterations
  int nIters_;

  //! Log of absolute residual norm
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> logResNorm_;

  //@}

};

  template <class ScalarType, class MV, class OP> 
  StatusTestLogResNorm<ScalarType,MV,OP>::StatusTestLogResNorm(int maxIters)
  {
    if (maxIters < 1)
      maxIters_ = 1;
    else
      maxIters_ = maxIters;

    logResNorm_.reserve( maxIters_ );
    
    nIters_ = 0;
  }
  
  template <class ScalarType, class MV, class OP>
  StatusType StatusTestLogResNorm<ScalarType,MV,OP>::checkStatus(Iteration<ScalarType,MV,OP> *iSolver )
  {
    // Check that this solve is a single-vector, single-block.
    const LinearProblem<ScalarType,MV,OP>& lp = iSolver->getProblem ();
    int blkSize = lp.getLSIndex().size();
    int numRHS = MVT::GetNumberVecs( *(lp.getRHS()) );

    int currIters = iSolver->getNumIters();

    if ( (numRHS==1) && (blkSize==1) && (currIters!=nIters_) )
    {
      std::vector<MagnitudeType> tmp_resvector( 1 );
      Teuchos::RCP<const MV> residMV = iSolver->getNativeResiduals (&tmp_resvector);
      if (! residMV.is_null ()) 
      {
        // We got a multivector back.  Compute the norms explicitly.
        MVT::MvNorm (*residMV, tmp_resvector, TwoNorm);
      }

      logResNorm_.push_back( tmp_resvector[0] );
      nIters_ = currIters;
    }
 
    return Undefined;
  }
  
  template <class ScalarType, class MV, class OP>
  void StatusTestLogResNorm<ScalarType,MV,OP>::reset()
  {
    nIters_ = 0;
    logResNorm_.clear();
    logResNorm_.reserve( maxIters_ );
  }    
    
  template <class ScalarType, class MV, class OP>
  void StatusTestLogResNorm<ScalarType,MV,OP>::print(std::ostream& os, int indent) const
  {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    printStatus(os, Undefined);
    os << "Logging Absolute Residual 2-Norm" << std::endl;
  }
 
  template <class ScalarType, class MV, class OP>
  void StatusTestLogResNorm<ScalarType,MV,OP>::printStatus(std::ostream& os, StatusType type) const 
  {
    os << std::left << std::setw(13) << std::setfill('.');
    os << "**";
    os << std::left << std::setfill(' ');
    return;
  } 

} // end Belos namespace

#endif /* BELOS_STATUS_TEST_LOGRESNORM_HPP */
