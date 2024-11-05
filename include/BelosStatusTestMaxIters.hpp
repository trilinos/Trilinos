// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
   
#ifndef BELOS_STATUS_TEST_MAXITERS_HPP
#define BELOS_STATUS_TEST_MAXITERS_HPP

/*!
  \file BelosStatusTestMaxIters.hpp
  \brief Belos::StatusTest class for specifying a maximum number of iterations.
*/

#include "BelosStatusTest.hpp"

/*! \class Belos::StatusTestMaxIters: 
    \brief A Belos::StatusTest class for specifying a maximum number of iterations.

    This implementation of the Belos::StatusTest base class tests the number of iterations performed
    against a maximum number allowed.
*/

namespace Belos {

template <class ScalarType, class MV, class OP>
class StatusTestMaxIters: public StatusTest<ScalarType,MV,OP> {

 public:

   //! @name Constructor/Destructor.
  //@{ 

  //! Constructor
  StatusTestMaxIters(int maxIters);

  //! Destructor
  virtual ~StatusTestMaxIters() {};
  //@}

  //! @name Status methods
  //@{ 

  //! Check convergence status of the iterative solver: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met using the current information from the 
    iterative solver.
  */
  StatusType checkStatus(Iteration<ScalarType,MV,OP> *iSolver );

  //! Return the result of the most recent CheckStatus call.
  StatusType getStatus() const {return(status_);}

  //@}

  //! @name Reset methods
  //@{ 

  //! Resets the status test to the initial internal state.
  void reset();

  //! Sets the maximum number of iterations allowed.
  void setMaxIters(int maxIters) { maxIters_ = maxIters; }

  //@}

  //! @name Accessor methods
  //@{ 

  //! Returns the maximum number of iterations set in the constructor.
  int getMaxIters() const { return(maxIters_); }

  //! Returns the current number of iterations from the most recent StatusTest call.
  int getNumIters() const { return(nIters_); }

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

  /** \brief Method to return description of the maximum iteration status test  */
  std::string description() const 
  {  
    std::ostringstream oss; 
    oss << "Belos::StatusTestMaxIters<>: [ " << getNumIters() << " / " << getMaxIters() << " ]"; 
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

  //! Status
  StatusType status_;
  //@}

};

  template <class ScalarType, class MV, class OP> 
  StatusTestMaxIters<ScalarType,MV,OP>::StatusTestMaxIters(int maxIters)
  {
    if (maxIters < 1)
      maxIters_ = 1;
    else
      maxIters_ = maxIters;
    
    nIters_ = 0;
    status_ = Undefined;
  }
  
  template <class ScalarType, class MV, class OP>
  StatusType StatusTestMaxIters<ScalarType,MV,OP>::checkStatus(Iteration<ScalarType,MV,OP> *iSolver )
  {
    status_ = Failed;
    nIters_ = iSolver->getNumIters();
    if (nIters_ >= maxIters_)
      status_ = Passed;
    return status_;
  }
  
  template <class ScalarType, class MV, class OP>
  void StatusTestMaxIters<ScalarType,MV,OP>::reset()
  {
    nIters_ = 0;
    status_ = Undefined;
  }    
    
  template <class ScalarType, class MV, class OP>
  void StatusTestMaxIters<ScalarType,MV,OP>::print(std::ostream& os, int indent) const
  {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    printStatus(os, status_);
    os << "Number of Iterations = ";
    os << nIters_;
    os << ((nIters_ < maxIters_) ? " < " : ((nIters_ == maxIters_) ? " == " : " > "));
    os << maxIters_;
    os << std::endl;
  }
 
  template <class ScalarType, class MV, class OP>
  void StatusTestMaxIters<ScalarType,MV,OP>::printStatus(std::ostream& os, StatusType type) const 
  {
    os << std::left << std::setw(13) << std::setfill('.');
    switch (type) {
    case  Passed:
      os << "Failed";
      break;
    case  Failed:
      os << "OK";
      break;
    case  Undefined:
    default:
      os << "**";
      break;
    }
    os << std::left << std::setfill(' ');
    return;
  } 

} // end Belos namespace

#endif /* BELOS_STATUS_TEST_MAXITERS_HPP */
