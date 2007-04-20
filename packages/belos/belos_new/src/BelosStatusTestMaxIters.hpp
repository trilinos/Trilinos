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
  StatusType getStatus() const {return(status_);};

  //@}

  //! @name Reset methods
  //@{ 

  //! Resets the status test to the initial internal state.
  void reset();

  //@}

  //! @name Accessor methods
  //@{ 

  //! Returns the maximum number of iterations set in the constructor.
  int getMaxIters() const { return(maxIters_); };

  //! Returns the current number of iterations from the most recent StatusTest call.
  int getNumIters() const { return(nIters_); };

  //@}

  //! @name Attribute methods
  //@{ 

  //! Indicates if residual vector is required by this convergence test: returns false for this class.
  bool residualVectorRequired() const { return(false); } ;
  //@}

  //! @name Print methods
  //@{ 

  //! Output formatted description of stopping test to output stream.
  void print(ostream& os, int indent = 0) const;

  //! Print message for each status specific to this stopping test.
  void printStatus(ostream& os, StatusType type) const;

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
  void StatusTestMaxIters<ScalarType,MV,OP>::print(ostream& os, int indent) const
  {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    printStatus(os, status_);
    os << "Number of Iterations = ";
    os << nIters_;
    os << ((nIters_ < maxIters_) ? " < " : ((nIters_ == maxIters_) ? " == " : " > "));
    os << maxIters_;
    os << endl;
  }
 
  template <class ScalarType, class MV, class OP>
  void StatusTestMaxIters<ScalarType,MV,OP>::printStatus(ostream& os, StatusType type) const 
  {
    os << setiosflags(ios::left) << setw(13) << setfill('.');
    switch (type) {
    case  Passed:
      os << "Failed";
    case  Failed:
      os << "OK";
      break;
    case  Undefined:
    default:
      os << "**";
      break;
    }
    os << setiosflags(ios::left) << setfill(' ');
    return;
  } 

} // end Belos namespace

#endif /* BELOS_STATUS_TEST_MAXITERS_HPP */
