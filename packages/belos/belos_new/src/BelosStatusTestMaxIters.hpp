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

template<class TYPE>
class StatusTestMaxIters: public StatusTest<TYPE> {

 public:

  //@{ \name Constructor/Destructor.

  //! Constructor
  StatusTestMaxIters(int maxIters);

  //! Destructor
  virtual ~StatusTestMaxIters() {};
  //@}

  //@{ \name Status methods

  //! Check convergence status of the iterative solver: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met using the current information from the 
    iterative solver.
  */
  StatusType CheckStatus(IterativeSolver<TYPE> *iSolver );

  //! Return the result of the most recent CheckStatus call.
  StatusType GetStatus() const {return(status_);};

  //@}

  //@{ \name Accessor methods

  //! Returns the maximum number of iterations set in the constructor.
  int GetMaxIters() const { return(maxIters_); };

  //! Returns the current number of iterations from the most recent StatusTest call.
  int GetNumIters() const { return(nIters_); };

  //@}

  //@{ \name Attribute methods

  //! Indicates if residual vector is required by this convergence test: returns false for this class.
  bool ResidualVectorRequired() const { return(false); } ;
  //@}

  //@{ \name Print methods

  //! Output formatted description of stopping test to output stream
  ostream& Print(ostream& os, int indent = 0) const;

  //@}
  

private:

  //@{ \name Private data members.
  //! Maximum number of iterations allowed
  int maxIters_;

  //! Current number of iterations
  int nIters_;

  //! Status
  StatusType status_;
  //@}

};

  template<class TYPE> 
  StatusTestMaxIters<TYPE>::StatusTestMaxIters(int maxIters)
  {
    if (maxIters < 1)
      maxIters_ = 1;
    else
      maxIters_ = maxIters;
    
    nIters_ = 0;
    status_ = Unchecked;
  }
  
  template<class TYPE>
  StatusType StatusTestMaxIters<TYPE>::CheckStatus(IterativeSolver<TYPE> *iSolver )
  {
    status_ = Unconverged;
    nIters_ = iSolver->GetNumIters();
    if (nIters_ >= maxIters_)
      status_ = Failed;
    return status_;
  }
  
  template<class TYPE>
  ostream& StatusTestMaxIters<TYPE>::Print(ostream& os, int indent) const
  {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    PrintStatus(os, status_);
    os << "Number of Iterations = ";
    os << nIters_;
    os << ((nIters_ < maxIters_) ? " < " : ((nIters_ == maxIters_) ? " = " : " > "));
    os << maxIters_;
    os << endl;
    return os;
  }
  
} // end Belos namespace

#endif /* BELOS_STATUS_TEST_MAXITERS_HPP */
