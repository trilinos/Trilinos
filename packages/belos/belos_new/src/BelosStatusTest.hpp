
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

#ifndef BELOS_STATUS_TEST_HPP
#define BELOS_STATUS_TEST_HPP

/*!
  \file BelosStatusTest.hpp
  \brief Pure virtual base class for defining the status testing capabilities of Belos.
*/

#include "BelosTypes.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosConfigDefs.hpp"

  /*! 
    \class Belos::StatusTest
    \brief A pure virtual class for defining the status tests for the Belos iterative solvers
    
    Belos::StatusTest is an interface that can be implemented to create convergence tests for
    all Belos solvers.  Almost any kind of test can be expressed using this mechanism, 
    including composite tests (see Belos::StatusTestCombo). 
  */

namespace Belos {

template <class TYPE, class OP, class MV>
class StatusTest {

 public:
  //@{ \name Constructors/destructors.

  //! Constructor
  StatusTest() {};

  //! Destructor
  virtual ~StatusTest() {};
  //@}

  //@{ \name Status methods
  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met.  The calling routine may pass in the
    current native residual vector (the one naturally produced as part of the iterative method) or a
    pre-computed estimate of the two-norm of the current residual, or both or neither.  The calling routine
    should also indicate if the solution of the linear problem has been updated to be compatible with
    the residual.  Some methods, such as GMRES do not update the solution at each iteration.

    \return Belos::StatusType: Unconverged, Converged or Failed.
  */
  virtual StatusType CheckStatus( IterativeSolver<TYPE,OP,MV>* iSolver ) = 0;

  //! Return the result of the most recent CheckStatus call.
  virtual StatusType GetStatus() const = 0;
  //@}

  //@{ \name Reset methods
  //! Informs the convergence test that it should reset its internal configuration to the initialized state.
  /*! This is necessary for the case when the status test is being reused by another solver or for another
    linear problem.  The status test may have information that pertains to a particular linear system.  The
    internal information will be reset back to the initialized state.  The user specified information that
    the convergence test uses will remain.
  */
  virtual void Reset() = 0;
  //@}

  //@{ \name Attribute methods

  //! Indicates if residual vector is required by this convergence test.
  /*! If this method returns true, then the ResidualVector argument to the Converged() method will
    defined.  If this method returns false, then the ResidualVector may not be defined when Converged() is
    called.  Some iterative methods do not explicitly construct the residual vector at each iteration.  Thus,
    if this vector is not required, this vector will not need to be constructed if this method returns false.
  */
  virtual bool ResidualVectorRequired() const = 0;
  //@}

  //@{ \name Print methods

  //! Output formatted description of stopping test to output stream.
  virtual ostream& Print(ostream& os, int indent = 0) const = 0;
 
  //! Output the result of the most recent CheckStatus call.
  virtual void PrintStatus(ostream& os, StatusType type) const {
    os << setiosflags(ios::left) << setw(13) << setfill('.');
    switch (type) {
    case  Failed:
      os << "Failed";
      break;
    case  Converged:
      os << "Converged";
      break;
    case  Unconverged:
    default:
      os << "**";
      break;
    }
    os << setiosflags(ios::left) << setfill(' ');
    return;
  };
  //@}

};

} // end of Belos namespace

#endif /* BELOS_STATUS_TEST_HPP */
