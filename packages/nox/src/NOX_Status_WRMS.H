// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_STATUS_WRMS_H
#define NOX_STATUS_WRMS_H

#include "NOX_Status_Test.H"	// base class

// Forward declarations
namespace NOX {
namespace Abstract {
class Vector;
}
}

namespace NOX {

namespace Status {

//! Convergence test based on the weighted root mean square norm used in the PVODE code.
/*!
  Convergence of the nonlinear solver is achieved if the following criteria is met:

  \f[
  ||E_n||_{rms,w} < 1.0
  \f]
 
  where \f$ E_n \f$ is the vector containing the local error estimate, \f$ C_q*(x_n-x_{n-1)} \f$, at nonlinear iteration \f$ n \f$:  

  \f[
  ||E_n||_{rms,w} =  \left[ \frac{1}{N} \sum_{i=1}^N \frac {E_{n,i}}{RTOL |x_i| + ATOL_i} \right]^\frac{1}{2} 
  \f]

  Here, \f$ E_{n,i} \f$ denotes the \f$ i \f$ -th component of \f$ E_n \f$, \f$RTOL \f$ is the relative error tolerance, and \f$ ATOL_i \f$ is the absolute error tolerance.  \f$ ATOL_i \f$ can be a scalar or a vector of length \f$ N \f$. \f$ C_q \f$ is a constant that depends on the time stepping algorithm (which BDF method is chosen).  

Reference: "PVODE, an ODE %Solver for Parallel Computers", Byrne, G. D., and Hindmarsh, A. C., Center for Applied Scientific Computing (CASC), UCRL-JC-132361, Rev. 1, May 1999.

*/
class WRMS : public Test {

public:

  //! Constructor where ATOL is a scalar
  WRMS(double rtol, double atol, double BDFMultiplier = 1.0);

  //! Constructor where ATOL is a vector
  WRMS(double rtol, Abstract::Vector& atol, double BDFMultiplier = 1.0);

  //! Destructor.
  virtual ~WRMS();

  virtual StatusType operator()(const Solver::Generic& problem);

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  //! Relative tolerance for convergence.
  double rtol;

  //! Absolute tolerance for convergence (scalar).
  double atol;

  //! flag to tell whether atol is a scalar or a vector.
  bool atolIsScalar;

  //! Absolute tolerance for convergence (vector).
  Abstract::Vector* atoli;

  //! Time integration method multiplier (BDF Multiplier).
  double factor;

  //! %Status
  StatusType status;

  //! Temporary vector used in computation
  Abstract::Vector* u;

  //! Temporary vector used in computation
  Abstract::Vector* v;

};

} // namespace Status
} // namespace NOX

#endif
