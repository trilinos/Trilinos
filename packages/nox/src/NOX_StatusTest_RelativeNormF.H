// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_STATUSTEST_RELATIVENORMF_H
#define NOX_STATUSTEST_RELATIVENORMF_H

#include "NOX_StatusTest_Generic.H"	// base class
#include "NOX_Abstract_Vector.H"        // for NormType
#include "NOX_Utils.H"                  // for cerr

// Forward declaration
namespace NOX {
namespace Abstract {
class Group;
}
}

namespace NOX {

namespace StatusTest {

/** \brief Requires reduction in residual norm compared to initial norm.

    At iteration zero, this computes \f$||F_0||\f$ and is unconverged.

    Returns a converged status if \f$||F_k|| < ||F_0|| *
    \mbox{tolerance}\f$, where \f$k\f$ is the nonlinear iteration number.

*/
class RelativeNormF : public Generic {

public:

  RelativeNormF(double tolerance,
		bool scale_norms_by_vector_length = false,
		const NOX::Utils* u = NULL);

  // derived
  virtual NOX::StatusTest::StatusType 
  checkStatus(const NOX::Solver::Generic& problem,
	      NOX::StatusTest::CheckType checkType);

  // derived
  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  NOX::StatusTest::StatusType status;

  //! Tolerance required for convergence
  double tolerance;

  //! Initial Norm of F
  double normF_0;

  //! Current Norm of F
  double normF;

  //! If set to true, the norms are scaled by the size of the vector
  bool scale_norms_by_vector_length;

  //! Ostream used to print errors
  NOX::Utils utils;

};

} // namespace Status
} // namespace NOX

#endif
