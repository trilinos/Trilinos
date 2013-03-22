// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef BRUSSELATORPROBLEMINTERFACE_H
#define BRUSSELATORPROBLEMINTERFACE_H

#include "Teuchos_RCP.hpp"

#include "LOCA_LAPACK_Interface.H" // base class
#include "NOX_LAPACK_Vector.H" // class data element
#include <fstream>

namespace LOCA {
  class GlobalData;
}

/*!
 *\brief LOCA-LAPACK problem interface for the Brussselator problem.
 */
/*!
 * %BrusselatorProblemInterface implements the LOCA::LAPACK::Interface for a 
 * 1D finite-difference discretization of the Brusselator problem:
 * \f[
 *     \frac{\partial T}{\partial t} = D_1\frac{\partial^2 T}{\partial x^2} + 
 *       \alpha - (\beta + 1)T + T^2 C \\
 *     \frac{\partial C}{\partial t} = D_2\frac{\partial^2 C}{\partial x^2} + 
 *       \beta T - T^2 C
 * \f]
 * subject to the boundar conditions \f$T(0) = T(1) = \alpha\f$, 
 * \f$C(0) = C(1) = \beta/\alpha\f$.  The parameters are \f$\alpha\f$, 
 * \f$\beta\f$, \f$D_1\f$, \f$D_2\f$ and \f$n\f$, the size of the 
 * discretization.  
 */
class BrusselatorProblemInterface : 
  public LOCA::LAPACK::Interface {

public:
 
  //! Constructor.  Also generates an initial guess
  BrusselatorProblemInterface(
		     const Teuchos::RCP<LOCA::GlobalData>& global_data,
		     int N, double a, double b, double d1, double d2,
		     std::ofstream& file);

  //! Destructor
  ~BrusselatorProblemInterface() {};

  //! Returns initial guess
  const NOX::LAPACK::Vector& getInitialGuess();

  //! Evaluates residual
  bool
  computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector &x);
  
  //! Evalues jacobian
  bool
  computeJacobian(NOX::LAPACK::Matrix<double>& J, 
		  const NOX::LAPACK::Vector & x);

  /*!
   * \brief Compute shifted matrix alpha*J + beta*M where J is the 
   * Jacobian matrix and M is the mass matrix.  Returns true if 
   * computation was successful.
   */
  virtual bool computeShiftedMatrix(double alpha, double beta,
				    const NOX::LAPACK::Vector& x,
				    NOX::LAPACK::Matrix<double>& A);

  //! Sets parameters
  void setParams(const LOCA::ParameterVector& p);

  //! Prints solution after successful step
  void printSolution(const NOX::LAPACK::Vector &x, const double conParam);

private:

   //! Global data
  Teuchos::RCP<LOCA::GlobalData> globalData;

  //! Stores the initial guess
  NOX::LAPACK::Vector initialGuess;

  //! Stores \f$\alpha\f$
  double alpha;

  //! Stores \f$\beta\f$
  double beta;

  //! Stores \f$D_1\f$
  double D1;

  //! Stores \f$D_2\f$
  double D2;

  //! Stores the size of the discretization
  int n;

  //! Stores reference to output file stream
  std::ofstream& outputFile;

};

#endif
