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

#ifndef LOCA_DERIVUTILS_H
#define LOCA_DERIVUTILS_H

#include "NOX_Abstract_Group.H" // For CopyType

namespace LOCA {
  namespace Abstract {
    class Group;
  }
}

namespace LOCA {

/*! 
  \brief %LOCA's generic derivative computation class to compute various
  derivatives via finite differencing.

  The %DerivUtils class provides generic methods to compute the following
  derivatives:
  \f[
      \frac{\partial F}{\partial p},\quad\frac{\partial Jn}{\partial p},\quad\frac{\partial Jn}{\partial x}a
  \f]
  where \f$J = \partial F/\partial x\f$ and \f$n\f$, \f$a\f$ are vectors.  These quantities are calculate by finite differencing.
  
*/

class DerivUtils {

public:
  
  /*! \brief Default constructor.  perturb is the relative perturbation size 
    used in differencing calculations. */
  DerivUtils(double perturb = 1.0e-6);

  //! Copy constructor.
  DerivUtils(const DerivUtils&);

  //! Destructor.
  virtual ~DerivUtils();

  /*! \brief Cloning function.  Creates a copy of the %DerivUtils object 
    of the same type. */
  virtual DerivUtils* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Compute derivative of f with respect to parameter, identified by param_id
  virtual NOX::Abstract::Group::ReturnType 
  computeDfDp(LOCA::Abstract::Group& grp, const int param_id,
	      NOX::Abstract::Vector& result) const;

  //! Compute derivative of Jn with respect to particular parameter param_id
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDp(LOCA::Abstract::Group& grp,
	       const NOX::Abstract::Vector& nullVector,
	       const int param_id, NOX::Abstract::Vector& result) const;

  //! Compute derivative of Jn with respect to parameter, given JnVector
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDp(LOCA::Abstract::Group& grp,
	       const NOX::Abstract::Vector& nullVector,
	       const int param_id,
	       const NOX::Abstract::Vector& JnVector,
	       NOX::Abstract::Vector& result) const;

  //! Compute vector (Jn)_{x}a  given vector a
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDxa(LOCA::Abstract::Group& grp,
		const NOX::Abstract::Vector& nullVector,
		const NOX::Abstract::Vector& aVector,
		NOX::Abstract::Vector& result) const;
  
  //! Compute vector (Jn)_{x}a  given vector a, given JnVector
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDxa(LOCA::Abstract::Group& grp,
		const NOX::Abstract::Vector& nullVector,
		const NOX::Abstract::Vector& aVector,
		const NOX::Abstract::Vector& JnVector,
		NOX::Abstract::Vector& result) const;

  //! Compute vector (Jn)_{x}a+(Jn)_p given vector a, param_id 
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDxaDp(LOCA::Abstract::Group& grp,
		  const NOX::Abstract::Vector& nullVector,
		  const NOX::Abstract::Vector& aVector,
		  const int param_id,
		  NOX::Abstract::Vector& result) const;

protected:

  //! Given 2 vectors and a perturbation, return differenced vector. 
  /*! The first argument (perturbedVector) is modified.  */
  virtual NOX::Abstract::Vector& doDifference(
                                NOX::Abstract::Vector& perturbedVector,
                                const NOX::Abstract::Vector& baseVector,
                                const double eps) const;

  //! Perturb parameter param_id in the group's parameter vector,
  //! and return perturbation size
  virtual double perturbParam(LOCA::Abstract::Group& grpPerturbed,
                              const int param_id) const;

  //! Perturb group's solution vector in direction of aVector,
  //! and return perturbation size
  virtual double perturbXVec(LOCA::Abstract::Group& grpPerturbed,
                             const NOX::Abstract::Vector& aVector) const;

  //! Return perturbation for finite differencing around given double variable
  virtual double epsScalar(const double p) const;

  //! Return perturbation for directional derivatives around first vector
  //! in direction of second vector
  virtual double epsVector(const NOX::Abstract::Vector& xVector,
                           const NOX::Abstract::Vector& aVector) const;

protected:

  //! Relative perturbation size using in differencing calculations.
  double perturb;

}; // class DerivUtils
} // namespace LOCA

#endif
