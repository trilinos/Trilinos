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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_DERIVUTILS_H
#define LOCA_DERIVUTILS_H

#include "NOX_Abstract_Group.H" // For CopyType

namespace NOX {
  namespace Abstract {
    class MultiVector;
  }
}

namespace LOCA {
  namespace Continuation {
    class AbstractGroup;
  }
  namespace Bifurcation {
    namespace HopfBord {
      class AbstractGroup;
    }
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
  computeDfDp(LOCA::Continuation::AbstractGroup& grp, int param_id,
	      NOX::Abstract::Vector& result) const;

  //! Compute derivative of f with respect to parameter, identified by param_id
  virtual NOX::Abstract::Group::ReturnType 
  computeDfDp(LOCA::Continuation::AbstractGroup& grp, 
	      const vector<int>& param_ids,
	      NOX::Abstract::MultiVector& result,
	      bool isValidF) const;

  //! Compute derivative of Jn with respect to particular parameter param_id
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDp(LOCA::Continuation::AbstractGroup& grp,
	       const NOX::Abstract::Vector& nullVector,
	       int param_id, NOX::Abstract::Vector& result) const;

  //! Compute derivative of Jn with respect to parameter, given JnVector
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDp(LOCA::Continuation::AbstractGroup& grp,
	       const NOX::Abstract::Vector& nullVector,
	       int param_id,
	       const NOX::Abstract::Vector& JnVector,
	       NOX::Abstract::Vector& result) const;

  //! Compute vector (Jn)_{x}a  given vector a
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
		const NOX::Abstract::Vector& nullVector,
		const NOX::Abstract::Vector& aVector,
		NOX::Abstract::Vector& result) const;
  
  //! Compute vector (Jn)_{x}a  given vector a, given JnVector
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
		const NOX::Abstract::Vector& nullVector,
		const NOX::Abstract::Vector& aVector,
		const NOX::Abstract::Vector& JnVector,
		NOX::Abstract::Vector& result) const;

  //! Compute vector (Jn)_{x}a  given multi-vector a
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
		const NOX::Abstract::Vector& nullVector,
		const NOX::Abstract::MultiVector& aVector,
		NOX::Abstract::MultiVector& result) const;
  
  //! Compute vector (Jn)_{x}a  given multi-vector a, given JnVector
  virtual NOX::Abstract::Group::ReturnType 
  computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
		const NOX::Abstract::Vector& nullVector,
		const NOX::Abstract::MultiVector& aVector,
		const NOX::Abstract::Vector& JnVector,
		NOX::Abstract::MultiVector& result) const;

  //! Compute derivative of (J+iwB)(y+iz) with respect to particular 
  //! parameter param_id
  virtual NOX::Abstract::Group::ReturnType 
  computeDCeDp(LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
	       const NOX::Abstract::Vector& yVector,
	       const NOX::Abstract::Vector& zVector,
	       double w,
	       const int param_id, 
	       NOX::Abstract::Vector& result_real,
	       NOX::Abstract::Vector& result_imag) const;

  //! Compute derivative of (J+iwB)(y+iz) with respect to parameter, 
  //! given (J+iwB)(y+iz) vectors
  virtual NOX::Abstract::Group::ReturnType 
  computeDCeDp(LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
	       const NOX::Abstract::Vector& yVector,
	       const NOX::Abstract::Vector& zVector,
	       double w,
	       const int param_id, 
	       const NOX::Abstract::Vector& Ce_real,
	       const NOX::Abstract::Vector& Ce_imag,
	       NOX::Abstract::Vector& result_real,
	       NOX::Abstract::Vector& result_imag) const;

  //! Compute vector (J+iwB)(y+iz))_{x}a  given vector a
  virtual NOX::Abstract::Group::ReturnType 
  computeDCeDxa(LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
		const NOX::Abstract::Vector& yVector,
		const NOX::Abstract::Vector& zVector,
		double w,
		const NOX::Abstract::Vector& aVector,
		NOX::Abstract::Vector& result_real,
		NOX::Abstract::Vector& result_imag) const;
  
  //! Compute vector (J+iwB)(y+iz))_{x}a  given vector a, given 
  //! (J+iwB)(y+iz) vector
  virtual NOX::Abstract::Group::ReturnType 
  computeDCeDxa(LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
		const NOX::Abstract::Vector& yVector,
		const NOX::Abstract::Vector& zVector,
		double w,
		const NOX::Abstract::Vector& aVector,
		const NOX::Abstract::Vector& Ce_real,
		const NOX::Abstract::Vector& Ce_imag,
		NOX::Abstract::Vector& result_real,
		NOX::Abstract::Vector& result_imag) const;

  //! Compute vector (J+iwB)(y+iz))_{x}a  given multi-vector a
  virtual NOX::Abstract::Group::ReturnType 
  computeDCeDxa(LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
		const NOX::Abstract::Vector& yVector,
		const NOX::Abstract::Vector& zVector,
		double w,
		const NOX::Abstract::MultiVector& aVector,
		NOX::Abstract::MultiVector& result_real,
		NOX::Abstract::MultiVector& result_imag) const;
  
  //! Compute vector (J+iwB)(y+iz))_{x}a  given multi-vector a, given 
  //! (J+iwB)(y+iz) vector
  virtual NOX::Abstract::Group::ReturnType 
  computeDCeDxa(LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
		const NOX::Abstract::Vector& yVector,
		const NOX::Abstract::Vector& zVector,
		double w,
		const NOX::Abstract::MultiVector& aVector,
		const NOX::Abstract::Vector& Ce_real,
		const NOX::Abstract::Vector& Ce_imag,
		NOX::Abstract::MultiVector& result_real,
		NOX::Abstract::MultiVector& result_imag) const;

protected:

  //! Perturb parameter param_id in the group's parameter vector,
  //! and return perturbation size
  virtual double perturbParam(LOCA::Continuation::AbstractGroup& grp,
			      double& param,
                              int param_id) const;

  //! Perturb group's solution vector in direction of aVector,
  //! and return perturbation size
  virtual double perturbXVec(LOCA::Continuation::AbstractGroup& grp,
			     const NOX::Abstract::Vector& xVector,
                             const NOX::Abstract::Vector& aVector) const;

  //! Return perturbation for finite differencing around given double variable
  virtual double epsScalar(double p) const;

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
