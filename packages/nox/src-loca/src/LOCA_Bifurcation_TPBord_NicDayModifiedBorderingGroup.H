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

#ifndef LOCA_BIFURCATION_TPBORD_NICDAYMODIFIEDBORDERINGGROUP_H
#define LOCA_BIFURCATION_TPBORD_NICDAYMODIFIEDBORDERINGGROUP_H

#include "LOCA_Bifurcation_TPBord_ExtendedGroup.H"         // base class

namespace LOCA {
  class ParameterVector;
}
namespace LOCA { 

  namespace Bifurcation { 
    
    namespace TPBord {

      /*! 
       * \brief Turning point group implementing a modified bordering 
       * algorithm.
       */
      /*!
       * The %NicDayModifiedBorderingGroup is derived from the 
       * LOCA::Bifurcation::TPBord::ExtendedGroup and reimplements the 
       * applyJacobianInverse() and applyJacobianInverseMulti() methods
       * using a modified bordering algorithm that avoids nearly-singular
       * matrix solves.  The algorithm is similar to the ModifedBorderingGroup
       * algorithm but uses 5 solves of the original Jacobian instead of
       * bordered solves.  See applyJacobianInverse() for a 
       * description of the algorithm.
       */
      class NicDayModifiedBorderingGroup 
	: public LOCA::Bifurcation::TPBord::ExtendedGroup {

      public:

	/*! 
	 * \brief Constructor with initial data passed through parameter lists.
	 */
	NicDayModifiedBorderingGroup(
			       LOCA::Bifurcation::TPBord::AbstractGroup& g, 
			       NOX::Parameter::List& bifParamList);
  
	/*! 
	 * \brief Constructor.  See LOCA::Bifurcation::TPBord::ExtendedGroup
	 * for a description.
	 */
	NicDayModifiedBorderingGroup(
		       LOCA::Bifurcation::TPBord::AbstractGroup& g, 
		       const NOX::Abstract::Vector& lenVec,
		       const NOX::Abstract::Vector& nullVec, int paramId);

	/*! 
	 * \brief Constructor.  See LOCA::Bifurcation::TPBord::ExtendedGroup
	 * for a description.
	 */
	NicDayModifiedBorderingGroup(
			 const LOCA::Bifurcation::TPBord::AbstractGroup& g, 
			 const NOX::Abstract::Vector& lenVec, 
			 const NOX::Abstract::Vector& nullVec, int paramId);

	//! Copy constructor.
	NicDayModifiedBorderingGroup(
			       const NicDayModifiedBorderingGroup& source, 
			       NOX::CopyType type = NOX::DeepCopy);

	//! Destructor.
	virtual ~NicDayModifiedBorderingGroup();

	//! Assignment operator
	virtual NicDayModifiedBorderingGroup& operator=(
				  const NicDayModifiedBorderingGroup& source);

	//! Assignment operator
	virtual NOX::Abstract::Group& 
	operator=(const NOX::Abstract::Group& source);

	//! Assignment operator
	virtual LOCA::Continuation::AbstractGroup& 
	operator=(const LOCA::Continuation::AbstractGroup& source);

	//! Assignment operator
	virtual LOCA::Extended::AbstractGroup& 
	operator=(const LOCA::Extended::AbstractGroup& source);

	//! Assignment operator
	virtual LOCA::Bifurcation::TPBord::ExtendedGroup& 
	operator=(const LOCA::Bifurcation::TPBord::ExtendedGroup& source);
     
	//! Cloning function
	virtual NOX::Abstract::Group* 
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	/*! 
	 * \brief Applies the inverse of the extended Jacobian matrix using 
	 * a modified bordering algorithm.
	 */
	/*!
	   Computes the solution to the turning point system
	   \f[
	      \begin{bmatrix}
	         J & 0 & \frac{\partial F}{\partial p} \\
		 \frac{\partial Jn}{\partial x} & J & \frac{\partial Jn}{\partial p} \\
		 0 & l^T & 0
	      \end{bmatrix}
	      \begin{bmatrix}
	         u \\
		 v \\
		 w
	      \end{bmatrix}
	      = 
	      \begin{bmatrix}
	         x \\
		 y \\
		 z
	      \end{bmatrix}
	  \f]
	  using the following modified bordering algorithm:
	  \f[
	       \begin{split}
	          \begin{bmatrix}
		     J & s \\
		     n^T & 0
		   \end{bmatrix}
		   \begin{bmatrix}
		      a \\
		      \bar{a}
		   \end{bmatrix} &= 
		   \begin{bmatrix} 
		      x \\
		      0
		   \end{bmatrix} \\
		   \begin{bmatrix}
		     J & s \\
		     n^T & 0
		   \end{bmatrix}
		   \begin{bmatrix}
		      b \\
		      \bar{b}
		   \end{bmatrix}&= 
		   \begin{bmatrix} 
		      \frac{\partial F}{\partial p} \\
		      0
		   \end{bmatrix} \\
		   \begin{bmatrix}
		     J & s \\
		     n^T & 0
		   \end{bmatrix}
		   \begin{bmatrix}
		      c \\
		      \bar{c}
		   \end{bmatrix}&= 
		   \begin{bmatrix} 
		      y - \frac{\partial Jn}{\partial x}a \\
		      0
		   \end{bmatrix} \\
		   \begin{bmatrix}
		     J & s \\
		     n^T & 0
		   \end{bmatrix}
		   \begin{bmatrix}
		      d \\
		      \bar{d}
		   \end{bmatrix}&= 
		   \begin{bmatrix} 
		      \frac{\partial Jn}{\partial p} - \frac{\partial Jn}{\partial x}a \\
		      0
		   \end{bmatrix} \\
		   \begin{bmatrix}
		     J & s \\
		     n^T & 0
		   \end{bmatrix}
		   \begin{bmatrix}
		      d \\
		      \bar{d}
		   \end{bmatrix}&= 
		   \begin{bmatrix} 
		      \frac{\partial Jn}{\partial x}n \\
		      0
		   \end{bmatrix} \\
		   \begin{bmatrix}
		      \sigma & \bar{e} & -l^Te \\
		      0 & \sigma & l^Tv \\
		      \bar{b} & \bar{d} & -l^Td
		   \end{bmatrix}
		   \begin{bmatrix}
                      \alpha \\
		      \beta \\
		      w
		   \end{bmatrix} &= 
		   \begin{bmatrix}
		      \bar{a} \\
		      \bar{c} \\
		      z - l^T c
		   \end{bmatrix} \\
		  u &= a - w b + \alpha v\\
		  v &= c - w d - \alpha e + \beta v
	       \end{split}
	  \f]
	  where \f$\sigma = \|Jn\|\f$ and \f$s = Jn/\sigma\f$.  The five
	  solves of the matrix
	  \f[
	      \begin{bmatrix}
                 J & s \\
		 n^T & 0
	      \end{bmatrix}
	  \f]
	  are computed via the applyBorderedJacobianInverse() method of the
	  underlying group.

	*/
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverse(NOX::Parameter::List& params, 
			     const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

	/*! 
	 * \brief Applies the inverse of the extended Jacobian matrix using 
	 * a modified bordering algorithm for multiple right-hand sides.
	 */
	/*!
	 * The current implemention solves each right hand side serially.
	 */
// 	virtual NOX::Abstract::Group::ReturnType 
// 	applyJacobianInverseMulti(
// 			     NOX::Parameter::List& params, 
// 			     const NOX::Abstract::Vector* const* inputs, 
// 			     NOX::Abstract::Vector** results, int nVecs) const;
      
      protected:
	
	//! Performs a deflated solve
	virtual NOX::Abstract::Group::ReturnType 
	deflatedJacobianInverse(NOX::Parameter::List& params, 
				const NOX::Abstract::Vector& u, 
				const NOX::Abstract::Vector& input, 
				NOX::Abstract::Vector& result,
				double& r) const;
      };
    } // namespace TPBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
