// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_BIFURCATION_TPBORD_ABSTRACTGROUP_H
#define LOCA_BIFURCATION_TPBORD_ABSTRACTGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class

namespace LOCA { 
  
  //! Groups and vectors for locating bifurcations
  namespace Bifurcation { 

    /*! 
     * \brief Groups and vectors for locating turning point bifurcations
     * using bordering techniques.
     */
    namespace TPBord {

      /*! 
       * \brief Interface to underlying groups for turning point calculations
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate turning points and pitchforks using the bordering
       * algorithms (see LOCA::Bifurcation::TPBord::ExtendedGroup and
       * LOCA::Bifurcation::PitchforkBord::ExtendedGroup for a description
       * of the governing equations and bordering algorithms).  
       *
       * This class
       * is derived from the LOCA::Continuation::AbstractGroup and declares
       * several pure virtual methods to solve the linear set of equations
       * \f$Jx=b\f$ for \f$J\f$ nearly singular, and to compute various 
       * derivatives of \f$Jn\f$ for a given \f$n\f$.  Default implementations
       * for the singular solves are encapsulated in the 
       * LOCA::Bifurcation::TPBord::SingularSolveGroup and for the derivatives
       * using finite differencing in the 
       * LOCA::Bifurcation::TPBord::FiniteDifferenceGroup.
       *
       * The last virtual method, innerProduct() has a default implementation
       * given by the dot product, but should be overloaded for any 
       * problem that has a different definition for the inner product (this
       * is only used in the LOCA::Bifurcation::PitchforkBord::ExtendedGroup
       * for the assymetry).
       */
      class AbstractGroup : 
	public virtual LOCA::Continuation::AbstractGroup {

      public:
  
	//! Default constructor.
	AbstractGroup() {}

	//! Destructor
	virtual ~AbstractGroup() {}

	/*! 
	 * @name Pure virtual methods 
	 * These methods must be defined by any concrete implementation
	 */
	//@{
  
	//! Assignment operator
	virtual AbstractGroup& operator=(const AbstractGroup& source) = 0;

	//! applyJacobianInverse for nearly singular systems
	/*!
	 * This method is to compute the solution to \f$Jx=b\f$ for nearly
	 * singular \f$J\f$.  \b approxNullVec contains an approximation to the
	 * null vector \f$n\f$ such that \f$Jn\approx0\f$, and 
	 * \b jacApproxNullVec
	 * contains \f$Jn\f$.  The idea is these quantities are supposed to
	 * help in obtaining an accurate solution to \f$Jx=b\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType
	applySingularJacobianInverse(
				NOX::Parameter::List& params,
				const NOX::Abstract::Vector& input,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) = 0;

	//! Version of applySingularJacobianInverse() for multiple RHS.
	virtual NOX::Abstract::Group::ReturnType 
	applySingularJacobianInverseMulti(
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector*const* inputs,
			       const NOX::Abstract::Vector& approxNullVec,
			       const NOX::Abstract::Vector& jacApproxNullVec,
			       NOX::Abstract::Vector** results,
			       int nVecs) = 0;

	/*! 
	 * \brief Computes the derivative \f$\partial Jn/\partial p\f$ 
	 * where \f$p\f$ is the parameter indexed by \b param_id.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDp(const NOX::Abstract::Vector& nullVector,
		     const int param_id, NOX::Abstract::Vector& result) = 0;

	/*! 
	 * \brief Computes the derivative \f$\partial Jn/\partial p\f$ 
	 * where \f$p\f$ is the parameter indexed by \b param_id.  The argument
	 * \b JnVector stores \f$Jn\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDp(const NOX::Abstract::Vector& nullVector,
		     const int param_id,
		     const NOX::Abstract::Vector& JnVector,
		     NOX::Abstract::Vector& result) = 0;

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::Vector& aVector,
		      NOX::Abstract::Vector& result) = 0;
  
	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. The argument \b JnVector stores \f$Jn\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::Vector& aVector,
		      const NOX::Abstract::Vector& JnVector,
		      NOX::Abstract::Vector& result) = 0;

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::MultiVector& aVector,
		      NOX::Abstract::MultiVector& result) = 0;
  
	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial Jn}{\partial x} a\f$ for the given direction 
	 * \f$a\f$. The argument \b JnVector stores \f$Jn\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDJnDxa(const NOX::Abstract::Vector& nullVector,
		      const NOX::Abstract::MultiVector& aVector,
		      const NOX::Abstract::Vector& JnVector,
		      NOX::Abstract::MultiVector& result) = 0;

	//@}

      
	/*! 
	 * @name Virtual methods with default implementations 
	 * These methods should be overloaded in a concrete implementation
	 * if more appropriate/efficient approaches are available.
	 */
	//@{

	//! Compute inner product between two vectors.  
	/*
	 * Default definition is x.dot(y).
	 */
	virtual double innerProduct(const NOX::Abstract::Vector& x,
				    const NOX::Abstract::Vector& y);

	//@}

	/*!
	 * @name Virtual methods specific to the modified bordering algorithm
	 */
	//@{

	/*!
	 * \brief Solve a bordered system of equations
	 */
	/*!
	 * Compute the solution to the bordered system of equations
	 * \f[
	 *     \begin{bmatrix}
	 *        J & a \\
	 *        b^T & 0
	 *     \end{bmatrix}
	 *     \begin{bmatrix}
	 *        \text{vResult} \\
	 *        \text{sResult}
	 *     \end{bmatrix} = 
	 *     \begin{bmatrix}
	 *        \text{vInput} \\
	 *        \text{sInput}
	 *     \end{bmatrix}
	 * \f]
	 */
	virtual NOX::Abstract::Group::ReturnType
	applyBorderedJacobianInverse(bool trans,
				     NOX::Parameter::List& params,
				     const NOX::Abstract::Vector& a,
				     const NOX::Abstract::Vector& b,
				     const NOX::Abstract::Vector& vInput,
				     double sInput,
				     NOX::Abstract::Vector& vResult,
				     double& sResult) const = 0;

	/*!
	 * \brief Solve a bordered system of equations with multiple right-
	 * hand-sides.
	 */
	/*!
	 * The default implementation here is to call 
	 * applyBorderedJacobianInverse() serially for each right-hand-side.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	applyBorderedJacobianInverseMulti(
			       bool trans,
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector& a,
			       const NOX::Abstract::Vector& b,
			       const NOX::Abstract::Vector*const* vInputs,
			       double *sInputs,
			       NOX::Abstract::Vector** vResults,
			       double *sResults,
			       int nVecs) const;

	//@}
      };
    } // namespace TPBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
