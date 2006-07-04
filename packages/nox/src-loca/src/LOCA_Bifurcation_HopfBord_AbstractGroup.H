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

#ifndef LOCA_BIFURCATION_HOPFBORD_ABSTRACTGROUP_H
#define LOCA_BIFURCATION_HOPFBORD_ABSTRACTGROUP_H

#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"    // base class
#include "LOCA_TimeDependent_AbstractGroup.H"    // base class

namespace LOCA { 
  
  namespace Bifurcation { 

    /*! 
     * \brief Groups and vectors for locating Hopf bifurcations
     * using bordering techniques.
     */
    namespace HopfBord {

      /*! 
       * \brief Interface to underlying groups for Hopf point calculations
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate Hopf bifurcations using the bordering
       * algorithms (see LOCA::Bifurcation::HopfBord::ExtendedGroup for a 
       * description of the governing equations and bordering algorithms).  
       *
       * This class
       * is derived from the LOCA::Bifurcation::TPBord::AbstractGroup and
       * LOCA::TimeDependent::AbstractGroup and declares several pure virtual
       * methods for applying, solving, and computing derivatives of the
       * complex matrix \f$J+i\omega B\f$ where \f$J\f$ is the Jacobian matrix,
       * \f$B\f$ is the mass matrix, and \f$\omega\f$ is a (real) scalar.
       *
       * The last virtual method, applyComplexInverseMulti() for solving 
       * systems with multiple right-hand-sides, has a default
       * implementation to serially call applyComplexInverse() for each
       * right-hand-side.  This should be overloaded for problems where
       * a block-complex solver is available.
       */
      class AbstractGroup : 
	public virtual LOCA::Bifurcation::TPBord::AbstractGroup,
	public virtual LOCA::TimeDependent::AbstractGroup {

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
	virtual AbstractGroup& 
	operator=(const AbstractGroup& source) = 0;

	//! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
	/*!
	 * The argument \b frequency stores \f$\omega\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType
	applyComplexInverse(Teuchos::ParameterList& params,
			    const NOX::Abstract::Vector& input_real,
			    const NOX::Abstract::Vector& input_imag,
			    double frequency,
			    NOX::Abstract::Vector& result_real,
			    NOX::Abstract::Vector& result_imag) const = 0;


	/*! 
	 * \brief Computes the derivative 
	 * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial p}\f$ 
	 * where \f$p\f$ is the parameter indexed by \b param_id.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDp(const NOX::Abstract::Vector& yVector,
		     const NOX::Abstract::Vector& zVector,
		     double w,
		     const int param_id, 
		     NOX::Abstract::Vector& result_real,
		     NOX::Abstract::Vector& result_imag) = 0;

	/*! 
	 * \brief Computes the derivative 
	 * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial p}\f$ 
	 * where \f$p\f$ is the parameter indexed by \b param_id.  The 
	 * arguments \b Ce_real and \b Ce_imag hold the real and imaginary
	 * components of \f$(J+i\omega B)(y+iz)\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDp(const NOX::Abstract::Vector& yVector,
		     const NOX::Abstract::Vector& zVector,
		     double w,
		     const int param_id, 
		     const NOX::Abstract::Vector& Ce_real,
		     const NOX::Abstract::Vector& Ce_imag,
		     NOX::Abstract::Vector& result_real,
		     NOX::Abstract::Vector& result_imag) = 0;

	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial x} a\f$ for the 
	 * given direction \f$a\f$. 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDxa(const NOX::Abstract::Vector& yVector,
		      const NOX::Abstract::Vector& zVector,
		      double w,
		      const NOX::Abstract::Vector& aVector,
		      NOX::Abstract::Vector& result_real,
		      NOX::Abstract::Vector& result_imag) = 0;
  
	/*! 
	 * \brief Computes the directional derivative 
	 * \f$\frac{\partial (J+i\omega B)(y+iz)}{\partial x} a\f$ for the 
	 * given direction \f$a\f$. The 
	 * arguments \b Ce_real and \b Ce_imag hold the real and imaginary
	 * components of \f$(J+i\omega B)(y+iz)\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDCeDxa(const NOX::Abstract::Vector& yVector,
		      const NOX::Abstract::Vector& zVector,
		      double w,
		      const NOX::Abstract::Vector& aVector,
		      const NOX::Abstract::Vector& Ce_real,
		      const NOX::Abstract::Vector& Ce_imag,
		      NOX::Abstract::Vector& result_real,
		      NOX::Abstract::Vector& result_imag) = 0;

	//@}

      
	/*! 
	 * @name Virtual methods with default implementations 
	 * These methods should be overloaded in a concrete implementation
	 * if more appropriate/efficient approaches are available.
	 */
	//@{

	//! Compute \f$(J+i\omega B)(y+iz)\f$
	/*! 
	 * The default implementation computes the complex product using
	 * the Jacobian and Mass matrix apply methods.
	 */
	virtual NOX::Abstract::Group::ReturnType
	applyComplex(const NOX::Abstract::Vector& input_real,
		     const NOX::Abstract::Vector& input_imag,
		     double frequency,
		     NOX::Abstract::Vector& result_real,
		     NOX::Abstract::Vector& result_imag) const;

	//! Solve Complex system for multiple RHS
	/*!
	 * The default implementation here is to call applyComplexInverse()
	 * for each right-hand-side.  This method should be overloaded if
	 * a block-complex solver is available.
	 */
	virtual NOX::Abstract::Group::ReturnType
	applyComplexInverseMulti(
			     Teuchos::ParameterList& params,
			     const NOX::Abstract::Vector* const* inputs_real,
			     const NOX::Abstract::Vector* const* inputs_imag,
			     double frequency,
			     NOX::Abstract::Vector** results_real,
			     NOX::Abstract::Vector** results_imag,
			     int nVecs) const;

	//@}

      };
    } // namespace HopfBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
