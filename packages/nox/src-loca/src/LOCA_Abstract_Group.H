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

#ifndef LOCA_ABSTRACT_GROUP_H
#define LOCA_ABSTRACT_GROUP_H

#include "LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.H"  // base class
#include "LOCA_Bifurcation_TPBord_SingularSolveGroup.H"       // base class
#include "LOCA_Homotopy_AbstractGroup.H"                      // base class

/*! \brief Library of continuation algorithms namespace.
 */
/*! All classes, functions, and global variables written for %LOCA reside
in the LOCA namespace.
*/
namespace LOCA { 

  //! %LOCA abstract interface namespace
  namespace Abstract { 

    /*!  \brief %LOCA abstract interface to a group, derived from the 
      NOX::Abstract::Group.  %LOCA %abstract groups provide the interface 
      necessary to perform continuation, i.e., compute families of solutions to
      \f$ F(x,p) = 0 \f$.
    */
    /*! 
      Concrete implemenations of this interface must provide implementations of
      all of the methods in the NOX::Abstract::Group interface as well as the
      parameter and scale vector setting/retrieval methods.
    */

    class Group : 
      public LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup,
      public LOCA::Bifurcation::TPBord::SingularSolveGroup,
      public virtual LOCA::Homotopy::AbstractGroup
    {

    public:

      //! Default constructor.
      Group(const DerivUtils& deriv = LOCA::DerivUtils());
  
      //! Default constructor.
      Group(NOX::Parameter::List& params,
	    const DerivUtils& deriv = LOCA::DerivUtils());
 
      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual LOCA::Continuation::AbstractGroup& 
      operator=(const LOCA::Continuation::AbstractGroup& source);

      //! Assignment operator
      virtual LOCA::Continuation::FiniteDifferenceGroup& 
      operator=(const LOCA::Continuation::FiniteDifferenceGroup& source);

      //! Assignment operator
      virtual LOCA::Bifurcation::TPBord::AbstractGroup& 
      operator=(const LOCA::Bifurcation::TPBord::AbstractGroup& source);

      //! Assignment operator
      virtual LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& 
      operator=(const LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& source);

      //! Assignment operator
      virtual LOCA::Bifurcation::TPBord::SingularSolveGroup& 
      operator=(const LOCA::Bifurcation::TPBord::SingularSolveGroup& source);

      //! Assignment operator
      virtual LOCA::TimeDependent::AbstractGroup& 
      operator=(const LOCA::TimeDependent::AbstractGroup& source);

      //! Assignment operator
      virtual LOCA::Bifurcation::HopfBord::AbstractGroup& 
      operator=(const LOCA::Bifurcation::HopfBord::AbstractGroup& source);

      //! Assignment operator
      virtual LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& 
      operator=(const LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& source);

      //! Assignment operator
      virtual LOCA::Homotopy::AbstractGroup& 
      operator=(const LOCA::Homotopy::AbstractGroup& source);
  
      //! Assignment operator
      virtual Group& operator=(const Group& source);

      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double conParamValue);

      virtual NOX::Abstract::Group::ReturnType
      computeMassMatrix();

      virtual NOX::Abstract::Group::ReturnType
      applyMassMatrix(const NOX::Abstract::Vector& input,
		      NOX::Abstract::Vector& result) const;

      virtual NOX::Abstract::Group::ReturnType
      applyComplexInverse(NOX::Parameter::List& params,
			  const NOX::Abstract::Vector& input_real,
			  const NOX::Abstract::Vector& input_imag,
			  double frequency,
			  NOX::Abstract::Vector& result_real,
			  NOX::Abstract::Vector& result_imag) const;
    };
  } // namespace Abstract
} // namespace LOCA

#endif
