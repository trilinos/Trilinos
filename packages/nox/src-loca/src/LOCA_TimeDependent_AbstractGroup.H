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

#ifndef LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H
#define LOCA_TIMEDEPENDENT_ABSTRACTGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class

namespace LOCA { 

  namespace TimeDependent {

    /*!  
     * \brief Interface to underlying groups for time dependent systems
     */
    /*!
     * This abstract class provides an interface for time dependent problems,
     * i.e., problems with a mass matrix (typically used in eignvalue or
     * Hopf calculations).  It provides two pure virtual methods, one to
     * compute the mass matrix and the second to apply the mass matrix to 
     * a vector.
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
      virtual AbstractGroup& 
      operator=(const AbstractGroup& source) = 0;
      
      //! Compute the mass matrix of the system
      virtual NOX::Abstract::Group::ReturnType
      computeMassMatrix() = 0;

      //! Multiply the mass matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyMassMatrix(const NOX::Abstract::Vector& input,
		      NOX::Abstract::Vector& result) const = 0;

      //! Is mass matrix valid
      virtual bool isMassMatrix() const = 0;

      //! Compute the shifted matrix, used by spectral transformations
      //! in linear stability analysis capability.
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(const double& shift)
      {  cout << "Error:  No implementation of computeShiftedMatrix in LOCA::TimeDependent::AbstractGroup" << endl;
         return NOX::Abstract::Group::BadDependency;};

      //! Multiply the shifted matrix by a vector
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
                         NOX::Abstract::Vector& result,
                         const double& shift, bool massMatrix) const
      { cout << "Error:  No implementation of applyShiftedMatrix in LOCA::TimeDependent::AbstractGroup" << endl;
        return NOX::Abstract::Group::BadDependency;};

      //! Apply the inverse of the shifted matrix by a vector, as needed by
      //! the shift-and-invert and generaized Cayley transformations for
      //! using the eigensolver as a linear stability analysis capability.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverse(const NOX::Abstract::Vector& input,
                                NOX::Abstract::Vector& result,
                                const double& shift, bool massMatrix,
                                NOX::Parameter::List& params) const
      { cout << "Error:  No implementation of applyShiftedMatrixInverse in LOCA::TimeDependent::AbstractGroup" << endl;
        return NOX::Abstract::Group::BadDependency;};


      //@}

    };
  } // namespace TimeDependent
} // namespace LOCA

#endif
