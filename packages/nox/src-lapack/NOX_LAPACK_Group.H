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

#ifndef NOX_LAPACK_GROUP_H
#define NOX_LAPACK_GROUP_H

#include "NOX_Abstract_Group.H"	// base class

#include "NOX_Common.H"             // class data element (string)
#include "NOX_LAPACK_Vector.H"	    // class data element
#include "NOX_LAPACK_Matrix.H"	    // class data element
#include "NOX_LAPACK_Interface.H"   // class data element

// Forward declares
namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace NOX {
  namespace LAPACK {

    //! A simple example of a group structure, based on BLAS/LAPACK.
    class Group : public virtual NOX::Abstract::Group {

    public:

      //! Constructor
      Group(Interface& i);

      //! Constructor with used and allocated dimensions of matrices
      Group(Interface& i, int m, int n);

      //! Copy constructor
      Group(const NOX::LAPACK::Group& source, NOX::CopyType type = DeepCopy);

      //! Destructor.
      ~Group();


      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      //! See above.
      NOX::Abstract::Group& operator=(const NOX::LAPACK::Group& source);

      /** @name "Compute" functions. */
      //@{

      void setX(const NOX::Abstract::Vector& y);
      //! See above
      void setX(const NOX::LAPACK::Vector& y);

      void computeX(const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step);
      //! See above.
      void computeX(const NOX::LAPACK::Group& grp, const NOX::LAPACK::Vector& d, double step);

      NOX::Abstract::Group::ReturnType computeF();

      NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeGradient();

      NOX::Abstract::Group::ReturnType computeNewton(NOX::Parameter::List& params);

      //@}

      /** @name Jacobian operations.
       *
       * Operations using the Jacobian matrix. These may not be defined in
       * matrix-free scenarios. */

      //@{
  
      NOX::Abstract::Group::ReturnType applyJacobian(const NOX::LAPACK::Vector& input, 
						     NOX::LAPACK::Vector& result) const;

      //! See above
      NOX::Abstract::Group::ReturnType applyJacobian(const NOX::Abstract::Vector& input, 
						     NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianTranspose(const NOX::LAPACK::Vector& input, 
							      NOX::LAPACK::Vector& result) const;

      //! See above
      NOX::Abstract::Group::ReturnType applyJacobianTranspose(const NOX::Abstract::Vector& input, 
							      NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianInverse(NOX::Parameter::List& params, 
							    const NOX::LAPACK::Vector& input, 
				Vector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianInverse(NOX::Parameter::List& params, 
							    const NOX::Abstract::Vector& input, 
							    NOX::Abstract::Vector& result) const;
  
      //@}

      /** @name "Is" functions
       *
       * Checks to see if various objects have been computed. Returns true
       * if the corresponding "compute" function has been called since the
       * last update to the solution vector (via instantiation or
       * computeX). */
      //@{

      bool isF() const;
      bool isJacobian() const;
      bool isGradient() const;
      bool isNewton() const;

      //@}

      /** @name "Get" functions 
       *
       * Note that these function do not check whether or not the vectors
       * are valid. Must use the "Is" functions for that purpose. */
      //@{

      const NOX::Abstract::Vector& getX() const;

      const NOX::Abstract::Vector& getF() const;
  
      double getNormF() const;

      const NOX::Abstract::Vector& getGradient() const;

      const NOX::Abstract::Vector& getNewton() const;

      //@}

      virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Print out the group
      void print() const;

    protected:

      //! resets the isValid flags to false
      void resetIsValid();

    protected:

      /** @name Vectors */
      //@{
      //! Solution vector.
      NOX::LAPACK::Vector xVector;
      //! Right-hand-side vector (function evaluation).
      NOX::LAPACK::Vector fVector;
      //! Newton direction vector.
      NOX::LAPACK::Vector newtonVector;
      //! Gradient vector (steepest descent vector).
      NOX::LAPACK::Vector gradientVector;
      //@}


      //! Jacobian Matrix
      NOX::LAPACK::Matrix jacobianMatrix;

      //! LU factorization of Jacobian Matrix
      mutable NOX::LAPACK::Matrix jacobianLUFact;

      //! Pivot array for LU factorization
      mutable vector<int> pivots;

      //! Problem interface
      NOX::LAPACK::Interface& problemInterface;

      /** @name IsValid flags 
       *  
       * True if the current solution is up-to-date with respect to the
       * currect xVector. */
      //@{
      bool isValidF;
      bool isValidJacobian;
      bool isValidGradient;
      bool isValidNewton;
      mutable bool isValidJacobianLUFact;
      //@}
  
      //! Norm of F
      double normF;


    };

  } // namespace LAPACK
} // namespace NOX


#endif
