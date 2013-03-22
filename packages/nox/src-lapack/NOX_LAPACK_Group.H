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

#ifndef NOX_LAPACK_GROUP_H
#define NOX_LAPACK_GROUP_H

#include "NOX_Abstract_Group.H"	// base class

#include "NOX_Common.H"             // class data element (std::string)
#include "NOX_LAPACK_Vector.H"	    // class data element
#include "NOX_LAPACK_LinearSolver.H"	    // class data element
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

      NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);

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

      NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& params, 
							    const NOX::LAPACK::Vector& input, 
				Vector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& params, 
							    const NOX::Abstract::Vector& input, 
							    NOX::Abstract::Vector& result) const;

      //! applyJacobianInverse for multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;
  
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

      Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;
  
      Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;

      //@}

      virtual Teuchos::RCP<NOX::Abstract::Group> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

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


      //! Jacobian Matrix & LAPACK solver
      mutable NOX::LAPACK::LinearSolver<double> jacSolver;

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
      //@}

    };

  } // namespace LAPACK
} // namespace NOX


#endif
