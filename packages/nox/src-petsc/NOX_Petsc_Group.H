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

#ifndef NOX_PETSC_GROUP_H
#define NOX_PETSC_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Petsc_Vector.H"	// class data element
#include "NOX_Common.H" // class data element (std::string)

#include "petscsnes.h"	
#include "petscksp.h"	

// Forward declares
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
namespace Petsc {
class Interface;
class SharedJacobian;
}
}

namespace NOX {
namespace Petsc {

//! Concrete implementation of NOX::Abstract::Group for %Petsc.
class Group : virtual public Abstract::Group {

public:

  //! Constructor.
  Group(Interface& i, Vec& x, Mat& J);

  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared Jacobian. */
  Group(const Group& source, CopyType type = DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual Abstract::Group& operator=(const NOX::Abstract::Group& source);
  //! See above.
  virtual Abstract::Group& operator=(const Group& source);

  /** @name "Compute" functions. */
  //@{

  virtual void setX(const Vector& y);

  //! See above
  virtual void setX(const Abstract::Vector& y);

  virtual void computeX(const Group& grp, 
			const Vector& d, 
			double step);
  //! See above.
  virtual void computeX(const Abstract::Group& grp,
			const Abstract::Vector& d, 
			double step);

  virtual Abstract::Group::ReturnType computeF();

  virtual Abstract::Group::ReturnType computeJacobian();

  virtual Abstract::Group::ReturnType computeGradient();

  virtual Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  virtual Abstract::Group::ReturnType 
  applyJacobian(const Vector& input, Vector& result) const;
  //! See above
  virtual Abstract::Group::ReturnType 
  applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const;

  virtual Abstract::Group::ReturnType 
  applyJacobianTranspose(const Vector& input, Vector& result) const;
  //! See above
  virtual Abstract::Group::ReturnType 
  applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  /*! 
      \brief Preconditions a vector using one of the available methods in %Petsc
 
      Currently, the preconditioning method corresponds to whatever
      is available in the user's %Petsc installation and is specified
      in the local \b .petscrc file.  For example, a direct solve can
      be accomplished by specifying \b -pc_type \b lu (serial only),
      whereas diagonal scaling can be achieved via \b -pc_type \b jacobi.
  */ 
  virtual Abstract::Group::ReturnType 
  applyRightPreconditioning(const Vector& input, Vector& result) const;
  //! See above
  virtual Abstract::Group::ReturnType 
  applyRightPreconditioning(Teuchos::ParameterList&, const Abstract::Vector& input, Abstract::Vector& result) const;
  
  //@}

  /** @name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */
  //@{

  virtual bool isF() const;
  virtual bool isJacobian() const;
  virtual bool isGradient() const;
  virtual bool isNewton() const;
  virtual bool isPreconditioner() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  virtual const Abstract::Vector& getX() const;

  virtual const Abstract::Vector& getF() const;
  
  virtual double getNormF() const;

  virtual const Abstract::Vector& getGradient() const;

  virtual const Abstract::Vector& getNewton() const;

  virtual Teuchos::RCP< const Abstract::Vector > getXPtr() const;

  virtual Teuchos::RCP< const Abstract::Vector > getFPtr() const;

  virtual Teuchos::RCP< const Abstract::Vector > getGradientPtr() const;

  virtual Teuchos::RCP< const Abstract::Vector > getNewtonPtr() const;

  //@}

  virtual Teuchos::RCP<NOX::Abstract::Group> 
  clone(CopyType type = DeepCopy) const;

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

protected:

  /** @name Vectors */
  //@{
  //! Solution vector.
  Vector xVector;
  //! Right-hand-side vector (function evaluation).
  Vector RHSVector;
  //! Gradient vector (steepest descent vector).
  Vector gradVector;
  //! Newton direction vector.
  Vector NewtonVector;
  //@}

  /** @name IsValid flags 
   *  
   * True if the current solution is up-to-date with respect to the
   * currect xVector. */
  //@{
  bool isValidRHS;
  bool isValidJacobian;
  bool isValidGrad;
  bool isValidNewton;
  bool isValidPreconditioner;
  //@}
  
  //! Norm of RHS
  double normRHS;

  //! Pointer to shared Jacobian matrix
  SharedJacobian* sharedJacobianPtr;

  //! Reference to shared Jacobian matrix 
  SharedJacobian& sharedJacobian;

  //! Specification of Jacobian Type
  std::string jacType;

  //! Reference to the user supplied interface functions
  Interface& userInterface;

};

} // namespace Petsc
} // namespace NOX



#endif
