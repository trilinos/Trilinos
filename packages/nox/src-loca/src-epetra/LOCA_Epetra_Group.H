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

#ifndef LOCA_EPETRA_GROUP_H
#define LOCA_EPETRA_GROUP_H

#include "LOCA_Abstract_Group.H"    // base class
#include "NOX_Epetra_Group.H"       // base class
#include "LOCA_Parameter_Vector.H"  // class data element
#include "NOX_Common.H"             // class data element (string)

// Forward declares
namespace NOX {
  namespace Epetra {
    class SharedOperator;
  }
  namespace Parameter {
    class List;
  }
}
namespace LOCA {
  namespace Epetra {
    class Interface;
  }
}

namespace LOCA {
namespace Epetra {

//! Concrete implementation of LOCA::Abstract::Group for %Trilinos/%Epetra.
/*! This group is set up to use the linear algebra services provided through the Trilinos/Epetra package with AztecOO for the linear %solver. 
 */

class Group : public NOX::Epetra::Group, public LOCA::Abstract::Group {

public:
  //! Constructor with Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
   * for preconditioning.  An Epetra_Operator must be supplied for 
   * the Jacobian even in Matrix-Free mode. linearSolverParams is 
   * the "Linear Solver" sublist of parameter list.
   */
  Group(const NOX::Parameter::List& linearSolverParams, Interface& i, 
 	const ParameterVector& p, NOX::Epetra::Vector& x, Epetra_Operator& J);

  //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
  //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
  Group(const NOX::Parameter::List& linearSolverParams, Interface& i, 
 	const ParameterVector& p, NOX::Epetra::Vector& x, Epetra_Operator& J, 
	Epetra_Operator& M);

  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared Jacobian and shared preconditioning matrix. */
  Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual Group& operator=(const Group& source);
  //! See above.
  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
  //! See above.
  virtual LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);

  /** @name "Compute" functions. */
  //@{

  //! Set the bifurcation/continuation parameter list in the group.
  virtual void setParams(const ParameterVector& p);

  //! Set the bifurcation/continuation parameter list in the group.
  virtual void computeParams(const ParameterVector& oldParams,
			     const ParameterVector& direction, 
			     double step);

  virtual NOX::Abstract::Group::ReturnType 
  computeF();

  virtual NOX::Abstract::Group::ReturnType 
  computeJacobian();

  virtual NOX::Abstract::Group::ReturnType 
  computeTangent(NOX::Parameter::List& params, int paramId);

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  virtual const ParameterVector& getParams() const;

  virtual const NOX::Abstract::Vector& getTangent() const;

  //@}

  virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Return the userInterface.
  virtual NOX::Epetra::Interface& getUserInterface();

  virtual void print() const {}

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

protected:

  //! Solution vector pointer.
  ParameterVector* pVectorPtr;

  NOX::Epetra::Vector* tangentVecPtr;

  NOX::Epetra::Vector& tangentVec;

  bool isValidTangent;

  /** @name Additional Shared Operators */
  //@{
  //! Pointer to shared mass matrix 
  //NOX::Epetra::SharedOperator* sharedMassMatrixPtr;

  //@}

  //! Reference to the user supplied interface functions
  Interface& userInterface;

};

} // namespace Epetra
} // namespace LOCA


#endif
