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
    /*! 
      This group is set up to use the linear algebra services provided 
      through the Trilinos/Epetra package with AztecOO for the linear %solver. 

      This class is derived both from the NOX::Epetra::Group and
      LOCA::Abstract::Group classes and therefore inherits the implementation
      of the NOX::Abstract::Group interface provided by NOX::Epetra::Group.
      
      This class additionally stores a ParameterVector params and a scaling
      vector (NOX::Epetra::Vector) scaleVec and implements the pure virtual
      methods of the LOCA::Abstract::Group class.  

      The computeF() and computeJacobian() methods of the NOX::Epetra::Group
      parent class are overloaded here.  They both set the entire contents
      of the parameter vector in the problem interface before calling the
      NOX::Epetra::Group computeF() and computeJacobian().
    */

    class Group : public NOX::Epetra::Group, public LOCA::Abstract::Group {

    public:
      //! Constructor with Jacobian Operator only.  
      /*! Either there is no preconditioning or the Jacobian will be used
       * for preconditioning.  An Epetra_Operator must be supplied for 
       * the Jacobian even in Matrix-Free mode. linearSolverParams is 
       * the "Linear Solver" sublist of parameter list.
       *
       * The scaling vector is set to a default of all one's
       */
      Group(NOX::Parameter::List& printParams, 
	    NOX::Parameter::List& linearSolverParams, Interface& i, 
	    const ParameterVector& p, NOX::Epetra::Vector& x, 
	    Epetra_Operator& J);

      //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
      /*! 
	linearSolverParams is the "Linear %Solver" sublist of parameter list.
	The scaling vector is set to a default of all one's
      */
      Group(NOX::Parameter::List& printParams, 
	    NOX::Parameter::List& linearSolverParams, Interface& i, 
	    const ParameterVector& p, NOX::Epetra::Vector& x, 
	    Epetra_Operator& J, 
	    Epetra_Operator& M);

      //! Constructor with Jacobian Operator only.  
      /*! Either there is no preconditioning or the Jacobian will be used
       * for preconditioning.  An Epetra_Operator must be supplied for 
       * the Jacobian even in Matrix-Free mode. linearSolverParams is 
       * the "Linear Solver" sublist of parameter list.
       */
      Group(NOX::Parameter::List& printParams, 
	    NOX::Parameter::List& linearSolverParams, Interface& i, 
	    const ParameterVector& p, NOX::Epetra::Vector& x, 
	    Epetra_Operator& J,
	    const NOX::Epetra::Vector& s);

      
      //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
      //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
      Group(NOX::Parameter::List& printParams, 
	    NOX::Parameter::List& linearSolverParams, Interface& i, 
	    const ParameterVector& p, NOX::Epetra::Vector& x, 
	    Epetra_Operator& J, 
	    Epetra_Operator& M, const NOX::Epetra::Vector& s);

      /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
	valid shared Jacobian and shared preconditioning matrix. */
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

      //! See above.
      virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

      //! See above.
      virtual LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);

      //! Cloning function
      virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Overloaded computeF()
      /*!
	Calls LOCA::Epetra::Interface::setParams before evalulating F.
      */
      virtual NOX::Abstract::Group::ReturnType 
      computeF();

      //! Overloaded computeJacobian()
      /*!
	Calls LOCA::Epetra::Interface::setParams before evalulating J.
      */
      virtual NOX::Abstract::Group::ReturnType 
      computeJacobian();

      //! Set the parameters
      virtual void setParams(const ParameterVector& p);

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);
  
      //! Set parameter indexed by paramID
      virtual void setParam(string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(string paramID) const;

      //! Return the userInterface.
      virtual NOX::Epetra::Interface& getUserInterface();

      virtual void print() const {}

      //! Set the scale vector
      virtual void setScaleVec(const NOX::Abstract::Vector& s);

      //! Set the scale vector
      void setScaleVec(const NOX::Epetra::Vector& s);

      //! Get the scale vector
      virtual const NOX::Abstract::Vector& getScaleVec() const;

    protected:
  
      //! Compute a default scale vector (vector of all one's)
      void computeScaleVec();

    protected:

      //! Parameter vector
      LOCA::ParameterVector params;

      /** @name Additional Shared Operators */
      //@{
      //! Pointer to shared mass matrix 
      //NOX::Epetra::SharedOperator* sharedMassMatrixPtr;

      //@}

      //! Reference to the user supplied interface functions
      Interface& userInterface;

      //! Scale vector
      NOX::Epetra::Vector scaleVec;

    };

  } // namespace Epetra
} // namespace LOCA


#endif
