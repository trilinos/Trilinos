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

  //! %LOCA %Epetra support
  namespace Epetra {

    //! Extension of the NOX::Epetra::Group to %LOCA.  
    /*!
      This class is derived both from the NOX::Epetra::Group and
      LOCA::Abstract::Group classes and therefore inherits the implementation
      of the NOX::Abstract::Group interface provided by NOX::Epetra::Group.

      This class provides implementations of %LOCA AbstractGroup virtual 
      methods specific to the %Epetra group.  It stores a parameter vector
      for setting/retrieving parameter values 
      (LOCA::Continuation::AbstractGroup), provides a facility for computing
      eigenvalues (LOCA::Continuation::AbstractGroup) using the Anasazi 
      package, and augements the Jacobian matrix for homotopy
      (LOCA::Homotopy::AbstractGroup), It currently has no support for
      storing a mass matrix and therefore cannot be used in Hopf or 
      generalized eigenvalue computations.  Since it is
      derived from the LOCA::Abstract::Group (which is in-turn derived
      from all FiniteDifference groups), this group implicitly uses the
      finite-difference implementations of parameter and second derivatives
      provided by the FiniteDifference groups.  

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

      /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
	valid shared Jacobian and shared preconditioning matrix. */
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

      /*! 
       * @name Overloaded NOX::Epetra::Group  methods.
       */
      //@{

      //! Assignment operator.
      virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

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

      //@}

      /*! 
       * @name Implementation of LOCA::Abstract::Group virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);

      //@}

      /*! 
       * @name Implementation of LOCA::Continuation::AbstractGroup virtual methods.
       */
      //@{
      

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

      //! Call the user interface print() routine, solution vector
      virtual void printSolution(const double conParam) const;

      //! Call the user interface print() routine, any vector
      virtual void printSolution(const NOX::Abstract::Vector& x,
                                 const double conParam) const;

      //! Computes leading eigenvalues/eigenvectors using Anasazi package
      /*!
       * To use Anasazi for computing eigenvalues, NOX/LOCA must be
       * configured with --enable-loca-anasazi
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(NOX::Parameter::List& params);

      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::TPBord::AbstractGroup virtual methods.
       */
      //@{

      //@}

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
				   double& sResult) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
       */
      //@{

      /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$J+pI\f$ where \f$I\f$ is the 
       * identity matrix and \f$p\f$ is a scalar (\b conParamValue).
       */
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double conParamValue);

      //@}

      //! Return the userInterface.
      virtual NOX::Epetra::Interface& getUserInterface();

      //! Call the user interface print() routine, any vector
      virtual void printSolution(const NOX::Epetra::Vector& x,
                                 const double conParam) const;

      protected:

      //! Sets aztec solver options for matrix-free bordered system
      void setBorderedAztecOptions(const NOX::Parameter::List& p, 
				   AztecOO& aztec) const;

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

      //! Extra vector needed for intermediate calculations of LOCA routines.
      /*! NOTE: there already is a tmpVectorPtr in the NOX::Epetra::Group.  This is a second temporary vector if that one extra isn't enough.
       */
      Epetra_Vector* tmpVectorPtr2;

      //! Step counter for option of computing eigenvalues every n continuation steps
      int eigenvalCounter;
    };

  } // namespace Epetra
} // namespace LOCA


#endif
