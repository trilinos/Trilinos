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

#ifndef LOCA_EPETRA_INTERFACE_H
#define LOCA_EPETRA_INTERFACE_H

#include "NOX_Common.H"
#include "NOX_Epetra_Interface.H"
#include "NOX_Epetra_Vector.H"
  
// Forward declarations
class Epetra_Vector;
class Epetra_RowMatrix;
class Epetra_Operator;

namespace LOCA {
  class ParameterVector;
}

namespace LOCA {
  namespace Epetra {

    /*! \brief 
      Used by LOCA::Epetra::Group to provide a link to the
      external code for RHS and Jacobian fills ()if needed. 
    */
    /*!
      This interface is derived from the NOX::Epetra::Interface and 
      and additionally provides a method for setting problem parameters.
    */
    class Interface : public NOX::Epetra::Interface {

    public:
 
      //! Constructor
      Interface() {};

      //! Destructor
      virtual ~Interface() {};

      //! Set the parameters in the user's app prior to calling one of the compute functions.
      virtual void setParameters(const ParameterVector& p) = 0;

      //! Call user's own print routine for vector-parameter pair
      virtual void printSolution(const Epetra_Vector& x_, const double conParam) {}

      //! Projects solution to a few scalars for multiparameter continuation
      /*! Default implementation is the max norm 
       */
      virtual void projectToDraw(const NOX::Epetra::Vector& x,
				 double *px) const {
	px[0] = x.norm(NOX::Abstract::Vector::MaxNorm);
      }

      //! Returns the dimension of the project to draw array
      virtual int projectToDrawDimension() const { return 1; }

    };
  } // namespace Epetra
} // namespace LOCA

#endif
