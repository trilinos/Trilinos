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

#ifndef NOX_EPETRANEW_LINEARSYSTEM_H
#define NOX_EPETRANEW_LINEARSYSTEM_H

#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Utils.H"          // class data element
#include "NOX_Common.H"         // class data element (string)

// Forward declares
namespace NOX {
  namespace EpetraNew {
    class Scaling;
  }
  namespace Parameter {
    class List;
  }
}
class Epetra_Vector;

namespace NOX {
namespace EpetraNew {

class LinearSystem {

protected:
 


public:
  //! Constructor.
  LinearSystem(){};

  //! Destructor.
  virtual ~LinearSystem(){};

  virtual bool applyJacobian(const NOX::Epetra::Vector& input, 
			     NOX::Epetra::Vector& result) const = 0;

  virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
				      NOX::Epetra::Vector& result) const = 0;
  
  virtual bool applyJacobianInverse(NOX::Parameter::List &params, 
				    const NOX::Epetra::Vector &input, 
				    NOX::Epetra::Vector &result) const = 0;

  virtual bool applyRightPreconditioning(bool useTranspose,
					 NOX::Parameter::List& params, 
					 const NOX::Epetra::Vector& input, 
					 NOX::Epetra::Vector& result) const = 0;

  //! Sets the diagonal scaling vector(s) used in scaling the linear system.  See NOX::EpetraNew::Scaling for details on how to specify scaling of the linear system.  
  virtual void resetScaling(NOX::EpetraNew::Scaling& scalingObject) = 0;

  virtual bool computeJacobian(Epetra_Vector& x) = 0;

  virtual bool createPreconditioner(Epetra_Vector& x, 
				    NOX::Parameter::List& p,
				    bool recomputeGraph) const = 0;

  virtual bool destroyPreconditioner() const = 0;
};

} // namespace EpetraNew
} // namespace NOX


#endif
