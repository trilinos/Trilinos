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

#ifndef _EPETRA_PRECONDITIONER_JACOBI_H_
#define _EPETRA_PRECONDITIONER_JACOBI_H_

//#include "Epetra_Comm.h"
//#include "Epetra_Operator.h"	// base class

#include "NOX_Common.H"
//#include "Epetra_Import.h"
//#include "Epetra_Vector.h"

#include "NOX_Epetra_Operator.H"

// Forward Declarations
//class Epetra_Comm;
//class Epetra_Map;

namespace NOX {

namespace Epetra {

/*! \brief Concrete implementation of the NOX::Epetra::Preconditioner class for Jacobi (Diagonal) Preconditioning.
 */
class JacobiPreconditioner : public Operator {

 public:

  //! Constructor
  JacobiPreconditioner(const Epetra_Vector& shape, const double value = 1.0e-6);

  //! Destructor
  virtual ~JacobiPreconditioner();

  virtual int SetUseTranspose(bool UseTranspose);

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual double NormInf() const;

  virtual const char* Label () const;

  virtual bool UseTranspose() const;

  virtual bool HasNormInf() const;

  virtual const Epetra_Comm & Comm() const;

  virtual const Epetra_Map& OperatorDomainMap () const;

  virtual const Epetra_Map& OperatorRangeMap () const;

  //! takes in the current solution vector and the Jacobian Epetra_Operator.  This function then computes the diagonal and stores it in an Epetra_Vector.
  virtual bool compute(const Epetra_Vector& x, const Epetra_Operator* Jacobian = 0);

private:

  //! Label for this operator
  string label;

  //! Pointer to the vector containing the diagonal of the Jacobian
  Epetra_Vector* diagonalVectorPtr;

  //! Reference to the vector containing the diagonal of the Jacobian
  Epetra_Vector& diagonalVector;

  //! The minimum diagonal value allowed for this preconditioner to be used.
  //! Defaults to 1.0e-6.
  double minValue;

};
} // namespace Epetra
} // namespace NOX

#endif /* _EPETRA_PRECONDITIONER_JACOBI_H_ */
