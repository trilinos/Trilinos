
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_PRECONDITIONER_JACOBI_H_
#define _EPETRA_PRECONDITIONER_JACOBI_H_

//#include "Epetra_Comm.h"
//#include "Epetra_Operator.h"	// base class

#include "NOX_Common.H"
//#include "Epetra_Import.h"
//#include "Epetra_Vector.h"

#include "NOX_Epetra_Preconditioner.H"

// Forward Declarations
//class Epetra_Comm;
//class Epetra_Map;

namespace NOX {

namespace Epetra {

/*! \brief Concrete implementation of the NOX::Epetra::Preconditioner class for Jacobi (Diagonal) Preconditioning.
 */
class JacobiPreconditioner : public Preconditioner {
      
 public:

  //! Constructor 
  JacobiPreconditioner(const Epetra_Vector& shape, const double value = 1.0e-6);

  //! Destructor
  virtual ~JacobiPreconditioner();

  virtual int SetUseTranspose(bool UseTranspose);

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  virtual double NormInf() const;
  
  virtual char* Label () const;
  
  virtual bool UseTranspose() const;

  virtual bool HasNormInf() const;

  virtual const Epetra_Comm & Comm() const;

  virtual const Epetra_BlockMap& DomainMap () const;

  virtual const Epetra_BlockMap& RangeMap () const;

  //! takes in the current solution vector and the Jacobian Epetra_Operator.  This function then computes the diagonal and stores it in an Epetra_Vector.
  virtual bool computePreconditioner(const Epetra_Vector& x, const Epetra_Operator* Jacobian = 0);

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
