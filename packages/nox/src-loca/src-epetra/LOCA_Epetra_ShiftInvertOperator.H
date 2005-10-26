
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef LOCA_EPETRA_SHIFTINVERTOPERATOR_H
#define LOCA_EPETRA_SHIFTINVERTOPERATOR_H

#include "Epetra_Operator.h"
#include "LOCA_Epetra_Group.H"
#include "LOCA_ErrorCheck.H"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_Interface_Jacobian.H"

class Epetra_MultiVector;
class Epetra_CrsMatrix;
class Epetra_Map;
class Epetra_Import;
class Epetra_BlockMap;
class Epetra_Comm;
class Epetra_Vector;

namespace LOCA {

  namespace Epetra {

    //! Interface for LOCA::Epetra::ShifterInvertOperator
    class ShiftInvertInterface : public NOX::Epetra::Interface::Jacobian {

    public:
      
      //! Constructor
      ShiftInvertInterface() {};
      
      //! Destructor
      ~ShiftInvertInterface() {};
      
      //! Compute Jacobian \f$J\f$.
      bool computeJacobian(const Epetra_Vector &x,
			   Epetra_Operator& Jac) {return true;}; 
    };

    //! Epetra operator for \f$(J-\sigma M)^{-1}\f$.
    class ShiftInvertOperator : public Epetra_Operator {
      
    public:

      //! Constructor
      ShiftInvertOperator(
		const Teuchos::RefCountPtr<const LOCA::Epetra::Group>& grp, 
		const Teuchos::RefCountPtr<const Epetra_Operator>& jac, 
		const double& shift, 
		bool hasMassMatrix);

      //! Destructor
      virtual ~ShiftInvertOperator();

      //! Set transpose
      virtual int SetUseTranspose(bool UseTranspose);
 
      //! Apply shifted operator
      virtual int Apply(const Epetra_MultiVector& X, 
			Epetra_MultiVector& Y) const;

      //! Apply shifted operator inverse
      virtual int ApplyInverse(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const; 
  
      //! Computing infinity norm
      virtual double NormInf() const;

      //! Label
      virtual const char * Label() const;

      //! Transpose
      virtual bool UseTranspose() const;

      //! Have norm-inf
      virtual bool HasNormInf() const;
  
      /*! 
       * \brief Returns a pointer to the Epetra_Comm communicator associated 
       * with this operator.
       */
      virtual const Epetra_Comm & Comm() const;
  
      /*! 
       * \brief Returns the Epetra_Map object associated with the domain of 
       * this operator.
       */
      virtual const Epetra_Map & OperatorDomainMap() const;
  
      /*! 
       * \brief Returns the Epetra_Map object associated with the range of 
       * this operator.
       */
      virtual const Epetra_Map & OperatorRangeMap() const;
 
    protected:  
  
      //! LOCA group
      Teuchos::RefCountPtr<const LOCA::Epetra::Group> locagrp; 

      //! Jacobian operator
      Teuchos::RefCountPtr<const Epetra_Operator> jacOper;

      //! Shift value
      const double& shift_;

      //! Do we have a mass matrix?
      bool massMatrix;

      //! Label
      const char * Label_;
    };
  } // namespace Epetra
} // namespace LOCA
#endif
