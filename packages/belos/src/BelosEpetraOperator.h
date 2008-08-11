// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

/*! \file BelosEpetraOperator.h
    \brief This file provides an Epetra_Operator interface so Belos can be integrated into
     other codes as an abstract operator.
*/

#ifndef BELOS_EPETRA_OPERATOR_H
#define BELOS_EPETRA_OPERATOR_H

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"

#include "BelosLinearProblem.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOutputManager.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"

#include "Teuchos_ParameterList.hpp"
using Teuchos::RCP;
using Teuchos::ParameterList;

/*! \class Belos::EpetraOperator
    \brief This class provides and interface to the Epetra_Operator class, so Belos can be 
    integrated into other codes as an abstract operator.
*/

namespace Belos {

///////////////////////////////////////////////////////////////
//-------- class BelosEpetraOperator --------------------
//
// This class will allow Belos to be called as an Epetra_Operator.
// Thus, it can use itself as a preconditioner if need be.  It can
// also be used as the inner iteration of Anasazi :)
//
///////////////////////////////////////////////////////////////

class EpetraOperator : public virtual Epetra_Operator {
public:
  //! @name Constructor / Destructor
  //@{ 
  
  //! Constructor
  EpetraOperator( const RCP<LinearProblem<double,Epetra_MultiVector,Epetra_Operator> >& lp, 
		  const RCP<ParameterList>& plist,
                  bool initSolnVec = false );
  
  //! Destructor
  virtual ~EpetraOperator() {}
  //@}
  
  //! @name Attribute methods
  //@{ 
  
  //! Set whether the operator or its inverse should be applied. [ This option is not implemented ]
  int SetUseTranspose( bool UseTranspose_in ) { return(-1); };
  //@}
  
  //! @name Operator application methods
  //@{ 
  
  //! Apply the operator.
  int Apply( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;
  
  //! Apply the operator's inverse.
  int ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const;
  //@}
  
  //! @name Norm methods
  //@{ 
  
  //! Compute the infinity norm of the operator. [ This option is not implemented ]
  double NormInf() const { return(0.0); };
  //@}
  
  //! @name Attribute access functions
  //@{ 
  
  //! Return the label of the operator.
  const char* Label() const { return(&Solver[0]); };
  
  //! Return whether the operator is using the transpose.
  bool UseTranspose() const { return(false); };
  
  //! Return whether the infinity norm is available for this operator.
  bool HasNormInf() const { return(false); };
  
  //! Return the communicator for this operator.
  const Epetra_Comm& Comm() const;
  
  //! Return the domain map for this operator.
  const Epetra_Map& OperatorDomainMap() const;
  
  //! Return the range map for this operator.
  const Epetra_Map& OperatorRangeMap() const;	
  //@}	   
private:

  RCP<SolverManager<double,Epetra_MultiVector,Epetra_Operator> > solver_;
  RCP<LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > lp_;
  RCP<ParameterList> plist_;

  std::vector<char> Solver;
  bool initSolnVec_;
};

} //end namespace Belos

// end of file BELOS_EPETRA_OPERATOR_H
#endif 

