//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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
#include "BelosPseudoBlockCGSolMgr.hpp"

#include "Teuchos_ParameterList.hpp"

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
  EpetraOperator( const Teuchos::RCP<LinearProblem<double,Epetra_MultiVector,Epetra_Operator> >& lp, 
		  const Teuchos::RCP<Teuchos::ParameterList>& plist,
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

  Teuchos::RCP<SolverManager<double,Epetra_MultiVector,Epetra_Operator> > solver_;
  Teuchos::RCP<LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > lp_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;

  std::vector<char> Solver;
  bool initSolnVec_;
};

} //end namespace Belos

// end of file BELOS_EPETRA_OPERATOR_H
#endif 

