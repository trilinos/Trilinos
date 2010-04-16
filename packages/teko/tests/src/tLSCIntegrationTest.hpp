/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER
*/

#ifndef __tLSCIntegrationTest_hpp__
#define __tLSCIntegrationTest_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tLSCIntegrationTest : public UnitTest {
public:
   virtual ~tLSCIntegrationTest() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return true; }

   bool test_withmassStable(int verbosity,std::ostream & os);
   bool test_nomassStable(int verbosity,std::ostream & os);
   bool test_plConstruction(int verbosity,std::ostream & os);

protected:
   void loadStableSystem();
   void solveList(Teuchos::ParameterList & paramList,int vcycles);

   double tolerance_;

   Teuchos::RCP<const Epetra_Map> velMap_;  // map of velocity space
   Teuchos::RCP<const Epetra_Map> prsMap_;  // map of pressure space
   Teuchos::RCP<const Epetra_Map> fullMap_; // map of pressure space

   // stable discretizations matrices
   Teuchos::RCP<const Epetra_CrsMatrix> sF_;
   Teuchos::RCP<const Epetra_CrsMatrix> sB_;
   Teuchos::RCP<const Epetra_CrsMatrix> sBt_;
   Teuchos::RCP<const Epetra_CrsMatrix> sQu_;
   Teuchos::RCP<Epetra_Operator> sA_;

   // stable rhs and IFISS solution
   Teuchos::RCP<Epetra_Vector> rhs_;
   Teuchos::RCP<const Epetra_Vector> sExact_;
};

} // end namespace Tests
} // end namespace Teko

#endif
