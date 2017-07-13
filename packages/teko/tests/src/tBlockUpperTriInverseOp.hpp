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

#ifndef __tBlockUpperTriInverseOp_hpp__
#define __tBlockUpperTriInverseOp_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Epetra_SerialComm.h"

#include <string>
#include <vector>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tBlockUpperTriInverseOp : public UnitTest {
public:
   virtual ~tBlockUpperTriInverseOp() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_apply(int verbosity,std::ostream & os);
   bool test_alphabeta(int verbosity,std::ostream & os);
   bool test_applyTranspose(int verbosity,std::ostream & os);

protected:
   double tolerance_;
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > A_;
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > invA_;

   std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > invDiag_;
};

} // end namespace Tests
} // end namespace Teko

#endif
