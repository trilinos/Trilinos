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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

// Teko-Package includes
#include "Teko_Utilities.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_SpmdVectorBase.hpp"

// Test-rig

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

TEUCHOS_UNIT_TEST(tUtilitiesTests, clipping)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   RCP<Thyra::VectorSpaceBase<double> > vs = Thyra::defaultSpmdVectorSpace<double>(10);
   RCP<Thyra::MultiVectorBase<double> > x = Thyra::createMembers<double>(vs,2);

   Thyra::assign(x.ptr(),1.0);

   // try to clip lower values
   Teko::clipLower(x,2.0);

   TEST_FLOATING_EQUALITY(Teko::norm_1(x,0),10.0*2.0,1e-16);
   TEST_FLOATING_EQUALITY(Teko::norm_1(x,1),10.0*2.0,1e-16);

   Teko::clipUpper(x,1.0);

   TEST_FLOATING_EQUALITY(Teko::norm_1(x,0),10.0*1.0,1e-16);
   TEST_FLOATING_EQUALITY(Teko::norm_1(x,1),10.0*1.0,1e-16);

   Teuchos::ArrayRCP<double> col0, col1; 
   rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(0))->getNonconstLocalData(Teuchos::ptrFromRef(col0));
   rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(1))->getNonconstLocalData(Teuchos::ptrFromRef(col1));

   TEST_EQUALITY(col0.size(),10);
   TEST_EQUALITY(col1.size(),10);

   for(int i=0;i<10;i++) {
      col0[i] = i-5;
      col1[i] = 5-i;
   }

   TEST_FLOATING_EQUALITY(Teko::norm_1(x,0),25.0,1e-16);
   TEST_FLOATING_EQUALITY(Teko::norm_1(x,1),25.0,1e-16);

   Teko::clipLower(x,0.0);

   TEST_FLOATING_EQUALITY(Teko::norm_1(x,0),10.0,1e-16);
   TEST_FLOATING_EQUALITY(Teko::norm_1(x,1),15.0,1e-16);

   for(int i=0;i<10;i++) {
      col0[i] = i-5;
      col1[i] = 5-i;
   }

   TEST_FLOATING_EQUALITY(Teko::norm_1(x,0),25.0,1e-16);
   TEST_FLOATING_EQUALITY(Teko::norm_1(x,1),25.0,1e-16);

   Teko::clipUpper(x,0.0);

   TEST_FLOATING_EQUALITY(Teko::norm_1(x,0),15.0,1e-16);
   TEST_FLOATING_EQUALITY(Teko::norm_1(x,1),10.0,1e-16);
}

TEUCHOS_UNIT_TEST(tUtilitiesTests, replaceValues)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   RCP<Thyra::VectorSpaceBase<double> > vs = Thyra::defaultSpmdVectorSpace<double>(10);
   RCP<Thyra::MultiVectorBase<double> > x = Thyra::createMembers<double>(vs,2);

   Teuchos::ArrayRCP<double> col0, col1; 
   rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(0))->getNonconstLocalData(Teuchos::ptrFromRef(col0));
   rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(1))->getNonconstLocalData(Teuchos::ptrFromRef(col1));

   for(int i=0;i<10;i++) {
      col0[i] =  i;
      col1[i] = -i;
   }

   Teko::replaceValue(x,0.0,99.0);

   TEST_EQUALITY(col0[0],99.0);
   TEST_EQUALITY(col1[0],99.0);
   for(int i=1;i<10;i++) {
      TEST_EQUALITY(col0[i],double(i));
      TEST_EQUALITY(col1[i],double(-i));
   }
}

TEUCHOS_UNIT_TEST(tUtilitiesTests, averages)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   RCP<Thyra::VectorSpaceBase<double> > vs = Thyra::defaultSpmdVectorSpace<double>(9);
   RCP<Thyra::MultiVectorBase<double> > x = Thyra::createMembers<double>(vs,2);

   Teuchos::ArrayRCP<double> col0, col1; 
   rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(0))->getNonconstLocalData(Teuchos::ptrFromRef(col0));
   rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(1))->getNonconstLocalData(Teuchos::ptrFromRef(col1));

   for(int i=0;i<9;i++) {
      col0[i] =  2.3+(i-4.0);
      col1[i] = -3.7+(i-4.0);
   }

   std::vector<double> averages;
   Teko::columnAverages(x,averages);

   TEST_EQUALITY(averages.size(),2);
   TEST_FLOATING_EQUALITY(averages[0], 2.3,1e-14);
   TEST_FLOATING_EQUALITY(averages[1],-3.7,1e-14);

   double avg = Teko::average(x);
   TEST_FLOATING_EQUALITY(avg,(9.0*2.3-9.0*3.7)/18.0,1e-14);
}
