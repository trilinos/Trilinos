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

#ifndef __test_utils_hpp__
#define __test_utils_hpp__

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Teuchos_Comm.hpp"

#include "Tpetra_Vector.hpp"
#include "Teko_Utilities.hpp"

#include <iostream>
#include <list>

namespace Teko {
namespace Test {

// build a 2x2 matrix...only in serial
const Teuchos::RCP<const Thyra::LinearOpBase<double> > build2x2(const Epetra_Comm & comm,double a,double b,double c,double d);
const Teuchos::RCP<const Thyra::LinearOpBase<ST> > build2x2(const Teuchos::RCP<const Teuchos::Comm<int> > comm,ST a,ST b,ST c,ST d);

// prints a vector, with string "s" as the name
void Print(std::ostream & os,const std::string & s,const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & v);

// prints a vector, with string "s" as the name
inline std::string Print(const std::string & s,const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & v)
{ std::stringstream ss; Print(ss,s,v); return ss.str(); }

// builds a single vector of the identity matrix
Teuchos::RCP<Thyra::VectorBase<double> > BuildIVector(int j,const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & vs);

// peforms an extraction of the matrix, in addition to printing it
// DANGER! This is a CPU intense routine, should only be done on small, in-core matricies
void HardExtract(std::ostream & os,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A);

// add two Thyra vectors
const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & x,
                                                   const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & y);
const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(double ax,const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & x,
                                                   double ay,const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & y);

// compute ||x-y||_2
double Difference(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & x,
                  const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & y);

// construct a diagonal matrix
const Teuchos::RCP<const Thyra::LinearOpBase<double> > DiagMatrix(int cnt,double * vec,std::string label="");

const Teuchos::RCP<const Thyra::LinearOpBase<ST> > DiagMatrix_tpetra(GO cnt,ST * vec,std::string label="");

// 2-Vector
const Teuchos::RCP<const Thyra::MultiVectorBase<double> > BlockVector(const Epetra_Vector & u, const Epetra_Vector & v,
        const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & vs);

const Teuchos::RCP<const Thyra::MultiVectorBase<ST> > BlockVector(const Tpetra::Vector<ST,LO,GO,NT> & u, const Tpetra::Vector<ST,LO,GO,NT> & v,
        const Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > & vs);

class UnitTest {
public:
   virtual ~UnitTest() {}
   virtual void initializeTest() = 0;
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun) = 0;
   virtual bool isParallel() const = 0;

   static void AddTest(const Teuchos::RCP<UnitTest> & ut,const std::string & name);
   static bool RunTests(int verbosity, std::ostream & stdstrm,std::ostream & failstrm);
   static bool RunTests_tpetra(int verbosity, std::ostream & stdstrm,std::ostream & failstrm);

   static Teuchos::RCP<const Epetra_Comm> GetComm();
   static Teuchos::RCP<const Teuchos::Comm<int> > GetComm_tpetra();
   static void SetComm(const Teuchos::RCP<const Epetra_Comm > & c);
   static void SetComm_tpetra(const Teuchos::RCP<const Teuchos::Comm<int> > & c);

   static void ClearTests();

protected:
   static std::list<std::pair<Teuchos::RCP<UnitTest>,std::string> > testList;
   static Teuchos::RCP<const Epetra_Comm > comm;
   static Teuchos::RCP<const Teuchos::Comm<int> > comm_tpetra;

   static bool CheckParallelBools(bool myBool,int & failPID);
   static bool CheckParallelBools_tpetra(bool myBool,int & failPID);
};

inline const std::string toString(bool status) { return status ? "PASSED" : "FAILED"; }

} // end namespace Tests
} // end namespace Teko

#define Teko_ADD_UNIT_TEST(str,name) Teko::Test::UnitTest::AddTest(Teuchos::rcp(new str()),#name)
#define Teko_TEST_MSG(os,level,msgp,msgf) {          \
    int failPID = -1;                              \
    status = UnitTest::CheckParallelBools(status,failPID); \
    if(verbosity>=level && status)                 \
       os << msgp << std::endl;                    \
    else if(verbosity>=level && not status)        \
       os << msgf << " ( PID = " << failPID << " )" << std::endl; \
    }

#define Teko_TEST_MSG_tpetra(os,level,msgp,msgf) {          \
    int failPID = -1;                              \
    status = UnitTest::CheckParallelBools_tpetra(status,failPID); \
    if(verbosity>=level && status)                 \
       os << msgp << std::endl;                    \
    else if(verbosity>=level && not status)        \
       os << msgf << " ( PID = " << failPID << " )" << std::endl; \
    }

#define TEST_EQUALITY(x,y,msg) \
   status = (x==y); \
   if(not status || verbosity>=10) { \
      os << msg << std::endl; \
   } \
   allPassed &= status;

#define TEST_NOT_EQUAL(x,y,msg) \
   status = (x!=y); \
   if(not status || verbosity>=10) { \
      os << msg << std::endl; \
   } \
   allPassed &= status;

#define TEST_MSG(msg) \
   if(verbosity>=10) { \
      os << msg << std::endl; \
   }

#define TEST_ASSERT(x,msg) \
   status = x; \
   if(not status || verbosity>=10) { \
      os << msg << std::endl; \
   } \
   allPassed &= status;

#endif
