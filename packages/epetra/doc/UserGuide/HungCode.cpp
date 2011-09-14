/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
*/

  // Scenario: We are solving a system of equations Ax = b.
  //           We know A and b, and have an approximation to x
  //           These two versions of code attempt compute the 
  //           2-norm of the residual r where r = b - Ax.



  // Code fragment 1:  THIS CODE WILL NOT WORK
  // Fails because the statement A.Multiply() typicall involves 
  // interprocessor communication but processor 0 is the only one executing
  // the Multiply() method.  The user program will probably stall.
  if (Comm().MyPID()==0) { // Only processor 0 will execute this code
      A.Multiply(false, q, r); // Compute Ax, store in r
      r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
      r.Norm2(&rnorm2);
	 cout << "2-norm of b - Ax = " << r << endl;
  } 



  // Code fragment 2:  THIS CODE WILL NOT WORK
  // The Multiply() method will work.  The Update() method will complete, but results 
  // will be incorrect.  Finally this code will stall on the call to Norm2() 
  // because it involved a collective operation where all processor must 
  // participate.
  A.Multiply(false, q, r); // Compute Ax, store in r
  if (Comm().MyPID()==0) { // Only processor 0 will execute this code
      r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
      r.Norm2(&rnorm2);
	 cout << "2-norm of b - Ax = " << r << endl;
  } 



  // Code fragment 3:  THIS CODE WILL NOT WORK
  // The Multiply() and Update() methods will work.  
  // As with the previous segment this code will stall on the call to Norm2() 
  // because it involved a collective operation where all processor must 
  // participate.
  A.Multiply(false, q, r); // Compute Ax, store in r
  r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
  if (Comm().MyPID()==0) { // Only processor 0 will execute this code
      r.Norm2(&rnorm2);
	 cout << "2-norm of b - Ax = " << r << endl;
  } 



  // Code fragment 4:  THIS CODE WILL WORK
  // All methods will work.  Only the output statement itself is restricted to
  // processor 0.
  A.Multiply(false, q, r); // Compute Ax, store in r
  r.Update(1.0, b, -1.0);  // r = b - r = b - Ax
      r.Norm2(&rnorm2);
  if (Comm().MyPID()==0)  // Only processor 0 will execute this code
	 cout << "2-norm of b - Ax = " << r << endl;
