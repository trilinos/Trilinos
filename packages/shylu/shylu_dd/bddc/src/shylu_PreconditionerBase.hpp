
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef BDDC_PRECONDITIONERBASE_H
#define BDDC_PRECONDITIONERBASE_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include <set>

namespace bddc {
  
template <class SX, class SM, class LO, class GO> 
  class PreconditionerBase
{
public:

  PreconditionerBase()
  {
  }

  virtual ~PreconditionerBase()
  {
  }

  virtual int InitialStaticCondensation(SX* rightHandSide1to1,
					SX* initialSolution1to1)
  {
    return 0;
  }

  virtual void StaticExpansion(SX* deltaSolution,
			       const bool interfaceValuesOnlyInput,
			       const bool useCurrentInteriorValues,
			       SX* initialSolution)
  {
  }

  virtual LO NumMyRows()=0;

  virtual LO NumMyRowsKrylov()=0;

  virtual void Apply(const SX* r, 
		     SX* Pr,
		     SX* APr,
		     const bool interfaceValuesOnly)=0;

  virtual void Apply(const SX* r, 
		     SX* Pr,
		     const bool interfaceValuesOnly)=0;

  virtual void ApplyFullOperator(SX* x, 
				 SX* Ax)=0;

  virtual void ReduceSum(SX* inputVec,
			 LO numTerms,
			 SX* summedVec)=0;

  virtual SM Norm2(SX* x, 
		   LO numTerm)=0;

  virtual SX DotProd(SX* x, 
		     SX* y, 
		     LO numTerm)=0;

  virtual int MyPID()
  {
    return 0;
  }

  virtual int numProc()
  {
    return 1;
  }

  virtual double GetTime()
  {
    return MPI_Wtime();
  }

private:

};

} // namespace bddc

#endif // BDDC_PRECONDITIONERBASE_H
  
