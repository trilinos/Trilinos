// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TabularOutputter.hpp"
#include "Domi_MDArray.hpp"

namespace
{

using Domi::MDArray;
using Domi::Ordinal;
using Teuchos::tuple;

int numLoops = 100;
int intPrec = 8;
int dblPrec = 6;
int dim1 = 10;
int dim2 = 8;
int dim3 = 6;
int dim4 = 4;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor & clp = Teuchos::UnitTestRepository::getCLP();
}

TEUCHOS_UNIT_TEST( MDArray, parenOperator1D )
{
  typedef Teuchos::TabularOutputter TO;

  out << std::endl << "TEUCHOS_ARRAY_BOUNDSCHECK is ";
  if (MDArray<double>::hasBoundsChecking()) out << "ON";
  else out << "OFF";
  out << std::endl << std::endl;

  TO outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT,    intPrec);

  outputter.pushFieldSpec("dim"         , TO::INT   );
  outputter.pushFieldSpec("num loops"   , TO::INT   );
  outputter.pushFieldSpec("paren access", TO::DOUBLE);

  outputter.outputHeader();

  int scale[4] = {1, 10, 100, 1000};
  for (int test_case_k = 0; test_case_k < 4; ++test_case_k)
  {
    Ordinal arrayDim = dim1 * scale[test_case_k];

    // dim
    outputter.outputField(arrayDim);

    // num loops
    outputter.outputField(numLoops);

    MDArray< double > mda(tuple(arrayDim));

    // paren access
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numLoops, arrayDim)
    {
      for (Ordinal ii=0; ii < arrayDim; ++ii)
        mda(ii) = 0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, time);

    outputter.nextRow();
  }
}

TEUCHOS_UNIT_TEST( MDArray, parenOperator2D )
{
  typedef Teuchos::TabularOutputter TO;

  TO outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT,    intPrec);

  outputter.pushFieldSpec("dim 1"       , TO::INT   );
  outputter.pushFieldSpec("dim 2"       , TO::INT   );
  outputter.pushFieldSpec("num loops"   , TO::INT   );
  outputter.pushFieldSpec("paren access", TO::DOUBLE);

  outputter.outputHeader();

  int scale[4] = {1, 2, 5, 10};
  for (int test_case_k = 0; test_case_k < 4; ++test_case_k)
  {
    Ordinal arrayDim1 = dim1 * scale[test_case_k];
    Ordinal arrayDim2 = dim2 * scale[test_case_k];

    // dim 1
    outputter.outputField(arrayDim1);

    // dim 2
    outputter.outputField(arrayDim2);

    // num loops
    outputter.outputField(numLoops);

    MDArray< double > mda(tuple(arrayDim1, arrayDim2));

    // paren access
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numLoops,
                                              arrayDim1*arrayDim2)
    {
      for (Ordinal jj=0; jj < arrayDim2; ++jj)
        for (Ordinal ii=0; ii < arrayDim1; ++ii)
          mda(ii,jj) = 0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, time);

    outputter.nextRow();
  }
}

TEUCHOS_UNIT_TEST( MDArray, parenOperator3D )
{
  typedef Teuchos::TabularOutputter TO;

  TO outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT,    intPrec);

  outputter.pushFieldSpec("dim 1"       , TO::INT   );
  outputter.pushFieldSpec("dim 2"       , TO::INT   );
  outputter.pushFieldSpec("dim 3"       , TO::INT   );
  outputter.pushFieldSpec("num loops"   , TO::INT   );
  outputter.pushFieldSpec("paren access", TO::DOUBLE);

  outputter.outputHeader();

  int scale[4] = {1, 2, 3, 4};
  for (int test_case_k = 0; test_case_k < 4; ++test_case_k)
  {
    Ordinal arrayDim1 = dim1 * scale[test_case_k];
    Ordinal arrayDim2 = dim2 * scale[test_case_k];
    Ordinal arrayDim3 = dim3 * scale[test_case_k];

    // dim 1
    outputter.outputField(arrayDim1);

    // dim 2
    outputter.outputField(arrayDim2);

    // dim 3
    outputter.outputField(arrayDim3);

    // num loops
    outputter.outputField(numLoops);

    MDArray< double > mda(tuple(arrayDim1, arrayDim2, arrayDim3));

    // paren access
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numLoops,
                                              arrayDim1*arrayDim2*arrayDim3)
    {
      for (Ordinal kk=0; kk < arrayDim3; ++kk)
        for (Ordinal jj=0; jj < arrayDim2; ++jj)
          for (Ordinal ii=0; ii < arrayDim1; ++ii)
            mda(ii,jj,kk) = 0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, time);

    outputter.nextRow();
  }
}

TEUCHOS_UNIT_TEST( MDArray, parenOperator4D )
{
  typedef Teuchos::TabularOutputter TO;

  TO outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT,    intPrec);

  outputter.pushFieldSpec("dim 1"       , TO::INT   );
  outputter.pushFieldSpec("dim 2"       , TO::INT   );
  outputter.pushFieldSpec("dim 3"       , TO::INT   );
  outputter.pushFieldSpec("dim 4"       , TO::INT   );
  outputter.pushFieldSpec("num loops"   , TO::INT   );
  outputter.pushFieldSpec("paren access", TO::DOUBLE);

  outputter.outputHeader();

  int scale[4] = {1, 2, 3};
  for (int test_case_k = 0; test_case_k < 3; ++test_case_k)
  {
    Ordinal arrayDim1 = dim1 * scale[test_case_k];
    Ordinal arrayDim2 = dim2 * scale[test_case_k];
    Ordinal arrayDim3 = dim3 * scale[test_case_k];
    Ordinal arrayDim4 = dim4 * scale[test_case_k];

    // dim 1
    outputter.outputField(arrayDim1);

    // dim 2
    outputter.outputField(arrayDim2);

    // dim 3
    outputter.outputField(arrayDim3);

    // dim 4
    outputter.outputField(arrayDim4);

    // num loops
    outputter.outputField(numLoops);

    MDArray< double > mda(tuple(arrayDim1, arrayDim2, arrayDim3, arrayDim4));

    // paren access
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numLoops,
                                              arrayDim1*arrayDim2*arrayDim3*arrayDim4)
    {
      for (Ordinal mm=0; mm < arrayDim4; ++mm)
        for (Ordinal kk=0; kk < arrayDim3; ++kk)
          for (Ordinal jj=0; jj < arrayDim2; ++jj)
            for (Ordinal ii=0; ii < arrayDim1; ++ii)
              mda(ii,jj,kk,mm) = 0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, time);

    outputter.nextRow();
  }
}

}
