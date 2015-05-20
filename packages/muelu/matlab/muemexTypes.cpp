// @HEADER
//
// ***********************************************************************
//
//                MueLu: A package for multigrid based preconditioning
//                                      Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                                        Jonathan Hu           (jhu@sandia.gov)
//                                        Andrey Prokopenko (aprokop@sandia.gov)
//                                        Ray Tuminaro          (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Tpetra_DefaultPlatform.hpp>
#include "MueLu_Utilities.hpp"
#include "muemexTypes_decl.hpp"
#include "mex.h"

//Debug this file
//#define VERBOSE_OUTPUT

using namespace std;
using namespace Teuchos;

/* Stuff for MATLAB R2006b vs. previous versions */
#if(defined(MX_API_VER) && MX_API_VER >= 0x07030000)
#else
typedef int mwIndex;
#endif

  
//Flag set to true if MATLAB's CSC matrix index type is not int (usually false)
bool rewrap_ints = sizeof(int) != sizeof(mwIndex);

int* mwIndex_to_int(int N, mwIndex* mwi_array)
{
  int* rv;
  rv = new int[N];
  for(int i = 0; i < N; i++)
    rv[i] = (int) mwi_array[i];
  return rv;
}




RCP<Epetra_CrsMatrix> epetraLoadMatrix(const mxArray* mxa)
{
  RCP<Epetra_CrsMatrix> matrix;
  try
    {
      int* colptr;
      int* rowind;
      double* vals = mxGetPr(mxa);
      int nr = mxGetM(mxa);
      int nc = mxGetN(mxa);
      if(rewrap_ints)
        {
          colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
          rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
        }
      else
        {
          rowind = (int*) mxGetIr(mxa);
          colptr = (int*) mxGetJc(mxa);
        }
      Epetra_SerialComm Comm;
      Epetra_Map RangeMap(nr, 0, Comm);
      Epetra_Map DomainMap(nc, 0, Comm);
      matrix = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RangeMap, DomainMap, 0));
      /* Do the matrix assembly */
      for(int i = 0; i < nc; i++)
        {
          for(int j = colptr[i]; j < colptr[i + 1]; j++)
            {
              //                global row, # of entries, value array, column indices array
              matrix->InsertGlobalValues(rowind[j], 1, &vals[j], &i);
            }
        }
      matrix->FillComplete(DomainMap, RangeMap);
      if(rewrap_ints)
        {
          delete [] rowind;
          delete [] colptr;
        }
    }
  catch(exception& e)
    {
      mexPrintf("An error occurred while setting up an Epetra matrix:\n");
      cout << e.what() << endl;
    }
  return matrix;
}



mxArray* createMatlabLOVector(RCP<Xpetra_ordinal_vector> vec)
{
  //this value might be a 64 bit int but it should never overflow a 32
  mwSize len = vec->getGlobalLength();
  //create a single column vector
  mwSize dimensions[] = {len, 1};
  return mxCreateNumericArray(2, dimensions, mxINT32_CLASS, mxREAL);
}


RCP<Epetra_MultiVector> loadEpetraMV(const mxArray* mxa)
{
  int nr = mxGetM(mxa);
  int nc = mxGetN(mxa);
  Epetra_SerialComm Comm;
  Epetra_BlockMap map(nr * nc, 1, 0, Comm);
  return rcp(new Epetra_MultiVector(Epetra_DataAccess::Copy, map, mxGetPr(mxa), nr, nc));
}

template<> mxArray* createMatlabMultiVector<double>(int numRows, int numCols)
{
  return mxCreateDoubleMatrix(numRows, numCols, mxREAL);
}

template<> mxArray* createMatlabMultiVector<complex_t>(int numRows, int numCols)
{
  return mxCreateDoubleMatrix(numRows, numCols, mxCOMPLEX);
}


mxArray* saveEpetraMV(RCP<Epetra_MultiVector> mv)
{
  mxArray* output = mxCreateDoubleMatrix(mv->GlobalLength(), mv->NumVectors(), mxREAL);
  double* dataPtr = mxGetPr(output);
  mv->ExtractCopy(dataPtr, mv->GlobalLength());
  return output;
}



int parseInt(const mxArray* mxa)
{
  mxClassID probIDtype = mxGetClassID(mxa);
  int rv;
  if(probIDtype == mxINT32_CLASS)
    {
      rv = *((int*) mxGetData(mxa));
    }
  else if(probIDtype == mxDOUBLE_CLASS)
    {
      rv = (int) *((double*) mxGetData(mxa));
    }
  else if(probIDtype == mxUINT32_CLASS)
    {
      rv = (int) *((unsigned int*) mxGetData(mxa));
    }
  else
    {
      rv = -1;
      throw runtime_error("Error: Unrecognized numerical type.");
    }
  return rv;
}

mxArray* saveEpetraMatrix(RCP<Epetra_CrsMatrix> mat)
{
  return saveMatrixToMatlab<double>(MueLu::EpetraCrs_To_XpetraMatrix<double, mm_LocalOrd, mm_GlobalOrd,mm_node_t>(mat));
}

RCP<Xpetra_ordinal_vector> loadLOVector(const mxArray* mxa)
{
  RCP<const Teuchos::Comm<int> > comm = rcp(new Teuchos::SerialComm<int>());
  const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
  RCP<const muemex_map_type> rowMap = rcp(new muemex_map_type(numGlobalIndices, (mm_GlobalOrd) 0, comm));
  if(mxGetClassID(mxa) != mxINT32_CLASS || mxGetN(mxa) != 1)
    throw runtime_error("Can only construct LOVector with int32 single vector.");
  int* array = (int*) mxGetData(mxa);
  ArrayView<int> dataView(array, mxGetM(mxa));
  RCP<Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loVec = rcp(new Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(rowMap, dataView));
  return Xpetra::toXpetra<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(loVec);
}

