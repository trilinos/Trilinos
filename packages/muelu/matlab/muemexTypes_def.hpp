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
#include "muemexTypes_decl.hpp"

#ifndef MUEMEX_TYPES_DEF_HPP
#define MUEMEX_TYPES_DEF_HPP



#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA) || !defined(HAVE_MUELU_TPETRA)
#error "Muemex types require MATLAB, Epetra and Tpetra."
#else
#include "mex.h"
#include <Tpetra_DefaultPlatform.hpp>
#include <stdexcept>

using Teuchos::RCP;
using Teuchos::rcp;

extern bool rewrap_ints;

/* ******************************* */
/* Specializations                 */
/* ******************************* */

template<>
RCP<Tpetra_CrsMatrix_double> tpetraLoadMatrix<double>(const mxArray* mxa)
{
  bool success = false;
  RCP<Tpetra_CrsMatrix_double> A;
  try
    {
      RCP<const Teuchos::Comm<int>> comm = rcp(new Teuchos::SerialComm<int>());
      //numGlobalIndices is just the number of rows in the matrix
      const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
      const mm_GlobalOrd indexBase = 0;
      RCP<const muemex_map_type> rowMap = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
      RCP<const muemex_map_type> domainMap = rcp(new muemex_map_type(mxGetN(mxa), indexBase, comm));
      A = Tpetra::createCrsMatrix<double, mm_GlobalOrd, mm_LocalOrd, mm_node_t>(rowMap);
      double* valueArray = mxGetPr(mxa);
      int* colptr;
      int* rowind;
      //int nr = mxGetM(mxa);
      int nc = mxGetN(mxa);
      if(rewrap_ints)
        {
          //mwIndex_to_int allocates memory so must delete[] later
          colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
          rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
        }
      else
        {
          rowind = (int*) mxGetIr(mxa);
          colptr = (int*) mxGetJc(mxa);
        }
      for(int i = 0; i < nc; i++)
        {
          for(int j = colptr[i]; j < colptr[i + 1]; j++)
            {
              //'array' of 1 element, containing column (in global matrix).
              Teuchos::ArrayView<mm_GlobalOrd> cols = Teuchos::ArrayView<mm_GlobalOrd>(&i, 1);
              //'array' of 1 element, containing value
              Teuchos::ArrayView<double> vals = Teuchos::ArrayView<double>(&valueArray[j], 1);
              A->insertGlobalValues(rowind[j], cols, vals);
            }
        }
      A->fillComplete(domainMap, rowMap);
      if(rewrap_ints)
        {
          delete[] rowind;
          delete[] colptr;
        }
      success = true;
    }
  catch(std::exception& e)
    {
      mexPrintf("Error while constructing Tpetra matrix:\n");
      std::cout << e.what() << std::endl;
    }
  if(!success)
    mexErrMsgTxt("An error occurred while setting up a Tpetra matrix.\n");
  return A;
}

template<>
RCP<Tpetra_CrsMatrix_complex> tpetraLoadMatrix<complex_t>(const mxArray* mxa)
{
  RCP<Tpetra_CrsMatrix_complex> A;
  //Create a map in order to create the matrix (taken from muelu basic example - complex)
  try
    {
      RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
      const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
      const mm_GlobalOrd indexBase = 0;
      RCP<const muemex_map_type> map = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
      A = rcp(new Tpetra_CrsMatrix_complex(map, 0));
      double* realArray = mxGetPr(mxa);
      double* imagArray = mxGetPi(mxa);
      int* colptr;
      int* rowind;
      int nc = mxGetN(mxa);
      if(rewrap_ints)
        {
          //mwIndex_to_int allocates memory so must delete[] later
          colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
          rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
        }
      else
        {
          rowind = (int*) mxGetIr(mxa);
          colptr = (int*) mxGetJc(mxa);
        }
      for(int i = 0; i < nc; i++)
        {
          for(int j = colptr[i]; j < colptr[i + 1]; j++)
            {
              //here assuming that complex_t will always be defined as std::complex<double>
              //use 'value' over and over again with Teuchos::ArrayViews to insert into matrix
              complex_t value = std::complex<double>(realArray[j], imagArray[j]);
              Teuchos::ArrayView<mm_GlobalOrd> cols = Teuchos::ArrayView<mm_GlobalOrd>(&i, 1);
              Teuchos::ArrayView<complex_t> vals = Teuchos::ArrayView<complex_t>(&value, 1);
              A->insertGlobalValues(rowind[j], cols, vals);
            }
        }
      A->fillComplete();
      if(rewrap_ints)
        {
          delete[] rowind;
          delete[] colptr;
        }
    }
  catch(std::exception& e)
    {
      mexPrintf("Error while constructing tpetra matrix:\n");
      std::cout << e.what() << std::endl;
    }
  return A;
}


template<>
RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadTpetraMV<double>(const mxArray* mxa)
{
  RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv;
  try
    {
      int nr = mxGetM(mxa);
      int nc = mxGetN(mxa);
      double* pr = mxGetPr(mxa);
      RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
      //numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
      RCP<const muemex_map_type> map = rcp(new muemex_map_type(nr, (mm_GlobalOrd) 0, comm));
      //Allocate a new array of complex values to use with the multivector
      Teuchos::ArrayView<const double> arrView(pr, nr * nc);
      mv = rcp(new Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(map, arrView, size_t(nr), size_t(nc)));
    }
  catch(std::exception& e)
    {
      mexPrintf("Error constructing Tpetra MultiVector.\n");
      std::cout << e.what() << std::endl;
    }
  return mv;
}

template<>
RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadTpetraMV<complex_t>(const mxArray* mxa)
{
  RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv;
  try
    {
      int nr = mxGetM(mxa);
      int nc = mxGetN(mxa);
      double* pr = mxGetPr(mxa);
      double* pi = mxGetPi(mxa);
      RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
      //numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
      RCP<const muemex_map_type> map = rcp(new muemex_map_type(nr, (mm_GlobalOrd) 0, comm));
      //Allocate a new array of complex values to use with the multivector
      complex_t* myArr = new complex_t[nr * nc];
      for(int n = 0; n < nc; n++)
        {
          for(int m = 0; m < nr; m++)
            {
              myArr[n * nr + m] = complex_t(pr[n * nr + m], pi[n * nr + m]);
            }
        }
      Teuchos::ArrayView<complex_t> arrView(myArr, nr * nc);
      mv = rcp(new Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(map, arrView, nr, nc));
    }
  catch(std::exception& e)
    {
      mexPrintf("Error constructing Tpetra MultiVector.\n");
      std::cout << e.what() << std::endl;
    }
  return mv;
}


/* ******************************* */
/* Begin MuemexData implementation */
/* ******************************* */
//Fully generic methods
template<typename T>
MuemexData<T>::MuemexData(T& dataToCopy, MUEMEX_TYPE dataType) : MuemexArg(dataType)
{
  data = dataToCopy;
}

template<typename T>
MuemexData<T>::~MuemexData() {}

template<typename T>
T& MuemexData<T>::getData()
{
  return data;
}

template<typename T>
void MuemexData<T>::setData(T& newData)
{
  this->data = data;
}

//string specializations
template<>
MuemexData<std::string>::MuemexData(const mxArray* mxa) : MuemexArg(STRING)
{
  data = "";
  if(!mxGetClassID(mxa) != mxCHAR_CLASS)
    {
      throw std::runtime_error("Can't construct string from anything but a char array.");
    }
  data = std::string(mxArrayToString(mxa));
}

template<>
mxArray* MuemexData<std::string>::convertToMatlab()
{
  return mxCreateString(data.c_str());
}

//int
template<>
MuemexData<int>::MuemexData(const mxArray* mxa) : MuemexArg(INT)
{
  data = parseInt(mxa);
}

template<>
mxArray* MuemexData<int>::convertToMatlab()
{
  mxArray* output = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
  int* ptr = (int*) mxGetData(output);
  *ptr = data;
  return output;
}

//double
template<>
MuemexData<double>::MuemexData(const mxArray* mxa) : MuemexArg(DOUBLE)
{
  data = *((double*) mxGetPr(mxa));
}

template<>
mxArray* MuemexData<double>::convertToMatlab()
{
  return mxCreateDoubleScalar(data);
}

//complex scalar
template<>
MuemexData<complex_t>::MuemexData(const mxArray* mxa) : MuemexArg(COMPLEX)
{
  double* realPart = mxGetPr(mxa);
  double* imagPart = mxGetPi(mxa);
  data = complex_t(*realPart, *imagPart);
}

template<>
mxArray* MuemexData<complex_t>::convertToMatlab()
{
  mxArray* output = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
  double* realPart = mxGetPr(output);
  double* imagPart = mxGetPi(output);
  *realPart = std::real<double>(data);
  *imagPart = std::imag<double>(data);
  return output;
}

//Epetra_Crs
template<>
MuemexData<RCP<Epetra_CrsMatrix>>::MuemexData(const mxArray* mxa) : MuemexArg(EPETRA_CRSMATRIX)
{
  data = epetraLoadMatrix(mxa);
}

template<>
mxArray* MuemexData<RCP<Epetra_CrsMatrix>>::convertToMatlab()
{
  return saveEpetraMatrix(data);
}

//Tpetra_Crs double
template<>
MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MATRIX_DOUBLE)
{
  data = tpetraLoadMatrix<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveMatrixToMatlab<double>(MueLu::TpetraCrs_To_XpetraMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(data));
}

//Tpetra_Crs complex
template<>
MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MATRIX_COMPLEX)
{
  data = tpetraLoadMatrix<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveMatrixToMatlab<complex_t>(MueLu::TpetraCrs_To_XpetraMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(data));
}

//Xpetra matrix double
template<>
MuemexData<RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MATRIX_DOUBLE)
{
  data = xpetraLoadMatrix<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveMatrixToMatlab<double>(data);
}

//Xpetra matrix complex
template<>
MuemexData<RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MATRIX_COMPLEX)
{
  data = xpetraLoadMatrix<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveMatrixToMatlab<complex_t>(data);
}

//Epetra MV
template<>
MuemexData<RCP<Epetra_MultiVector>>::MuemexData(const mxArray* mxa) : MuemexArg(EPETRA_MULTIVECTOR)
{
  data = loadEpetraMV(mxa);
}

template<>
mxArray* MuemexData<RCP<Epetra_MultiVector>>::convertToMatlab()
{
  return saveEpetraMV(data);
}

//Tpetra MV double
template<>
MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MULTIVECTOR_DOUBLE)
{
  data = loadTpetraMV<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveTpetraMV<double>(data);
}

//Tpetra MV complex
template<>
MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(TPETRA_MULTIVECTOR_COMPLEX)
{
  data = loadTpetraMV<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveTpetraMV<complex_t>(data);
}

//Xpetra ordinal vector
template<>
MuemexData<RCP<Xpetra_ordinal_vector>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_ORDINAL_VECTOR)
{
  data = loadLOVector(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra_ordinal_vector>>::convertToMatlab()
{
  return createMatlabLOVector(data);
}

//Xpetra multivector double
template<>
MuemexData<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MULTIVECTOR_DOUBLE)
{
  data = loadXpetraMV<double>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveMultiVectorToMatlab<double>(data);
}

//Xpetra multivector complex
template<>
MuemexData<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::MuemexData(const mxArray* mxa) : MuemexArg(XPETRA_MULTIVECTOR_COMPLEX)
{
  data = loadXpetraMV<complex_t>(mxa);
}

template<>
mxArray* MuemexData<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>::convertToMatlab()
{
  return saveMultiVectorToMatlab<complex_t>(data);
}

/* ***************************** */
/* End MuemexData implementation */
/* ***************************** */

template<typename Scalar = double>
mxArray* saveMatrixToMatlab(RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mat)
{
  int nr = mat->getGlobalNumRows();
  int nc = mat->getGlobalNumCols();
  int nnz = mat->getGlobalNumEntries();
#ifdef VERBOSE_OUTPUT
  RCP<Teuchos::FancyOStream> fancyStream = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  mat->describe(*fancyStream, Teuchos::VERB_EXTREME);
#endif
  mxArray* mxa = createMatlabSparse<Scalar>(nr, nc, nnz);
  mwIndex* ir = mxGetIr(mxa);
  mwIndex* jc = mxGetJc(mxa);
  for(int i = 0; i < nc + 1; i++)
    {
      jc[i] = 0;
    }
  size_t maxEntriesPerRow = mat->getGlobalMaxNumRowEntries();
  int* rowProgress = new int[nc];
  //The array that will be copied to Pr and (if complex) Pi later
  Scalar* sparseVals = new Scalar[nnz];
  size_t numEntries;
  if(mat->isLocallyIndexed())
    {
      Scalar* rowValArray = new Scalar[maxEntriesPerRow];
      Teuchos::ArrayView<Scalar> rowVals(rowValArray, maxEntriesPerRow);
      mm_LocalOrd* rowIndicesArray = new mm_LocalOrd[maxEntriesPerRow];
      Teuchos::ArrayView<mm_LocalOrd> rowIndices(rowIndicesArray, maxEntriesPerRow);
      for(mm_LocalOrd m = 0; m < nr; m++)       //All rows in the Xpetra matrix
        {
          mat->getLocalRowCopy(m, rowIndices, rowVals, numEntries);     //Get the row
          for(mm_LocalOrd entry = 0; entry < int(numEntries); entry++)  //All entries in row
            {
              jc[rowIndices[entry] + 1]++; //for each entry, increase jc for the entry's column
            }
        }
      //now jc holds the number of elements in each column, but needs cumulative sum over all previous columns also
      int entriesAccum = 0;
      for(int n = 0; n <= nc; n++)
        {
          int temp = entriesAccum;
          entriesAccum += jc[n];
          jc[n] += temp;
        }
      //Jc now populated with colptrs
      for(int i = 0; i < nc; i++)
        {
          rowProgress[i] = 0;
        }
      //Row progress values like jc but keep track as the MATLAB matrix is being filled in
      for(mm_LocalOrd m = 0; m < nr; m++)       //rows
        {
          mat->getLocalRowCopy(m, rowIndices, rowVals, numEntries);
          for(mm_LocalOrd i = 0; i < int(numEntries); i++)      //entries in row m (NOT columns)
            {
              //row is m, col is rowIndices[i], val is rowVals[i]
              mm_LocalOrd col = rowIndices[i];
              sparseVals[jc[col] + rowProgress[col]] = rowVals[i];      //Set value
              ir[jc[col] + rowProgress[col]] = m;                                               //Set row at which value occurs
              rowProgress[col]++;
            }
        }
      delete[] rowIndicesArray;
    }
  else
    {
      Teuchos::ArrayView<const mm_GlobalOrd> rowIndices;
      Teuchos::ArrayView<const Scalar> rowVals;
      for(mm_GlobalOrd m = 0; m < nr; m++)
        {
          mat->getGlobalRowView(m, rowIndices, rowVals);
          for(mm_GlobalOrd n = 0; n < rowIndices.size(); n++)
            {
              jc[rowIndices[n] + 1]++;
            }
        }
      //Last element of jc is just nnz
      jc[nc] = nnz;
      //Jc now populated with colptrs
      for(int i = 0; i < nc; i++)
        {
          rowProgress[i] = 0;
        }
      int entriesAccum = 0;
      for(int n = 0; n <= nc; n++)
        {
          int temp = entriesAccum;
          entriesAccum += jc[n];
          jc[n] += temp;
        }
      //Row progress values like jc but keep track as the MATLAB matrix is being filled in
      for(mm_GlobalOrd m = 0; m < nr; m++)      //rows
        {
          mat->getGlobalRowView(m, rowIndices, rowVals);
          for(mm_LocalOrd i = 0; i < rowIndices.size(); i++)    //entries in row m (NOT == columns)
            {
              //row is m, col is rowIndices[i], val is rowVals[i]
              mm_GlobalOrd col = rowIndices[i];
              sparseVals[jc[col] + rowProgress[col]] = rowVals[i];      //Set value
              ir[jc[col] + rowProgress[col]] = m;                                               //Set row at which value occurs
              rowProgress[col]++;
            }
        }
    }
  //finally, copy sparseVals into pr (and pi, if complex)
  fillMatlabArray<Scalar>(sparseVals, mxa, nnz);
  delete[] sparseVals;
  delete[] rowProgress;
  return mxa;
}

template<> mxArray* createMatlabSparse<double>(int numRows, int numCols, int nnz)
{
  return mxCreateSparse(numRows, numCols, nnz, mxREAL);
}

template<> mxArray* createMatlabSparse<complex_t>(int numRows, int numCols, int nnz)
{
  return mxCreateSparse(numRows, numCols, nnz, mxCOMPLEX);
}

template<> void fillMatlabArray<double>(double* array, const mxArray* mxa, int n)
{
  memcpy(mxGetPr(mxa), array, n * sizeof(double));
}

template<> void fillMatlabArray<complex_t>(complex_t* array, const mxArray* mxa, int n)
{
  double* pr = mxGetPr(mxa);
  double* pi = mxGetPi(mxa);
  for(int i = 0; i < n; i++)
    {
      pr[i] = std::real<double>(array[i]);
      pi[i] = std::imag<double>(array[i]);
    }
}

template<typename Scalar>
mxArray* saveMultiVectorToMatlab(RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv)
{
  //Precondition: Memory has already been allocated by MATLAB for the array.
  int nr = mv->getGlobalLength();
  int nc = mv->getNumVectors();
  mxArray* output = createMatlabMultiVector<Scalar>(nr, nc);
  Scalar* data = new Scalar[nr * nc];
  for(int col = 0; col < nc; col++)
    {
      Teuchos::ArrayRCP<const Scalar> colData = mv->getData(col);
      for(int row = 0; row < nr; row++)
        {
          data[col * nr + row] = colData[row];
        }
    }
  fillMatlabArray<Scalar>(data, output, nc * nr);
  return output;
}

template<typename Scalar>
mxArray* saveTpetraMV(RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv)
{
  //Precondition: Memory has already been allocated by MATLAB for the array.
  int nr = mv->getGlobalLength();
  int nc = mv->getNumVectors();
  mxArray* output = createMatlabMultiVector<Scalar>(nr, nc);
  Scalar* data = new Scalar[nr * nc];
  for(int col = 0; col < nc; col++)
    {
      Teuchos::ArrayRCP<const Scalar> colData = mv->getData(col);
      for(int row = 0; row < nr; row++)
        {
          data[col * nr + row] = colData[row];
        }
    }
  fillMatlabArray<Scalar>(data, output, nc * nr);
  return output;
}

template<typename Scalar>
RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadXpetraMV(const mxArray* mxa)
{
  RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > tmv = loadTpetraMV<Scalar>(mxa);
  return MueLu::TpetraMultiVector_To_XpetraMultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tmv);
}

#endif //HAVE_MUELU_MATLAB error handler
#endif //MUEMEX_TYPES_DEF_HPP guard
