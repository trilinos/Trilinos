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
#include "MueLu_MatlabUtils_def.hpp"
#include <stdexcept>

#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA) || !defined(HAVE_MUELU_TPETRA)
#error "Muemex types require MATLAB, Epetra and Tpetra."
#else

/* Stuff for MATLAB R2006b vs. previous versions */
#if(defined(MX_API_VER) && MX_API_VER >= 0x07030000)
#else
typedef int mwIndex;
#endif

using namespace std;
using namespace Teuchos;

namespace MueLu {

/* Explicit instantiation of Muemexdata */
  template class MuemexData<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >; 
  template class MuemexData<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
  template class MuemexData<RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >; 
  template class MuemexData<RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
  template class MuemexData<RCP<MAggregates>>;
  template class MuemexData<RCP<MAmalInfo>>;
  /*  
  template class MuemexData<int>;

template<> class MuemexData<complex_t>;	  
  template<> class MuemexData<std::string>; 
  template<> class MuemexData<double>;					
  template<> class MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >; 
  template<> class MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >; 
  template<> class MuemexData<RCP<Epetra_MultiVector> >;			
  template<> class MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >; 
  template<> class MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >; 

  */

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

/* ******************************* */
/* Specializations                 */
/* ******************************* */

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

/******************************/
/* callback functions         */
/******************************/

void callMatlabNoArgs(std::string function)
{
  int result = mexEvalString(function.c_str());
  if(result != 0)
    mexPrintf("An error occurred while running a MATLAB command.");\
}

std::vector<RCP<MuemexArg>> callMatlab(std::string function, int numOutputs, std::vector<RCP<MuemexArg>> args)
{
  using Teuchos::rcp_static_cast;
  mxArray** matlabArgs = new mxArray* [args.size()];
  mxArray** matlabOutput = new mxArray* [numOutputs];
  std::vector<RCP<MuemexArg>> output;
  for(int i = 0; i < int(args.size()); i++)
    {
      try
        {
          switch(args[i]->type)
            {
            case INT:
              matlabArgs[i] = rcp_static_cast<MuemexData<int>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<double>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case STRING:
              matlabArgs[i] = rcp_static_cast<MuemexData<std::string>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<complex_t>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_ORDINAL_VECTOR:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_ordinal_vector>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MULTIVECTOR_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MULTIVECTOR_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MATRIX_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MATRIX_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MATRIX_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_Matrix_double>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MATRIX_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_Matrix_complex>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MULTIVECTOR_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_MultiVector_double>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MULTIVECTOR_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_MultiVector_complex>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case EPETRA_CRSMATRIX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Epetra_CrsMatrix>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case EPETRA_MULTIVECTOR:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Epetra_MultiVector>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case AGGREGATES:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<MAggregates>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case AMALGAMATION_INFO:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<MAmalInfo>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            }
        }
      catch (std::exception& e)
        {
          mexPrintf("An error occurred while converting arg #%d to MATLAB:\n", i);
	  std::cout << e.what() << std::endl;
          mexPrintf("Passing 0 instead.\n");
          matlabArgs[i] = mxCreateDoubleScalar(0);
        }
    }
  //now matlabArgs is populated with MATLAB data types
  int result = mexCallMATLAB(numOutputs, matlabOutput, args.size(), matlabArgs, function.c_str());
  if(result != 0)
    mexPrintf("Matlab encountered an error while running command through muemexCallbacks.\n");
  //now, if all went well, matlabOutput contains all the output to return to user
  for(int i = 0; i < numOutputs; i++)
    {
      try
        {
          //Identify the type of each output, and put into output vector
          mxArray* item = matlabOutput[i];
          switch(mxGetClassID(item))
            {
            case mxCHAR_CLASS:
              //string
              output.push_back(rcp(new MuemexData<std::string>(item)));
              break;
            case mxINT32_CLASS:
              if(mxGetM(item) == 1 && mxGetN(item) == 1)
                //individual integer
                output.push_back(rcp(new MuemexData<int>(item)));
              else if(mxGetM(item) != 1 || mxGetN(item) != 1)
                //ordinal vector
                output.push_back(rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(item)));
              else
                throw std::runtime_error("Error: Don't know what to do with integer array.\n");
              break;
            case mxDOUBLE_CLASS:
              if(mxGetM(item) == 1 && mxGetN(item) == 1)
                {
                  if(mxIsComplex(item))
                    //single double (scalar, real)
                    output.push_back(rcp(new MuemexData<double>(item)));
                  else
                    //single complex scalar
                    output.push_back(rcp(new MuemexData<complex_t>(item)));
              }
              else if(mxIsSparse(item))
              {
                  //Default to Tpetra matrix for this
                  if(mxIsComplex(item))
                    //complex Tpetra matrix (sparse)
                    output.push_back(rcp(new MuemexData<RCP<Tpetra_CrsMatrix_double>>(item)));
                  else
                    //real Tpetra matrix
                    output.push_back(rcp(new MuemexData<RCP<Tpetra_CrsMatrix_complex>>(item)));
              }
              else
              {
                  //Default to Tpetra multivector for this case
                  if(mxIsComplex(item))
                    output.push_back(rcp(new MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(item)));
                  else
                    output.push_back(rcp(new MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(item)));
              }
              break;
            case mxSTRUCT_CLASS:
            {
              //the only thing that should get here currently is an Aggregates struct
              //verify that is has the correct fields with the correct types
              //also assume that aggregates data will not be stored in an array of more than 1 element.
              if(isValidMatlabAggregates(item))
                output.push_back(rcp(new MuemexData<RCP<MAggregates>>(item)));
              else
                throw runtime_error("Invalid aggregates struct passed in from MATLAB.");
              break;
            }
            default:
	        throw std::runtime_error("MATLAB returned an unsupported type as a function output.\n");
            }
        }
	catch(std::exception& e)
        {
          mexPrintf("An error occurred while converting output #%d from MATLAB:\n", i);
	  std::cout << e.what() << std::endl;
        }
    }
  return output;
}

/******************************/
/* More utility functions     */
/******************************/

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
              //global row, # of entries, value array, column indices array
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
  catch(std::exception& e)
    {
      mexPrintf("An error occurred while setting up an Epetra matrix:\n");
      std::cout << e.what() << std::endl;
    }
  return matrix;
}

mxArray* createMatlabLOVector(RCP<Xpetra_ordinal_vector>& vec)
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

mxArray* saveEpetraMV(RCP<Epetra_MultiVector>& mv)
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
      throw std::runtime_error("Error: Unrecognized numerical type.");
    }
  return rv;
}

mxArray* saveEpetraMatrix(RCP<Epetra_CrsMatrix>& mat)
{
  RCP<Xpetra_Matrix_double> xmat = MueLu::EpetraCrs_To_XpetraMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(mat);
  return saveMatrixToMatlab<double>(xmat);
}

RCP<Xpetra_ordinal_vector> loadLOVector(const mxArray* mxa)
{
  RCP<const Teuchos::Comm<int> > comm = rcp(new Teuchos::SerialComm<int>());
  const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
  RCP<const muemex_map_type> rowMap = rcp(new muemex_map_type(numGlobalIndices, (mm_GlobalOrd) 0, comm));
  if(mxGetClassID(mxa) != mxINT32_CLASS || mxGetN(mxa) != 1)
    throw std::runtime_error("Can only construct LOVector with int32 single vector.");
  int* array = (int*) mxGetData(mxa);
  Teuchos::ArrayView<int> dataView(array, mxGetM(mxa));
  RCP<Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loVec = rcp(new Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(rowMap, dataView));
  return Xpetra::toXpetra<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(loVec);
}

template<typename Scalar = double>
RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> xpetraLoadMatrix(const mxArray* mxa)
{
  return MueLu::TpetraCrs_To_XpetraMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tpetraLoadMatrix<Scalar>(mxa));
}

//Aggregates structs in MATLAB either have to be constructed with constructAggregates() or at least have the same fields
RCP<MAggregates> loadAggregates(const mxArray* mxa)
{
  if(mxGetNumberOfElements(mxa) != 1)
    throw runtime_error("Aggregates must be individual structs in MATLAB.");
  if(!mxIsStruct(mxa))
    throw runtime_error("Trying to pull aggregates from non-struct MATLAB object.");
  //assume that in matlab aggregate structs will only be stored in a 1x1 array
  //mxa must have the same fields as the ones declared in constructAggregates function in muelu.m for this to work
  const int correctNumFields = 5; //change if more fields are added to the aggregates representation in constructAggregates in muelu.m
  if(mxGetNumberOfFields(mxa) != correctNumFields)
    throw runtime_error("Aggregates structure has wrong number of fields.");
  //Pull MuemexData types back out
  int nVert;
  int nAgg;
  RCP<Xpetra_ordinal_vector> vertToAggIdLOV;
  RCP<Xpetra_ordinal_vector> aggSizesLOV;
  RCP<Xpetra_ordinal_vector> rootNodesIn; //roots are stored in MATLAB struct as list of ints (representing nodeids)
  { //scope this so variable names can be simpler, and extra MuemexData wrappers are disposed of immediately
    RCP<MuemexData<int>> mmData = rcp(new MuemexData<int>(mxGetField(mxa, 0, "nVertices")));
    nVert = mmData->getData();
  }
  {
    RCP<MuemexData<int>> mmData = rcp(new MuemexData<int>(mxGetField(mxa, 0, "nAggregates")));
    nAgg = mmData->getData();
  }
  {
    RCP<MuemexData<RCP<Xpetra_ordinal_vector>>> mmData = rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(mxGetField(mxa, 0, "vertexToAggID")));
    vertToAggIdLOV = mmData->getData();
  }
  {
    RCP<MuemexData<RCP<Xpetra_ordinal_vector>>> mmData = rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(mxGetField(mxa, 0, "aggSizes")));
    aggSizesLOV = mmData->getData();
  }
  {
    RCP<MuemexData<RCP<Xpetra_ordinal_vector>>> mmData = rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(mxGetField(mxa, 0, "rootNodes")));
    rootNodesIn = mmData->getData();
  }
  //Now have all the data needed to fully reconstruct the aggregate
  //Use similar approach as UserAggregationFactory (which is written for >1 thread but will just be serial here)
  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
  int myRank = comm->getRank();
  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
  RCP<Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t>> map = Xpetra::MapFactory<mm_LocalOrd, mm_GlobalOrd, mm_node_t>::Build(lib, nVert, 0, comm);
  RCP<MAggregates> agg = rcp(new MAggregates(map));
  agg->SetNumAggregates(nAgg);
  //Get handles for the vertex2AggId and procwinner arrays in reconstituted aggregates object
  //this is serial so all procwinner values will be same (0)
  ArrayRCP<mm_LocalOrd> vertex2AggId = agg->GetVertex2AggId()->getDataNonConst(0);  //the '0' means first (and only) column of multivector, since is just vector
  ArrayRCP<mm_LocalOrd> procWinner = agg->GetProcWinner()->getDataNonConst(0);
  //Also get const ArrayRCP views into vert2AggId and aggSizes vectors (of LOVectors extracted from matlab struct)
  ArrayRCP<const mm_LocalOrd> vertex2AggIDToCopy = vertToAggIdLOV->getData(0);
  ArrayRCP<const mm_LocalOrd> aggSizesToCopy = aggSizesLOV->getData(0);
  ArrayRCP<const mm_LocalOrd> rootNodesToCopy = rootNodesIn->getData(0);
  //mm_LocalOrd and int are equivalent, so is ok to talk about aggSize with just 'int'
  //Deep copy the entire vertex2AggID and isRoot arrays, which are both nVert items long
  //At the same time, set ProcWinner
  for(int i = 0; i < nVert; i++)
  {
    vertex2AggId[i] = vertex2AggIDToCopy[i];
    procWinner[i] = myRank; //all nodes are going to be on the same proc
    agg->SetIsRoot(i, false); //the ones that are root will be set in next loop
  }
  for(int i = 0; i < nAgg; i++) //rootNodesToCopy is an array of node IDs which are the roots of their aggs
  {
    agg->SetIsRoot(rootNodesToCopy[i], true);
  }
  //Now recompute the aggSize array and cache the results in the object
  agg->ComputeAggregateSizes(true, true);
  agg->AggregatesCrossProcessors(false);
  return agg;
}

mxArray* saveAggregates(RCP<MAggregates>& agg)
{
  //Grab data from agg
  mm_LocalOrd nAgg = agg->GetNumAggregates();
  RCP<Xpetra_ordinal_vector> vertexToAggID = (RCP<Xpetra_ordinal_vector>) agg->GetVertex2AggId();
  mm_LocalOrd nVert = vertexToAggID->getLocalLength();
  //Construct Matlab objects to use in struct
  vector<RCP<MuemexArg>> structMembers; //input to the matlab callback
  vector<RCP<MuemexArg>> builtStruct;   //output
  mwSize arrLen[] = {1}; //array of one element to hold dimensions of 1-D array for mxCreateNumArr
  mxArray* nVertices = mxCreateNumericArray(1, arrLen, mxINT32_CLASS, mxREAL);
  mxArray* nAggregates = mxCreateNumericArray(1, arrLen, mxINT32_CLASS, mxREAL);
  *((int*) mxGetData(nVertices)) = (int) nVert;
  *((int*) mxGetData(nAggregates)) = (int) nAgg;
  //agg sizes
  ArrayRCP<mm_LocalOrd> rawAggSizes = agg->ComputeAggregateSizes(false, true);
  mm_LocalOrd* rawAggSizeArray = rawAggSizes.getRawPtr();
  ArrayView<mm_LocalOrd> aggSizeArrayView(rawAggSizeArray, nAgg);
  RCP<const Teuchos::Comm<int>> comm = rcp(new Teuchos::SerialComm<int>());
  const Tpetra::global_size_t numGlobalIndices = nAgg;
  RCP<const muemex_map_type> rowMap = rcp(new muemex_map_type(numGlobalIndices, (mm_GlobalOrd) 0, comm));
  RCP<Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> aggSizeTempLOVec = rcp(new Tpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(rowMap, aggSizeArrayView));
  RCP<Xpetra_ordinal_vector> aggSizes = Xpetra::toXpetra<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(aggSizeTempLOVec);
  //Place matlab objects in MuemexData objects
  RCP<MuemexData<int>> nVertArg = rcp(new MuemexData<int>(nVert));
  RCP<MuemexData<int>> nAggArg = rcp(new MuemexData<int>(nAgg));
  RCP<Xpetra::TpetraMap<mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetraVecMap = rcp(new Xpetra::TpetraMap<mm_LocalOrd, mm_GlobalOrd, mm_node_t>(nAgg, (mm_GlobalOrd) 0, comm));
  RCP<Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t>> vecMap = rcp_implicit_cast<Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t>, Xpetra::TpetraMap<mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(tpetraVecMap);
  //now have the Xpetra::Map that we need to construct rootnode ordinal vector
  mm_LocalOrd* rawRootNodes = (mm_LocalOrd*) malloc(sizeof(mm_LocalOrd) * nAgg);
  int iter = 0;
  for(int i = 0; i < nVert; i++)
  {
    if(iter > nAgg)
      throw runtime_error("Number of root nodes and number of aggregates don't match!");
    if(agg->IsRoot(i))
    {
      rawRootNodes[iter] = i;
      iter++;
    }
  }
  ArrayView<mm_LocalOrd> rootNodeArrayView(rawRootNodes, nAgg);
  RCP<Xpetra::TpetraVector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetraRootNodeVec = rcp(new Xpetra::TpetraVector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(vecMap, rootNodeArrayView));
  if(rawRootNodes == NULL)
    throw runtime_error("Malloc failed!");
  RCP<Xpetra_ordinal_vector> rootNodes = rcp_implicit_cast<Xpetra_ordinal_vector, Xpetra::TpetraVector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(tpetraRootNodeVec);
  RCP<MuemexData<RCP<Xpetra_ordinal_vector>>> vertexToAggIDArg = rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(vertexToAggID));
  RCP<MuemexData<RCP<Xpetra_ordinal_vector>>> rootNodesArg = rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(rootNodes));
  RCP<MuemexData<RCP<Xpetra_ordinal_vector>>> aggSizesArg = rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(aggSizes));
  //Finally, upcast the inputs to MuemexArg and put in arg vector (order matters here!)
  structMembers.push_back(rcp_implicit_cast<MuemexArg, MuemexData<int>>(nVertArg));
  structMembers.push_back(rcp_implicit_cast<MuemexArg, MuemexData<int>>(nAggArg));
  structMembers.push_back(rcp_implicit_cast<MuemexArg, MuemexData<RCP<Xpetra_ordinal_vector>>>(vertexToAggIDArg));
  structMembers.push_back(rcp_implicit_cast<MuemexArg, MuemexData<RCP<Xpetra_ordinal_vector>>>(rootNodesArg));
  structMembers.push_back(rcp_implicit_cast<MuemexArg, MuemexData<RCP<Xpetra_ordinal_vector>>>(aggSizesArg));
  builtStruct = callMatlab("constructAggregates", 1, structMembers);
  //builtStruct[0] should be a MuemexArg wrapping the MAggregates object. Get that, convert to matlab struct and return
  //avoid buffer overrun, throw if builtStruct[0] doesn't exist
  if(builtStruct.size() == 0)
    throw runtime_error("Matlab aggregates struct builder didn't return the struct.");
  free(rawRootNodes); //have to free this, the ArrayView doesn't automatically
  RCP<MuemexData<RCP<MAggregates>>> outData = rcp_static_cast<MuemexData<RCP<MAggregates>>, MuemexArg>(builtStruct[0]);
  return outData->convertToMatlab();
}

RCP<MAmalInfo> loadAmalInfo(const mxArray* mxa)
{
  RCP<MAmalInfo> amal;
  throw runtime_error("AmalgamationInfo not supported in Muemex yet.");
  return amal;
}

mxArray* saveAmalInfo(RCP<MAmalInfo>& amalInfo)
{
  throw runtime_error("AmalgamationInfo not supported in MueMex yet.");
  return mxCreateDoubleScalar(0);
}

bool isValidMatlabAggregates(const mxArray* mxa)
{
  bool isValidAggregates = true;
  if(!mxIsStruct(mxa))
    return false;
  int numFields = mxGetNumberOfFields(mxa); //check that struct has correct # of fields
  if(numFields != 5)
    isValidAggregates = false;
  if(isValidAggregates)
  {
    const char* mem1 = mxGetFieldNameByNumber(mxa, 0);
    if(mem1 == NULL || strcmp(mem1, "nVertices") != 0)
      isValidAggregates = false;
    const char* mem2 = mxGetFieldNameByNumber(mxa, 1);
    if(mem2 == NULL || strcmp(mem2, "nAggregates") != 0)
      isValidAggregates = false;
    const char* mem3 = mxGetFieldNameByNumber(mxa, 2);
    if(mem3 == NULL || strcmp(mem3, "vertexToAggID") != 0)
      isValidAggregates = false;
    const char* mem4 = mxGetFieldNameByNumber(mxa, 3);
    if(mem3 == NULL || strcmp(mem4, "rootNodes") != 0)
      isValidAggregates = false;
    const char* mem5 = mxGetFieldNameByNumber(mxa, 4);
    if(mem4 == NULL || strcmp(mem5, "aggSizes") != 0)
      isValidAggregates = false;
  }
  return isValidAggregates;
}

}//end namespace
#endif // HAVE_MUELU_MATLAB
