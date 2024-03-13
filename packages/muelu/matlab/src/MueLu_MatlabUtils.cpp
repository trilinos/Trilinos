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

#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA)
#error "Muemex types require MATLAB, Epetra and Tpetra."
#else

/* Stuff for MATLAB R2006b vs. previous versions */
#if (defined(MX_API_VER) && MX_API_VER >= 0x07030000)
#else
typedef int mwIndex;
#endif

using namespace std;
using namespace Teuchos;

namespace MueLu {

/* Explicit instantiation of MuemexData variants */
template class MuemexData<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<MAggregates> >;
template class MuemexData<RCP<MAmalInfo> >;
template class MuemexData<int>;
template class MuemexData<bool>;
template class MuemexData<complex_t>;
template class MuemexData<string>;
template class MuemexData<double>;
template class MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<Epetra_MultiVector> >;
template class MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;
template class MuemexData<RCP<Xpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >;

// Flag set to true if MATLAB's CSC matrix index type is not int (usually false)
bool rewrap_ints = sizeof(int) != sizeof(mwIndex);

int* mwIndex_to_int(int N, mwIndex* mwi_array) {
  // int* rv = (int*) malloc(N * sizeof(int));
  int* rv = new int[N];  // not really better but may avoid confusion for valgrind
  for (int i = 0; i < N; i++)
    rv[i] = (int)mwi_array[i];
  return rv;
}

/* ******************************* */
/* Specializations                 */
/* ******************************* */

template <>
mxArray* createMatlabSparse<double>(int numRows, int numCols, int nnz) {
  return mxCreateSparse(numRows, numCols, nnz, mxREAL);
}

template <>
mxArray* createMatlabSparse<complex_t>(int numRows, int numCols, int nnz) {
  return mxCreateSparse(numRows, numCols, nnz, mxCOMPLEX);
}

template <>
void fillMatlabArray<double>(double* array, const mxArray* mxa, int n) {
  memcpy(mxGetPr(mxa), array, n * sizeof(double));
}

template <>
void fillMatlabArray<complex_t>(complex_t* array, const mxArray* mxa, int n) {
  double* pr = mxGetPr(mxa);
  double* pi = mxGetPi(mxa);
  for (int i = 0; i < n; i++) {
    pr[i] = std::real<double>(array[i]);
    pi[i] = std::imag<double>(array[i]);
  }
}

/******************************/
/* Callback Functions         */
/******************************/

void callMatlabNoArgs(std::string function) {
  int result = mexEvalString(function.c_str());
  if (result != 0)
    mexPrintf("An error occurred while running a MATLAB command.");
}

std::vector<RCP<MuemexArg> > callMatlab(std::string function, int numOutputs, std::vector<RCP<MuemexArg> > args) {
  using Teuchos::rcp_static_cast;
  mxArray** matlabArgs   = new mxArray*[args.size()];
  mxArray** matlabOutput = new mxArray*[numOutputs];
  std::vector<RCP<MuemexArg> > output;

  for (int i = 0; i < int(args.size()); i++) {
    try {
      switch (args[i]->type) {
        case BOOL:
          matlabArgs[i] = rcp_static_cast<MuemexData<bool>, MuemexArg>(args[i])->convertToMatlab();
          break;
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
        case XPETRA_MAP:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_map> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case XPETRA_ORDINAL_VECTOR:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_ordinal_vector> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case TPETRA_MULTIVECTOR_DOUBLE:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case TPETRA_MULTIVECTOR_COMPLEX:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case TPETRA_MATRIX_DOUBLE:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case TPETRA_MATRIX_COMPLEX:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case XPETRA_MATRIX_DOUBLE:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_Matrix_double> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case XPETRA_MATRIX_COMPLEX:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_Matrix_complex> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case XPETRA_MULTIVECTOR_DOUBLE:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_MultiVector_double> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case XPETRA_MULTIVECTOR_COMPLEX:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_MultiVector_complex> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case EPETRA_CRSMATRIX:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Epetra_CrsMatrix> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case EPETRA_MULTIVECTOR:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Epetra_MultiVector> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case AGGREGATES:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<MAggregates> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case AMALGAMATION_INFO:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<MAmalInfo> >, MuemexArg>(args[i])->convertToMatlab();
          break;
        case GRAPH:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<MGraph> >, MuemexArg>(args[i])->convertToMatlab();
#ifdef HAVE_MUELU_INTREPID2
        case FIELDCONTAINER_ORDINAL:
          matlabArgs[i] = rcp_static_cast<MuemexData<RCP<FieldContainer_ordinal> >, MuemexArg>(args[i])->convertToMatlab();
          break;
#endif
      }
    } catch (std::exception& e) {
      mexPrintf("An error occurred while converting arg #%d to MATLAB:\n", i);
      std::cout << e.what() << std::endl;
      mexPrintf("Passing 0 instead.\n");
      matlabArgs[i] = mxCreateDoubleScalar(0);
    }
  }
  // now matlabArgs is populated with MATLAB data types
  int result = mexCallMATLAB(numOutputs, matlabOutput, args.size(), matlabArgs, function.c_str());
  if (result != 0)
    mexPrintf("Matlab encountered an error while running command through muemexCallbacks.\n");
  // now, if all went well, matlabOutput contains all the output to return to user
  for (int i = 0; i < numOutputs; i++) {
    try {
      output.push_back(convertMatlabVar(matlabOutput[i]));
    } catch (std::exception& e) {
      mexPrintf("An error occurred while converting output #%d from MATLAB:\n", i);
      std::cout << e.what() << std::endl;
    }
  }
  delete[] matlabOutput;
  delete[] matlabArgs;
  return output;
}

/******************************/
/* More utility functions     */
/******************************/

template <>
mxArray* createMatlabMultiVector<double>(int numRows, int numCols) {
  return mxCreateDoubleMatrix(numRows, numCols, mxREAL);
}

template <>
mxArray* createMatlabMultiVector<complex_t>(int numRows, int numCols) {
  return mxCreateDoubleMatrix(numRows, numCols, mxCOMPLEX);
}

mxArray* saveAmalInfo(RCP<MAmalInfo>& amalInfo) {
  throw runtime_error("AmalgamationInfo not supported in MueMex yet.");
  return mxCreateDoubleScalar(0);
}

bool isValidMatlabAggregates(const mxArray* mxa) {
  bool isValidAggregates = true;
  if (!mxIsStruct(mxa))
    return false;
  int numFields = mxGetNumberOfFields(mxa);  // check that struct has correct # of fields
  if (numFields != 5)
    isValidAggregates = false;
  if (isValidAggregates) {
    const char* mem1 = mxGetFieldNameByNumber(mxa, 0);
    if (mem1 == NULL || strcmp(mem1, "nVertices") != 0)
      isValidAggregates = false;
    const char* mem2 = mxGetFieldNameByNumber(mxa, 1);
    if (mem2 == NULL || strcmp(mem2, "nAggregates") != 0)
      isValidAggregates = false;
    const char* mem3 = mxGetFieldNameByNumber(mxa, 2);
    if (mem3 == NULL || strcmp(mem3, "vertexToAggID") != 0)
      isValidAggregates = false;
    const char* mem4 = mxGetFieldNameByNumber(mxa, 3);
    if (mem3 == NULL || strcmp(mem4, "rootNodes") != 0)
      isValidAggregates = false;
    const char* mem5 = mxGetFieldNameByNumber(mxa, 4);
    if (mem4 == NULL || strcmp(mem5, "aggSizes") != 0)
      isValidAggregates = false;
  }
  return isValidAggregates;
}

bool isValidMatlabGraph(const mxArray* mxa) {
  bool isValidGraph = true;
  if (!mxIsStruct(mxa))
    return false;
  int numFields = mxGetNumberOfFields(mxa);  // check that struct has correct # of fields
  if (numFields != 2)
    isValidGraph = false;
  if (isValidGraph) {
    const char* mem1 = mxGetFieldNameByNumber(mxa, 0);
    if (mem1 == NULL || strcmp(mem1, "edges") != 0)
      isValidGraph = false;
    const char* mem2 = mxGetFieldNameByNumber(mxa, 1);
    if (mem2 == NULL || strcmp(mem2, "boundaryNodes") != 0)
      isValidGraph = false;
  }
  return isValidGraph;
}

std::vector<std::string> tokenizeList(const std::string& params) {
  using namespace std;
  vector<string> rlist;
  const char* delims = ",";
  char* copy         = (char*)malloc(params.length() + 1);
  strcpy(copy, params.c_str());
  char* mark = (char*)strtok(copy, delims);
  while (mark != NULL) {
    // Remove leading and trailing whitespace in token
    char* tail = mark + strlen(mark) - 1;
    while (*mark == ' ')
      mark++;
    while (*tail == ' ' && tail > mark)
      tail--;
    tail++;
    *tail = 0;
    string tok(mark);  // copies the characters to string object
    rlist.push_back(tok);
    mark = strtok(NULL, delims);
  }
  free(copy);
  return rlist;
}

Teuchos::RCP<Teuchos::ParameterList> getInputParamList() {
  using namespace Teuchos;
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory for the matrix A.");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Factory for the prolongator.");
  validParamList->set<RCP<const FactoryBase> >("R", Teuchos::null, "Factory for the restrictor.");
  validParamList->set<RCP<const FactoryBase> >("Ptent", Teuchos::null, "Factory for the tentative (unsmoothed) prolongator.");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for the node coordinates.");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Factory for the nullspace.");
  validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Factory for the aggregates.");
  validParamList->set<RCP<const FactoryBase> >("UnamalgamationInfo", Teuchos::null, "Factory for amalgamation.");
#ifdef HAVE_MUELU_INTREPID2
  validParamList->set<RCP<const FactoryBase> >("pcoarsen: element to node map", Teuchos::null, "Generating factory of the element to node map");
#endif
  return validParamList;
}

Teuchos::RCP<MuemexArg> convertMatlabVar(const mxArray* mxa) {
  switch (mxGetClassID(mxa)) {
    case mxCHAR_CLASS:
      // string
      return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<std::string>(mxa)));
      break;
    case mxLOGICAL_CLASS:
      // boolean
      return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<bool>(mxa)));
      break;
    case mxINT32_CLASS:
      if (mxGetM(mxa) == 1 && mxGetN(mxa) == 1)
        // individual integer
        return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<int>(mxa)));
      else if (mxGetM(mxa) != 1 || mxGetN(mxa) != 1)
        // ordinal vector
        return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<Xpetra_ordinal_vector> >(mxa)));
      else
        throw std::runtime_error("Error: Don't know what to do with integer array.\n");
      break;
    case mxDOUBLE_CLASS:
      if (mxGetM(mxa) == 1 && mxGetN(mxa) == 1) {
        if (mxIsComplex(mxa))
          // single double (scalar, real)
          return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<complex_t>(mxa)));
        else
          // single complex scalar
          return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<double>(mxa)));
      } else if (mxIsSparse(mxa))  // use a CRS matrix
      {
        // Default to Tpetra matrix for this
        if (mxIsComplex(mxa))
          // complex matrix
          return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<Xpetra_Matrix_complex> >(mxa)));
        else
          // real-valued matrix
          return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<Xpetra_Matrix_double> >(mxa)));
      } else {
        // Default to Xpetra multivector for this case
        if (mxIsComplex(mxa))
          return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(mxa)));
        else
          return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(mxa)));
      }
      break;
    case mxSTRUCT_CLASS: {
      // the only thing that should get here currently is an Aggregates struct or Graph struct
      // verify that it has the correct fields with the correct types
      // also assume that aggregates data will not be stored in an array of more than 1 element.
      if (isValidMatlabAggregates(mxa)) {
        return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<MAggregates> >(mxa)));
      } else if (isValidMatlabGraph(mxa)) {
        return rcp_implicit_cast<MuemexArg>(rcp(new MuemexData<RCP<MGraph> >(mxa)));
      } else {
        throw runtime_error("Invalid aggregates or graph struct passed in from MATLAB.");
        return Teuchos::null;
      }
      break;
    }
    default:
      throw std::runtime_error("MATLAB returned an unsupported type as a function output.\n");
      return Teuchos::null;
  }
}

/******************************/
/* Explicit Instantiations    */
/******************************/

template bool loadDataFromMatlab<bool>(const mxArray* mxa);
template int loadDataFromMatlab<int>(const mxArray* mxa);
template double loadDataFromMatlab<double>(const mxArray* mxa);
template complex_t loadDataFromMatlab<complex_t>(const mxArray* mxa);
template string loadDataFromMatlab<string>(const mxArray* mxa);
template RCP<Xpetra_ordinal_vector> loadDataFromMatlab<RCP<Xpetra_ordinal_vector> >(const mxArray* mxa);
template RCP<Tpetra_MultiVector_double> loadDataFromMatlab<RCP<Tpetra_MultiVector_double> >(const mxArray* mxa);
template RCP<Tpetra_MultiVector_complex> loadDataFromMatlab<RCP<Tpetra_MultiVector_complex> >(const mxArray* mxa);
template RCP<Tpetra_CrsMatrix_double> loadDataFromMatlab<RCP<Tpetra_CrsMatrix_double> >(const mxArray* mxa);
template RCP<Tpetra_CrsMatrix_complex> loadDataFromMatlab<RCP<Tpetra_CrsMatrix_complex> >(const mxArray* mxa);
template RCP<Xpetra_Matrix_double> loadDataFromMatlab<RCP<Xpetra_Matrix_double> >(const mxArray* mxa);
template RCP<Xpetra_Matrix_complex> loadDataFromMatlab<RCP<Xpetra_Matrix_complex> >(const mxArray* mxa);
template RCP<Xpetra_MultiVector_double> loadDataFromMatlab<RCP<Xpetra_MultiVector_double> >(const mxArray* mxa);
template RCP<Xpetra_MultiVector_complex> loadDataFromMatlab<RCP<Xpetra_MultiVector_complex> >(const mxArray* mxa);
template RCP<Epetra_CrsMatrix> loadDataFromMatlab<RCP<Epetra_CrsMatrix> >(const mxArray* mxa);
template RCP<Epetra_MultiVector> loadDataFromMatlab<RCP<Epetra_MultiVector> >(const mxArray* mxa);
template RCP<MAggregates> loadDataFromMatlab<RCP<MAggregates> >(const mxArray* mxa);
template RCP<MAmalInfo> loadDataFromMatlab<RCP<MAmalInfo> >(const mxArray* mxa);

template mxArray* saveDataToMatlab(bool& data);
template mxArray* saveDataToMatlab(int& data);
template mxArray* saveDataToMatlab(double& data);
template mxArray* saveDataToMatlab(complex_t& data);
template mxArray* saveDataToMatlab(string& data);
template mxArray* saveDataToMatlab(RCP<Xpetra_ordinal_vector>& data);
template mxArray* saveDataToMatlab(RCP<Tpetra_MultiVector_double>& data);
template mxArray* saveDataToMatlab(RCP<Tpetra_MultiVector_complex>& data);
template mxArray* saveDataToMatlab(RCP<Tpetra_CrsMatrix_double>& data);
template mxArray* saveDataToMatlab(RCP<Tpetra_CrsMatrix_complex>& data);
template mxArray* saveDataToMatlab(RCP<Xpetra_Matrix_double>& data);
template mxArray* saveDataToMatlab(RCP<Xpetra_Matrix_complex>& data);
template mxArray* saveDataToMatlab(RCP<Xpetra_MultiVector_double>& data);
template mxArray* saveDataToMatlab(RCP<Xpetra_MultiVector_complex>& data);
template mxArray* saveDataToMatlab(RCP<Epetra_CrsMatrix>& data);
template mxArray* saveDataToMatlab(RCP<Epetra_MultiVector>& data);
template mxArray* saveDataToMatlab(RCP<MAggregates>& data);
template mxArray* saveDataToMatlab(RCP<MAmalInfo>& data);

template vector<RCP<MuemexArg> > processNeeds<double>(const Factory* factory, string& needsParam, Level& lvl);
template vector<RCP<MuemexArg> > processNeeds<complex_t>(const Factory* factory, string& needsParam, Level& lvl);
template void processProvides<double>(vector<RCP<MuemexArg> >& mexOutput, const Factory* factory, string& providesParam, Level& lvl);
template void processProvides<complex_t>(vector<RCP<MuemexArg> >& mexOutput, const Factory* factory, string& providesParam, Level& lvl);

}  // namespace MueLu
#endif  // HAVE_MUELU_MATLAB
