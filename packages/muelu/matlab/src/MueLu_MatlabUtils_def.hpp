// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MATLABUTILS_DEF_HPP
#define MUELU_MATLABUTILS_DEF_HPP

#include "MueLu_MatlabUtils_decl.hpp"

#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA)
#error "Muemex types require MATLAB, Epetra and Tpetra."
#else

using Teuchos::RCP;
using Teuchos::rcp;
using namespace std;

namespace MueLu {

extern bool rewrap_ints;

/* ******************************* */
/* getMuemexType                   */
/* ******************************* */

template <typename T>
MuemexType getMuemexType(const T& data) { throw std::runtime_error("Unknown Type"); }

template <>
MuemexType getMuemexType(const int& data) { return INT; }
template <>
MuemexType getMuemexType<int>() { return INT; }
template <>
MuemexType getMuemexType<bool>() { return BOOL; }

template <>
MuemexType getMuemexType(const double& data) { return DOUBLE; }
template <>
MuemexType getMuemexType<double>() { return DOUBLE; }

template <>
MuemexType getMuemexType(const std::string& data) { return STRING; }
template <>
MuemexType getMuemexType<string>() { return STRING; }

template <>
MuemexType getMuemexType(const complex_t& data) { return COMPLEX; }
template <>
MuemexType getMuemexType<complex_t>() { return COMPLEX; }

template <>
MuemexType getMuemexType(const RCP<Xpetra_map>& data) { return XPETRA_MAP; }
template <>
MuemexType getMuemexType<RCP<Xpetra_map> >() { return XPETRA_MAP; }

template <>
MuemexType getMuemexType(const RCP<Xpetra_ordinal_vector>& data) { return XPETRA_ORDINAL_VECTOR; }
template <>
MuemexType getMuemexType<RCP<Xpetra_ordinal_vector> >() { return XPETRA_ORDINAL_VECTOR; }

template <>
MuemexType getMuemexType(const RCP<Tpetra_MultiVector_double>& data) { return TPETRA_MULTIVECTOR_DOUBLE; }
template <>
MuemexType getMuemexType<RCP<Tpetra_MultiVector_double> >() { return TPETRA_MULTIVECTOR_DOUBLE; }

template <>
MuemexType getMuemexType(const RCP<Tpetra_MultiVector_complex>& data) { return TPETRA_MULTIVECTOR_COMPLEX; }
template <>
MuemexType getMuemexType<RCP<Tpetra_MultiVector_complex> >() { return TPETRA_MULTIVECTOR_COMPLEX; }

template <>
MuemexType getMuemexType(const RCP<Tpetra_CrsMatrix_double>& data) { return TPETRA_MATRIX_DOUBLE; }
template <>
MuemexType getMuemexType<RCP<Tpetra_CrsMatrix_double> >() { return TPETRA_MATRIX_DOUBLE; }

template <>
MuemexType getMuemexType(const RCP<Tpetra_CrsMatrix_complex>& data) { return TPETRA_MATRIX_COMPLEX; }
template <>
MuemexType getMuemexType<RCP<Tpetra_CrsMatrix_complex> >() { return TPETRA_MATRIX_COMPLEX; }

template <>
MuemexType getMuemexType(const RCP<Xpetra_MultiVector_double>& data) { return XPETRA_MULTIVECTOR_DOUBLE; }
template <>
MuemexType getMuemexType<RCP<Xpetra_MultiVector_double> >() { return XPETRA_MULTIVECTOR_DOUBLE; }

template <>
MuemexType getMuemexType(const RCP<Xpetra_MultiVector_complex>& data) { return XPETRA_MULTIVECTOR_COMPLEX; }
template <>
MuemexType getMuemexType<RCP<Xpetra_MultiVector_complex> >() { return XPETRA_MULTIVECTOR_COMPLEX; }

template <>
MuemexType getMuemexType(const RCP<Xpetra_Matrix_double>& data) { return XPETRA_MATRIX_DOUBLE; }
template <>
MuemexType getMuemexType<RCP<Xpetra_Matrix_double> >() { return XPETRA_MATRIX_DOUBLE; }

template <>
MuemexType getMuemexType(const RCP<Xpetra_Matrix_complex>& data) { return XPETRA_MATRIX_COMPLEX; }
template <>
MuemexType getMuemexType<RCP<Xpetra_Matrix_complex> >() { return XPETRA_MATRIX_COMPLEX; }

template <>
MuemexType getMuemexType(const RCP<Epetra_CrsMatrix>& data) { return EPETRA_CRSMATRIX; }
template <>
MuemexType getMuemexType<RCP<Epetra_CrsMatrix> >() { return EPETRA_CRSMATRIX; }

template <>
MuemexType getMuemexType(const RCP<Epetra_MultiVector>& data) { return EPETRA_MULTIVECTOR; }
template <>
MuemexType getMuemexType<RCP<Epetra_MultiVector> >() { return EPETRA_MULTIVECTOR; }

template <>
MuemexType getMuemexType(const RCP<MAggregates>& data) { return AGGREGATES; }
template <>
MuemexType getMuemexType<RCP<MAggregates> >() { return AGGREGATES; }

template <>
MuemexType getMuemexType(const RCP<MAmalInfo>& data) { return AMALGAMATION_INFO; }
template <>
MuemexType getMuemexType<RCP<MAmalInfo> >() { return AMALGAMATION_INFO; }

template <>
MuemexType getMuemexType(const RCP<MGraph>& data) { return GRAPH; }
template <>
MuemexType getMuemexType<RCP<MGraph> >() { return GRAPH; }

#ifdef HAVE_MUELU_INTREPID2
template <>
MuemexType getMuemexType(const RCP<FieldContainer_ordinal>& data) { return FIELDCONTAINER_ORDINAL; }
template <>
MuemexType getMuemexType<RCP<FieldContainer_ordinal> >() { return FIELDCONTAINER_ORDINAL; }
#endif

/* "prototypes" for specialized functions used in other specialized functions */

template <>
mxArray* createMatlabSparse<double>(int numRows, int numCols, int nnz);
template <>
mxArray* createMatlabSparse<complex_t>(int numRows, int numCols, int nnz);
template <>
mxArray* createMatlabMultiVector<double>(int numRows, int numCols);
template <>
mxArray* createMatlabMultiVector<complex_t>(int numRows, int numCols);
template <>
void fillMatlabArray<double>(double* array, const mxArray* mxa, int n);
template <>
void fillMatlabArray<complex_t>(complex_t* array, const mxArray* mxa, int n);
template <>
mxArray* saveDataToMatlab(RCP<Xpetra_MultiVector_double>& data);
template <>
mxArray* saveDataToMatlab(RCP<Xpetra_MultiVector_complex>& data);
template <>
mxArray* saveDataToMatlab(RCP<Xpetra_Matrix_double>& data);
template <>
mxArray* saveDataToMatlab(RCP<Xpetra_Matrix_complex>& data);

/* ******************************* */
/* loadDataFromMatlab              */
/* ******************************* */

template <>
int loadDataFromMatlab<int>(const mxArray* mxa) {
  mxClassID probIDtype = mxGetClassID(mxa);
  int rv;
  if (probIDtype == mxINT32_CLASS) {
    rv = *((int*)mxGetData(mxa));
  } else if (probIDtype == mxLOGICAL_CLASS) {
    rv = (int)*((bool*)mxGetData(mxa));
  } else if (probIDtype == mxDOUBLE_CLASS) {
    rv = (int)*((double*)mxGetData(mxa));
  } else if (probIDtype == mxUINT32_CLASS) {
    rv = (int)*((unsigned int*)mxGetData(mxa));
  } else {
    rv = -1;
    throw std::runtime_error("Error: Unrecognized numerical type.");
  }
  return rv;
}

template <>
bool loadDataFromMatlab<bool>(const mxArray* mxa) {
  return *((bool*)mxGetData(mxa));
}

template <>
double loadDataFromMatlab<double>(const mxArray* mxa) {
  return *((double*)mxGetPr(mxa));
}

template <>
complex_t loadDataFromMatlab<complex_t>(const mxArray* mxa) {
  double realpart = real<double>(*((double*)mxGetPr(mxa)));
  double imagpart = imag<double>(*((double*)mxGetPi(mxa)));
  return complex_t(realpart, imagpart);
}

template <>
string loadDataFromMatlab<string>(const mxArray* mxa) {
  string rv = "";
  if (mxGetClassID(mxa) != mxCHAR_CLASS) {
    throw runtime_error("Can't construct string from anything but a char array.");
  }
  rv = string(mxArrayToString(mxa));
  return rv;
}

template <>
RCP<Xpetra_map> loadDataFromMatlab<RCP<Xpetra_map> >(const mxArray* mxa) {
  RCP<const Teuchos::Comm<int> > comm = rcp(new Teuchos::SerialComm<int>());
  int nr                              = mxGetM(mxa);
  int nc                              = mxGetN(mxa);
  if (nr != 1)
    throw std::runtime_error("A Xpetra::Map representation from MATLAB must be a single row vector.");
  double* pr                    = mxGetPr(mxa);
  mm_GlobalOrd numGlobalIndices = nc;

  std::vector<mm_GlobalOrd> localGIDs(numGlobalIndices);
  for (int i = 0; i < int(numGlobalIndices); i++) {
    localGIDs[i] = Teuchos::as<mm_GlobalOrd>(pr[i]);
  }

  const Teuchos::ArrayView<const mm_GlobalOrd> localGIDs_view(&localGIDs[0], localGIDs.size());
  RCP<Xpetra_map> map =
      Xpetra::MapFactory<mm_LocalOrd, mm_GlobalOrd, mm_node_t>::Build(
          Xpetra::UseTpetra,
          Teuchos::OrdinalTraits<mm_GlobalOrd>::invalid(),
          localGIDs_view,
          0, comm);

  if (map.is_null())
    throw runtime_error("Failed to create Xpetra::Map.");
  return map;
}

template <>
RCP<Xpetra_ordinal_vector> loadDataFromMatlab<RCP<Xpetra_ordinal_vector> >(const mxArray* mxa) {
  RCP<const Teuchos::Comm<int> > comm = rcp(new Teuchos::SerialComm<int>());
  if (mxGetN(mxa) != 1 && mxGetM(mxa) != 1)
    throw std::runtime_error("An OrdinalVector from MATLAB must be a single row or column vector.");
  mm_GlobalOrd numGlobalIndices                               = mxGetM(mxa) * mxGetN(mxa);
  RCP<Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t> > map = Xpetra::MapFactory<mm_LocalOrd, mm_GlobalOrd, mm_node_t>::Build(Xpetra::UseTpetra, numGlobalIndices, 0, comm);
  if (mxGetClassID(mxa) != mxINT32_CLASS)
    throw std::runtime_error("Can only construct LOVector with int32 data.");
  int* array = (int*)mxGetData(mxa);
  if (map.is_null())
    throw runtime_error("Failed to create map for Xpetra ordinal vector.");
  RCP<Xpetra_ordinal_vector> loVec = Xpetra::VectorFactory<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t>::Build(map, false);
  if (loVec.is_null())
    throw runtime_error("Failed to create ordinal vector with Xpetra::VectorFactory.");
  for (int i = 0; i < int(numGlobalIndices); i++) {
    loVec->replaceGlobalValue(i, 0, array[i]);
  }
  return loVec;
}

template <>
RCP<Tpetra_MultiVector_double> loadDataFromMatlab<RCP<Tpetra_MultiVector_double> >(const mxArray* mxa) {
  RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > mv;
  try {
    int nr                              = mxGetM(mxa);
    int nc                              = mxGetN(mxa);
    double* pr                          = mxGetPr(mxa);
    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    // numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
    RCP<const muemex_map_type> map = rcp(new muemex_map_type(nr, (mm_GlobalOrd)0, comm));
    // Allocate a new array of complex values to use with the multivector
    Teuchos::ArrayView<const double> arrView(pr, nr * nc);
    mv = rcp(new Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(map, arrView, size_t(nr), size_t(nc)));
  } catch (std::exception& e) {
    mexPrintf("Error constructing Tpetra MultiVector.\n");
    std::cout << e.what() << std::endl;
  }
  return mv;
}

template <>
RCP<Tpetra_MultiVector_complex> loadDataFromMatlab<RCP<Tpetra_MultiVector_complex> >(const mxArray* mxa) {
  RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > mv;
  try {
    int nr                              = mxGetM(mxa);
    int nc                              = mxGetN(mxa);
    double* pr                          = mxGetPr(mxa);
    double* pi                          = mxGetPi(mxa);
    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    // numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
    RCP<const muemex_map_type> map = rcp(new muemex_map_type(nr, (mm_GlobalOrd)0, comm));
    // Allocate a new array of complex values to use with the multivector
    complex_t* myArr = new complex_t[nr * nc];
    for (int n = 0; n < nc; n++) {
      for (int m = 0; m < nr; m++) {
        myArr[n * nr + m] = complex_t(pr[n * nr + m], pi[n * nr + m]);
      }
    }
    Teuchos::ArrayView<complex_t> arrView(myArr, nr * nc);
    mv = rcp(new Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(map, arrView, nr, nc));
  } catch (std::exception& e) {
    mexPrintf("Error constructing Tpetra MultiVector.\n");
    std::cout << e.what() << std::endl;
  }
  return mv;
}

template <>
RCP<Tpetra_CrsMatrix_double> loadDataFromMatlab<RCP<Tpetra_CrsMatrix_double> >(const mxArray* mxa) {
  bool success = false;
  RCP<Tpetra_CrsMatrix_double> A;

  int* colptr = NULL;
  int* rowind = NULL;

  try {
    RCP<const Teuchos::Comm<int> > comm = rcp(new Teuchos::SerialComm<int>());
    // numGlobalIndices is just the number of rows in the matrix
    const size_t numGlobalIndices        = mxGetM(mxa);
    RCP<const muemex_map_type> rowMap    = rcp(new muemex_map_type(numGlobalIndices, 0, comm));
    RCP<const muemex_map_type> domainMap = rcp(new muemex_map_type(mxGetN(mxa), 0, comm));
    double* valueArray                   = mxGetPr(mxa);
    int nc                               = mxGetN(mxa);
    if (rewrap_ints) {
      // mwIndex_to_int allocates memory so must delete[] later
      colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
      rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
    } else {
      rowind = (int*)mxGetIr(mxa);
      colptr = (int*)mxGetJc(mxa);
    }
    // Need this to convert CSC colptrs to CRS row counts
    Teuchos::Array<size_t> rowCounts(numGlobalIndices);
    for (int i = 0; i < nc; i++) {
      for (int j = colptr[i]; j < colptr[i + 1]; j++) {
        rowCounts[rowind[j]]++;
      }
    }
    A = rcp(new Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(rowMap, rowCounts()));
    for (int i = 0; i < nc; i++) {
      for (int j = colptr[i]; j < colptr[i + 1]; j++) {
        //'array' of 1 element, containing column (in global matrix).
        Teuchos::ArrayView<mm_GlobalOrd> cols = Teuchos::ArrayView<mm_GlobalOrd>(&i, 1);
        //'array' of 1 element, containing value
        Teuchos::ArrayView<double> vals = Teuchos::ArrayView<double>(&valueArray[j], 1);
        A->insertGlobalValues(rowind[j], cols, vals);
      }
    }
    A->fillComplete(domainMap, rowMap);
    if (rewrap_ints) {
      delete[] rowind;
      rowind = NULL;
      delete[] colptr;
      colptr = NULL;
    }
    success = true;
  } catch (std::exception& e) {
    if (rewrap_ints) {
      if (rowind != NULL) delete[] rowind;
      if (colptr != NULL) delete[] colptr;
      rowind = NULL;
      colptr = NULL;
    }
    mexPrintf("Error while constructing Tpetra matrix:\n");
    std::cout << e.what() << std::endl;
  }
  if (!success)
    mexErrMsgTxt("An error occurred while setting up a Tpetra matrix.\n");
  return A;
}

template <>
RCP<Tpetra_CrsMatrix_complex> loadDataFromMatlab<RCP<Tpetra_CrsMatrix_complex> >(const mxArray* mxa) {
  RCP<Tpetra_CrsMatrix_complex> A;
  // Create a map in order to create the matrix (taken from muelu basic example - complex)
  try {
    RCP<const Teuchos::Comm<int> > comm          = Tpetra::getDefaultComm();
    const Tpetra::global_size_t numGlobalIndices = mxGetM(mxa);
    const mm_GlobalOrd indexBase                 = 0;
    RCP<const muemex_map_type> rowMap            = rcp(new muemex_map_type(numGlobalIndices, indexBase, comm));
    RCP<const muemex_map_type> domainMap         = rcp(new muemex_map_type(mxGetN(mxa), indexBase, comm));
    double* realArray                            = mxGetPr(mxa);
    double* imagArray                            = mxGetPi(mxa);
    int* colptr;
    int* rowind;
    int nc = mxGetN(mxa);
    if (rewrap_ints) {
      // mwIndex_to_int allocates memory so must delete[] later
      colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
      rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
    } else {
      rowind = (int*)mxGetIr(mxa);
      colptr = (int*)mxGetJc(mxa);
    }
    // Need this to convert CSC colptrs to CRS row counts
    Teuchos::Array<size_t> rowCounts(numGlobalIndices);
    for (int i = 0; i < nc; i++) {
      for (int j = colptr[i]; j < colptr[i + 1]; j++) {
        rowCounts[rowind[j]]++;
      }
    }
    A = rcp(new Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(rowMap, rowCounts()));
    for (int i = 0; i < nc; i++) {
      for (int j = colptr[i]; j < colptr[i + 1]; j++) {
        // here assuming that complex_t will always be defined as std::complex<double>
        // use 'value' over and over again with Teuchos::ArrayViews to insert into matrix
        complex_t value                       = std::complex<double>(realArray[j], imagArray[j]);
        Teuchos::ArrayView<mm_GlobalOrd> cols = Teuchos::ArrayView<mm_GlobalOrd>(&i, 1);
        Teuchos::ArrayView<complex_t> vals    = Teuchos::ArrayView<complex_t>(&value, 1);
        A->insertGlobalValues(rowind[j], cols, vals);
      }
    }
    A->fillComplete(domainMap, rowMap);
    if (rewrap_ints) {
      delete[] rowind;
      delete[] colptr;
    }
  } catch (std::exception& e) {
    mexPrintf("Error while constructing tpetra matrix:\n");
    std::cout << e.what() << std::endl;
  }
  return A;
}

template <>
RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > loadDataFromMatlab<RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(const mxArray* mxa) {
  RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > tmat = loadDataFromMatlab<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(mxa);
  return MueLu::TpetraCrs_To_XpetraMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tmat);
}

template <>
RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > loadDataFromMatlab<RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(const mxArray* mxa) {
  RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > tmat = loadDataFromMatlab<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(mxa);
  return MueLu::TpetraCrs_To_XpetraMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tmat);
}

template <>
RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > loadDataFromMatlab<RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(const mxArray* mxa) {
  RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > tpetraMV = loadDataFromMatlab<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(mxa);
  return MueLu::TpetraMultiVector_To_XpetraMultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tpetraMV);
}

template <>
RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > loadDataFromMatlab<RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(const mxArray* mxa) {
  RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > tpetraMV = loadDataFromMatlab<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > >(mxa);
  return MueLu::TpetraMultiVector_To_XpetraMultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tpetraMV);
}

template <>
RCP<Epetra_CrsMatrix> loadDataFromMatlab<RCP<Epetra_CrsMatrix> >(const mxArray* mxa) {
  RCP<Epetra_CrsMatrix> matrix;
  try {
    int* colptr;
    int* rowind;
    double* vals = mxGetPr(mxa);
    int nr       = mxGetM(mxa);
    int nc       = mxGetN(mxa);
    if (rewrap_ints) {
      colptr = mwIndex_to_int(nc + 1, mxGetJc(mxa));
      rowind = mwIndex_to_int(colptr[nc], mxGetIr(mxa));
    } else {
      rowind = (int*)mxGetIr(mxa);
      colptr = (int*)mxGetJc(mxa);
    }
    Epetra_SerialComm Comm;
    Epetra_Map RangeMap(nr, 0, Comm);
    Epetra_Map DomainMap(nc, 0, Comm);
    matrix = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RangeMap, DomainMap, 0));
    /* Do the matrix assembly */
    for (int i = 0; i < nc; i++) {
      for (int j = colptr[i]; j < colptr[i + 1]; j++) {
        // global row, # of entries, value array, column indices array
        matrix->InsertGlobalValues(rowind[j], 1, &vals[j], &i);
      }
    }
    matrix->FillComplete(DomainMap, RangeMap);
    if (rewrap_ints) {
      delete[] rowind;
      delete[] colptr;
    }
  } catch (std::exception& e) {
    mexPrintf("An error occurred while setting up an Epetra matrix:\n");
    std::cout << e.what() << std::endl;
  }
  return matrix;
}

template <>
RCP<Epetra_MultiVector> loadDataFromMatlab<RCP<Epetra_MultiVector> >(const mxArray* mxa) {
  int nr = mxGetM(mxa);
  int nc = mxGetN(mxa);
  Epetra_SerialComm Comm;
  Epetra_BlockMap map(nr * nc, 1, 0, Comm);
  return rcp(new Epetra_MultiVector(Epetra_DataAccess::Copy, map, mxGetPr(mxa), nr, nc));
}

template <>
RCP<MAggregates> loadDataFromMatlab<RCP<MAggregates> >(const mxArray* mxa) {
  if (mxGetNumberOfElements(mxa) != 1)
    throw runtime_error("Aggregates must be individual structs in MATLAB.");
  if (!mxIsStruct(mxa))
    throw runtime_error("Trying to pull aggregates from non-struct MATLAB object.");
  // assume that in matlab aggregate structs will only be stored in a 1x1 array
  // mxa must have the same fields as the ones declared in constructAggregates function in muelu.m for this to work
  const int correctNumFields = 5;  // change if more fields are added to the aggregates representation in constructAggregates in muelu.m
  if (mxGetNumberOfFields(mxa) != correctNumFields)
    throw runtime_error("Aggregates structure has wrong number of fields.");
  // Pull MuemexData types back out
  int nVert = *(int*)mxGetData(mxGetField(mxa, 0, "nVertices"));
  int nAgg  = *(int*)mxGetData(mxGetField(mxa, 0, "nAggregates"));
  // Now have all the data needed to fully reconstruct the aggregate
  // Use similar approach as UserAggregationFactory (which is written for >1 thread but will just be serial here)
  RCP<const Teuchos::Comm<int> > comm                         = Teuchos::DefaultComm<int>::getComm();
  int myRank                                                  = comm->getRank();
  Xpetra::UnderlyingLib lib                                   = Xpetra::UseTpetra;
  RCP<Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t> > map = Xpetra::MapFactory<mm_LocalOrd, mm_GlobalOrd, mm_node_t>::Build(lib, nVert, 0, comm);
  RCP<MAggregates> agg                                        = rcp(new MAggregates(map));
  agg->SetNumAggregates(nAgg);
  // Get handles for the vertex2AggId and procwinner arrays in reconstituted aggregates object
  // this is serial so all procwinner values will be same (0)
  ArrayRCP<mm_LocalOrd> vertex2AggId = agg->GetVertex2AggId()->getDataNonConst(0);  // the '0' means first (and only) column of multivector, since is just vector
  ArrayRCP<mm_LocalOrd> procWinner   = agg->GetProcWinner()->getDataNonConst(0);
  // mm_LocalOrd and int are equivalent, so is ok to talk about aggSize with just 'int'
  // Deep copy the entire vertex2AggID and isRoot arrays, which are both nVert items long
  // At the same time, set ProcWinner
  mxArray* vertToAggID_in  = mxGetField(mxa, 0, "vertexToAggID");
  int* vertToAggID_inArray = (int*)mxGetData(vertToAggID_in);
  mxArray* rootNodes_in    = mxGetField(mxa, 0, "rootNodes");
  int* rootNodes_inArray   = (int*)mxGetData(rootNodes_in);
  for (int i = 0; i < nVert; i++) {
    vertex2AggId[i] = vertToAggID_inArray[i];
    procWinner[i]   = myRank;  // all nodes are going to be on the same proc
    agg->SetIsRoot(i, false);  // the ones that are root will be set in next loop
  }
  for (int i = 0; i < nAgg; i++)  // rootNodesToCopy is an array of node IDs which are the roots of their aggs
  {
    agg->SetIsRoot(rootNodes_inArray[i], true);
  }
  // Now recompute the aggSize array the results in the object
  agg->ComputeAggregateSizes(true);
  agg->AggregatesCrossProcessors(false);
  return agg;
}

template <>
RCP<MAmalInfo> loadDataFromMatlab<RCP<MAmalInfo> >(const mxArray* mxa) {
  RCP<MAmalInfo> amal;
  throw runtime_error("AmalgamationInfo not supported in Muemex yet.");
  return amal;
}

template <>
RCP<MGraph> loadDataFromMatlab<RCP<MGraph> >(const mxArray* mxa) {
  // mxa must be struct with logical sparse matrix called 'edges' and Nx1 int32 array 'boundaryNodes'
  mxArray* edges         = mxGetField(mxa, 0, "edges");
  mxArray* boundaryNodes = mxGetField(mxa, 0, "boundaryNodes");
  if (edges == NULL)
    throw runtime_error("Graph structure in MATLAB must have a field called 'edges' (logical sparse matrix)");
  if (boundaryNodes == NULL)
    throw runtime_error("Graph structure in MATLAB must have a field called 'boundaryNodes' (int32 array containing list of boundary nodes)");
  int* boundaryList = (int*)mxGetData(boundaryNodes);
  if (!mxIsSparse(edges) || mxGetClassID(edges) != mxLOGICAL_CLASS)
    throw runtime_error("Graph edges must be stored as a logical sparse matrix.");
  // Note that Matlab stores sparse matrices in column major format.
  mwIndex* matlabColPtrs    = mxGetJc(edges);
  mwIndex* matlabRowIndices = mxGetIr(edges);
  mm_GlobalOrd nRows        = (mm_GlobalOrd)mxGetM(edges);

  // Create and populate row-major CRS data structures for Xpetra::TpetraCrsGraph.

  // calculate number of nonzeros in each row
  Teuchos::Array<int> entriesPerRow(nRows);
  int nnz = matlabColPtrs[mxGetN(edges)];  // last entry in matlabColPtrs
  for (int i = 0; i < nnz; i++)
    entriesPerRow[matlabRowIndices[i]]++;
  // Populate usual row index array.  We don't need this for the Xpetra Graph ctor, but
  // it's convenient for building up the column index array, which the ctor does need.
  Teuchos::Array<int> rows(nRows + 1);
  rows[0] = 0;
  for (int i = 0; i < nRows; i++)
    rows[i + 1] = rows[i] + entriesPerRow[i];
  Teuchos::Array<int> cols(nnz);                   // column index array
  Teuchos::Array<int> insertionsPerRow(nRows, 0);  // track of #insertions done per row
  int ncols = mxGetN(edges);
  for (int colNum = 0; colNum < ncols; ++colNum) {
    int ci = matlabColPtrs[colNum];
    for (int j = ci; j < Teuchos::as<int>(matlabColPtrs[colNum + 1]); ++j) {
      int rowNum                                    = matlabRowIndices[j];
      cols[rows[rowNum] + insertionsPerRow[rowNum]] = colNum;
      insertionsPerRow[rowNum]++;
    }
  }
  // Find maximum
  int maxNzPerRow = 0;
  for (int i = 0; i < nRows; i++) {
    if (maxNzPerRow < entriesPerRow[i])
      maxNzPerRow = entriesPerRow[i];
  }

  RCP<const Teuchos::Comm<mm_GlobalOrd> > comm = rcp(new Teuchos::SerialComm<mm_GlobalOrd>());
  typedef Xpetra::TpetraMap<mm_LocalOrd, mm_GlobalOrd, mm_node_t> MMap;
  RCP<MMap> map = rcp(new MMap(nRows, 0, comm));
  typedef Xpetra::TpetraCrsGraph<mm_LocalOrd, mm_GlobalOrd, mm_node_t> TpetraGraph;
  RCP<TpetraGraph> tgraph = rcp(new TpetraGraph(map, (size_t)maxNzPerRow));
  // Populate tgraph in compressed-row format. Must get each row individually...
  for (int i = 0; i < nRows; ++i) {
    tgraph->insertGlobalIndices((mm_GlobalOrd)i, cols(rows[i], entriesPerRow[i]));
  }
  tgraph->fillComplete(map, map);
  RCP<MGraph> mgraph = rcp(new MueLu::Graph<mm_LocalOrd, mm_GlobalOrd, mm_node_t>(tgraph));
  // Set boundary nodes
  int numBoundaryNodes = mxGetNumberOfElements(boundaryNodes);
  bool* boundaryFlags  = new bool[nRows];
  for (int i = 0; i < nRows; i++) {
    boundaryFlags[i] = false;
  }
  for (int i = 0; i < numBoundaryNodes; i++) {
    boundaryFlags[boundaryList[i]] = true;
  }
  ArrayRCP<bool> boundaryNodesInput(boundaryFlags, 0, nRows, true);
  mgraph->SetBoundaryNodeMap(boundaryNodesInput);
  return mgraph;
}

#ifdef HAVE_MUELU_INTREPID2
template <>
RCP<FieldContainer_ordinal> loadDataFromMatlab<RCP<FieldContainer_ordinal> >(const mxArray* mxa) {
  if (mxGetClassID(mxa) != mxINT32_CLASS)
    throw runtime_error("FieldContainer must have integer storage entries");

  int* data = (int*)mxGetData(mxa);
  int nr    = mxGetM(mxa);
  int nc    = mxGetN(mxa);

  RCP<FieldContainer_ordinal> fc = rcp(new FieldContainer_ordinal("FC from Matlab", nr, nc));
  for (int col = 0; col < nc; col++) {
    for (int row = 0; row < nr; row++) {
      (*fc)(row, col) = data[col * nr + row];
    }
  }
  return fc;
}
#endif

/* ******************************* */
/* saveDataToMatlab                */
/* ******************************* */

template <>
mxArray* saveDataToMatlab(int& data) {
  mwSize dims[]           = {1, 1};
  mxArray* mxa            = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  *((int*)mxGetData(mxa)) = data;
  return mxa;
}

template <>
mxArray* saveDataToMatlab(bool& data) {
  mwSize dims[]            = {1, 1};
  mxArray* mxa             = mxCreateLogicalArray(2, dims);
  *((bool*)mxGetData(mxa)) = data;
  return mxa;
}

template <>
mxArray* saveDataToMatlab(double& data) {
  return mxCreateDoubleScalar(data);
}

template <>
mxArray* saveDataToMatlab(complex_t& data) {
  mwSize dims[]            = {1, 1};
  mxArray* mxa             = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  *((double*)mxGetPr(mxa)) = real<double>(data);
  *((double*)mxGetPi(mxa)) = imag<double>(data);
  return mxa;
}

template <>
mxArray* saveDataToMatlab(string& data) {
  return mxCreateString(data.c_str());
}

template <>
mxArray* saveDataToMatlab(RCP<Xpetra_map>& data) {
  // Precondition: Memory has already been allocated by MATLAB for the array.
  int nc          = data->getGlobalNumElements();
  int nr          = 1;
  mxArray* output = createMatlabMultiVector<double>(nr, nc);
  double* array   = (double*)malloc(sizeof(double) * nr * nc);
  for (int col = 0; col < nc; col++) {
    mm_GlobalOrd gid = data->getGlobalElement(col);
    array[col]       = Teuchos::as<double>(gid);
  }
  fillMatlabArray<double>(array, output, nc * nr);
  free(array);
  return output;
}

template <>
mxArray* saveDataToMatlab(RCP<Xpetra_ordinal_vector>& data) {
  mwSize len = data->getGlobalLength();
  // create a single column vector
  mwSize dimensions[]             = {len, 1};
  mxArray* rv                     = mxCreateNumericArray(2, dimensions, mxINT32_CLASS, mxREAL);
  int* dataPtr                    = (int*)mxGetData(rv);
  ArrayRCP<const mm_LocalOrd> arr = data->getData(0);
  for (int i = 0; i < int(data->getGlobalLength()); i++) {
    dataPtr[i] = arr[i];
  }
  return rv;
}

template <>
mxArray* saveDataToMatlab(RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> >& data) {
  RCP<Xpetra_MultiVector_double> xmv = MueLu::TpetraMultiVector_To_XpetraMultiVector(data);
  return saveDataToMatlab(xmv);
}

template <>
mxArray* saveDataToMatlab(RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> >& data) {
  RCP<Xpetra_MultiVector_complex> xmv = MueLu::TpetraMultiVector_To_XpetraMultiVector(data);
  return saveDataToMatlab(xmv);
}

template <>
mxArray* saveDataToMatlab(RCP<Tpetra_CrsMatrix_double>& data) {
  RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > xmat = TpetraCrs_To_XpetraMatrix(data);
  return saveDataToMatlab(xmat);
}

template <>
mxArray* saveDataToMatlab(RCP<Tpetra_CrsMatrix_complex>& data) {
  RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > xmat = TpetraCrs_To_XpetraMatrix(data);
  return saveDataToMatlab(xmat);
}

template <>
mxArray* saveDataToMatlab(RCP<Xpetra_Matrix_double>& data) {
  typedef double Scalar;
  // Compute global constants, if we need them
  Teuchos::rcp_const_cast<Xpetra_CrsGraph>(data->getCrsGraph())->computeGlobalConstants();

  int nr                                 = data->getGlobalNumRows();
  int nc                                 = data->getGlobalNumCols();
  int nnz                                = data->getGlobalNumEntries();

#ifdef VERBOSE_OUTPUT
  RCP<Teuchos::FancyOStream> fancyStream = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  mat->describe(*fancyStream, Teuchos::VERB_EXTREME);
#endif
  mxArray* mxa = createMatlabSparse<Scalar>(nr, nc, nnz);
  mwIndex* ir  = mxGetIr(mxa);
  mwIndex* jc  = mxGetJc(mxa);
  for (int i = 0; i < nc + 1; i++) {
    jc[i] = 0;
  }

  size_t maxEntriesPerRow = data->getGlobalMaxNumRowEntries();
  if (maxEntriesPerRow == Teuchos::OrdinalTraits<size_t>::invalid() || maxEntriesPerRow == 0) maxEntriesPerRow = data->getLocalMaxNumRowEntries();

  int* rowProgress = new int[nc];
  // The array that will be copied to Pr and (if complex) Pi later
  Scalar* sparseVals = new Scalar[nnz];
  size_t numEntries;
  if (data->isLocallyIndexed()) {
    Scalar* rowValArray = new Scalar[maxEntriesPerRow];
    Teuchos::ArrayView<Scalar> rowVals(rowValArray, maxEntriesPerRow);
    mm_LocalOrd* rowIndicesArray = new mm_LocalOrd[maxEntriesPerRow];
    Teuchos::ArrayView<mm_LocalOrd> rowIndices(rowIndicesArray, maxEntriesPerRow);
    for (mm_LocalOrd m = 0; m < nr; m++)  // All rows in the Xpetra matrix
    {
      data->getLocalRowCopy(m, rowIndices, rowVals, numEntries);     // Get the row
      for (mm_LocalOrd entry = 0; entry < int(numEntries); entry++)  // All entries in row
      {
        jc[rowIndices[entry] + 1]++;  // for each entry, increase jc for the entry's column
      }
    }

    // now jc holds the number of elements in each column, but needs cumulative sum over all previous columns also
    int entriesAccum = 0;
    for (int n = 0; n <= nc; n++) {
      int temp = entriesAccum;
      entriesAccum += jc[n];
      jc[n] += temp;
    }
    // Jc now populated with colptrs
    for (int i = 0; i < nc; i++) {
      rowProgress[i] = 0;
    }
    // Row progress values like jc but keep track as the MATLAB matrix is being filled in
    for (mm_LocalOrd m = 0; m < nr; m++)  // rows
    {
      data->getLocalRowCopy(m, rowIndices, rowVals, numEntries);
      for (mm_LocalOrd i = 0; i < int(numEntries); i++)  // entries in row m (NOT columns)
      {
        // row is m, col is rowIndices[i], val is rowVals[i]
        mm_LocalOrd col                        = rowIndices[i];
        sparseVals[jc[col] + rowProgress[col]] = rowVals[i];  // Set value
        ir[jc[col] + rowProgress[col]]         = m;           // Set row at which value occurs
        rowProgress[col]++;
      }
    }
    delete[] rowIndicesArray;
  } else {
    Teuchos::ArrayView<const mm_GlobalOrd> rowIndices;
    Teuchos::ArrayView<const Scalar> rowVals;
    for (mm_GlobalOrd m = 0; m < nr; m++) {
      data->getGlobalRowView(m, rowIndices, rowVals);
      for (mm_GlobalOrd n = 0; n < rowIndices.size(); n++) {
        jc[rowIndices[n] + 1]++;
      }
    }
    // Last element of jc is just nnz
    jc[nc] = nnz;
    // Jc now populated with colptrs
    for (int i = 0; i < nc; i++) {
      rowProgress[i] = 0;
    }
    int entriesAccum = 0;
    for (int n = 0; n <= nc; n++) {
      int temp = entriesAccum;
      entriesAccum += jc[n];
      jc[n] += temp;
    }
    // Row progress values like jc but keep track as the MATLAB matrix is being filled in
    for (mm_GlobalOrd m = 0; m < nr; m++)  // rows
    {
      data->getGlobalRowView(m, rowIndices, rowVals);
      for (mm_LocalOrd i = 0; i < rowIndices.size(); i++)  // entries in row m
      {
        mm_GlobalOrd col                       = rowIndices[i];  // row is m, col is rowIndices[i], val is rowVals[i]
        sparseVals[jc[col] + rowProgress[col]] = rowVals[i];     // Set value
        ir[jc[col] + rowProgress[col]]         = m;              // Set row at which value occurs
        rowProgress[col]++;
      }
    }
  }
  // finally, copy sparseVals into pr (and pi, if complex)
  fillMatlabArray<Scalar>(sparseVals, mxa, nnz);
  delete[] sparseVals;
  delete[] rowProgress;
  return mxa;
}

template <>
mxArray* saveDataToMatlab(RCP<Xpetra_Matrix_complex>& data) {
  typedef complex_t Scalar;

  // Compute global constants, if we need them
  Teuchos::rcp_const_cast<Xpetra_CrsGraph>(data->getCrsGraph())->computeGlobalConstants();

  int nr                                 = data->getGlobalNumRows();
  int nc                                 = data->getGlobalNumCols();
  int nnz                                = data->getGlobalNumEntries();
#ifdef VERBOSE_OUTPUT
  RCP<Teuchos::FancyOStream> fancyStream = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  mat->describe(*fancyStream, Teuchos::VERB_EXTREME);
#endif
  mxArray* mxa = createMatlabSparse<Scalar>(nr, nc, nnz);
  mwIndex* ir  = mxGetIr(mxa);
  mwIndex* jc  = mxGetJc(mxa);
  for (int i = 0; i < nc + 1; i++) {
    jc[i] = 0;
  }
  size_t maxEntriesPerRow = data->getGlobalMaxNumRowEntries();
  int* rowProgress        = new int[nc];
  // The array that will be copied to Pr and (if complex) Pi later
  Scalar* sparseVals = new Scalar[nnz];
  size_t numEntries;
  if (data->isLocallyIndexed()) {
    Scalar* rowValArray = new Scalar[maxEntriesPerRow];
    Teuchos::ArrayView<Scalar> rowVals(rowValArray, maxEntriesPerRow);
    mm_LocalOrd* rowIndicesArray = new mm_LocalOrd[maxEntriesPerRow];
    Teuchos::ArrayView<mm_LocalOrd> rowIndices(rowIndicesArray, maxEntriesPerRow);
    for (mm_LocalOrd m = 0; m < nr; m++)  // All rows in the Xpetra matrix
    {
      data->getLocalRowCopy(m, rowIndices, rowVals, numEntries);     // Get the row
      for (mm_LocalOrd entry = 0; entry < int(numEntries); entry++)  // All entries in row
      {
        jc[rowIndices[entry] + 1]++;  // for each entry, increase jc for the entry's column
      }
    }
    // now jc holds the number of elements in each column, but needs cumulative sum over all previous columns also
    int entriesAccum = 0;
    for (int n = 0; n <= nc; n++) {
      int temp = entriesAccum;
      entriesAccum += jc[n];
      jc[n] += temp;
    }
    // Jc now populated with colptrs
    for (int i = 0; i < nc; i++) {
      rowProgress[i] = 0;
    }
    // Row progress values like jc but keep track as the MATLAB matrix is being filled in
    for (mm_LocalOrd m = 0; m < nr; m++)  // rows
    {
      data->getLocalRowCopy(m, rowIndices, rowVals, numEntries);
      for (mm_LocalOrd i = 0; i < int(numEntries); i++)  // entries in row m (NOT columns)
      {
        // row is m, col is rowIndices[i], val is rowVals[i]
        mm_LocalOrd col                        = rowIndices[i];
        sparseVals[jc[col] + rowProgress[col]] = rowVals[i];  // Set value
        ir[jc[col] + rowProgress[col]]         = m;           // Set row at which value occurs
        rowProgress[col]++;
      }
    }
    delete[] rowIndicesArray;
  } else {
    Teuchos::ArrayView<const mm_GlobalOrd> rowIndices;
    Teuchos::ArrayView<const Scalar> rowVals;
    for (mm_GlobalOrd m = 0; m < nr; m++) {
      data->getGlobalRowView(m, rowIndices, rowVals);
      for (mm_GlobalOrd n = 0; n < rowIndices.size(); n++) {
        jc[rowIndices[n] + 1]++;
      }
    }
    // Last element of jc is just nnz
    jc[nc] = nnz;
    // Jc now populated with colptrs
    for (int i = 0; i < nc; i++) {
      rowProgress[i] = 0;
    }
    int entriesAccum = 0;
    for (int n = 0; n <= nc; n++) {
      int temp = entriesAccum;
      entriesAccum += jc[n];
      jc[n] += temp;
    }
    // Row progress values like jc but keep track as the MATLAB matrix is being filled in
    for (mm_GlobalOrd m = 0; m < nr; m++)  // rows
    {
      data->getGlobalRowView(m, rowIndices, rowVals);
      for (mm_LocalOrd i = 0; i < rowIndices.size(); i++)  // entries in row m
      {
        mm_GlobalOrd col                       = rowIndices[i];  // row is m, col is rowIndices[i], val is rowVals[i]
        sparseVals[jc[col] + rowProgress[col]] = rowVals[i];     // Set value
        ir[jc[col] + rowProgress[col]]         = m;              // Set row at which value occurs
        rowProgress[col]++;
      }
    }
  }
  // finally, copy sparseVals into pr (and pi, if complex)
  fillMatlabArray<Scalar>(sparseVals, mxa, nnz);
  delete[] sparseVals;
  delete[] rowProgress;
  return mxa;
}

/*
template<>
mxArray* saveDataToMatlab(RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>& data)
{
  //Precondition: Memory has already been allocated by MATLAB for the array.
  int nr = data->getGlobalLength();
  int nc = data->getNumVectors();
  mxArray* output = createMatlabMultiVector<Scalar>(nr, nc);
  Scalar* array = (Scalar*) malloc(sizeof(Scalar) * nr * nc);
  for(int col = 0; col < nc; col++)
  {
    Teuchos::ArrayRCP<const Scalar> colData = data->getData(col);
    for(int row = 0; row < nr; row++)
    {
      array[col * nr + row] = colData[row];
    }
  }
  fillMatlabArray<Scalar>(array, output, nc * nr);
  free(array);
  return output;
}
*/

template <>
mxArray* saveDataToMatlab(RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> >& data) {
  // Precondition: Memory has already been allocated by MATLAB for the array.
  int nr          = data->getGlobalLength();
  int nc          = data->getNumVectors();
  mxArray* output = createMatlabMultiVector<double>(nr, nc);
  double* array   = (double*)malloc(sizeof(double) * nr * nc);
  for (int col = 0; col < nc; col++) {
    Teuchos::ArrayRCP<const double> colData = data->getData(col);
    for (int row = 0; row < nr; row++) {
      array[col * nr + row] = colData[row];
    }
  }
  fillMatlabArray<double>(array, output, nc * nr);
  free(array);
  return output;
}

template <>
mxArray* saveDataToMatlab(RCP<Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> >& data) {
  // Precondition: Memory has already been allocated by MATLAB for the array.
  int nr           = data->getGlobalLength();
  int nc           = data->getNumVectors();
  mxArray* output  = createMatlabMultiVector<complex_t>(nr, nc);
  complex_t* array = (complex_t*)malloc(sizeof(complex_t) * nr * nc);
  for (int col = 0; col < nc; col++) {
    Teuchos::ArrayRCP<const complex_t> colData = data->getData(col);
    for (int row = 0; row < nr; row++) {
      array[col * nr + row] = colData[row];
    }
  }
  fillMatlabArray<complex_t>(array, output, nc * nr);
  free(array);
  return output;
}

template <>
mxArray* saveDataToMatlab(RCP<Epetra_CrsMatrix>& data) {
  RCP<Xpetra_Matrix_double> xmat = EpetraCrs_To_XpetraMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(data);
  return saveDataToMatlab(xmat);
}

template <>
mxArray* saveDataToMatlab(RCP<Epetra_MultiVector>& data) {
  mxArray* output = mxCreateDoubleMatrix(data->GlobalLength(), data->NumVectors(), mxREAL);
  double* dataPtr = mxGetPr(output);
  data->ExtractCopy(dataPtr, data->GlobalLength());
  return output;
}

template <>
mxArray* saveDataToMatlab(RCP<MAggregates>& data) {
  // Set up array of inputs for matlab constructAggregates
  int numNodes = data->GetVertex2AggId()->getData(0).size();
  int numAggs  = data->GetNumAggregates();
  mxArray* dataIn[5];
  mwSize singleton[]                        = {1, 1};
  dataIn[0]                                 = mxCreateNumericArray(2, singleton, mxINT32_CLASS, mxREAL);
  *((int*)mxGetData(dataIn[0]))             = numNodes;
  dataIn[1]                                 = mxCreateNumericArray(2, singleton, mxINT32_CLASS, mxREAL);
  *((int*)mxGetData(dataIn[1]))             = numAggs;
  mwSize nodeArrayDims[]                    = {(mwSize)numNodes, 1};  // dimensions for Nx1 array, where N is number of nodes (vert2Agg)
  dataIn[2]                                 = mxCreateNumericArray(2, nodeArrayDims, mxINT32_CLASS, mxREAL);
  int* vtaid                                = (int*)mxGetData(dataIn[2]);
  ArrayRCP<const mm_LocalOrd> vertexToAggID = data->GetVertex2AggId()->getData(0);
  for (int i = 0; i < numNodes; i++) {
    vtaid[i] = vertexToAggID[i];
  }
  mwSize aggArrayDims[] = {(mwSize)numAggs, 1};  // dims for Nx1 array, where N is number of aggregates (rootNodes, aggSizes)
  dataIn[3]             = mxCreateNumericArray(2, aggArrayDims, mxINT32_CLASS, mxREAL);
  // First, find out if the aggregates even have 1 root node per aggregate. If not, assume roots are invalid and assign ourselves
  int totalRoots = 0;
  for (int i = 0; i < numNodes; i++) {
    if (data->IsRoot(i))
      totalRoots++;
  }
  bool reassignRoots = false;
  if (totalRoots != numAggs) {
    cout << endl
         << "Warning: Number of root nodes and number of aggregates do not match." << endl;
    cout << "Will reassign root nodes when writing aggregates to matlab." << endl
         << endl;
    reassignRoots = true;
  }
  int* rn = (int*)mxGetData(dataIn[3]);  // list of root nodes (in no particular order)
  {
    if (reassignRoots) {
      // For each aggregate, just pick the first node we see in it and set it as root
      int lastFoundNode = 0;  // heuristic for speed, a node in aggregate N+1 is likely to come very soon after a node in agg N
      for (int i = 0; i < numAggs; i++) {
        rn[i] = -1;
        for (int j = lastFoundNode; j < lastFoundNode + numNodes; j++) {
          int index = j % numNodes;
          if (vertexToAggID[index] == i) {
            rn[i]         = index;
            lastFoundNode = index;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(rn[i] == -1, runtime_error, "Invalid aggregates: Couldn't find any node in aggregate #" << i << ".");
      }
    } else {
      int i = 0;  // iterates over aggregate IDs
      for (int j = 0; j < numNodes; j++) {
        if (data->IsRoot(j)) {
          if (i == numAggs)
            throw runtime_error("Cannot store invalid aggregates in MATLAB - more root nodes than aggregates.");
          rn[i] = j;  // now we know this won't go out of bounds (rn's underlying matlab array is numAggs in length)
          i++;
        }
      }
      if (i + 1 < numAggs)
        throw runtime_error("Cannot store invalid aggregates in MATLAB - fewer root nodes than aggregates.");
    }
  }
  dataIn[4]                      = mxCreateNumericArray(1, aggArrayDims, mxINT32_CLASS, mxREAL);
  int* as                        = (int*)mxGetData(dataIn[4]);  // list of aggregate sizes
  ArrayRCP<mm_LocalOrd> aggSizes = data->ComputeAggregateSizes();
  for (int i = 0; i < numAggs; i++) {
    as[i] = aggSizes[i];
  }
  mxArray* matlabAggs[1];
  int result = mexCallMATLAB(1, matlabAggs, 5, dataIn, "constructAggregates");
  if (result != 0)
    throw runtime_error("Matlab encountered an error while constructing aggregates struct.");
  return matlabAggs[0];
}

template <>
mxArray* saveDataToMatlab(RCP<MAmalInfo>& data) {
  throw runtime_error("AmalgamationInfo not supported in MueMex yet.");
  return mxCreateDoubleScalar(0);
}

template <>
mxArray* saveDataToMatlab(RCP<MGraph>& data) {
  int numEntries        = (int)data->GetGlobalNumEdges();
  int numRows           = (int)data->GetDomainMap()->getGlobalNumElements();  // assume numRows == numCols
  mxArray* mat          = mxCreateSparseLogicalMatrix(numRows, numRows, numEntries);
  mxLogical* outData    = (mxLogical*)mxGetData(mat);
  mwIndex* rowInds      = mxGetIr(mat);
  mwIndex* colPtrs      = mxGetJc(mat);
  mm_LocalOrd* dataCopy = new mm_LocalOrd[numEntries];
  mm_LocalOrd* iter     = dataCopy;
  int* entriesPerRow    = new int[numRows];
  int* entriesPerCol    = new int[numRows];
  for (int i = 0; i < numRows; i++) {
    entriesPerRow[i] = 0;
    entriesPerCol[i] = 0;
  }
  for (int i = 0; i < numRows; i++) {
    ArrayView<const mm_LocalOrd> neighbors = data->getNeighborVertices(i);  // neighbors has the column indices for row i
    memcpy(iter, neighbors.getRawPtr(), sizeof(mm_LocalOrd) * neighbors.size());
    entriesPerRow[i] = neighbors.size();
    for (int j = 0; j < neighbors.size(); j++) {
      entriesPerCol[neighbors[j]]++;
    }
    iter += neighbors.size();
  }
  mwIndex** rowIndsByColumn  = new mwIndex*[numRows];  // rowIndsByColumn[0] points to array of row indices in column 1
  mxLogical** valuesByColumn = new mxLogical*[numRows];
  int* numEnteredPerCol      = new int[numRows];
  int accum                  = 0;
  for (int i = 0; i < numRows; i++) {
    rowIndsByColumn[i] = &rowInds[accum];
    // cout << "Entries in column " << i << " start at offset " << accum << endl;
    valuesByColumn[i] = &outData[accum];
    accum += entriesPerCol[i];
    if (accum > numEntries)
      throw runtime_error("potato");
  }
  for (int i = 0; i < numRows; i++) {
    numEnteredPerCol[i] = 0;  // rowIndsByColumn[n][numEnteredPerCol[n]] gives the next place to put a row index
  }
  // entriesPerCol now has Jc information (col offsets)
  accum = 0;  // keep track of cumulative index in dataCopy
  for (int row = 0; row < numRows; row++) {
    for (int entryInRow = 0; entryInRow < entriesPerRow[row]; entryInRow++) {
      int col = dataCopy[accum];
      accum++;
      rowIndsByColumn[col][numEnteredPerCol[col]] = row;
      valuesByColumn[col][numEnteredPerCol[col]]  = (mxLogical)1;
      numEnteredPerCol[col]++;
    }
  }
  accum = 0;  // keep track of total entries over all columns
  for (int col = 0; col < numRows; col++) {
    colPtrs[col] = accum;
    accum += entriesPerCol[col];
  }
  colPtrs[numRows] = accum;  // the last entry in jc, which is equivalent to numEntries
  delete[] numEnteredPerCol;
  delete[] rowIndsByColumn;
  delete[] valuesByColumn;
  delete[] dataCopy;
  delete[] entriesPerRow;
  delete[] entriesPerCol;
  // Construct list of boundary nodes
  const ArrayRCP<const bool> boundaryFlags = data->GetBoundaryNodeMap();
  int numBoundaryNodes                     = 0;
  for (int i = 0; i < boundaryFlags.size(); i++) {
    if (boundaryFlags[i])
      numBoundaryNodes++;
  }
  cout << "Graph has " << numBoundaryNodes << " Dirichlet boundary nodes." << endl;
  mwSize dims[]         = {(mwSize)numBoundaryNodes, 1};
  mxArray* boundaryList = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  int* dest             = (int*)mxGetData(boundaryList);
  int* destIter         = dest;
  for (int i = 0; i < boundaryFlags.size(); i++) {
    if (boundaryFlags[i]) {
      *destIter = i;
      destIter++;
    }
  }
  mxArray* constructArgs[] = {mat, boundaryList};
  mxArray* out[1];
  mexCallMATLAB(1, out, 2, constructArgs, "constructGraph");
  return out[0];
}

#ifdef HAVE_MUELU_INTREPID2
template <>
mxArray* saveDataToMatlab(RCP<FieldContainer_ordinal>& data) {
  int rank = data->rank();
  // NOTE: Only supports rank 2 arrays
  if (rank != 2)
    throw std::runtime_error("Error: Only rank two FieldContainers are supported.");

  int nr = data->extent(0);
  int nc = data->extent(1);

  mwSize dims[] = {(mwSize)nr, (mwSize)nc};
  mxArray* mxa  = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  int* array    = (int*)mxGetData(mxa);

  for (int col = 0; col < nc; col++) {
    for (int row = 0; row < nr; row++) {
      array[col * nr + row] = (*data)(row, col);
    }
  }
  return mxa;
}
#endif

template <typename T>
MuemexData<T>::MuemexData(const mxArray* mxa)
  : MuemexArg(getMuemexType<T>()) {
  data = loadDataFromMatlab<T>(mxa);
}

template <typename T>
mxArray* MuemexData<T>::convertToMatlab() {
  return saveDataToMatlab<T>(data);
}

template <typename T>
MuemexData<T>::MuemexData(T& dataToCopy, MuemexType dataType)
  : MuemexArg(dataType) {
  data = dataToCopy;
}

template <typename T>
MuemexData<T>::MuemexData(T& dataToCopy)
  : MuemexData(dataToCopy, getMuemexType(dataToCopy)) {}

template <typename T>
T& MuemexData<T>::getData() {
  return data;
}

template <typename T>
void MuemexData<T>::setData(T& newData) {
  this->data = newData;
}

/* ***************************** */
/* More Template Functions       */
/* ***************************** */

template <typename T>
void addLevelVariable(const T& data, std::string& name, Level& lvl, const Factory* fact) {
  lvl.AddKeepFlag(name, fact, MueLu::UserData);
  lvl.Set<T>(name, data, fact);
}

template <typename T>
const T& getLevelVariable(std::string& name, Level& lvl) {
  try {
    return lvl.Get<T>(name);
  } catch (std::exception& e) {
    throw std::runtime_error("Requested custom variable " + name + " is not in the level.");
  }
}

// Functions used to put data through matlab factories - first arg is "this" pointer of matlab factory
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
std::vector<Teuchos::RCP<MuemexArg> > processNeeds(const Factory* factory, std::string& needsParam, Level& lvl) {
  using namespace std;
  using namespace Teuchos;
  typedef RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Matrix_t;
  typedef RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MultiVector_t;
  typedef RCP<Aggregates<LocalOrdinal, GlobalOrdinal, Node> > Aggregates_t;
  typedef RCP<AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node> > AmalgamationInfo_t;
  typedef RCP<MGraph> Graph_t;
  vector<string> needsList = tokenizeList(needsParam);
  vector<RCP<MuemexArg> > args;
  for (size_t i = 0; i < needsList.size(); i++) {
    if (needsList[i] == "A" || needsList[i] == "P" || needsList[i] == "R" || needsList[i] == "Ptent") {
      Matrix_t mydata = lvl.Get<Matrix_t>(needsList[i], factory->GetFactory(needsList[i]).get());
      args.push_back(rcp(new MuemexData<Matrix_t>(mydata)));
    } else if (needsList[i] == "Nullspace" || needsList[i] == "Coordinates") {
      MultiVector_t mydata = lvl.Get<MultiVector_t>(needsList[i], factory->GetFactory(needsList[i]).get());
      args.push_back(rcp(new MuemexData<MultiVector_t>(mydata)));
    } else if (needsList[i] == "Aggregates") {
      Aggregates_t mydata = lvl.Get<Aggregates_t>(needsList[i], factory->GetFactory(needsList[i]).get());
      args.push_back(rcp(new MuemexData<Aggregates_t>(mydata)));
    } else if (needsList[i] == "UnAmalgamationInfo") {
      AmalgamationInfo_t mydata = lvl.Get<AmalgamationInfo_t>(needsList[i], factory->GetFactory(needsList[i]).get());
      args.push_back(rcp(new MuemexData<AmalgamationInfo_t>(mydata)));
    } else if (needsList[i] == "Level") {
      int levelNum = lvl.GetLevelID();
      args.push_back(rcp(new MuemexData<int>(levelNum)));
    } else if (needsList[i] == "Graph") {
      Graph_t mydata = lvl.Get<Graph_t>(needsList[i], factory->GetFactory(needsList[i]).get());
      args.push_back(rcp(new MuemexData<Graph_t>(mydata)));
    } else {
      vector<string> words;
      string badNameMsg = "Custom Muemex variables in \"Needs\" list require a type and a name, e.g. \"double myVal\". \n Leading and trailing spaces are OK. \n Don't know how to handle \"" + needsList[i] + "\".\n";
      // compare type without case sensitivity
      char* buf = (char*)malloc(needsList[i].size() + 1);
      strcpy(buf, needsList[i].c_str());
      for (char* iter = buf; *iter != ' '; iter++) {
        if (*iter == 0) {
          free(buf);
          throw runtime_error(badNameMsg);
        }
        *iter = (char)tolower(*iter);
      }
      const char* wordDelim = " ";
      char* mark            = strtok(buf, wordDelim);
      while (mark != NULL) {
        string wordStr(mark);
        words.push_back(wordStr);
        mark = strtok(NULL, wordDelim);
      }
      if (words.size() != 2) {
        free(buf);
        throw runtime_error(badNameMsg);
      }
      char* typeStr = (char*)words[0].c_str();
      if (strstr(typeStr, "ordinalvector")) {
        typedef RCP<Xpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > LOVector_t;
        LOVector_t mydata = getLevelVariable<LOVector_t>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<LOVector_t>(mydata)));
      } else if (strstr(typeStr, "map")) {
        typedef RCP<Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t> > Map_t;
        Map_t mydata = getLevelVariable<Map_t>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<Map_t>(mydata)));
      } else if (strstr(typeStr, "scalar")) {
        Scalar mydata = getLevelVariable<Scalar>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<Scalar>(mydata)));
      } else if (strstr(typeStr, "double")) {
        double mydata = getLevelVariable<double>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<double>(mydata)));
      } else if (strstr(typeStr, "complex")) {
        complex_t mydata = getLevelVariable<complex_t>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<complex_t>(mydata)));
      } else if (strstr(typeStr, "matrix")) {
        Matrix_t mydata = getLevelVariable<Matrix_t>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<Matrix_t>(mydata)));
      } else if (strstr(typeStr, "multivector")) {
        MultiVector_t mydata = getLevelVariable<MultiVector_t>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<MultiVector_t>(mydata)));
      } else if (strstr(typeStr, "int")) {
        int mydata = getLevelVariable<int>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<int>(mydata)));
      } else if (strstr(typeStr, "string")) {
        string mydata = getLevelVariable<string>(needsList[i], lvl);
        args.push_back(rcp(new MuemexData<string>(mydata)));
      } else {
        free(buf);
        throw std::runtime_error(words[0] + " is not a known variable type.");
      }
      free(buf);
    }
  }
  return args;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void processProvides(std::vector<Teuchos::RCP<MuemexArg> >& mexOutput, const Factory* factory, std::string& providesParam, Level& lvl) {
  using namespace std;
  using namespace Teuchos;
  typedef RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Matrix_t;
  typedef RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > MultiVector_t;
  typedef RCP<Aggregates<LocalOrdinal, GlobalOrdinal, Node> > Aggregates_t;
  typedef RCP<AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node> > AmalgamationInfo_t;
  typedef RCP<MGraph> Graph_t;
  vector<string> provides = tokenizeList(providesParam);
  for (size_t i = 0; i < size_t(provides.size()); i++) {
    if (provides[i] == "A" || provides[i] == "P" || provides[i] == "R" || provides[i] == "Ptent") {
      RCP<MuemexData<Matrix_t> > mydata = Teuchos::rcp_static_cast<MuemexData<Matrix_t> >(mexOutput[i]);
      lvl.Set(provides[i], mydata->getData(), factory);
    } else if (provides[i] == "Nullspace" || provides[i] == "Coordinates") {
      RCP<MuemexData<MultiVector_t> > mydata = Teuchos::rcp_static_cast<MuemexData<MultiVector_t> >(mexOutput[i]);
      lvl.Set(provides[i], mydata->getData(), factory);
    } else if (provides[i] == "Aggregates") {
      RCP<MuemexData<Aggregates_t> > mydata = Teuchos::rcp_static_cast<MuemexData<Aggregates_t> >(mexOutput[i]);
      lvl.Set(provides[i], mydata->getData(), factory);
    } else if (provides[i] == "UnAmalgamationInfo") {
      RCP<MuemexData<AmalgamationInfo_t> > mydata = Teuchos::rcp_static_cast<MuemexData<AmalgamationInfo_t> >(mexOutput[i]);
      lvl.Set(provides[i], mydata->getData(), factory);
    } else if (provides[i] == "Graph") {
      RCP<MuemexData<Graph_t> > mydata = Teuchos::rcp_static_cast<MuemexData<Graph_t> >(mexOutput[i]);
      lvl.Set(provides[i], mydata->getData(), factory);
    } else {
      vector<string> words;
      string badNameMsg = "Custom Muemex variables in \"Provides\" list require a type and a name, e.g. \"double myVal\". \n Leading and trailing spaces are OK. \n Don't know how to handle \"" + provides[i] + "\".\n";
      // compare type without case sensitivity
      char* buf = (char*)malloc(provides[i].size() + 1);
      strcpy(buf, provides[i].c_str());
      for (char* iter = buf; *iter != ' '; iter++) {
        if (*iter == 0) {
          free(buf);
          throw runtime_error(badNameMsg);
        }
        *iter = (char)tolower(*iter);
      }
      const char* wordDelim = " ";
      char* mark            = strtok(buf, wordDelim);
      while (mark != NULL) {
        string wordStr(mark);
        words.push_back(wordStr);
        mark = strtok(NULL, wordDelim);
      }
      if (words.size() != 2) {
        free(buf);
        throw runtime_error(badNameMsg);
      }
      const char* typeStr = words[0].c_str();
      if (strstr(typeStr, "ordinalvector")) {
        typedef RCP<Xpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t> > LOVector_t;
        RCP<MuemexData<LOVector_t> > mydata = Teuchos::rcp_static_cast<MuemexData<LOVector_t> >(mexOutput[i]);
        addLevelVariable<LOVector_t>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "map")) {
        typedef RCP<Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t> > Map_t;
        RCP<MuemexData<Map_t> > mydata = Teuchos::rcp_static_cast<MuemexData<Map_t> >(mexOutput[i]);
        addLevelVariable<Map_t>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "scalar")) {
        RCP<MuemexData<Scalar> > mydata = Teuchos::rcp_static_cast<MuemexData<Scalar> >(mexOutput[i]);
        addLevelVariable<Scalar>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "double")) {
        RCP<MuemexData<double> > mydata = Teuchos::rcp_static_cast<MuemexData<double> >(mexOutput[i]);
        addLevelVariable<double>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "complex")) {
        RCP<MuemexData<complex_t> > mydata = Teuchos::rcp_static_cast<MuemexData<complex_t> >(mexOutput[i]);
        addLevelVariable<complex_t>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "matrix")) {
        RCP<MuemexData<Matrix_t> > mydata = Teuchos::rcp_static_cast<MuemexData<Matrix_t> >(mexOutput[i]);
        addLevelVariable<Matrix_t>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "multivector")) {
        RCP<MuemexData<MultiVector_t> > mydata = Teuchos::rcp_static_cast<MuemexData<MultiVector_t> >(mexOutput[i]);
        addLevelVariable<MultiVector_t>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "int")) {
        RCP<MuemexData<int> > mydata = Teuchos::rcp_static_cast<MuemexData<int> >(mexOutput[i]);
        addLevelVariable<int>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "bool")) {
        RCP<MuemexData<bool> > mydata = Teuchos::rcp_static_cast<MuemexData<bool> >(mexOutput[i]);
        addLevelVariable<bool>(mydata->getData(), words[1], lvl, factory);
      } else if (strstr(typeStr, "string")) {
        RCP<MuemexData<string> > mydata = Teuchos::rcp_static_cast<MuemexData<string> >(mexOutput[i]);
        addLevelVariable<string>(mydata->getData(), words[1], lvl, factory);
      } else {
        free(buf);
        throw std::runtime_error(words[0] + " is not a known variable type.");
      }
      free(buf);
    }
  }
}

// Throwable Stubs for long long

template <>
std::vector<Teuchos::RCP<MuemexArg> > processNeeds<double, mm_LocalOrd, long long, mm_node_t>(const Factory* factory, std::string& needsParam, Level& lvl) {
  throw std::runtime_error("Muemex does not support long long for global indices");
}

template <>
std::vector<Teuchos::RCP<MuemexArg> > processNeeds<complex_t, mm_LocalOrd, long long, mm_node_t>(const Factory* factory, std::string& needsParam, Level& lvl) {
  throw std::runtime_error("Muemex does not support long long for global indices");
}

template <>
void processProvides<double, mm_LocalOrd, long long, mm_node_t>(std::vector<Teuchos::RCP<MuemexArg> >& mexOutput, const Factory* factory, std::string& providesParam, Level& lvl) {
  throw std::runtime_error("Muemex does not support long long for global indices");
}

template <>
void processProvides<complex_t, mm_LocalOrd, long long, mm_node_t>(std::vector<Teuchos::RCP<MuemexArg> >& mexOutput, const Factory* factory, std::string& providesParam, Level& lvl) {
  throw std::runtime_error("Muemex does not support long long for global indices");
}

}  // namespace MueLu
#endif  // HAVE_MUELU_MATLAB error handler
#endif  // MUELU_MATLABUTILS_DEF_HPP guard
