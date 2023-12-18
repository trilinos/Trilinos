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

#ifndef MUELU_MATLABUTILS_DECL_HPP
#define MUELU_MATLABUTILS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA)
#error "Muemex requires MATLAB, Epetra and Tpetra."
#else

#include "mex.h"
#include <string>
#include <complex>
#include <stdexcept>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "MueLu_Factory.hpp"
#include "MueLu_Hierarchy_decl.hpp"
#include "MueLu_Aggregates_decl.hpp"
#include "MueLu_AmalgamationInfo_decl.hpp"
#include "MueLu_Utilities_decl.hpp"
#include "MueLu_Graph_decl.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_VectorFactory.hpp"
#include <Tpetra_Core.hpp>

#include "Kokkos_DynRankView.hpp"

namespace MueLu {

enum MuemexType {
  INT,
  BOOL,
  DOUBLE,
  COMPLEX,
  STRING,
  XPETRA_MAP,
  XPETRA_ORDINAL_VECTOR,
  TPETRA_MULTIVECTOR_DOUBLE,
  TPETRA_MULTIVECTOR_COMPLEX,
  TPETRA_MATRIX_DOUBLE,
  TPETRA_MATRIX_COMPLEX,
  XPETRA_MATRIX_DOUBLE,
  XPETRA_MATRIX_COMPLEX,
  XPETRA_MULTIVECTOR_DOUBLE,
  XPETRA_MULTIVECTOR_COMPLEX,
  EPETRA_CRSMATRIX,
  EPETRA_MULTIVECTOR,
  AGGREGATES,
  AMALGAMATION_INFO,
  GRAPH
#ifdef HAVE_MUELU_INTREPID2
  ,
  FIELDCONTAINER_ORDINAL
#endif
};

typedef Tpetra::KokkosCompat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> mm_node_t;
typedef typename Tpetra::Map<>::local_ordinal_type mm_LocalOrd;  // these are used for LocalOrdinal and GlobalOrdinal of all xpetra/tpetra templated types
typedef typename Tpetra::Map<>::global_ordinal_type mm_GlobalOrd;
typedef std::complex<double> complex_t;
typedef Tpetra::Map<> muemex_map_type;
typedef Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_double;
typedef Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_complex;
typedef Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector_double;
typedef Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector_complex;
typedef Xpetra::Map<mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_map;
typedef Xpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_ordinal_vector;
typedef Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_double;
typedef Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_complex;
typedef Xpetra::CrsGraph<mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_CrsGraph;
typedef Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_double;
typedef Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_complex;
typedef MueLu::Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_double;
typedef MueLu::Hierarchy<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_complex;
typedef MueLu::Aggregates<mm_LocalOrd, mm_GlobalOrd, mm_node_t> MAggregates;
typedef MueLu::AmalgamationInfo<mm_LocalOrd, mm_GlobalOrd, mm_node_t> MAmalInfo;
typedef MueLu::GraphBase<mm_LocalOrd, mm_GlobalOrd, mm_node_t> MGraph;

#ifdef HAVE_MUELU_INTREPID2
typedef Kokkos::DynRankView<mm_LocalOrd, typename mm_node_t::device_type> FieldContainer_ordinal;
#endif

class MuemexArg {
 public:
  MuemexArg(MuemexType dataType) { type = dataType; }
  MuemexType type;
};

template <typename T>
MuemexType getMuemexType(const T& data);

template <typename T>
class MuemexData : public MuemexArg {
 public:
  MuemexData(T& data);                   // Construct from pre-existing data, to pass to MATLAB.
  MuemexData(T& data, MuemexType type);  // Construct from pre-existing data, to pass to MATLAB.
  MuemexData(const mxArray* mxa);        // Construct from MATLAB array, to get from MATLAB.
  mxArray* convertToMatlab();            // Create a MATLAB object and copy this data to it
  T& getData();                          // Set and get methods
  void setData(T& data);

 private:
  T data;
};

template <typename T>
MuemexType getMuemexType(const T& data);

template <typename T>
MuemexType getMuemexType();

template <typename T>
T loadDataFromMatlab(const mxArray* mxa);

template <typename T>
mxArray* saveDataToMatlab(T& data);

// Add data to level. Set the keep flag on the data to "user-provided" so it's not deleted.
template <typename T>
void addLevelVariable(const T& data, std::string& name, Level& lvl, const FactoryBase* fact = NoFactory::get());

template <typename T>
const T& getLevelVariable(std::string& name, Level& lvl);

// Functions used to put data through matlab factories - first arg is "this" pointer of matlab factory
template <typename Scalar = double, typename LocalOrdinal = mm_LocalOrd, typename GlobalOrdinal = mm_GlobalOrd, typename Node = mm_node_t>
std::vector<Teuchos::RCP<MuemexArg> > processNeeds(const Factory* factory, std::string& needsParam, Level& lvl);

template <typename Scalar = double, typename LocalOrdinal = mm_LocalOrd, typename GlobalOrdinal = mm_GlobalOrd, typename Node = mm_node_t>
void processProvides(std::vector<Teuchos::RCP<MuemexArg> >& mexOutput, const Factory* factory, std::string& providesParam, Level& lvl);

// create a sparse array in Matlab
template <typename Scalar>
mxArray* createMatlabSparse(int numRows, int numCols, int nnz);
template <typename Scalar>
mxArray* createMatlabMultiVector(int numRows, int numCols);
template <typename Scalar>
void fillMatlabArray(Scalar* array, const mxArray* mxa, int n);
int* mwIndex_to_int(int N, mwIndex* mwi_array);
bool isValidMatlabAggregates(const mxArray* mxa);
bool isValidMatlabGraph(const mxArray* mxa);
std::vector<std::string> tokenizeList(const std::string& param);
// The two callback functions that MueLu can call to run anything in MATLAB
void callMatlabNoArgs(std::string function);
std::vector<Teuchos::RCP<MuemexArg> > callMatlab(std::string function, int numOutputs, std::vector<Teuchos::RCP<MuemexArg> > args);
Teuchos::RCP<Teuchos::ParameterList> getInputParamList();
Teuchos::RCP<MuemexArg> convertMatlabVar(const mxArray* mxa);

// trim from start
static inline std::string& ltrim(std::string& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline std::string& rtrim(std::string& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline std::string& trim(std::string& s) {
  return ltrim(rtrim(s));
}

}  // namespace MueLu

#endif  // HAVE_MUELU_MATLAB error handler
#endif  // MUELU_MATLABUTILS_DECL_HPP guard
