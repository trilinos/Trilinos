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

#include "mex.h"
#include <string>
#include <complex>
#include <stdexcept>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Factory.hpp"
#include "MueLu_Hierarchy_decl.hpp"
#include "MueLu_Aggregates_decl.hpp"
#include "MueLu_AmalgamationInfo_decl.hpp"
#include "MueLu_Utilities_decl.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_VectorFactory.hpp"
#include <Tpetra_DefaultPlatform.hpp>

#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA) || !defined(HAVE_MUELU_TPETRA)
#error "Muemex requires MATLAB, Epetra and Tpetra."
#else

namespace MueLu
{

enum MUEMEX_TYPE
{
  INT,
  DOUBLE,
  COMPLEX,
  STRING,
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
};

typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> mm_node_t;
typedef int mm_LocalOrd;  //these are used for LocalOrdinal and GlobalOrdinal of all xpetra/tpetra templated types
typedef int mm_GlobalOrd;
typedef std::complex<double> complex_t;
typedef Tpetra::Map<> muemex_map_type;
typedef Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_double;
typedef Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_complex;
typedef Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector_double;
typedef Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector_complex;
typedef Xpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_ordinal_vector;
typedef Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_double;
typedef Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_complex;
typedef Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_double;
typedef Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_complex;
typedef MueLu::Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_double;
typedef MueLu::Hierarchy<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_complex;
typedef MueLu::Aggregates<mm_LocalOrd, mm_GlobalOrd, mm_node_t> MAggregates;
typedef MueLu::AmalgamationInfo<mm_LocalOrd, mm_GlobalOrd, mm_node_t> MAmalInfo;
typedef MueLu::Graph<mm_LocalOrd, mm_GlobalOrd, mm_node_t> MGraph;

class MuemexArg
{
  public:
    MuemexArg(MUEMEX_TYPE dataType) {type = dataType;}
    MUEMEX_TYPE type;
};

template<typename T> 
MUEMEX_TYPE getMuemexType(const T & data);

template<typename T>
class MuemexData : public MuemexArg
{
  public:
    MuemexData(T& data); //Construct from pre-existing data, to pass to MATLAB.
    MuemexData(T& data, MUEMEX_TYPE type);        //Construct from pre-existing data, to pass to MATLAB.
    MuemexData(const mxArray* mxa); //Construct from MATLAB array, to get from MATLAB.
    mxArray* convertToMatlab(); //Create a MATLAB object and copy this data to it
    T& getData();                         //Set and get methods
    void setData(T& data);
  private:
    T data;
};

template<typename T> 
MUEMEX_TYPE getMuemexType(const T & data);

template<typename T>
MUEMEX_TYPE getMuemexType();

template<typename T>
T loadDataFromMatlab(const mxArray* mxa);

template<typename T>
mxArray* saveDataToMatlab(T& data);

//Add data to level. Set the keep flag on the data to "user-provided" so it's not deleted.
template<typename T>
void addLevelVariable(const T& data, std::string& name, Level& lvl);

template<typename T>
const T& getLevelVariable(std::string& name, Level& lvl);

//Functions used to put data through matlab factories - first arg is "this" pointer of matlab factory
template<typename Scalar = double, typename LocalOrdinal = mm_LocalOrd, typename GlobalOrdinal = mm_GlobalOrd, typename Node = mm_node_t>
std::vector<Teuchos::RCP<MuemexArg>> processNeeds(const Factory* factory, std::string& needsParam, Level& lvl);

template<typename Scalar = double, typename LocalOrdinal = mm_LocalOrd, typename GlobalOrdinal = mm_GlobalOrd, typename Node = mm_node_t>
void processProvides(std::vector<Teuchos::RCP<MuemexArg>>& mexOutput, const Factory* factory, std::string& providesParam, Level& lvl);

//create a sparse array in Matlab
template<typename Scalar> mxArray* createMatlabSparse(int numRows, int numCols, int nnz);
template<typename Scalar> mxArray* createMatlabMultiVector(int numRows, int numCols);
template<typename Scalar> void fillMatlabArray(Scalar* array, const mxArray* mxa, int n);
int* mwIndex_to_int(int N, mwIndex* mwi_array);
bool isValidMatlabAggregates(const mxArray* mxa);
std::vector<std::string> tokenizeList(const std::string& param);
//The two callback functions that MueLu can call to run anything in MATLAB
void callMatlabNoArgs(std::string function);
std::vector<Teuchos::RCP<MuemexArg>> callMatlab(std::string function, int numOutputs, std::vector<Teuchos::RCP<MuemexArg>> args);
Teuchos::RCP<Teuchos::ParameterList> getInputParamList();
Teuchos::RCP<MuemexArg> convertMatlabVar(const mxArray* mxa);

}//end namespace

#endif //HAVE_MUELU_MATLAB error handler
#endif //MUELU_MATLABUTILS_DECL_HPP guard
