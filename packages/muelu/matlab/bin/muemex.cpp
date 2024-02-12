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

#include "muemex.h"

#define IS_FALSE 0
#define IS_TRUE 1
#define MUEMEX_ERROR -1

// Do not compile MueMex if any of these aren't available
#if !defined HAVE_MUELU_EPETRA || !defined HAVE_MUELU_MATLAB
#error "MueMex requires Epetra, Tpetra and MATLAB."
#endif

#include <Tpetra_Core.hpp>
#include "MueLu_MatlabUtils.hpp"
#include "MueLu_TwoLevelMatlabFactory.hpp"
#include "MueLu_SingleLevelMatlabFactory.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"

using namespace std;
using namespace Teuchos;

extern void _main();

/* MUEMEX Teuchos Parameters*/
#define MUEMEX_INTERFACE "Linear Algebra"

/* Default values */
#define MUEMEX_DEFAULT_LEVELS 10
#define MUEMEX_DEFAULT_NUMPDES 1
#define MUEMEX_DEFAULT_ADAPTIVEVECS 0
#define MUEMEX_DEFAULT_USEDEFAULTNS true
#define MMABS(x) ((x) > 0 ? (x) : (-(x)))
#define MMISINT(x) ((x) == 0 ? (((x - (int)(x)) < 1e-15) ? true : false) : (((x - (int)(x)) < 1e-15 * MMABS(x)) ? true : false))

/* Debugging */
//#define VERBOSE_OUTPUT

namespace MueLu {

// Need subclass of Hierarchy that gives public access to list of FactoryManagers
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class OpenHierarchy : public Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  const RCP<const FactoryManagerBase>& GetFactoryManager(const int levelID) const;
};

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
const RCP<const FactoryManagerBase>& OpenHierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetFactoryManager(const int levelID) const {
  TEUCHOS_TEST_FOR_EXCEPTION(levelID < 0 || levelID > this->GetNumLevels(), Exceptions::RuntimeError,
                             "MueLu::Hierarchy::GetFactoryManager(): invalid input parameter value: LevelID = " << levelID);
  return this->GetLevelManager(levelID);
}

// Declare and call default constructor for data_pack_list vector (starts empty)
vector<RCP<MuemexSystem>> MuemexSystemList::list;
int MuemexSystemList::nextID = 0;

// Need a global flag to keep track of Epetra vs. Tpetra for constructing multivectors for param lists
bool useEpetra = false;

// Parse a string to get each bit of Belos verbosity
int strToMsgType(const char* str) {
  if (str == NULL)
    return Belos::Errors;
  else if (strstr(str, "Warnings") != NULL)
    return Belos::Warnings;
  else if (strstr(str, "IterationDetails") != NULL)
    return Belos::IterationDetails;
  else if (strstr(str, "OrthoDetails") != NULL)
    return Belos::OrthoDetails;
  else if (strstr(str, "FinalSummary") != NULL)
    return Belos::FinalSummary;
  else if (strstr(str, "TimingDetails") != NULL)
    return Belos::TimingDetails;
  else if (strstr(str, "StatusTestDetails") != NULL)
    return Belos::StatusTestDetails;
  else if (strstr(str, "Debug") != NULL)
    return Belos::Debug;
  // This has no effect when added/OR'd with flags
  return Belos::Errors;
}

MuemexType strToDataType(const char* str, char* typeName, bool complexFlag = false) {
  std::string temp(str);
  std::string myStr = trim(temp);
  MuemexType matrixType, multivectorType, scalarType;
  if (!complexFlag) {
    matrixType      = XPETRA_MATRIX_DOUBLE;
    multivectorType = XPETRA_MULTIVECTOR_DOUBLE;
    scalarType      = DOUBLE;
  } else {
    matrixType      = XPETRA_MATRIX_COMPLEX;
    multivectorType = XPETRA_MULTIVECTOR_COMPLEX;
    scalarType      = COMPLEX;
  }
  size_t npos = string::npos;
  if (myStr == "A" ||
      myStr == "P" || myStr == "Ptent" ||
      myStr == "R")
    return matrixType;
  if (myStr == "Nullspace")
    return multivectorType;
  if (myStr == "Aggregates")
    return AGGREGATES;
  if (myStr == "Graph")
    return GRAPH;
  if (myStr == "Coordinates")
    return XPETRA_MULTIVECTOR_DOUBLE;
#ifdef HAVE_MUELU_INTREPID2
  if (myStr == "pcoarsen: element to node map")
    return FIELDCONTAINER_ORDINAL;
#endif
  // Check for custom variable
  size_t firstWordStart = myStr.find_first_not_of(' ');
  size_t firstWordEnd   = myStr.find(' ', firstWordStart);
  std::string firstWord = myStr.substr(firstWordStart, firstWordEnd - firstWordStart);
  if (firstWord.length() > 0) {
    temp                   = myStr.substr(firstWordEnd, myStr.length() - firstWordEnd);
    std::string secondWord = trim(temp);
    if (secondWord.length() > 0) {
      // make first word lowercase
      std::transform(firstWord.begin(), firstWord.end(), firstWord.begin(), ::tolower);
      // compare first word with possible values
      if (firstWord.find("matrix") != npos)
        return matrixType;
      if (firstWord.find("multivector") != npos)
        return multivectorType;
      if (firstWord.find("map") != npos)
        return XPETRA_MAP;
      if (firstWord.find("ordinalvector") != npos)
        return XPETRA_ORDINAL_VECTOR;
      if (firstWord.find("int") != npos)
        return INT;
      if (firstWord.find("scalar") != npos)
        return scalarType;
      if (firstWord.find("double") != npos)
        return DOUBLE;
      if (firstWord.find("complex") != npos)
        return COMPLEX;
    }
  }
  if (typeName) {
    std::string typeString(typeName);
    std::transform(typeString.begin(), typeString.end(), typeString.begin(), ::tolower);
    if (typeString.find("matrix") != npos)
      return matrixType;
    if (typeString.find("multivector") != npos)
      return multivectorType;
    if (typeString.find("map") != npos)
      return XPETRA_MAP;
    if (typeString.find("ordinalvector") != npos)
      return XPETRA_ORDINAL_VECTOR;
    if (typeString.find("int") != npos)
      return INT;
    if (typeString.find("scalar") != npos)
      return scalarType;
    if (typeString.find("double") != npos)
      return DOUBLE;
    if (typeString.find("complex") != npos)
      return COMPLEX;
    string errMsg = typeString + " is not a valid type.";
    throw runtime_error(errMsg);
  }
  throw runtime_error("Could not determine type of data.");
}

// Parse a string to get Belos output style (Brief is default)
int strToOutputStyle(const char* str) {
  if (strstr("General", str) != NULL)
    return Belos::General;
  else
    return Belos::Brief;
}

// Get Belos "Verbosity" setting for its ParameterList
int getBelosVerbosity(const char* input) {
  int result = 0;
  char* str  = (char*)input;
  char* pch;
  pch = strtok(str, " +,");
  if (pch == NULL)
    return result;
  result |= strToMsgType(pch);
  while (pch != NULL) {
    pch = strtok(NULL, " +,");
    if (pch == NULL)
      return result;
    result |= strToMsgType(pch);
  }
  return result;
}

template <typename Scalar>
mxArray* TpetraSystem<Scalar>::solve(RCP<ParameterList> params, RCP<Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> matrix, const mxArray* b, int& iters) {
  mxArray* output;
  try {
    int matSize = A->getGlobalNumRows();
    // Define Tpetra vector/multivector types for convenience
    typedef Tpetra::Vector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Vector;
    typedef Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector;
    typedef Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_Operator;
    RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
    // numGlobalIndices for map constructor is the number of rows in matrix/vectors, right?
    RCP<const muemex_map_type> map = rcp(new muemex_map_type(matSize, (mm_GlobalOrd)0, comm));
    RCP<Tpetra_MultiVector> rhs    = loadDataFromMatlab<RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(b);
    RCP<Tpetra_MultiVector> lhs    = rcp(new Tpetra_MultiVector(map, rhs->getNumVectors()));
    // rhs is initialized, lhs is not
    //  Default params
    params->get("Output Frequency", 1);
    params->get("Output Style", Belos::Brief);

#ifdef VERBOSE_OUTPUT
    params->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
#else
    params->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::IterationDetails | Belos::Warnings | Belos::StatusTestDetails);
#endif
    // register all possible solvers
    auto problem = rcp(new Belos::LinearProblem<Scalar, Tpetra_MultiVector, Tpetra_Operator>(matrix, lhs, rhs));
    problem->setRightPrec(prec);
    if (!problem->setProblem()) {
      throw std::runtime_error("ERROR: failed to set up Belos problem.");
    }
    std::string solverName = "GMRES";
    if (params->isParameter("solver")) {
      solverName = params->template get<std::string>("solver");
    }
    // Convert from basic MueMex solver names to the official Belos names.
    // At the same time, check that solverName is in the valid set.
    std::string belosSolverName;
    if (solverName == "GMRES") {
      belosSolverName = "PseudoBlock GMRES";
    } else if (solverName == "CG") {
      belosSolverName = "PseudoBlock CG";
    } else {
      std::string msg = std::string("ERROR: requested solver \"") + solverName + "\" not supported. Currently supported solvers: CG, GMRES";
      mexPrintf("%s\n", msg.c_str());
      output = mxCreateDoubleScalar(0);
      return output;
    }
    Teuchos::RCP<Belos::SolverManager<Scalar, Tpetra_MultiVector, Tpetra_Operator>> solver;
    Belos::SolverFactory<Scalar, Tpetra_MultiVector, Tpetra_Operator> factory;
    try {
      // Just use the default parameters for the solver
      solver = factory.create(belosSolverName, params);
    } catch (std::exception& e) {
      mexPrintf("%s\n", e.what());
      output = mxCreateDoubleScalar(0);
      return output;
    }
    solver->setProblem(problem);
    Belos::ReturnType ret = solver->solve();
    iters                 = solver->getNumIters();
    if (ret == Belos::Converged) {
      mexPrintf("Success, Belos converged!\n");
      output = saveDataToMatlab(lhs);
    } else {
      mexPrintf("Belos failed to converge.\n");
      iters  = 0;
      output = mxCreateDoubleScalar(0);
    }
  } catch (exception& e) {
    mexPrintf("Error occurred while running Belos solver:\n");
    cout << e.what() << endl;
    output = mxCreateDoubleScalar(0);
  }
  return output;
}

template <typename Scalar>
mxArray* TpetraSystem<Scalar>::apply(const mxArray* r) {
  typedef Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_MultiVector;
  RCP<Tpetra_MultiVector> rhs = loadDataFromMatlab<RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(r);
  RCP<Tpetra_MultiVector> lhs = rcp(new Tpetra_MultiVector(rhs->getMap(), rhs->getNumVectors()));
  try {
    this->prec->apply(*rhs, *lhs);
  } catch (exception& e) {
    mexPrintf("Error occurred while applying MueLu-Tpetra preconditioner:\n");
    cout << e.what() << endl;
  }
  return saveDataToMatlab(lhs);
}

template <>
RCP<Hierarchy_double> getDatapackHierarchy<double>(MuemexSystem* dp) {
  RCP<MueLu::Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier;
  switch (dp->type) {
    case EPETRA: {
      EpetraSystem* pack = (EpetraSystem*)dp;
      hier               = pack->getHierarchy();
      break;
    }
    case TPETRA: {
      TpetraSystem<double>* pack = (TpetraSystem<double>*)dp;
      hier                       = pack->getHierarchy();
      break;
    }
    default: {
      throw runtime_error("Got unexpected linear system type for real-valued functions.");
    }
  }
  return hier;
}

#ifdef HAVE_COMPLEX_SCALARS
template <>
RCP<Hierarchy_complex> getDatapackHierarchy<complex_t>(MuemexSystem* dp) {
  return ((TpetraSystem<complex_t>*)dp)->getHierarchy();
}
#endif

template <typename Scalar, typename T>
void setHierarchyData(MuemexSystem* problem, int levelID, T& data, string& dataName) {
  RCP<Level> level;
  if (problem->type == EPETRA) {
    RCP<Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier = ((EpetraSystem*)problem)->getHierarchy();
    level                                                             = hier->GetLevel(levelID);
  } else if (problem->type == TPETRA) {
    RCP<Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier = ((TpetraSystem<double>*)problem)->getHierarchy();
    level                                                             = hier->GetLevel(levelID);
  } else if (problem->type == TPETRA_COMPLEX) {
#ifdef HAVE_COMPLEX_SCALARS
    RCP<Hierarchy<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier = ((TpetraSystem<complex_t>*)problem)->getHierarchy();
    level                                                                = hier->GetLevel(levelID);
#else
    throw std::runtime_error("setHierarchyData(): complex scalars not supported.");
#endif
  }
  if (level.is_null())
    throw runtime_error("Error getting level when setting custom level data.");
  level->Set(dataName, data);
  level->AddKeepFlag(dataName, NoFactory::get(), UserData);
}

// data pack base class implementation

MuemexSystem::MuemexSystem(DataPackType probType)
  : id(MUEMEX_ERROR)
  , type(probType) {}
MuemexSystem::~MuemexSystem() {}

mxArray* MuemexSystem::getHierarchyData(string dataName, MuemexType dataType, int levelID) {
  mxArray* output = NULL;
  try {
    // First, get Level, which doesn't depend on Epetra vs. Tpetra
    RCP<MueLu::Level> level;
    RCP<const FactoryManagerBase> fmb;
    if (this->type == TPETRA) {
      TpetraSystem<double>* tsys = (TpetraSystem<double>*)this;
      if (tsys->keepAll)
        fmb = tsys->systemManagers[levelID];
    } else if (this->type == TPETRA_COMPLEX) {
#ifdef HAVE_COMPLEX_SCALARS
      TpetraSystem<complex_t>* tsys = (TpetraSystem<complex_t>*)this;
      if (tsys->keepAll)
        fmb = tsys->systemManagers[levelID];
#else
      throw std::runtime_error("getHierarchyData(): complex scalars not supported.");
#endif
    }
    const FactoryBase* factory = NoFactory::get();  //(ptr to constant)
    bool needFMB               = true;
    if (dataName == "A" || dataName == "P")  // these are kept by default, don't use actual factory pointer
      // Otherwise would break getting A and P when 'keep' is off
      needFMB = false;
    switch (this->type) {
      case EPETRA:
      case TPETRA: {
        RCP<OpenHierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier = rcp_static_cast<OpenHierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(getDatapackHierarchy<double>(this));
        level                                                                 = hier->GetLevel(levelID);
        if (needFMB) {
          if (fmb.is_null())
            fmb = (RCP<const FactoryManagerBase>)hier->GetFactoryManager(levelID);
          if (!fmb.is_null()) {
            try {
              factory = fmb->GetFactory(dataName).get();
            } catch (exception& e) {
            }  // forced to try using NoFactory (which will work with default keeps A, P)
          }
        }
        break;
      }
      case TPETRA_COMPLEX: {
#ifdef HAVE_COMPLEX_SCALARS
        RCP<OpenHierarchy<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier = rcp_static_cast<OpenHierarchy<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(getDatapackHierarchy<complex_t>(this));
        level                                                                    = hier->GetLevel(levelID);
        if (needFMB) {
          if (fmb.is_null())
            fmb = (RCP<const FactoryManagerBase>)hier->GetFactoryManager(levelID);
          if (!fmb.is_null()) {
            try {
              factory = fmb->GetFactory(dataName).get();
            } catch (exception& e) {
            }  // attempt to use NoFactory
          }
        }
        break;
#else
        throw std::runtime_error("Complex scalars not supported");
#endif
      }
    }
    if (level.is_null())
      throw runtime_error("Can't get level data because level is null.");
    bool dataIsAvailable = level->IsAvailable(dataName, factory);
    if (!dataIsAvailable) {
      // Give the level the FactoryManager again so it can provide the data (by re-creating it, if necessary)
      level->SetFactoryManager(fmb);
    }
    // Given the dataName and factory pointer, all data in the level should now be accessible
    switch (dataType) {
      case XPETRA_MATRIX_DOUBLE:
        return saveDataToMatlab(level->Get<RCP<Xpetra_Matrix_double>>(dataName, factory));
      case XPETRA_MATRIX_COMPLEX: {
#ifdef HAVE_COMPLEX_SCALARS
        return saveDataToMatlab(level->Get<RCP<Xpetra_Matrix_complex>>(dataName, factory));
#endif
      }
      case XPETRA_MULTIVECTOR_DOUBLE:
        if (dataName == "Coordinates") {
          // Coordinates is special because it's always user-provided on level 0, not always provided at all, not always kept in the level (only kept if doing agg viz, etc), and is always MV<double> regardless of problem scalar type
          double errReturn = -1;
          if (level->GetLevelID() == 0) {
            // Try to get coordinates as if it's user data, but don't be surprised if it's not there at all.
            try {
              RCP<Xpetra_MultiVector_double> coords = level->Get<RCP<Xpetra_MultiVector_double>>(dataName, NoFactory::get());
              if (coords.is_null())
                throw runtime_error("Coordinates were not available (Level 0).");  // just print the message below and return -1
              return saveDataToMatlab(coords);
            } catch (exception& e) {
              cout << endl
                   << "Coordinates were not available on Level 0." << endl;
              cout << "They must be provided by the user and aren't generated or kept by default (even in MueMex 'keep' mode)." << endl;
              cout << "User-provided coordinates for Level 0 will be kept and passed to other levels if something requires them:" << endl;
              cout << "aggregate visualization, brick aggregation, repartitioning or distance laplacian filtering." << endl
                   << endl;
              return saveDataToMatlab(errReturn);
            }
          } else {
            // If coords are provided & kept, they are produced by CoordinatesTransferFactory in levels > 0.
            try {
              RCP<Xpetra_MultiVector_double> coords = level->Get<RCP<Xpetra_MultiVector_double>>(dataName, factory);
              if (coords.is_null())
                throw runtime_error("Coordinates were not available (Level > 0).");
              return saveDataToMatlab(coords);
            } catch (exception& e) {
              cout << "Coordinates must be provided by the user and aren't generated or kept by default (even in MueMex 'keep' mode)." << endl;
              cout << "User-provided coordinates for Level 0 will be kept and passed to other levels if something requires them:" << endl;
              cout << "aggregate visualization, brick aggregation, repartitioning or distance laplacian filtering." << endl
                   << endl;
              return saveDataToMatlab(errReturn);
            }
          }
        } else {
          return saveDataToMatlab(level->Get<RCP<Xpetra_MultiVector_double>>(dataName, factory));
        }
      case XPETRA_MULTIVECTOR_COMPLEX:
        return saveDataToMatlab(level->Get<RCP<Xpetra_MultiVector_complex>>(dataName, factory));
      case XPETRA_MAP:
        return saveDataToMatlab(level->Get<RCP<Xpetra_map>>(dataName, factory));
      case XPETRA_ORDINAL_VECTOR:
        return saveDataToMatlab(level->Get<RCP<Xpetra_ordinal_vector>>(dataName, factory));
      case DOUBLE:
        return saveDataToMatlab(level->Get<double>(dataName, factory));
      case COMPLEX:
        return saveDataToMatlab(level->Get<complex_t>(dataName, factory));
      case INT:
        return saveDataToMatlab(level->Get<int>(dataName, factory));
      case AGGREGATES:
        return saveDataToMatlab(level->Get<RCP<MAggregates>>(dataName, factory));
      case GRAPH:
        return saveDataToMatlab(level->Get<RCP<MGraph>>(dataName, factory));
#ifdef HAVE_MUELU_INTREPID2
        return saveDataToMatlab(level->Get<RCP<FieldContainer_ordinal>>(dataName, factory));
#endif
      default:
        throw runtime_error("Invalid MuemexType for getting hierarchy data.");
    }
  } catch (exception& e) {
    mexPrintf("Error occurred while getting hierarchy data.\n");
    cout << e.what() << endl;
  }
  return output;
}

// EpetraSystem impl

EpetraSystem::EpetraSystem()
  : MuemexSystem(EPETRA) {}
EpetraSystem::~EpetraSystem() {}

int EpetraSystem::status() {
  mexPrintf("**** Problem ID %d [MueLu_Epetra] ****\n", id);
  if (!A.is_null())
    mexPrintf("Matrix: %dx%d w/ %d nnz\n", A->NumGlobalRows(), A->NumGlobalCols(), A->NumMyNonzeros());
  mexPrintf("Operator Complexity: %f\n", operatorComplexity);
  if (!List.is_null()) {
    mexPrintf("Parameter List:\n");
    List->print();
  }
  mexPrintf("\n");
  return IS_TRUE;
} /*end status*/

int EpetraSystem::setup(const mxArray* matlabA, bool haveCoords, const mxArray* matlabCoords) {
  bool success = false;
  try {
    /* Matrix Fill */
    A = loadDataFromMatlab<RCP<Epetra_CrsMatrix>>(matlabA);
    if (haveCoords) {
      // Create 'user data' sublist if it doesn't already exist
      auto userData = Teuchos::sublist(List, "user data");
      userData->set("Coordinates", loadDataFromMatlab<RCP<Epetra_MultiVector>>(matlabCoords));
    }
    prec = MueLu::CreateEpetraPreconditioner(A, *List);
    // underlying the Epetra_Operator prec is a MueLu::EpetraOperator
    RCP<MueLu::EpetraOperator> meo = rcp_static_cast<MueLu::EpetraOperator, Epetra_Operator>(prec);
    operatorComplexity             = meo->GetHierarchy()->GetOperatorComplexity();
    success                        = true;
  } catch (exception& e) {
    mexPrintf("Error occurred while setting up epetra problem:\n");
    cout << e.what() << endl;
  }
  return success ? IS_TRUE : IS_FALSE;
} /*end setup*/

/* EpetraSystem::solve - Given two Teuchos lists, one in the EpetraSystem, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   TPL     - Teuchos list of solve-time options [I]
   A       - The matrix to solve with (may not be the one the preconditioned was used for)
   b       - RHS vector [I]
   x       - solution vector [O]
   iters   - number of iterations taken [O]
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
mxArray* EpetraSystem::solve(RCP<ParameterList> TPL, RCP<Epetra_CrsMatrix> matrix, const mxArray* b, int& iters) {
  mxArray* output;
  try {
    // Set up X and B
    Epetra_Map map              = matrix->DomainMap();
    RCP<Epetra_MultiVector> rhs = loadDataFromMatlab<RCP<Epetra_MultiVector>>(b);
    RCP<Epetra_MultiVector> lhs = rcp(new Epetra_MultiVector(map, rhs->NumVectors(), true));
    // Default params
    TPL->get("Output Frequency", 1);
    TPL->get("Output Style", Belos::Brief);
#ifdef VERBOSE_OUTPUT
    TPL->get("Verbosity", Belos::Errors | Belos::Warnings | Belos::Debug | Belos::FinalSummary | Belos::IterationDetails | Belos::OrthoDetails | Belos::TimingDetails | Belos::StatusTestDetails);
#else
    TPL->get("Verbosity", Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::Warnings + Belos::StatusTestDetails);
#endif
    RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>> problem = rcp(new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(matrix, lhs, rhs));
    RCP<Belos::EpetraPrecOp> epo                                                   = rcp(new Belos::EpetraPrecOp(prec));
    problem->setRightPrec(epo);
    bool set = problem->setProblem();
    TEUCHOS_TEST_FOR_EXCEPTION(!set, runtime_error, "Linear Problem failed to set up correctly!");
    Belos::SolverFactory<double, Epetra_MultiVector, Epetra_Operator> factory;
    // Get the solver name from the parameter list, default to PseudoBlockGmres if none specified by user
    string solverName                                                             = TPL->get("solver", "GMRES");
    RCP<Belos::SolverManager<double, Epetra_MultiVector, Epetra_Operator>> solver = factory.create(solverName, TPL);
    solver->setProblem(problem);
    Belos::ReturnType ret = solver->solve();
    if (ret == Belos::Converged) {
      mexPrintf("Success, Belos converged!\n");
      iters  = solver->getNumIters();
      output = saveDataToMatlab(lhs);
    } else {
      mexPrintf("Belos failed to converge.\n");
      iters  = 0;
      output = mxCreateDoubleScalar(0);
    }
    output = saveDataToMatlab(lhs);
  } catch (exception& e) {
    mexPrintf("Error occurred during Belos solve:\n");
    cout << e.what() << endl;
    output = mxCreateDoubleScalar(0);
  }
  return output;
}

mxArray* EpetraSystem::apply(const mxArray* r) {
  RCP<Epetra_MultiVector> rhs = loadDataFromMatlab<RCP<Epetra_MultiVector>>(r);
  Epetra_SerialComm Comm;
  Epetra_Map map(rhs->GlobalLength(), 0, Comm);
  RCP<Epetra_MultiVector> lhs = rcp(new Epetra_MultiVector(map, rhs->NumVectors(), true));
  try {
    this->prec->Apply(*rhs, *lhs);
  } catch (exception& e) {
    mexPrintf("Error occurred while applying MueLu-Epetra preconditioner:\n");
    cout << e.what() << endl;
  }
  return saveDataToMatlab(lhs);
}

RCP<Hierarchy_double> EpetraSystem::getHierarchy() {
  RCP<MueLu::EpetraOperator> meo                                           = rcp_static_cast<MueLu::EpetraOperator, Epetra_Operator>(prec);
  RCP<MueLu::Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> hier = meo->GetHierarchy();
  if (hier.is_null())
    throw runtime_error("Hierarchy from Epetra problem was null.");
  return hier;
}

// tpetra_double_data_pack implementation

template <>
TpetraSystem<double>::TpetraSystem()
  : MuemexSystem(TPETRA) {}
template <>
TpetraSystem<double>::~TpetraSystem() {}

template <typename Scalar>
int TpetraSystem<Scalar>::setup(const mxArray* matlabA, bool haveCoords, const mxArray* matlabCoords) {
  // decide whether do do default or custom setup
  bool doCustomSetup = List->isParameter("keep") && List->isType<bool>("keep") && List->get<bool>("keep");
  List->remove("keep", false);  //"keep" would cause Plist validation to fail if left in
  if (doCustomSetup) {
    try {
      customSetup(matlabA, haveCoords, matlabCoords);
    } catch (exception& e) {
      cout << "An error occurred during Tpetra custom problem setup:" << endl;
      cout << e.what() << endl;
      return IS_FALSE;
    }
  } else {
    try {
      normalSetup(matlabA, haveCoords, matlabCoords);
    } catch (exception& e) {
      cout << "An error occurred during Tpetra preconditioner setup:" << endl;
      cout << e.what();
      return IS_FALSE;
    }
  }
  return IS_TRUE;
}

template <typename Scalar>
void TpetraSystem<Scalar>::normalSetup(const mxArray* matlabA, bool haveCoords, const mxArray* matlabCoords) {
  keepAll = false;
  A       = loadDataFromMatlab<RCP<Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(matlabA);
  RCP<Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> opA(A);
  RCP<MueLu::TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mop;
  if (haveCoords) {
    auto userData = Teuchos::sublist(List, "user data");
    userData->set("Coordinates", loadDataFromMatlab<RCP<Tpetra_MultiVector_double>>(matlabCoords));
  }
  // Create the nullspace if not already set by user through XML
  if (!(List->isSublist("level 0") && List->sublist("level 0", true).isParameter("Nullspace")) && !(List->isSublist("user data") && List->sublist("user data", true).isParameter("Nullspace"))) {
    int nPDE = MasterList::getDefault<int>("number of equations");
    if (List->isSublist("Matrix")) {
      // Factory style parameter list
      const Teuchos::ParameterList& operatorList = List->sublist("Matrix");
      if (operatorList.isParameter("PDE equations"))
        nPDE = operatorList.get<int>("PDE equations");
    } else if (List->isParameter("number of equations")) {
      // Easy style parameter list
      nPDE = List->get<int>("number of equations");
    }
    mexPrintf("** Constructing nullspace for %d PDEs\n", nPDE);
    auto domainMap = A->getDomainMap();
    auto nullspace = rcp(new Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(domainMap, nPDE));
    if (nPDE == 1) {
      nullspace->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    } else {
      typedef typename Teuchos::ArrayRCP<Scalar>::size_type arrayRCPSizeType;
      for (int i = 0; i < nPDE; i++) {
        Teuchos::ArrayRCP<Scalar> nsData = nullspace->getDataNonConst(i);
        for (arrayRCPSizeType j = 0; j < nsData.size(); j++) {
          mm_GlobalOrd GID = domainMap->getGlobalElement(j) - domainMap->getIndexBase();
          if ((GID - i) % nPDE == 0)
            nsData[j] = Teuchos::ScalarTraits<Scalar>::one();
        }
      }
    }
    auto userData = Teuchos::sublist(List, "user data");
    userData->set("Nullspace", nullspace);
  }
  mop  = MueLu::CreateTpetraPreconditioner<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(opA, *List);
  prec = rcp_implicit_cast<Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(mop);

  // print data??
  // mop->GetHierarchy()->GetLevel(0)->print(std::cout, MueLu::Debug);

  operatorComplexity = mop->GetHierarchy()->GetOperatorComplexity();
}

template <typename Scalar>
void TpetraSystem<Scalar>::customSetup(const mxArray* matlabA, bool haveCoords, const mxArray* matlabCoords) {
  keepAll = true;
  A       = loadDataFromMatlab<RCP<Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(matlabA);
  RCP<MueLu::TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mop;
  // Now modify CreateTpetraPreconditioner to set keep flags on all factories
  typedef Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> MultiVector;
  typedef Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Matrix;
  typedef Hierarchy<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy;
  typedef HierarchyManager<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> HierarchyManager;
  RCP<HierarchyManager> mueluFactory = rcp(new ParameterListInterpreter<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(*List, A->getComm()));
  RCP<Hierarchy> H                   = mueluFactory->CreateHierarchy();
  H->setlib(Xpetra::UseTpetra);
  RCP<Matrix> xA = TpetraCrs_To_XpetraMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(A);
  H->GetLevel(0)->Set("A", xA);
  if (haveCoords) {
    RCP<Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> coords = loadDataFromMatlab<RCP<Xpetra_MultiVector_double>>(matlabCoords);
    H->GetLevel(0)->Set("Coordinates", coords);
  }
  // Decide whether user passed level 0 Nullspace in parameter list. If not, make it here.
  if (!List->isSublist("level 0") || !List->sublist("level 0", true).isParameter("Nullspace")) {
    int nPDE = MasterList::getDefault<int>("number of equations");
    if (List->isSublist("Matrix")) {
      // Factory style parameter list
      const Teuchos::ParameterList& operatorList = List->sublist("Matrix");
      if (operatorList.isParameter("PDE equations"))
        nPDE = operatorList.get<int>("PDE equations");
    } else if (List->isParameter("number of equations")) {
      // Easy style parameter list
      nPDE = List->get<int>("number of equations");
    }
    RCP<MultiVector> nullspace = Xpetra::MultiVectorFactory<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>::Build(xA->getDomainMap(), nPDE);
    if (nPDE == 1) {
      nullspace->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    } else {
      typedef typename Teuchos::ArrayRCP<Scalar>::size_type arrayRCPSizeType;
      for (int i = 0; i < nPDE; i++) {
        Teuchos::ArrayRCP<Scalar> nsData = nullspace->getDataNonConst(i);
        for (arrayRCPSizeType j = 0; j < nsData.size(); j++) {
          // TODO optimizations:
          // TODO This can be optimized by getting the domain map and index base outside the loop.
          // TODO Also, the whole local-to-global lookup table can be fetched one time, instead of repeatedly
          // TODO calling getGlobalElement.
          mm_GlobalOrd GID = A->getDomainMap()->getGlobalElement(j) - A->getDomainMap()->getIndexBase();
          if ((GID - i) % nPDE == 0)
            nsData[j] = Teuchos::ScalarTraits<Scalar>::one();
        }
      }
    }
    H->GetLevel(0)->Set("Nullspace", nullspace);
  }
  Teuchos::ParameterList nonSerialList, dummyList;
  ExtractNonSerializableData(*List, dummyList, nonSerialList);
  HierarchyUtils<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>::AddNonSerializableDataToHierarchy(*mueluFactory, *H, nonSerialList);
  // Set up dummy levels in hierarchy
  for (int i = 0; i < 5; i++) {
    RCP<Level> l = rcp(new Level());
    H->AddLevel(l);
  }
  // Set keep flags on ALL factories in ALL levels
  // We have access to H's list of FactoryManagers so we know how to get factory pointer given the name of the factory.
  // We don't know which names are in the FactoryManagers though so a brute force approach is needed...
  vector<string> keepItems                                               = {"A", "P", "R", "Ptent", "Aggregates", "Coordinates", "UnAmalgamationInfo", "Smoother", "PreSmoother", "PostSmoother", "CoarseSolver", "Graph", "CoarseMap", "Nullspace", "Ppattern", "Constraint", "CoarseNumZLayers", "LineDetection_Layers", "LineDetection_VertLineIds", "Partition", "Importer", "DofsPerNode",
                                                                            "Filtering"
                                                                                                                          "pcoarsen: element to node map"};
  RCP<OpenHierarchy<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> openH = rcp_static_cast<OpenHierarchy<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>, Hierarchy>(H);
  if (openH.is_null())
    throw runtime_error("Could not cast RCP<Hierarchy> to subclass.");
  for (int lvl = 0; lvl < H->GetNumLevels(); lvl++) {
    RCP<const FactoryManagerBase> fman = (RCP<const FactoryManagerBase>)mueluFactory->GetFactoryManager(lvl);
    systemManagers.push_back(fman);
    for (auto s : keepItems) {
      try {
        const RCP<const FactoryBase> fact = fman->GetFactory(s);  // will throw if factory doesn't exist, ignore in that case
        if (!fact.is_null()) {
          FactoryBase* factPtr = (FactoryBase*)fact.get();
          // Add keep flag to level
          H->GetLevel(lvl)->Keep(s, factPtr);
        }
      } catch (exception& e) {
      }
    }
  }
  mueluFactory->SetupHierarchy(*H);
  operatorComplexity = H->GetOperatorComplexity();
  prec               = rcp(new TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>(H));
}

template <>
int TpetraSystem<double>::status() {
  mexPrintf("**** Problem ID %d [MueLu_Tpetra] ****\n", id);
  if (!A.is_null())
    mexPrintf("Matrix: %dx%d w/ %d nnz\n", A->getGlobalNumRows(), A->getGlobalNumCols(), A->getGlobalNumEntries());
  mexPrintf("Operator Complexity: %f\n", operatorComplexity);
  if (!List.is_null()) {
    mexPrintf("Parameter List:\n");
    List->print();
  }
  mexPrintf("\n");
  return IS_TRUE;
}

template <typename Scalar>
RCP<MueLu::Hierarchy<Scalar, mm_GlobalOrd, mm_LocalOrd, mm_node_t>> TpetraSystem<Scalar>::getHierarchy() {
  RCP<MueLu::TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mueluOp = rcp_static_cast<MueLu::TpetraOperator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>, Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>(prec);
  if (mueluOp.is_null())
    throw runtime_error("Tpetra precondition operator was null.");
  RCP<MueLu::Hierarchy<Scalar, mm_GlobalOrd, mm_LocalOrd, mm_node_t>> hier = mueluOp->GetHierarchy();
  if (hier.is_null())
    throw runtime_error("Hierarchy from Tpetra problem was null.");
  return hier;
}

// tpetra_complex_data_pack implementation

#ifdef HAVE_COMPLEX_SCALARS
template <>
TpetraSystem<complex_t>::TpetraSystem()
  : MuemexSystem(TPETRA_COMPLEX) {}
template <>
TpetraSystem<complex_t>::~TpetraSystem() {}

template <>
int TpetraSystem<complex_t>::status() {
  mexPrintf("**** Problem ID %d [MueLu_Tpetra (Complex Scalars)] ****\n", id);
  if (!A.is_null())
    mexPrintf("Matrix: %dx%d w/ %d nnz\n", A->getGlobalNumRows(), A->getGlobalNumCols(), A->getGlobalNumEntries());
  mexPrintf("Operator Complexity: %f\n", operatorComplexity);
  if (!List.is_null()) {
    mexPrintf("Parameter List:\n");
    List->print();
  }
  mexPrintf("\n");
  return IS_TRUE;
}
#endif  // HAVE_COMPLEX_SCALARS

// MuemexSystemList namespace implementation

void MuemexSystemList::clearAll() {
  // When items are cleared, RCPs will auto-delete the datapacks
  list.clear();
}

/* add - Adds an MuemexSystem to the list.
   Parameters:
   D       - The MuemexSystem. [I]
   Returns: problem id number of D
*/
int MuemexSystemList::add(RCP<MuemexSystem> D) {
  TEUCHOS_ASSERT(!D.is_null());
  D->id = nextID;
  nextID++;
  list.push_back(D);
  return D->id;
} /*end add*/

/* find - Finds problem by id
   Parameters:
   id      - ID number [I]
   Returns: pointer to MuemexSystem matching 'id', if found, NULL if not
   found.
*/
RCP<MuemexSystem> MuemexSystemList::find(int id) {
  if (isInList(id)) {
    for (auto problem : list) {
      if (problem->id == id)
        return problem;
    }
  }
  RCP<MuemexSystem> notFound;
  return notFound;
} /*end find*/

/* remove - Removes problem by id
   Parameters:
   id      - ID number [I]
   Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
*/
int MuemexSystemList::remove(int id) {
  int index = -1;
  for (int i = 0; i < int(list.size()); i++) {
    if (list[i]->id == id) {
      index = i;
      break;
    }
  }
  if (index == -1) {
    mexErrMsgTxt("Error: Tried to clean up a problem that doesn't exist.");
    return IS_FALSE;
  }
  list.erase(list.begin() + index);
  return IS_TRUE;
} /*end remove*/

/* size - Number of stored problems */
int MuemexSystemList::size() {
  return list.size();
}

/* Returns the status of all members of the list
   Returns IS_TRUE
*/
int MuemexSystemList::status_all() {
  // This prints all the existing problems in ascending order by ID
  for (int i = 0; i < nextID; i++) {
    for (auto problem : list) {
      if (problem->id == i) {
        problem->status();
        break;
      }
    }
  }
  return IS_TRUE;
} /*end status_all */

bool MuemexSystemList::isInList(int id) {
  bool rv = false;
  for (auto problem : list) {
    if (problem->id == id) {
      rv = true;
      break;
    }
  }
  return rv;
}

/* sanity_check - sanity checks the first couple of arguements and returns the
   program mode.
   Parameters:
   nrhs    - Number of program inputs [I]
   prhs    - The problem inputs [I]
   Return value: Which mode to run the program in.
*/

MODE_TYPE sanity_check(int nrhs, const mxArray* prhs[]) {
  MODE_TYPE rv = MODE_ERROR;
  /* Check for mode */
  if (nrhs == 0)
    mexErrMsgTxt("Error: muelu() expects at least one argument\n");
  /* Pull mode data from 1st Input */
  MODE_TYPE mode = (MODE_TYPE)loadDataFromMatlab<int>(prhs[0]);
  switch (mode) {
    case MODE_SETUP:
      if (nrhs > 1 && mxIsSparse(prhs[1])) {
        if (nrhs > 3 && mxIsSparse(prhs[2]) && mxIsSparse(prhs[3]))
          rv = MODE_ERROR;
        else
          rv = MODE_SETUP;
      } else {
        mexErrMsgTxt("Error: Invalid input for setup\n");
      }
      break;
    case MODE_SOLVE:
      // problem ID and matrix or rhs must be numeric
      if (nrhs >= 2 && mxIsNumeric(prhs[1]) && mxIsNumeric(prhs[2]))
        rv = MODE_SOLVE;
      else
        mexErrMsgTxt("Error: Invalid input for solve\n");
      break;
    case MODE_APPLY:
      // problem ID and RHS must be numeric
      if (nrhs == 3 && mxIsNumeric(prhs[1]) && mxIsNumeric(prhs[2]))
        rv = MODE_APPLY;
      else
        mexErrMsgTxt("Error: Invalid input for apply\n");
      break;
    case MODE_CLEANUP:
      if (nrhs == 1 || nrhs == 2)
        rv = MODE_CLEANUP;
      else
        mexErrMsgTxt("Error: Extraneous args for cleanup\n");
      break;
    case MODE_STATUS:
      if (nrhs == 1 || nrhs == 2)
        rv = MODE_STATUS;
      else
        mexErrMsgTxt("Error: Extraneous args for status\n");
      break;
    case MODE_AGGREGATE:
      if (nrhs > 1 && mxIsSparse(prhs[1]))
        // Uncomment the next line and remove one after when implementing aggregate mode
        // rv = MODE_AGGREGATE;
        rv = MODE_ERROR;
      else
        mexErrMsgTxt("Error: Invalid input for aggregate\n");
      break;
    case MODE_GET:
      if (nrhs < 4 || nrhs > 5)
        mexErrMsgTxt("Error: Wrong number of args for get\n");
      else
        rv = MODE_GET;
      break;
    default:
      printf("Mode number = %d\n", (int)mode);
      mexErrMsgTxt("Error: Invalid input mode\n");
  };
  return rv;
}
/*end sanity_check*/

void csc_print(int n, int* rowind, int* colptr, double* vals) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = colptr[i]; j < colptr[i + 1]; j++) {
      mexPrintf("%d %d %20.16e\n", rowind[j], i, vals[j]);
    }
  }
}

void parse_list_item(RCP<ParameterList> List, char* option_name, const mxArray* prhs) {
  // List shouldn't be NULL but if it is, initialize here
  if (List.is_null()) {
    List = rcp(new ParameterList);
  }
  mxClassID cid;
  int i, M, N, *opt_int;
  char* opt_char;
  double* opt_float;
  string opt_str;
  RCP<ParameterList> sublist = rcp(new ParameterList);
  mxArray *cell1, *cell2;
  /* Pull relevant info the the option value */
  cid = mxGetClassID(prhs);
  M   = mxGetM(prhs);
  N   = mxGetN(prhs);
  /* Add to the Teuchos list */

  // extract potential typeStr. The code is based on the assumption that
  // the typedefinition is the first word. It is necessary to distinguish
  // between "map" type (representing a Xpetra::Map) and multivector (default)
  vector<string> typestring = tokenizeList(option_name);
  std::transform(typestring[0].begin(), typestring[0].end(), typestring[0].begin(), ::tolower);
  size_t WordStart    = typestring[0].find_first_not_of(' ');
  size_t WordEnd      = typestring[0].find(' ', WordStart);
  std::string typeStr = typestring[0].substr(WordStart, WordEnd - WordStart);

  /////

  switch (cid) {
    case mxCHAR_CLASS:
      // String
      opt_char = mxArrayToString(prhs);
      opt_str  = opt_char;
      List->set(option_name, opt_str);
      if (strcmp(option_name, MUEMEX_INTERFACE) == 0) {
        if (strcmp(opt_str.c_str(), "epetra") == 0)
          useEpetra = true;
        else if (strcmp(opt_str.c_str(), "tpetra") == 0)
          useEpetra = false;
      }
      mxFree(opt_char);
      break;
    case mxDOUBLE_CLASS:
    case mxSINGLE_CLASS:
      // Single or double, real or complex
      if (mxIsComplex(prhs)) {
#ifndef HAVE_COMPLEX_SCALARS
        opt_float              = mxGetPr(prhs);
        double* opt_float_imag = mxGetPi(prhs);
        // assuming user wants std::complex<double> here...
        if (M == 1 && N == 1) {
          List->set(option_name, complex_t(*opt_float, *opt_float_imag));
        } else if (M == 0 || N == 0) {
          List->set(option_name, (complex_t*)NULL);
        } else {
          if (mxIsSparse(prhs))
            List->set(option_name, loadDataFromMatlab<RCP<Xpetra_Matrix_complex>>(prhs));
          else
            List->set(option_name, loadDataFromMatlab<RCP<Xpetra_MultiVector_complex>>(prhs));
        }
#else
        std::cerr << "Error: cannot load argument \"" << option_name << "\" because complex is not instantiated in this build.\n";
        throw std::invalid_argument("Complex not supported");
#endif
      } else {
        opt_float = mxGetPr(prhs);
        if (M == 1 && N == 1 && MMISINT(opt_float[0])) {
          List->set(option_name, (int)opt_float[0]);
        } else if (M == 1 && N == 1) {
          List->set(option_name, opt_float[0]);
        } else if (M == 0 || N == 0) {
          List->set(option_name, (double*)NULL);
        } else {
          if (mxIsSparse(prhs))
            List->set(option_name, loadDataFromMatlab<RCP<Xpetra_Matrix_double>>(prhs));
          else {
            if (typeStr == "map")  // data stored as Xpetra::Map type
              List->set(option_name, loadDataFromMatlab<RCP<Xpetra_map>>(prhs));
            else  // data stored as MultiVector
              List->set(option_name, loadDataFromMatlab<RCP<Xpetra_MultiVector_double>>(prhs));
          }
        }
      }
      break;
    case mxLOGICAL_CLASS:
      // Bool
      if (M == 1 && N == 1)
        List->set(option_name, mxIsLogicalScalarTrue(prhs));
      else
        List->set(option_name, mxGetLogicals(prhs));
      // NTS: The else probably doesn't work.
      break;
    case mxINT8_CLASS:
    case mxUINT8_CLASS:
    case mxINT16_CLASS:
    case mxUINT16_CLASS:
    case mxINT32_CLASS:
    case mxUINT32_CLASS:
      // Integer
      opt_int = (int*)mxGetData(prhs);
      if (M == 1 && N == 1)
        List->set(option_name, *opt_int);
#ifdef HAVE_MUELU_INTREPID2
      else if (strcmp(option_name, "pcoarsen: element to node map") == 0)
        List->set(option_name, loadDataFromMatlab<RCP<FieldContainer_ordinal>>(prhs));
#endif
      else
        List->set(option_name, loadDataFromMatlab<RCP<Xpetra_ordinal_vector>>(prhs));
      break;
      // NTS: 64-bit ints will break on a 32-bit machine.        We
      // should probably detect machine type, or somthing, but that would
      // involve a non-trivial quantity of autoconf kung fu.
    case mxCELL_CLASS:
      // Interpret a cell list as a nested teuchos list.
      // NTS: Assuming that it's a 1D row ordered array
      for (i = 0; i < N; i += 2) {
        cell1 = mxGetCell(prhs, i);
        cell2 = mxGetCell(prhs, i + 1);
        if (!mxIsChar(cell1))
          mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");
        opt_char = mxArrayToString(cell1);
        parse_list_item(sublist, opt_char, cell2);
        List->set(option_name, *sublist);
        mxFree(opt_char);
      }
      break;
    case mxINT64_CLASS:
    case mxUINT64_CLASS:
    case mxFUNCTION_CLASS:
    case mxUNKNOWN_CLASS:
    case mxSTRUCT_CLASS:
      // Currently Graph and Aggregates are stored as structures
      if (isValidMatlabAggregates(prhs)) {
        try {
          List->set(option_name, loadDataFromMatlab<RCP<MAggregates>>(prhs));
          break;
        } catch (exception& e) {
          cout << e.what();
          throw runtime_error("Parsing aggregates in parameter list failed.");
        }
      } else if (isValidMatlabGraph(prhs)) {
        try {
          List->set(option_name, loadDataFromMatlab<RCP<MGraph>>(prhs));
          break;
        } catch (exception& e) {
          cout << e.what();
          throw runtime_error("Parsing graph in parameter list failed.");
        }
      }
    default:
      mexPrintf("Error parsing input option: %s [type=%d]\n", option_name, cid);
      mexErrMsgTxt("Error: An input option is invalid!\n");
  };
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
/* build_teuchos_list - takes the inputs (barring the solver mode and
   matrix/rhs) and turns them into a Teuchos list.
   Parameters:
   nrhs    - Number of program inputs [I]
   prhs    - The problem inputs [I]
   Return value: Teuchos list containing all parameters passed in by the user.
*/
RCP<ParameterList> build_teuchos_list(int nrhs, const mxArray* prhs[]) {
  RCP<ParameterList> TPL = rcp(new ParameterList);
  char* option_name;
  for (int i = 0; i < nrhs; i += 2) {
    if (i == nrhs - 1 || !mxIsChar(prhs[i]))
      mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");
    /* What option are we setting? */
    option_name = mxArrayToString(prhs[i]);
    /* Parse */
    parse_list_item(TPL, option_name, prhs[i + 1]);
    /* Free memory */
    mxFree(option_name);
  }
  TPL->print(cout);
  return TPL;
}
/*end build_teuchos_list*/

}  // namespace MueLu

using namespace MueLu;  //...but give mexFunction access to all MueLu members defined above

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // Lazily initialize Tpetra
  if (!Tpetra::isInitialized()) {
    int argc    = 0;
    char** argv = NULL;
    Tpetra::initialize(&argc, &argv);
  }
  double* id;
  int rv;
  // Arrays representing vectors
  string intf;
  RCP<ParameterList> List;
  MODE_TYPE mode;
  RCP<MuemexSystem> D;
  /* Sanity Check Input */
  mode = sanity_check(nrhs, prhs);

#define FORCE_FACTORIES_TO_COMPILE
#ifdef FORCE_FACTORIES_TO_COMPILE
  {
    // debug
    MueLu::TwoLevelMatlabFactory<double, int, int> f1;
    MueLu::SingleLevelMatlabFactory<double, int, int> f2;
  }
#endif

  switch (mode) {
    case MODE_SETUP: {
      try {
        double oc       = 0;
        bool haveCoords = false;
        if (nrhs == 2)
          List = rcp(new ParameterList);
        else if (nrhs == 3) {
          if (!mxIsNumeric(prhs[2]) || mxIsSparse(prhs[2]) || mxIsComplex(prhs[2]))
            throw runtime_error("Expected real-valued, dense Coordinates array as the third muelu argument");
          else {
            haveCoords = true;
            List       = rcp(new ParameterList);
          }
        } else {
          if (mxIsNumeric(prhs[2]) && !mxIsSparse(prhs[2]) && !mxIsComplex(prhs[2])) {
            List       = build_teuchos_list(nrhs - 3, &(prhs[3]));
            haveCoords = true;
          } else {
            // assume that the parameters start at third argument if it doesn't seem like coords are there
            List = build_teuchos_list(nrhs - 2, &(prhs[2]));
          }
        }
        // Combine xml and easy parameter lists
        if (List->isParameter("xml parameter file")) {
          RCP<ParameterList> xmlParams = Teuchos::getParametersFromXmlFile(List->get<string>("xml parameter file"));
          List->remove("xml parameter file");
          List->setParametersNotAlreadySet(*xmlParams);
        }
        if (mxIsComplex(prhs[1])) {
          // Abort if input is complex but complex isn't supported
#ifndef HAVE_COMPLEX_SCALARS
          mexPrintf("Error: Complex scalars unsupported by this build of Trilinos.\n");
          throw runtime_error("Complex scalars not supported.");
#endif
        }
        intf = List->get(MUEMEX_INTERFACE, "tpetra");
        List->remove(MUEMEX_INTERFACE);  // no longer need this parameter
        if (intf == "epetra") {
          if (mxIsComplex(prhs[1])) {
            mexPrintf("Error: Attempting to use complex-valued matrix with Epetra, which is unsupported.\n");
            mexPrintf("Use Tpetra with complex matrices instead.\n");
            throw runtime_error("Tried to use complex matrix with Epetra");
          }
          RCP<EpetraSystem> dp = rcp(new EpetraSystem());
          dp->List             = List;
          dp->setup(prhs[1], haveCoords, haveCoords ? prhs[2] : (mxArray*)NULL);
          oc = dp->operatorComplexity;
          D  = rcp_implicit_cast<MuemexSystem>(dp);
        } else if (intf == "tpetra") {
          // infer scalar type from prhs (can be double or complex<double>)
          if (mxIsComplex(prhs[1])) {
#ifdef HAVE_COMPLEX_SCALARS
            RCP<TpetraSystem<complex_t>> dp = rcp(new TpetraSystem<complex_t>());
            dp->List                        = List;
            dp->setup(prhs[1], haveCoords, haveCoords ? prhs[2] : (mxArray*)NULL);
            oc = dp->operatorComplexity;
            D  = rcp_implicit_cast<MuemexSystem>(dp);
#else
            throw runtime_error("Complex scalars not supported.");
#endif
          } else {
            RCP<TpetraSystem<double>> dp = rcp(new TpetraSystem<double>());
            dp->List                     = List;
            dp->setup(prhs[1], haveCoords, haveCoords ? prhs[2] : (mxArray*)NULL);
            oc = dp->operatorComplexity;
            D  = rcp_implicit_cast<MuemexSystem>(dp);
          }
        }
        rv = MuemexSystemList::add(D);
        mexPrintf("Set up problem #%d\n", rv);
        if (nlhs > 0) {
          plhs[0]                     = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
          *((int*)mxGetData(plhs[0])) = rv;
          // output OC also
          if (nlhs > 1)
            plhs[1] = mxCreateDoubleScalar(oc);
        }
        mexLock();
      } catch (exception& e) {
        mexPrintf("An error occurred during setup routine:\n");
        cout << e.what() << endl;
        if (nlhs > 0) {
          plhs[0]                     = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
          *((int*)mxGetData(plhs[0])) = -1;
        }
        if (nlhs > 1)
          plhs[1] = mxCreateDoubleScalar(0);
      }
      break;
    }
    case MODE_SOLVE: {
      int iters;
      try {
        bool reuse;
        // MODE_SOLVE, probID, matrix, bVec, params OR
        // MODE_SOLVE, probID, bVec, params
        // prhs[0] holds the MODE_SOLVE enum value
        // if reusing, either nrhs == 3 or prhs[3] is not a string
        //(because bVec can't be a string)
        if (MuemexSystemList::size() == 0)
          throw runtime_error("No linear systems are set up.");
        if (nrhs == 3 || (nrhs > 3 && mxGetClassID(prhs[3]) == mxCHAR_CLASS)) {
          // No matrix supplied as argument, use one from setup
          if (nrhs > 3)
            List = build_teuchos_list(nrhs - 3, &prhs[3]);
          else
            List = rcp(new ParameterList);
          reuse = true;
        } else {
          if (nrhs > 4)
            List = build_teuchos_list(nrhs - 4, &prhs[4]);
          else
            List = rcp(new ParameterList);
          reuse = false;
        }
        if (List->isType<string>("Output Style")) {
          int type = strToOutputStyle(List->get("Output Style", "Belos::Brief").c_str());
          List->remove("Output Style");
          // Reset the ParameterList entry to be of type int instead of string
          List->set("Output Style", type);
        }
        // Convert Belos msg type string to int in List
        // Note: if the parameter value is already an int, don't touch it.
        if (List->isType<string>("Verbosity")) {
          // Errors + Warnings is already the default Belos verbosity setting
          int verb = getBelosVerbosity(List->get("Verbosity", "Belos::Errors + Belos::Warnings").c_str());
          List->remove("Verbosity");
          List->set("Verbosity", verb);
        }
        int probID           = loadDataFromMatlab<int>(prhs[1]);
        RCP<MuemexSystem> dp = MuemexSystemList::find(probID);
        if (dp.is_null())
          throw runtime_error("Problem handle not allocated.");
        // get pointer to MATLAB array that will be "B" or "rhs" multivector
        const mxArray* rhs = reuse ? prhs[2] : prhs[3];
        switch (dp->type) {
          case EPETRA: {
            RCP<EpetraSystem> esys = rcp_static_cast<EpetraSystem, MuemexSystem>(dp);
            RCP<Epetra_CrsMatrix> matrix;
            if (reuse)
              matrix = esys->GetMatrix();
            else
              matrix = loadDataFromMatlab<RCP<Epetra_CrsMatrix>>(prhs[2]);
            plhs[0] = esys->solve(List, matrix, rhs, iters);
            break;
          }
          case TPETRA: {
            RCP<TpetraSystem<double>> tsys = rcp_static_cast<TpetraSystem<double>, MuemexSystem>(dp);
            RCP<Tpetra_CrsMatrix_double> matrix;
            if (reuse)
              matrix = tsys->GetMatrix();
            else
              matrix = loadDataFromMatlab<RCP<Tpetra_CrsMatrix_double>>(prhs[2]);
            plhs[0] = tsys->solve(List, matrix, rhs, iters);
            break;
          }
          case TPETRA_COMPLEX: {
#ifdef HAVE_COMPLEX_SCALARS
            RCP<TpetraSystem<complex_t>> tsys = rcp_static_cast<TpetraSystem<complex_t>, MuemexSystem>(dp);
            RCP<Tpetra_CrsMatrix_complex> matrix;
            if (reuse)
              matrix = tsys->GetMatrix();
            else
              matrix = loadDataFromMatlab<RCP<Tpetra_CrsMatrix_complex>>(prhs[2]);
            plhs[0] = tsys->solve(List, matrix, rhs, iters);
            break;
#else
            std::cerr << "Cannot solve complex-valued system because complex is not enabled in this build.\n";
            throw std::invalid_argument("Complex not supported");
#endif
          }
        }
        if (nlhs > 1)
          plhs[1] = saveDataToMatlab<int>(iters);
      } catch (exception& e) {
        mexPrintf("An error occurred during the solve routine:\n");
        cout << e.what() << endl;
      }
      break;
    }
    case MODE_APPLY: {
      try {
        // MODE_APPLY, probID, rhsVec
        // prhs[0] holds the MODE_APPLY enum value
        if (MuemexSystemList::size() == 0)
          throw runtime_error("No linear systems are set up.");
        int probID           = loadDataFromMatlab<int>(prhs[1]);
        RCP<MuemexSystem> dp = MuemexSystemList::find(probID);
        if (dp.is_null())
          throw runtime_error("Problem handle not allocated.");
        // get pointer to MATLAB array that will be "B" or "rhs" multivector
        const mxArray* rhs = prhs[2];
        switch (dp->type) {
          case EPETRA: {
            RCP<EpetraSystem> esys = rcp_static_cast<EpetraSystem, MuemexSystem>(dp);
            plhs[0]                = esys->apply(rhs);
            break;
          }
          case TPETRA: {
            RCP<TpetraSystem<double>> tsys = rcp_static_cast<TpetraSystem<double>, MuemexSystem>(dp);
            plhs[0]                        = tsys->apply(rhs);
            break;
          }
          case TPETRA_COMPLEX: {
#ifdef HAVE_COMPLEX_SCALARS
            RCP<TpetraSystem<complex_t>> tsys = rcp_static_cast<TpetraSystem<complex_t>, MuemexSystem>(dp);
            plhs[0]                           = tsys->apply(rhs);
            break;
#else
            std::cerr << "Cannot solve complex-valued system because complex is not enabled in this build.\n";
            throw std::invalid_argument("Complex not supported");
#endif
          }
        }
      } catch (exception& e) {
        mexPrintf("An error occurred during the apply routine:\n");
        cout << e.what() << endl;
      }
      break;
    }
    case MODE_CLEANUP: {
      try {
        mexPrintf("MueMex in cleanup mode.\n");
        if (MuemexSystemList::size() > 0 && nrhs == 1) {
          /* Cleanup all problems */
          for (int i = 0; i < MuemexSystemList::size(); i++)
            mexUnlock();
          MuemexSystemList::clearAll();
          rv = 1;
        } else if (MuemexSystemList::size() > 0 && nrhs == 2) {
          /* Cleanup one problem */
          int probID = loadDataFromMatlab<int>(prhs[1]);
          mexPrintf("Cleaning up problem #%d\n", probID);
          rv = MuemexSystemList::remove(probID);
          if (rv)
            mexUnlock();
        } /*end elseif*/
        else {
          rv = 0;
        }
        /* Set return value */
        plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        id      = (double*)mxGetData(plhs[0]);
        *id     = double(rv);
      } catch (exception& e) {
        mexPrintf("An error occurred during the cleanup routine:\n");
        cout << e.what() << endl;
        if (nlhs > 0) {
          plhs[0]                   = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
          *(int*)mxGetData(plhs[0]) = 0;
        }
      }
      break;
    }
    case MODE_STATUS: {
      try {
        // mexPrintf("MueMex in status checking mode.\n");
        if (MuemexSystemList::size() > 0 && nrhs == 1) {
          /* Status check on all problems */
          rv = MuemexSystemList::status_all();
        } /*end if*/
        else if (MuemexSystemList::size() > 0 && nrhs == 2) {
          /* Status check one problem */
          int probID = loadDataFromMatlab<int>(prhs[1]);
          D          = MuemexSystemList::find(probID);
          if (D.is_null())
            throw runtime_error("Error: Problem handle not allocated.\n");
          rv = D->status();
        } /*end elseif*/
        else
          mexPrintf("No problems set up.\n");
        if (nlhs > 0) {
          int outVal = 0;
          plhs[0]    = saveDataToMatlab<int>(outVal);
        }
      } catch (exception& e) {
        mexPrintf("An error occurred during the status routine:\n");
        cout << e.what() << endl;
        int outVal = -1;
        if (nlhs > 0)
          plhs[0] = saveDataToMatlab<int>(outVal);
      }
      break;
    }
    case MODE_GET: {
      try {
        int probID            = loadDataFromMatlab<int>(prhs[1]);
        int levelID           = loadDataFromMatlab<int>(prhs[2]);
        char* dataName        = mxArrayToString(prhs[3]);
        MuemexType outputType = INT;
        RCP<MuemexSystem> dp  = MuemexSystemList::find(probID);
        if (dp.is_null()) {
          throw runtime_error("Problem handle not allocated.");
        }

        // See if typeHint was given
        char* paramTypeName = NULL;
        if (nrhs > 4) {
          paramTypeName = mxArrayToString(prhs[4]);
          mexPrintf("paramTypeName %s", paramTypeName);
        }
        bool complexFlag = dp->type == TPETRA_COMPLEX;
        // std::cout << "before strToDataType dataName=" << dataName << " paramTypeName " << paramTypeName << std::endl;

        outputType = strToDataType(dataName, paramTypeName, complexFlag);
        plhs[0]    = dp->getHierarchyData(string(dataName), outputType, levelID);
      } catch (exception& e) {
        mexPrintf("An error occurred during the get routine:\n");
        cout << e.what() << endl;
        plhs[0] = mxCreateDoubleScalar(0);
      }
      break;
    }
    case MODE_ERROR:
      mexPrintf("MueMex error.");
      break;
    case MODE_AGGREGATE:
    default:
      mexPrintf("Mode not supported yet.");
  }
}
