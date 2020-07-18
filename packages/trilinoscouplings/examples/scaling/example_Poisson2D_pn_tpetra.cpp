// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   example_Poisson.cpp
    \brief  Example solution of a Poisson equation on a quad mesh using
    nodal (Hgrad) elements.

    This example uses the following Trilinos packages:
    \li     Pamgen to generate a Quad mesh.
    \li     Sacado to form the source term from user-specified manufactured solution.
    \li     Intrepid to build the discretization matrix and right-hand side.
    \li     Tpetra to handle the global matrix and vector.
    \li     Isorropia to partition the matrix. (Optional)
    \li     ML to solve the linear system.


    \verbatim

    Poisson system:

    div A grad u = f in Omega
    u = g on Gamma

    where
    A is a symmetric, positive definite material tensor
    f is a given source term

    Corresponding discrete linear system for nodal coefficients(x):

    Kx = b

    K - HGrad stiffness matrix
    b - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson, D. Hensinger, C. Siefert.

    \remark Usage:
    \code   ./TrilinosCouplings_examples_scaling_example_Poisson.exe \endcode

    \remark Example driver requires input file named Poisson2D.xml with Pamgen
    formatted mesh description and settings for Isorropia (a version
    is included in the Trilinos repository with this driver).

    \remark The exact solution (u) and material tensor (A) are set in the
    functions "exactSolution" and "materialTensor" and may be
    modified by the user.

*/

/*** Uncomment if you would like output data for plotting ***/
//#define DUMP_DATA

/**************************************************************/
/*                          Includes                          */
/**************************************************************/
#include <stdlib.h>

// TrilinosCouplings includes
#include "TrilinosCouplings_config.h"
#include "TrilinosCouplings_Pamgen_Utils.hpp"

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Intrepid2 Includes
#include "Kokkos_DynRankView.hpp"

// Tpetra includes
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// Pamgen includes
#include "create_inline_mesh.h"
#include "pamgen_im_exodusII_l.h"
#include "pamgen_im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// AztecOO includes
#include "AztecOO.h"

// MueLu Includes
#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
#  include "MueLu.hpp"
#  include "MueLu_ParameterListInterpreter.hpp"
#  include "MueLu_TpetraOperator.hpp"
#  include "MueLu_Utilities.hpp"
#  include "MueLu_HierarchyManager.hpp"
#  include "MueLu_FactoryManagerBase.hpp"
#  include "MueLu_CreateTpetraPreconditioner.hpp"
// MueLu/Avatar Includes
#ifdef HAVE_TRILINOSCOUPLINGS_AVATAR
#  include "MueLu_AvatarInterface.hpp"
#endif

#endif // HAVE_TRILINOSCOUPLINGS_MUELU

#ifdef HAVE_INTREPID_KOKKOSCORE
#include "Sacado.hpp"
#else
// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#endif

//#if defined(HAVE_TRINOSCOUPLINGS_BELOS) && defined(HAVE_TRILINOSCOUPLINGS_MUELU)
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosFixedPointSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
//#endif

using namespace std;
using namespace Intrepid;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::TimeMonitor;
using Teuchos::ParameterList;

// for debugging
#define zwrap(x) (std::abs(x)<1e-10 ? 0.0 : x)

// Statistics for machine learning
const int NUM_STATISTICS = 8;
std::vector<double> local_stat_max(NUM_STATISTICS);
std::vector<double> local_stat_min(NUM_STATISTICS);
std::vector<double> local_stat_sum(NUM_STATISTICS);


// Disable the full set of problem statistics
const bool use_new_problem_stats = true;

struct fecomp{
  bool operator () ( topo_entity* x,  topo_entity*  y )const
  {
    if(x->sorted_local_node_ids < y->sorted_local_node_ids)return true;
    return false;
  }
};

template<class FC>
double distance2(const FC & coord, int n1, int n2) {
  double dist = 0.0;
  for(int i=0; i<coord.dimension(1); i++)
    dist += (coord(n2,i) -coord(n1,i)) * (coord(n2,i) -coord(n1,i));
  return sqrt(dist);
}


/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef double scalar_type;
typedef Teuchos::ScalarTraits<scalar_type> ScalarTraits;
using local_ordinal_type = Tpetra::Map<>::local_ordinal_type;
using global_ordinal_type = Tpetra::Map<>::global_ordinal_type;
typedef KokkosClassic::DefaultNode::DefaultNodeType NO;
typedef Sacado::Fad::SFad<double,2>      Fad2; //# ind. vars fixed at 2
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;

typedef Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,NO>    operator_type;
typedef Tpetra::CrsGraph<local_ordinal_type,global_ordinal_type,NO>   crs_graph_type;
typedef Tpetra::CrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>   crs_matrix_type;
typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,NO>      vector_type;
typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,NO> multivector_type;
typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,NO>            driver_map_type;
Tpetra::global_size_t INVALID_GO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
Tpetra::global_size_t INVALID_LO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();


typedef Xpetra::Matrix<scalar_type,local_ordinal_type,global_ordinal_type,NO> xpetra_crs_matrix_type;
typedef Xpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,NO> xpetra_multivector_type;

typedef Belos::LinearProblem<scalar_type, multivector_type, operator_type> linear_problem_type;

typedef MueLu::TpetraOperator<scalar_type,local_ordinal_type,global_ordinal_type,NO> muelu_tpetra_operator;

//typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;

#define TC_sumAll(rcpComm, in, out) \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out))
#define TC_minAll(rcpComm, in, out) \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out))
#define TC_maxAll(rcpComm, in, out) \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out))


double myDistance2(const multivector_type &v, int i0, int i1) {
 const size_t numVectors = v.getNumVectors();
 double distance = 0.0;
 for (size_t j=0; j<numVectors; j++) {
   distance += (v.getData(j)[i0]-v.getData(j)[i1])*(v.getData(j)[i0]-v.getData(j)[i1]);
 }
 return distance;
}

// forward declarations
void PromoteMesh_Pn_Kirby(const int degree,  const EPointType & pointType, const FieldContainer<int> & P1_elemToNode, const FieldContainer<double> & P1_nodeCoord, const FieldContainer<double> & P1_edgeCoord,  const FieldContainer<int> & P1_elemToEdge,  const FieldContainer<int> & P1_elemToEdgeOrient, const FieldContainer<int> & P1_nodeOnBoundary,
                          FieldContainer<int> & Pn_elemToNode, FieldContainer<double> & Pn_nodeCoord,FieldContainer<int> & Pn_nodeOnBoundary, std::vector<int> &Pn_edgeNodes, std::vector<int> &Pn_cellNodes);

void GenerateEdgeEnumeration(const FieldContainer<int> & elemToNode, const FieldContainer<double> & nodeCoord, FieldContainer<int> & elemToEdge, FieldContainer<int> & elemToEdgeOrient, FieldContainer<int> & edgeToNode, FieldContainer<double> & edgeCoord);

void CreateP1MeshFromPnMesh(const int degree, const FieldContainer<int> & P2_elemToNode, FieldContainer<int> &aux_P1_elemToNode);

void CreateLinearSystem(int numWorkSets,
                        int desiredWorksetSize,
                        FieldContainer<int> const &elemToNode,
                        FieldContainer<double> const &nodeCoord,
                        FieldContainer<double> const &cubPoints,
                        FieldContainer<double> const &cubWeights,
                        RCP<Basis<double,FieldContainer<double> > > &myBasis_rcp,
                        FieldContainer<double> const &HGBGrads,
                        FieldContainer<double> const &HGBValues,
                        std::vector<global_ordinal_type>       const &globalNodeIds,
                        crs_matrix_type &StiffMatrix,
                        RCP<multivector_type> &rhsVector,
                        std::string &msg
                        );


void GenerateIdentityCoarsening_pn_to_p1(const FieldContainer<int> & P2_elemToNode,
                      //const FieldContainer<int> & P1_elemToNode,
                      RCP<const driver_map_type> const & P1_map_aux, RCP<const driver_map_type> const &P2_map,
                      RCP<crs_matrix_type> & Interpolation,
                      RCP<crs_matrix_type> & Restriction);

//GenerateIdentityCoarsening_pn_to_p1(elemToNode, StiffMatrix_aux.getDomainMap(), StiffMatrix.getRangeMap(), P_identity, R_identity);

void Apply_Dirichlet_BCs(std::vector<int> &BCNodes, crs_matrix_type & A, multivector_type & x, multivector_type & b,
                         const multivector_type & soln);

/**********************************************************************************/
/***************** FUNCTION DECLARATION FOR ML PRECONDITIONER *********************/
/**********************************************************************************/

/** \brief  ML Preconditioner

    \param  ProblemType        [in]    problem type
    \param  MLList             [in]    ML parameter list
    \param  A                  [in]    discrete operator matrix
    \param  nCoord             [in]    Nodal Coordinates
    \param  xexact             [in]    exact solution
    \param  b                  [in]    right-hand-side vector
    \param  maxits             [in]    max iterations
    \param  tol                [in]    solver tolerance
    \param  uh                 [in/out]   solution vector
    \param  TotalErrorResidual [out]   error residual
    \param  TotalErrorExactSol [out]   error in uh

*/
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
                                 ParameterList   & AMGList,
                                 RCP<crs_matrix_type>   const & A,
                                 RCP<multivector_type> & nCoord,
                                 RCP<multivector_type> const & xexact,
                                 RCP<multivector_type> & b,
				 int maxits,
				 double tol,	
                                 RCP<multivector_type> & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol,
				 std::string &amgType,
				  std::string &solveType);
					


/**********************************************************************************/
/******** FUNCTION DECLARATIONS FOR EXACT SOLUTION AND SOURCE TERMS ***************/
/**********************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point

    \return Value of the exact solution at (x,y)
*/
template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y);

/** \brief  User-defined material tensor.

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y).
*/
template<typename Scalar>
void materialTensor(Scalar material[][2], const Scalar&  x, const Scalar&  y);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
*/
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[2], const Scalar& x, const Scalar& y);


/** \brief Computes source term: f = -div(A.grad u).  Requires user-defined exact solution
    and material tensor.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y)
*/
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y);


/** \brief Computation of the material tensor at array of points in physical space.

    \param worksetMaterialValues      [out]     Rank-2, 3 or 4 array with dimensions (C,P), (C,P,D) or (C,P,D,D)
    with the values of the material tensor
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        worksetMaterialValues,
                            const ArrayIn &   evaluationPoints);


/** \brief Computation of the source term at array of points in physical space.

    \param sourceTermValues           [out]     Rank-2 (C,P) array with the values of the source term
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints);

/** \brief Computation of the exact solution at array of points in physical space.

    \param exactSolutionValues        [out]     Rank-2 (C,P) array with the values of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints);


/** \brief Computation of the gradient of the exact solution at array of points in physical space.

    \param exactSolutionGradValues    [out]     Rank-3 (C,P,D) array with the values of the gradient of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints);



// Copy field containers
template<class FC1, class FC2>
void CopyFieldContainer2D(const FC1 & c1, FC2 & c2) {
  Kokkos::resize(c2,c1.dimension(0),c1.dimension(1));
  for(size_t i=0; i<(size_t)c1.dimension(0); i++)
    for(size_t j=0; j<(size_t)c1.dimension(1); j++)
      c2(i,j) = c1(i,j);
}



/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/

int main(int argc, char *argv[]) {

  int error = 0;
  int numProcs=1;
  int rank=0;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
  int MyPID = Comm->getRank();

  Teuchos::CommandLineProcessor clp(false);
  Teuchos::ParameterList problemStatistics;
  
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Pamgen Setup")));

  std::string optMeshFile = "Poisson2D.xml";
  clp.setOption("mesh",  &optMeshFile, "xml file with description of Pamgen mesh");
  std::string optSolverFile  = "";
  clp.setOption("solver",   &optSolverFile,  "xml file containing linear solver options");
  bool optPrintTimings = false;
  clp.setOption("timings", "notimings",  &optPrintTimings,      "print timer summary");

  // If matrixFilename is nonempty, dump the matrix to that file
  // in MatrixMarket format.
  std::string matrixFilename;
  clp.setOption ("matrixFilename", &matrixFilename, "If nonempty, dump the "
		  "generated matrix to that file in MatrixMarket format.");

  // If coordsFilename is nonempty, dump the rhs to that file
  // in MatrixMarket format.
  std::string rhsFilename;
  clp.setOption ("rhsFilename", &rhsFilename, "If nonempty, dump the "
		  "generated rhs to that file in MatrixMarket format.");
  
  // If coordsFilename is nonempty, dump the coords to that file
  // in MatrixMarket format.
  std::string coordsFilename;
  clp.setOption ("coordsFilename", &coordsFilename, "If nonempty, dump the "
		  "generated coordinates to that file in MatrixMarket format.");
  
  // Random number seed
  int randomSeed=24601;
  clp.setOption ("seed", &randomSeed, "Random Seed.");
  
  
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }
  
  // Initialize RNG
  srand(randomSeed);

  if (MyPID == 0){
    std::cout \
      << "===============================================================================\n" \
      << "|                                                                             |\n" \
      << "|          Example: Solve Poisson Equation on Hexahedral Mesh                 |\n" \
      << "|                                                                             |\n" \
      << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
      << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
      << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
      << "|                                                                             |\n" \
      << "|  Intrepid's website:   http://trilinos.sandia.gov/packages/intrepid         |\n" \
      << "|  Pamgen's website:     http://trilinos.sandia.gov/packages/pamgen           |\n" \
      << "|  ML's website:         http://trilinos.sandia.gov/packages/ml               |\n" \
      << "|  Isorropia's website:  http://trilinos.sandia.gov/packages/isorropia        |\n" \
      << "|  Trilinos website:     http://trilinos.sandia.gov                           |\n" \
      << "|                                                                             |\n" \
      << "===============================================================================\n";
  }


#ifdef HAVE_MPI
  if (MyPID == 0) {
    std::cout << "PARALLEL executable \n";
  }
#else
  if (MyPID == 0) {
    std::cout << "SERIAL executable \n";
  }
#endif

  /**********************************************************************************/
  /********************************** GET XML INPUTS ********************************/
  /**********************************************************************************/

  // Command line for xml file, otherwise use default
  std::string   xmlMeshInFileName(optMeshFile), xmlSolverInFileName(optSolverFile);

  // Read xml file into parameter list
  ParameterList inputMeshList;
  ParameterList inputSolverList;

  if(xmlMeshInFileName.length()) {
    if (MyPID == 0) {
      std::cout << "\nReading parameter list from the XML file \""<<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile (xmlMeshInFileName, Teuchos::ptr (&inputMeshList));
    if (MyPID == 0) {
      inputMeshList.print(std::cout,2,true,true);
      std::cout << "\n";
    }
  }
  else
    {
      std::cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
      return 0;
    }

  if(xmlSolverInFileName.length()) {
    if (MyPID == 0)
      std::cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n\n";
    Teuchos::updateParametersFromXmlFile(xmlSolverInFileName, Teuchos::inoutArg(inputSolverList));
  } else if (MyPID == 0) std::cout << "Using default solver values ..." << std::endl;

  // Get pamgen mesh definition
  std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,"meshInput");
  int degree = inputMeshList.get("degree",2);

  // Get Isorropia and Zoltan parameters.
  ParameterList iso_paramlist = inputMeshList.sublist
    ("Isorropia Input") ;
  if (MyPID == 0) {
    std::cout << "Isorropia/Zoltan parameters" << std::endl;
    iso_paramlist.print(std::cout,2,true,true);
  }


  /**********************************************************************************/
  /***************************** GET CELL TOPOLOGY **********************************/
  /**********************************************************************************/

  // Get cell topology for base hexahedron
  shards::CellTopology P1_cellType(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  shards::CellTopology Pn_cellType(shards::getCellTopologyData<shards::Quadrilateral<> >() );
  assert(P1_cellType.getDimension() == Pn_cellType.getDimension());

  // Get dimensions
  int P1_numNodesPerElem = P1_cellType.getNodeCount();
  int P1_numEdgesPerElem = P1_cellType.getEdgeCount();
  int P1_numNodesPerEdge = 2; // for any rational universe
  int spaceDim = P1_cellType.getDimension();
  int dim = 2;


  // Build reference element edge to node map
  FieldContainer<int> refEdgeToNode(P1_numEdgesPerElem, P1_numNodesPerEdge);
  for (int i=0; i<P1_numEdgesPerElem; i++){
    refEdgeToNode(i,0)=P1_cellType.getNodeMap(1, i, 0);
    refEdgeToNode(i,1)=P1_cellType.getNodeMap(1, i, 1);
  }

  /**********************************************************************************/
  /******************************* GENERATE MESH ************************************/
  /**********************************************************************************/

  if (MyPID == 0) {
    std::cout << "Generating mesh ... \n\n";
  }

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  long long cr_result = Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
  TrilinosCouplings::pamgen_error_check(std::cout,cr_result);

  string msg("Poisson: ");
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Mesh queries")));

  // Get mesh size info
  char title[100];
  long long numDim;
  long long numNodes;
  long long numElems;
  long long numElemBlk;
  long long numNodeSets;
  long long numSideSets;
  int id = 0;

  im_ex_get_init_l(id, title, &numDim, &numNodes,
                   &numElems, &numElemBlk, &numNodeSets,
                   &numSideSets);

  long long numNodesGlobal;
  long long numElemsGlobal;
  long long numElemBlkGlobal;
  long long numNodeSetsGlobal;
  long long numSideSetsGlobal;

  im_ne_get_init_global_l(id, &numNodesGlobal, &numElemsGlobal,
                          &numElemBlkGlobal, &numNodeSetsGlobal,
                          &numSideSetsGlobal);

  // Print mesh information
  if (MyPID == 0){
    std::cout << " Number of Global Elements: " << numElemsGlobal << " \n";
    std::cout << " Number of Global Nodes: " << numNodesGlobal << " \n\n";
  }

  long long * block_ids = new long long [numElemBlk];
  error += im_ex_get_elem_blk_ids_l(id, block_ids);


  long long  *nodes_per_element   = new long long [numElemBlk];
  long long  *element_attributes  = new long long [numElemBlk];
  long long  *elements            = new long long [numElemBlk];
  char      **element_types       = new char * [numElemBlk];
  long long **elmt_node_linkage   = new long long * [numElemBlk];


  for(long long i = 0; i < numElemBlk; i ++){
    element_types[i] = new char [MAX_STR_LENGTH + 1];
    error += im_ex_get_elem_block_l(id,
                                    block_ids[i],
                                    element_types[i],
                                    (long long*)&(elements[i]),
                                    (long long*)&(nodes_per_element[i]),
                                    (long long*)&(element_attributes[i]));
  }

  /*connectivity*/
  for(long long b = 0; b < numElemBlk; b++){
    elmt_node_linkage[b] =  new long long [nodes_per_element[b]*elements[b]];
    error += im_ex_get_elem_conn_l(id,block_ids[b],elmt_node_linkage[b]);
  }

  // Get node-element connectivity
  int telct = 0;
  FieldContainer<int> P1_elemToNode(numElems,P1_numNodesPerElem);
  for(long long b = 0; b < numElemBlk; b++){
    for(long long el = 0; el < elements[b]; el++){
      for (int j=0; j<P1_numNodesPerElem; j++) {
        P1_elemToNode(telct,j) = elmt_node_linkage[b][el*P1_numNodesPerElem + j]-1;
      }
      telct ++;
    }
  }

  // Read node coordinates and place in field container
  FieldContainer<double> P1_nodeCoord(numNodes,dim);
  double * nodeCoordx = new double [numNodes];
  double * nodeCoordy = new double [numNodes];
  im_ex_get_coord_l(id,nodeCoordx,nodeCoordy,0);
  for (int i=0; i<numNodes; i++) {
    P1_nodeCoord(i,0)=nodeCoordx[i];
    P1_nodeCoord(i,1)=nodeCoordy[i];
  }
  /*parallel info*/
  long long num_internal_nodes;
  long long num_border_nodes;
  long long num_external_nodes;
  long long num_internal_elems;
  long long num_border_elems;
  long long num_node_comm_maps;
  long long num_elem_comm_maps;
  im_ne_get_loadbal_param_l( id,
                             &num_internal_nodes,
                             &num_border_nodes,
                             &num_external_nodes,
                             &num_internal_elems,
                             &num_border_elems,
                             &num_node_comm_maps,
                             &num_elem_comm_maps,
                             0/*unused*/ );

  if(num_node_comm_maps > 0){
    node_comm_proc_ids   = new long long  [num_node_comm_maps];
    node_cmap_node_cnts  = new long long  [num_node_comm_maps];
    node_cmap_ids        = new long long  [num_node_comm_maps];
    comm_node_ids        = new long long* [num_node_comm_maps];
    comm_node_proc_ids   = new long long* [num_node_comm_maps];

    long long *  elem_cmap_ids        = new long long [num_elem_comm_maps];
    long long *  elem_cmap_elem_cnts  = new long long [num_elem_comm_maps];


    if ( im_ne_get_cmap_params_l( id,
                                  node_cmap_ids,
                                  (long long*)node_cmap_node_cnts,
                                  elem_cmap_ids,
                                  (long long*)elem_cmap_elem_cnts,
                                  0/*not used proc_id*/ ) < 0 )++error;

    for(long long j = 0; j < num_node_comm_maps; j++) {
      comm_node_ids[j]       = new long long [node_cmap_node_cnts[j]];
      comm_node_proc_ids[j]  = new long long [node_cmap_node_cnts[j]];
      if ( im_ne_get_node_cmap_l( id,
                                  node_cmap_ids[j],
                                  comm_node_ids[j],
                                  comm_node_proc_ids[j],
                                  0/*not used proc_id*/ ) < 0 )++error;
      node_comm_proc_ids[j] = comm_node_proc_ids[j][0];
    }

    delete [] elem_cmap_ids;
    delete [] elem_cmap_elem_cnts;
  }

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Global Node Nums")));

  //Calculate global node ids
  long long * P1_globalNodeIds = new long long[numNodes];
  bool * P1_nodeIsOwned = new bool[numNodes];

  calc_global_node_ids(P1_globalNodeIds,
                       P1_nodeIsOwned,
                       numNodes,
                       num_node_comm_maps,
                       node_cmap_node_cnts,
                       node_comm_proc_ids,
                       comm_node_ids,
                       rank);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Boundary Conds")));

  // Container indicating whether a node is on the boundary (1-yes 0-no)
  FieldContainer<int> P1_nodeOnBoundary(numNodes);

  // Get boundary (side set) information
  long long * sideSetIds = new long long [numSideSets];
  long long numSidesInSet;
  long long numDFinSet;
  im_ex_get_side_set_ids_l(id,sideSetIds);
  for (int i=0; i<numSideSets; i++) {
    im_ex_get_side_set_param_l(id,sideSetIds[i],&numSidesInSet,&numDFinSet);
    if (numSidesInSet > 0){
      long long * sideSetElemList = new long long [numSidesInSet];
      long long * sideSetSideList = new long long [numSidesInSet];
      im_ex_get_side_set_l(id,sideSetIds[i],sideSetElemList,sideSetSideList);
      for (int j=0; j<numSidesInSet; j++) {

        int sideNode0 = P1_cellType.getNodeMap(1,sideSetSideList[j]-1,0);
        int sideNode1 = P1_cellType.getNodeMap(1,sideSetSideList[j]-1,1);

        P1_nodeOnBoundary(P1_elemToNode(sideSetElemList[j]-1,sideNode0))=1;
        P1_nodeOnBoundary(P1_elemToNode(sideSetElemList[j]-1,sideNode1))=1;
      }
      delete [] sideSetElemList;
      delete [] sideSetSideList;
    }
  }
  delete [] sideSetIds;


 // Enumerate edges 
  // NOTE: Only correct in serial
  FieldContainer<int> P1_elemToEdge(numElems,4);// Because quads
  FieldContainer<int> P1_elemToEdgeOrient(numElems,4);
  FieldContainer<double> P1_edgeCoord(1,dim);//will be resized  
  FieldContainer<int> P1_edgeToNode(1,2);//will be resized
  GenerateEdgeEnumeration(P1_elemToNode, P1_nodeCoord, P1_elemToEdge,P1_elemToEdgeOrient,P1_edgeToNode,P1_edgeCoord);
 

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Mesh Statistics")));

/**********************************************************************************/
/****************************** STATISTICS (Part I) *******************************/
/**********************************************************************************/
  // Statistics: Compute max / min of sigma parameter, mesh information
  // Put the current statistics into their own functions so this looks less awful
  double edge_length_max = 0;
  double edge_length_min = 0;
  double diag_length_max = 0;
  double diag_length_min = 0;
  double maxmin_ratio    = 0;
  double diag_ratio      = 0;
  double stretch         = 0;
  double taper           = 0;
  double skew            = 0;

  double dist1           = 0.0;
  double dist2           = 0.0;
  double dist3           = 0.0;
  double dist4           = 0.0;
 
  int edge = P1_elemToEdge(0, 0);
  int node1 = P1_edgeToNode(edge, 0);
  int node2 = P1_edgeToNode(edge, 1);
  double dist = 0;
  static const int NUM_NODE_PAIRS = 2;
  int diag_nodes1[] = {0, 1};
  int diag_nodes2[] = {2, 3};

  double x0 = 0.0;
  double x1 = 0.0;
  double x2 = 0.0;
  double x3 = 0.0;
  double y0 = 0.0;
  double y1 = 0.0;
  double y2 = 0.0;
  double y3 = 0.0;
  double e3 = 0.0;
  double f3 = 0.0;

  // Intialize
  for(int i=0; i<NUM_STATISTICS; i++) {
    local_stat_max[i] = 0.0;
    local_stat_min[i] = std::numeric_limits<double>::max();
    local_stat_sum[i] = 0.0;
  }

  for(int i=0; i<numElems; i++) {
    // 0 - Material property
    local_stat_max[0] = 1.0; 
    local_stat_min[0] = 1.0; 
    local_stat_sum[0] += 1.0;
    
    edge = P1_elemToEdge(i, 0);
    node1 = P1_edgeToNode(edge, 0);
    node2 = P1_edgeToNode(edge, 1);
    
    // 1 - Max/min edge - ratio of max to min edge length
    edge_length_max = distance2(P1_nodeCoord, node1, node2);
    edge_length_min = edge_length_max;
    for (int j=0; j<P1_numEdgesPerElem; j++) {
      edge = P1_elemToEdge(i,j);
      node1 = P1_edgeToNode(edge,0);
      node2 = P1_edgeToNode(edge,1);
      dist = distance2(P1_nodeCoord,node1,node2);
      edge_length_max = std::max(edge_length_max,dist);
      edge_length_min = std::min(edge_length_min,dist);
    }
    maxmin_ratio = edge_length_max / edge_length_min;
    local_stat_max[1] = std::max(local_stat_max[1],maxmin_ratio);
    local_stat_min[1] = std::min(local_stat_min[1],maxmin_ratio);
    local_stat_sum[1] += maxmin_ratio;
   
    // 2 - det of cell Jacobian (later)

    // 3 - Stretch
    diag_length_max = distance2(P1_nodeCoord, P1_elemToNode(i, diag_nodes1[0]),
                                           P1_elemToNode(i, diag_nodes2[0]));
    diag_length_min = distance2(P1_nodeCoord, P1_elemToNode(i, diag_nodes1[0]),
                                           P1_elemToNode(i, diag_nodes2[0]));
    for (int j=0; j<NUM_NODE_PAIRS; j++) {
      node1 = P1_elemToNode(i, diag_nodes1[j]);
      node2 = P1_elemToNode(i, diag_nodes2[j]);
      dist = distance2(P1_nodeCoord, node1, node2);
      diag_length_max = std::max(diag_length_max, dist);
      diag_length_min = std::min(diag_length_min, dist);
    }
    stretch = sqrt(3) * edge_length_min / diag_length_max;
    diag_ratio = diag_length_min / diag_length_max;
    local_stat_max[3] = std::max(local_stat_max[3], stretch);
    local_stat_min[3] = std::min(local_stat_min[3], stretch);
    local_stat_sum[3] += stretch;

    // 4 - Diagonal Ratio
    local_stat_max[4] = std::max(local_stat_max[4], diag_ratio);
    local_stat_min[4] = std::min(local_stat_min[4], diag_ratio);
    local_stat_sum[4] += diag_ratio;

    // 5 - Inverse Taper
    int opposite_edges1[2][2] = {{0, 1}, {0, 3}};
    int opposite_edges2[2][2] = {{2, 3}, {1, 2}};
    dist1 = distance2(P1_nodeCoord, P1_elemToNode(i, opposite_edges1[0][0]), P1_elemToNode(i, opposite_edges1[0][1]));
    dist2 = distance2(P1_nodeCoord, P1_elemToNode(i, opposite_edges2[0][0]), P1_elemToNode(i, opposite_edges2[0][1]));
    dist3 = distance2(P1_nodeCoord, P1_elemToNode(i, opposite_edges1[1][0]), P1_elemToNode(i, opposite_edges1[1][1]));
    dist4 = distance2(P1_nodeCoord, P1_elemToNode(i, opposite_edges2[1][0]), P1_elemToNode(i, opposite_edges2[1][1]));
    taper = std::max(dist1/dist2, dist2/dist1);
    taper = std::max(taper, dist3/dist4);
    taper = std::max(taper, dist4/dist3);
    local_stat_max[5] = std::max(local_stat_max[5], taper);
    local_stat_min[5] = std::min(local_stat_min[5], taper);
    local_stat_sum[5] += taper;


    // 6 - Skew
    // See "Quadrilateral and hexagonal shape parameters" - Robinson, J. 1994
    x0 = P1_nodeCoord(P1_elemToNode(i, 0), 0);
    x1 = P1_nodeCoord(P1_elemToNode(i, 1), 0);
    x2 = P1_nodeCoord(P1_elemToNode(i, 2), 0);
    x3 = P1_nodeCoord(P1_elemToNode(i, 3), 0);

    y0 = P1_nodeCoord(P1_elemToNode(i, 0), 1);
    y1 = P1_nodeCoord(P1_elemToNode(i, 1), 1);
    y2 = P1_nodeCoord(P1_elemToNode(i, 2), 1);
    y3 = P1_nodeCoord(P1_elemToNode(i, 3), 1);

    e3 = 0.25 * (x0-x1+x2-x3);
    f3 = 0.25 * (-y0-y1+y2+y3);
    skew = e3 / f3;
    local_stat_max[6] = std::max(local_stat_max[4], skew);
    local_stat_min[6] = std::min(local_stat_max[4], skew);
    local_stat_sum[6] += skew;


  }


/**********************************************************************************/
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting cubature")));

  


  // Generate higher order mesh
  // NOTE: Only correct in serial
  int Pn_numNodes = numNodes + (degree-1)*P1_edgeCoord.dimension(0) + (degree-1)*(degree-1)*numElems;
  int Pn_numNodesperElem = (degree+1)*(degree+1); // Quads
  FieldContainer<int>    elemToNode(numElems,Pn_numNodesperElem); 
  FieldContainer<double> nodeCoord(Pn_numNodes,dim);
  FieldContainer<int>   nodeOnBoundary(Pn_numNodes);
  std::vector<int> Pn_edgeNodes;
  std::vector<int> Pn_cellNodes;


  PromoteMesh_Pn_Kirby(degree,POINTTYPE_EQUISPACED,P1_elemToNode,P1_nodeCoord,P1_edgeCoord,P1_elemToEdge,P1_elemToEdgeOrient,
                       P1_nodeOnBoundary, elemToNode, nodeCoord, nodeOnBoundary, Pn_edgeNodes, Pn_cellNodes);


  long long numElems_aux = numElems*degree*degree;  //degree^2 P1 elements per Pn element in auxiliary mesh
  FieldContainer<int> aux_P1_elemToNode(numElems_aux,P1_numNodesPerElem); //4 P1 elements per Pn element
  CreateP1MeshFromPnMesh(degree, elemToNode, aux_P1_elemToNode);

  // Only works in serial
  std::vector<bool>Pn_nodeIsOwned(Pn_numNodes,true);
  std::vector<global_ordinal_type>Pn_globalNodeIds(Pn_numNodes);
  for(int i=0; i<Pn_numNodes; i++)
    Pn_globalNodeIds[i]=static_cast<global_ordinal_type>(i);

  std::vector<double> Pn_nodeCoordx(Pn_numNodes);
  std::vector<double> Pn_nodeCoordy(Pn_numNodes);
  for (int i=0; i<Pn_numNodes; i++) {
    Pn_nodeCoordx[i] = nodeCoord(i,0);
    Pn_nodeCoordy[i] = nodeCoord(i,1);
  }

  // Reset constants
  int P1_numNodes =numNodes;
  numNodes = Pn_numNodes;
  numNodesGlobal = numNodes;

  // Print mesh information
  if (MyPID == 0){
    std::cout << " Number of Pn Global Elements: " << numElemsGlobal << " \n";
    std::cout << " Number of Pn Global Nodes: " << numNodesGlobal << " \n";
    std::cout << " Number of faux P1 Global Elements: " << aux_P1_elemToNode.dimension(0) << " \n\n";
  }


  /**********************************************************************************/
  /********************************* GET CUBATURE ***********************************/
  /**********************************************************************************/

  // Get numerical integration points and weights
  DefaultCubatureFactory<double>  cubFactory;
  int cubDegree = 2*degree;
  RCP<Cubature<double> > myCub = cubFactory.create(Pn_cellType, cubDegree);

  int cubDim       = myCub->getDimension();
  int numCubPoints = myCub->getNumPoints();

  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);

  myCub->getCubature(cubPoints, cubWeights);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting cubature")));

  /**********************************************************************************/
  /*********************************** GET BASIS ************************************/
  /**********************************************************************************/

  // Define basis
  Basis_HGRAD_QUAD_Cn_FEM<double, FieldContainer<double> > myHGradBasis(degree,POINTTYPE_EQUISPACED);
  RCP<Basis<double,FieldContainer<double> > > myHGradBasis_rcp(&myHGradBasis, false);
  RCP<Basis_HGRAD_QUAD_Cn_FEM<double,FieldContainer<double> > > myHGradBasisWithDofCoords_rcp(&myHGradBasis, false);

  // Auxillary p=1 basis
  Basis_HGRAD_QUAD_C1_FEM<double, FieldContainer<double> > myHGradBasis_aux;
  RCP<Basis<double,FieldContainer<double> > > myHGradBasis_aux_rcp(&myHGradBasis_aux, false);


  int numFieldsG = myHGradBasis.getCardinality();

  FieldContainer<double> HGBValues(numFieldsG, numCubPoints);
  FieldContainer<double> HGBGrads(numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
  myHGradBasis.getValues(HGBValues, cubPoints, OPERATOR_VALUE);
  myHGradBasis.getValues(HGBGrads, cubPoints, OPERATOR_GRAD);


  //#define OUTPUT_REFERENCE_FUNCTIONS
#ifdef OUTPUT_REFERENCE_FUNCTIONS
  // Evaluate basis values and gradients at DOF points 
  FieldContainer<double> DofCoords(numFieldsG,spaceDim);
  myHGradBasis.getDofCoords(DofCoords);
  FieldContainer<double> HGBValues_at_Dofs(numFieldsG,numFieldsG);
  myHGradBasis.getValues(HGBValues_at_Dofs, DofCoords, OPERATOR_VALUE);


  printf("*** Dof Coords ***\n");
  for(int j=0; j<numFieldsG; j++)
    printf("[%d] %10.2e %10.2e\n",j,DofCoords(j,0),DofCoords(j,1));
  
  printf("*** CubPoints & Weights ***\n");
  for(int j=0; j<numCubPoints; j++)
    printf("[%d] %10.2e %10.2e (%10.2e)\n",j,cubPoints(j,0),cubPoints(j,1),cubWeights(j));
  
  printf("*** HGBValues @ Dofs ***\n");
  for(int i=0; i<numFieldsG; i++) {
    printf("[%d] ",i);
    for(int j=0; j<numFieldsG; j++) {
      printf("%10.2e ",zwrap(HGBValues_at_Dofs(i,j)));      
    }
    printf("\n");
  }
  
  printf("*** HGBValues ***\n");
  for(int i=0; i<numFieldsG; i++) {
    printf("[%d] ",i);
    for(int j=0; j<numCubPoints; j++) {
      printf("%10.2e ",zwrap(HGBValues(i,j)));
      
    }
    printf("\n");
  }
  
  printf("*** HGBGrad ***\n");
  for(int i=0; i<numFieldsG; i++) {
    printf("[%d] ",i);
    for(int j=0; j<numCubPoints; j++) {
      printf("(%10.2e %10.2e)",zwrap(HGBGrads(i,j,0)),zwrap(HGBGrads(i,j,1)));
      printf("\n");
    }
  }    

  printf("*** Pn_nodeIsOwned ***\n");
  for(int i=0; i<Pn_numNodes; i++)
    printf("%d ",(int)Pn_nodeIsOwned[i]);
  printf("\n");

  printf("*** Pn_nodeOnBoundary ***\n");
  for(int i=0; i<Pn_numNodes; i++)
    printf("%d ",(int)nodeOnBoundary[i]);
  printf("\n");

#endif

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Build global maps")));

  /**********************************************************************************/
  /********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
  /**********************************************************************************/
  // Count owned nodes (Pn)
  int Pn_ownedNodes=0;
  for(int i=0;i<numNodes;i++)
    if(Pn_nodeIsOwned[i]) Pn_ownedNodes++;

  // Build a list of the OWNED global ids...
  // NTS: will need to switch back to long long
  std::vector<global_ordinal_type> Pn_ownedGIDs(Pn_ownedNodes);
  int oidx=0;
  for(int i=0;i<numNodes;i++)
    if(Pn_nodeIsOwned[i]){
      Pn_ownedGIDs[oidx]=(global_ordinal_type)Pn_globalNodeIds[i];
      oidx++;
    }

  // Count owned nodes (P1)
  int P1_ownedNodes=0;
  for(int i=0;i<P1_numNodes;i++)
    if(P1_nodeIsOwned[i]) P1_ownedNodes++;
  std::vector<global_ordinal_type> P1_ownedGIDs(P1_ownedNodes);
  oidx=0;
  for(int i=0;i<P1_numNodes;i++)
    if(P1_nodeIsOwned[i]){
      P1_ownedGIDs[oidx]=(global_ordinal_type)P1_globalNodeIds[i];
      oidx++;
    }
  //TODO JJH 8-June-2016 need to populate edge and elem seed nodes

  //seed points for block relaxation
  ArrayRCP<int> nodeSeeds(Pn_numNodes,INVALID_LO);
  //unknowns at mesh nodes
  oidx=0;
  for(int i=0;i<P1_numNodes;i++) {
    if(P1_nodeIsOwned[i]){
      nodeSeeds[(int)P1_globalNodeIds[i]] = oidx;
      oidx++;
    }
  }
  int numNodeSeeds = oidx;

  //unknowns on edges
  ArrayRCP<int> edgeSeeds(Pn_numNodes,INVALID_LO);
  for (size_t i=0; i<Pn_edgeNodes.size(); ++i)
    edgeSeeds[Pn_edgeNodes[i]] = i;
  int numEdgeSeeds = Pn_edgeNodes.size();

  //unknowns in cell interiors
  ArrayRCP<int> cellSeeds(Pn_numNodes,INVALID_LO);
  for (size_t i=0; i<Pn_cellNodes.size(); ++i)
    cellSeeds[Pn_cellNodes[i]] = i;
  int numCellSeeds = Pn_cellNodes.size();
       
  // Generate map for nodes
  RCP<driver_map_type> globalMapG = rcp(new driver_map_type(INVALID_GO,&Pn_ownedGIDs[0],Pn_ownedNodes,0,Comm));
    
  // Generate p1 map
  RCP<driver_map_type> P1_globalMap = rcp(new driver_map_type(INVALID_GO,&P1_ownedGIDs[0],P1_ownedNodes,0,Comm));

  // Genetrate Pn-to-P1 coarsening.
  Kokkos::DynRankView<local_ordinal_type,typename NO::device_type>  elemToNodeI2;

  CopyFieldContainer2D(elemToNode,elemToNodeI2);
  
  if (inputSolverList.isParameter("aux P1") && inputSolverList.isParameter("linear P1"))
    throw std::runtime_error("Can only specify \"aux P1\" or \"linear P1\", not both.");
  if (inputSolverList.isParameter("linear P1")) {
    printf("Activating Linear scheduled p-coarsening...\n");
    Teuchos::ParameterList & mymuelu = inputSolverList.sublist("MueLu");
    Teuchos::ParameterList & level0  = mymuelu.sublist("level 0");
    level0.set("pcoarsen: element to node map",rcp(&elemToNodeI2,false));
    inputSolverList.remove("linear P1"); //even though LevelWrap happily accepts this parameter
  }



  // Global arrays in Tpetra format
  crs_matrix_type StiffMatrix(globalMapG, 20*numFieldsG);
  RCP<multivector_type> rhsVector = rcp(new multivector_type(globalMapG,1));
  RCP<multivector_type> femCoefficients = rcp(new multivector_type(globalMapG,1));

  tm = Teuchos::null;

  /**********************************************************************************/
  /**** COOORDINATES FOR DISTANCE LAPLACIAN                                       ***/
  /**********************************************************************************/

  // Put coordinates in multivector for output
  RCP<multivector_type> nCoord = rcp(new multivector_type(globalMapG,dim));

  int indOwned = 0;
  for (int inode=0; inode<Pn_numNodes; inode++) {
    if (Pn_nodeIsOwned[inode]) {
      nCoord->getDataNonConst(0)[indOwned]=nodeCoord(inode,0);
      nCoord->getDataNonConst(1)[indOwned]=nodeCoord(inode,1);
      indOwned++;
    }
  }

  /**********************************************************************************/
  /************************** DIRICHLET BC SETUP ************************************/
  /**********************************************************************************/
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Get Dirichlet boundary values")));

  int numBCNodes = 0;
  for (int inode = 0; inode < numNodes; inode++){
    if (nodeOnBoundary(inode) && Pn_nodeIsOwned[inode]){
      numBCNodes++;
    }
  }

  // Vector for use in applying BCs
  multivector_type v(globalMapG,true);
  ArrayRCP<scalar_type> vdata = v.getDataNonConst(0);

  // Set v to boundary values on Dirichlet nodes
  std::vector<int> BCNodes(numBCNodes);
  int indbc=0;
  int iOwned=0;
  for (int inode=0; inode<numNodes; inode++){
    if (Pn_nodeIsOwned[inode]){
      if (nodeOnBoundary(inode)){
        BCNodes[indbc]=iOwned;
        indbc++;
        double x  = nodeCoord(inode, 0);
        double y  = nodeCoord(inode, 1);
        vdata[iOwned]=exactSolution(x, y);
      }
      iOwned++;
    }
  }
    
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Global assembly")));
  
  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/

  // Define desired workset size and count how many worksets there are on this processor's mesh block
  int desiredWorksetSize = numElems;                      // change to desired workset size!
  int numWorksets        = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksets*desiredWorksetSize < numElems) numWorksets += 1;

  if (MyPID == 0) {
    std::cout << "Building discretization matrix and right hand side... \n\n";
    std::cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
    std::cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
  }

  // Create Pn matrix and RHS
  CreateLinearSystem(numWorksets,
                     desiredWorksetSize,
                     elemToNode,
                     nodeCoord,
                     cubPoints,
                     cubWeights,
                     myHGradBasis_rcp,
                     HGBGrads,
                     HGBValues,
                     Pn_globalNodeIds,
                     StiffMatrix,
                     rhsVector,
                     msg
                     );

  /**********************************************************************************/
  /********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
  /**********************************************************************************/

  printf("Example: StiffMatrix (ra,ro,co,do) = (%d,%d,%d,%d)\n",
	 (int)StiffMatrix.getRangeMap()->getGlobalNumElements(),
	 (int)StiffMatrix.getRowMap()->getGlobalNumElements(),
	 (int)StiffMatrix.getColMap()->getGlobalNumElements(),
	 (int)StiffMatrix.getRangeMap()->getGlobalNumElements());

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting cubature for auxiliary P1 mesh")));

  /////////////////////////////////////////////////////////////////////
  // Create P1 matrix and RHS to be used in preconditioner
  /////////////////////////////////////////////////////////////////////

  //Cubature

  // Get numerical integration points and weights
  int cubDegree_aux = 2; //TODO  This was 3 for P=2.  I think this should be 2 now .... Ask Chris.
  RCP<Cubature<double> > myCub_aux = cubFactory.create(P1_cellType, cubDegree_aux);

  int cubDim_aux       = myCub_aux->getDimension();
  int numCubPoints_aux = myCub_aux->getNumPoints();

  FieldContainer<double> cubPoints_aux(numCubPoints_aux, cubDim_aux);
  FieldContainer<double> cubWeights_aux(numCubPoints_aux);
  myCub_aux->getCubature(cubPoints_aux, cubWeights_aux);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting basis for auxiliary P1 mesh")));

  //Basis

  // Define basis
  int numFieldsG_aux = myHGradBasis_aux.getCardinality();
  FieldContainer<double> HGBValues_aux(numFieldsG_aux, numCubPoints_aux);
  FieldContainer<double> HGBGrads_aux(numFieldsG_aux, numCubPoints_aux, spaceDim);

  // Evaluate basis values and gradients at cubature points
  myHGradBasis_aux.getValues(HGBValues_aux, cubPoints_aux, OPERATOR_VALUE);
  myHGradBasis_aux.getValues(HGBGrads_aux, cubPoints_aux, OPERATOR_GRAD);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Global assembly (auxiliary system)")));

  crs_matrix_type StiffMatrix_aux(globalMapG, 20*numFieldsG_aux);
  RCP<multivector_type> rhsVector_aux = rcp(new multivector_type(globalMapG,1));

  desiredWorksetSize = numElems_aux;
  numWorksets        = numElems_aux/desiredWorksetSize;
  if(numWorksets*desiredWorksetSize < numElems_aux) numWorksets += 1;
  CreateLinearSystem(numWorksets,
                     desiredWorksetSize,
                     aux_P1_elemToNode,
                     nodeCoord,
                     cubPoints_aux,
                     cubWeights_aux,
                     myHGradBasis_aux_rcp,
                     HGBGrads_aux,
                     HGBValues_aux,
                     Pn_globalNodeIds,
                     StiffMatrix_aux,
                     rhsVector_aux,
                     msg
                     );


/**********************************************************************************/
/***************************** STATISTICS (Part III) ******************************/
/**********************************************************************************/
  std::vector<double> global_stat_max(NUM_STATISTICS),global_stat_min(NUM_STATISTICS), global_stat_sum(NUM_STATISTICS);
  Teuchos::reduceAll(*Comm,Teuchos::REDUCE_MIN,NUM_STATISTICS,local_stat_min.data(),global_stat_min.data());
  Teuchos::reduceAll(*Comm,Teuchos::REDUCE_MAX,NUM_STATISTICS,local_stat_max.data(),global_stat_max.data());
  Teuchos::reduceAll(*Comm,Teuchos::REDUCE_SUM,NUM_STATISTICS,local_stat_sum.data(),global_stat_sum.data());
  // NOTE: All output properties should be unitless if we want to compare across problems.
  // NOTE: Should the mean be weighted by cell volume?  That is not currently done.

  // 0 - Material property
  problemStatistics.set("sigma: min/mean",global_stat_min[0]/global_stat_sum[0]*numElemsGlobal);
  problemStatistics.set("sigma: max/mean",global_stat_max[0]/global_stat_sum[0]*numElemsGlobal);

  // 1 - Max/min edge ratio
  problemStatistics.set("element edge ratio: min",global_stat_min[1]);
  problemStatistics.set("element edge ratio: max",global_stat_max[1]);
  problemStatistics.set("element edge ratio: mean",global_stat_sum[1] / numElemsGlobal);

  // 2 - det of cell Jacobian (later)
  problemStatistics.set("element det jacobian: min/mean",global_stat_min[2]/global_stat_sum[2]*numElemsGlobal);
  problemStatistics.set("element det jacobian: max/mean",global_stat_max[2]/global_stat_sum[2]*numElemsGlobal);

  if(use_new_problem_stats) {
    // 3 - Stretch
    problemStatistics.set("Stretch max", global_stat_max[3]);
    problemStatistics.set("Stretch min", global_stat_min[3]);
    problemStatistics.set("Stretch mean", global_stat_sum[3] / numElemsGlobal);

    // 4 - Diagonal Ratio
    problemStatistics.set("Diagonal Ratio max", global_stat_max[4]);
    problemStatistics.set("Diagonal Ratio min", global_stat_min[4]);
    problemStatistics.set("Diagonal Ratio mean", global_stat_sum[4]/numElemsGlobal);

    // 5 - Inverse Taper
    problemStatistics.set("Inverse Taper max", global_stat_max[5]);
    problemStatistics.set("Inverse Taper min", global_stat_min[5]);
    problemStatistics.set("Inverse Taper mean", global_stat_sum[5] / numElemsGlobal);
    
    // 6 - Skew
    problemStatistics.set("Skew max", global_stat_max[6]);
    problemStatistics.set("Skew min", global_stat_min[6]);
    problemStatistics.set("Skew mean", global_stat_sum[6] / numElemsGlobal);

    // 7 - Lapl Diag
    problemStatistics.set("Lapl Diag max", global_stat_max[7]);
    problemStatistics.set("Lapl Diag min", global_stat_min[7]);
    problemStatistics.set("Lapl Diag mean", global_stat_sum[7] / StiffMatrix.getGlobalNumEntries());
  }

  // Print Problem Statistics
  std::cout<<"*** Problem Statistics ***\n"<<problemStatistics<<std::endl;


  /**********************************************************************************/
  /********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
  /**********************************************************************************/
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs")));

  // Generate Pn-to-P1 identity coarsening (base mesh to auxiliary mesh).
  RCP<crs_matrix_type> P_identity, R_identity;

  if (inputSolverList.isParameter("aux P1")) {
    printf("Generating Identity Pn-to-P1 coarsening...\n");
    GenerateIdentityCoarsening_pn_to_p1(elemToNode, StiffMatrix_aux.getDomainMap(), StiffMatrix.getRangeMap(), P_identity, R_identity);
    inputSolverList.remove("aux P1"); //even though LevelWrap happily accepts this parameter
  }
  
/**********************************************************************************/
/******************************* ADJUST MATRIX DUE TO BC **************************/
/**********************************************************************************/

  // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  Apply_Dirichlet_BCs(BCNodes,StiffMatrix,*femCoefficients,*rhsVector,v);

  ///////////////////////////////////////////
  // Apply BCs to auxiliary P1 matrix
  ///////////////////////////////////////////
  // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  std::cout << "numBCNodes = " << numBCNodes << std::endl;
  std::cout << "globalMapG #elts = " << globalMapG->getNodeNumElements() << std::endl;
  Apply_Dirichlet_BCs(BCNodes,StiffMatrix_aux,*rhsVector_aux,*rhsVector_aux,v);

  tm = Teuchos::null;

  // Optionally dump the matrix and/or its coords to files.
  {
    typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
    if (matrixFilename != "") {
      writer_type::writeSparseFile (matrixFilename, StiffMatrix);
    }
    if (rhsFilename != "") {
      writer_type::writeDenseFile (rhsFilename, rhsVector);
    }
    if (coordsFilename != "") {
      writer_type::writeDenseFile (coordsFilename, nCoord);
    }
  }


  /**********************************************************************************/
  /*********************************** SOLVE ****************************************/
  /**********************************************************************************/
  // Which Solver?
  std::string solveType("cg");
  if(inputSolverList.isParameter("solver")) {
    solveType = inputSolverList.get<std::string>("solver");
  }


  // Run the solver
  std::string amgType("MueLu");

  ParameterList amgList;
  std::string seedType = inputSolverList.get("seed type","node");
  inputSolverList.remove("seed type");
  if (inputSolverList.isSublist("MueLu"))
    amgList = inputSolverList.sublist("MueLu");
  else
    amgList = inputSolverList;
  std::string lev0List = "level 0";
  if (amgList.isSublist(lev0List)) {
    std::cout << "found \"" << lev0List << "\" sublist" << std::endl;
    ParameterList &sl = amgList.sublist(lev0List);
    std::string smooType = sl.get<std::string>("smoother: type");
    if ( (smooType == "SPARSE BLOCK RELAXATION" || smooType == "BLOCK RELAXATION") && sl.isParameter("smoother: params")) {
      ParameterList &ssl = sl.sublist("smoother: params");
      std::cout << "found \"smoother: params\" for block relaxation" << std::endl;
      int numLocalParts = -1;
      std::cout << "setting \"partitioner: map\"" << std::endl;
      if (seedType == "node") {
        ssl.set("partitioner: map", nodeSeeds);
        numLocalParts = numNodeSeeds;
      } else if (seedType == "edge") {
        ssl.set("partitioner: map", edgeSeeds);
        numLocalParts = numEdgeSeeds;
      } else if (seedType == "cell") {
        ssl.set("partitioner: map", cellSeeds);
        numLocalParts = numCellSeeds;
      }
      if(numLocalParts!=-1){
	std::cout << "setting \"partitioner: local parts\" = " << numLocalParts << std::endl;
	ssl.set("partitioner: local parts", numLocalParts);
      }
    }
    if(sl.isParameter("coarse: type")) {
	std::string coarseType = sl.get<std::string>("coarse: type");
	if ((coarseType == "SPARSE BLOCK RELAXATION"|| coarseType == "BLOCK RELAXATION") && sl.isParameter("coarse: params")) {
	  ParameterList &ssl = sl.sublist("coarse: params");
	  std::cout << "found \"smoother: params\" for block relaxation" << std::endl;
	  int numLocalParts=-1;
	  std::cout << "setting \"partitioner: map\"" << std::endl;
	  if (seedType == "node") {
	    ssl.set("partitioner: map", nodeSeeds);
	    numLocalParts = numNodeSeeds;
	  } else if (seedType == "edge") {
	    ssl.set("partitioner: map", edgeSeeds);
	    numLocalParts = numEdgeSeeds;
	  } else if (seedType == "cell") {
	    ssl.set("partitioner: map", cellSeeds);
	    numLocalParts = numCellSeeds;
	  }
	  if(numLocalParts!=-1){
	    std::cout << "setting \"partitioner: local parts\" = " << numLocalParts << std::endl;
	    ssl.set("partitioner: local parts", numLocalParts);
	  }
	}
    }

  }

  // /////////////////////////////////////////////////////////////////////// //

  // Shove coordinates in multivector for saving to file.
  multivector_type coordinates(StiffMatrix.getRowMap(),2);
  ArrayRCP<scalar_type> xdat = coordinates.getDataNonConst(0);
  ArrayRCP<scalar_type> ydat = coordinates.getDataNonConst(1);
  for (size_t i=0; i<Pn_nodeCoordx.size(); ++i) {
    xdat[i] = Pn_nodeCoordx[i];
    ydat[i] = Pn_nodeCoordy[i];
  }
  xdat = Teuchos::null; ydat = Teuchos::null;
  //writer_type::writeDenseFile("coords.m", coordinates);

  multivector_type P1coordinates(P1_globalMap,2);
  xdat = P1coordinates.getDataNonConst(0);
  ydat = P1coordinates.getDataNonConst(1);
  for (size_t i=0; i<P1coordinates.getLocalLength(); i++) {
    xdat[i] = nodeCoordx[i];
    ydat[i] = nodeCoordy[i];
  }
  xdat = Teuchos::null; ydat = Teuchos::null;
  //writer_type::writeDenseFile("P1coords.m", P1coordinates);

  multivector_type nBound(globalMapG,true);

  RCP<multivector_type> exactNodalVals = rcp(new multivector_type(globalMapG,1));
  double TotalErrorResidual = 0.0;
  double TotalErrorExactSol = 0.0;

  // Get exact solution at nodes
  for (int i = 0; i<numNodes; i++) {
    if (Pn_nodeIsOwned[i]){
      double x = nodeCoord(i,0);
      double y = nodeCoord(i,1);
      double exactu = exactSolution(x, y);

      int rowindex=Pn_globalNodeIds[i];
      exactNodalVals->sumIntoGlobalValue(rowindex, 0, exactu);
    }
  }

  char probType[10] = "laplace";

  RCP<crs_matrix_type> interpolationMatrix, restrictionMatrix;
  if (P_identity != Teuchos::null) {
    Teuchos::ParameterList & level1 = amgList.sublist("level 1");
    RCP<xpetra_crs_matrix_type> xA1 = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>(rcpFromRef(StiffMatrix_aux));
    RCP<xpetra_crs_matrix_type> xP = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>(P_identity);
    level1.set("A",xA1);
    level1.set("P",xP);
    amgList.set("transpose: use implicit",true);
    // create nullspace for Level 1 (first coarse AMG level)
    RCP<multivector_type> nullspace = rcp(new multivector_type(exactNodalVals->getMap(),1));
    ArrayRCP<scalar_type> data = nullspace->getDataNonConst(0);
    for (int i=0; i<data.size(); ++i)
      data[i] = 1.0;
    RCP<xpetra_multivector_type> xnullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector<scalar_type,local_ordinal_type,global_ordinal_type,NO>(nullspace);
    level1.set("Nullspace",xnullspace);
  }

  int maxits = inputSolverList.get("Maximum Iterations",(int)100);
  double tol = inputSolverList.get("Convergence Tolerance",(double)1e-10);

  if (amgList.isParameter("solve Ae=0")) {
    rhsVector->scale(0.0);
    femCoefficients->randomize();
  }


  TestMultiLevelPreconditionerLaplace(probType, amgList,
                                      rcpFromRef(StiffMatrix),
				      nCoord,
				      exactNodalVals,
                                      rhsVector,            
				      maxits,
				      tol,
				      femCoefficients,
                                      TotalErrorResidual,   TotalErrorExactSol,
                                      amgType,
				      solveType);

  /**********************************************************************************/
  /**************************** CALCULATE ERROR *************************************/
  /**********************************************************************************/

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Calculate error")));

  double L2err = 0.0;
  double L2errTot = 0.0;
  double H1err = 0.0;
  double H1errTot = 0.0;
  double Linferr = 0.0;
  double LinferrTot = 0.0;

  // Import solution onto current processor
  //int numNodesGlobal = globalMapG.NumGlobalElements();
  RCP<driver_map_type>  solnMap = rcp(new driver_map_type(static_cast<Tpetra::global_size_t>(numNodesGlobal), static_cast<size_t>(numNodesGlobal), 0, Comm));
  Tpetra::Import<local_ordinal_type, global_ordinal_type, NO> solnImporter(globalMapG,solnMap);
  multivector_type  uCoeff(solnMap,1);
  uCoeff.doImport(*femCoefficients, solnImporter, Tpetra::INSERT);

  // Define desired workset size
  desiredWorksetSize = numElems;
  int numWorksetsErr    = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksetsErr*desiredWorksetSize < numElems) numWorksetsErr += 1;

  // Get cubature points and weights for error calc (may be different from previous)
  Intrepid::DefaultCubatureFactory<double>  cubFactoryErr;
  int cubDegErr = 3*degree;
  RCP<Intrepid::Cubature<double> > cellCubatureErr = cubFactoryErr.create(Pn_cellType, cubDegErr);
  int cubDimErr       = cellCubatureErr->getDimension();
  int numCubPointsErr = cellCubatureErr->getNumPoints();
  Intrepid::FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
  Intrepid::FieldContainer<double> cubWeightsErr(numCubPointsErr);
  cellCubatureErr->getCubature(cubPointsErr, cubWeightsErr);

  // Evaluate basis values and gradients at cubature points
  Intrepid::FieldContainer<double> uhGVals(numFieldsG, numCubPointsErr);
  Intrepid::FieldContainer<double> uhGrads(numFieldsG, numCubPointsErr, spaceDim);
  myHGradBasis.getValues(uhGVals, cubPointsErr, Intrepid::OPERATOR_VALUE);
  myHGradBasis.getValues(uhGrads, cubPointsErr, Intrepid::OPERATOR_GRAD);

  // Loop over worksets
  for(int workset = 0; workset < numWorksetsErr; workset++){

    // compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // when numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
    worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // now we know the actual workset size and can allocate the array for the cell nodes
    worksetSize  = worksetEnd - worksetBegin;
    Intrepid::FieldContainer<double> cellWorksetEr(worksetSize, numFieldsG, spaceDim);
    Intrepid::FieldContainer<double> worksetApproxSolnCoef(worksetSize, numFieldsG);

    // loop over cells to fill arrays with coordinates and discrete solution coefficient
    int cellCounter = 0;
    ArrayRCP<const scalar_type> uCoeffData = uCoeff.getData(0);
    ArrayRCP<const scalar_type> femCoeffData = femCoefficients->getData(0);
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      for (int node = 0; node < numFieldsG; node++) {
        cellWorksetEr(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
        cellWorksetEr(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);

        int rowIndex  = Pn_globalNodeIds[elemToNode(cell, node)];
#ifdef HAVE_MPI
        worksetApproxSolnCoef(cellCounter, node) = uCoeffData[rowIndex];
#else
        worksetApproxSolnCoef(cellCounter, node) = femCoeffData[rowIndex];
#endif
      }

      cellCounter++;

    } // end cell loop

      // Containers for Jacobian
    Intrepid::FieldContainer<double> worksetJacobianE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
    Intrepid::FieldContainer<double> worksetJacobInvE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
    Intrepid::FieldContainer<double> worksetJacobDetE(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetCubWeightsE(worksetSize, numCubPointsErr);

    // Containers for basis values and gradients in physical space
    Intrepid::FieldContainer<double> uhGValsTrans(worksetSize,numFieldsG, numCubPointsErr);
    Intrepid::FieldContainer<double> uhGradsTrans(worksetSize, numFieldsG, numCubPointsErr, spaceDim);

    // compute cell Jacobians, their inverses and their determinants
    IntrepidCTools::setJacobian(worksetJacobianE, cubPointsErr, cellWorksetEr,  myHGradBasis_rcp);
    IntrepidCTools::setJacobianInv(worksetJacobInvE, worksetJacobianE );
    IntrepidCTools::setJacobianDet(worksetJacobDetE, worksetJacobianE );

    // map cubature points to physical frame
    Intrepid::FieldContainer<double> worksetCubPoints(worksetSize, numCubPointsErr, cubDimErr);
    IntrepidCTools::mapToPhysicalFrame(worksetCubPoints, cubPointsErr, cellWorksetEr, myHGradBasis_rcp);

    // evaluate exact solution and gradient at cubature points
    Intrepid::FieldContainer<double> worksetExactSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetExactSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    evaluateExactSolution(worksetExactSoln, worksetCubPoints);
    evaluateExactSolutionGrad(worksetExactSolnGrad, worksetCubPoints);

    // transform basis values to physical coordinates
    IntrepidFSTools::HGRADtransformVALUE<double>(uhGValsTrans, uhGVals);
    IntrepidFSTools::HGRADtransformGRAD<double>(uhGradsTrans, worksetJacobInvE, uhGrads);

    // compute weighted measure
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeightsE, worksetJacobDetE, cubWeightsErr);

    // evaluate the approximate solution and gradient at cubature points
    Intrepid::FieldContainer<double> worksetApproxSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetApproxSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    IntrepidFSTools::evaluate<double>(worksetApproxSoln, worksetApproxSolnCoef, uhGValsTrans);
    IntrepidFSTools::evaluate<double>(worksetApproxSolnGrad, worksetApproxSolnCoef, uhGradsTrans);

    // get difference between approximate and exact solutions
    Intrepid::FieldContainer<double> worksetDeltaSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetDeltaSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    IntrepidRSTools::subtract(worksetDeltaSoln, worksetApproxSoln, worksetExactSoln);
    IntrepidRSTools::subtract(worksetDeltaSolnGrad, worksetApproxSolnGrad, worksetExactSolnGrad);

    // take absolute values
    IntrepidRSTools::absval(worksetDeltaSoln);
    IntrepidRSTools::absval(worksetDeltaSolnGrad);
    // apply cubature weights to differences in values and grads for use in integration
    Intrepid::FieldContainer<double> worksetDeltaSolnWeighted(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetDeltaSolnGradWeighted(worksetSize, numCubPointsErr, spaceDim);
    IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnWeighted,
                                                    worksetCubWeightsE, worksetDeltaSoln);
    IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnGradWeighted,
                                                    worksetCubWeightsE, worksetDeltaSolnGrad);

    // integrate to get errors on each element
    Intrepid::FieldContainer<double> worksetL2err(worksetSize);
    Intrepid::FieldContainer<double> worksetH1err(worksetSize);
    IntrepidFSTools::integrate<double>(worksetL2err, worksetDeltaSoln,
                                       worksetDeltaSolnWeighted, Intrepid::COMP_BLAS);
    IntrepidFSTools::integrate<double>(worksetH1err, worksetDeltaSolnGrad,
                                       worksetDeltaSolnGradWeighted, Intrepid::COMP_BLAS);

    // loop over cells to get errors for total workset
    cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // loop over cubature points
      for(int nPt = 0; nPt < numCubPointsErr; nPt++){

        Linferr = std::max(Linferr, worksetDeltaSoln(cellCounter,nPt));

      }

      L2err += worksetL2err(cellCounter);
      H1err += worksetH1err(cellCounter);

      cellCounter++;

    } // end cell loop

  } // end loop over worksets

#ifdef HAVE_MPI
  // sum over all processors
  TC_sumAll(Comm,L2err,L2errTot);
  TC_sumAll(Comm,H1err,H1errTot);
  TC_maxAll(Comm,Linferr,LinferrTot);
#else
  L2errTot = L2err;
  H1errTot = H1err;
  LinferrTot = Linferr;
#endif


  if (MyPID == 0) {
    std::cout << "\n" << "L2 Error:  " << sqrt(L2errTot) <<"\n";
    std::cout << "H1 Error:  " << sqrt(H1errTot) <<"\n";
    std::cout << "LInf Error:  " << LinferrTot <<"\n\n";
  }

  const bool alwaysWriteLocal = false;
  const bool writeGlobalStats = true;
  const bool writeZeroTimers  = false;
  const bool ignoreZeroTimers = true;
  const std::string filter    = "";
  if (optPrintTimings) {
    tm.reset();
    TimeMonitor::summarize(Comm.ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                           writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
  }

  // Cleanup
  for(long long b = 0; b < numElemBlk; b++){
    delete [] elmt_node_linkage[b];
    delete [] element_types[b];
  }
  delete [] block_ids;
  delete [] nodes_per_element;
  delete [] element_attributes;
  delete [] element_types;
  delete [] elmt_node_linkage;
  delete [] elements;
  delete [] P1_globalNodeIds;
  delete [] P1_nodeIsOwned;
  if(num_node_comm_maps > 0){
    delete [] node_comm_proc_ids;
    delete [] node_cmap_node_cnts;
    delete [] node_cmap_ids;
    for(long long i=0;i<num_node_comm_maps;i++){
      delete [] comm_node_ids[i];
      delete [] comm_node_proc_ids[i];
    }

    delete [] comm_node_ids;
    delete [] comm_node_proc_ids;
  }

  delete [] nodeCoordx;
  delete [] nodeCoordy;

  // delete mesh
  Delete_Pamgen_Mesh();

  return 0;

}
/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/

/**********************************************************************************/
/************ USER DEFINED FUNCTIONS FOR EXACT SOLUTION ***************************/
/**********************************************************************************/

template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y) {

  // Patch test: bi-linear function is in the FE space and should be recovered
    return 1. + x + y + x*y;

  // Analytic solution with homogeneous Dirichlet boundary data
  //  return sin(M_PI*x)*sin(M_PI*y)**exp(x+y);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  //  return exp(x + y )/(1. + x*y);
}


template<typename Scalar>
void materialTensor(Scalar material[][2], const Scalar& x, const Scalar& y) {

  material[0][0] = 1.;
  material[0][1] = 0.;
  //
  material[1][0] = 0.;
  material[1][1] = 1.;
}

/**********************************************************************************/
/************** AUXILIARY FUNCTIONS FROM EXACT SOLUTION ***************************/
/**********************************************************************************/

/************ Grad of Exact Solution ****************/
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[2], const Scalar& x, const Scalar& y) {

  // To enable derivatives of the gradient (i.e., 2nd derivatives of the exact solution) need 2 levels of fad types
  Sacado::Fad::SFad<Scalar,2> fad_x = x;
  Sacado::Fad::SFad<Scalar,2> fad_y = y;
  Sacado::Fad::SFad<Scalar,2> u;

  // Indicate the independent variables
  fad_x.diff(0,2);
  fad_y.diff(1,2);

  u = exactSolution(fad_x, fad_y);

  gradExact[0] = u.dx(0);
  gradExact[1] = u.dx(1);
}

/************ Source Term (RHS) ****************/
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y){

  Scalar u;
  Scalar grad_u[2];
  Scalar flux[2] = {0.0, 0.0};
  Scalar material[2][2];
  Scalar f = 0.;

  // Indicate the independent variables
  x.diff(0,2);
  y.diff(1,2);

  // Get exact solution and its gradient
  u = exactSolution(x, y);
  exactSolutionGrad(grad_u, x, y);

  // Get material tensor
  materialTensor<Scalar>(material, x, y);

  // Compute total flux = (A.grad u)
  for(int i = 0; i < 2; i++){

    // Add diffusive flux
    for(int j = 0; j < 2; j++){
      flux[i] += material[i][j]*grad_u[j];
    }
  }

  // Compute source term (right hand side): f = -div(A.grad u)
  f = -(flux[0].dx(0) + flux[1].dx(1));

  return f;
}

/**********************************************************************************/
/*************************** EVALUATION METHODS ***********************************/
/**********************************************************************************/

/************ Material Tensor ****************/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        matTensorValues,
                            const ArrayIn &   evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints        = evaluationPoints.dimension(1);
  int spaceDim         = evaluationPoints.dimension(2);

  double material[2][2];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);

      materialTensor<double>(material, x, y);

      for(int row = 0; row < spaceDim; row++){
        for(int col = 0; col < spaceDim; col++){
          matTensorValues(cell, pt, row, col) = material[row][col];
        }
      }
    }
  }
}
/************ Source Term (RHS) ****************/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      Sacado::Fad::SFad<double,2> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<double,2> y = evaluationPoints(cell, pt, 1);

      sourceTermValues(cell, pt) = sourceTerm<Sacado::Fad::SFad<double,2> >(x, y).val();
    }
  }
}

/************ Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);

      exactSolutionValues(cell, pt) = exactSolution<double>(x, y);
    }
  }
}
/************ Grad of Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);
  int spaceDim  = evaluationPoints.dimension(2);

  double gradient[2];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);

      exactSolutionGrad<double>(gradient, x, y);

      for(int row = 0; row < spaceDim; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}

/**********************************************************************************/
/******************************* TEST ML ******************************************/
/**********************************************************************************/

// Test ML
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
                                 ParameterList   & amgList,
                                 RCP<crs_matrix_type> const &A0,
                                 RCP<multivector_type> & nCoord,
                                 RCP<multivector_type> const & xexact,
                                 RCP<multivector_type> & b,
				 int maxIts,
				 double tol,	
                                 RCP<multivector_type> & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol,
				 std::string &amgType,
				 std::string &solveType)
{

  //  int mypid = A0->getComm()->getRank();
  RCP<multivector_type>x = uh;

  linear_problem_type Problem(linear_problem_type(A0, x, b));
  RCP<multivector_type> lhs = Problem.getLHS();
  RCP<const multivector_type> rhs = Problem.getRHS();

  RCP<TimeMonitor>
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Total solve time")));

  // =================== //
  // call ML and Belos   //
  // =================== //

  if (amgType != "ML" && amgType != "MueLu") {
    std::ostringstream errStr;
    errStr << "Error: \"amgType\" has value \"" << amgType << "\". Valid values are \"ML\" or \"MueLu\".";
    throw std::runtime_error(errStr.str());
  }

  int numIterations=0;
  if (amgType == "ML") {
    throw std::runtime_error("Error: ML does not support Tpetra objects");
  } else if (amgType == "MueLu") {

    // Multigrid Hierarchy, the easy way  
    RCP<operator_type> A0op = A0;
    amgList.sublist("user data").set("Coordinates",nCoord);
    Teuchos::RCP<muelu_tpetra_operator> M = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,NO>(A0op, amgList);
    Problem.setRightPrec(M);

    bool set = Problem.setProblem();
    if (set == false) {
      std::cout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    ParameterList belosList;
    belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
    belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency",      1);
    belosList.set("Output Style",          Belos::Brief);
    bool scaleResidualHist = true;
    if (!scaleResidualHist)
      belosList.set("Implicit Residual Scaling", "None");

    // Create an iterative solver manager
    RCP< Belos::SolverManager<scalar_type, multivector_type, operator_type> > solver;
    if(solveType == "cg")
      solver = rcp(new Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type>(rcpFromRef(Problem), rcp(&belosList, false)));
    else if (solveType == "gmres") { 
      solver = rcp(new Belos::PseudoBlockGmresSolMgr<scalar_type, multivector_type, operator_type>(rcpFromRef(Problem), rcp(&belosList, false)));
    }
    else if (solveType == "fixed point" || solveType == "fixed-point") {
      solver = rcp(new Belos::FixedPointSolMgr<scalar_type, multivector_type, operator_type>(rcpFromRef(Problem), rcp(&belosList, false)));   
    }
    else {
      std::cout << "\nERROR:  Invalid solver '"<<solveType<<"'" << std::endl;
      return EXIT_FAILURE;
    }


    // Perform solve
    solver->solve();
    //writer_type::writeDenseFile("sol.m", x);
    numIterations = solver->getNumIters();
  }

  ArrayRCP<const scalar_type> lhsdata = lhs->getData(0);
  ArrayRCP<const scalar_type> xexactdata = xexact->getData(0);

  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  double d = 0.0, d_tot = 0.0 , s =0.0, s_tot=0.0;
  for( size_t i=0 ; i<lhs->getMap()->getNodeNumElements() ; ++i ) {
    d += (lhsdata[i] - xexactdata[i]) * (lhsdata[i] - xexactdata[i]);
    s +=  xexactdata[i]* xexactdata[i];
  }
  lhsdata = Teuchos::null;
  xexactdata = Teuchos::null;

  TC_sumAll(A0->getComm(),d,d_tot);
  TC_sumAll(A0->getComm(),s,s_tot);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  Array<typename ScalarTraits::magnitudeType> Norm(1);
  vector_type Ax(rhs->getMap());
  A0->apply(*lhs, Ax);
  Ax.update(1.0, *rhs, -1.0);
  Ax.norm2(Norm);

  string msg = ProblemType;

  int numProc = A0->getComm()->getSize();
  if (A0->getComm()->getRank() == 0) {
    cout << msg << endl << "......Using " << numProc << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm[0] << endl;
    cout << msg << "......||x_exact - x||_2/||x_exact||_2 = " << sqrt(d_tot/s_tot) << endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm[0];

  return(numIterations);

}

/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


void GenerateEdgeEnumeration(const FieldContainer<int> & elemToNode, const FieldContainer<double> & nodeCoord, FieldContainer<int> & elemToEdge, FieldContainer<int> & elemToEdgeOrient, FieldContainer<int> & edgeToNode, FieldContainer<double> & edgeCoord) {
  // Not especially efficient, but effective... at least in serial!
  
  int numElems        = elemToNode.dimension(0);
  int numNodesperElem = elemToNode.dimension(1);
  int dim             = nodeCoord.dimension(1);

  // Sanity checks
  if(numNodesperElem !=4) throw std::runtime_error("Error: GenerateEdgeEnumeration only works on Quads!");
  if(elemToEdge.dimension(0)!=numElems || elemToEdge.dimension(1)!=4 || elemToEdge.dimension(0)!=elemToEdgeOrient.dimension(0) || elemToEdge.dimension(1)!=elemToEdgeOrient.dimension(1)) 
    throw std::runtime_error("Error: GenerateEdgeEnumeration array size mismatch");

  int edge_node0_id[4]={0,1,2,3};
  int edge_node1_id[4]={1,2,3,0};
  
  // Run over all the elements and start enumerating edges
  typedef std::map<std::pair<int,int>,int> en_map_type;
  en_map_type e2n;

  int num_edges=0;
  en_map_type edge_map;
  for(int i=0; i<numElems; i++) {

    // Generate edge pairs, based on the global orientation of "edges go from low ID to high ID"
    for(int j=0; j<4; j++) {
      int lo = std::min(elemToNode(i,edge_node0_id[j]),elemToNode(i,edge_node1_id[j]));
      int hi = std::max(elemToNode(i,edge_node0_id[j]),elemToNode(i,edge_node1_id[j]));

      std::pair<int,int> ep(lo,hi);

      int edge_id;
      en_map_type::iterator iter = edge_map.find(ep);      
      if(iter==edge_map.end()) {
        edge_map[ep] = num_edges;
        edge_id = num_edges;
        num_edges++;
      }
      else
        edge_id = (*iter).second;
      
      elemToEdge(i,j) = edge_id;
      elemToEdgeOrient(i,j) = (lo==elemToNode(i,edge_node0_id[j]))?1:-1;
    }
  }      

  // Fill out the edge centers (clobbering data if needed)
  edgeCoord.resize(num_edges,dim);
  for(int i=0; i<numElems; i++) {
    for(int j=0; j<4; j++) {      
      int n0 = elemToNode(i,edge_node0_id[j]);
      int n1 = elemToNode(i,edge_node1_id[j]);
      for(int k=0; k<dim; k++)
        edgeCoord(elemToEdge(i,j),k) = (nodeCoord(n0,k)+nodeCoord(n1,k))/2.0;
    }
  }

  // Edge to Node connectivity
  edgeToNode.resize(num_edges,2);
  for(int i=0; i<numElems; i++) {
    for(int j=0; j<4; j++) {      
      int lo = std::min(elemToNode(i,edge_node0_id[j]),elemToNode(i,edge_node1_id[j]));
      int hi = std::max(elemToNode(i,edge_node0_id[j]),elemToNode(i,edge_node1_id[j]));
      edgeToNode(elemToEdge(i,j),0) = lo;
      edgeToNode(elemToEdge(i,j),1) = hi;
    }
  }

  
  //#define DEBUG_EDGE_ENUMERATION
#ifdef DEBUG_EDGE_ENUMERATION
  printf("**** Edge coordinates ***\n");
  for(int i=0; i<num_edges; i++)
    printf("[%2d] %10.2f %10.2f\n",i,edgeCoord(i,0),edgeCoord(i,1));
#endif

}



/*********************************************************************************************************/
void PromoteMesh_Pn_Kirby(const int degree, const EPointType & pointType,
                 const FieldContainer<int>    & P1_elemToNode,
                 const FieldContainer<double> & P1_nodeCoord,
                 const FieldContainer<double> & P1_edgeCoord,
                 const FieldContainer<int>    & P1_elemToEdge,
                 const FieldContainer<int>    & P1_elemToEdgeOrient,
                 const FieldContainer<int>    & P1_nodeOnBoundary,
                 FieldContainer<int>          & Pn_elemToNode,
                 FieldContainer<double>       & Pn_nodeCoord,
                 FieldContainer<int>          & Pn_nodeOnBoundary,
                 std::vector<int> & Pn_edgeNodes,
                 std::vector<int> & Pn_cellNodes) {
//#define DEBUG_PROMOTE_MESH
  int numElems           = P1_elemToNode.dimension(0);
  int P1_numNodesperElem = P1_elemToNode.dimension(1);
  int P1_numEdgesperElem = P1_elemToEdge.dimension(1);
  int Pn_numNodesperElem = Pn_elemToNode.dimension(1);
  int P1_numNodes        = P1_nodeCoord.dimension(0);
  int P1_numEdges        = P1_edgeCoord.dimension(0);
  int dim                = P1_nodeCoord.dimension(1);


#ifdef DEBUG_PROMOTE_MESH
  int Pn_numNodes        = Pn_nodeCoord.dimension(0);
#endif

  int Pn_ExpectedNodesperElem = P1_numNodesperElem + (degree-1)*P1_numEdgesperElem + (degree-1)*(degree-1);
  int Pn_ExpectedNumNodes     = P1_numNodes + (degree-1)*P1_numEdges + (degree-1)*(degree-1)*numElems;  

  // Sanity checks
  if(P1_numNodesperElem !=4 || Pn_numNodesperElem !=Pn_ExpectedNodesperElem ) throw std::runtime_error("Error: PromoteMesh_Pn_Kirby only works on Quads!");
  if(P1_elemToEdge.dimension(0)!=numElems || P1_elemToEdge.dimension(1)!=4 || P1_elemToEdge.dimension(0)!=P1_elemToEdgeOrient.dimension(0) || P1_elemToEdge.dimension(1)!=P1_elemToEdgeOrient.dimension(1) ||
     Pn_elemToNode.dimension(0)!=numElems || Pn_nodeCoord.dimension(0) != Pn_ExpectedNumNodes)
    throw std::runtime_error("Error: PromoteMesh_Pn_Kirby array size mismatch");

  const CellTopologyData &cellTopoData = *shards::getCellTopologyData<shards::Quadrilateral<4> >();                                                  
  shards::CellTopology cellTopo(&cellTopoData);

  // Kirby elements are ordered strictly in lex ordering from the bottom left to the top right
  /*
    Kirby Quad-9 Layout (p=2)
    inode6 -- inode7 -- inode8    
    |                   |
    inode3    inode4    inode5
    |                   |       
    inode0 -- inode1 -- inode2


    Kirby Quad-16 Layout (p=3)
    inode12-- inode13-- inode14-- inode15
    |                             |       
    inode8    inode9    inode10   inode11
    |                             |       
    inode4    inode5    inode6    inode7
    |                             |       
    inode0 -- inode1 -- inode2 -- inode3
  */

  // We still number the global nodes in the exodus-style (node,edge,face) ordering, *by global orientation* but
  // the elem2node map has to account for Kirby-style elements
  int p1_node_in_pn[4] = {0,degree, (degree+1)*(degree+1)-1, degree*(degree+1)};
  int edge_node0_id[4]={0,1,2,3};
  int edge_node1_id[4]={1,2,3,0};

  // As you advance along each edge, we'll skip by the following amount to get to the next Kirby dof
  // edge_skip is in *local* orientation, so we'll multiply it by the elemToEdgeOrient value
  int edge_skip[4] = {1,degree+1,-1,-(degree+1)};
  int center_root = degree+2;

  // Make the new el2node array
  for(int i=0; i<numElems; i++)  {    
    // P1 nodes
    for(int j=0; j<P1_numNodesperElem; j++)
      Pn_elemToNode(i,p1_node_in_pn[j]) = P1_elemToNode(i,j);
  
    // P1 edges
    for(int j=0; j<P1_numEdgesperElem; j++){
      int orient   = P1_elemToEdgeOrient(i,j);
      int base_id = (orient==1) ? p1_node_in_pn[edge_node0_id[j]] : p1_node_in_pn[edge_node1_id[j]];
      int skip     =  orient*edge_skip[j];
      for(int k=0; k<degree-1; k++) {
        int node = P1_numNodes+P1_elemToEdge(i,j)*(degree-1)+k;
        Pn_elemToNode(i,base_id+(k+1)*skip) = node;
        Pn_edgeNodes.push_back(node);
      }
    }

    // P1 cells
    for(int j=0; j<degree-1; j++) 
      for(int k=0; k<degree-1; k++) {
        int node = P1_numNodes+P1_numEdges*(degree-1)+i*(degree-1)*(degree-1) +  j*(degree-1)+k;
        Pn_elemToNode(i,center_root+j*(degree+1)+k) = node;
        Pn_cellNodes.push_back(node);
      }

  }

  // Get the p=n node locations
  FieldContainer<double> RefNodeCoords(Pn_numNodesperElem,dim);  
  FieldContainer<double> RefNodeCoords2(1,Pn_numNodesperElem,dim);  
  FieldContainer<double> PhysNodeCoords(1,Pn_numNodesperElem,dim);  
  Basis_HGRAD_QUAD_Cn_FEM<double, FieldContainer<double> > BasisPn(degree,pointType);
  BasisPn.getDofCoords(RefNodeCoords);

#ifdef DEBUG_PROMOTE_MESH
  printf("**** P=%d reference nodal coordinates ***\n",degree);
  for(int i=0; i<Pn_numNodesperElem; i++)
    printf("[%2d] %10.2f %10.2f\n",i,RefNodeCoords(i,0),RefNodeCoords(i,1));
#endif

  for(int i=0; i<Pn_numNodesperElem; i++)
    for(int k=0; k<dim; k++)
      RefNodeCoords2(0,i,k) = RefNodeCoords(i,k);

  // Make the new coordinates (inefficient, but correct)
 for(int i=0; i<numElems; i++)  {    
   // Get Element coords
   FieldContainer<double> my_p1_nodes(1,P1_numNodesperElem,dim);
   for(int j=0; j<P1_numNodesperElem; j++) {
     int jj = P1_elemToNode(i,j);
     for(int k=0; k<dim; k++)
       my_p1_nodes(0,j,k) = P1_nodeCoord(jj,k);
   }

   // Map the reference node location to physical   
   Intrepid::CellTools<double>::mapToPhysicalFrame(PhysNodeCoords,RefNodeCoords2,my_p1_nodes,cellTopo);


#ifdef DEBUG_PROMOTE_MESH
   printf("[%2d] PhysNodes  : ",i);
   for(int j=0; j<Pn_numNodesperElem; j++)
     printf("(%10.2f %10.2f) ",PhysNodeCoords(0,j,0),PhysNodeCoords(0,j,1));
   printf("\n");
#endif


   // Copy the physical node locations to the Pn_nodeCoord array
   for(int j=0; j<Pn_numNodesperElem; j++) {
     int jj = Pn_elemToNode(i,j);
     for(int k=0; k<dim; k++)
       Pn_nodeCoord(jj,k) = PhysNodeCoords(0,j,k);
   }

 }

#ifdef DEBUG_PROMOTE_MESH
 for(int i=0; i<numElems; i++)  { 
   printf("[%2d] Effective  : ",i);
   for(int j=0; j<Pn_numNodesperElem; j++)
     printf("(%10.2f %10.2f) ",Pn_nodeCoord(Pn_elemToNode(i,j),0),Pn_nodeCoord(Pn_elemToNode(i,j),1));
   printf("\n");
 }
#endif

  // Update the boundary conditions

  // P1 nodes
  for(int i=0; i<P1_numNodes; i++) 
    Pn_nodeOnBoundary(i) = P1_nodeOnBoundary(i);
  
  // P1 edges
  // Not all that efficient
  for(int i=0; i<numElems; i++) {
    for(int j=0; j<P1_numEdgesperElem; j++) {
      int n0 = P1_elemToNode(i,edge_node0_id[j]);
      int n1 = P1_elemToNode(i,edge_node1_id[j]);
      if(P1_nodeOnBoundary(n0) && P1_nodeOnBoundary(n1)) {
        for(int k=0; k<degree-1; k++) 
          Pn_nodeOnBoundary(P1_numNodes+P1_elemToEdge(i,j)*(degree-1)+k) =1;
      }
    }
  }

#ifdef DEBUG_PROMOTE_MESH
  // debug
  printf("**** P=%d elem2node  ***\n",degree);
  for(int i=0; i<numElems; i++) {
    printf("[%2d] ",i);
    for(int j=0; j<Pn_numNodesperElem; j++)
      printf("%2d ",Pn_elemToNode(i,j));
    printf("\n");
  }


  printf("**** P=%d nodes  ***\n",degree);
  for(int i=0; i<Pn_numNodes; i++) {
    printf("[%2d] %10.2e %10.2e (%d)\n",i,Pn_nodeCoord(i,0),Pn_nodeCoord(i,1),(int)Pn_nodeOnBoundary(i));   
  }
  
  printf("**** P=1 nodes  ***\n");
  for(int i=0; i<P1_numNodes; i++) {
    printf("[%2d] %10.2e %10.2e (%d)\n",i,P1_nodeCoord(i,0),P1_nodeCoord(i,1),(int)P1_nodeOnBoundary(i));   
  }
#endif


}


/*********************************************************************************************************/

void CreateP1MeshFromPnMesh(int degree,
                            FieldContainer<int> const    & Pn_elemToNode,
                            FieldContainer<int>          & P1_elemToNode)
{

  /*
    Main idea:
    For each P=n element
    Create n^2 P1 elements
    For each of those P1 elements
    Create element to node map
  */

  int Pn_numElems        = Pn_elemToNode.dimension(0);

  int p1ElemCtr=0;
  for (int i=0; i<Pn_numElems; ++i) {

    // Kirby elements are ordered strictly in lex ordering from the bottom left to the top right
    /*
      Kirby Quad-9 Layout (p=2)
      inode6 -- inode7 -- inode8    
      |                   |
      inode3    inode4    inode5
      |                   |       
      inode0 -- inode1 -- inode2


      Kirby Quad-16 Layout (p=3)
      inode12-- inode13-- inode14-- inode15
      |                             |       
      inode8    inode9    inode10   inode11
      |                             |       
      inode4    inode5    inode6    inode7
      |                             |       
      inode0 -- inode1 -- inode2 -- inode3
    */

    /*
      How "P1" elements are traversed in a Kirby P2 element

      inode6 -- inode7 -- inode8

      |   2    |     3    |

      inode3 -- inode4 -- inode5

      |   0    |     1    |       

      inode0 -- inode1 -- inode2
    */

    int dp1 = degree+1;
    for (int j=0; j<degree; ++j) {  //"y" direction
      for (int k=0; k<degree; ++k) {  //"x" direction
      
        int lowerLeftNode = dp1*j + k;
        //Exodus ordering (counter-clockwise)
        P1_elemToNode(p1ElemCtr,0) = Pn_elemToNode(i,lowerLeftNode);
        P1_elemToNode(p1ElemCtr,1) = Pn_elemToNode(i,lowerLeftNode+1);
        P1_elemToNode(p1ElemCtr,2) = Pn_elemToNode(i,lowerLeftNode+dp1+1);
        P1_elemToNode(p1ElemCtr,3) = Pn_elemToNode(i,lowerLeftNode+dp1);
        p1ElemCtr++;

      }
    }

  }
} //CreateP1MeshFromPnMesh

/*********************************************************************************************************/

void CreateLinearSystem(int numWorksets,
                        int desiredWorksetSize,
                        FieldContainer<int>    const &elemToNode,
                        FieldContainer<double> const &nodeCoord,
                        FieldContainer<double> const &cubPoints,
                        FieldContainer<double> const &cubWeights,
                        RCP<Basis<double,FieldContainer<double> > >&myBasis_rcp,
                        FieldContainer<double> const &HGBGrads,
                        FieldContainer<double> const &HGBValues,
                        std::vector<global_ordinal_type>       const &globalNodeIds,
                        crs_matrix_type &StiffMatrix,
                        RCP<multivector_type> &rhsVector,
                        std::string &msg
                        )
{
  int numCubPoints = cubPoints.dimension(0);
  int cubDim = cubPoints.dimension(1);
  int spaceDim = nodeCoord.dimension(1);
  int numFieldsG = HGBGrads.dimension(0);
  long long numElems = elemToNode.dimension(0);
  int numNodesPerElem = elemToNode.dimension(1);


  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Allocate arrays")));
  
  std::cout << "CreateLinearSystem:" << std::endl;
  std::cout << "     numCubPoints = " << numCubPoints << std::endl;
  std::cout << "     cubDim = " << cubDim << std::endl;
  std::cout << "     spaceDim = " << spaceDim << std::endl;
  std::cout << "     numFieldsG = " << numFieldsG << std::endl;
  std::cout << "     numElems = " << numElems << std::endl;
  std::cout << "     numNodesPerElem = " << numNodesPerElem << std::endl;
  std::cout << "     length(globalNodeIds) = " << globalNodeIds.size() << std::endl;
  std::cout << "     length(nodeCoord) = " << nodeCoord.dimension(0) << std::endl;

  if(nodeCoord.dimension(0) != Teuchos::as<int>(globalNodeIds.size())) {
    std::ostringstream errStr;
    errStr << "Error: CreateLinearSystem: length of coordinates != #nodes ("
           << nodeCoord.dimension(0) << "!=" << globalNodeIds.size() << ")";
    throw std::runtime_error(errStr.str());
  }

  for(int workset = 0; workset < numWorksets; workset++){

    // Compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
    worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // Now we know the actual workset size and can allocate the array for the cell nodes
    worksetSize  = worksetEnd - worksetBegin;
    //std::cout << "     cellWorkset: begin=" << worksetBegin << ", end=" << worksetEnd
    //          << ", size=" << worksetSize << std::endl;
    FieldContainer<double> cellWorkset(worksetSize, numNodesPerElem, spaceDim);

    // Copy coordinates into cell workset
    int cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){
      for (int node = 0; node < numNodesPerElem; node++) {
        cellWorkset(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
        cellWorkset(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);
      }
      cellCounter++;
    }

    /**********************************************************************************/
    /*                                Allocate arrays                                 */
    /**********************************************************************************/

    // Containers for Jacobians, integration measure & cubature points in workset cells
    FieldContainer<double> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet  (worksetSize, numCubPoints);
    FieldContainer<double> worksetCubWeights(worksetSize, numCubPoints);
    FieldContainer<double> worksetCubPoints (worksetSize, numCubPoints, cubDim);

    // Containers for basis values transformed to workset cells and them multiplied by cubature weights
    FieldContainer<double> worksetHGBValues        (worksetSize, numFieldsG, numCubPoints);
    FieldContainer<double> worksetHGBValuesWeighted(worksetSize, numFieldsG, numCubPoints);
    FieldContainer<double> worksetHGBGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
    FieldContainer<double> worksetHGBGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for diffusive & advective fluxes & non-conservative adv. term and reactive terms
    FieldContainer<double> worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for material values and source term. Require user-defined functions
    FieldContainer<double> worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetSourceTerm   (worksetSize, numCubPoints);

    // Containers for workset contributions to the discretization matrix and the right hand side
    FieldContainer<double> worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
    FieldContainer<double> worksetRHS         (worksetSize, numFieldsG);

    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Calculate Jacobians")));

    /**********************************************************************************/
    /*                                Calculate Jacobians                             */
    /**********************************************************************************/
    // Due to the way setJacobian works, we need to use the basis form here
    //    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, myBasis_rcp);
    IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

    //#define DEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
      printf("*** cubPoints ***\n");
      for(int j=0; j<numCubPoints; j++)
        printf("(%d) [%10.2e %10.2e]\n",j,cubPoints(j,0),cubPoints(j,1));

      printf("*** cellWorkset ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numNodesPerElem; j++)
          printf("(%d,%d) [%10.2e %10.2e]\n",i,j,cellWorkset(i,j,0),cellWorkset(i,j,1));
      
      printf("*** Jacobian ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) [%10.2e %10.2e; %10.2e %10.2e]\n",i,j,zwrap(worksetJacobian(i,j,0,0)),zwrap(worksetJacobian(i,j,0,1)),zwrap(worksetJacobian(i,j,1,0)),zwrap(worksetJacobian(i,j,1,1)));

      printf("*** det J ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) %10.2e\n",i,j,zwrap(worksetJacobDet(i,j)));

#endif
    

    /**********************************************************************************/
    /*          Cubature Points to Physical Frame and Compute Data                    */
    /**********************************************************************************/

    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Map to physical frame and get source term")));

    // map cubature points to physical frame
       IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, myBasis_rcp);

    // get A at cubature points
    evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

    // get source term at cubature points
    evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);

    /**********************************************************************************/
    /*                         Compute Stiffness Matrix                               */
    /**********************************************************************************/

    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Compute stiffness matrix")));

    // Transform basis gradients to physical frame:
    IntrepidFSTools::HGRADtransformGRAD<double>(worksetHGBGrads,                // DF^{-T}(grad u)
                                                worksetJacobInv,   HGBGrads);

    // Compute integration measure for workset cells:
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeights,              // Det(DF)*w = J*w
                                                worksetJacobDet, cubWeights);


    // Multiply transformed (workset) gradients with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBGradsWeighted,           // DF^{-T}(grad u)*J*w
                                             worksetCubWeights, worksetHGBGrads);


    // Compute the diffusive flux:
    IntrepidFSTools::tensorMultiplyDataField<double>(worksetDiffusiveFlux,      //  A*(DF^{-T}(grad u)
                                                     worksetMaterialVals,
                                                     worksetHGBGrads);

    // Integrate to compute workset diffusion contribution to global matrix:
    IntrepidFSTools::integrate<double>(worksetStiffMatrix,                    //(DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
                                       worksetHGBGradsWeighted,
                                       worksetDiffusiveFlux, COMP_BLAS);

    /**********************************************************************************/
    /*                                   Compute RHS                                  */
    /**********************************************************************************/

    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Compute right-hand side")));

    // Transform basis values to physical frame:
    IntrepidFSTools::HGRADtransformVALUE<double>(worksetHGBValues,              // clones basis values (u)
                                                 HGBValues);

    // Multiply transformed (workset) values with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBValuesWeighted,          // (u)*J*w
                                             worksetCubWeights, worksetHGBValues);

    // Integrate worksetSourceTerm against weighted basis function set
    IntrepidFSTools::integrate<double>(worksetRHS,                             // f.(u)*J*w
                                       worksetSourceTerm,
                                       worksetHGBValuesWeighted,  COMP_BLAS);

    /**********************************************************************************/
    /***************************** STATISTICS (Part II) ******************************/
    /**********************************************************************************/
    for(int i=0; i<worksetSize; i++) {
      // 0 - Material property
      // 1 - Max/min edge - ratio of max to min edge length
      // 2 - det of cell Jacobian (later)
      double elementdetJ = 0.0, elementWeight=0.0;
      for(int j=0; j<numCubPoints; j++) {
        elementdetJ   += worksetJacobDet(i,j) * worksetCubWeights(i,j);
        elementWeight += worksetCubWeights(i,j);
      }
      double detJ = elementdetJ / elementWeight;
      local_stat_max[2] = std::max(local_stat_max[2],detJ);
      local_stat_min[2] = std::min(local_stat_min[2],detJ);
      local_stat_sum[2] += detJ;
    }


#ifdef DEBUG_OUTPUT      
      printf("*** worksetSourceTerm ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) %10.2e\n",i,j,worksetSourceTerm(i,j));


      printf("*** worksetCubWeights ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) %10.2e\n",i,j,worksetCubWeights(i,j));


      printf("*** worksetHGBValues ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numFieldsG; j++) {
          printf("(%d,%d) ",i,j);
          for(int k=0; k<numCubPoints; k++) {
            printf("%10.2e ",zwrap(worksetHGBValues(i,j,k)));       
          }
          printf("\n");
        }
      
      printf("*** worksetHGBValuesWeighted ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numFieldsG; j++) {
          printf("(%d,%d) ",i,j);
          for(int k=0; k<numCubPoints; k++) {
            printf("%10.2e ",zwrap(worksetHGBValues(i,j,k)));       
          }
          printf("\n");
        }     

      printf("*** worksetRHS ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numFieldsG; j++)
          printf("(%d,%d) %10.2e\n",i,j,worksetRHS(i,j));

      printf("*** worksetStiffMatrix ***\n");
      for(int i=0; i<worksetSize; i++) {
        printf("(%d) [ ",i);
        for(int j=0; j<numFieldsG; j++) {
          for(int k=0; k<numFieldsG; k++) {
            printf("%10.2e ",zwrap(worksetStiffMatrix(i,j,k)));     
          }
          if(j!=numFieldsG-1) printf(";");
        }
        printf("]\n");
      }     



#endif

    /**********************************************************************************/
    /*                         Assemble into Global Matrix and RHS                    */
    /**********************************************************************************/

    //"WORKSET CELL" loop: local cell ordinal is relative to numElems

    Array<scalar_type> vals1(1);
    Array<global_ordinal_type> cols1(1);
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;

      // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
      for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

        int localRow  = elemToNode(cell, cellRow);
        int globalRow = globalNodeIds[localRow];
        double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);


        rhsVector->sumIntoGlobalValue(globalRow, 0, sourceTermContribution);

        // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
        for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

          int localCol  = elemToNode(cell, cellCol);
          int globalCol = globalNodeIds[localCol];
          double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);
          cols1[0] = globalCol;
          vals1[0] = operatorMatrixContribution;
          StiffMatrix.insertGlobalValues(globalRow, cols1(), vals1());
	  //	  printf("Inserting A(%d,%d) = %6.4e\n",(int)globalCol,(int)cols1[0],vals1[0]);
        }// *** cell col loop ***
      }// *** cell row loop ***
    }// *** workset cell loop **

  }// *** workset loop ***

  StiffMatrix.fillComplete();

/**********************************************************************************/
/***************************** STATISTICS (Part IIb) ******************************/
/**********************************************************************************/
  Teuchos::RCP<const crs_graph_type> gl_StiffGraph = StiffMatrix.getCrsGraph();
  Teuchos::RCP<const driver_map_type> rowMap = gl_StiffGraph->getRowMap();
  Teuchos::RCP<multivector_type> coordsOwnedPlusShared;

  RCP<multivector_type> nCoord= rcp(new multivector_type(rowMap,2));
  {
    auto nCoordD = nCoord->get2dViewNonConst();
    for (int inode=0; inode<nodeCoord.dimension(0); inode++) {
      nCoordD[0][inode]=nodeCoord(inode,0);
      nCoordD[1][inode]=nodeCoord(inode,1);
    }  
  }

  if (!(gl_StiffGraph->getImporter().is_null())) {
    coordsOwnedPlusShared = rcp(new multivector_type(gl_StiffGraph->getColMap(), 3, true));
    coordsOwnedPlusShared->doImport(*nCoord, *gl_StiffGraph->getImporter(), Tpetra::CombineMode::ADD, false);
  }
  else {
    coordsOwnedPlusShared = nCoord;
  }

  vector_type laplDiagOwned(rowMap, true);
  Teuchos::ArrayView<const local_ordinal_type> indices;
  Teuchos::ArrayView<const scalar_type> values;
  size_t numOwnedRows = rowMap->getNodeNumElements();
  for (size_t row=0; row<numOwnedRows; row++) {
    StiffMatrix.getLocalRowView(row, indices, values);
    size_t numIndices = indices.size();
    for (size_t j=0; j<numIndices; j++) {
      size_t col = indices[j];
      if (row == col) continue;
      laplDiagOwned.sumIntoLocalValue(row, 1/myDistance2(*coordsOwnedPlusShared, row, col));
    }
  }
  Teuchos::RCP<vector_type> laplDiagOwnedPlusShared;
  if (!gl_StiffGraph->getImporter().is_null()) {
    laplDiagOwnedPlusShared = rcp(new vector_type(gl_StiffGraph->getColMap(), true));
    laplDiagOwnedPlusShared->doImport(laplDiagOwned, *gl_StiffGraph->getImporter(), Tpetra::CombineMode::ADD, false);
  }
  else {
    laplDiagOwnedPlusShared = rcp(&laplDiagOwned, false);
  }
  for (size_t row=0; row<numOwnedRows; row++) {
    StiffMatrix.getLocalRowView(row, indices, values);
    size_t numIndices = indices.size();
    for(size_t j=0; j<numIndices; j++) {
      size_t col = indices[j];
      if (row==col) continue;
      double laplVal = 1.0 / myDistance2(*coordsOwnedPlusShared, row, col);
      double aiiajj = std::abs(laplDiagOwnedPlusShared->getData()[row]*laplDiagOwnedPlusShared->getData()[col]);
      double aij = laplVal * laplVal;
      double ratio = sqrt(aij / aiiajj);
      local_stat_max[7] = std::max(local_stat_max[7], ratio);
      local_stat_min[7] = std::min(local_stat_min[7], ratio);
      local_stat_sum[7] += ratio;
    }
  }

} //CreateLinearSystem

/*********************************************************************************************************/
void GenerateIdentityCoarsening_pn_to_p1(const FieldContainer<int> & Pn_elemToNode,
                      RCP<const driver_map_type> const & P1_map_aux, RCP<const driver_map_type> const &Pn_map,
                      RCP<crs_matrix_type> & P,
                      RCP<crs_matrix_type> & R) {

  // Generate prolongator matrix that will be used to transfer from P1 on auxiliary mesh to Pn on base mesh.
  // It's just the identity matrix.
  // By construction, the node numbering on the P1 auxiliary mesh is the same as on the Pn base mesh
  // (see CreateP1MeshFromPnMesh).  The P2 map is the range map, the P1 auxiliary map is the domain map.
  
  double one = 1.0;
  P = rcp(new crs_matrix_type(Pn_map,1));

  //We must keep track of the nodes already encountered.  Inserting more than once will cause
  //the values to be summed.  Using a hashtable would work -- we abuse std::map for this purpose.
  std::map<int,int> hashTable;
  int Nelem=Pn_elemToNode.dimension(0);
  Array<scalar_type> vals1(1);
  vals1[0] = one;
  Array<global_ordinal_type> cols1(1);
  for(int i=0; i<Nelem; i++) {
    for(int j=0; j<Pn_elemToNode.dimension(1); j++) {
      int row = Pn_elemToNode(i,j);
      if (hashTable.find(row) == hashTable.end()) {
        //not found
        cols1[0] = row;
        P->insertGlobalValues(row,cols1(),vals1());
        hashTable[row] = 1;
      }
    }
  }
  P->fillComplete(P1_map_aux,Pn_map);

  /*
  //JJH FIXME no need to generate R ... it's the identity, stupid
  hashTable.clear();
  R = rcp(new crs_matrix_type(P1_map_aux,1));
  Nelem = P1_elemToNode.dimension(0);
  int nodesPerElem = P1_elemToNode.dimension(1);
  for(int i=0; i<Nelem; ++i) {
    for(int j=0; j<nodesPerElem; ++j) {
      int row = P1_elemToNode(i,j);
      if (hashTable.find(row) == hashTable.end()) {
        //not found
        cols1[0] = row;
        R->insertGlobalValues(row,cols1(),vals1());
        hashTable[row] = 1;
      }
    }
  }
  R->fillComplete(Pn_map,P1_map_aux);
  */

}


/*********************************************************************************************************/
void Apply_Dirichlet_BCs(std::vector<int> &BCNodes, crs_matrix_type & A, multivector_type & x, multivector_type & b,
                         const multivector_type & soln) {
  int N=(int)BCNodes.size();
  ArrayRCP<scalar_type> xdata = x.getDataNonConst(0);
  ArrayRCP<scalar_type> bdata = b.getDataNonConst(0);
  ArrayRCP<const scalar_type> solndata = soln.getData(0);

  A.resumeFill();

  for(int i=0; i<N; i++) {
    local_ordinal_type lrid = BCNodes[i];

    xdata[lrid]=bdata[lrid] = solndata[lrid];

    size_t numEntriesInRow = A.getNumEntriesInLocalRow(lrid);
    Array<local_ordinal_type> cols(numEntriesInRow);
    Array<scalar_type> vals(numEntriesInRow);
    A.getLocalRowCopy(lrid, cols(), vals(), numEntriesInRow);
    
    for(int j=0; j<vals.size(); j++)
      vals[j] = (cols[j] == lrid) ? 1.0 : 0.0;

    A.replaceLocalValues(lrid, cols(), vals());
  }

  A.fillComplete();

}
