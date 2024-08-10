// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   example_StabilizedADR.cpp
    \brief  Example solution of a steady-state advection-diffusion-reaction equation with Dirichlet
            boundary conditon on a hexahedral mesh using nodal (Hgrad) elements and stabilization.

            This example requires the following Trilinos packages:
    \li     Pamgen to generate a Hexahedral mesh;
    \li     Sacado to form the source term from user-specified manufactured solution, and
            diffusion, advection and reaction functions.
    \li     Intrepid to build the discretization matrix and right-hand side
    \li     ML to solve the linear system.

    \verbatim
     Steady-state advection-diffusion boundary value problem in conservative form:

            -div (A.grad(u) - b.u) + c.u = f  in Omega
                                       u = 0  on Gamma
     where
            A   is a symmetric, positive definite diffusivity tensor
            b   is advective velocity vector which is not required to be divergence free
            c   is a reaction coefficient
            f   is a given source term

     The user must provide definitions for these quantities and definition for the exact solution


     Corresponding discrete linear system for nodal coefficients(x):

                 Kx = d

            K - HGrad stiffness matrix
            d - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson, D. Hensinger, C. Siefert.

    \remark Usage:
    \code   ./TrilinosCouplings_examples_scaling_example_StabilizedADR.exe \endcode

    \remark Example requires Pamgen formatted mesh input file named ADR.xml.
*/

//#define DUMP_DATA
// TrilinosCouplings includes
#include "TrilinosCouplings_config.h"
#include "TrilinosCouplings_Pamgen_Utils.hpp"

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

#ifdef HAVE_INTREPID_KOKKOS
#include "Sacado.hpp"
#else
// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#endif

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Import.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "pamgen_im_exodusII_l.h"
#include "pamgen_im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

// STL includes
#include "map"


// Typedefs: FAD with # of ind. vars fixed at 3
typedef Sacado::Fad::SFad<double,3>   Fad3;
typedef shards::CellTopology          ShardsCellTopology;
typedef Intrepid::FunctionSpaceTools  IntrepidFSTools;
typedef Intrepid::CellTools<double>   IntrepidCTools;

typedef long long                     int128;

using namespace std;
using namespace Intrepid;

// Global physics data
double g_advection[3]={0.,0.,0.};
double g_reaction;


/**************************************************************************************************
 *                                                                                                *
 *    Problem specification: Declaration of user-defined functions for                            *
 *    exact solution, diffusion, advection and reaction terms                                     *
 *                                                                                                *
 **************************************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \return Value of the exact solution at (x,y,z)
 */
template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z);



/** \brief  User-defined diffusion tensor.

    \param  diffusion   [out]   3 x 3 diffusivity tensor evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y,z).
*/
template<typename Scalar>
void diffusionTensor(Scalar diffusion[][3], const Scalar&  x, const Scalar&  y, const Scalar&  z);



/** \brief  User-defined advective vector.

    \param  advection   [out]   3-dimensional advective vector evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \remark Advective vector is not required to be divergence-free
*/
template<typename Scalar>
void advectiveVector(Scalar advection[3], const Scalar& x, const Scalar& y, const Scalar& z);



/** \brief  User-defined reaction coefficient.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \return Value of the reaction coefficient function at (x,y,z)
*/
template<typename Scalar>
const Scalar reactionTerm(const Scalar& x, const Scalar& y, const Scalar& z);


/**************************************************************************************************
 *                                                                                                *
 *    Declarations of auxiliary functions that require user-defined functions for                 *
 *    exact solution, diffusion, advection and reaction terms                                     *
 *                                                                                                *
 **************************************************************************************************/


/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y,z)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point
 */
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z);



/** \brief Computes source term: f = -div(K.grad u - b.u) + c.u. Requires user-defined exact solution,
           diffusion tensor, advective vector and reaction coefficient.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution, diffusivity, advection and
            reaction coefficients, evaluated at (x,y,z)
 */
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z);


/**************************************************************************************************
 *                                                                                                *
 *    Evaluation methods accepting multidimensional arrays of points                              *
 *                                                                                                *
 **************************************************************************************************/


/** \brief Computation of the diffusion tensor at array of points in physical space.

    \param worksetDiffusionValues      [out]     Rank-2, 3 or 4 array with dimensions (C,P), (C,P,D) or (C,P,D,D)
                                                with the values of the diffusion tensor
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateDiffusionTensor(ArrayOut &        worksetDiffusionValues,
                             const ArrayIn &   evaluationPoints);



/** \brief Computation of the advective vector at array of points in physical space..

    \param worksetAdvectionValues      [out]     Rank-3 (C,P,D) array with the values of the advective vector
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateAdvectiveVector(ArrayOut &       worksetAdvectionValues,
                             const ArrayIn &  evaluationPoints);



/** \brief Computation of the reaction coefficient at array of points in physical space..

    \param worksetReactionValues  [out]     Rank-2 (C,P) array with the values of the reaction coefficient
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateReactionCoefficient(ArrayOut &       worksetReactionValues,
                                 const ArrayIn &  evaluationPoints);



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


/**************************************************************************************************
 *                                                                                                *
 *    ML related declarations                                                                     *
 *                                                                                                *
 **************************************************************************************************/
int TestMultiLevelPreconditioner(char                        ProblemType[],
                                 Teuchos::ParameterList &    MLList,
                                 Epetra_CrsMatrix &          A,
                                 const Epetra_MultiVector &  xexact,
                                 Epetra_MultiVector &        b,
                                 Epetra_MultiVector &        uh,
                                 double &                    TotalErrorResidual,
                                 double &                    TotalErrorExactSol);



/**************************************************************************************************
 *                                                                                                *
 *    PAMGEN wrapper (UNSTABLE - IN CONSTRUCTION)                                                 *
 *                                                                                                *
 **************************************************************************************************/
/** \brief Pamgen wrapper

    \param localNodeCoordsFC          [out]     Rank-2 (LN,D) array with the local nodes
    \param localCellToNodeFC          [out]     Rank-2 (C,N) array with the local nodes connectivity
    \param nodeOnBoundaryFC           [out]     Rank-1 (LN) array tells which local nodes are on the boundary
    \param nodeIsOwnedFC              [out]     Rank-1 (LN) array tells which local nodes are owned
    \param globalNodeIdsFC            [out]     Rank-1 (LN) array with the global node Ids of the local nodes
    \param meshInput                  [in]      The list of parameters for Create_Pamgen_Mesh, extracted from the mesh file
    \param procRank                   [in]      Rank of the processor
    \param numProcs                   [in]      Number of processors
    \param Comm                       [in]      Communicator object
    \param Time                       [in]      Time object used for timings of various mesh tasks in the wrapper
    \param message                    [in]      String that is printed as part of the various diagnostic outputs
    \param verbose                    [in]      Forces diagnostic output if set to 1, default is 0
*/
template<class Scalar>
void getPamgenMesh(FieldContainer<Scalar>    & localNodeCoordsFC,
                   FieldContainer<long long> & localCellToNodeFC,
                   FieldContainer<int> &       nodeOnBoundaryFC,
                   FieldContainer<bool> &      nodeIsOwnedFC,
                   FieldContainer<long long> & globalNodeIdsFC,
                   const std::string &         meshInput,
                   const int &                 procRank,
                   const int &                 numProcs,
                   const Epetra_Comm &         Comm,
                   Epetra_Time &               Time,
                   const string &              message,
                   const int                   verbose = 0);



/**************************************************************************************************
 *                                                                                                *
 *    Argument input (UNSTABLE - IN CONSTRUCTION)                                                 *
 *                                                                                                *
 **************************************************************************************************/
void getInputArguments(Teuchos::ParameterList &  inputMeshList,
                       std::string            &  meshInput,
                       Teuchos::ParameterList &  inputSolverList,
                       int &                     argc,
                       char *                    argv[],
                       const Epetra_Comm &       Comm,
                       const std::string &       inputMeshFile = "ADR.xml",
                       const std::string &       intputSolverFile = "ML_nonsym.xml");



/**************************************************************************************************
 *                                                                                                *
 *    DRIVER: Advection-Diffusion-Reaction Boundary Value Problem                                 *
 *                                                                                                *
 **************************************************************************************************/

int main(int argc, char *argv[]) {
  int    numProcs = 1;
  int    procRank = 0;
  string msg("Advection-diffusion-reaction driver: ");

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  procRank = mpiSession.getRank();
  numProcs = mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  Epetra_Time Time(Comm);

 if (MyPID == 0){
  std::cout \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|  Example: Solve Advection-diffusion-reaction Equation on Hexahedral Mesh    |\n" \
    << "|                                                                             |\n" \
    << "|  Requires user-defined functions for:                                       |\n" \
    << "|                                                                             |\n" \
    << "|    the exact solution         -->   exactSolution(...)                      |\n" \
    << "|    the diffusion tensor       -->   diffusionTensor(...)                    |\n" \
    << "|    the advective vector       -->   advectiveVector(...)                    |\n" \
    << "|    the reaction coefficient   -->   reactionTerm(...)                       |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Pamgen's website:   http://trilinos.sandia.gov/packages/pamgen             |\n" \
    << "|  Sacado's website:   http://trilinos.sandia.gov/packages/sacado             |\n" \
    << "|  ML's website:       http://trilinos.sandia.gov/packages/ml                 |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  }

#ifdef HAVE_MPI
  if (MyPID == 0) {
    std::cout << msg << "PARALLEL executable \n";
  }
#else
  if (MyPID == 0) {
    std::cout << msg << "SERIAL executable \n";
  }
#endif

  /*************************************************************************************************
    *                                                                                              *
    *  DISCRETIZATION SETUP: get inputs, generate mesh, choose cell topology, cubature and basis   *
    *                                                                                              *
    ************************************************************************************************/

  // *********************** Read XML file with mesh specifications *******************************

  Teuchos::ParameterList   inputMeshList;
  std::string              meshInput;
  Teuchos::ParameterList   inputSolverList;
  int verbose = 1;

  getInputArguments(inputMeshList, meshInput, inputSolverList, argc, argv, Comm);


  // ************************************** GENERATE MESH *****************************************

  // Arrays for mesh data: will be resized to required dimensions by the Pamgen wrapper
  FieldContainer<int128>  globalNodeOrdinalsFC;
  FieldContainer<int128>  localCellToNodeFC;
  FieldContainer<double>  localNodeCoordsFC;
  FieldContainer<int>     localNodeOnBoundaryFC;
  FieldContainer<bool>    localNodeIsMineFC;

  // Generate Pamgen mesh and return mesh data in Field Containers
  getPamgenMesh<double>(localNodeCoordsFC,    localCellToNodeFC,      localNodeOnBoundaryFC,
                        localNodeIsMineFC,    globalNodeOrdinalsFC,   meshInput,
                        procRank,             numProcs,               Comm,
                        Time,                 msg,                    verbose);

  // ********************************* GET PHYSICS INFORMATION ************************************
  // Read physics info from mesh xml file if it's there.  Otherwise do nothing.
  if(inputMeshList.isSublist("Physics Input")){
    Teuchos::ParameterList & physicsList=inputMeshList.sublist("Physics Input");
    g_advection[0]=physicsList.get("x advection",0.0);
    g_advection[1]=physicsList.get("y advection",0.0);
    g_advection[2]=physicsList.get("z advection",0.0);
  }

  // *********************************** DEFINE CELL TOPOLOGY *************************************

  // All cells in the Pamgen grid are of the Hexahedron<8> type.
  ShardsCellTopology cellType(shards::getCellTopologyData<shards::Hexahedron<8> >() );

  int numNodesPerCell = cellType.getNodeCount();
  int spaceDim = cellType.getDimension();


  // ************************************* DEFINE CUBATURE ****************************************

  if (MyPID == 0) std::cout << "Getting cubature ... \n\n";

  // Define cubature of the specified degree for the cellType
  DefaultCubatureFactory<double>  cubFactory;
  int cubDegree = 2;
  Teuchos::RCP<Cubature<double> > cellCubature = cubFactory.create(cellType, cubDegree);

  int cubDim       = cellCubature -> getDimension();
  int numCubPoints = cellCubature -> getNumPoints();

  // Get numerical integration points and weights
  FieldContainer<double> cubPoints (numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);

  cellCubature -> getCubature(cubPoints, cubWeights);


  // *************************************** DEFINE BASIS ****************************************

  if (MyPID == 0) std::cout << "Getting basis ... \n\n";

  // Define basis
  Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > HGB;
  int numHGBFields = HGB.getCardinality();

  FieldContainer<double> HGBValues(numHGBFields, numCubPoints);
  FieldContainer<double> HGBGrads (numHGBFields, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points on the reference Hex cell
  HGB.getValues(HGBValues, cubPoints, OPERATOR_VALUE);
  HGB.getValues(HGBGrads,  cubPoints, OPERATOR_GRAD);


  /*************************************************************************************************
    *                                                                                              *
    *  GLOBAL MATRIX PROBLEM SETUP: generate Epetra map & define global Epetra arrays              *
    *                                                                                              *
    ************************************************************************************************/

  // Count number of myDofs: for the lowest-order nodal FEM = number of myNodes
  int numMyDofs     = 0;
  int numMyNodes    = 0;
  int numLocalNodes = localNodeCoordsFC.dimension(0);

  for(int localNode = 0; localNode < numLocalNodes; localNode++) {
    if(localNodeIsMineFC[localNode]) numMyNodes++;
  }
  numMyDofs = numMyNodes;

  // Build a list of the global Dof ordinals of myDofs: for the lowest-order nodal FEM  = global node ordinals
  int *myGlobalDofs       = new int[numMyDofs];
  int  myGlobalDofOrdinal = 0;
  for(int localNode = 0; localNode < numLocalNodes; localNode++){
    if(localNodeIsMineFC[localNode]){
      myGlobalDofs[myGlobalDofOrdinal] = (int)globalNodeOrdinalsFC[localNode];    // NTS: will need to switch back to int128
      myGlobalDofOrdinal++;
    }
  }

  // Define global Epetra arrays
  Epetra_Map          globalMapG            (-1, numMyDofs, myGlobalDofs, 0, Comm);
  Epetra_FECrsMatrix  globalOperatorMatrix  (Copy, globalMapG, 5*numHGBFields);
  Epetra_FEVector     globalSourceTermVector(globalMapG);


  /*************************************************************************************************
    *                                                                                              *
    *  DIRICHLET BOUNDARY CONDITIONS SETUP: count myDofs that are BV and set BVs into a map        *
    *                                                                                              *
    ************************************************************************************************/

  // Map, keyed by local node ordinal, with the bdry values at all local nodes that are on the bdry
  std::map<int, double> localBoundaryValuesMap;

  // Number of myDofs that are also boundary values: for lowest-order FEM = number of myNodes that are on the bdry
  int numMyBoundaryValueDofs = 0;


  // Set the map with the local boundary values & count myDofs that are boundary values
  for(int localNode = 0; localNode < numLocalNodes; localNode++){

    // If the localNode is on the boundary compute and set the exact solution value into the map
    if (localNodeOnBoundaryFC(localNode) ){
      double x      = localNodeCoordsFC(localNode, 0);
      double y      = localNodeCoordsFC(localNode, 1);
      double z      = localNodeCoordsFC(localNode, 2);
      localBoundaryValuesMap[localNode] = exactSolution(x, y, z);

      // If this local node is also "mine", i.e., owned, increment counter for myDofs that are bdry values
      if(localNodeIsMineFC[localNode]){
        numMyBoundaryValueDofs++;
      }
    }
  }

  // Lookup table for ML_Epetra::Apply_OAZToMatrix method:  myBoundaryValueDofs[myBvOrdinal] = myDofOrdinal
  int * myBoundaryValueDofs = new int [numMyBoundaryValueDofs];
  int   myBvOrdinal  = 0;
  int   myDofOrdinal = 0;

  // Count the bdry values and set their associated myDofOrdinals into the lookup table
  for (int localNode = 0; localNode < numLocalNodes; localNode++){
    if (localNodeIsMineFC[localNode]){
      if (localNodeOnBoundaryFC(localNode)){
        myBoundaryValueDofs[myBvOrdinal] = myDofOrdinal;
        myBvOrdinal++;
      }
      myDofOrdinal++;
    }
  }


  /*************************************************************************************************
    *                                                                                              *
    *  ASSEMBLY OF WORKSET CONTRIBUTIONS TO MATRIX AND RIGHT HAND SIDE: workset loop               *
    *                                                                                              *
    ************************************************************************************************/

  // Define desired workset size and count how many worksets there are on this processor's mesh block
  int numLocalGridCells  = localCellToNodeFC.dimension(0);
  int desiredWorksetSize = numLocalGridCells;                      // change to desired workset size!
  //int desiredWorksetSize = 100;                      // change to desired workset size!
  int numWorksets        = numLocalGridCells/desiredWorksetSize;

  // When numLocalGridCells is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksets*desiredWorksetSize < numLocalGridCells) numWorksets += 1;

  if (MyPID == 0) {
    std::cout << "Building discretization matrix and right hand side... \n\n";
    std::cout << "\tDesired workset size: " << desiredWorksetSize <<"\n";
    std::cout << "\tNumber of worksets:   " << numWorksets <<"\n\n";
    Time.ResetStartTime();
  }


  // ************************************** Loop over worksets *************************************

  for(int workset = 0; workset < numWorksets; workset++){

    // Compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // When numLocalGridCells is not divisible by desiredWorksetSize, the last workset ends at numLocalGridCells
    worksetEnd   = (worksetEnd <= numLocalGridCells) ? worksetEnd : numLocalGridCells;

    // Now we know the actual workset size and can allocate the array for the cell nodes
    worksetSize  = worksetEnd - worksetBegin;
    FieldContainer<double> cellWorkset(worksetSize, numNodesPerCell, spaceDim);

    // Copy cell nodes to the cell workset
    int cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){
      for (int node = 0; node < numNodesPerCell; node++) {
        cellWorkset(cellCounter, node, 0) = localNodeCoordsFC( localCellToNodeFC(cell, node), 0);
        cellWorkset(cellCounter, node, 1) = localNodeCoordsFC( localCellToNodeFC(cell, node), 1);
        cellWorkset(cellCounter, node, 2) = localNodeCoordsFC( localCellToNodeFC(cell, node), 2);
      }
      cellCounter++;
    }


    /*************************************************************************************************
      *     Allocate arrays                                                                          *
      ************************************************************************************************/

    // Containers for Jacobians, integration measure & cubature points in workset cells
    FieldContainer<double> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet  (worksetSize, numCubPoints);
    FieldContainer<double> worksetCubWeights(worksetSize, numCubPoints);
    FieldContainer<double> worksetCubPoints (worksetSize, numCubPoints, cubDim);

    // Containers for basis values transformed to workset cells and them multiplied by cubature weights
    FieldContainer<double> worksetHGBValues        (worksetSize, numHGBFields, numCubPoints);
    FieldContainer<double> worksetHGBValuesWeighted(worksetSize, numHGBFields, numCubPoints);
    FieldContainer<double> worksetHGBGrads         (worksetSize, numHGBFields, numCubPoints, spaceDim);
    FieldContainer<double> worksetHGBGradsWeighted (worksetSize, numHGBFields, numCubPoints, spaceDim);

    // Containers for diffusive & advective fluxes & non-conservative adv. term and reactive terms
    FieldContainer<double> worksetDiffusiveFlux(worksetSize, numHGBFields, numCubPoints, spaceDim);
    FieldContainer<double> worksetAdvectiveFlux(worksetSize, numHGBFields, numCubPoints, spaceDim);
    FieldContainer<double> worksetAdvectiveTerm(worksetSize, numHGBFields, numCubPoints);
    FieldContainer<double> worksetReactiveTerm (worksetSize, numHGBFields, numCubPoints);

    // Containers for diffusion, advection and reaction values. Require user-defined functions
    FieldContainer<double> worksetDiffusion (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetAdvection (worksetSize, numCubPoints, spaceDim);
    FieldContainer<double> worksetReaction  (worksetSize, numCubPoints);
    FieldContainer<double> worksetSourceTerm(worksetSize, numCubPoints);

    // Containers for workset contributions to the discretization matrix and the right hand side
    FieldContainer<double> worksetDiffusionMatrix (worksetSize, numHGBFields, numHGBFields);
    FieldContainer<double> worksetAdvectionMatrix (worksetSize, numHGBFields, numHGBFields);
    FieldContainer<double> worksetReactionMatrix  (worksetSize, numHGBFields, numHGBFields);
    FieldContainer<double> worksetSourceTermVector(worksetSize, numHGBFields);


    if(MyPID == 0) {
      std::cout << msg << "Allocate arrays  = " << Time.ElapsedTime() << " sec."<< endl;
      Time.ResetStartTime();
    }




    /*************************************************************************************************
      *     Reference frame computation: Jacobian data for each cell in the workset                  *
      ************************************************************************************************/

    // Compute cell Jacobians, their inverses and their determinants
    IntrepidCTools::setJacobian   (worksetJacobian, cubPoints, cellWorkset, cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );


    /*************************************************************************************************
      *     Physical frame computation: user-defined diffusion, advection, reaction and source terms *
      ************************************************************************************************/

    // Map cubature points to workset cells & compute diffusion, adection, reaction and source terms
    IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);
    evaluateDiffusionTensor            (worksetDiffusion, worksetCubPoints);
    evaluateAdvectiveVector            (worksetAdvection, worksetCubPoints);
    evaluateReactionCoefficient        (worksetReaction,  worksetCubPoints);

    if(MyPID == 0) {
      std::cout << msg << "tabulate problem data " << Time.ElapsedTime() << " sec."<< endl;
      Time.ResetStartTime();
    }


    evaluateSourceTerm                 (worksetSourceTerm, worksetCubPoints);


    if(MyPID == 0) {
      std::cout << msg << "Compute RHS " << Time.ElapsedTime() << " sec."<< endl;
      Time.ResetStartTime();
    }



    /*************************************************************************************************
      *     Transform basis functions and their gradients and compute integration measure            *
      ************************************************************************************************/

    // Transform basis values to physical frame:
    IntrepidFSTools::HGRADtransformVALUE<double>(worksetHGBValues,              // clones basis values (u)
                                                 HGBValues);

    // Transform basis gradients to physical frame:
    IntrepidFSTools::HGRADtransformGRAD<double>(worksetHGBGrads,                // DF^{-T}(grad u)
                                                worksetJacobInv,   HGBGrads);

    // Compute integration measure for workset cells:
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeights,              // Det(DF)*w = J*w
                                                worksetJacobDet, cubWeights);

    // Multiply transformed (workset) values with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBValuesWeighted,          // (u)*J*w
                                             worksetCubWeights, worksetHGBValues);

    // Multiply transformed (workset) gradients with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBGradsWeighted,           // DF^{-T}(grad u)*J*w
                                             worksetCubWeights, worksetHGBGrads);


    if(MyPID == 0) {
      std::cout << msg << "Transform time  = " << Time.ElapsedTime() << " sec."<< endl;
      Time.ResetStartTime();
    }


    /*************************************************************************************************
      *   Diffusion term:   DF^{-T}(grad u)*K*DF^{-T}(grad u)*J*w                                    *
      ************************************************************************************************/

    // Compute the diffusive flux:
    IntrepidFSTools::tensorMultiplyDataField<double>(worksetDiffusiveFlux,      //  K*(DF^{-T}(grad u)
                                                     worksetDiffusion,
                                                     worksetHGBGrads);

    // Integrate to compute workset diffusion contribution to global matrix:
    IntrepidFSTools::integrate<double>(worksetDiffusionMatrix,                  // (DF^{-T}(grad u)*J*w)*(K*DF^{-T}(grad u))
                                       worksetHGBGradsWeighted,
                                       worksetDiffusiveFlux, COMP_BLAS);



    if(MyPID == 0) {
      std::cout << msg << "diffusion time  = " << Time.ElapsedTime() << " sec."<< endl;
      Time.ResetStartTime();
    }



  /*************************************************************************************************
    *     Advection term:   b*DF^{-T}(grad u)*(u)*J*w  requires div b = 0                     *
    ************************************************************************************************/

    // Compute the advective term in non-conservative form: requires div b = 0!
    IntrepidFSTools::dotMultiplyDataField<double>(worksetAdvectiveTerm,         // b*(DF^{-T}*(grad u)
                                                  worksetAdvection,
                                                  worksetHGBGrads);

    // Integrate to compute workset advection contribution to global matrix:
    IntrepidFSTools::integrate<double>(worksetAdvectionMatrix,                  // b*(DF^{-T}*(grad u)(u)*J*w
                                       worksetAdvectiveTerm,
                                       worksetHGBValuesWeighted, COMP_BLAS);



    if(MyPID == 0) {
      std::cout << msg << "advection time  = " << Time.ElapsedTime() << " sec."<< endl;
      Time.ResetStartTime();
    }



   /*************************************************************************************************
     *     Reaction term:  C.(u).(u)*J*w                                                            *
     ************************************************************************************************/

    // Compute the reaction term:
    IntrepidFSTools::scalarMultiplyDataField<double> (worksetReactiveTerm,      // C.(u)
                                                      worksetReaction,
                                                      worksetHGBValues);

    // Integrate to compute workset reaction contribution to global matrix:
    IntrepidFSTools::integrate<double>(worksetReactionMatrix,                    // C.(u).(u)*J*w
                                       worksetReactiveTerm,
                                       worksetHGBValuesWeighted, COMP_BLAS);

///////////////////**************************** TEST ***********************/////////////////////////
/*
    // Compute the advective term in non-conservative form: requires div b = 0!
    IntrepidFSTools::dotMultiplyDataField<double>(worksetAdvectiveTerm,         // b*(DF^{-T}*(grad u)
                                                  worksetAdvection,
                                                  worksetHGBGrads);

    // Compute the reaction term:
    IntrepidFSTools::scalarMultiplyDataField<double> (worksetReactiveTerm,      // C.(u)
                                                      worksetReaction,
                                                      worksetHGBValues);

    // Add both terms
    Intrepid::RealSpaceTools<double>::add(worksetAdvectiveTerm, worksetReactiveTerm);

    // Integrate to compute workset advection contribution to global matrix:
    IntrepidFSTools::integrate<double>(worksetAdvectionMatrix,                  // b*(DF^{-T}*(grad u)(u)*J*w
                                       worksetAdvectiveTerm,
                                       worksetHGBValuesWeighted, COMP_BLAS);
*/

    if(MyPID == 0) {
      std::cout << msg << "reaction time  = " << Time.ElapsedTime() << " sec."<< endl;
      Time.ResetStartTime();
    }


  /*************************************************************************************************
    *      Source term (right hand side): f.(u)*J*w                                                *
    ************************************************************************************************/

    // Integrate worksetSourceTerm against weighted basis function set
    IntrepidFSTools::integrate<double>(worksetSourceTermVector,                  // f.(u)*J*w
                                       worksetSourceTerm,
                                       worksetHGBValuesWeighted,  COMP_BLAS);


  /*************************************************************************************************
    *     Assemble into global matrix and RHS vector & adjust RHS vector for inhomogeneous BC      *
    ************************************************************************************************/

    //"WORKSET CELL" loop: local cell ordinal is relative to numLocalGridCells
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;

      // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
      for (int cellRow = 0; cellRow < numHGBFields; cellRow++){

        int localRow  = localCellToNodeFC(cell, cellRow);
        int globalRow = globalNodeOrdinalsFC[localRow];
        double sourceTermContribution =  worksetSourceTermVector(worksetCellOrdinal, cellRow);

        globalSourceTermVector.SumIntoGlobalValues(1, &globalRow, &sourceTermContribution);

        // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
        for (int cellCol = 0; cellCol < numHGBFields; cellCol++){

          int localCol  = localCellToNodeFC(cell, cellCol);
          int globalCol = globalNodeOrdinalsFC[localCol];
          double operatorMatrixContribution =
            worksetDiffusionMatrix(worksetCellOrdinal, cellRow, cellCol) +
            worksetAdvectionMatrix(worksetCellOrdinal, cellRow, cellCol) +
            worksetReactionMatrix (worksetCellOrdinal, cellRow, cellCol);

          globalOperatorMatrix.InsertGlobalValues(1, &globalRow, 1, &globalCol, &operatorMatrixContribution);

          // Adjust global source term for inhomogeneous boundary conditions:
          if(localNodeOnBoundaryFC(localCol) ){
            map<int,double>::iterator localBvIterator  = localBoundaryValuesMap.find(localCol);

            if(localBvIterator != localBoundaryValuesMap.end() ) {
              double boundaryValue          = localBvIterator -> second;
              double globalSourceAdjustment = -1.0*boundaryValue*operatorMatrixContribution;

              globalSourceTermVector.SumIntoGlobalValues(1, &globalRow, &globalSourceAdjustment);
            }
            else{
              std::cout << " Error: requested boundary value does not exist...exiting!\n"; exit(1);
            }
          }// *** adjust source ***
        }// *** cell col loop ***
      }// *** cell row loop ***
    }// *** workset cell loop **
  }// *** workset loop ***


  /*************************************************************************************************
    *                                                                                              *
    *  ASSEMBLE OVER MULTIPLE PROCESSORS & CONVERT TO LOCAL INDICES                                *
    *                                                                                              *
    ************************************************************************************************/

  globalOperatorMatrix.GlobalAssemble();
  globalOperatorMatrix.FillComplete();
  globalSourceTermVector.GlobalAssemble();

  if(MyPID == 0) {
    std::cout << msg << "Matrix & Right Hand Side Assembly  = " << Time.ElapsedTime() << " sec."<< endl;
    Time.ResetStartTime();
  }


  /*************************************************************************************************
    *                                                                                              *
    *  FINAL ADJUSTMENTS TO GLOBAL MATRIX AND RHS VECTOR FOR INHOMOGENEOUS BOUNDARY CONDITIONS     *
    *                                                                                              *
    ************************************************************************************************/

  /*************************************************************************************************
    *  Step 1 - Adjust RHS vector: set appropriate global RHS entries to the exact boundary values        *
    ************************************************************************************************/

  // Reset myDof counter
  myDofOrdinal = 0;
  for (int localNode = 0; localNode < numLocalNodes; localNode++){
    if (localNodeIsMineFC[localNode]){
      if (localNodeOnBoundaryFC(localNode)){
        map<int,double>::iterator localBvIterator  = localBoundaryValuesMap.find(localNode);

        // Replace the global source term value by the exact boundary value
        if(localBvIterator != localBoundaryValuesMap.end() ) {
          globalSourceTermVector[0][myDofOrdinal] = localBvIterator -> second;
        }
        else{
          std::cout << " Error: requested boundary value does not exist...exiting!\n"; exit(1);
        }
      }
      myDofOrdinal++;
    }
  }
   // Alternative access method to globalSourceTermVector[0][myDofOrdinal]
   // int globalRow = globalNodeIdsFC[localNode];
   // globalSourceTermVector.ReplaceGlobalValue(globalRow, 0, newValue);


  /*************************************************************************************************
    * Step 2 - Adjust MATRIX: zero out Dirichlet rows and columns and add 1 to the diagonal       *
    ************************************************************************************************/

  ML_Epetra::Apply_OAZToMatrix(myBoundaryValueDofs, numMyBoundaryValueDofs, globalOperatorMatrix);
  delete [] myBoundaryValueDofs;


#ifdef DUMP_DATA
  // Dump matrices to disk
  EpetraExt::RowMatrixToMatlabFile("ADR_matrix.dat",globalOperatorMatrix);
  EpetraExt::MultiVectorToMatrixMarketFile("ADR_rhs_vector.dat",globalSourceTermVector,0,0,false);
#endif


  /*************************************************************************************************
    *                                                                                              *
    *  SOLVE THE LINEAR SYSTEM                                                                     *
    *                                                                                              *
    ************************************************************************************************/

  // Run the solver
  Teuchos::ParameterList MLList = inputSolverList;
  ML_Epetra::SetDefaults("SA", MLList, 0, 0, false);
  Epetra_FEVector exactNodalVals(globalMapG);
  Epetra_FEVector femCoefficients(globalMapG);
  double TotalErrorResidual = 0.0;
  double TotalErrorExactSol = 0.0;

  // Get exact solution at nodes
  for (int i = 0; i<numLocalNodes; i++) {
    if (localNodeIsMineFC[i]){
      double x = localNodeCoordsFC(i,0);
      double y = localNodeCoordsFC(i,1);
      double z = localNodeCoordsFC(i,2);
      double exactu = exactSolution(x, y, z);

      int rowindex=globalNodeOrdinalsFC[i];
      exactNodalVals.SumIntoGlobalValues(1, &rowindex, &exactu);
    }
  }
  exactNodalVals.GlobalAssemble();

  char probType[10] = "conv-diff";

  TestMultiLevelPreconditioner(probType,             MLList,
                               globalOperatorMatrix,     exactNodalVals,
                               globalSourceTermVector,       femCoefficients,
                               TotalErrorResidual,   TotalErrorExactSol);


  /*************************************************************************************************
    *                                                                                              *
    *  COMPUTE FINITE ELEMENT ERROR                                                                *
    *                                                                                              *
    ************************************************************************************************/

  double L2err = 0.0;
  double L2errTot = 0.0;
  double H1err = 0.0;
  double H1errTot = 0.0;
  double Linferr = 0.0;
  double LinferrTot = 0.0;

#ifdef HAVE_MPI
  // Import solution onto current processor
  int numNodesGlobal = globalMapG.NumGlobalElements();
  Epetra_Map     solnMap(numNodesGlobal, numNodesGlobal, 0, Comm);
  Epetra_Import  solnImporter(solnMap, globalMapG);
  Epetra_Vector  uCoeff(solnMap);
  uCoeff.Import(femCoefficients, solnImporter, Insert);
#endif

  // Set workset size to 1
  desiredWorksetSize = 1;

  // Get cubature points and weights for error calc (may be different from previous)
  DefaultCubatureFactory<double>  cubFactoryErr;
  int cubDegErr = 3;
  Teuchos::RCP<Cubature<double> > cellCubatureErr = cubFactoryErr.create(cellType, cubDegErr);
  int cubDimErr       = cellCubatureErr->getDimension();
  int numCubPointsErr = cellCubatureErr->getNumPoints();
  FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
  FieldContainer<double> cubWeightsErr(numCubPointsErr);
  cellCubatureErr->getCubature(cubPointsErr, cubWeightsErr);

  // Containers for Jacobian
  FieldContainer<double> worksetJacobianE(desiredWorksetSize, numCubPointsErr, spaceDim, spaceDim);
  FieldContainer<double> worksetJacobInvE(desiredWorksetSize, numCubPointsErr, spaceDim, spaceDim);
  FieldContainer<double> worksetJacobDetE(desiredWorksetSize, numCubPointsErr);
  FieldContainer<double> worksetCubWeightsE(desiredWorksetSize, numCubPointsErr);

  // Evaluate basis values and gradients at cubature points
  FieldContainer<double> uhGVals(numHGBFields, numCubPointsErr);
  FieldContainer<double> uhGValsTrans(desiredWorksetSize,numHGBFields, numCubPointsErr);
  FieldContainer<double> uhGrads(numHGBFields, numCubPointsErr, spaceDim);
  FieldContainer<double> uhGradsTrans(desiredWorksetSize, numHGBFields, numCubPointsErr, spaceDim);
  HGB.getValues(uhGVals, cubPointsErr, OPERATOR_VALUE);
  HGB.getValues(uhGrads, cubPointsErr, OPERATOR_GRAD);


  // Loop over elements
  for (int k=0; k<numLocalGridCells; k++){

    double L2errElem = 0.0;
    double H1errElem = 0.0;
    double uExact;
    double graduExact1, graduExact2, graduExact3;

    // physical cell coordinates
    FieldContainer<double> cellWorksetEr(desiredWorksetSize,numNodesPerCell, spaceDim);

    for (int i=0; i<numNodesPerCell; i++) {
      cellWorksetEr(0,i,0) = localNodeCoordsFC(localCellToNodeFC(k,i),0);
      cellWorksetEr(0,i,1) = localNodeCoordsFC(localCellToNodeFC(k,i),1);
      cellWorksetEr(0,i,2) = localNodeCoordsFC(localCellToNodeFC(k,i),2);
    }

    // compute cell Jacobians, their inverses and their determinants
    IntrepidCTools::setJacobian(worksetJacobianE, cubPointsErr, cellWorksetEr, cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInvE, worksetJacobianE );
    IntrepidCTools::setJacobianDet(worksetJacobDetE, worksetJacobianE );

    // transform integration points to physical points
    FieldContainer<double> worksetCubPoints(desiredWorksetSize,numCubPointsErr, cubDimErr);
    IntrepidCTools::mapToPhysicalFrame(worksetCubPoints, cubPointsErr, cellWorksetEr, cellType);

    // transform basis values to physical coordinates
    IntrepidFSTools::HGRADtransformVALUE<double>(uhGValsTrans, uhGVals);
    IntrepidFSTools::HGRADtransformGRAD<double>(uhGradsTrans, worksetJacobInvE, uhGrads);

    // compute weighted measure
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeightsE, worksetJacobDetE, cubWeightsErr);

    // loop over cubature points
    for (int nPt = 0; nPt < numCubPointsErr; nPt++){

      // get exact solution and gradients using AD
      Fad3 x = worksetCubPoints(0,nPt,0);
      Fad3 y = worksetCubPoints(0,nPt,1);
      Fad3 z = worksetCubPoints(0,nPt,2);
      Fad3 grad_u[3];

      uExact = exactSolution(x, y, z).val();
      exactSolutionGrad(grad_u, x, y, z);

      graduExact1 = grad_u[0].val();
      graduExact2 = grad_u[1].val();
      graduExact3 = grad_u[2].val();

      // calculate approximate solution and gradients
      double uApprox = 0.0;
      double graduApprox1 = 0.0;
      double graduApprox2= 0.0;
      double graduApprox3 = 0.0;
      for (int i = 0; i < numHGBFields; i++){
        int rowIndex = globalNodeOrdinalsFC[localCellToNodeFC(k,i)];
#ifdef HAVE_MPI
        double uh1 = uCoeff.Values()[rowIndex];
#else
        double uh1 = femCoefficients.Values()[rowIndex];
#endif
        uApprox += uh1*uhGValsTrans(0,i,nPt);
        graduApprox1 += uh1*uhGradsTrans(0,i,nPt,0);
        graduApprox2 += uh1*uhGradsTrans(0,i,nPt,1);
        graduApprox3 += uh1*uhGradsTrans(0,i,nPt,2);
      }

      // REMOVE
      //std::cout << " exact solution = " << uExact << "\n";
      //std::cout << " FEM solution   = " << uApprox << "\n\n";



      // evaluate the error at cubature points
      Linferr = max(Linferr, abs(uExact - uApprox));

      L2errElem+=(uExact - uApprox)*(uExact - uApprox)*worksetCubWeightsE(0,nPt);
      H1errElem+=((graduExact1 - graduApprox1)*(graduExact1 - graduApprox1))
        *worksetCubWeightsE(0,nPt);
      H1errElem+=((graduExact2 - graduApprox2)*(graduExact2 - graduApprox2))
        *worksetCubWeightsE(0,nPt);
      H1errElem+=((graduExact3 - graduApprox3)*(graduExact3 - graduApprox3))
        *worksetCubWeightsE(0,nPt);
    }

    L2err+=L2errElem;
    H1err+=H1errElem;
  }


#ifdef HAVE_MPI
  // sum over all processors
  Comm.SumAll(&L2err,&L2errTot,1);
  Comm.SumAll(&H1err,&H1errTot,1);
  Comm.MaxAll(&Linferr,&LinferrTot,1);
#else
  L2errTot = L2err;
  H1errTot = H1err;
  LinferrTot = Linferr;
#endif

  if (MyPID == 0) {
    std::cout << "\n" << "L2 Error:  " << sqrt(L2errTot) <<"\n";
    std::cout << "H1 Error:  " << sqrt(H1errTot) <<"\n";
    std::cout << "LInf Error:  " << LinferrTot <<"\n";
  }

  // Cleanup
  delete [] myGlobalDofs;

  // reset format state of std::cout
  //   std::cout.copyfmt(oldFormatState);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  exit(0);

}


/**************************************************************************************************
 *                                                                                                *
 *    Problem specification:                                                                      *
 *    Begin user-defined functions for exact solution, diffusion, advection and reaction terms    *
 *                                                                                                *
 **************************************************************************************************/

template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z) {

  // Patch test: tri-linear function is in the FE space and should be recovered
  return 1. + x + y + z + x*y + x*z + y*z + x*y*z;

  // Analytic solution with homogeneous Dirichlet boundary data
  // return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  // return exp(x + y + z)/(1. + x*y + y*z + x*y*z);
}



template<typename Scalar>
void diffusionTensor(Scalar diffusion[][3], const Scalar& x, const Scalar& y, const Scalar& z) {
  diffusion[0][0] = 10.;
  diffusion[0][1] = 0.;
  diffusion[0][2] = 0.;
  //
  diffusion[1][0] = 0.;
  diffusion[1][1] = 10.;
  diffusion[1][2] = 0.;
  //
  diffusion[2][0] = 0.;
  diffusion[2][1] = 0.;
  diffusion[2][2] = 10.;
}



template<typename Scalar>
void advectiveVector(Scalar advection[3], const Scalar& x, const Scalar& y, const Scalar& z){

  advection[0] = g_advection[0];
  advection[1] = g_advection[1];
  advection[2] = g_advection[2];
}


template<typename Scalar>
const Scalar reactionTerm(const Scalar& x, const Scalar& y, const Scalar& z){
  return g_reaction;
}



/**************************************************************************************************
 *                                                                                                *
 *    Definitions of auxiliary functions that require user-defined functions for                  *
 *    exact solution, diffusion, advection and reaction terms                                     *
 *                                                                                                *
 **************************************************************************************************/

template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z) {

  // To enable derivatives of the gradient (i.e., 2nd derivatives of the exact solution) need 2 levels of fad types
  Sacado::Fad::SFad<Scalar,3> fad_x = x;
  Sacado::Fad::SFad<Scalar,3> fad_y = y;
  Sacado::Fad::SFad<Scalar,3> fad_z = z;
  Sacado::Fad::SFad<Scalar,3> u;

  // Indicate the independent variables
  fad_x.diff(0,3);
  fad_y.diff(1,3);
  fad_z.diff(2,3);

  u = exactSolution(fad_x, fad_y, fad_z);

  gradExact[0] = u.dx(0);
  gradExact[1] = u.dx(1);
  gradExact[2] = u.dx(2);
}



template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z){

  Scalar u;
  Scalar grad_u[3];
  Scalar flux[3] = {0.0, 0.0, 0.0};
  Scalar diffusion[3][3];
  Scalar advection[3];
  Scalar f = 0.;

  // Indicate the independent variables
  x.diff(0,3);
  y.diff(1,3);
  z.diff(2,3);

  // Get exact solution and its gradient
  u = exactSolution(x, y, z);
  exactSolutionGrad(grad_u, x, y, z);

  // Get diffusion tensor and advective vector
  diffusionTensor<Scalar>(diffusion, x, y, z);
  advectiveVector<Scalar>(advection, x, y, z);

  // Compute total flux = (K.grad u - b.u)
  for(int i = 0; i < 3; i++){

    // Initialize by the advective flux
    flux[i] = -advection[i]*u;
    //flux[i] = 0.;

    // Add diffusive flux
    for(int j = 0; j < 3; j++){
      flux[i] += diffusion[i][j]*grad_u[j];
    }
  }

  // Compute source term (right hand side): f = -div(flux) + reactionCoeff*u
  f = -(flux[0].dx(0) + flux[1].dx(1) + flux[2].dx(2)) +  reactionTerm<Scalar>(x, y, z)*u;

  return f;
}



/**************************************************************************************************
*                                                                                                *
*    Evaluation methods accepting multidimensional arrays of points                              *
*                                                                                                *
**************************************************************************************************/

template<class ArrayOut, class ArrayIn>
void evaluateDiffusionTensor(ArrayOut &        diffTensorValues,
                             const ArrayIn &   evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints        = evaluationPoints.dimension(1);
  int spaceDim         = evaluationPoints.dimension(2);

  double diffusion[3][3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      diffusionTensor<double>(diffusion, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        for(int col = 0; col < spaceDim; col++){
          diffTensorValues(cell, pt, row, col) = diffusion[row][col];
        }
      }
    }
  }
}



template<class ArrayOut, class ArrayIn>
void evaluateAdvectiveVector(ArrayOut &       worksetAdvectionValues,
                             const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);
  int spaceDim  = evaluationPoints.dimension(2);

  double advection[3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      advectiveVector<double>(advection, x, y, z);

      for(int row = 0; row < spaceDim; row++){
          worksetAdvectionValues(cell, pt, row) = advection[row];
      }
    }
  }
}



template<class ArrayOut, class ArrayIn>
void evaluateReactionCoefficient(ArrayOut &       worksetReactionValues,
                                 const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      worksetReactionValues(cell, pt) = reactionTerm<double>(x, y, z);
    }
  }
}



template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      Sacado::Fad::SFad<double,3> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<double,3> y = evaluationPoints(cell, pt, 1);
      Sacado::Fad::SFad<double,3> z = evaluationPoints(cell, pt, 2);

      sourceTermValues(cell, pt) = sourceTerm<Sacado::Fad::SFad<double,3> >(x, y, z).val();
    }
  }
}



template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionValues(cell, pt) = exactSolution<double>(x, y, z);
    }
  }
}



template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);
  int spaceDim  = evaluationPoints.dimension(2);

  double gradient[3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionGrad<double>(gradient, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}


/**************************************************************************************************
 *                                                                                                *
 *    Solver part                                                                                 *
 *                                                                                                *
 **************************************************************************************************/

// Test ML
int TestMultiLevelPreconditioner(char ProblemType[],
                                 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol)
{
  Epetra_MultiVector x(xexact);
  x.PutScalar(0.0);

  Epetra_LinearProblem Problem(&A,&x,&b);
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();

  Epetra_Time Time(A.Comm());

  // =================== //
  // call ML and AztecOO //
  // =================== //

  AztecOO solver(Problem);
  ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 10);

  solver.Iterate(500, 1e-10);

  delete MLPrec;

  uh = *lhs;

  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  double d = 0.0, d_tot = 0.0;
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - xexact[0][i]) * ((*lhs)[0][i] - xexact[0][i]);

  A.Comm().SumAll(&d,&d_tot,1);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A.Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);

  string msg = ProblemType;

  if (A.Comm().MyPID() == 0) {
    cout << msg << endl << "......Using " << A.Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;

  return( solver.NumIters() );

}


/***************************************************************************************************
  *                                                                                                *
  *    Pamgen mesh                                                                                 *
  *                                                                                                *
  **************************************************************************************************/

template<class Scalar>
void getPamgenMesh(FieldContainer<Scalar>    & localNodeCoordsFC,
                   FieldContainer<long long> & localCellToNodeFC,
                   FieldContainer<int> &       nodeOnBoundaryFC,
                   FieldContainer<bool> &      nodeIsOwnedFC,
                   FieldContainer<long long> & globalNodeIdsFC,
                   const std::string &         meshInput,
                   const int &                 procRank,
                   const int &                 numProcs,
                   const Epetra_Comm &         Comm,
                   Epetra_Time &               Time,
                   const string &              message,
                   const int                   verbose) {

  typedef long long                     int128;

  // All cells in a Pamgen grid have Hexahedron<8> type.
  ShardsCellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );
  int spaceDim = hex_8.getDimension();

  if (!Comm.MyPID() && verbose) {
    std::cout << message << "Generating Pamgen mesh ... \n\n";  Time.ResetStartTime();
  }

  // Error flag
  int pamgenError = 0;

  // Generate mesh with Pamgen
  int128 maxInt = 9223372036854775807LL;
  int128 *  node_comm_proc_ids   = NULL;
  int128 *  node_cmap_node_cnts  = NULL;
  int128 *  node_cmap_ids        = NULL;
  int128 ** comm_node_ids        = NULL;
  int128 ** comm_node_proc_ids   = NULL;


  long long cr_result = Create_Pamgen_Mesh(meshInput.c_str(), spaceDim, procRank, numProcs, maxInt);
  TrilinosCouplings::pamgen_error_check(std::cout,cr_result);

  if(!Comm.MyPID() && verbose ) {
    std::cout << message << " Timing of Pamgen tasks: \n\n" ;
    std::cout << "                  Pamgen Setup: " << Time.ElapsedTime() << " sec.\n";
    Time.ResetStartTime();
  }


/***************************************************************************************************
  *                                                                                                *
  *      Process local element blocks: local to the rank processor for which mesh was created      *
  *                                                                                                *
  **************************************************************************************************/

  // Get local mesh data: limited to entities local to the rank processor
  char title[100];
  int128 numDim;
  int128 numNodesLocal;
  int128 numGridCellsLocal;
  int128 numElemBlkLocal;
  int128 numNodeSetsLocal;
  int128 numSideSetsLocal;
  int id = 0;

  pamgenError += im_ex_get_init_l(id, title,
                                  &numDim,          &numNodesLocal,     &numGridCellsLocal,
                                  &numElemBlkLocal, &numNodeSetsLocal,  &numSideSetsLocal);

  // Get local element block data: retrieve element node count, connectivity, etc on this processor
  int128  * block_ids           = new int128   [numElemBlkLocal];
  int128  * nodes_per_element   = new int128   [numElemBlkLocal];
  int128  * element_attributes  = new int128   [numElemBlkLocal];
  int128  * elements            = new int128   [numElemBlkLocal];
  char   ** element_types       = new char   * [numElemBlkLocal];
  int128 ** elmt_node_linkage   = new int128 * [numElemBlkLocal];

  // First get the ids of the element blocks on this processor
  pamgenError += im_ex_get_elem_blk_ids_l(id, block_ids);

  // Then get the mesh data for each one of the blocks on this processor
  for(int128 i = 0; i < numElemBlkLocal; i ++){
    element_types[i] = new char [MAX_STR_LENGTH + 1];
    pamgenError += im_ex_get_elem_block_l(id,
                                          block_ids[i],
                                          element_types[i],
                                          (int128*) & (elements[i]),
                                          (int128*) & (nodes_per_element[i]),
                                          (int128*) & (element_attributes[i]));
  }

  // Get connectivity of the elements in each block
  for(int128 b = 0; b < numElemBlkLocal; b++){
    elmt_node_linkage[b] =  new int128 [nodes_per_element[b]*elements[b]];
    pamgenError += im_ex_get_elem_conn_l(id,
                                         block_ids[b],
                                         elmt_node_linkage[b]);
  }

  // Copy element-to-node connectivity to multi-dimensional array: collect connectivity from all element
  // blocks. This requires all blocks to have the same number of nodes per element!
  int telct = 0;
  int numNodesPerElem = nodes_per_element[0]; // This will not work if blocks have elements with diff. node counts!
  localCellToNodeFC.resize(numGridCellsLocal,numNodesPerElem);

  for(int128 b = 0; b < numElemBlkLocal; b++){
    for(int128 el = 0; el < elements[b]; el++){
      for (int j = 0; j < numNodesPerElem; j++) {
        localCellToNodeFC(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
      }
      telct ++;
    }
  }

  // Read node coordinates for this processor and place in field container
  localNodeCoordsFC.resize(numNodesLocal,spaceDim);

  double * nodeCoordx = new double [numNodesLocal];
  double * nodeCoordy = new double [numNodesLocal];
  double * nodeCoordz = new double [numNodesLocal];

  pamgenError += im_ex_get_coord_l(id,nodeCoordx,nodeCoordy,nodeCoordz);

  for (int i=0; i<numNodesLocal; i++) {
    localNodeCoordsFC(i,0) = nodeCoordx[i];
    localNodeCoordsFC(i,1) = nodeCoordy[i];
    localNodeCoordsFC(i,2) = nodeCoordz[i];
  }
  delete [] nodeCoordx;
  delete [] nodeCoordy;
  delete [] nodeCoordz;


/***************************************************************************************************
  *                                                                                                *
  *      Parallel info: get global mesh data                                                       *
  *                                                                                                *
  **************************************************************************************************/

  // Get global mesh data: retrieve mesh sizing  information for the complete (global) mesh data base
  int128 numNodesGlobal;
  int128 numGridCellsGlobal;
  int128 numElemBlkGlobal;
  int128 numNodeSetsGlobal;
  int128 numSideSetsGlobal;

  pamgenError += im_ne_get_init_global_l(id,
                                         &numNodesGlobal,    &numGridCellsGlobal, &numElemBlkGlobal,
                                         &numNodeSetsGlobal, &numSideSetsGlobal);

  int128 num_internal_nodes;
  int128 num_border_nodes;
  int128 num_external_nodes;
  int128 num_internal_elems;
  int128 num_border_elems;
  int128 num_node_comm_maps;
  int128 num_elem_comm_maps;

  // Get sizing information for processor communication data. This is the  first step in gathering
  // all the information required to construct communication protocols between adjacent regions of decomposed mesh.
  pamgenError += im_ne_get_loadbal_param_l( id,
                                            &num_internal_nodes,
                                            &num_border_nodes,
                                            &num_external_nodes,
                                            &num_internal_elems,
                                            &num_border_elems,
                                            &num_node_comm_maps,
                                            &num_elem_comm_maps,
                                            0/*unused*/ );

  if(num_node_comm_maps > 0){
    node_comm_proc_ids   = new int128  [num_node_comm_maps];
    node_cmap_node_cnts  = new int128  [num_node_comm_maps];
    node_cmap_ids        = new int128  [num_node_comm_maps];
    comm_node_ids        = new int128* [num_node_comm_maps];
    comm_node_proc_ids   = new int128* [num_node_comm_maps];

    int128 *  elem_cmap_ids        = new int128 [num_elem_comm_maps];
    int128 *  elem_cmap_elem_cnts  = new int128 [num_elem_comm_maps];

    // Get communication maps and ids
    if( im_ne_get_cmap_params_l(id,
                                node_cmap_ids,
                                (int128*)node_cmap_node_cnts,
                                elem_cmap_ids,
                                (int128*)elem_cmap_elem_cnts,
                                0/*not used proc_id*/ ) < 0 ) ++pamgenError;

    for(int128 j = 0; j < num_node_comm_maps; j++) {
      comm_node_ids[j]       = new int128 [node_cmap_node_cnts[j]];
      comm_node_proc_ids[j]  = new int128 [node_cmap_node_cnts[j]];

      // Get node communication map
      if( im_ne_get_node_cmap_l(id,
                                node_cmap_ids[j],
                                comm_node_ids[j],
                                comm_node_proc_ids[j],
                                0/*not used proc_id*/ ) < 0 ) ++pamgenError;
      node_comm_proc_ids[j] = comm_node_proc_ids[j][0];
    }
    delete [] elem_cmap_ids;
    delete [] elem_cmap_elem_cnts;
  }

  if(!Comm.MyPID() && verbose) {
    std::cout << "                  Mesh Queries: " << Time.ElapsedTime() << " sec.\n";
    Time.ResetStartTime();
  }


  //Calculate global node ids
  int128 * globalNodeIds = new int128[numNodesLocal];
  bool *   nodeIsOwned   = new bool[numNodesLocal];

  calc_global_node_ids(globalNodeIds,
                       nodeIsOwned,
                       numNodesLocal,
                       num_node_comm_maps,
                       node_cmap_node_cnts,
                       node_comm_proc_ids,
                       comm_node_ids,
                       procRank);

  if(!Comm.MyPID() && verbose) {
    std::cout << "  Calculating Global Node Nums: " << Time.ElapsedTime() << " sec.\n";
    Time.ResetStartTime();
  }


  // Container indicating whether a node is owned or ghost (1-yes 0-no)
  nodeIsOwnedFC.resize(numNodesLocal);
  globalNodeIdsFC.resize(numNodesLocal);

  for(int128 i = 0; i < numNodesLocal; i++){
    nodeIsOwnedFC(i) = nodeIsOwned[i];
    globalNodeIdsFC(i) = globalNodeIds[i];
  }

  if(!Comm.MyPID() && verbose) {
    std::cout << "        Processing Owned Nodes: " << Time.ElapsedTime() << " sec.\n";
    Time.ResetStartTime();
  }


  // Get boundary (side set) information
  int128 * sideSetIds = new int128 [numSideSetsLocal];
  int128 numSidesInSet;
  int128 numDFinSet;

  pamgenError += im_ex_get_side_set_ids_l(id,sideSetIds);


  // Container indicating whether a node is on the boundary (1-yes 0-no).
  nodeOnBoundaryFC.resize(numNodesLocal);

  for (int i = 0; i < numSideSetsLocal; i++) {
    pamgenError += im_ex_get_side_set_param_l(id, sideSetIds[i], &numSidesInSet, &numDFinSet);
    if (numSidesInSet > 0){
      int128 * sideSetElemList = new int128 [numSidesInSet];
      int128 * sideSetSideList = new int128 [numSidesInSet];
      im_ex_get_side_set_l(id,sideSetIds[i],sideSetElemList,sideSetSideList);
      for (int j=0; j<numSidesInSet; j++) {

        int sideNode0 = hex_8.getNodeMap(2, sideSetSideList[j]-1, 0);
        int sideNode1 = hex_8.getNodeMap(2, sideSetSideList[j]-1, 1);
        int sideNode2 = hex_8.getNodeMap(2, sideSetSideList[j]-1, 2);
        int sideNode3 = hex_8.getNodeMap(2, sideSetSideList[j]-1, 3);

        nodeOnBoundaryFC(localCellToNodeFC(sideSetElemList[j]-1,sideNode0))=1;
        nodeOnBoundaryFC(localCellToNodeFC(sideSetElemList[j]-1,sideNode1))=1;
        nodeOnBoundaryFC(localCellToNodeFC(sideSetElemList[j]-1,sideNode2))=1;
        nodeOnBoundaryFC(localCellToNodeFC(sideSetElemList[j]-1,sideNode3))=1;
      }
      delete [] sideSetElemList;
      delete [] sideSetSideList;
    }
  }
  delete [] sideSetIds;

  if(!Comm.MyPID() && verbose) {
    std::cout << "Processing Boundary Conditions: " << Time.ElapsedTime() << " sec.\n";
    Time.ResetStartTime();
  }

  // Print mesh information
  if (!Comm.MyPID()){
    std::cout << "     Number of Global Elements: " << numGridCellsGlobal << " \n";
    std::cout << "        Number of Global Nodes: " << numNodesGlobal << " \n\n";
  }


#ifdef DUMP_DATA
  // Print coords

  std::stringstream fname;
  fname << "nodeCoords";
  fname << Comm.MyPID() << ".dat";
  FILE *f=fopen(fname.str().c_str(),"w");
  for (int i=0; i<numNodesLocal; i++) {
    if (nodeIsOwned[i]) {
      fprintf(f,"%22.16e %22.16e %22.16e\n", localNodeCoordsFC(i,0),localNodeCoordsFC(i,1),localNodeCoordsFC(i,2));
    }
  }
  fclose(f);

  // Output element to node connectivity
  std::stringstream efname;
  efname << "elem2node";
  efname << Comm.MyPID() << ".dat";
  ofstream el2nout(efname.str().c_str());
  for (int i=0; i<numGridCellsLocal; i++) {
    for (int m=0; m<numNodesPerElem; m++) {
      el2nout << globalNodeIds[localCellToNodeFC(i,m)] << "  ";
    }
    el2nout << "\n";
  }
  el2nout.close();
#endif


  // Delete mesh & clean up
  for(int128 b = 0; b < numElemBlkLocal; b++){
    delete [] elmt_node_linkage[b];
    delete [] element_types[b];
  }
  delete [] nodeIsOwned;
  delete [] block_ids;
  delete [] nodes_per_element;
  delete [] element_attributes;
  delete [] element_types;
  delete [] elmt_node_linkage;
  delete [] elements;
  delete [] globalNodeIds;

  if(num_node_comm_maps > 0){
    delete [] node_comm_proc_ids;
    delete [] node_cmap_node_cnts;
    delete [] node_cmap_ids;
    for(int128 i=0;i<num_node_comm_maps;i++){
      delete [] comm_node_ids[i];
      delete [] comm_node_proc_ids[i];
    }
    delete [] comm_node_ids;
    delete [] comm_node_proc_ids;
  }

  Delete_Pamgen_Mesh();

  TEUCHOS_TEST_FOR_EXCEPTION( !(pamgenError == 0), std::runtime_error, " Pamgen error... Exiting!");
}// *** PAMGEN wrapper ***

/**************************************************************************************************
 *                                                                                                *
 *    Argument input                                                                              *
 *                                                                                                *
 **************************************************************************************************/
void getInputArguments(Teuchos::ParameterList &  inputMeshList,
                       std::string            &  meshInput,
                       Teuchos::ParameterList &  inputSolverList,
                       int &                     argc,
                       char *                    argv[],
                       const Epetra_Comm &       Comm,
                       const std::string &       inputMeshFile,
                       const std::string &       inputSolverFile){

  // Command line for xml file, otherwise use default
  std::string   xmlMeshInFileName, xmlSolverInFileName;

  //Check number of arguments
  if (argc > 3) {
    std::cout << "\n>>> ERROR: Invalid number of arguments.\n\n";
    std::cout << "Usage:\n\n";
    std::cout << "  ./TrilinosCouplings_examples_scaling_Example_StabilizedADR.exe [meshfile.xml] [solver.xml]\n\n";
    std::cout << "   meshfile.xml(optional) - xml file with description of Pamgen mesh\n\n";
    std::cout << "   solver.xml(optional)   - xml file with ML solver options\n\n";
    TEUCHOS_TEST_FOR_EXCEPTION( argc == 3, std::invalid_argument, " ...exiting!");
  }

  // Solver is specified
  if(argc >= 3) xmlSolverInFileName = string(argv[2]);
  else          xmlSolverInFileName = inputSolverFile;

  // Custom mesh file is specified
  if(argc >= 2) xmlMeshInFileName   = string(argv[1]);
  else          xmlMeshInFileName   = inputMeshFile;

  if(xmlMeshInFileName.length()) {
    if (!Comm.MyPID()) {
      std::cout << "\nReading parameter list from the XML file \""<< xmlMeshInFileName << "\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile (xmlMeshInFileName, Teuchos::ptr (&inputMeshList));

    if (!Comm.MyPID()) {
      inputMeshList.print(std::cout, 2, true, true);
      std::cout << "\n";
    }
  }
  else{
    std::cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, " ...exiting!");
  }

  if(xmlSolverInFileName.length()) {
    if (!Comm.MyPID()){
      std::cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile (xmlSolverInFileName, Teuchos::ptr (&inputSolverList));
  }
  else if (!Comm.MyPID()) {
    std::cout << "Using default solver values ..." << std::endl;
  }

  // Get pamgen mesh definition
  meshInput = Teuchos::getParameter<std::string>(inputMeshList,"meshInput");

}















