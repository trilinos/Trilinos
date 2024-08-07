// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/********************************************************************************/
/****************** Solution of Advection Diffusion Equation ********************/
/********************** Using CVFEM with Multi-dimensional **********************/
/************************ Scharfetter-Gummel Updwinding *************************/
/********************************************************************************/

/** \file   example_CVFEM.cpp
    \brief  Example solution of an Advection Diffusion equation on a quadrilateral
            or triangular mesh using the CVFEM.

    \verbatim

     Advection diffusion system:

          - div (epsilon grad phi - u phi) = f in Omega
                                       phi = g on Gamma

       where
             u is the advection velocity
             epsilon is the diffusion coefficient
             f is a given source term


     Corresponding discrete linear system for nodal coefficients(x):

                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector

    \endverbatim

    \author K. Peterson

    \remark Usage:
    \code   ./TrilinosCouplings_Example_CVFEM.exe <input.xml>  \endcode

    \remark Example requires an xml input file with mesh and solver settings.


    NOTE: This problem does not currently work on more than one core (csiefer@sandia.gov 3/17/16)
*/

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HCURL_TRI_I1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CubatureControlVolume.hpp"
#include "Intrepid_CubatureControlVolumeSide.hpp"
#include "Intrepid_Utils.hpp"

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

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// AztecOO includes
#include "AztecOO.h"

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

// Sacado includes
#include "Sacado.hpp"

/*** Uncomment if you would like output data for plotting ***/
//#define DUMP_DATA

/**********************************************************************************/
/**************************** FUNCTION DECLARATIONS *******************************/
/**********************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point

    \return Value of the exact solution at (x,y)
 */
template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y);

/** \brief  User-defined advection velocity.

    \param  advVel      [out]   advection velocity evaluated at (x,y)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
 */
template<typename Scalar>
void advectionVelocity(Scalar advVel[2], const Scalar& x, const Scalar& y, const std::string problem);

/** \brief  User-defined diffusivity

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  epsilon     [in]    diffusion coefficient from xml imput file

    \return Value of the diffusivity at (x,y)
 */
template<typename Scalar1, typename Scalar2>
const Scalar1 diffusivity(const Scalar1& x, const Scalar1& y, const Scalar2& epsilon, const bool variableEpsilon);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
 */
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[2], const Scalar& x, const Scalar& y);

/** \brief Computes source term: f = -div(J_n).  Requires user-defined exact solution
           and material parameters.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  epsilon    [in]    diffusion coefficient

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y)
 */
template<typename Scalar1, typename Scalar2>
const Scalar1 sourceTerm(Scalar1& x,    Scalar1& y, Scalar2& epsilon,
                        const bool variableEpsilon, const std::string problem);


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

/** \brief Create quadrilateral or triangle mesh on the domain [0,1]x[0,1]

    \param elemToNode                [out]      Array with element to node mapping
    \param nodeCoords                [out]      Array with nodal coordinates
    \param bcLeftId                  [out]      Array with ids of left boundary nodes,
    \param bcRightId                 [out]      Array with ids of right boundary nodes
    \param bcTopId                   [out]      Array with ids of top boundary nodes
    \param bcBotId                   [out]      Array with ids of bottom boundary nodes
    \param meshType                   [in]      Mesh type (quad or tri)
    \param meshSize                   [in]      Number of elements in each direction
*/
template<class ArrayOut1, class ArrayOut2, class Scalar>
void createMesh(ArrayOut1 &  elemToNode, ArrayOut2 & nodeCoords,
                ArrayOut1 &  bcLeftId, ArrayOut1 & bcRightId,
                ArrayOut1 &  bcTopId, ArrayOut1 & bcBotId,
                const std::string &  meshType, const Scalar & meshSize);

/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  int numProcs=1;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  numProcs=mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  Epetra_Time Time(Comm);

  bool showTiming = 0;

#ifdef HAVE_MPI
  if (numProcs > 1) {
     if (MyPID == 0) {
       std::cout <<"\n>>> ERROR: Example will only run on a single processor for now. \n\n";
     }
     exit(1);
  }
#endif

/**********************************************************************************/
/*********************************** XML INPUT ************************************/
/**********************************************************************************/

   // get xml file from command line if provided, otherwise use default
    std::string  xmlInFileName;
    if(argc>=2) xmlInFileName=std::string(argv[1]);
    else xmlInFileName="CVFEM.xml";

   // default values
    std::string meshType = "quad";
    std::string meshMotion = "none";
    int meshSize = 10;
    std::string problem = "manufacturedSoln";
    std::string stabilization = "type1";
    double epsilon = 0.001;
    bool variableEpsilon = 0;

   // read xml file into parameter list
    Teuchos::ParameterList inputList;

     if(xmlInFileName.length()) {

       if (MyPID == 0)
         std::cout << "\nReading parameter list from \""<<xmlInFileName<<"\" \n\n";

       Teuchos::updateParametersFromXmlFile(xmlInFileName, Teuchos::ptr (&inputList));

       // get mesh input
       Teuchos::ParameterList & meshList=inputList.sublist("Mesh Input");
       meshType   = Teuchos::getParameter<std::string>(meshList,"meshType");
       meshMotion = Teuchos::getParameter<std::string>(meshList,"meshMotion");
       meshSize   = Teuchos::getParameter<int>(meshList,"meshSize");
       meshMotion = Teuchos::getParameter<std::string>(meshList,"meshMotion");

       // get problem definition
       problem = Teuchos::getParameter<std::string>(inputList,"problem");

       // get stabilization type
       stabilization = Teuchos::getParameter<std::string>(inputList,"stabilization");

       // get solver inputs
       Teuchos::ParameterList inputSolverList = inputList.sublist("ML Input") ;

       // get physical parameter values
       Teuchos::ParameterList & physicsList=inputList.sublist("Physics Input");
       epsilon = physicsList.get("epsilon",0.001);
       variableEpsilon = Teuchos::getParameter<bool>(physicsList,"variableEpsilon");

        std::cout <<"       meshType         = " << meshType << " \n";
        std::cout <<"       meshSize         = " << meshSize << " \n";
        std::cout <<"       meshMotion       = " << meshMotion << " \n";
        std::cout <<"       problem          = " << problem << " \n";
        std::cout <<"       stabilization    = " << stabilization<< " \n";
        std::cout <<"       epsilon          = " << epsilon << " \n";
        std::cout <<"       variableEpsilon  = " << variableEpsilon << " \n\n";

     }
     else {

       if (MyPID == 0){
        std::cout << "\nSetting default parameters: \n";
        std::cout <<"       meshType         = " << meshType << " \n";
        std::cout <<"       meshSize         = " << meshSize << " \n";
        std::cout <<"       meshMotion       = " << meshMotion << " \n";
        std::cout <<"       problem          = " << problem << " \n";
        std::cout <<"       stabilization    = " << stabilization<< " \n";
        std::cout <<"       epsilon          = " << epsilon << " \n";
        std::cout <<"       variableEpsilon  = " << variableEpsilon << " \n\n";
       }

    }

      TEUCHOS_TEST_FOR_EXCEPTION( ( (stabilization != "type1")             &&
                               (stabilization != "type2")  && (stabilization != "type3")),
                                std::invalid_argument,
                               "Unknown stabilization type./n/n");

      TEUCHOS_TEST_FOR_EXCEPTION( ( (problem != "manufacturedSoln") && (problem != "skewAdvection") &&
                               (problem != "horizAdvection")  && (problem != "doubleGlazing")),
                                std::invalid_argument,
                               "Unknown problem type./n/n");

/**********************************************************************************/
/*********************************** BUILD MESH ***********************************/
/**********************************************************************************/

   // 2-D meshes only for now
    int spaceDim = 2;

   // For selected mesh type define element to node mapping and element node coordinates
    Intrepid::FieldContainer<int> elemToNode;
    Intrepid::FieldContainer<double> nodeCoords;

   // Also define arrays for use in defining boundary conditions
    Intrepid::FieldContainer<int> bcLeftId;
    Intrepid::FieldContainer<int> bcRightId;
    Intrepid::FieldContainer<int> bcTopId;
    Intrepid::FieldContainer<int> bcBotId;

   // Create mesh on square domain [0,1]x[0,1]
    createMesh(elemToNode, nodeCoords, bcLeftId,
               bcRightId, bcTopId, bcBotId,
               meshType, meshSize);

    int numNodes = nodeCoords.dimension(0);
    int numElems = elemToNode.dimension(0);

    Intrepid::FieldContainer<int> nodeOnBoundary(numNodes);
    for (int i = 0; i<bcLeftId.dimension(0); i++){
       nodeOnBoundary(bcLeftId(i)) = 1;
    }
    for (int i = 0; i<bcRightId.dimension(0); i++){
       nodeOnBoundary(bcRightId(i)) = 1;
    }
    for (int i = 0; i<bcTopId.dimension(0); i++){
       nodeOnBoundary(bcTopId(i)) = 1;
    }
    for (int i = 0; i<bcBotId.dimension(0); i++){
       nodeOnBoundary(bcBotId(i)) = 1;
    }

   // For now assume regular mesh
     double hx = nodeCoords(elemToNode(0,1),0) - nodeCoords(elemToNode(0,0),0);
     double hy = nodeCoords(elemToNode(0,2),1) - nodeCoords(elemToNode(0,1),1);

   // Modify nodal coordinates, if necessary
      for (int inode = 0; inode < numNodes; inode++) {

         double x = nodeCoords(inode,0);
         double y = nodeCoords(inode,1);

         // randomly perturbed mesh
         if (meshMotion == "random"){
           if (!nodeOnBoundary(inode)) {
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            nodeCoords(inode,0) = x + 0.25 * hx * rx;
            nodeCoords(inode,1) = y + 0.25 * hy * ry;
           } // end if !node on boundary
         }

         // Wavy grid motion 1
         else if (meshMotion == "wave"){
           double alpha = 0.5/5.0;
           nodeCoords(inode,0) = x + alpha*sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
           nodeCoords(inode,1) = y + alpha*sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
         }

         // Tensor product grid
         else if (meshMotion == "tensor"){
           double gamma = 0.1;
           double alpha = sin(4.0*M_PI*gamma)/2.0;
           nodeCoords(inode,0) = (1.0-alpha)*x + alpha*x*x*x;
           nodeCoords(inode,1) = (1.0-alpha)*y + alpha*y*y*y;
         }

       } // end node loop


    if(showTiming) {std::cout << "Read/Create mesh                            "
                   << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

/**********************************************************************************/
/********************************SET CELL TOPOLOGY ********************************/
/**********************************************************************************/

    // For now set cell topology based on number of nodes per element.
    //  WARNING: only works for 2-D 1st degree elements.

    Teuchos::RCP<shards::CellTopology> cellType;

    int numNodesPerElem = elemToNode.dimension(1);

    TEUCHOS_TEST_FOR_EXCEPTION( ( (numNodesPerElem != 4) && (numNodesPerElem != 3) ),
                                std::invalid_argument,
                               "Unknown number of nodes per element for cell topology selection. Please use 3 for Triangle or 4 for Quadrilateral.");

    const CellTopologyData &myCellData =
            (numNodesPerElem == 4) ? *shards::getCellTopologyData<shards::Quadrilateral<4> >() : *shards::getCellTopologyData<shards::Triangle<3> >();

//    shards::CellTopology cellType(&myCellData);

    cellType = Teuchos::rcp(new shards::CellTopology(&myCellData));

    int numEdgesPerElem = cellType->getEdgeCount();

    if(showTiming) {std::cout << "Get cell topology                           "
                   << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

/**********************************************************************************/
/******************* GET CONTROL VOLUME INTEGRATION POINTS ************************/
/**********************************************************************************/

    //  Loop over cells to fill cell coordinate array
    Intrepid::FieldContainer<double> cellCoords(numElems, numNodesPerElem, spaceDim);

    // loop over elements
    for (int ielem = 0; ielem < numElems; ielem++) {

      // loop over nodes and fill cell coordinates array
       for (int inode = 0; inode < numNodesPerElem; inode++) {

         cellCoords(ielem, inode, 0) = nodeCoords(elemToNode(ielem,inode),0);
         cellCoords(ielem, inode, 1) = nodeCoords(elemToNode(ielem,inode),1);

       } // end node loop

    } // end cell loop

    Teuchos::RCP<Intrepid::Cubature<double,Intrepid::FieldContainer<double>  > > controlVolCub;
    Teuchos::RCP<Intrepid::Cubature<double,Intrepid::FieldContainer<double>  > > controlVolSideCub;

    // cubature rule for subcontrol volume points
    controlVolCub  = Teuchos::rcp(new Intrepid::CubatureControlVolume<double,Intrepid::FieldContainer<double>,Intrepid::FieldContainer<double> >(cellType));
    int numVolPoints = controlVolCub->getNumPoints();

    // cubature rule for points along control volume edges
    controlVolSideCub  = Teuchos::rcp(new Intrepid::CubatureControlVolumeSide<double,Intrepid::FieldContainer<double>,Intrepid::FieldContainer<double> >(cellType));
    int numSidePoints = controlVolSideCub->getNumPoints();

    // get volume cubature points and weights - points are in physical coordinates!
    Intrepid::FieldContainer<double> volCubPoints(numElems,numVolPoints,spaceDim);
    Intrepid::FieldContainer<double> volCubWeights(numElems,numVolPoints);
    controlVolCub->getCubature(volCubPoints,volCubWeights,cellCoords);

    // get side cubature points and weights - points are in physical coordinates!
    // also note that the weights are vectors, it is the side normal weighted by the side length!
    Intrepid::FieldContainer<double> sideCubPoints(numElems,numSidePoints,spaceDim);
    Intrepid::FieldContainer<double> sideCubWeights(numElems,numSidePoints,spaceDim);
    controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);

    if(showTiming) {std::cout << "Get control volume cubature                 "
                   << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

/**********************************************************************************/
/************************ GET JACOBIANS AT CUBATURE POINTS ************************/
/**********************************************************************************/

    // get cubature points in reference space
    Intrepid::FieldContainer<double> refVolPoints(numElems,numVolPoints,spaceDim);
    Intrepid::FieldContainer<double> refSidePoints(numElems,numSidePoints,spaceDim);
    Intrepid::CellTools<double>::mapToReferenceFrame(refVolPoints, volCubPoints, cellCoords,*cellType);
    Intrepid::CellTools<double>::mapToReferenceFrame(refSidePoints, sideCubPoints, cellCoords,*cellType);

    // get Jacobian
    Intrepid::FieldContainer<double> cellJacobianVol(numElems,numVolPoints,spaceDim,spaceDim);
    Intrepid::FieldContainer<double> cellJacobianSide(numElems,numSidePoints,spaceDim,spaceDim);
    Intrepid::FieldContainer<double> cellJacobInvVol(numElems,numVolPoints,spaceDim,spaceDim);
    Intrepid::FieldContainer<double> cellJacobInvSide(numElems,numSidePoints,spaceDim,spaceDim);
    Intrepid::CellTools<double>::setJacobian(cellJacobianVol, refVolPoints, cellCoords,*cellType,-1);
    Intrepid::CellTools<double>::setJacobianInv(cellJacobInvVol, cellJacobianVol );
    Intrepid::CellTools<double>::setJacobian(cellJacobianSide, refSidePoints, cellCoords,*cellType,-1);
    Intrepid::CellTools<double>::setJacobianInv(cellJacobInvSide, cellJacobianSide );


    if(showTiming) {std::cout << "Get cell jacobians                          "
                   << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

/**********************************************************************************/
/*********************************** GET BASIS ************************************/
/**********************************************************************************/

   // Choose appropriate HGRAD and HCURL bases based on cell topology
   // Evaluate basis at cubature points in cell loop because points are in physical space
    Teuchos::RCP<Intrepid::Basis<double, Intrepid::FieldContainer<double> > >  HGradBasis;
    Teuchos::RCP<Intrepid::Basis<double, Intrepid::FieldContainer<double> > >  HCurlBasis;

    switch (cellType->getKey()) {

       case shards::Triangle<3>::key:
         HGradBasis = Teuchos::rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, Intrepid::FieldContainer<double> > );
         HCurlBasis = Teuchos::rcp(new Intrepid::Basis_HCURL_TRI_I1_FEM<double, Intrepid::FieldContainer<double> > );
         break;

       case shards::Quadrilateral<4>::key:
         HGradBasis = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, Intrepid::FieldContainer<double> > );
         HCurlBasis = Teuchos::rcp(new Intrepid::Basis_HCURL_QUAD_I1_FEM<double, Intrepid::FieldContainer<double> > );
         break;

       default:
         TEUCHOS_TEST_FOR_EXCEPTION( ( (cellType->getKey() != shards::Quadrilateral<4>::key)             &&
                               (cellType->getKey() != shards::Triangle<3>::key) ),
                                std::invalid_argument,
                               "Unknown cell topology for basis selection. Please use Triangle_3 or Quadrilateral_4.");
     }

    int numFields = HGradBasis->getCardinality();
    int numFieldsEdge = HCurlBasis->getCardinality();


    if(showTiming) {std::cout << "Getting basis                               "
                 << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}


/**********************************************************************************/
/********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
/**********************************************************************************/

  // For now assume local nodes = global nodes (and assume DOFs = nodes)

    Epetra_Map   globalMapNodes(numNodes, 0, Comm);
    Epetra_FECrsMatrix OpMatrix(Copy, globalMapNodes, numFields);
    Epetra_FEVector rhsVector(globalMapNodes);

    if(showTiming) {std::cout <<  "Build global maps                           "
                 << Time.ElapsedTime() << " sec \n";  Time.ResetStartTime();}


/**********************************************************************************/
/**** PUT COORDINATES AND NODAL VALUES IN ARRAYS FOR OUTPUT (FOR PLOTTING ONLY) ***/
/**********************************************************************************/

  // Put coordinates in multivector for output
     Epetra_MultiVector nCoord(globalMapNodes,2);

  // Put element to node mapping in multivector for output
     Epetra_Map   globalMapElem(numElems, 0, Comm);
     Epetra_MultiVector elem2node(globalMapElem,numNodesPerElem);

    // Loop over elements
     for (int j = 0; j < numElems; j++) {

        for (int i = 0; i < numNodesPerElem; i++) {
          int inode = elemToNode(j,i);
          elem2node[i][j] = inode;
          nCoord[0][inode] = nodeCoords(inode,0);
          nCoord[1][inode] = nodeCoords(inode,1);
        }

     } // end loop over elements

    // output multivectors
#ifdef DUMP_DATA
     EpetraExt::MultiVectorToMatrixMarketFile("elem2node.dat",elem2node,0,0,false);
     EpetraExt::MultiVectorToMatrixMarketFile("coords.dat",nCoord,0,0,false);
#endif

    if(showTiming) {Time.ResetStartTime();}

/**********************************************************************************/
/************************** DIRICHLET BC SETUP ************************************/
/**********************************************************************************/

 // Vector for use in applying BCs
    Epetra_MultiVector v(globalMapNodes,true);
    v.PutScalar(0.0);

    int numBCNodes =  bcLeftId.dimension(0) + bcTopId.dimension(0)
                     + bcRightId.dimension(0) + bcBotId.dimension(0);

    int * bcNodeVec = new int [numBCNodes];

     int bcNodeCount = 0;

    // Loop over right boundary
     for (int i = 0; i < bcRightId.dimension(0); i++) {

        int bcNodeId = bcRightId(i);
        bcNodeVec[bcNodeCount] = bcNodeId;

        // get value of exact solution on boundary
        double x  = nodeCoords(bcNodeId,0);
        double y  = nodeCoords(bcNodeId,1);

        v[0][bcRightId(i)]=1.0;

        if (problem == "manufacturedSoln")
             v[0][bcNodeId]=exactSolution(x,y);
        if (problem == "horizAdvection")
             v[0][bcNodeId]=0.0;

        bcNodeCount ++;

     } // end loop over right boundary nodes

    // Loop over left boundary
     for (int i = 0; i < bcLeftId.dimension(0); i++) {

        int bcNodeId = bcLeftId(i);

        bcNodeVec[bcNodeCount] = bcNodeId;

        // get value of exact solution on boundary
        double x  = nodeCoords(bcNodeId,0);
        double y  = nodeCoords(bcNodeId,1);

        v[0][bcNodeId]=0.0;

        if (problem == "manufacturedSoln")
             v[0][bcNodeId]=exactSolution(x,y);
        if (problem == "horizAdvection")
             v[0][bcNodeId]=1.0;

        bcNodeCount ++;

     } // end loop over left boundary nodes

    // Loop over top boundary
     for (int i = 0; i < bcTopId.dimension(0); i++) {

        int bcNodeId = bcTopId(i);
        bcNodeVec[bcNodeCount] = bcNodeId;

        // get value of exact solution on boundary
        double x  = nodeCoords(bcNodeId,0);
        double y  = nodeCoords(bcNodeId,1);

        v[0][bcNodeId]=0.0;

        if (problem == "manufacturedSoln")
             v[0][bcNodeId]=exactSolution(x,y);
        if (problem == "horizAdvection")
             v[0][bcNodeId]=1.0;

        bcNodeCount ++;

     } // end loop over top boundary nodes

    // Loop over bottom boundary
     for (int i = 0; i < bcBotId.dimension(0); i++) {

        int bcNodeId = bcBotId(i);
        bcNodeVec[bcNodeCount] = bcNodeId;

        // get value of exact solution on boundary
        double x  = nodeCoords(bcNodeId,0);
        double y  = nodeCoords(bcNodeId,1);

        v[0][bcNodeId]=0.0;

        if (problem == "manufacturedSoln")
             v[0][bcNodeId]=exactSolution(x,y);
        if (problem == "horizAdvection")
             v[0][bcNodeId]=1.0;
        if (problem == "skewAdvection" && x > 0.5)
             v[0][bcNodeId]=1.0;

        bcNodeCount ++;

     } // end loop over bottom boundary nodes

   if(showTiming) {std::cout << "Get Dirichlet boundary values               "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


/**********************************************************************************/
/*********************** ASSEMBLE MATRIX AND VECTOR *******************************/
/**********************************************************************************/

  // Define desired workset size
  int desiredWorksetSize = 50;
  int numWorksets        = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksets*desiredWorksetSize < numElems) numWorksets += 1;

  // Loop over worksets
  for(int workset = 0; workset < numWorksets; workset++){

    // Compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
     worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;
     worksetSize  = worksetEnd - worksetBegin;


 /**********************************************************************************/
 /*                                   Compute RHS                                  */
 /**********************************************************************************/

     Intrepid::FieldContainer<double> worksetRHS (worksetSize, numNodesPerElem);

     int cellCounter = 0;
     for(int icell = worksetBegin; icell < worksetEnd; icell++){

        // loop over nodes to compute contribution from subcontrol volumes
        // (one subcontrol volume associated with each primary cell node)
         for (int inode = 0; inode < numNodesPerElem; inode++) {
            double x = volCubPoints(icell,inode,0);
            double y = volCubPoints(icell,inode,1);
            double source =  0;
            if (problem == "manufacturedSoln")
            {
                Sacado::Fad::SFad<double,2> xfad = x;
                Sacado::Fad::SFad<double,2> yfad = y;
                source = sourceTerm<Sacado::Fad::SFad<double,2> >(xfad, yfad, epsilon, variableEpsilon, problem).val();
            }
            worksetRHS(cellCounter,inode) = source * volCubWeights(icell,inode);
         }
         cellCounter++;

      } // end cell workset loop

//      if(showTiming) {std::cout << "Compute right-hand side                     "
//                  << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


 /**********************************************************************************/
 /*                   Compute Diffusion, Advection, Stabilization                  */
 /**********************************************************************************/

    // Loop over cells to assemble diffusion, advection, and stabilization terms
     Intrepid::FieldContainer<double> worksetDiffOp(worksetSize, numFields, numFields);
     Intrepid::FieldContainer<double> worksetAdvOp(worksetSize, numFields, numFields);
     Intrepid::FieldContainer<double> worksetStabOp(worksetSize, numFields, numFields);

     cellCounter = 0;
     for(int icell = worksetBegin; icell < worksetEnd; icell++){

       // get subcontrol volume edge midpoints in reference frame
       // needed for evaluation of basis function
        Intrepid::FieldContainer<double> thisCellCoords(1,numNodesPerElem, spaceDim);
        Intrepid::FieldContainer<double> thisRefSidePoint(numSidePoints, spaceDim);
        Intrepid::FieldContainer<double> thiscellJacobInv(1, numSidePoints, spaceDim, spaceDim);
         for (int idim = 0; idim < spaceDim; idim++){
             for (int inode = 0; inode < numNodesPerElem; inode++){
               thisCellCoords(0,inode,idim) = cellCoords(icell,inode,idim);
             }
             for (int ipt = 0; ipt < numSidePoints; ipt++){
               thisRefSidePoint(ipt,idim) = refSidePoints(icell,ipt,idim);
                 for (int jdim = 0; jdim < spaceDim; jdim++){
                    thiscellJacobInv(0,ipt,idim,jdim) = cellJacobInvSide(icell,ipt,idim,jdim);
                 }
             }
          }

       // evaluate nodal basis function gradients at subcontrol volume side midpoints
        Intrepid::FieldContainer<double> nodeBasisGrads(numFields, numSidePoints, spaceDim);
        HGradBasis->getValues(nodeBasisGrads, thisRefSidePoint, Intrepid::OPERATOR_GRAD);

       // evaluate nodal basis function values at subcontrol volume side midpoints
        Intrepid::FieldContainer<double> nodeBasisVals(numFields, numSidePoints);
        HGradBasis->getValues(nodeBasisVals, thisRefSidePoint, Intrepid::OPERATOR_VALUE);

       // evaluate edge basis functions at subcontrol volume side midpoints
        Intrepid::FieldContainer<double> edgeBasisVals(numFieldsEdge, numSidePoints, spaceDim);
        HCurlBasis->getValues(edgeBasisVals, thisRefSidePoint, Intrepid::OPERATOR_VALUE);

       // transform node basis function values to physical coordinates
        Intrepid::FieldContainer<double> nodeBasisValsTrans(1, numFields, numSidePoints);
        Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(nodeBasisValsTrans, nodeBasisVals);

       // transform node basis function gradients to physical coordinates
        Intrepid::FieldContainer<double> nodeBasisGradsTrans(1, numFields, numSidePoints, spaceDim);
        Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>(nodeBasisGradsTrans, thiscellJacobInv, nodeBasisGrads);

       // transform edge basis functions to physical coordinates
        Intrepid::FieldContainer<double> edgeBasisValsTrans(1, numFieldsEdge, numSidePoints, spaceDim);
        Intrepid::FunctionSpaceTools::HCURLtransformVALUE<double>(edgeBasisValsTrans, thiscellJacobInv,
                                             edgeBasisVals);

        // Get velocity and diffusivity at subcontrol volume side midpoints
         Intrepid::FieldContainer<double> sideAdvVel(numEdgesPerElem,spaceDim);
         Intrepid::FieldContainer<double> sideDiffusivity(numEdgesPerElem);
         for (int iedge = 0; iedge < numEdgesPerElem; iedge++) {

            double x = sideCubPoints(icell,iedge,0);
            double y = sideCubPoints(icell,iedge,1);
            double advVel[2];
            advectionVelocity(advVel,x,y,problem);
            sideAdvVel(iedge,0) = advVel[0];
            sideAdvVel(iedge,1) = advVel[1];
            sideDiffusivity(iedge) = diffusivity(x,y,epsilon,variableEpsilon);

         }

         // Assemble cell diffusion operator
         // Number of side cubature points is equal to number of primary cell edges
         for (int iedge = 0; iedge < numSidePoints; iedge++) {

           // get local node Ids for primary cell edge
            int rowNodeId0 = cellType->getNodeMap(1,iedge,0);
            int rowNodeId1 = cellType->getNodeMap(1,iedge,1);

            for (int inode = 0; inode < numFields; inode++){

                int colNodeId = inode;

                double operatorContrib0 = 0.0;
                for (int idim = 0; idim < spaceDim; idim++){
                   operatorContrib0 += nodeBasisGradsTrans(0,inode,iedge,idim) * sideCubWeights(icell,iedge,idim);
                }

                operatorContrib0 *= sideDiffusivity(iedge);
                double operatorContrib1 = -operatorContrib0;

                worksetDiffOp(cellCounter,colNodeId,rowNodeId0) += operatorContrib1;
                worksetDiffOp(cellCounter,colNodeId,rowNodeId1) += operatorContrib0;

             }
         }

         // Assemble cell advection operator
         for (int iedge = 0; iedge < numFieldsEdge; iedge++) {

           // get local element node Ids
            int rowNodeId0 = cellType->getNodeMap(1,iedge,0);
            int rowNodeId1 = cellType->getNodeMap(1,iedge,1);

            for (int inode = 0; inode < numFields; inode++){

                int colNodeId = inode;

                double operatorContrib0 = 0.0;
                for (int idim = 0; idim < spaceDim; idim++){
                   operatorContrib0 += nodeBasisValsTrans(0,inode,iedge) * sideAdvVel(iedge,idim)
                                       * sideCubWeights(icell,iedge,idim);
                }

                double operatorContrib1 = -operatorContrib0;

                worksetAdvOp(cellCounter,colNodeId,rowNodeId0) += operatorContrib0;
                worksetAdvOp(cellCounter,colNodeId,rowNodeId1) += operatorContrib1;

             }
         }

        // Get edge coefficients
         Intrepid::FieldContainer<double> edgeCoef(numFieldsEdge);
         Intrepid::FieldContainer<double> edgeCoefOther(numFieldsEdge);
         for (int iedge = 0; iedge < numFieldsEdge; iedge++) {

           // get local element node Ids
            int rowNodeId0 = cellType->getNodeMap(1,iedge,0);
            int rowNodeId1 = cellType->getNodeMap(1,iedge,1);

           // get coords of nodes
             double x0 = cellCoords(icell,rowNodeId0,0); double y0 = cellCoords(icell,rowNodeId0,1);
             double x1 = cellCoords(icell,rowNodeId1,0); double y1 = cellCoords(icell,rowNodeId1,1);

           // get diffusion coefficient at nodes and average
             double D0 = diffusivity(x0,y0,epsilon,variableEpsilon);
             double D1 = diffusivity(x1,y1,epsilon,variableEpsilon);
             double edgeDiffCoef = (D0 + D1)/2.0;

           // get edge length and tangent
             Intrepid::FieldContainer<double> edgeTan(spaceDim);
             double xlen = x1-x0;
             double ylen = y1-y0;
             double edgeLen = sqrt(xlen*xlen + ylen*ylen);
             edgeTan(0) = xlen/edgeLen;
             edgeTan(1) = ylen/edgeLen;

            // get average velocity and dot with tangent
             Intrepid::FieldContainer<double> edgeVel(spaceDim);
             double advVel0[2]; double advVel1[2];
             advectionVelocity(advVel0,x0,y0,problem);
             advectionVelocity(advVel1,x1,y1,problem);
             edgeVel(0) = (advVel0[0] + advVel1[0])/2.0;
             edgeVel(1) = (advVel0[1] + advVel1[1])/2.0;
             double barEdgeVel = edgeVel(0)*edgeTan(0) + edgeVel(1)*edgeTan(1);

            // compute edge Peclet number
             //double edgeAlpha = barEdgeVel * edgeLen / (2.0*edgeDiffCoef);
             double edgeAlpha = barEdgeVel * edgeLen / (edgeDiffCoef);

            // compute nodal coefficients for this edge
             double edgeCoef0 = 0.0;
             double edgeCoef1 = edgeDiffCoef;
             double edgeCoef2 = edgeDiffCoef;
             double edgeCoef3 = edgeDiffCoef;


             if (abs(edgeAlpha) > 1.0e-10) {
                double coth = 1.0/tanh(edgeAlpha/2.0);
                edgeCoef0 = barEdgeVel*edgeLen/2.0 * coth - edgeDiffCoef;
                edgeCoef1 = barEdgeVel*edgeLen/2.0 * (coth);
                edgeCoef2 = barEdgeVel*edgeLen/2.0 * (coth + 1.0);
                edgeCoef3 = barEdgeVel*edgeLen/2.0 * (coth - 1.0);
             }


            if (stabilization == "type1")
               edgeCoef(iedge) = edgeCoef0;
            if (stabilization == "type2")
               edgeCoef(iedge) = edgeCoef1;
            if (stabilization == "type3")
            {
               edgeCoef(iedge) = edgeCoef2;
               edgeCoefOther(iedge) = -edgeCoef3;
            }

        } // end edge loop

          // TESTING THIS MAY NOT WORK

        // loop over primary cell edges
	for (int iedge = 0; iedge < numEdgesPerElem; iedge++) {

          // get node Ids - these are local element node ids
            int rowNodeId0 = cellType->getNodeMap(1,iedge,0);
            int rowNodeId1 = cellType->getNodeMap(1,iedge,1);

           // loop over subcontrol volume sides
	    for (int jedge = 0; jedge < numEdgesPerElem; jedge++) {

                int colNodeId0 = cellType->getNodeMap(1,jedge,0);
                int colNodeId1 = cellType->getNodeMap(1,jedge,1);

               // get operator contributions
                 double Wdotn = 0.0;
                 for (int idim = 0; idim < spaceDim; idim++) {
                     Wdotn += edgeBasisValsTrans(0,iedge,jedge,idim)*sideCubWeights(icell,jedge,idim);
                 }
                 double operatorContrib = edgeCoef(iedge)*Wdotn;
                 double operatorContrib2 = edgeCoefOther(iedge)*Wdotn;

               // sum into operator matrix
               if (stabilization == "type1" || stabilization == "type2"){

                  worksetStabOp(cellCounter,colNodeId0,rowNodeId0) += operatorContrib;
                  worksetStabOp(cellCounter,colNodeId0,rowNodeId1) += -operatorContrib;
                  worksetStabOp(cellCounter,colNodeId1,rowNodeId0) += -operatorContrib;
                  worksetStabOp(cellCounter,colNodeId1,rowNodeId1) += operatorContrib;
               }
               if (stabilization == "type3"){

                  worksetStabOp(cellCounter,colNodeId0,rowNodeId0) += operatorContrib;
                  worksetStabOp(cellCounter,colNodeId0,rowNodeId1) += operatorContrib2;
                  worksetStabOp(cellCounter,colNodeId1,rowNodeId0) += -operatorContrib;
                  worksetStabOp(cellCounter,colNodeId1,rowNodeId1) += -operatorContrib2;
               }

             } // end loop over edges

           } // end loop over edges

         cellCounter++;

      } // end workset cell loop

 /**********************************************************************************/
 /*                         Assemble into Global Matrix                            */
 /**********************************************************************************/
     cellCounter = 0;
     for(int icell = worksetBegin; icell < worksetEnd; icell++){

         //std::cout << "cellCounter = " << cellCounter <<",  icell = " << icell <<"\n";
         // Loop over matrix row
          for (int cellRow = 0; cellRow < numFields; cellRow++){

             int globalRow = elemToNode(icell,cellRow);
             double sourceTermContribution =  worksetRHS(cellCounter, cellRow);
             rhsVector.SumIntoGlobalValues(1, &globalRow, &sourceTermContribution);

             // Loop over matrix column
             for (int cellCol = 0; cellCol < numFields; cellCol++){

               int globalCol = elemToNode(icell,cellCol);

// TEMP - no stab.
//               double operatorMatrixContribution = worksetDiffOp(cellCounter, cellRow, cellCol)
//                                                     - worksetAdvOp(cellCounter, cellRow, cellCol);



               double operatorMatrixContribution = worksetStabOp(cellCounter, cellRow, cellCol);

               if (stabilization == "type2" || stabilization == "type1" ){
                     operatorMatrixContribution += -worksetAdvOp(cellCounter, cellRow, cellCol);
               }

               if (stabilization == "type1" ){
                     operatorMatrixContribution += worksetDiffOp(cellCounter, cellRow, cellCol);
               }



               OpMatrix.InsertGlobalValues(1, &globalRow, 1, &globalCol, &operatorMatrixContribution);

            }// end cell col loop

         }// end cell row loop

         cellCounter++;

      } // end workset cell loop

    } // end loop over worksets


 /**********************************************************************************/
 /*                              Assemble Global Matrix                            */
 /**********************************************************************************/

    OpMatrix.GlobalAssemble();
    OpMatrix.FillComplete();
    rhsVector.GlobalAssemble();

    if(showTiming) {std::cout << "Global assembly                             "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


/**********************************************************************************/
/************************ ADJUST MATRIX AND RHS FOR BCs ***************************/
/**********************************************************************************/

    // Apply stiffness matrix to v
     Epetra_MultiVector rhsDir(globalMapNodes,true);
     OpMatrix.Apply(v,rhsDir);

    // Update right-hand side
     rhsVector.Update(-1.0,rhsDir,1.0);

    // Loop over boundary nodes and replace rhs values with boundary values
     for (int i = 0; i < numBCNodes; i++) {

        int bcNodeId = bcNodeVec[i];

        rhsVector[0][bcNodeId]=v[0][bcNodeId];

     } // end loop over boundary nodes

     // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
     //  and add one to diagonal.
      ML_Epetra::Apply_OAZToMatrix(bcNodeVec, numBCNodes, OpMatrix);

     //delete [] bcNodeVec;

     if(showTiming) {std::cout << "Adjust global matrix and rhs due to BCs     " << Time.ElapsedTime()
                  << " sec \n"; Time.ResetStartTime();}


#ifdef DUMP_DATA
    // Dump matrices to disk
     EpetraExt::RowMatrixToMatlabFile("op_matrix_cvfem.dat",OpMatrix);
     EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector_cvfem.dat",rhsVector,0,0,false);
#endif

/**********************************************************************************/
/*********************************** SOLVE ****************************************/
/**********************************************************************************/

    // Build ML preconditioner
     Teuchos::ParameterList MLList;
     ML_Epetra::SetDefaults("SA",MLList);
     MLList.set("ML output",1);
     //MLList.set("max levels",10);
     MLList.set("max levels",1);
     MLList.set("coarse: type","Amesos-KLU");
     MLList.set("smoother: type","ILU");
     MLList.set("smoother: sweeps",2);
     MLList.set("coarse: max size",1000);
     ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(OpMatrix, MLList, true);

    // Solve
     Epetra_MultiVector x(rhsVector);
     x.PutScalar(0.0);

     Epetra_LinearProblem Problem(&OpMatrix,&x,&rhsVector);
     Epetra_MultiVector* lhs = Problem.GetLHS();

     AztecOO solver(Problem);
     solver.SetPrecOperator(MLPrec);
     solver.SetAztecOption(AZ_solver, AZ_gmres);
     solver.SetAztecOption(AZ_output, 1);
     solver.Iterate(100, 1e-10);

#ifdef DUMP_DATA
     EpetraExt::MultiVectorToMatrixMarketFile("soln_vector_cvfem.dat",*lhs,0,0,false);
#endif

     delete MLPrec;

     if(showTiming) {std::cout << "Solve linear system                         " << Time.ElapsedTime()
                  << " sec \n"; Time.ResetStartTime();}


/**********************************************************************************/
/**************************** CALCULATE ERROR *************************************/
/**********************************************************************************/

   if (problem == "manufacturedSoln") // perform error calculation
   {

     if (showTiming) {Time.ResetStartTime();}

     double L2err = 0.0;
     double L2errTot = 0.0;
     double H1err = 0.0;
     double H1errTot = 0.0;
     double Linferr = 0.0;
     double LinferrTot = 0.0;

/*
    // NOTE: If using multiple processors need to import solution for error calculation !!
    // Import solution onto current processor
    // FIXME - not efficient
     int numNodesGlobal = globalMapNodes.NumGlobalElements();
     Epetra_Map     solnMap(numNodesGlobal, numNodesGlobal, 0, Comm);
     Epetra_Import  solnImporter(solnMap, globalMapNodes);
     Epetra_Vector  uCoeff(solnMap);
     uCoeff.Import(*lhs, solnImporter, Insert);
*/
    // Get cubature points and weights for error calc (may be different from previous)
     Intrepid::DefaultCubatureFactory<double>  cubFactoryErr;
     int cubDegErr = 3;
     Teuchos::RCP<Intrepid::Cubature<double> > cellCubatureErr = cubFactoryErr.create(*cellType, cubDegErr);
     int cubDimErr       = cellCubatureErr->getDimension();
     int numCubPointsErr = cellCubatureErr->getNumPoints();
     Intrepid::FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
     Intrepid::FieldContainer<double> cubWeightsErr(numCubPointsErr);
     cellCubatureErr->getCubature(cubPointsErr, cubWeightsErr);

    // Evaluate basis values and gradients at cubature points
     Intrepid::FieldContainer<double> basisValsErr(numFields, numCubPointsErr);
     Intrepid::FieldContainer<double> basisGradsErr(numFields, numCubPointsErr, spaceDim);
     HGradBasis->getValues(basisValsErr, cubPointsErr, Intrepid::OPERATOR_VALUE);
     HGradBasis->getValues(basisGradsErr, cubPointsErr, Intrepid::OPERATOR_GRAD);

    // Loop over worksets
     for (int workset = 0; workset < numWorksets; workset++){

      // compute cell numbers where the workset starts and ends
       int worksetSize  = 0;
       int worksetBegin = (workset + 0)*desiredWorksetSize;
       int worksetEnd   = (workset + 1)*desiredWorksetSize;

      // when numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
        worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;
        worksetSize  = worksetEnd - worksetBegin;

      // loop over cells to fill arrays with coordinates and discrete solution coefficient
        Intrepid::FieldContainer<double> approxSolnCoef(worksetSize, numNodesPerElem);
        Intrepid::FieldContainer<double> cellWorksetErr(worksetSize, numNodesPerElem, spaceDim);

        int cellCounter = 0;
       // loop over workset cells and fill arrays with approximate solution
        for(int icell = worksetBegin; icell < worksetEnd; icell++){

           for (int inode = 0; inode < numNodesPerElem; inode++) {

             int rowIndex = elemToNode(icell,inode);
             approxSolnCoef(cellCounter, inode) = (*lhs).Values()[rowIndex];
             cellWorksetErr(cellCounter, inode, 0) = nodeCoords(elemToNode(icell,inode),0);
             cellWorksetErr(cellCounter, inode, 1) = nodeCoords(elemToNode(icell,inode),1);
           }
            cellCounter++;
        }

       // Containers for Jacobian
        Intrepid::FieldContainer<double> worksetJacobianE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
        Intrepid::FieldContainer<double> worksetJacobInvE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
        Intrepid::FieldContainer<double> worksetJacobDetE(worksetSize, numCubPointsErr);
        Intrepid::FieldContainer<double> worksetCubWeightsE(worksetSize, numCubPointsErr);

       // Containers for basis values and gradients in physical space
        Intrepid::FieldContainer<double> basisValsErrTrans(worksetSize, numFields, numCubPointsErr);
        Intrepid::FieldContainer<double> basisGradsErrTrans(worksetSize, numFields, numCubPointsErr, spaceDim);

       // compute cell Jacobians, their inverses and their determinants
        Intrepid::CellTools<double>::setJacobian(worksetJacobianE, cubPointsErr, cellWorksetErr, *cellType);
        Intrepid::CellTools<double>::setJacobianInv(worksetJacobInvE, worksetJacobianE );
        Intrepid::CellTools<double>::setJacobianDet(worksetJacobDetE, worksetJacobianE );

       // map cubature points to physical frame
        Intrepid::FieldContainer<double> worksetCubPointsErr(worksetSize, numCubPointsErr, cubDimErr);
        Intrepid::CellTools<double>::mapToPhysicalFrame(worksetCubPointsErr, cubPointsErr, cellWorksetErr, *cellType);

       // evaluate exact solution and gradient at cubature points
        Intrepid::FieldContainer<double> worksetExactSoln(worksetSize, numCubPointsErr);
        Intrepid::FieldContainer<double> worksetExactSolnGrad(worksetSize, numCubPointsErr, spaceDim);
        evaluateExactSolution(worksetExactSoln, worksetCubPointsErr);
        evaluateExactSolutionGrad(worksetExactSolnGrad, worksetCubPointsErr);

       // transform basis values to physical coordinates
        Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(basisValsErrTrans, basisValsErr);
        Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>(basisGradsErrTrans, worksetJacobInvE, basisGradsErr);

       // compute weighted measure
        Intrepid::FunctionSpaceTools::computeCellMeasure<double>(worksetCubWeightsE, worksetJacobDetE, cubWeightsErr);

       // evaluate the approximate solution and gradient at cubature points
        Intrepid::FieldContainer<double> worksetApproxSoln(worksetSize, numCubPointsErr);
        Intrepid::FieldContainer<double> worksetApproxSolnGrad(worksetSize, numCubPointsErr, spaceDim);
        Intrepid::FunctionSpaceTools::evaluate<double>(worksetApproxSoln, approxSolnCoef, basisValsErrTrans);
        Intrepid::FunctionSpaceTools::evaluate<double>(worksetApproxSolnGrad, approxSolnCoef, basisGradsErrTrans);

       // get difference between approximate and exact solutions
        Intrepid::FieldContainer<double> worksetDeltaSoln(worksetSize, numCubPointsErr);
        Intrepid::FieldContainer<double> worksetDeltaSolnGrad(worksetSize, numCubPointsErr, spaceDim);
        Intrepid::RealSpaceTools<double>::subtract(worksetDeltaSoln, worksetApproxSoln, worksetExactSoln);
        Intrepid::RealSpaceTools<double>::subtract(worksetDeltaSolnGrad, worksetApproxSolnGrad, worksetExactSolnGrad);

       // take absolute values
        Intrepid::RealSpaceTools<double>::absval(worksetDeltaSoln);
        Intrepid::RealSpaceTools<double>::absval(worksetDeltaSolnGrad);

       // apply cubature weights to differences in values and grads for use in integration
        Intrepid::FieldContainer<double> worksetDeltaSolnWeighted(worksetSize, numCubPointsErr);
        Intrepid::FieldContainer<double> worksetDeltaSolnGradWeighted(worksetSize, numCubPointsErr, spaceDim);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(worksetDeltaSolnWeighted,
                                                 worksetCubWeightsE, worksetDeltaSoln);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataData<double>(worksetDeltaSolnGradWeighted,
                                                worksetCubWeightsE, worksetDeltaSolnGrad);

       // integrate to get errors on each element
        Intrepid::FieldContainer<double> worksetL2err(worksetSize);
        Intrepid::FieldContainer<double> worksetH1err(worksetSize);
        Intrepid::FunctionSpaceTools::integrate<double>(worksetL2err, worksetDeltaSoln,
                                           worksetDeltaSolnWeighted, Intrepid::COMP_BLAS);
        Intrepid::FunctionSpaceTools::integrate<double>(worksetH1err, worksetDeltaSolnGrad,
                                           worksetDeltaSolnGradWeighted, Intrepid::COMP_BLAS);

        // loop over cells to get errors for total workset
        cellCounter = 0;
        for(int icell = worksetBegin; icell < worksetEnd; icell++){

           // loop over cubature points
            for(int nPt = 0; nPt < numCubPointsErr; nPt++){
                Linferr = std::max(Linferr, worksetDeltaSoln(cellCounter,nPt));
            }

           L2err += worksetL2err(cellCounter);
           H1err += worksetH1err(cellCounter);

           cellCounter++;

         } // end workset cell loop

     } // end loop over worksets



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
       std::cout << "\n";
       std::cout << "     L2 Error:   " << sqrt(L2errTot) <<"\n";
       std::cout << "     H1 Error:   " << sqrt(H1errTot) <<"\n";
       std::cout << "   LInf Error:   " << LinferrTot <<"\n\n";
     }


    if(showTiming) {std::cout << "Calculate error                             "
                  << Time.ElapsedTime() << " s \n"; Time.ResetStartTime();}

  } // end if problem == manufacturedSoln


 return 0;

}
/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/


/**********************************************************************************/
/********************** USER DEFINED FUNCTIONS FOR TEST ***************************/
/**********************************************************************************/

template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y) {

   //return  x + y;
   return  x*x*x - y*y;

}

template<typename Scalar>
void advectionVelocity(Scalar advVel[2], const Scalar& x, const Scalar& y, const std::string problem) {

     advVel[0] = 0;
     advVel[1] = 0;

     if (problem == "horizAdvection"){
        advVel[0] = 1.0;
        advVel[1] = 0.0;
     }

     if (problem == "skewAdvection" || problem == "manufacturedSoln"){
        advVel[0] = -sin(M_PI/6);
        advVel[1] = cos(M_PI/6);
     }

     if (problem == "doubleGlazing"){
        advVel[0] = 2.0*(2.0*y-1.0)*(1.0-(2.0*x - 1.0)*(2.0*x - 1.0));
        advVel[1] = -2.0*(2.0*x-1.0)*(1.0-(2.0*y - 1.0)*(2.0*y - 1.0));
     }

}

template<typename Scalar1, typename Scalar2>
const Scalar1 diffusivity(const Scalar1& x, const Scalar1& y, const Scalar2& epsilon, const bool variableEpsilon) {

  Scalar1 D = epsilon;
  if (variableEpsilon) {

     // Variable diffusion coefficient
       D = 0.5*epsilon*sin(M_PI*x)*cos(M_PI*y) + epsilon;

   }

   return D;

}


/**********************************************************************************/
/***************** AUXILIARY FUNCTIONS FROM USER DEFINED **************************/
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
template<typename Scalar1, typename Scalar2>
const Scalar1 sourceTerm(Scalar1& x, Scalar1& y, Scalar2& epsilon,
                         bool variableEpsilon, std::string problem)  {

  Scalar1 phi;
  Scalar1 D;
  Scalar1 u[2];
  Scalar1 grad_phi[2];
  Scalar1 flux[2];
  Scalar1 f = 0.;

  // Indicate the independent variables
   x.diff(0,2);
   y.diff(1,2);

  // Get exact solution and its gradient
   phi = exactSolution(x, y);
   exactSolutionGrad(grad_phi, x, y);

  // Get diffusivity
   D = diffusivity(x,y,epsilon,variableEpsilon);

  // Get advection velocity
   advectionVelocity(u, x, y, problem);

  // Compute total flux = J_n = epsilon grad phi - u phi
  for(int i = 0; i < 2; i++){
    flux[i] = D * grad_phi[i] - phi * u[i];
  }

  // Compute source term (right hand side): f = -div J_n
   f = -(flux[0].dx(0) + flux[1].dx(1));

  return f;
}

/**********************************************************************************/
/*************************** EVALUATION METHODS ***********************************/
/**********************************************************************************/

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
/****************************** CREATE MESH ***************************************/
/**********************************************************************************/

template<class ArrayOut1, class ArrayOut2, class Scalar>
void createMesh(ArrayOut1 &  elemToNode, ArrayOut2 & nodeCoords,
                ArrayOut1 &  bcLeftId, ArrayOut1 & bcRightId,
                ArrayOut1 &  bcTopId, ArrayOut1 & bcBotId,
                const std::string &  meshType, const Scalar & meshSize){

   int spaceDim = 2;

   // Create regular quadrilateral mesh on square domain [0,1]x[0,1]
   if (meshType == "quad")
   {
      int numElems = meshSize*meshSize;
      int numNodes = (meshSize+1)*(meshSize+1);
      int numNodesPerElem = 4;

      nodeCoords.resize(numNodes,spaceDim);
      elemToNode.resize(numElems,numNodesPerElem);

      double hx = 1.0/((double)meshSize);
      int nx = meshSize;

      int ielem = 0;
      for (int j=0; j<nx; j++){
         for (int i=0; i<nx; i++){

            int inode0 = (nx+1)*(j)   + i;
            int inode1 = (nx+1)*(j)   + i+1;
            int inode2 = (nx+1)*(j+1) + i+1;
            int inode3 = (nx+1)*(j+1) + i;

            nodeCoords(inode0,0) = (double)i*hx;
            nodeCoords(inode0,1) = (double)j*hx;
            nodeCoords(inode1,0) = (double)(i+1)*hx;
            nodeCoords(inode1,1) = (double)j*hx;
            nodeCoords(inode2,0) = (double)(i+1)*hx;
            nodeCoords(inode2,1) = (double)(j+1)*hx;
            nodeCoords(inode3,0) = (double)i*hx;
            nodeCoords(inode3,1) = (double)(j+1)*hx;

            elemToNode(ielem,0) = inode0;
            elemToNode(ielem,1) = inode1;
            elemToNode(ielem,2) = inode2;
            elemToNode(ielem,3) = inode3;

            ielem++;
         }
      }

   }
   // Create regular triangular mesh on square domain [0,1]x[0,1]
   else if (meshType == "tri")
   {
      int numElems = 2*meshSize*meshSize;
      int numNodes = (meshSize+1)*(meshSize+1);
      int numNodesPerElem = 3;

      nodeCoords.resize(numNodes,spaceDim);
      elemToNode.resize(numElems,numNodesPerElem);

      double hx = 1.0/((double)meshSize);
      int nx = meshSize;

      // loop over quads and create two triangles for each quad
      //  _________
      // |        /|
      // |      /  |
      // |    /    |
      // |  /      |
      // |/________|
      //
      int ielem = 0;
      for (int j=0; j<nx; j++){
         for (int i=0; i<nx; i++){

            int inode0 = (nx+1)*(j)   + i;
            int inode1 = (nx+1)*(j)   + i+1;
            int inode2 = (nx+1)*(j+1) + i;

            int inode3 = (nx+1)*(j+1) + i;
            int inode4 = (nx+1)*(j)   + i+1;
            int inode5 = (nx+1)*(j+1) + i+1;

            nodeCoords(inode0,0) = (double)i*hx;
            nodeCoords(inode0,1) = (double)j*hx;
            nodeCoords(inode1,0) = (double)(i+1)*hx;
            nodeCoords(inode1,1) = (double)j*hx;
            nodeCoords(inode2,0) = (double)i*hx;
            nodeCoords(inode2,1) = (double)(j+1)*hx;

            nodeCoords(inode3,0) = (double)i*hx;
            nodeCoords(inode3,1) = (double)(j+1)*hx;
            nodeCoords(inode4,0) = (double)(i+1)*hx;
            nodeCoords(inode4,1) = (double)j*hx;
            nodeCoords(inode5,0) = (double)(i+1)*hx;
            nodeCoords(inode5,1) = (double)(j+1)*hx;

            elemToNode(ielem,0) = (nx+1)*j + i;
            elemToNode(ielem,1) = (nx+1)*j + i+1;
            elemToNode(ielem,2) = (nx+1)*(j+1) + i;

            elemToNode(ielem+1,0) = (nx+1)*(j+1) + i;
            elemToNode(ielem+1,1) = (nx+1)*j + i+1;
            elemToNode(ielem+1,2) = (nx+1)*(j+1) + i+1;

            ielem+=2;
         }
      }

   }


   // Fill boundary arrays
    bcLeftId.resize(meshSize);
    bcRightId.resize(meshSize);
    bcTopId.resize(meshSize);
    bcBotId.resize(meshSize);

   // For this square grid each boundary has meshSize nodes
   //  (the left corner of each boundary belongs to that boundary)
    int nx = meshSize;
    for (int ibcnode = 0; ibcnode < nx; ibcnode++){

       bcBotId(ibcnode) = ibcnode;
       bcTopId(ibcnode) = (nx+1)*nx + ibcnode + 1;
       bcLeftId(ibcnode) = (nx+1)*(ibcnode+1);
       bcRightId(ibcnode) = (nx+1)*(ibcnode) + nx;

    }

}

/**********************************************************************************/
/************************************** END ***************************************/
/**********************************************************************************/
