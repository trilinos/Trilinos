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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   example_Poisson.cpp
    \brief  Example solution of a Poisson equation on a hexahedral mesh using
            nodal (Hgrad) elements.

            This example uses Pamgen to generate a hexahedral mesh, Intrepid to
            build the stiffness matrix and right-hand side, and ML to solve.

    \verbatim

     Poisson system:
 
            div grad u = f in Omega
                     u = 0 on Gamma

     Corresponding discrete linear system for nodal coefficients(x):

                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson, D. Hensinger, C. Siefert.

    \remark Usage:
    \code   ./TrilinosCouplings_examples_scaling_example_Poisson.exe \endcode

    \remark Example requires Pamgen formatted mesh input file named PoissonMesh.in.
*/

//#define DUMP_DATA

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
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

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "im_exodusII_l.h"
#include "im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"


int TestMultiLevelPreconditionerLaplace(char ProblemType[],
				 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol);


using namespace std;
using namespace Intrepid;


// Functions to evaluate exact solution and derivatives
double evalu(double & x, double & y, double & z);
int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3);
double evalDivGradu(double & x, double & y, double & z);

int main(int argc, char *argv[]) {
  int error = 0;
  int numProcs=1;
  int rank=0;
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  rank=mpiSession.getRank();
  numProcs=mpiSession.getNProc();
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
    << "|          Example: Solve Poisson Equation on Hexahedral Mesh                 |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Pamgen's website:   http://trilinos.sandia.gov/packages/pamgen             |\n" \
    << "|  ML's website:       http://trilinos.sandia.gov/packages/ml                 |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  }

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

#ifdef HAVE_MPI
  if (MyPID == 0) {
    std::cout << "PARALLEL executable \n"; 
  }
#else
  if (MyPID == 0) {
    std::cout << "SERIAL executable \n";
  }
#endif


// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int spaceDim = hex_8.getDimension();
    int dim = 3;

// *********************************** GENERATE MESH ************************************

  if (MyPID == 0) {
    std::cout << "Generating mesh ... \n\n";
  }

  // Read in Pamgen mesh file
    string  meshInput;
    string  tmp;

    ifstream finput;
    finput.open("PoissonMesh.in");
    if (finput.is_open()){
      while(getline(finput,tmp)){
        meshInput += tmp;
        meshInput += "\n";
      }
    }
    else {
       std::cout << "Cannot open mesh file: PoissonMesh.in" <<"\n";
       return 0;
    }
    finput.close();

    if (MyPID == 0) {
      std::cout << meshInput <<"\n";
    }

   // Generate mesh with Pamgen
    long long maxInt = 9223372036854775807LL;
    Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);

    string msg("Poisson: ");
    if(MyPID == 0) {cout << msg << "Pamgen Setup     = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}
    
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
       std::cout << "    Number of Global Nodes: " << numNodesGlobal << " \n\n";
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
    FieldContainer<int> elemToNode(numElems,numNodesPerElem);
    for(long long b = 0; b < numElemBlk; b++){
      for(long long el = 0; el < elements[b]; el++){
	for (int j=0; j<numNodesPerElem; j++) {
	  elemToNode(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
	}
	telct ++;
      }
    }
 

   // Read node coordinates and place in field container
    FieldContainer<double> nodeCoord(numNodes,dim);
    double * nodeCoordx = new double [numNodes];
    double * nodeCoordy = new double [numNodes];
    double * nodeCoordz = new double [numNodes];
    im_ex_get_coord_l(id,nodeCoordx,nodeCoordy,nodeCoordz);
    for (int i=0; i<numNodes; i++) {          
      nodeCoord(i,0)=nodeCoordx[i];
      nodeCoord(i,1)=nodeCoordy[i];
      nodeCoord(i,2)=nodeCoordz[i];
    }
    delete [] nodeCoordx;
    delete [] nodeCoordy;


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

    if(!Comm.MyPID()) {cout << msg << "Mesh Queries     = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}

    //Calculate global node ids
    long long * globalNodeIds = new long long[numNodes];
    bool * nodeIsOwned = new bool[numNodes];

    calc_global_node_ids(globalNodeIds,
			 nodeIsOwned,
			 numNodes,
			 num_node_comm_maps,
			 node_cmap_node_cnts,
			 node_comm_proc_ids,
			 comm_node_ids,
			 rank);    

    if(MyPID==0) {cout << msg << "Global Node Nums = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}

#ifdef DUMP_DATA
    // Print coords
    std::stringstream fname;
      fname << "coords";
      fname << MyPID << ".dat";
    FILE *f=fopen(fname.str().c_str(),"w");
    for (int i=0; i<numNodes; i++) {
      if (nodeIsOwned[i]) {
       fprintf(f,"%22.16e %22.16e %22.16e\n",nodeCoord(i,0),nodeCoord(i,1),nodeCoord(i,2));
      }
    }
    fclose(f);

  // Output element to node connectivity
    std::stringstream efname;
      efname << "elem2node";
      efname << MyPID << ".dat";
    ofstream el2nout(efname.str().c_str());
    for (int i=0; i<numElems; i++) {
      for (int m=0; m<numNodesPerElem; m++) {
        el2nout << globalNodeIds[elemToNode(i,m)] << "  ";
      }
      el2nout << "\n";
    }
    el2nout.close();
#endif
  
   // Container indicating whether a node is on the boundary (1-yes 0-no)
    FieldContainer<int> nodeOnBoundary(numNodes);

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
             
             int sideNode0 = hex_8.getNodeMap(2,sideSetSideList[j]-1,0);
             int sideNode1 = hex_8.getNodeMap(2,sideSetSideList[j]-1,1);
             int sideNode2 = hex_8.getNodeMap(2,sideSetSideList[j]-1,2);
             int sideNode3 = hex_8.getNodeMap(2,sideSetSideList[j]-1,3);
             
             nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode0))=1;
             nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode1))=1;
             nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode2))=1;
             nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode3))=1;
          }
          delete [] sideSetElemList;
          delete [] sideSetSideList;
       }
    }
    delete [] sideSetIds;

   if(MyPID ==0) {cout << msg << "Boundary Conds   = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}

// ************************************ CUBATURE **************************************

  if (MyPID == 0) {
    std::cout << "Getting cubature ... \n\n";
  }

   // Get numerical integration points and weights
    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);


// ************************************** BASIS ***************************************

  if (MyPID == 0) {
     std::cout << "Getting basis ... \n\n";
  }

   // Define basis 
     Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;
     int numFieldsG = hexHGradBasis.getCardinality();
     FieldContainer<double> hexGVals(numFieldsG, numCubPoints); 
     FieldContainer<double> hexGrads(numFieldsG, numCubPoints, spaceDim); 

  // Evaluate basis values and gradients at cubature points
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);
     hexHGradBasis.getValues(hexGrads, cubPoints, OPERATOR_GRAD);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL STIFFNESS MATRIX *************

  if (MyPID == 0) {
    std::cout << "Building stiffness matrix and right hand side ... \n\n";
  }

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    int numCells = 1; 

   // Container for nodes
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HGRAD stiffness matrix
    FieldContainer<double> localStiffMatrix(numCells, numFieldsG, numFieldsG);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> hexGradsTransformed(numCells, numFieldsG, numCubPoints, spaceDim);
    FieldContainer<double> hexGradsTransformedWeighted(numCells, numFieldsG, numCubPoints, spaceDim);
   // Containers for right hand side vectors
    FieldContainer<double> rhsData(numCells, numCubPoints);
    FieldContainer<double> localRHS(numCells, numFieldsG);
    FieldContainer<double> hexGValsTransformed(numCells, numFieldsG, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    // Count owned nodes
    int ownedNodes=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]) ownedNodes++;

    // Build a list of the OWNED global ids...
    // NTS: will need to switch back to long long
    int *ownedGIDs=new int[ownedNodes];    
    int oidx=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]){
        ownedGIDs[oidx]=(int)globalNodeIds[i];
        oidx++;
      }
    // Generate epetra map    
    Epetra_Map globalMapG(-1,ownedNodes,ownedGIDs,0,Comm);


    // Global arrays in Epetra format
    Epetra_FECrsMatrix StiffMatrix(Copy, globalMapG, numFieldsG);
    Epetra_FEVector rhs(globalMapG);

    
 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );

// ************************** Compute element HGrad stiffness matrices *******************************
  
     // transform to physical coordinates 
      fst::HGRADtransformGRAD<double>(hexGradsTransformed, hexJacobInv, hexGrads);
      
     // compute weighted measure
      fst::computeMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGradsTransformedWeighted,
                                   weightedMeasure, hexGradsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(localStiffMatrix,
                             hexGradsTransformed, hexGradsTransformedWeighted, COMP_CPP);

      // assemble into global matrix
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = globalNodeIds[elemToNode(k,row)];
            int colIndex = globalNodeIds[elemToNode(k,col)];
            double val = localStiffMatrix(0,row,col);
            StiffMatrix.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side functions at physical points
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);

          rhsData(0,nPt) = evalDivGradu(x, y, z);
       }

     // transform basis values to physical coordinates 
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);
      
     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasure, hexGValsTransformed);

     // integrate rhs term
      fst::integrate<double>(localRHS, rhsData, hexGValsTransformedWeighted,
                             COMP_CPP);

    // assemble into global vector
     for (int row = 0; row < numFieldsG; row++){
           int rowIndex = globalNodeIds[elemToNode(k,row)];
           double val = -localRHS(0,row);
           rhs.SumIntoGlobalValues(1, &rowIndex, &val);
      }
 
     
 } // *** end element loop ***


  // Assemble over multiple processors
   StiffMatrix.GlobalAssemble(); StiffMatrix.FillComplete();
   rhs.GlobalAssemble();

   if(MyPID == 0) {cout << msg << "Matrix Assembly  = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}

  // Adjust matrix and rhs due to Dirichlet boundary conditions
    int numBCNodes=0;
    for (int i=0; i<numNodes; i++){
        if (nodeOnBoundary(i) && nodeIsOwned[i]){
            numBCNodes++;
        }
    }
    int * BCNodes = new int [numBCNodes];
    int indbc=0;
    int indOwned=0;
    for (int i=0; i<numNodes; i++){
       if (nodeIsOwned[i]){
          if (nodeOnBoundary(i)){
             BCNodes[indbc]=indOwned;
             indbc++;
             rhs[0][indOwned]=0;
          }
          indOwned++;
       }
    }
   // ML routine that zeroes out Dirichlet rows and columns and adds 1 to diagonal 
    ML_Epetra::Apply_OAZToMatrix(BCNodes, numBCNodes, StiffMatrix);

    delete [] BCNodes;


#ifdef DUMP_DATA   
  // Dump matrices to disk
    EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
    EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhs,0,0,false);
#endif


   // Run the solver
   Teuchos::ParameterList MLList;
   ML_Epetra::SetDefaults("SA",MLList);
   Epetra_FEVector xexact(globalMapG);
   Epetra_FEVector uh(globalMapG);
   double TotalErrorResidual=0, TotalErrorExactSol=0;

   // Get exact solution at nodes
    for (int i = 0; i<numNodes; i++) {
       if (nodeIsOwned[i]){
          double x = nodeCoord(i,0);
          double y = nodeCoord(i,1);
          double z = nodeCoord(i,2);
          double exactu = evalu(x, y, z);
          int rowindex=globalNodeIds[i];
          xexact.SumIntoGlobalValues(1, &rowindex, &exactu);
       }
    }
    xexact.GlobalAssemble();
       
   char probType[10] = "laplace";
   
    TestMultiLevelPreconditionerLaplace(probType,MLList,StiffMatrix,xexact,rhs,uh,
                                       TotalErrorResidual, TotalErrorExactSol);

   // ********  Calculate Error in Solution *************** 
     double L2err = 0.0;
     double L2errTot = 0.0;
     double H1err = 0.0;
     double H1errTot = 0.0;
     double Linferr = 0.0;
     double LinferrTot = 0.0;

#ifdef HAVE_MPI
   // Import solution onto current processor
     Epetra_Map  solnMap(numNodesGlobal, numNodesGlobal, 0, Comm);
     Epetra_Import  solnImporter(solnMap, globalMapG);
     Epetra_Vector  uCoeff(solnMap);
     uCoeff.Import(uh, solnImporter, Insert);
#endif

   // Get cubature points and weights for error calc (may be different from previous)
     DefaultCubatureFactory<double>  cubFactoryErr;                                   
     int cubDegErr = 3;
     Teuchos::RCP<Cubature<double> > hexCubErr = cubFactoryErr.create(hex_8, cubDegErr); 
     int cubDimErr       = hexCubErr->getDimension();
     int numCubPointsErr = hexCubErr->getNumPoints();
     FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
     FieldContainer<double> cubWeightsErr(numCubPointsErr);
     hexCubErr->getCubature(cubPointsErr, cubWeightsErr);

   // Containers for Jacobian
     FieldContainer<double> hexJacobianE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobInvE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobDetE(numCells, numCubPointsErr);
     FieldContainer<double> weightedMeasureE(numCells, numCubPointsErr);

  // Evaluate basis values and gradients at cubature points
     FieldContainer<double> uhGVals(numFieldsG, numCubPointsErr); 
     FieldContainer<double> uhGValsTrans(numCells,numFieldsG, numCubPointsErr); 
     FieldContainer<double> uhGrads(numFieldsG, numCubPointsErr, spaceDim); 
     FieldContainer<double> uhGradsTrans(numCells, numFieldsG, numCubPointsErr, spaceDim); 
     hexHGradBasis.getValues(uhGVals, cubPointsErr, OPERATOR_VALUE);
     hexHGradBasis.getValues(uhGrads, cubPointsErr, OPERATOR_GRAD);


   // Loop over elements
    for (int k=0; k<numElems; k++){

      double L2errElem = 0.0;
      double H1errElem = 0.0;
      double uExact; 
      double graduExact1, graduExact2, graduExact3;

     // physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

    // compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobianE, cubPointsErr, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInvE, hexJacobianE );
       CellTools::setJacobianDet(hexJacobDetE, hexJacobianE );

      // transform integration points to physical points
       FieldContainer<double> physCubPoints(numCells,numCubPointsErr, cubDimErr);
       CellTools::mapToPhysicalFrame(physCubPoints, cubPointsErr, hexNodes, hex_8);

      // transform basis values to physical coordinates 
       fst::HGRADtransformVALUE<double>(uhGValsTrans, uhGVals);
       fst::HGRADtransformGRAD<double>(uhGradsTrans, hexJacobInvE, uhGrads);

      // compute weighted measure
       fst::computeMeasure<double>(weightedMeasureE, hexJacobDetE, cubWeightsErr);

      // loop over cubature points
       for (int nPt = 0; nPt < numCubPointsErr; nPt++){

         // get exact solution and gradients
          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          uExact = evalu(x, y, z);
          evalGradu(x, y, z, graduExact1, graduExact2, graduExact3);

         // calculate approximate solution and gradients
          double uApprox = 0.0;
          double graduApprox1 = 0.0;
          double graduApprox2= 0.0;
          double graduApprox3 = 0.0;
          for (int i = 0; i < numFieldsG; i++){
             int rowIndex = globalNodeIds[elemToNode(k,i)];
#ifdef HAVE_MPI
             double uh1 = uCoeff.Values()[rowIndex];
#else
             double uh1 = uh.Values()[rowIndex];
#endif
             uApprox += uh1*uhGValsTrans(0,i,nPt); 
             graduApprox1 += uh1*uhGradsTrans(0,i,nPt,0); 
             graduApprox2 += uh1*uhGradsTrans(0,i,nPt,1); 
             graduApprox3 += uh1*uhGradsTrans(0,i,nPt,2); 
          }

         // evaluate the error at cubature points
          Linferr = max(Linferr, abs(uExact - uApprox));

          L2errElem+=(uExact - uApprox)*(uExact - uApprox)*weightedMeasureE(0,nPt);
          H1errElem+=((graduExact1 - graduApprox1)*(graduExact1 - graduApprox1))
                     *weightedMeasureE(0,nPt);
          H1errElem+=((graduExact2 - graduApprox2)*(graduExact2 - graduApprox2))
                     *weightedMeasureE(0,nPt);
          H1errElem+=((graduExact3 - graduApprox3)*(graduExact3 - graduApprox3))
                     *weightedMeasureE(0,nPt);
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
   for(long long b = 0; b < numElemBlk; b++){     
     delete [] elmt_node_linkage[b];
     delete [] element_types[b];
   }
   delete [] block_ids;
   delete [] nodes_per_element;
   delete [] element_attributes;
   delete [] element_types;
   delete [] elmt_node_linkage;
   delete [] ownedGIDs;
   delete [] elements;
   delete [] globalNodeIds;
   delete [] nodeIsOwned;
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

   // delete mesh
   Delete_Pamgen_Mesh();
   
   // reset format state of std::cout
//   std::cout.copyfmt(oldFormatState);
   
#ifdef HAVE_MPI
   MPI_Finalize();
#endif
 
   exit(0);

}


// Calculates value of exact solution u
 double evalu(double & x, double & y, double & z)
 {
   // u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)
  // double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);

  // or

   // u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x+y+z)
   double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

   return exactu;
 }

// Calculates gradient of exact solution u
 int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3)
 {
   //  for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)
   /*
        gradu1 = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
        gradu2 = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
        gradu3 = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
   */
  

  // or

   // for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x+y+z)
     gradu1 = (M_PI*cos(M_PI*x)+sin(M_PI*x))
                  *sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);
     gradu2 = (M_PI*cos(M_PI*y)+sin(M_PI*y))
                  *sin(M_PI*x)*sin(M_PI*z)*exp(x+y+z);
     gradu3 = (M_PI*cos(M_PI*z)+sin(M_PI*z))
                  *sin(M_PI*x)*sin(M_PI*y)*exp(x+y+z);
  
   return 0;
 }

// Calculates Laplacian of exact solution u
 double evalDivGradu(double & x, double & y, double & z)
 {
   //  for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)
   //double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);

  // or

   // for u(x,y,z)=sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x+y+z)
   double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*exp(x+y+z)
                    + 3.0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);
  
  
   
   return divGradu;
 }

#ifdef HAVE_MPI

#endif


// Test ML
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
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
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 32);

  solver.Iterate(200, 1e-10);
  
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


