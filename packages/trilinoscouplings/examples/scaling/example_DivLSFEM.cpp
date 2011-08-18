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

/** \file   example_DivLSFEM.cpp
    \brief  Example solution of a div-curl system on a hexahedral mesh using
            div-conforming (face) elements.

          This example uses the following Trilinos packages:
    \li        Pamgen to generate a Hexahedral mesh.
    \li        Intrepid to build the discretization matrices and right-hand side.
    \li        Epetra to handle the global matrix and vector.
    \li        ML to solve the linear system.


         For more details on the formulation see

    \li     P. Bochev, K. Peterson, C. Siefert, "Analysis and Computation
            of Compatible Least-Squares Methods for Div-Curl Equations",
            SIAM J. Numerical Analysis, vol 49, pp 159-191, 2011.

    \verbatim

            Div-Curl System:

                       curl u = g  in Omega
                        div u = h  in Omega
                          u.n = 0  on Gamma

            Corresponding discrete linear system for face element coeficients (x):

                      (Kd + Md*Dc*McInv*Dc'*Md)x = b

                      Kd    - Hdiv stiffness matrix
                      Md    - Hdiv mass matrix
                      Dc    - Edge to Face incidence matrix
                      McInv - Hgrad mass matrix inverse
                      b     - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson, D. Hensinger, C. Siefert.

     \remark Usage
     \verbatim

     ./TrilinosCouplings_examples_scaling_Example_DivLSFEM.exe  inputfile.xml

        inputfile.xml (optional)  -  xml input file containing Pamgen mesh description
                                     and material parameters for each Pamgen block,
                                     if not present code attempts to read DivLSFEMin.xml.


     \endverbatim
 **/


// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HDIV_HEX_I1_FEM.hpp"
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
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "im_exodusII_l.h"
#include "im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_epetra_utils.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_FaceMatrixFreePreconditioner.h"
#include "ml_RefMaxwell.h"

#define ABS(x) ((x)>0?(x):-(x))

/*** Uncomment if you would like output data for plotting ***/
//#define DUMP_DATA

using namespace std;
using namespace Intrepid;


/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;


struct fecomp{
  bool operator () ( topo_entity* x,  topo_entity*  y )const
  {
    if(x->sorted_local_node_ids < y->sorted_local_node_ids)return true;    
    return false;
  }
};

/**********************************************************************************/
/***************** FUNCTION DECLARATIONS FOR ML PRECONDITIONER ********************/
/**********************************************************************************/

/** \brief Multiplies abs(A)x = y, where all non-zero entries of A are replaced with their absolute values value  

    \param  A          [in]    matrix
    \param  x          [in]    vector
    \param  y          [out]   vector 
 */
int Multiply_Abs(const Epetra_CrsMatrix &A,const Epetra_Vector &x,Epetra_Vector &y);


/** \brief  ML Preconditioner

    \param  ProblemType        [in]    problem type
    \param  MLList             [in]    ML parameter list
    \param  GradDiv            [in]    H(curl) stiffness matrix
    \param  D0clean            [in]    Edge to node incidence matrix
    \param  D1clean            [in]    Face to edge incidence matrix
    \param  FaceNode           [in]    Face to node incidence matrix
    \param  M1inv              [in]    H(curl) mass matrix inverse
    \param  M2                 [in]    H(div) mass matrix
    \param  xh                 [out]   solution vector
    \param  b                  [in]    right-hand-side vector
    \param  TotalErrorResidual [out]   error residual
    \param  TotalErrorExactSol [out]   error in xh

 */
void TestMultiLevelPreconditioner_DivLSFEM(char ProblemType[],
                                           Teuchos::ParameterList   & MLList,
                                           Epetra_CrsMatrix   & GradDiv,
                                           Epetra_CrsMatrix   & D0clean,
                                           Epetra_CrsMatrix   & D1clean,
                                           Epetra_CrsMatrix   & FaceNode,
                                           Epetra_CrsMatrix   & M1,
                                           Epetra_CrsMatrix   & M1inv,
                                           Epetra_CrsMatrix   & M2,
                                           Epetra_MultiVector & xh,
                                           Epetra_MultiVector & b,
                                           double & TotalErrorResidual,
                                           double & TotalErrorExactSol);

/**********************************************************************************/
/******** FUNCTION DECLARATIONS FOR EXACT SOLUTION AND SOURCE TERMS ***************/
/**********************************************************************************/

/** \brief  Exact solution evaluation.

    \param  uExact0            [out]   first component of exact solution at (x,y,z)
    \param  uExact1            [out]   second component of exact solution at (x,y,z)
    \param  uExact2            [out]   third component of exact solution at (x,y,z)
    \param  x                  [in]    x coordinate
    \param  y                  [in]    y coordinate
    \param  z                  [in]    z coordinate
 */
int evalu(double & uExact0, 
          double & uExact1, 
          double & uExact2, 
          double & x, 
          double & y, 
          double & z);

/** \brief  Divergence of exact solution.

    \param  x                  [in]    x coordinate
    \param  y                  [in]    y coordinate
    \param  z                  [in]    z coordinate

    \return Value of the divergence of exact solution at (x,y,z)
 */
double evalDivu(double & x, 
                double & y, 
                double & z);

/** \brief  Curl of exact solution.

    \param  curlu0            [out]   first component of curl of exact solution
    \param  curlu1            [out]   second component of curl of exact solution
    \param  curlu2            [out]   third component of curl of exact solution
    \param  x                  [in]    x coordinate
    \param  y                  [in]    y coordinate
    \param  z                  [in]    z coordinate
    \param  mu                 [in]    material parameter
 */
int evalCurlu(double & curlu0, 
              double & curlu1, 
              double & curlu2, 
              double & x, 
              double & y, 
              double & z, 
              double & mu);

/** \brief  Curl of curl of exact solution.

    \param  curlCurlu0         [out]   first component of curl curl of exact solution
    \param  curlCurlu1         [out]   second component of curl curl of exact solution
    \param  curlCurlu2         [out]   third component of curl curl of exact solution
    \param  x                  [in]    x coordinate
    \param  y                  [in]    y coordinate
    \param  z                  [in]    z coordinate
    \param  mu                 [in]    material parameter
 */
int evalCurlCurlu(double & curlCurlu0, 
                  double & curlCurlu1, 
                  double & curlCurlu2, 
                  double & x, 
                  double & y, 
                  double & z, 
                  double & mu);
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/



/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/
int main(int argc, char *argv[]) {

  int error = 0;
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  int rank=mpiSession.getRank();
  int numProcs=mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  int MyPID = Comm.MyPID();
#else
  int rank=0;
  int numProcs=1;
  int MyPID = 0;
  Epetra_SerialComm Comm;
#endif
  Epetra_Time Time(Comm);

   //Check number of arguments
  if (argc > 2) {
      std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
      std::cout <<"Usage:\n\n";
      std::cout <<"  ./TrilinosCouplings_examples_scaling_Example_DivLSFEM.exe [inputfile.xml] \n\n";
      std::cout <<"   inputfile.xml(optional) - xml file with description of Pamgen mesh \n";
      std::cout <<"                             and material parameters for each block \nn";
      exit(1);
   }
  
 if (MyPID == 0) {
  std::cout \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|          Example: Div-Curl System on Hexahedral Mesh                        |\n" \
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

  std::set < topo_entity * , fecomp > edge_set;
  std::set < topo_entity * , fecomp > face_set;

  std::vector < topo_entity * > edge_vector;
  std::vector < topo_entity * > face_vector;

  std::vector < int > edge_comm_procs;

  int dim = 3;

/**********************************************************************************/
/********************************** GET XML INPUTS ********************************/
/**********************************************************************************/

  // Command line for xml file, otherwise use default
    std::string   xmlInFileName;
    if(argc>=2) xmlInFileName=string(argv[1]);
    else xmlInFileName="DivLSFEMin.xml";

  // Read xml file into parameter list
    Teuchos::ParameterList inputList;

   if(xmlInFileName.length()) {
     if (MyPID == 0) {
      std::cout << "\nReading parameter list from the XML file \""<<xmlInFileName<<"\" ...\n\n";
     }
      Teuchos::updateParametersFromXmlFile(xmlInFileName,&inputList);
     if (MyPID == 0) {
      inputList.print(std::cout,2,true,true);
      std::cout << "\n";
     }
    }
    else
    {
      std::cout << "Cannot read input file: " << xmlInFileName << "\n";
      return 0;
    }

  // Get pamgen mesh definition
    std::string meshInput = Teuchos::getParameter<std::string>(inputList,"meshInput");
 

/**********************************************************************************/
/***************************** GET CELL TOPOLOGY **********************************/
/**********************************************************************************/

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology cellType(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = cellType.getNodeCount();
    int numEdgesPerElem = cellType.getEdgeCount();
    int numFacesPerElem = cellType.getSideCount();
    int numNodesPerEdge = 2;
    int numNodesPerFace = 4;
    int numEdgesPerFace = 4;
    int spaceDim = cellType.getDimension();

   // Build reference element edge to node map
    FieldContainer<int> refEdgeToNode(numEdgesPerElem,numNodesPerEdge);
    for (int i=0; i<numEdgesPerElem; i++){
        refEdgeToNode(i,0)=cellType.getNodeMap(1, i, 0);
        refEdgeToNode(i,1)=cellType.getNodeMap(1, i, 1);
    }

   // Build reference element face to node map
    FieldContainer<int> refFaceToNode(numFacesPerElem,numNodesPerFace);
    for (int i=0; i<numFacesPerElem; i++){
        refFaceToNode(i,0)=cellType.getNodeMap(2, i, 0);
        refFaceToNode(i,1)=cellType.getNodeMap(2, i, 1);
        refFaceToNode(i,2)=cellType.getNodeMap(2, i, 2);
        refFaceToNode(i,3)=cellType.getNodeMap(2, i, 3);
    }

   // Build reference element face to edge map (Hardcoded for now)
    FieldContainer<int> refFaceToEdge(numFacesPerElem,numEdgesPerFace);
        refFaceToEdge(0,0)=0; refFaceToEdge(0,1)=9;
        refFaceToEdge(0,2)=4; refFaceToEdge(0,3)=8;
        refFaceToEdge(1,0)=1; refFaceToEdge(1,1)=10;
        refFaceToEdge(1,2)=5; refFaceToEdge(1,3)=9;
        refFaceToEdge(2,0)=2; refFaceToEdge(2,1)=11;
        refFaceToEdge(2,2)=6; refFaceToEdge(2,3)=10;
        refFaceToEdge(3,0)=3; refFaceToEdge(3,1)=8;
        refFaceToEdge(3,2)=7; refFaceToEdge(3,3)=11;
        refFaceToEdge(4,0)=0; refFaceToEdge(4,1)=1;
        refFaceToEdge(4,2)=2; refFaceToEdge(4,3)=3;
        refFaceToEdge(5,0)=4; refFaceToEdge(5,1)=5;
        refFaceToEdge(5,2)=6; refFaceToEdge(5,3)=7;

/**********************************************************************************/
/******************************* GENERATE MESH ************************************/
/**********************************************************************************/

   // Generate mesh with Pamgen

    long long maxInt = 9223372036854775807LL;
    Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
    
   // Get local mesh size info
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
  // Get global mesh size info
    long long numNodesGlobal;
    long long numElemsGlobal;
    long long numElemBlkGlobal;
    long long numNodeSetsGlobal;
    long long numSideSetsGlobal;

    im_ne_get_init_global_l(id, &numNodesGlobal, &numElemsGlobal,
                         &numElemBlkGlobal, &numNodeSetsGlobal,
                         &numSideSetsGlobal);


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

   // Get mu value for each block of elements from parameter list
    double  *mu = new double [numElemBlk];
    for(int b = 0; b < numElemBlk; b++){
       stringstream muBlock;
       muBlock.clear();
       muBlock << "mu" << b;
       mu[b] = inputList.get(muBlock.str(),1.0);
    }

  // Get node-element connectivity and set element mu value
    int telct = 0;
    FieldContainer<int> elemToNode(numElems,numNodesPerElem);
    FieldContainer<double> muVal(numElems);
    for(long long b = 0; b < numElemBlk; b++){
      for(long long el = 0; el < elements[b]; el++){
	for (int j=0; j<numNodesPerElem; j++) {
	  elemToNode(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
	}
        muVal(telct) = mu[b];     
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
    }


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
  // Count owned nodes
    int ownedNodes=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]) ownedNodes++;

   // Build a list of the OWNED global ids...
    int *ownedGIDs=new int [ownedNodes];
    int oidx=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]){
        ownedGIDs[oidx]=(int)globalNodeIds[i];
        oidx++;
      }


    FieldContainer<int> elemToEdge(numElems,numEdgesPerElem);
    FieldContainer<int> elemToFace(numElems,numFacesPerElem);

   // calculate edge and face ids
    int elct = 0;
    for(long long b = 0; b < numElemBlk; b++){
      if(nodes_per_element[b] == 4){
      }
      else if (nodes_per_element[b] == 8){
	//loop over all elements and push their edges onto a set if they are not there already
	for(long long el = 0; el < elements[b]; el++){
	  std::set< topo_entity *, fecomp > ::iterator fit;
	  for (int i=0; i < numEdgesPerElem; i++){
	    topo_entity * teof = new topo_entity;
	    for(int j = 0; j < numNodesPerEdge;j++){
	      teof->add_node(elmt_node_linkage[b][el*numNodesPerElem + refEdgeToNode(i,j)],globalNodeIds);
	    }
	    teof->sort();
	    fit = edge_set.find(teof);
	    if(fit == edge_set.end()){
	      teof->local_id = edge_vector.size();
	      edge_set.insert(teof);
	      elemToEdge(elct,i)= edge_vector.size();
	      edge_vector.push_back(teof);
	    }
	    else{
	      elemToEdge(elct,i) = (*fit)->local_id;
	      delete teof;
	    }
	  }
	  for (int i=0; i < numFacesPerElem; i++){
	    topo_entity * teof = new topo_entity;
	    for(int j = 0; j < numNodesPerFace;j++){
	      teof->add_node(elmt_node_linkage[b][el*numNodesPerElem + refFaceToNode(i,j)],globalNodeIds);
	    }
	    teof->sort();
	    fit = face_set.find(teof);
	    if(fit == face_set.end()){
	      teof->local_id = face_vector.size();
	      face_set.insert(teof);
	      elemToFace(elct,i)= face_vector.size();
	      face_vector.push_back(teof);
	    }
	    else{
	      elemToFace(elct,i) = (*fit)->local_id;
	      delete teof;
	    }
	  }
	  elct ++;
	}
      }
    }
    
   // Edge to Node connectivity
    FieldContainer<int> edgeToNode(edge_vector.size(), numNodesPerEdge);
    for(unsigned ect = 0; ect != edge_vector.size(); ect++){
      std::list<long long>::iterator elit;
      int nct = 0;
      for(elit  = edge_vector[ect]->local_node_ids.begin();
	  elit != edge_vector[ect]->local_node_ids.end();
	  elit ++){
	edgeToNode(ect,nct) = *elit-1;
	nct++;
      }
    }

   // Face to Node connectivity
    FieldContainer<int> faceToNode(face_vector.size(), numNodesPerFace);
    for(unsigned fct = 0; fct != face_vector.size(); fct++){
      std::list<long long>::iterator flit;
      int nct = 0;
      for(flit  = face_vector[fct]->local_node_ids.begin();
	  flit != face_vector[fct]->local_node_ids.end();
	  flit ++){
	faceToNode(fct,nct) = *flit-1;
	nct++;
      }
    }


   // Face to Edge connectivity
    FieldContainer<int> faceToEdge(face_vector.size(), numEdgesPerFace);
    FieldContainer<bool> faceDone(face_vector.size());
    for (int ielem = 0; ielem < numElems; ielem++){
       for (int iface = 0; iface < numFacesPerElem; iface++){
         if (!faceDone(elemToFace(ielem,iface))){
           for (int iedge = 0; iedge < numEdgesPerFace; iedge++){
              faceToEdge(elemToFace(ielem,iface),iedge) = 
                           elemToEdge(ielem,refFaceToEdge(iface,iedge));
              faceDone(elemToFace(ielem,iface))=1;
           }
         }
       }
    }   

    int numEdges = edge_vector.size();
    int numFaces = face_vector.size();

   // Calculate global edge and face numbering
    std::string doing_type;
    doing_type = "EDGES";
    calc_global_ids(edge_vector,
	       comm_node_ids,
	       node_comm_proc_ids, 
	       node_cmap_node_cnts,
	       num_node_comm_maps,
	       rank,
	       doing_type);


    doing_type = "FACES";
    calc_global_ids(face_vector,
	       comm_node_ids,
	       node_comm_proc_ids, 
	       node_cmap_node_cnts,
	       num_node_comm_maps,
	       rank,
	       doing_type);

  // Build list of owned global edge ids
    long long * globalEdgeIds = new long long[numEdges];
    bool * edgeIsOwned = new bool[numEdges];
    int numOwnedEdges=0;
    for (int i=0; i<numEdges; i++) {
        edgeIsOwned[i] = edge_vector[i]->owned;
        globalEdgeIds[i] = edge_vector[i]->global_id;
        if (edgeIsOwned[i]){
           numOwnedEdges++;
        }
     }
    int * ownedEdgeIds = new int[numOwnedEdges];
    int nedge=0;
    for (int i=0; i<numEdges; i++) {
        if (edgeIsOwned[i]){
           ownedEdgeIds[nedge]=(int)globalEdgeIds[i];
           nedge++;
        }
     }

  // Build list of owned global face ids
    long long * globalFaceIds = new long long[numFaces];
    bool * faceIsOwned = new bool[numFaces];
    int numOwnedFaces=0;
    for (int i=0; i<numFaces; i++) {
        faceIsOwned[i] = face_vector[i]->owned;
        globalFaceIds[i] = face_vector[i]->global_id;
        if (faceIsOwned[i]){
           numOwnedFaces++;
        }
     }
    int * ownedFaceIds = new int[numOwnedFaces];
    int nface=0;
    for (int i=0; i<numFaces; i++) {
        if (faceIsOwned[i]){
           ownedFaceIds[nface]=(int)globalFaceIds[i];
           nface++;
        }
     }

  // Calculate number of global edges and faces
    int numEdgesGlobal;
    int numFacesGlobal;
#ifdef HAVE_MPI
    Comm.SumAll(&numOwnedEdges,&numEdgesGlobal,1);
    Comm.SumAll(&numOwnedFaces,&numFacesGlobal,1);
#else
    numEdgesGlobal = numEdges;
    numFacesGlobal = numFaces;
#endif

 // Print mesh size information
  if (MyPID == 0) {
    std::cout << " Number of Elements: " << numElemsGlobal << " \n";
    std::cout << "    Number of Nodes: " << numNodesGlobal << " \n";
    std::cout << "    Number of Edges: " << numEdgesGlobal << " \n";
    std::cout << "    Number of Faces: " << numFacesGlobal << " \n\n";
  }

   // Container indicating whether a face is on the boundary (1-yes 0-no)
    FieldContainer<int> edgeOnBoundary(numEdges);
    FieldContainer<int> faceOnBoundary(numFaces);

 // Get boundary (side set) information
    long long * sideSetIds = new long long [numSideSets];
    long long numSidesInSet;
    long long numDFinSet;
    int numBndyFaces=0;
    im_ex_get_side_set_ids_l(id,sideSetIds);
    for (int i=0; i<numSideSets; i++) {
        im_ex_get_side_set_param_l(id,sideSetIds[i],&numSidesInSet,&numDFinSet);
        if (numSidesInSet > 0){
          long long * sideSetElemList = new long long [numSidesInSet];
          long long * sideSetSideList = new long long [numSidesInSet];
          im_ex_get_side_set_l(id,sideSetIds[i],sideSetElemList,sideSetSideList);
          for (int j=0; j<numSidesInSet; j++) {
             int iface = sideSetSideList[j]-1;
             faceOnBoundary(elemToFace(sideSetElemList[j]-1,iface))=1;
             numBndyFaces++;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),0))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),1))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),2))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),3))=1;
          }
          delete [] sideSetElemList;
          delete [] sideSetSideList;
       }
    }

    delete [] sideSetIds;

#ifdef RANDOM_GRID
/***********************************************************************************/
/* This block of code will create a randomly perturbed mesh by jiggling the node
   coordinates up to 1/4 of an edge length away from their initial location.
         !!!NOTE: THIS WILL ONLY WORK CORRECTLY ON A SINGLE PROCESSOR!!!           */
/***********************************************************************************/

   // Container indicating whether a node is on the boundary (1-yes 0-no)
    FieldContainer<int> nodeOnBoundary(numNodes);
     for (int i=0; i<numEdges; i++) {
        if (edgeOnBoundary(i)){
           nodeOnBoundary(edgeToNode(i,0))=1;
           nodeOnBoundary(edgeToNode(i,1))=1;
        }
     }

     // Side length assuming an initially regular grid
         double hx = nodeCoord(elemToNode(0,1),0)-nodeCoord(elemToNode(0,0),0);
         double hy = nodeCoord(elemToNode(0,3),1)-nodeCoord(elemToNode(0,0),1);
         double hz = nodeCoord(elemToNode(0,4),2)-nodeCoord(elemToNode(0,0),2);

     // Loop over nodes
       for (int inode = 0; inode < numNodes; inode++){
         if (!nodeOnBoundary(inode) & nodeIsOwned[inode]) {
           // random numbers between -1.0 and 1.0
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double rz = 2.0 * (double)rand()/RAND_MAX - 1.0;
           // limit variation to 1/4 edge length
            nodeCoord(inode,0) = nodeCoord(inode,0) + 0.125 * hx * rx;
            nodeCoord(inode,1) = nodeCoord(inode,1) + 0.125 * hy * ry;
            nodeCoord(inode,2) = nodeCoord(inode,2) + 0.125 * hz * rz;
            nodeCoordx[inode] = nodeCoord(inode,0);
            nodeCoordy[inode] = nodeCoord(inode,1);
            nodeCoordz[inode] = nodeCoord(inode,2);
         }
       }
#endif

/**********************************************************************************/
/********************************* GET CUBATURE ***********************************/
/**********************************************************************************/

   // Get numerical integration points and weights
    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(cellType, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);

 /**********************************************************************************/
 /*     Get numerical integration points and weights for hexahedron face           */
 /*                  (needed for rhs boundary term)                                */
 /**********************************************************************************/

    // Define topology of the face parametrization domain as [-1,1]x[-1,1]
    CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Define cubature
    DefaultCubatureFactory<double>  cubFactoryFace;
    Teuchos::RCP<Cubature<double> > hexFaceCubature = cubFactoryFace.create(paramQuadFace, 3);
    int cubFaceDim    = hexFaceCubature -> getDimension();
    int numFacePoints = hexFaceCubature -> getNumPoints();

    // Define storage for cubature points and weights on [-1,1]x[-1,1]
    FieldContainer<double> paramFaceWeights(numFacePoints);
    FieldContainer<double> paramFacePoints(numFacePoints,cubFaceDim);

    // Define storage for cubature points on workset faces
    hexFaceCubature -> getCubature(paramFacePoints, paramFaceWeights);

  if(MyPID==0) {std::cout << "Getting cubature                            "
                 << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}


/**********************************************************************************/
/*********************************** GET BASIS ************************************/
/**********************************************************************************/

   // Define basis 
    Basis_HCURL_HEX_I1_FEM<double, FieldContainer<double> > hexHCurlBasis;
    Basis_HDIV_HEX_I1_FEM<double, FieldContainer<double> > hexHDivBasis;
    Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;

    int numFieldsC = hexHCurlBasis.getCardinality();
    int numFieldsD = hexHDivBasis.getCardinality();
    int numFieldsG = hexHGradBasis.getCardinality();

  // Evaluate basis at cubature points
     FieldContainer<double> HCVals(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> HDVals(numFieldsD, numCubPoints, spaceDim); 
     FieldContainer<double> HDivs(numFieldsD, numCubPoints); 
     FieldContainer<double> HGVals(numFieldsG, numCubPoints); 

     hexHCurlBasis.getValues(HCVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(HDVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(HDivs, cubPoints, OPERATOR_DIV);
     hexHGradBasis.getValues(HGVals, cubPoints, OPERATOR_VALUE);

/**********************************************************************************/
/********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
/**********************************************************************************/

   // Define global epetra maps
    Epetra_Map globalMapG(-1,ownedNodes,ownedGIDs,0,Comm);
    Epetra_Map globalMapC(-1,numOwnedEdges,ownedEdgeIds,0,Comm);
    Epetra_Map globalMapD(-1,numOwnedFaces,ownedFaceIds,0,Comm);

   // Global arrays in Epetra format
    Epetra_FECrsMatrix MassMatrixC(Copy, globalMapC, numFieldsC);
    Epetra_FECrsMatrix MassMatrixD(Copy, globalMapD, numFieldsD);
    Epetra_FECrsMatrix MassMatrixG(Copy, globalMapG, numFieldsG);
    Epetra_FECrsMatrix StiffMatrixD(Copy, globalMapD, numFieldsD);
    Epetra_FEVector rhsVector(globalMapD);

 if(MyPID==0) {std::cout << "Build global maps                           "
                 << Time.ElapsedTime() << " sec \n";  Time.ResetStartTime();}


/**********************************************************************************/
/************************** OUTPUT CONNECTIVITY (FOR PLOTTING) ********************/
/**********************************************************************************/

  // Build the coordinate vectors for ML solver (including owned nodes only)
    Epetra_Vector Nx(globalMapG), Ny(globalMapG),Nz(globalMapG);
    for(int i=0,nlid=0;i<numNodes;i++)
      if(nodeIsOwned[i]) {
        Nx[nlid]=nodeCoordx[i];
        Ny[nlid]=nodeCoordy[i];
        Nz[nlid]=nodeCoordz[i];
        nlid++;
      }

#ifdef DUMP_DATA
    // Put element to node mapping in multivector for output
     Epetra_Map   globalMapElem(numElemsGlobal, numElems, 0, Comm);
     Epetra_MultiVector elem2node(globalMapElem, numNodesPerElem);
     for (int ielem=0; ielem<numElems; ielem++) {
        for (int inode=0; inode<numNodesPerElem; inode++) {
          elem2node[inode][ielem]=globalNodeIds[elemToNode(ielem,inode)];
        }
      }
     EpetraExt::MultiVectorToMatrixMarketFile("elem2node.dat",elem2node,0,0,false);

    // Put element to face mapping in multivector for output
     Epetra_MultiVector elem2face(globalMapElem, numFacesPerElem);
     for (int ielem=0; ielem<numElems; ielem++) {
        for (int iface=0; iface<numFacesPerElem; iface++) {
          elem2face[iface][ielem]=globalFaceIds[elemToFace(ielem,iface)];
        }
      }
     EpetraExt::MultiVectorToMatrixMarketFile("elem2face.dat",elem2face,0,0,false);

    // Put face to edge mapping in multivector for output
     Epetra_MultiVector face2edge(globalMapD, numEdgesPerFace);
     int ownedFace = 0;
     for (int iface=0; iface<numFaces; iface++) {
       if (faceIsOwned[iface]) {
         for (int iedge=0; iedge<numEdgesPerFace; iedge++) {
           face2edge[iedge][ownedFace]=globalNodeIds[faceToEdge(iface,iedge)];
         }
         ownedFace++;
        }
      }
     EpetraExt::MultiVectorToMatrixMarketFile("face2edge.dat",face2edge,0,0,false);

    // Put face to node mapping in multivector for output
     Epetra_MultiVector face2node(globalMapD, numNodesPerFace);
     int iFace = 0;
     for (int iface=0; iface<numFaces; iface++) {
       if (faceIsOwned[iface]) {
         for (int inode=0; inode<numNodesPerFace; inode++) {
           face2node[inode][iFace]=globalNodeIds[faceToNode(iface,inode)];
         }
         iFace++;
        }
      }
     EpetraExt::MultiVectorToMatrixMarketFile("face2node.dat",face2node,0,0,false);

    // Define multi-vector for cell face and edge signs (fill during cell loop)
     Epetra_MultiVector faceSign(globalMapElem, numFacesPerElem);
     Epetra_MultiVector edgeSign(globalMapElem, numEdgesPerElem);

  if(MyPID==0) {Time.ResetStartTime();}
#endif

/**********************************************************************************/
/*************************BUILD INCIDENCE MATRIX***********************************/
/**********************************************************************************/

  // Edge to node (DGrad) and face to edge (DCurl) incidence matrices
    Epetra_FECrsMatrix DCurl(Copy, globalMapD, 4);
    Epetra_FECrsMatrix DGrad(Copy, globalMapC, 2);

   // Edge to node incidence matrix
    double edgevals[2];
    edgevals[0]=-1.0; edgevals[1]=1.0;
    for (int j=0; j<numEdges; j++){
      if (edgeIsOwned[j]){
        int rowNum = globalEdgeIds[j];
        int colNum[2];
        colNum[0] = globalNodeIds[edgeToNode(j,0)];
        colNum[1] = globalNodeIds[edgeToNode(j,1)];
        DGrad.InsertGlobalValues(1, &rowNum, 2, colNum, edgevals);
      }
    }

   // Edge to face incidence matrix
   FieldContainer<bool> faceDone2(numFaces);
   for (int i=0; i<numElems; i++){
      for (int k=0; k<numFacesPerElem; k++){
         int iface=elemToFace(i,k);
         if (faceIsOwned[iface] && !faceDone2(iface)){
             double vals[4];
             int rowNum = globalFaceIds[iface];
             int colNum[4];

            for (int m=0; m<numEdgesPerFace; m++){
              colNum[m] = globalEdgeIds[faceToEdge(iface,m)];
              int indm = m+1;
              if (indm >= numEdgesPerFace) indm=0;
              if (edgeToNode(faceToEdge(iface,m),1) == edgeToNode(faceToEdge(iface,indm),0) ||
                 edgeToNode(faceToEdge(iface,m),1) == edgeToNode(faceToEdge(iface,indm),1)){
                 vals[m]=1.0;}
              else vals[m]=-1.0;

            // This is a convoluted way to account for edge orientations that 
            // may be incorrect on the local processor because the edge is
            // not owned by the local processor.
             int edgeIndex = -1;
             if (!edgeIsOwned[faceToEdge(iface,m)]){
                 for (int l=0; l<numEdgesPerElem; l++){
                    if (faceToEdge(iface,m)==elemToEdge(i,l)) 
                       edgeIndex=l;
                }
             }
               if (edgeIndex != -1 && edgeIndex < 8){
                 if (edgeIndex < 4 && faceIsOwned[elemToFace(i,4)]){
                   vals[m]=-1.0*vals[m];                 
                 }
                 else if (edgeIndex > 3 && faceIsOwned[elemToFace(i,5)]){
                   vals[m]=-1.0*vals[m];
                 }
               }
            } // end loop over face edges

           DCurl.InsertGlobalValues(1, &rowNum, 4, colNum, vals);
           faceDone2(iface)=1;

       } // end if face is owned and face not done
    } // end loop over element faces
  } // end loop over elements


    DGrad.GlobalAssemble(globalMapG,globalMapC);
    DGrad.FillComplete(MassMatrixG.RowMap(),MassMatrixC.RowMap());

    DCurl.GlobalAssemble(globalMapC,globalMapD);
    DCurl.FillComplete(MassMatrixC.RowMap(),MassMatrixD.RowMap());

  if(MyPID==0) {std::cout << "Building incidence matrices                 "
                 << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}


/**********************************************************************************/
/******************** INHOMOGENEOUS BOUNDARY CONDITIONS ***************************/
/**********************************************************************************/


    FieldContainer<double> bndyFaceVal(numBndyFaces);
    FieldContainer<int>    bndyFaceToFace(numFaces);
    FieldContainer<double> refFacePoints(numFacePoints,spaceDim);
    FieldContainer<double> bndyFacePoints(1,numFacePoints,spaceDim);
    FieldContainer<double> bndyFaceJacobians(1,numFacePoints,spaceDim,spaceDim);
    FieldContainer<double> faceNorm(1,numFacePoints,spaceDim);
    FieldContainer<double> uDotNormal(1,numFacePoints);
    FieldContainer<double> uFace(numFacePoints,spaceDim);
    FieldContainer<double> nodes(1, numNodesPerElem, spaceDim);

    int ibface=0;
    // Evaluate normal at face quadrature points
    for (int ielem=0; ielem<numElems; ielem++) {
       for (int inode=0; inode<numNodesPerElem; inode++) {
         nodes(0,inode,0) = nodeCoord(elemToNode(ielem,inode),0);
         nodes(0,inode,1) = nodeCoord(elemToNode(ielem,inode),1);
         nodes(0,inode,2) = nodeCoord(elemToNode(ielem,inode),2);
       }
       for (int iface=0; iface<numFacesPerElem; iface++){
          if(faceOnBoundary(elemToFace(ielem,iface))){

          // map evaluation points from reference face to reference cell
             IntrepidCTools::mapToReferenceSubcell(refFacePoints,
                                   paramFacePoints,
                                   2, iface, cellType);

          // calculate Jacobian
             IntrepidCTools::setJacobian(bndyFaceJacobians, refFacePoints,
                         nodes, cellType);

          // map evaluation points from reference cell to physical cell
             IntrepidCTools::mapToPhysicalFrame(bndyFacePoints,
                                refFacePoints,
                                nodes, cellType);

          // Compute face normals
             IntrepidCTools::getPhysicalFaceNormals(faceNorm,
                                              bndyFaceJacobians,
                                              iface, cellType);

          // evaluate exact solution and dot with normal
           for(int nPt = 0; nPt < numFacePoints; nPt++){

             double x = bndyFacePoints(0, nPt, 0);
             double y = bndyFacePoints(0, nPt, 1);
             double z = bndyFacePoints(0, nPt, 2);

             evalu(uFace(nPt,0), uFace(nPt,1), uFace(nPt,2), x, y, z);
             uDotNormal(0,nPt)=(uFace(nPt,0)*faceNorm(0,nPt,0)+uFace(nPt,1)*faceNorm(0,nPt,1)+uFace(nPt,2)*faceNorm(0,nPt,2));
           }

          // integrate
           for(int nPt = 0; nPt < numFacePoints; nPt++){
             bndyFaceVal(ibface)=bndyFaceVal(ibface)+uDotNormal(0,nPt)*paramFaceWeights(nPt);
           }
           bndyFaceToFace(elemToFace(ielem,iface))=ibface;
           ibface++;
         } // end if face on boundary

       } // end loop over element faces

    } // end loop over elements


   // Count boundary faces
     int numBCFaces=0;
     for (int i=0; i<numFaces; i++){
         if (faceOnBoundary(i) && faceIsOwned[i]){
             numBCFaces++;
         }
      }

    // Vector for use in applying BCs
     Epetra_MultiVector v(globalMapD,true);
     v.PutScalar(0.0);

   // Set v to boundary values on Dirichlet faces
     int * BCFaces = new int [numBCFaces];
     int indbc=0;
     int indOwned=0;
     for (int i=0; i<numFaces; i++){
       if (faceIsOwned[i]){
        if (faceOnBoundary(i)){
           BCFaces[indbc]=indOwned;
           indbc++;
           v[0][indOwned]=bndyFaceVal(bndyFaceToFace(i));
          }
         indOwned++;
        }
     }

  if(MyPID==0) {std::cout << "Boundary Condition Setup                    "
                 << Time.ElapsedTime() << " sec \n\n"; Time.ResetStartTime();}

/**********************************************************************************/
/******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
/**********************************************************************************/

// Define desired workset size and count how many worksets there are on this processor's mesh block
  int desiredWorksetSize = numElems;                      // change to desired workset size!
  //int desiredWorksetSize = 100;                      // change to desired workset size!
  int numWorksets        = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksets*desiredWorksetSize < numElems) numWorksets += 1;

 if (MyPID == 0) {
    std::cout << "Building discretization matrix and right hand side... \n\n";
    std::cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
    std::cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
    Time.ResetStartTime();
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
     FieldContainer<double> cellWorkset(worksetSize, numNodesPerElem, spaceDim);
     FieldContainer<double> worksetEdgeSigns(worksetSize, numEdgesPerElem);
     FieldContainer<double> worksetFaceSigns(worksetSize, numFacesPerElem);

   // Copy coordinates into cell workset
    int cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Physical cell coordinates
       for (int inode=0; inode<numNodesPerElem; inode++) {
         cellWorkset(cellCounter,inode,0) = nodeCoord(elemToNode(cell,inode),0);
         cellWorkset(cellCounter,inode,1) = nodeCoord(elemToNode(cell,inode),1);
         cellWorkset(cellCounter,inode,2) = nodeCoord(elemToNode(cell,inode),2);
       }

     // Face signs
      for (int iface=0; iface<numFacesPerElem; iface++) {
         worksetFaceSigns(cellCounter,iface) = -1.0;
         for (int i=0; i<numNodesPerFace; i++) {
           int indf=i+1;
           if (indf >= numNodesPerFace) indf=0;
           if (elemToNode(cell,refFaceToNode(iface,0))==faceToNode(elemToFace(cell,iface),i) &&
               elemToNode(cell,refFaceToNode(iface,1))==faceToNode(elemToFace(cell,iface),indf))
                worksetFaceSigns(cellCounter,iface) = 1.0;
          }
           if (!faceIsOwned[elemToFace(cell,iface)]){
              worksetFaceSigns(cellCounter,iface)=-1.0*worksetFaceSigns(cellCounter,iface);
           }
       } 

#ifdef DUMP_DATA
       for (int iface=0; iface<numFacesPerElem; iface++) {
          faceSign[iface][cell] = worksetFaceSigns(cellCounter,iface);
       }
#endif

      // Edge signs
       for (int iedge=0; iedge<numEdgesPerElem; iedge++) {
          if (elemToNode(cell,refEdgeToNode(iedge,0))==edgeToNode(elemToEdge(cell,iedge),0) &&
              elemToNode(cell,refEdgeToNode(iedge,1))==edgeToNode(elemToEdge(cell,iedge),1))
              worksetEdgeSigns(cellCounter,iedge) = 1.0;
          else
              worksetEdgeSigns(cellCounter,iedge) = -1.0;
        }

      // modify signs for edges that are owned by another processor
      // (Note: this is particular to the numbering of edges used in Pamgen!!)
       if (!faceIsOwned[elemToFace(cell,0)]) {
            worksetEdgeSigns(cellCounter,0)=-1.0*worksetEdgeSigns(cellCounter,0);
            worksetEdgeSigns(cellCounter,4)=-1.0*worksetEdgeSigns(cellCounter,4);
        }
       if (!faceIsOwned[elemToFace(cell,1)]) {
            worksetEdgeSigns(cellCounter,1)=-1.0*worksetEdgeSigns(cellCounter,1);
            worksetEdgeSigns(cellCounter,5)=-1.0*worksetEdgeSigns(cellCounter,5);
        }
       if (!faceIsOwned[elemToFace(cell,2)]) {
            worksetEdgeSigns(cellCounter,2)=-1.0*worksetEdgeSigns(cellCounter,2);
            worksetEdgeSigns(cellCounter,6)=-1.0*worksetEdgeSigns(cellCounter,6);
        }
       if (!faceIsOwned[elemToFace(cell,3)]) {
            worksetEdgeSigns(cellCounter,3)=-1.0*worksetEdgeSigns(cellCounter,3);
            worksetEdgeSigns(cellCounter,7)=-1.0*worksetEdgeSigns(cellCounter,7);
        }

#ifdef DUMP_DATA
       for (int iedge=0; iedge<numEdgesPerElem; iedge++) {
          edgeSign[iedge][cell] = worksetEdgeSigns(cellCounter,iedge);
       }
#endif
        cellCounter++;

     } // end loop over workset cells

 /**********************************************************************************/
 /*                                Allocate arrays                                 */
 /**********************************************************************************/

   // Containers for Jacobian
    FieldContainer<double> worksetJacobian (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobInv (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet (worksetSize, numCubPoints);

   // Container for cubature points in physical space
    FieldContainer<double> worksetCubPoints (worksetSize,numCubPoints, cubDim);

   // Containers for element HCURL mass matrix
    FieldContainer<double> massMatrixHCurl           (worksetSize, numFieldsC, numFieldsC);
    FieldContainer<double> weightedMeasure           (worksetSize, numCubPoints);
    FieldContainer<double> weightedMeasureMu         (worksetSize, numCubPoints);
    FieldContainer<double> HCValsTransformed         (worksetSize, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> HCValsTransformedWeighted (worksetSize, numFieldsC, numCubPoints, spaceDim);

   // Containers for element HDIV mass matrix
    FieldContainer<double> massMatrixHDiv            (worksetSize, numFieldsD, numFieldsD);
    FieldContainer<double> HDValsTransformed         (worksetSize, numFieldsD, numCubPoints, spaceDim);
    FieldContainer<double> HDValsTransformedWeighted (worksetSize, numFieldsD, numCubPoints, spaceDim);

   // Containers for element HDIV stiffness matrix
    FieldContainer<double> stiffMatrixHDiv           (worksetSize, numFieldsD, numFieldsD);
    FieldContainer<double> HDivsTransformed          (worksetSize, numFieldsD, numCubPoints);
    FieldContainer<double> HDivsTransformedWeighted  (worksetSize, numFieldsD, numCubPoints);

   // Containers for element HGRAD mass matrix
    FieldContainer<double> massMatrixHGrad           (worksetSize, numFieldsG, numFieldsG);
    FieldContainer<double> HGValsTransformed         (worksetSize, numFieldsG, numCubPoints);
    FieldContainer<double> HGValsTransformedWeighted (worksetSize, numFieldsG, numCubPoints);

   // Containers for right hand side vectors
    FieldContainer<double> rhsDatag           (worksetSize, numCubPoints, cubDim);
    FieldContainer<double> rhsDatah           (worksetSize, numCubPoints);
    FieldContainer<double> gD                 (worksetSize, numFieldsD);
    FieldContainer<double> hD                 (worksetSize, numFieldsD);

  // Containers for right hand side boundary term
    FieldContainer<double> gDBoundary         (1, numFieldsD);
    FieldContainer<double> refFacePoints      (numFacePoints,spaceDim);
    FieldContainer<double> cellNodes          (1, numNodesPerElem, spaceDim);
    FieldContainer<double> worksetFacePoints  (1, numFacePoints, spaceDim);
    FieldContainer<double> faceJacobians      (1, numFacePoints, spaceDim, spaceDim);
    FieldContainer<double> faceJacobInv       (1, numFacePoints, spaceDim, spaceDim);
    FieldContainer<double> faceJacobDet       (1, numFacePoints);
    FieldContainer<double> faceNormal         (1, numFacePoints, spaceDim);
    FieldContainer<double> bcFaceDVals        (numFieldsD, numFacePoints, spaceDim);
    FieldContainer<double> bcDValsTransformed (1, numFieldsD, numFacePoints, spaceDim);
    FieldContainer<double> curluFace          (1, numFacePoints, spaceDim);
    FieldContainer<double> bcDataCrossField   (1, numFieldsD, numFacePoints, spaceDim);
    FieldContainer<double> bcFaceSigns        (1, numFieldsD);


 /**********************************************************************************/
 /*                                Calculate Jacobians                             */
 /**********************************************************************************/

      IntrepidCTools::setJacobian   (worksetJacobian, cubPoints, cellWorkset, cellType);
      IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
      IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

   if(MyPID==0) {std::cout << "Calculate Jacobians                         "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


 /**********************************************************************************/
 /*                          Compute HCURL Mass Matrix                             */
 /**********************************************************************************/

     // transform to physical coordinates 
      IntrepidFSTools::HCURLtransformVALUE<double>(HCValsTransformed, worksetJacobInv, 
                                   HCVals);

     // compute weighted measure
      IntrepidFSTools::computeCellMeasure<double>(weightedMeasure, worksetJacobDet, cubWeights);

     // combine mu value with weighted measure

      cellCounter = 0;
      for(int cell = worksetBegin; cell < worksetEnd; cell++){
         for (int nPt = 0; nPt < numCubPoints; nPt++){
          weightedMeasureMu(cellCounter,nPt) = weightedMeasure(cellCounter,nPt) * muVal(cell);
        }
        cellCounter++;
      }

     // multiply by weighted measure
      IntrepidFSTools::multiplyMeasure<double>(HCValsTransformedWeighted,
                                               weightedMeasureMu, HCValsTransformed);

     // integrate to compute element mass matrix
      IntrepidFSTools::integrate<double>(massMatrixHCurl,
                                         HCValsTransformed, HCValsTransformedWeighted,
                                         COMP_BLAS);

     // apply edge signs
      IntrepidFSTools::applyLeftFieldSigns<double>(massMatrixHCurl, worksetEdgeSigns);
      IntrepidFSTools::applyRightFieldSigns<double>(massMatrixHCurl, worksetEdgeSigns);

   if(MyPID==0) {std::cout << "Compute HCURL Mass Matrix                   "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

 /**********************************************************************************/
 /*                         Compute HDiv Mass Matrix                               */
 /**********************************************************************************/

     // transform to physical coordinates 
      IntrepidFSTools::HDIVtransformVALUE<double>(HDValsTransformed, worksetJacobian, 
                                                  worksetJacobDet, HDVals);

     // multiply by weighted measure
      IntrepidFSTools::multiplyMeasure<double>(HDValsTransformedWeighted,
                                               weightedMeasure, HDValsTransformed);

     // integrate to compute element mass matrix
      IntrepidFSTools::integrate<double>(massMatrixHDiv,
                                         HDValsTransformed, HDValsTransformedWeighted,
                                         COMP_BLAS);

     // apply face signs
      IntrepidFSTools::applyLeftFieldSigns<double>(massMatrixHDiv, worksetFaceSigns);
      IntrepidFSTools::applyRightFieldSigns<double>(massMatrixHDiv, worksetFaceSigns);


   if(MyPID==0) {std::cout << "Compute HDIV Mass Matrix                    "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

 /**********************************************************************************/
 /*                         Compute HDiv Stiffness Matrix                          */
 /**********************************************************************************/

      // transform to physical coordinates 
      IntrepidFSTools::HDIVtransformDIV<double>(HDivsTransformed, worksetJacobDet,
                                                HDivs);

     // multiply by weighted measure
      IntrepidFSTools::multiplyMeasure<double>(HDivsTransformedWeighted,
                                               weightedMeasure, HDivsTransformed);

     // integrate to compute element stiffness matrix
      IntrepidFSTools::integrate<double>(stiffMatrixHDiv,
                                         HDivsTransformed, HDivsTransformedWeighted,
                                         COMP_BLAS);

     // apply face signs
      IntrepidFSTools::applyLeftFieldSigns<double>(stiffMatrixHDiv, worksetFaceSigns);
      IntrepidFSTools::applyRightFieldSigns<double>(stiffMatrixHDiv, worksetFaceSigns);

  if(MyPID==0) {std::cout << "Compute HDiv Stiffness Matrix               "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

 /**********************************************************************************/
 /*                         Compute HGrad Mass Matrix                              */
 /**********************************************************************************/

     // transform to physical coordinates
      IntrepidFSTools::HGRADtransformVALUE<double>(HGValsTransformed, HGVals);

     // multiply values with weighted measure
      IntrepidFSTools::multiplyMeasure<double>(HGValsTransformedWeighted,
                                               weightedMeasure, HGValsTransformed);

     // integrate to compute element mass matrix
      IntrepidFSTools::integrate<double>(massMatrixHGrad,
                                         HGValsTransformed, HGValsTransformedWeighted, 
                                         COMP_BLAS);

  if(MyPID==0) {std::cout << "Compute HGRAD Mass Matrix                   "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

 /**********************************************************************************/
 /*                          Build Right Hand Side                                 */
 /**********************************************************************************/

      // transform integration points to physical points
       IntrepidCTools::mapToPhysicalFrame(worksetCubPoints, cubPoints, cellWorkset, cellType);

      // evaluate right hand side functions at physical points
       for(int cell = worksetBegin; cell < worksetEnd; cell++){

        // Compute cell ordinal relative to the current workset
         int worksetCellOrdinal = cell - worksetBegin;

           for (int nPt = 0; nPt < numCubPoints; nPt++){

              double x = worksetCubPoints(worksetCellOrdinal,nPt,0);
              double y = worksetCubPoints(worksetCellOrdinal,nPt,1);
              double z = worksetCubPoints(worksetCellOrdinal,nPt,2);
              double du1, du2, du3;

              evalCurlCurlu(du1, du2, du3, x, y, z, muVal(cell));
              rhsDatag(worksetCellOrdinal,nPt,0) = du1;
              rhsDatag(worksetCellOrdinal,nPt,1) = du2;
              rhsDatag(worksetCellOrdinal,nPt,2) = du3;
         
              rhsDatah(worksetCellOrdinal,nPt) = evalDivu(x, y, z);
           }

         } // end workset cell loop

        // integrate (g,curl w) term
         IntrepidFSTools::integrate<double>(gD, rhsDatag, HDValsTransformedWeighted,
                                            COMP_BLAS);

        // integrate (h,div w) term
         IntrepidFSTools::integrate<double>(hD, rhsDatah, HDivsTransformedWeighted,
                                            COMP_BLAS);

        // apply signs
         IntrepidFSTools::applyFieldSigns<double>(gD, worksetFaceSigns);
         IntrepidFSTools::applyFieldSigns<double>(hD, worksetFaceSigns);


        // calculate RHS boundary term
        for(int cell = worksetBegin; cell < worksetEnd; cell++){

        // Compute cell ordinal relative to the current workset
          int worksetCellOrdinal = cell - worksetBegin;

        // evaluate boundary term
          for (int iface = 0; iface < numFacesPerElem; iface++){
             if (faceOnBoundary(elemToFace(cell,iface))){

               // cell nodal coordinates
                for (int inode =0; inode < numNodesPerElem; inode++){
                   cellNodes(0,inode,0) = cellWorkset(worksetCellOrdinal,inode,0);
                   cellNodes(0,inode,1) = cellWorkset(worksetCellOrdinal,inode,1);
                   cellNodes(0,inode,2) = cellWorkset(worksetCellOrdinal,inode,2);
                }

               // cell face signs
                for (int jface =0; jface < numFacesPerElem; jface++){
                   bcFaceSigns(0,jface) = worksetFaceSigns(worksetCellOrdinal,jface);
                }

              // map Gauss points on quad to reference face: paramFacePoints -> refFacePoints
                IntrepidCTools::mapToReferenceSubcell(refFacePoints,
                                   paramFacePoints,
                                   2, iface, cellType);

              // get basis values at points on reference cell
                hexHDivBasis.getValues(bcFaceDVals, refFacePoints, OPERATOR_VALUE);

              // compute Jacobians at Gauss pts. on reference face for all parent cells
                IntrepidCTools::setJacobian(faceJacobians, refFacePoints,
                                            cellNodes, cellType);
                IntrepidCTools::setJacobianDet(faceJacobDet, faceJacobians);

              // transform to physical coordinates
                IntrepidFSTools::HDIVtransformVALUE<double>(bcDValsTransformed, faceJacobians,
                                                            faceJacobDet, bcFaceDVals);

              // map Gauss points on quad from ref. face to face workset: refFacePoints -> worksetFacePoints
                IntrepidCTools::mapToPhysicalFrame(worksetFacePoints,
                                                   refFacePoints,
                                                   cellNodes, cellType);

               // Compute face normals
                IntrepidCTools::getPhysicalFaceNormals(faceNormal,
                                                       faceJacobians,
                                                       iface, cellType);

               // evaluate curl u at face points
                for(int nPt = 0; nPt < numFacePoints; nPt++){

                   double x = worksetFacePoints(0, nPt, 0);
                   double y = worksetFacePoints(0, nPt, 1);
                   double z = worksetFacePoints(0, nPt, 2);

                   evalCurlu(curluFace(0,nPt,0), curluFace(0,nPt,1), 
                             curluFace(0,nPt,2), x, y, z, muVal(cell));
                 }

               // compute the cross product of curluFace with basis and multiply by weights
                for (int nF = 0; nF < numFieldsD; nF++){
                   for(int nPt = 0; nPt < numFacePoints; nPt++){
                     bcDataCrossField(0,nF,nPt,0) = (curluFace(0,nPt,1)*bcDValsTransformed(0,nF,nPt,2)
                                 - curluFace(0,nPt,2)*bcDValsTransformed(0,nF,nPt,1))
                                  * paramFaceWeights(nPt);
                     bcDataCrossField(0,nF,nPt,1) = (curluFace(0,nPt,2)*bcDValsTransformed(0,nF,nPt,0)
                                 - curluFace(0,nPt,0)*bcDValsTransformed(0,nF,nPt,2))
                                  * paramFaceWeights(nPt);
                     bcDataCrossField(0,nF,nPt,2) = (curluFace(0,nPt,0)*bcDValsTransformed(0,nF,nPt,1)
                                 - curluFace(0,nPt,1)*bcDValsTransformed(0,nF,nPt,0))
                                  *paramFaceWeights(nPt);
                    } //nPt
                 } //nF

                // integrate
                 IntrepidFSTools::integrate<double>(gDBoundary, faceNormal, bcDataCrossField,
                                                    COMP_CPP);

                // apply signs
                 IntrepidFSTools::applyFieldSigns<double>(gDBoundary, bcFaceSigns);

                // add into  gD term
                 for (int nF = 0; nF < numFieldsD; nF++){
                   gD(worksetCellOrdinal,nF) = gD(worksetCellOrdinal,nF) - gDBoundary(0,nF);
                 }

              } // if faceOnBoundary

          } // numFaces

       }// *** workset cell loop **



  if(MyPID==0) {std::cout << "Compute right-hand side                     "
                  << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

 /**********************************************************************************/
 /*                         Assemble into Global Matrix                            */
 /**********************************************************************************/

   // Loop over workset cells
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;


      /*** Assemble H(grad) mass matrix ***/

      // loop over nodes for matrix row
      for (int cellNodeRow = 0; cellNodeRow < numFieldsG; cellNodeRow++){

        int localNodeRow  = elemToNode(cell, cellNodeRow);
        int globalNodeRow = globalNodeIds[localNodeRow];

       // loop over nodes for matrix column
        for (int cellNodeCol = 0; cellNodeCol < numFieldsG; cellNodeCol++){

          int localNodeCol  = elemToNode(cell, cellNodeCol);
          int globalNodeCol = globalNodeIds[localNodeCol];
          double massGContribution = massMatrixHGrad(worksetCellOrdinal, cellNodeRow, cellNodeCol);

          MassMatrixG.InsertGlobalValues(1, &globalNodeRow, 1, &globalNodeCol, &massGContribution);

        }// *** cell node col loop ***
      }// *** cell node row loop ***


      /*** Assemble H(curl) mass matrix ***/

      // loop over edges for matrix row
      for (int cellEdgeRow = 0; cellEdgeRow < numFieldsC; cellEdgeRow++){

        int localEdgeRow  = elemToEdge(cell, cellEdgeRow);
        int globalEdgeRow = globalEdgeIds[localEdgeRow];

       // loop over edges for matrix column
        for (int cellEdgeCol = 0; cellEdgeCol < numFieldsC; cellEdgeCol++){

          int localEdgeCol  = elemToEdge(cell, cellEdgeCol);
          int globalEdgeCol = globalEdgeIds[localEdgeCol];

          double massCContribution  = massMatrixHCurl (worksetCellOrdinal, cellEdgeRow, cellEdgeCol);

          MassMatrixC.InsertGlobalValues (1, &globalEdgeRow, 1, &globalEdgeCol, &massCContribution);

        }// *** cell edge col loop ***
      }// *** cell edge row loop ***


      /*** Assemble H(div) mass matrix, stiffness matrix and right-hand side ***/

      // loop over faces for matrix row
      for (int cellFaceRow = 0; cellFaceRow < numFieldsD; cellFaceRow++){

        int localFaceRow  = elemToFace(cell, cellFaceRow);
        int globalFaceRow = globalFaceIds[localFaceRow];
        double rhsContribution = gD(worksetCellOrdinal, cellFaceRow) + hD(worksetCellOrdinal, cellFaceRow);

        rhsVector.SumIntoGlobalValues(1, &globalFaceRow, &rhsContribution);

       // loop over faces for matrix column
        for (int cellFaceCol = 0; cellFaceCol < numFieldsD; cellFaceCol++){

          int localFaceCol  = elemToFace(cell, cellFaceCol);
          int globalFaceCol = globalFaceIds[localFaceCol];

          double massDContribution  = massMatrixHDiv (worksetCellOrdinal, cellFaceRow, cellFaceCol);
          double stiffDContribution = stiffMatrixHDiv(worksetCellOrdinal, cellFaceRow, cellFaceCol);

          MassMatrixD.InsertGlobalValues (1, &globalFaceRow, 1, &globalFaceCol, &massDContribution);
          StiffMatrixD.InsertGlobalValues(1, &globalFaceRow, 1, &globalFaceCol, &stiffDContribution);


        }// *** cell face col loop ***
      }// *** cell face row loop ***

    }// *** workset cell loop **
  }// *** workset loop ***

  if(MyPID==0) {std::cout << "Assemble Matrices                           "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


/**********************************************************************************/
/********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
/**********************************************************************************/

    // Assemble over multiple processors, if necessary
    MassMatrixG.GlobalAssemble();  MassMatrixG.FillComplete();
    MassMatrixC.GlobalAssemble();  MassMatrixC.FillComplete();
    MassMatrixD.GlobalAssemble();  MassMatrixD.FillComplete();
    StiffMatrixD.GlobalAssemble(); StiffMatrixD.FillComplete();
    rhsVector.GlobalAssemble();


  if(MyPID==0) {std::cout << "Global assembly                             "
                 << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


#ifdef DUMP_DATA
    // Node Coordinates
    EpetraExt::VectorToMatrixMarketFile("coordx.dat",Nx,0,0,false);
    EpetraExt::VectorToMatrixMarketFile("coordy.dat",Ny,0,0,false);
    EpetraExt::VectorToMatrixMarketFile("coordz.dat",Nz,0,0,false);

    // Edge signs
    EpetraExt::MultiVectorToMatrixMarketFile("edge_signs.dat",edgeSign,0,0,false);

    // Face signs
    EpetraExt::MultiVectorToMatrixMarketFile("face_signs.dat",faceSign,0,0,false);
#endif

/**********************************************************************************/
/*********************** ADJUST MATRICES AND RHS FOR BCs **************************/
/**********************************************************************************/

  // Build the inverse diagonal for MassMatrixC
  // We do this by computing the absolute rowsum and inverting that.
    Epetra_CrsMatrix MassMatrixCinv(Copy,MassMatrixC.RowMap(),MassMatrixC.RowMap(),1);
    Epetra_Vector DiagC(MassMatrixC.RowMap());
    Epetra_Vector temp(MassMatrixC.RowMap());
    temp.PutScalar(1.0);
    Multiply_Abs(MassMatrixC,temp,DiagC);

    for(int i=0; i<DiagC.MyLength(); i++) {
      DiagC[i]=1.0/DiagC[i];
    }
    for(int i=0; i<DiagC.MyLength(); i++) {
      int CID=MassMatrixC.GCID(i);
      MassMatrixCinv.InsertGlobalValues(MassMatrixC.GRID(i),1,&(DiagC[i]),&CID);
    }
    MassMatrixCinv.FillComplete();

   // Set value to zero on diagonal that corresponds to boundary edge
    for(int i=0;i<numEdges;i++) {
       if (edgeOnBoundary(i)){
          double val=0.0;
          int index = globalEdgeIds[i];
          MassMatrixCinv.ReplaceGlobalValues(index,1,&val,&index);
       }
    }

   // Get the full operator matrix
   ML_Epetra::ML_RefMaxwell_11_Operator MatrixD(StiffMatrixD,DCurl,MassMatrixCinv,MassMatrixD);

  // Apply it to v
   Epetra_MultiVector rhsDir(globalMapD,true);
   MatrixD.Apply(v,rhsDir);

   // Update right-hand side
   rhsVector.Update(-1.0,rhsDir,1.0);

   // Update rhs values on Dirichlet edges
     indbc=0;
     int iOwned=0;
     for (int i=0; i<numFaces; i++){
       if (faceIsOwned[i]){
        if (faceOnBoundary(i)){
           indbc++;
           rhsVector[0][iOwned]=bndyFaceVal(bndyFaceToFace(i));
          }
         iOwned++;
        }
     }

   // Zero out rows and columns of stiffness and mass matrix corresponding to Dirichlet faces
   //  and add one to diagonal.
     ML_Epetra::Apply_OAZToMatrix(BCFaces, numBCFaces, StiffMatrixD);
     ML_Epetra::Apply_OAZToMatrix(BCFaces, numBCFaces, MassMatrixD);

     delete [] BCFaces;

   // Build Face-Node Incidence Matrix
     Epetra_CrsMatrix DGrad1(DGrad);
     Epetra_CrsMatrix DCurl1(DCurl);
     Epetra_CrsMatrix FaceNode(Copy,globalMapD,0);
     DGrad1.PutScalar(1.0);
     DCurl1.PutScalar(1.0);
     EpetraExt::MatrixMatrix::Multiply(DCurl1,false,DGrad1,false,FaceNode);
     FaceNode.PutScalar(1.0);

 if(MyPID==0) {std::cout << "Adjust global matrix and rhs due to BCs     " << Time.ElapsedTime()
                  << " sec \n"; Time.ResetStartTime();}

#ifdef DUMP_DATA
  // Dump matrices to disk
   EpetraExt::RowMatrixToMatlabFile("mag_m0_matrix.dat",MassMatrixG);
   EpetraExt::RowMatrixToMatlabFile("mag_m1_matrix.dat",MassMatrixC);
   EpetraExt::RowMatrixToMatlabFile("mag_m1inv_matrix.dat",MassMatrixCinv);
   EpetraExt::RowMatrixToMatlabFile("mag_m2_matrix.dat",MassMatrixD);
   EpetraExt::RowMatrixToMatlabFile("mag_k2_matrix.dat",StiffMatrixD);
   EpetraExt::RowMatrixToMatlabFile("mag_t0_matrix.dat",DGrad);
   EpetraExt::RowMatrixToMatlabFile("mag_t1_matrix.dat",DCurl);
   EpetraExt::RowMatrixToMatlabFile("mag_fn_matrix.dat",FaceNode);
   EpetraExt::MultiVectorToMatrixMarketFile("mag_rhs2.dat",rhsVector,0,0,false);
   EpetraExt::VectorToMatrixMarketFile("diagc.dat",DiagC,0,0,false);
#endif


/**********************************************************************************/
/*********************************** SOLVE ****************************************/
/**********************************************************************************/

   // Parameter list for ML
   Teuchos::ParameterList MLList,dummy,dummy2;
   double TotalErrorResidual=0, TotalErrorExactSol=0;   
   ML_Epetra::SetDefaultsRefMaxwell(MLList);
   Teuchos::ParameterList MLList2=MLList.get("refmaxwell: 11list",MLList);
   MLList2.set("aggregation: type","Uncoupled-MIS");
   MLList2.set("x-coordinates",&Nx[0]);
   MLList2.set("y-coordinates",&Ny[0]);
   MLList2.set("z-coordinates",&Nz[0]);   
   MLList2.set("ML output",10);
   MLList2.set("smoother: sweeps (level 0)",3);
   MLList2.set("smoother: sweeps",3);
   MLList2.set("smoother: type","Chebyshev");
   MLList2.set("eigen-analysis: type", "power-method");
   MLList2.get("edge matrix free: coarse",dummy);
   ML_Epetra::SetDefaults("SA",dummy,0,0,false);
   dummy.set("PDE equations",3);
   dummy.set("ML output",10);
   dummy.set("smoother: sweeps",3);
   dummy.set("smoother: type","Chebyshev");
   dummy.set("aggregation: type","Uncoupled-MIS");
   dummy.set("smoother: pre or post","both");
   dummy.set("max levels",10);
   dummy.set("coarse: type","Amesos-KLU");
   dummy.set("repartition: enable",1);
   dummy.set("repartition: min per proc",1000);
   dummy.set("repartition: max min ratio",1.4);
   dummy.set("repartition: Zoltan dimensions",3);
   dummy.set("x-coordinates",&Nx[0]);
   dummy.set("y-coordinates",&Ny[0]);
   dummy.set("z-coordinates",&Nz[0]);   
   MLList2.set("face matrix free: coarse",dummy);
   MLList2.set("edge matrix free: coarse","disabled");
   
  if (MyPID == 0) {
   cout<<MLList2<<endl;
  }

   Epetra_FEVector xh(rhsVector);  
   MassMatrixG.SetLabel("M0");
   MassMatrixC.SetLabel("M1");
   MassMatrixD.SetLabel("M2");
   StiffMatrixD.SetLabel("D2");
   DGrad.SetLabel("D0");
   DCurl.SetLabel("D1");
   MassMatrixCinv.SetLabel("M1^{-1}");
   
   char probType[12] = "div_lsfem";

   TestMultiLevelPreconditioner_DivLSFEM(probType, MLList2, StiffMatrixD,
					 DGrad, DCurl, FaceNode,
					 MassMatrixC, MassMatrixCinv,
					 MassMatrixD, xh, rhsVector,
					 TotalErrorResidual, TotalErrorExactSol);

/**********************************************************************************/
/**************************** CALCULATE ERROR *************************************/
/**********************************************************************************/

  if (MyPID == 0) {Time.ResetStartTime();}

     double L2err = 0.0;
     double HDiverr = 0.0;
     double Linferr = 0.0;
     double L2errTot = 0.0;
     double HDiverrTot = 0.0;
     double LinferrTot = 0.0;

#ifdef HAVE_MPI
   // Import solution onto current processor
     Epetra_Map  solnMap(numFacesGlobal, numFacesGlobal, 0, Comm);
     Epetra_Import  solnImporter(solnMap, globalMapD);
     Epetra_FEVector  uCoeff(solnMap);
     uCoeff.Import(xh, solnImporter, Insert);
#endif

     int numCells = 1;
     FieldContainer<double> hexFaceSigns(numCells, numFieldsD);
     FieldContainer<double> hexNodes(numCells, numFieldsG, spaceDim);

   // Get cubature points and weights for error calc (may be different from previous)
     DefaultCubatureFactory<double>  cubFactoryErr;
     int cubDegErr = 3;
     Teuchos::RCP<Cubature<double> > hexCubErr = cubFactoryErr.create(cellType, cubDegErr);
     int cubDimErr       = hexCubErr->getDimension();
     int numCubPointsErr = hexCubErr->getNumPoints();
     FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
     FieldContainer<double> cubWeightsErr(numCubPointsErr);
     hexCubErr->getCubature(cubPointsErr, cubWeightsErr);
     FieldContainer<double> physCubPointsE(numCells,numCubPointsErr, cubDimErr);

   // Containers for Jacobian
     FieldContainer<double> hexJacobianE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobInvE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobDetE(numCells, numCubPointsErr);
     FieldContainer<double> weightedMeasureE(numCells, numCubPointsErr);

 // Evaluate basis values and curls at cubature points
     FieldContainer<double> uhDVals(numFieldsD, numCubPointsErr, spaceDim);
     FieldContainer<double> uhDValsTrans(numCells,numFieldsD, numCubPointsErr, spaceDim);
     FieldContainer<double> uhDivs(numFieldsD, numCubPointsErr);
     FieldContainer<double> uhDivsTrans(numCells, numFieldsD, numCubPointsErr);
     hexHDivBasis.getValues(uhDVals, cubPointsErr, OPERATOR_VALUE);
     hexHDivBasis.getValues(uhDivs, cubPointsErr, OPERATOR_DIV);


   // Loop over elements
    for (int k=0; k<numElems; k++){

      double L2errElem = 0.0;
      double HDiverrElem = 0.0;
      double uExact1, uExact2, uExact3;
      double divuExact;

     // physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }
     // Face signs
      for (int j=0; j<numFacesPerElem; j++) {
         hexFaceSigns(0,j) = -1.0;
         for (int i=0; i<numNodesPerFace; i++) {
           int indf=i+1;
           if (indf > numNodesPerFace) indf=0;
           if (elemToNode(k,refFaceToNode(j,0))==faceToNode(elemToFace(k,j),i) &&
               elemToNode(k,refFaceToNode(j,1))==faceToNode(elemToFace(k,j),indf))
                hexFaceSigns(0,j) = 1.0;
           }
           if (!faceIsOwned[elemToFace(k,j)]){
              hexFaceSigns(0,j)=-1.0*hexFaceSigns(0,j);
           }
       }

    // compute cell Jacobians, their inverses and their determinants
       IntrepidCTools::setJacobian(hexJacobianE, cubPointsErr, hexNodes, cellType);
       IntrepidCTools::setJacobianInv(hexJacobInvE, hexJacobianE );
       IntrepidCTools::setJacobianDet(hexJacobDetE, hexJacobianE );

      // transform integration points to physical points
       IntrepidCTools::mapToPhysicalFrame(physCubPointsE, cubPointsErr, hexNodes, cellType);

      // transform basis values to physical coordinates
       IntrepidFSTools::HDIVtransformVALUE<double>(uhDValsTrans, hexJacobianE, hexJacobDetE, uhDVals);
       IntrepidFSTools::HDIVtransformDIV<double>(uhDivsTrans, hexJacobDetE, uhDivs);

      // compute weighted measure
       IntrepidFSTools::computeCellMeasure<double>(weightedMeasureE, hexJacobDetE, cubWeightsErr);

     // loop over cubature points
       for (int nPt = 0; nPt < numCubPointsErr; nPt++){

         // get exact solution and divs
          double x = physCubPointsE(0,nPt,0);
          double y = physCubPointsE(0,nPt,1);
          double z = physCubPointsE(0,nPt,2);
          evalu(uExact1, uExact2, uExact3, x, y, z);
          divuExact = evalDivu(x, y, z);

         // calculate approximate solution and divs
          double uApprox1 = 0.0;
          double uApprox2 = 0.0;
          double uApprox3 = 0.0;
          double divuApprox = 0.0;
          for (int i = 0; i < numFieldsD; i++){
             int rowIndex = globalFaceIds[elemToFace(k,i)];
#ifdef HAVE_MPI
             double uh1 = uCoeff.Values()[rowIndex];
#else
             double uh1 = xh.Values()[rowIndex];
#endif
             uApprox1 += uh1*uhDValsTrans(0,i,nPt,0)*hexFaceSigns(0,i);
             uApprox2 += uh1*uhDValsTrans(0,i,nPt,1)*hexFaceSigns(0,i);
             uApprox3 += uh1*uhDValsTrans(0,i,nPt,2)*hexFaceSigns(0,i);
             divuApprox += uh1*uhDivsTrans(0,i,nPt)*hexFaceSigns(0,i);
          }

         // evaluate the error at cubature points
          Linferr = max(Linferr, abs(uExact1 - uApprox1));
          Linferr = max(Linferr, abs(uExact2 - uApprox2));
          Linferr = max(Linferr, abs(uExact3 - uApprox3));
          L2errElem+=(uExact1 - uApprox1)*(uExact1 - uApprox1)*weightedMeasureE(0,nPt);
          L2errElem+=(uExact2 - uApprox2)*(uExact2 - uApprox2)*weightedMeasureE(0,nPt);
          L2errElem+=(uExact3 - uApprox3)*(uExact3 - uApprox3)*weightedMeasureE(0,nPt);
          HDiverrElem+=((divuExact - divuApprox)*(divuExact - divuApprox))
                     *weightedMeasureE(0,nPt);
        }

       L2err+=L2errElem;
       HDiverr+=HDiverrElem;
     }

#ifdef HAVE_MPI
   // sum over all processors
    Comm.SumAll(&L2err,&L2errTot,1);
    Comm.SumAll(&HDiverr,&HDiverrTot,1);
    Comm.MaxAll(&Linferr,&LinferrTot,1);

#else
    L2errTot = L2err;
    HDiverrTot = HDiverr;
    LinferrTot = Linferr;
#endif


  if (MyPID == 0) {
    std::cout << "\n" << "L2 Error:  " << sqrt(L2errTot) <<"\n";
    std::cout << "HDiv Error:  " << sqrt(HDiverrTot) <<"\n";
    std::cout << "LInf Error:  " << LinferrTot <<"\n\n";
  }


 // delete mesh
 Delete_Pamgen_Mesh();

 //clean up
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
   delete [] nodeCoordx;
   delete [] nodeCoordy;
   delete [] nodeCoordz;

   delete [] globalNodeIds;
   delete [] nodeIsOwned;
   delete [] globalEdgeIds;
   delete [] edgeIsOwned;
   delete [] globalFaceIds;
   delete [] faceIsOwned;
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


   return 0;
}
/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/


/*************************************************************************************/
/****************************** Multiply Ones ****************************************/
/*************************************************************************************/
/** \brief Multiplies Ax = y, where all non-zero entries of A are replaced with the value 1.0

    \param  A                [in]    matrix
    \param  x                [in]    vector
    \param  y                [in]    vector
 */
int Multiply_Ones(const Epetra_CrsMatrix &A,const Epetra_Vector &x,Epetra_Vector &y){
  if(!A.Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();
  const Epetra_Import* Importer_=A.Importer();
  const Epetra_Export* Exporter_=A.Exporter();
  Epetra_Vector *xcopy=0, *ImportVector_=0, *ExportVector_=0;

  if (&x==&y && Importer_==0 && Exporter_==0) {
    xcopy = new Epetra_Vector(x);
    xp = (double *) xcopy->Values();
  }
  else if (Importer_)
    ImportVector_ = new Epetra_Vector(Importer_->TargetMap());
  else if (Exporter_)
    ExportVector_ = new Epetra_Vector(Exporter_->SourceMap());
  

  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if(Importer_ != 0) {
    EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer_, Insert));
    xp = (double*) ImportVector_->Values();
    }
  
  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if(Exporter_ != 0)  yp = (double*) ExportVector_->Values();
  
  // Do actual computation
  for(int i = 0; i < A.NumMyRows(); i++) {
    int NumEntries,*RowIndices;
    A.Graph().ExtractMyRowView(i,NumEntries,RowIndices);
    double sum = 0.0;
    for(int j = 0; j < NumEntries; j++) 
      sum += xp[*RowIndices++];    
    yp[i] = sum;    
  }
  
  if(Exporter_ != 0) {
    y.PutScalar(0.0); // Make sure target is zero
    EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter_, Add)); // Fill y with Values from export vector
  }
  // Handle case of rangemap being a local replicated map
  if (!A.Graph().RangeMap().DistributedGlobal() && A.Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());

  delete xcopy;
  delete ImportVector_;
  delete ExportVector_;


  return(0);
}

/*************************************************************************************/
/****************************** Multiply Abs *****************************************/
/*************************************************************************************/
// Multiplies abs(A)x = y, where all non-zero entries of A are replaced with their absolute values value 
int Multiply_Abs(const Epetra_CrsMatrix &A,const Epetra_Vector &x,Epetra_Vector &y){
  if(!A.Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();
  const Epetra_Import* Importer_=A.Importer();
  const Epetra_Export* Exporter_=A.Exporter();
  Epetra_Vector *xcopy=0, *ImportVector_=0, *ExportVector_=0;

  if (&x==&y && Importer_==0 && Exporter_==0) {
    xcopy = new Epetra_Vector(x);
    xp = (double *) xcopy->Values();
  }
  else if (Importer_)
    ImportVector_ = new Epetra_Vector(Importer_->TargetMap());
  else if (Exporter_)
    ExportVector_ = new Epetra_Vector(Exporter_->SourceMap());
  

  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if(Importer_ != 0) {
    EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer_, Insert));
    xp = (double*) ImportVector_->Values();
    }
  
  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if(Exporter_ != 0)  yp = (double*) ExportVector_->Values();
  
  // Do actual computation
  for(int i = 0; i < A.NumMyRows(); i++) {
    int NumEntries,*RowIndices;
    double *RowValues;
    A.ExtractMyRowView(i,NumEntries,RowValues,RowIndices);
    double sum = 0.0;
    for(int j = 0; j < NumEntries; j++){
      double v=*RowValues++;
      sum += ABS(v) * xp[*RowIndices++];	
    }
    yp[i] = sum;    
  }
  
  if(Exporter_ != 0) {
    y.PutScalar(0.0); // Make sure target is zero
    EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter_, Add)); // Fill y with Values from export vector
  }
  // Handle case of rangemap being a local replicated map
  if (!A.Graph().RangeMap().DistributedGlobal() && A.Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());

  delete xcopy;
  delete ImportVector_;
  delete ExportVector_;


  return(0);
}

/*************************************************************************************/
/**************************** GET ML RESIDUAL ****************************************/
/*************************************************************************************/
/** \brief Compute ML solution residual

    \param  A                  [in]    discrete operator
    \param  lhs                [in]    solution vector
    \param  rhs                [in]    right hand side vector
    \param  Time               [in]    elapsed time for output
    \param  TotalErrorResidual [out]   error residual
    \param  TotalErrorExactSol [out]   error in xh (not an appropriate measure
                                                    for H(div) basis functions)
 */
void solution_test(string msg, const Epetra_Operator &A,const Epetra_MultiVector &lhs,const Epetra_MultiVector &rhs,const Epetra_MultiVector &xexact,Epetra_Time & Time, double & TotalErrorExactSol, double& TotalErrorResidual){
  // ==================================================== //  
  // compute difference between exact solution and ML one //
  // ==================================================== //  
  double d = 0.0, d_tot = 0.0;  
  for( int i=0 ; i<lhs.Map().NumMyElements() ; ++i )
    d += (lhs[0][i] - xexact[0][i]) * (lhs[0][i] - xexact[0][i]);
  
  A.Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs.Map());
  A.Apply(lhs, Ax);
  Ax.Update(1.0, rhs, -1.0);
  Ax.Norm2(&Norm);
  
  if (A.Comm().MyPID() == 0) {
    cout << msg << "......Using " << A.Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
//    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
}


/*************************************************************************************/
/*************************** ML PRECONDITIONER****************************************/
/*************************************************************************************/
void TestMultiLevelPreconditioner_DivLSFEM(char ProblemType[],
                                           Teuchos::ParameterList   & MLList,
                                           Epetra_CrsMatrix   & GradDiv,
                                           Epetra_CrsMatrix   & D0clean,
                                           Epetra_CrsMatrix   & D1clean,
                                           Epetra_CrsMatrix   & FaceNode,
                                           Epetra_CrsMatrix   & M1,
                                           Epetra_CrsMatrix   & M1inv,
                                           Epetra_CrsMatrix   & M2,
                                           Epetra_MultiVector & xh,
                                           Epetra_MultiVector & b,
                                           double & TotalErrorResidual,
                                           double & TotalErrorExactSol){
  
  /* Nuke M1 for D0, OAZ*/
  Epetra_CrsMatrix D1(D1clean);
  ML_Epetra::Apply_BCsToGradient(GradDiv,D1);  
  
  /* Get the BC faces*/
  int numBCfaces;  
  int* BCfaces=ML_Epetra::FindLocalDiricheltRowsFromOnesAndZeros(GradDiv,numBCfaces);  
  ML_Epetra::Apply_OAZToMatrix(BCfaces,numBCfaces,M2);

  if(!GradDiv.Comm().MyPID())
    cout<<"Total number of rows = "<<GradDiv.NumGlobalRows()<<endl;

  /* Build the (1,1) Block Operator */
  ML_Epetra::ML_RefMaxwell_11_Operator Operator11(GradDiv,D1,M1inv,M2);
  
  /* Build the AztecOO stuff */
  Epetra_MultiVector x(xh);
  x.PutScalar(0.0);
  
  Epetra_LinearProblem Problem(&Operator11,&x,&b); 
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  
  Epetra_Time Time(GradDiv.Comm());
    
  /* Build the aggregation guide matrix */
  Epetra_CrsMatrix *TMT_Agg_Matrix;
  ML_Epetra::ML_Epetra_PtAP(M1,D0clean,TMT_Agg_Matrix,false);
  
  /* Approximate the diagonal for EMFP: 2a^2 b guy */
  Epetra_Vector Diagonal(GradDiv.DomainMap());
  Epetra_Vector EdgeDiagonal(D1.DomainMap());
  Epetra_Vector FaceDiagonal(M2.DomainMap());
  Epetra_Vector DivDiagonal(GradDiv.DomainMap());

  M1inv.ExtractDiagonalCopy(EdgeDiagonal);
  M2.ExtractDiagonalCopy(FaceDiagonal);
  GradDiv.ExtractDiagonalCopy(DivDiagonal);

  Multiply_Ones(D1,EdgeDiagonal,Diagonal);

  for(int i=0;i<GradDiv.NumMyRows();i++) {
    Diagonal[i]=Diagonal[i]*(FaceDiagonal[i]*FaceDiagonal[i]) + DivDiagonal[i];
    if(ABS(Diagonal[i])<1e-12) Diagonal[i]=1.0;
  }

#ifdef DUMP_DATA
  // Dump matrices to disk
   EpetraExt::VectorToMatrixMarketFile("mag_est_diag.dat",Diagonal,0,0,false);
   EpetraExt::RowMatrixToMatlabFile("mag_tmt_matrix.dat",*TMT_Agg_Matrix);
#endif


  /* Build the EMFP Preconditioner */  
  ML_Epetra::FaceMatrixFreePreconditioner FMFP(Operator11,Diagonal,D1,FaceNode,*TMT_Agg_Matrix,BCfaces,numBCfaces,MLList);

  /* Solve! */
  AztecOO solver(Problem);  
  solver.SetPrecOperator(&FMFP);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-10);
  //  solver.Iterate(1, 1e-10);

  Epetra_MultiVector xexact(xh);
  xexact.PutScalar(0.0);
  
  // accuracy check
  string msg = ProblemType;
  solution_test(msg,Operator11,*lhs,*rhs,xexact,Time,TotalErrorExactSol,TotalErrorResidual);

#ifdef DUMP_DATA
   EpetraExt::MultiVectorToMatrixMarketFile("mag_lhs2.dat",*lhs,0,0,false);
#endif


  xh = *lhs;
  
  // Cleanup
  delete TMT_Agg_Matrix;
  delete [] BCfaces;
}

/**********************************************************************************/
/************ USER DEFINED FUNCTIONS FOR EXACT SOLUTION ***************************/
/**********************************************************************************/

// Calculates value of exact solution u
 int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z)
 {

/*
   // Exact Solution 1 - homogeneous boundary conditions, nonzero curl
    uExact0 = exp(y+z)*(x+1.0)*(x-1.0);
    uExact1 = exp(x+z)*(y+1.0)*(y-1.0);
    uExact2 = exp(x+y)*(z+1.0)*(z-1.0);
*/
  
/*
   // Exact Solution 2 - homogeneous boundary conditions, nonzero curl
    uExact0 = cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    uExact1 = cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    uExact2 = cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
*/

 /*
   // Exact Solution 3 - homogeneous boundary conditions, zero curl
    uExact0 = x*x-1.0;
    uExact1 = y*y-1.0;
    uExact2 = z*z-1.0;
 */  

   // Exact solution 4 - patch test with inhomogeneous boundary conditions,
   //                    zero curl, linear field should be recovered
    uExact0 = 1.0 + 2.0*x;
    uExact1 = 3.0 + 4.0*y;
    uExact2 = 5.0 + 6.0*z;

   return 0;
 }

// Calculates divergence of exact solution u
 double evalDivu(double & x, double & y, double & z)
 {
   
   // Exact Solution 1 - homogeneous boundary conditions, nonzero curl
   // double divu = 2.0*x*exp(y+z)+2.0*y*exp(x+z)+2.0*z*exp(x+y);

   // Exact Solution 2 - homogeneous boundary conditions, nonzero curl
   //double divu = 2.0*x*cos(M_PI*y)*cos(M_PI*z) + 2.0*y*cos(M_PI*x)*cos(M_PI*z)
   //               + 2.0*z*cos(M_PI*x)*cos(M_PI*y);
   
   // Exact Solution 3 - homogeneous boundary conditions, zero curl
   // double divu = 2.0*(x + y + z);

   // Exact solution 4 - patch test with inhomogeneous boundary conditions,
   //                    zero curl, linear field should be recovered
     double divu = 12.0;

   return divu;
 }


// Calculates curl of exact solution u
 int evalCurlu(double & curlu0, double & curlu1, double & curlu2, 
                double & x, double & y, double & z, double & mu)
 {
  
  /*
   // Exact Solution 1 - homogeneous boundary conditions, nonzero curl
    double duxdy = exp(y+z)*(x+1.0)*(x-1.0);
    double duxdz = exp(y+z)*(x+1.0)*(x-1.0);
    double duydx = exp(x+z)*(y+1.0)*(y-1.0);
    double duydz = exp(x+z)*(y+1.0)*(y-1.0);
    double duzdx = exp(x+y)*(z+1.0)*(z-1.0);
    double duzdy = exp(x+y)*(z+1.0)*(z-1.0);
  */
 

  /*
   // Exact Solution 2 - homogeneous boundary conditions, nonzero curl
    double duxdy = -M_PI*sin(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    double duxdz = -M_PI*sin(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0);
    double duydx = -M_PI*sin(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    double duydz = -M_PI*sin(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double duzdx = -M_PI*sin(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double duzdy = -M_PI*sin(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0);
  */

  /*
    curlu0 = (duzdy - duydz)/mu;
    curlu1 = (duxdz - duzdx)/mu;
    curlu2 = (duydx - duxdy)/mu;
  */
 
  /*
   // Exact Solution 3 - homogeneous boundary conditions, zero curl
    curlu0 = 0;
    curlu1 = 0;
    curlu2 = 0;
  */

   // Exact solution 4 - patch test with inhomogeneous boundary conditions,
   //                    zero curl, linear field should be recovered
    curlu0 = 0;
    curlu1 = 0;
    curlu2 = 0;
  
   return 0;
 }

// Calculates curl of the curl of exact solution u
 int evalCurlCurlu(double & curlCurlu0, double & curlCurlu1, double & curlCurlu2, 
                    double & x, double & y, double & z, double & mu)
{
   
 /*
   // Exact Solution 1 - homogeneous boundary conditions, nonzero curl
    double dcurlu0dy = exp(x+y)*(z+1.0)*(z-1.0) - 2.0*y*exp(x+z);
    double dcurlu0dz = 2.0*z*exp(x+y) - exp(x+z)*(y+1.0)*(y-1.0); 
    double dcurlu1dx = 2.0*x*exp(y+z) - exp(x+y)*(z+1.0)*(z-1.0); 
    double dcurlu1dz = exp(y+z)*(x+1.0)*(x-1.0) - 2.0*z*exp(x+y);
    double dcurlu2dx = exp(x+z)*(y+1.0)*(y-1.0) - 2.0*x*exp(y+z);
    double dcurlu2dy = 2.0*y*exp(x+z) - exp(y+z)*(x+1.0)*(x-1.0);
  */
                       

  /*
   // Exact Solution 2 - homogeneous boundary conditions, nonzero curl
    double dcurlu0dy = -M_PI*M_PI*cos(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0)
                           + 2.0*y*M_PI*sin(M_PI*z)*cos(M_PI*x);
    double dcurlu0dz = -2.0*z*M_PI*sin(M_PI*y)*cos(M_PI*x)
                          + M_PI*M_PI*cos(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double dcurlu1dx = -2.0*x*M_PI*sin(M_PI*z)*cos(M_PI*y)
                          + M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double dcurlu1dz = -M_PI*M_PI*cos(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0)
                           + 2.0*z*M_PI*sin(M_PI*x)*cos(M_PI*y);
    double dcurlu2dx = -M_PI*M_PI*cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0)
                           + 2.0*x*M_PI*sin(M_PI*y)*cos(M_PI*z);
    double dcurlu2dy = -2.0*y*M_PI*sin(M_PI*x)*cos(M_PI*z)
                          + M_PI*M_PI*cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
  */
                       
  /*
    curlCurlu0 = (dcurlu2dy - dcurlu1dz)/mu;
    curlCurlu1 = (dcurlu0dz - dcurlu2dx)/mu;
    curlCurlu2 = (dcurlu1dx - dcurlu0dy)/mu;
  */
 
 /*
   // Exact Solution 3 - homogeneous boundary conditions, zero curl
    curlCurlu0 = 0.0;
    curlCurlu1 = 0.0;
    curlCurlu2 = 0.0;
 */

   // Exact solution 4 - patch test with inhomogeneous boundary conditions,
   //                    zero curl, linear field should be recovered
    curlCurlu0 = 0.0;
    curlCurlu1 = 0.0;
    curlCurlu2 = 0.0;

    return 0;
}

