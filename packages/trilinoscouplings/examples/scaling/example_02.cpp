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

/** \file   example_01.cpp
    \brief  Example creation of mass and stiffness matrices for div-curl system on a hexadedral mesh using curl-conforming elements.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
 
    \remark Sample command line
    \code   ./example_01.exe 10 10 10 false 1.0 10.0 0.0 1.0 -1.0 1.0 -1.0 1.0 \endcode
*/

#undef DEBUG_PRINTING
// #define DEBUG_PRINTING

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_SerialComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

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
#include "../mesh_spec_lt/im_exodusII_l.h"
#include "../mesh_spec_lt/im_ne_nemesisI_l.h"

using namespace std;
using namespace Intrepid;
using namespace shards;

class topo_entity{
public:
  topo_entity(){
    local_id = -1;
    global_id = -1;
    owned = true;
  };
  void add_node(long long the_val,long long * global_nids){
    local_node_ids.push_back(the_val);
    sorted_local_node_ids.push_back(the_val);
    sorted_global_node_ids.push_back(global_nids[the_val-1]);
  }
  void sort(){
    sorted_local_node_ids.sort();
    sorted_global_node_ids.sort();
  }
  ~topo_entity(){};
  std::list <long long > local_node_ids;
  std::list <long long > sorted_local_node_ids;
  std::list <long long > sorted_global_node_ids;
  long long local_id;
  long long global_id;
  bool owned;
};

/*******************************************************************************/
inline bool compare_sorted_global_node_ids ( topo_entity* const x,  topo_entity* const y )
/*******************************************************************************/
{
  assert(x->sorted_global_node_ids.size() == y->sorted_global_node_ids.size());
  if(x->sorted_global_node_ids < y->sorted_global_node_ids)return true;    
  return false;
}

struct fecomp{
  bool operator () ( topo_entity* x,  topo_entity*  y )const
  {
    if(x->sorted_local_node_ids < y->sorted_local_node_ids)return true;    
    return false;
  }
};

void  Conform_Boundary_IDS_topo_entity(std::vector < std:: vector < topo_entity * > > & topo_entities,
					long long * proc_ids, 
					long long  rank);

void  Conform_Boundary_IDS(long long ** comm_entities,
			   long long * entity_counts,
			   long long * proc_ids, 
			   long long * data_array,
			   long long num_comm_pairs,
			   long long  rank);

/*******************************************************************************/
void calc_global_node_ids(long long * globalNodeIds,
			  long long numNodes,
			  long long num_node_comm_maps,
			  long long * node_cmap_node_cnts,
			  long long * node_comm_proc_ids,
			  long long * * comm_node_ids,
			  int rank)
/*******************************************************************************/
{
  for(long long i = 0; i < numNodes; i ++)globalNodeIds[i] = 1l;
  for(long long j = 0; j < num_node_comm_maps; j++) {
    for(long long k = 0; k < node_cmap_node_cnts[j] ; k ++){
      if(node_comm_proc_ids[j] < rank)globalNodeIds[comm_node_ids[j][k]-1] = -1;	
    }
  }
  long long num_unique_nodes = 0;
  for(long long  i = 0 ; i < numNodes; i ++)if(globalNodeIds[i] == 1l)num_unique_nodes ++;
  long long start_id = 0;
  MPI_Scan(&num_unique_nodes,&start_id,1,
	   MPI_LONG_LONG_INT,
	   MPI_SUM,
	   MPI_COMM_WORLD);
  start_id -= num_unique_nodes;

  int num_assigned = 0;
  for(long long  i = 0 ; i < numNodes; i ++)if(globalNodeIds[i] == 1l){
    globalNodeIds[i] = num_assigned + start_id;
    num_assigned ++;
  }

  //Conforms global nodal ids
  Conform_Boundary_IDS(comm_node_ids,
		       node_cmap_node_cnts,
		       node_comm_proc_ids, 
		       globalNodeIds,
		       num_node_comm_maps,
		       rank);

}


/*******************************************************************************/
void calc_global_ids(std::vector < topo_entity * > eof_vec,
		long long **comm_node_ids,
		long long * node_comm_proc_ids,
		long long * node_cmap_node_cnts,
		int num_node_comm_maps,
		int rank,
		std::string fname_string)
/*******************************************************************************/
{
  std::vector < std:: vector < topo_entity *> > topo_entities;
  int nncm = num_node_comm_maps;
  // make a vector of sets of comm nodes
  std::vector < std::set < long long> > comm_vector_set;
  for(int i = 0; i < nncm; i ++){
    std::vector < topo_entity * >v;
    topo_entities.push_back(v);
    std::set < long long > as;
    for(int j = 0; j < node_cmap_node_cnts[i];j ++){
      as.insert(comm_node_ids[i][j]);
    }
    comm_vector_set.push_back(as);
  }
/*run over all edges, faces*/

  for(unsigned tec = 0;tec != eof_vec.size();tec ++){
    topo_entity * teof = eof_vec[tec];

    for(int i = 0; i < nncm; i ++){
      bool found = true;
      std::list <long long > :: iterator lit;
      for( lit = teof->local_node_ids.begin();
	   lit != teof->local_node_ids.end() && found == true;
	   lit ++){
	if(comm_vector_set[i].find(*lit) == comm_vector_set[i].end())found = false;
      }
      //if all component nodes found then add face,edge to comm lists
      if(found){
	topo_entities[i].push_back(teof);
	if(node_comm_proc_ids[i] < rank)teof->owned = false;//not owned if found on lower processor
      }
      else{
      }
    }
  }

  //need to sort the edges_or_face vectors by their sorted global node ids
  for(unsigned i = 0; i < topo_entities.size(); i ++){
    if(!topo_entities[i].empty()){
      std::sort(topo_entities[i].begin(),topo_entities[i].end(),compare_sorted_global_node_ids);
    }
  }
#ifdef DEBUG_PRINTING
 std::stringstream aname;
  aname << fname_string;
  aname << rank;
  aname << ".txt";
  ofstream fout(aname.str().c_str());
#endif

  //need to sort the edges_or_face vectors by their sorted global node ids
  for(unsigned i = 0; i < topo_entities.size(); i ++){
    if(!topo_entities[i].empty()){
      std::sort(topo_entities[i].begin(),topo_entities[i].end(),compare_sorted_global_node_ids);
    }
#ifdef DEBUG_PRINTING
    fout << " from proc rank " << rank 
	 << " to proc rank " << node_comm_proc_ids[i] 
	 << " has " << topo_entities[i].size() 
	 << " entries " << std::endl;

    if(!topo_entities[i].empty()){
      for(unsigned j = 0; j < topo_entities[i].size();j ++){
	topo_entity * eof = topo_entities[i][j];
	for(std::list < long long > :: iterator klit = eof->sorted_global_node_ids.begin();
	    klit != eof->sorted_global_node_ids.end(); klit ++){
	  fout << (*klit) << " " ;
	}
	fout << endl;
      }
    }
#endif
  }

  //count the number of entities owned;
  long long owned_entities = 0;
  for(unsigned ict = 0; ict < eof_vec.size(); ict ++){
    if(eof_vec[ict]->owned)owned_entities ++;
  }
#ifdef DEBUG_PRINTING
  fout << " proc " << rank << " owns " << owned_entities << " edges " << std::endl;
#endif

  long long start_id = 0;
  MPI_Scan(&owned_entities,&start_id,1,
	   MPI_LONG_LONG_INT,
	   MPI_SUM,
	   MPI_COMM_WORLD);
  start_id -= owned_entities;
#ifdef DEBUG_PRINTING
  fout << " proc " << rank << " start_id " << start_id << std::endl;
#endif
  //DMH
  long long num_assigned = 0;
  for(unsigned ict = 0; ict < eof_vec.size(); ict ++){
    if(eof_vec[ict]->owned){
      eof_vec[ict]->global_id = num_assigned + start_id;
      start_id ++;
    }
  }

  Conform_Boundary_IDS_topo_entity(topo_entities,
				    node_comm_proc_ids, 
				    rank);

#ifdef DEBUG_PRINTING
  for(unsigned ict = 0; ict < eof_vec.size(); ict ++){
    fout << "on proc " << rank << " entity " << ict << " has id " << eof_vec[ict]->global_id << std::endl;
  }

  fout.close();
#endif

}

// Functions to evaluate exact solution and derivatives
int evalu(double & uExact0, 
	  double & uExact1, 
	  double & uExact2, 
	  double & x, 
	  double & y, 
	  double & z);

double evalDivu(double & x, 
		double & y, 
		double & z);

int evalCurlu(double & curlu0, 
	      double & curlu1, 
	      double & curlu2, 
	      double & x, 
	      double & y, 
	      double & z);

int evalGradDivu(double & gradDivu0, 
		 double & gradDivu1, 
		 double & gradDivu2, 
		 double & x, 
		 double & y, 
		 double & z);




int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int error = 0;
  int rank=mpiSession.getRank();
  int numProcs=mpiSession.getNProc();
   //Check number of arguments
    TEST_FOR_EXCEPTION( ( argc < 13 ),
                      std::invalid_argument,
                      ">>> ERROR (example_01): Invalid number of arguments. See code listing for requirements.");
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 12)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|          Example: Div-Curl System on Hexahedral Mesh                        |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";


  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;
  
  std::set < topo_entity * , fecomp > edge_set;
  std::set < topo_entity * , fecomp > face_set;

  std::vector < topo_entity * > node_vector;
  std::vector < topo_entity * > edge_vector;
  std::vector < topo_entity * > face_vector;

  std::vector < int > edge_comm_procs;

  int dim = 3;
  int spaceDim = 3;


// ************************************ GET INPUTS **************************************

  /* In the implementation for discontinuous material properties only the boundaries for
     region 1, associated with mu1, are input. The remainder of the grid is assumed to use mu2.
     Note that the material properties are assigned using the undeformed grid. */

    int NX            = atoi(argv[1]);  // num intervals in x direction (assumed box domain, -1,1)
    int NY            = atoi(argv[2]);  // num intervals in y direction (assumed box domain, -1,1)
    int NZ            = atoi(argv[3]);  // num intervals in z direction (assumed box domain, -1,1)
    int randomMesh    = atoi(argv[4]);  // 1 if mesh randomizer is to be used 0 if not
    double mu1        = atof(argv[5]);  // material property value for region 1
    double mu2        = atof(argv[6]);  // material property value for region 2
    double mu1LeftX   = atof(argv[7]);  // left X boundary for region 1
    double mu1RightX  = atof(argv[8]);  // right X boundary for region 1
    double mu1LeftY   = atof(argv[9]);  // left Y boundary for region 1
    double mu1RightY  = atof(argv[10]); // right Y boundary for region 1
    double mu1LeftZ   = atof(argv[11]); // left Z boundary for region 1
    double mu1RightZ  = atof(argv[12]); // right Z boundary for region 1

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8;
    if(dim == 3){
    hex_8 = CellTopology(shards::getCellTopologyData<Hexahedron<8> >() );
    }
    else{//only handling 3D currently

      TEST_FOR_EXCEPTION( ( dim != 3 ),
			  std::invalid_argument,
			  ">>> ERROR (example_02): only 3D supported in this executable");
      
    }

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int numEdgesPerElem = hex_8.getEdgeCount();
    int numFacesPerElem = hex_8.getSideCount();
    int numNodesPerFace = 4;
    int numNodesPerEdge = 2;

   // Build reference element edge to node map
    FieldContainer<int> refEdgeToNode(numEdgesPerElem,numNodesPerEdge);
    for (int i=0; i<numEdgesPerElem; i++){
      for(int j = 0; j < numNodesPerEdge;j++){
        refEdgeToNode(i,j)=hex_8.getNodeMap(1, i, j);
      }
    }

    FieldContainer<int> refFaceToNode(numFacesPerElem,numNodesPerFace);
    for (int i=0; i<numFacesPerElem; i++){
      for(int j = 0; j < numNodesPerFace;j++){
	refFaceToNode(i,j)=hex_8.getNodeMap(2, i, j);
      }
    }

// *********************************** GENERATE MESH ************************************

    std::cout << "Generating mesh ... \n\n";

    std::cout << "    NX" << "   NY" << "   NZ\n";
    std::cout << std::setw(5) << NX <<
                 std::setw(5) << NY <<
                 std::setw(5) << NZ << "\n\n";

   // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    double leftZ = -1.0, rightZ = 1.0;

   // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

   // Create Pamgen input file
    stringstream ss;
    ss.clear();
    ss << "mesh \n";
    ss << "  rectilinear \n"; 
    ss << "     nx = " << NX << "\n";
    ss << "     ny = " << NY << "\n"; 
    ss << "     nz = " << NZ << "\n"; 
    ss << "     bx = 1\n";
    ss << "     by = 1\n"; 
    ss << "     bz = 1\n"; 
    ss << "     gmin = " << leftX << " " << leftY << " " << leftZ << "\n";
    ss << "     gmax = " << rightX << " " << rightY << " " << rightZ << "\n";
    ss << "  end \n";
    ss << "  set assign \n";
    ss << "     sideset, ilo, 1\n"; 
    ss << "     sideset, jlo, 2\n"; 
    ss << "     sideset, klo, 3\n"; 
    ss << "     sideset, ihi, 4\n"; 
    ss << "     sideset, jhi, 5\n"; 
    ss << "     sideset, khi, 6\n"; 
    ss << "  end \n";
    ss << "end \n";

    string meshInput = ss.str();
    std::cout << meshInput <<"\n";

   // Generate mesh with Pamgen

    long int maxInt = 9223372036854775807LL;
    Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
    
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
    long long * block_ids = new long long [numElemBlk];
    error += im_ex_get_elem_blk_ids_l(id, block_ids);

    int edgePerElem = 4;
    if(dim == 3)edgePerElem = 12;
    int facePerElem = 0;
    if(dim == 3)facePerElem = 6;
    FieldContainer<int> elemToEdge(numElems,edgePerElem);
    FieldContainer<int> elemToFace(numElems,facePerElem);



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
    delete [] nodeCoordz;

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

    calc_global_node_ids(globalNodeIds,
			 numNodes,
			 num_node_comm_maps,
			 node_cmap_node_cnts,
			 node_comm_proc_ids,
			 comm_node_ids,
			 rank);    

    //create edges and calculate edge ids
    /*connectivity*/
    int elct = 0;
    for(long long b = 0; b < numElemBlk; b++){
      if(nodes_per_element[b] == 4){
      }
      else if (nodes_per_element[b] == 8){
	//loop over all elements and push their edges onto a set if they are not there already
	for(long long el = 0; el < elements[b]; el++){
	  std::set< topo_entity * > ::iterator fit;
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
    
    FieldContainer<int> edgeToNode(edge_vector.size(), numNodesPerEdge);
    for(unsigned ect = 0; ect != edge_vector.size(); ect++){
      std::list<long long>::iterator elit;
      int nct = 0;
      for(elit  = edge_vector[ect]->local_node_ids.begin();
	  elit != edge_vector[ect]->local_node_ids.end();
	  elit ++){
	edgeToNode(ect,nct) = *elit;
	nct++;
      }
    }
    int numEdges = edge_vector.size();
    int numFaces = face_vector.size();
    std::cout << " Number of Elements: " << numElems << " \n";
    std::cout << "    Number of Nodes: " << numNodes << " \n";
    std::cout << "    Number of Edges: " << numEdges << " \n";
    std::cout << "    Number of Faces: " << numFaces << " \n\n";
   
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
  

    MPI_Finalize();
    exit(0);

   // Get boundary (side set) information
   // Side set 1 - left,  Side set 2 - front, Side set 3 - bottom,
   // Side set 4 - right, Side set 5 - back,  Side set 6 - top
    long long * sideSetIds = new long long [numSideSets];
    FieldContainer<int> numElemsOnBoundary(numSideSets);
    long long numSidesinSet;
    long long numDFinSet;
    long long maxNumSidesinSet = 0;
    im_ex_get_side_set_ids_l(id,sideSetIds);
    for (int i=0; i<numSideSets; i++) {
        im_ex_get_side_set_param_l(id,sideSetIds[i],&numSidesinSet,&numDFinSet);
        numElemsOnBoundary(i)=numSidesinSet;
        if (numSidesinSet > maxNumSidesinSet)
           maxNumSidesinSet = numSidesinSet; 
     }
   // Container for global element numbers of boundary elements for side sets 1-6
    FieldContainer<int> elemsOnBoundary(numSideSets,maxNumSidesinSet);
   // Container for local face numbers corresponding to global elements in previous array
    FieldContainer<int> facesOnBoundary(numSideSets,maxNumSidesinSet);
    for (int i=0; i<numSideSets; i++) {
        numSidesinSet=numElemsOnBoundary(i);
        long long * sideSetElemList = new long long [numSidesinSet];
        long long * sideSetSideList = new long long [numSidesinSet];
        im_ex_get_side_set_l(id,sideSetIds[i],sideSetElemList,sideSetSideList);
        for (int j=0; j<numSidesinSet; j++) {
          elemsOnBoundary(i,j)=sideSetElemList[j] - 1;
          facesOnBoundary(i,j)=sideSetSideList[j] - 1;
        }
        delete [] sideSetElemList;
        delete [] sideSetSideList;
     }
    delete [] sideSetIds;

   // Print mesh information  

  
 
   // Output element to edge connectivity
    ofstream fout("elem2edge.dat");
   
    fout.close();


   // Set material properties using undeformed grid assuming each element has only one value of mu
    FieldContainer<double> muVal(numElems);
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          int ielem = i + j * NX + k * NX * NY;
          double midElemX = nodeCoord(elemToNode(ielem,0),0) + hx/2.0;
          double midElemY = nodeCoord(elemToNode(ielem,0),1) + hy/2.0;
          double midElemZ = nodeCoord(elemToNode(ielem,0),2) + hz/2.0;
          if ( (midElemX > mu1LeftX) && (midElemY > mu1LeftY) && (midElemZ > mu1LeftZ) &&
               (midElemX <= mu1RightX) && (midElemY <= mu1RightY) && (midElemZ <= mu1RightZ) ){
             muVal(ielem) = mu1;
          }
           else {
             muVal(ielem) = mu2;
          }
        }
      }
    }

   // Perturb mesh coordinates (only interior nodes)
    if (randomMesh){
      for (int k=1; k<NZ; k++) {
        for (int j=1; j<NY; j++) {
          for (int i=1; i<NX; i++) {
            int inode = i + j * (NX + 1) + k * (NX + 1) * (NY + 1);
           // random numbers between -1.0 and 1.0
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double rz = 2.0 * (double)rand()/RAND_MAX - 1.0; 
           // limit variation to 1/4 edge length
            nodeCoord(inode,0) = nodeCoord(inode,0) + 0.125 * hx * rx;
            nodeCoord(inode,1) = nodeCoord(inode,1) + 0.125 * hy * ry;
            nodeCoord(inode,2) = nodeCoord(inode,2) + 0.125 * hz * rz;
          }
        }
      }
    }




// **************************** INCIDENCE MATRIX **************************************

   // Node to edge incidence matrix
    std::cout << "Building incidence matrix ... \n\n";

    Epetra_SerialComm Comm;
    Epetra_Map globalMapC(numEdges, 0, Comm);
    Epetra_Map globalMapG(numNodes, 0, Comm);
    Epetra_FECrsMatrix DGrad(Copy, globalMapC, globalMapG, 2);

    double vals[2];
    vals[0]=-1; vals[1]=1;
    for (int j=0; j<numEdges; j++){
        int rowNum = j;
        int colNum[2];
        colNum[0] = edgeToNode(j,0);
        colNum[1] = edgeToNode(j,1);
        DGrad.InsertGlobalValues(1, &rowNum, 2, colNum, vals);
    }


// ************************************ CUBATURE **************************************

   // Get numerical integration points and weights
    std::cout << "Getting cubature ... \n\n";

    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);


// ************************************** BASIS ***************************************

   // Define basis 
    std::cout << "Getting basis ... \n\n";
    Basis_HCURL_HEX_I1_FEM<double, FieldContainer<double> > hexHCurlBasis;
    Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;

    int numFieldsC = hexHCurlBasis.getCardinality();
    int numFieldsG = hexHGradBasis.getCardinality();

  // Evaluate basis at cubature points
     FieldContainer<double> hexGVals(numFieldsG, numCubPoints); 
     FieldContainer<double> hexCVals(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> hexCurls(numFieldsC, numCubPoints, spaceDim); 

     hexHCurlBasis.getValues(hexCVals, cubPoints, OPERATOR_VALUE);
     hexHCurlBasis.getValues(hexCurls, cubPoints, OPERATOR_CURL);
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL MASS and STIFFNESS MATRICES *************


    std::cout << "Building mass and stiffness matrices ... \n\n";

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    typedef ArrayTools art;
    int numCells = 1; 

   // Containers for nodes and edge signs 
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
    FieldContainer<double> hexEdgeSigns(numCells, numFieldsC);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HGRAD mass matrix
    FieldContainer<double> massMatrixG(numCells, numFieldsG, numFieldsG);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> weightedMeasureMuInv(numCells, numCubPoints);
    FieldContainer<double> hexGValsTransformed(numCells, numFieldsG, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
   // Containers for element HCURL mass matrix
    FieldContainer<double> massMatrixC(numCells, numFieldsC, numFieldsC);
    FieldContainer<double> hexCValsTransformed(numCells, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> hexCValsTransformedWeighted(numCells, numFieldsC, numCubPoints, spaceDim);
   // Containers for element HCURL stiffness matrix
    FieldContainer<double> stiffMatrixC(numCells, numFieldsC, numFieldsC);
    FieldContainer<double> weightedMeasureMu(numCells, numCubPoints);    
    FieldContainer<double> hexCurlsTransformed(numCells, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> hexCurlsTransformedWeighted(numCells, numFieldsC, numCubPoints, spaceDim);
   // Containers for right hand side vectors
    FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
    FieldContainer<double> rhsDatah(numCells, numCubPoints, cubDim);
    FieldContainer<double> gC(numCells, numFieldsC);
    FieldContainer<double> hC(numCells, numFieldsC);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    
   // Global arrays in Epetra format
    Epetra_FECrsMatrix MassG(Copy, globalMapG, numFieldsG);
    Epetra_FECrsMatrix MassC(Copy, globalMapC, numFieldsC);
    Epetra_FECrsMatrix StiffC(Copy, globalMapC, numFieldsC);
    Epetra_FEVector rhsC(globalMapC);

    ofstream fSignsout("edgeSigns.dat");

 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

     // Edge signs
      for (int j=0; j<numEdgesPerElem; j++) {
          if (elemToNode(k,refEdgeToNode(j,0))==edgeToNode(elemToEdge(k,j),0) &&
              elemToNode(k,refEdgeToNode(j,1))==edgeToNode(elemToEdge(k,j),1))
              hexEdgeSigns(0,j) = 1.0;
          else 
              hexEdgeSigns(0,j) = -1.0;

         fSignsout << hexEdgeSigns(0,j) << "  ";
       }
       fSignsout << "\n";

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );

// ************************** Compute element HGrad mass matrices *******************************
  
     // transform to physical coordinates 
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);
      
     // compute weighted measure
      fst::computeMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

      // combine mu value with weighted measure
      for (int nC = 0; nC < numCells; nC++){
        for (int nPt = 0; nPt < numCubPoints; nPt++){
          weightedMeasureMuInv(nC,nPt) = weightedMeasure(nC,nPt) / muVal(k);
        }
      }
      
     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasureMuInv, hexGValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixG,
                             hexGValsTransformed, hexGValsTransformedWeighted, COMP_CPP);

      // assemble into global matrix
      int err = 0;
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = elemToNode(k,row);
            int colIndex = elemToNode(k,col);
            double val = massMatrixG(0,row,col);
     //       err = MassG.SumIntoGlobalValues(1, &rowIndex, 1, &colIndex, &val);
     //       if (err > 0) {
                MassG.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
    //        }
         }
      }

// ************************** Compute element HCurl mass matrices *******************************

     // transform to physical coordinates 
      fst::HCURLtransformVALUE<double>(hexCValsTransformed, hexJacobInv, hexEdgeSigns,
                                   hexCVals);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexCValsTransformedWeighted,
                                   weightedMeasure, hexCValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixC,
                             hexCValsTransformed, hexCValsTransformedWeighted,
                             COMP_CPP);

     // assemble into global matrix
      err = 0;
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = massMatrixC(0,row,col);
            err = MassC.SumIntoGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            if (err > 0) {
                MassC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            }
         }
      }

// ************************ Compute element HCurl stiffness matrices *****************************

      // transform to physical coordinates 
      fst::HCURLtransformCURL<double>(hexCurlsTransformed, hexJacobian, hexJacobDet, 
                                   hexEdgeSigns, hexCurls);

      // combine mu value with weighted measure
      for (int nC = 0; nC < numCells; nC++){
        for (int nPt = 0; nPt < numCubPoints; nPt++){
          weightedMeasureMu(nC,nPt) = weightedMeasure(nC,nPt) / muVal(k);
         }
      }

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexCurlsTransformedWeighted,
                                   weightedMeasureMu, hexCurlsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(stiffMatrixC,
                             hexCurlsTransformed, hexCurlsTransformedWeighted,
                             COMP_CPP);

     // assemble into global matrix
      err = 0;
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = stiffMatrixC(0,row,col);
            err = StiffC.SumIntoGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            if (err > 0) {
                StiffC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            }
         }
      }

// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side functions at physical points
       FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
       FieldContainer<double> rhsDatah(numCells, numCubPoints, cubDim);
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          double du1, du2, du3;

          evalCurlu(du1, du2, du3, x, y, z);
          rhsDatag(0,nPt,0) = du1;
          rhsDatag(0,nPt,1) = du2;
          rhsDatag(0,nPt,2) = du3;
         
          evalGradDivu(du1, du2, du3,  x, y, z);
          rhsDatah(0,nPt,0) = du1;
          rhsDatah(0,nPt,1) = du2;
          rhsDatah(0,nPt,2) = du3;
       }

     // integrate (g,curl w) term
      fst::integrate<double>(gC, rhsDatag, hexCurlsTransformedWeighted,
                             COMP_CPP);

     // integrate (h,div w) term
      fst::integrate<double>(hC, rhsDatah, hexCValsTransformedWeighted,
                             COMP_CPP);

    // assemble into global vector
     for (int row = 0; row < numFieldsC; row++){
           int rowIndex = elemToEdge(k,row);
           double val = gC(0,row)-hC(0,row);
           err = rhsC.SumIntoGlobalValues(1, &rowIndex, &val);
     }
 
     
 } // *** end element loop ***

  // Assemble over multiple processors, if necessary
   DGrad.GlobalAssemble();  DGrad.FillComplete();    
   MassG.GlobalAssemble();  MassG.FillComplete();
   MassC.GlobalAssemble();  MassC.FillComplete();
   StiffC.GlobalAssemble(); StiffC.FillComplete();
   rhsC.GlobalAssemble();
   
  // Dump matrices to disk
   EpetraExt::RowMatrixToMatlabFile("mag_m0_matrix.dat",MassG);
   EpetraExt::RowMatrixToMatlabFile("mag_m1_matrix.dat",MassC);
   EpetraExt::RowMatrixToMatlabFile("mag_k1_matrix.dat",StiffC);
   EpetraExt::RowMatrixToMatlabFile("mag_t_matrix.dat",DGrad);
   EpetraExt::MultiVectorToMatlabFile("rhs1_vector.dat",rhsC);

   fSignsout.close();

 // delete mesh
 Delete_Pamgen_Mesh();

 // reset format state of std::cout
 std::cout.copyfmt(oldFormatState);
 
 return 0;
}

// Calculates value of exact solution u
 int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z)
 {
    uExact0 = cos(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0);
    uExact1 = cos(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0);
    uExact2 = cos(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);

   return 0;
 }

// Calculates divergence of exact solution u
 double evalDivu(double & x, double & y, double & z)
 {
   double divu = -M_PI*sin(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0)
                 -M_PI*sin(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0)
                 -M_PI*sin(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);
   return divu;
 }


// Calculates curl of exact solution u
 int evalCurlu(double & curlu0, double & curlu1, double & curlu2, double & x, double & y, double & z)
 {
   double duxdy = cos(M_PI*x)*exp(y*z)*(z+1.0)*(z-1.0)*(z*(y+1.0)*(y-1.0) + 2.0*y);
   double duxdz = cos(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(y*(z+1.0)*(z-1.0) + 2.0*z);
   double duydx = cos(M_PI*y)*exp(x*z)*(z+1.0)*(z-1.0)*(z*(x+1.0)*(x-1.0) + 2.0*x);
   double duydz = cos(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(x*(z+1.0)*(z-1.0) + 2.0*z);
   double duzdx = cos(M_PI*z)*exp(x*y)*(y+1.0)*(y-1.0)*(y*(x+1.0)*(x-1.0) + 2.0*x);
   double duzdy = cos(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(x*(y+1.0)*(y-1.0) + 2.0*y);

   curlu0 = duzdy - duydz;
   curlu1 = duxdz - duzdx;
   curlu2 = duydx - duxdy;

   return 0;
 }

// Calculates gradient of the divergence of exact solution u
 int evalGradDivu(double & gradDivu0, double & gradDivu1, double & gradDivu2, double & x, double & y, double & z)
{
    gradDivu0 = -M_PI*M_PI*cos(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0)
                  -M_PI*sin(M_PI*y)*exp(x*z)*(z+1.0)*(z-1.0)*(z*(x+1.0)*(x-1.0)+2.0*x)
                  -M_PI*sin(M_PI*z)*exp(x*y)*(y+1.0)*(y-1.0)*(y*(x+1.0)*(x-1.0)+2.0*x);
    gradDivu1 = -M_PI*sin(M_PI*x)*exp(y*z)*(z+1.0)*(z-1.0)*(z*(y+1.0)*(y-1.0)+2.0*y)
                  -M_PI*M_PI*cos(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0)
                  -M_PI*sin(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(x*(y+1.0)*(y-1.0)+2.0*y);
    gradDivu2 = -M_PI*sin(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(y*(z+1.0)*(z-1.0)+2.0*z)
                  -M_PI*sin(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(x*(z+1.0)*(z-1.0)+2.0*z)
                  -M_PI*M_PI*cos(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);
    return 0;
}



/*****************************************************************************/
void  Conform_Boundary_IDS(long long ** comm_entities,
			   long long * entity_counts,
			   long long * proc_ids, 
			   long long * data_array,
			   long long num_comm_pairs,
			   long long  rank)
/*****************************************************************************/
{
  //relies on data array having a -1 for unassigned values
  unsigned nncm = num_comm_pairs;

  //Will load an ownership flag along with the data to allow the conform
  MPI_Request * req = new MPI_Request[nncm];
  
  long long ** send_buffer = new long long * [nncm];
  long long ** receive_buffer = new long long * [nncm];
  for(unsigned i = 0; i < nncm; i ++)send_buffer[i]    = new long long [entity_counts[i]];
  for(unsigned i = 0; i < nncm; i ++)receive_buffer[i] = new long long [entity_counts[i]];

  // load up send buffer
  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < entity_counts[i];j++){
      send_buffer[i][j] = data_array[comm_entities[i][j]-1];
    }
  }

  //communicate

  for(unsigned i = 0; i < nncm ;i ++){
    int size = entity_counts[i];
    int proc = proc_ids[i];
    MPI_Irecv(receive_buffer[i],size, MPI_LONG_LONG_INT, proc, 1, MPI_COMM_WORLD, req + i);
  }

  for(unsigned i = 0; i < nncm ;i ++){
    int size = entity_counts[i];
    int proc = proc_ids[i];
    MPI_Send(send_buffer[i], size, MPI_LONG_LONG_INT, proc, 1,MPI_COMM_WORLD);
  }

  for(unsigned i = 0; i < nncm ;i ++){
    MPI_Status stat;
    MPI_Wait(req + i, &stat);
  }

  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < entity_counts[i];j++){
      if(receive_buffer[i][j] >= 0)data_array[comm_entities[i][j]-1] = receive_buffer[i][j];
    }
  }
  
  for(unsigned i = 0; i < nncm; i ++){
    if(send_buffer[i])   delete [] send_buffer[i];
    if(receive_buffer[i])delete [] receive_buffer[i];
  }
  delete [] send_buffer;
  delete [] receive_buffer;
}


/*****************************************************************************/
void  Conform_Boundary_IDS_topo_entity(std::vector < std:: vector < topo_entity * > > & topo_entities,
				       long long * proc_ids, 
				       long long  rank)
/*****************************************************************************/
{
  //relies on data array having a -1 for unassigned values
  unsigned nncm = topo_entities.size();

  //Will load an ownership flag along with the data to allow the conform
  MPI_Request * req = new MPI_Request[nncm];
  
  long long ** send_buffer = new long long * [nncm];
  long long ** receive_buffer = new long long * [nncm];
  for(unsigned i = 0; i < nncm; i ++){
    if(topo_entities[i].size() > 0){
      send_buffer[i]    = new long long [topo_entities[i].size()];
      receive_buffer[i] = new long long [topo_entities[i].size()];
    }
    else{
      send_buffer[i] = NULL;
      receive_buffer[i] = NULL;
    }
  }

  // load up send buffer
  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < topo_entities[i].size();j++){
      send_buffer[i][j] = topo_entities[i][j]->global_id;
    }
  }

  //communicate

  for(unsigned i = 0; i < nncm ;i ++){
    int size = topo_entities[i].size();
    if(size > 0){
      int proc = proc_ids[i];
      MPI_Irecv(receive_buffer[i],size, MPI_LONG_LONG_INT, proc, 1, MPI_COMM_WORLD, req + i);
    }
  }

  for(unsigned i = 0; i < nncm ;i ++){
    int size = topo_entities[i].size();
    if(size > 0){
      int proc = proc_ids[i];
      MPI_Send(send_buffer[i], size, MPI_LONG_LONG_INT, proc, 1,MPI_COMM_WORLD);
    }
  }

  for(unsigned i = 0; i < nncm ;i ++){
    MPI_Status stat;
    int size = topo_entities[i].size();
    if(size > 0){
      MPI_Wait(req + i, &stat);
    }
  }

  for(unsigned i = 0; i < nncm; i ++){
    for(unsigned j = 0; j < topo_entities[i].size();j++){
      if(receive_buffer[i][j] >= 0)topo_entities[i][j]->global_id = receive_buffer[i][j];
    }
  }

  for(unsigned i = 0; i < nncm; i ++){
    if(send_buffer[i])   delete [] send_buffer[i];
    if(receive_buffer[i])delete [] receive_buffer[i];
  }
  delete [] send_buffer;
  delete [] receive_buffer;
}
