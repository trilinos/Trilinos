// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Basic testing of Zoltan2::APFMeshAdapter

#include <Zoltan2_APFMeshAdapter.hpp>

#ifdef HAVE_ZOLTAN2_PARMA
#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <gmi_mesh.h>
#endif

// Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace std;
using Teuchos::RCP;

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
//Tpetra typedefs
typedef Tpetra::DefaultPlatform::DefaultPlatformType            Platform;
typedef Tpetra::MultiVector<double, int, int>     tMVector_t;

//Topology type helper function
enum Zoltan2::EntityTopologyType topologyAPFtoZ2(enum apf::Mesh::Type ttype) {
  if (ttype==apf::Mesh::VERTEX)
    return Zoltan2::POINT;
  else if (ttype==apf::Mesh::EDGE)
    return Zoltan2::LINE_SEGMENT;
  else if (ttype==apf::Mesh::TRIANGLE)
    return Zoltan2::TRIANGLE;
  else if (ttype==apf::Mesh::QUAD)
    return Zoltan2::QUADRILATERAL;
  else if (ttype==apf::Mesh::TET)
    return Zoltan2::TETRAHEDRON;
  else if (ttype==apf::Mesh::HEX)
    return Zoltan2::HEXAHEDRON;
  else if (ttype==apf::Mesh::PRISM)
    return Zoltan2::PRISM;
  else if (ttype==apf::Mesh::PYRAMID)
    return Zoltan2::PYRAMID;
  else 
    throw "No such APF topology type";
    
}

/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int narg, char *arg[]) {

  Teuchos::GlobalMPISession mpiSession(&narg, &arg,0);
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > CommT = platform.getComm();

#ifdef HAVE_ZOLTAN2_PARMA
  //Open up PCU for communication required by all APF operations
  PCU_Comm_Init();
  
  //Open up the mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("pumiTri14/plate.dmg","pumiTri14/2/");
  apf::verify(m);

  int dim = m->getDimension();

  //Contruct the MeshAdapter
  typedef Zoltan2::APFMeshAdapter<apf::Mesh2*> Adapter;
  typedef Adapter::lno_t lno_t;
  typedef Adapter::zgid_t zgid_t;
  typedef Adapter::scalar_t scalar_t;

  std::string pri = "face";
  std::string adj = "edge";
  if (dim==3) {
    adj=pri;
    pri="region";
  }
  
  Zoltan2::APFMeshAdapter<apf::Mesh2*> ia(*CommT,m,pri,adj,true);

  Zoltan2::MeshEntityType ents[4] = {Zoltan2::MESH_VERTEX,
                                     Zoltan2::MESH_EDGE,
                                     Zoltan2::MESH_FACE,
                                     Zoltan2::MESH_REGION};
  int* numberGlobal = new int[dim];
  //Check the local number of each entity
  for (int i=0;i<=dim;i++) {
    if (ia.getLocalNumOf(ents[i])!=m->count(i)) {      
      std::cerr<<"Local number of entities does not match in dimension "<<i<<"\n";
      return 1;
    }
    numberGlobal[i] = countOwned(m,i);
  }

  //Check the global number of each entity
  PCU_Add_Ints(numberGlobal,dim+1);
  for (int i=0;i<=dim;i++)
    if (ia.getGlobalNumOf(ents[i])!=static_cast<unsigned int>(numberGlobal[i])) {
      std::cerr<<"Global number of entities does not match in dimension "<<i<<"\n";
      return 1;
    }
  delete [] numberGlobal;

  //Check the coordinate dimension
  apf::GlobalNumbering** gnums = new apf::GlobalNumbering*[dim];
  apf::Numbering** lnums = new apf::Numbering*[dim];
  for (int i=0;i<=dim;i++) {
    gnums[i] = m->getGlobalNumbering(i);
    lnums[i] = m->getNumbering(i);
  }

  for (int i=0;i<=dim;i++) {
    const zgid_t* gids;
    const Zoltan2::EntityTopologyType* topTypes;
    const scalar_t* x_coords;
    const scalar_t* y_coords;
    const scalar_t* z_coords;
    int x_stride;
    int y_stride;
    int z_stride;
    const lno_t** offsets = new const lno_t*[dim];
    const zgid_t** adj_gids = new const zgid_t*[dim];

    ia.getIDsViewOf(ents[i],gids);
    ia.getTopologyViewOf(ents[i],topTypes);
    ia.getCoordinatesViewOf(ents[i],x_coords,x_stride,0);
    ia.getCoordinatesViewOf(ents[i],y_coords,y_stride,1);
    ia.getCoordinatesViewOf(ents[i],z_coords,z_stride,2);
    //Check availability of First Adjacency
    for (int j=0;j<=dim;j++) {
      if (ia.availAdjs(ents[i],ents[j])!=(i!=j)) {
        std::cerr<<"First Adjacency does not exist from "<<i<<" to "<< j<<"\n";
        return 5;
      }
      ia.getAdjsView(ents[i],ents[j],offsets[j],adj_gids[j]);
    }
    int j=0;
    apf::MeshIterator* itr = m->begin(i);
    apf::MeshEntity* ent;
    size_t* numAdjs = new size_t[dim+1];
    for (int k=0;k<=dim;k++)
      numAdjs[k]=0;
    while ((ent=m->iterate(itr))) {
      //Check Local ID numbers
      if (apf::getNumber(lnums[i],ent,0,0)!=j) {
        std::cerr<<"Local numbering does not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 2;
      } 

      //Check Global Id numbers
      if (apf::getNumber(gnums[i],apf::Node(ent,0))!=gids[j]) {
        std::cerr<<"Global numbering does not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 2;
      }

      //Check Topology Types
      if (topologyAPFtoZ2(m->getType(ent))!=topTypes[j]) {
        std::cerr<<"Topology types do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 3;
      }
      
      //Check the coordinates
      apf::Vector3 pnt;
      if (i==0)
        m->getPoint(ent,0,pnt);
      else
        pnt = apf::getLinearCentroid(m,ent);
      float eps=.00001;
      if (fabs(pnt[0] - x_coords[j*x_stride])>eps) {
        std::cerr<<"X coordinate do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 4;
      }
      if (fabs(pnt[1] - y_coords[j*y_stride])>eps) {
        std::cerr<<"Y coordinate do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 4;
      }
      if (fabs(pnt[2] - z_coords[j*z_stride])>eps) {
        std::cerr<<"Z coordinate do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 4;
      }  
      
      //Check first adjacencies
      for (int k=0;k<=dim;k++) {
        if (i==k) {
          //Check first adjacency to self is set to NULL
          if (offsets[k]!=NULL || adj_gids[k]!=NULL) {
            std::cerr<<"[WARNING] First adjacency to self is not set to NULL in dimension "<<i
                     <<" to dimension "<<k<<"\n";
          }
            
          continue;
        }
        apf::Adjacent adjs;
        m->getAdjacent(ent,k,adjs);
        lno_t ind = offsets[k][j];
        for (unsigned int l=0;l<adjs.getSize();l++) {
          if (apf::getNumber(gnums[k],apf::Node(adjs[l],0))!=adj_gids[k][ind]) {
            std::cerr<<"First adjacency does not match in dimension " <<i<<" to dimension "<<k
                     <<" on entity "<<j<<"\n";
            return 7;
          }
          ind++;
        }
        if (ind!=offsets[k][j+1]) {
          std::cerr<<"First adjacency length does not match in dimension "<<i<<" to dimension "
                   <<k<<" on entity "<<j<<"\n";
          return 8;
        }
        numAdjs[k]+=adjs.getSize();
        
      }
      j++;
    }
    m->end(itr);
    delete [] offsets;
    delete [] adj_gids;
    //Check the number of first adjacency
    for (int k=0;k<=dim;k++) {
      if (ia.getLocalNumAdjs(ents[i],ents[k])!=numAdjs[k]) {
        std::cerr<<"Local number of first adjacencies do not match in dimension "<<i
                 <<" through dimension "<<k<<"\n";
        return 6;
      }
    }
    delete [] numAdjs;
   
  }
  delete [] lnums;
  delete [] gnums;

  
  //Delete the MeshAdapter
  ia.destroy();

  //Delete the APF Mesh
  m->destroyNative();
  apf::destroyMesh(m);
  
  //End PCU communications
  PCU_Comm_Free();
#endif
  std::cout<<"PASS\n";
  /***************************************************************************/
  /***************************** GENERATE MESH *******************************/
  /***************************************************************************/
  /*
  if (me == 0) cout << "Generating mesh ... \n\n";

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh(meshInput.c_str(), dim, me, numProcs, maxInt);

  // Creating mesh adapter
  if (me == 0) cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::PamgenMeshAdapter<tMVector_t> inputAdapter_t;

  inputAdapter_t ia(*CommT, "region");
  inputAdapter_t::zgid_t const *adjacencyIds=NULL;
  inputAdapter_t::lno_t const *offsets=NULL;
  ia.print(me);
  Zoltan2::MeshEntityType primaryEType = ia.getPrimaryEntityType();
  Zoltan2::MeshEntityType adjEType = ia.getAdjacencyEntityType();

  if (ia.availAdjs(primaryEType, adjEType)) {
    ia.getAdjsView(primaryEType, adjEType, offsets, adjacencyIds);
    int dimension, num_nodes, num_elem;
    int error = 0;
    char title[100];
    int exoid = 0;
    int num_elem_blk, num_node_sets, num_side_sets;
    error += im_ex_get_init(exoid, title, &dimension, &num_nodes, &num_elem,
			    &num_elem_blk, &num_node_sets, &num_side_sets);

    if ((int)ia.getLocalNumOf(primaryEType) != num_elem) {
      cout << "Number of elements do not match\n";
      return 2;
    }

    int *element_num_map = new int [num_elem];
    error += im_ex_get_elem_num_map(exoid, element_num_map);

    inputAdapter_t::zgid_t *node_num_map = new int [num_nodes];
    error += im_ex_get_node_num_map(exoid, node_num_map);

    int *elem_blk_ids = new int [num_elem_blk];
    error += im_ex_get_elem_blk_ids(exoid, elem_blk_ids);

    int *num_nodes_per_elem = new int [num_elem_blk];
    int *num_attr           = new int [num_elem_blk];
    int *num_elem_this_blk  = new int [num_elem_blk];
    char **elem_type        = new char * [num_elem_blk];
    int **connect           = new int * [num_elem_blk];

    for(int i = 0; i < num_elem_blk; i++){
      elem_type[i] = new char [MAX_STR_LENGTH + 1];
      error += im_ex_get_elem_block(exoid, elem_blk_ids[i], elem_type[i],
				    (int*)&(num_elem_this_blk[i]),
				    (int*)&(num_nodes_per_elem[i]),
				    (int*)&(num_attr[i]));
      delete[] elem_type[i];
    }

    delete[] elem_type;
    elem_type = NULL;
    delete[] num_attr;
    num_attr = NULL;

    for(int b = 0; b < num_elem_blk; b++) {
      connect[b] = new int [num_nodes_per_elem[b]*num_elem_this_blk[b]];
      error += im_ex_get_elem_conn(exoid, elem_blk_ids[b], connect[b]);
    }

    delete[] elem_blk_ids;
    elem_blk_ids = NULL;
    int telct = 0;

    for (int b = 0; b < num_elem_blk; b++) {
      for (int i = 0; i < num_elem_this_blk[b]; i++) {
	if (offsets[telct + 1] - offsets[telct] != num_nodes_per_elem[b]) {
	  std::cout << "Number of adjacencies do not match" << std::endl;
	  return 3;
	}

	for (int j = 0; j < num_nodes_per_elem[b]; j++) {
	  ssize_t in_list = -1;

	  for(inputAdapter_t::lno_t k=offsets[telct];k<offsets[telct+1];k++) {
	    if(adjacencyIds[k] ==
	       node_num_map[connect[b][i*num_nodes_per_elem[b]+j]-1]) {
	      in_list = k;
	      break;
	    }
	  }

	  if (in_list < 0) {
	    std::cout << "Adjacency missing" << std::endl;
	    return 4;
	  }
	}

	++telct;
      }
    }

    if (telct != num_elem) {
      cout << "Number of elements do not match\n";
      return 2;
    }
  }
  else{
    std::cout << "Adjacencies not available" << std::endl;
    return 1;
  }

  // delete mesh
  if (me == 0) cout << "Deleting the mesh ... \n\n";

  Delete_Pamgen_Mesh();

  if (me == 0)
    std::cout << "PASS" << std::endl;
  */
  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
