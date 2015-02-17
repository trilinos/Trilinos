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

/*! \file Zoltan2_PamgenMeshAdapter.hpp
    \brief Defines the PamgenMeshAdapter class.
*/

#ifndef _ZOLTAN2_PAMGENMESHADAPTER_HPP_
#define _ZOLTAN2_PAMGENMESHADAPTER_HPP_

#include <Zoltan2_MeshAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <vector>

#include <pamgen_im_exodusII.h>
#include "pamgen_im_ne_nemesisI.h"

namespace Zoltan2 {

/*! \brief This class represents a mesh.
 *
 *  A mesh can be a collection of global Identifiers
 *           and their associated weights, if any.
 *
 *  The user supplies the identifiers and weights by way of pointers
 *    to arrays.  
 *
    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent coordinates, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

 */

template <typename User>
  class PamgenMeshAdapter: public MeshAdapter<User> {

public:

  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::zgid_t    zgid_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef MeshAdapter<User>       base_adapter_t;
  typedef User user_t;

  /*! \brief Constructor for mesh with identifiers but no coordinates or edges
   *  \param etype is the mesh entity type of the identifiers
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */

  PamgenMeshAdapter(const Comm<int> &comm, std::string typestr="region");

  void print(int);

  ////////////////////////////////////////////////////////////////
  // The MeshAdapter interface.
  // This is the interface that would be called by a model or a problem .
  ////////////////////////////////////////////////////////////////

  size_t getGlobalNumOf(MeshEntityType etype) const
  {

    if ((MESH_REGION == etype && 3 == dimension_) ||
	(MESH_FACE == etype && 2 == dimension_)) {
      return num_elems_global_;
    }

    if (MESH_VERTEX == etype) {
      return num_nodes_global_;
    }

    return 0;
  }

  size_t getLocalNumOf(MeshEntityType etype) const
  {
    if ((MESH_REGION == etype && 3 == dimension_) ||
	(MESH_FACE == etype && 2 == dimension_)) {
      return num_elem_;
    }

    if (MESH_VERTEX == etype) {
      return num_nodes_;
    }

    return 0;
  }
   
  void getIDsViewOf(MeshEntityType etype, const zgid_t *&Ids) const
  {
    if ((MESH_REGION == etype && 3 == dimension_) ||
	(MESH_FACE == etype && 2 == dimension_)) {
      Ids = element_num_map_;
    }

    else if (MESH_VERTEX == etype) {
      Ids = node_num_map_;
    }

    else Ids = NULL;
  }

  void getWeightsViewOf(MeshEntityType etype, const scalar_t *&weights,
			int &stride, int idx = 0) const
  {
    weights = NULL;
    stride = 0;
  }

  int getDimension() const { return dimension_; }

  void getCoordinatesViewOf(MeshEntityType etype, const scalar_t *&coords,
			    int &stride, int dim) const {
    if ((MESH_REGION == etype && 3 == dimension_) ||
	       (MESH_FACE == etype && 2 == dimension_)) {
      if (dim == 0) {
	coords = Acoords_;
      } else if (dim == 1) {
	coords = Acoords_ + num_elem_;
      } else if (dim == 2) {
	coords = Acoords_ + 2 * num_elem_;
      }
      stride = 1;
    } else if (MESH_REGION == etype && 2 == dimension_) {
      coords = NULL;
      stride = 0;
    } else if (MESH_VERTEX == etype) {
      if (dim == 0) {
	coords = coords_;
      } else if (dim == 1) {
	coords = coords_ + num_nodes_;
      } else if (dim == 2) {
	coords = coords_ + 2 * num_nodes_;
      }
      stride = 1;
    } else {
      coords = NULL;
      stride = 0;
      Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
    }
  }

  bool availAdjs(MeshEntityType source, MeshEntityType target) const {
    if ((MESH_REGION == source && MESH_VERTEX == target && 3 == dimension_) ||
	(MESH_FACE == source && MESH_VERTEX == target && 2 == dimension_)) {
      return TRUE;
    }

    return false;
  }

  size_t getLocalNumAdjs(MeshEntityType source, MeshEntityType target) const
  {
    if (availAdjs(source, target)) {
      return tnoct_;
    }
    
    return 0;
  }

  void getAdjsView(MeshEntityType source, MeshEntityType target,
		   const lno_t *&offsets, const zgid_t *& adjacencyIds) const
  {
    if ((MESH_REGION == source && MESH_VERTEX == target && 3 == dimension_) ||
	(MESH_FACE == source && MESH_VERTEX == target && 2 == dimension_)) {
      offsets = elemOffsets_;
      adjacencyIds = elemToNode_;
    } else if (MESH_REGION == source && 2 == dimension_) {
      offsets = NULL;
      adjacencyIds = NULL;
    } else {
      offsets = NULL;
      adjacencyIds = NULL;
      Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
    }
  }

#define USE_MESH_ADAPTER
#ifndef USE_MESH_ADAPTER
  bool avail2ndAdjs(MeshEntityType sourcetarget, MeshEntityType through) const
  {
    if (through == MESH_VERTEX) {
      if (sourcetarget == MESH_REGION && dimension_ == 3) return true;
      if (sourcetarget == MESH_FACE && dimension_ == 2) return true;
    }
    return false;
  }

  size_t getLocalNum2ndAdjs(MeshEntityType sourcetarget, 
			    MeshEntityType through) const
  {
    if (avail2ndAdjs(sourcetarget, through)) {
      return nadj_;
    }

    return 0;
  }

  void get2ndAdjsView(MeshEntityType sourcetarget, MeshEntityType through, 
		      const lno_t *&offsets, const zgid_t *&adjacencyIds) const
  {
    if (avail2ndAdjs(sourcetarget, through)) {
      offsets = start_;
      adjacencyIds = adj_;
    } else {
      offsets = NULL;
      adjacencyIds = NULL;
      Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
    }
  }
#endif

private:
  int dimension_, num_nodes_global_, num_elems_global_, num_nodes_, num_elem_;
  zgid_t *element_num_map_, *node_num_map_;
  int *elemToNode_, tnoct_, *elemOffsets_;
  double *coords_, *Acoords_;
  lno_t *start_;
  zgid_t *adj_;
  size_t nadj_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

static
ssize_t in_list(const int value, size_t count, int *vector)
{
  for(size_t i=0; i < count; i++) {
    if(vector[i] == value)
      return i;
  }
  return -1;
}

template <typename User>
PamgenMeshAdapter<User>::PamgenMeshAdapter(const Comm<int> &comm,
					   std::string typestr):
  dimension_(0)
{
  this->setEntityTypes(typestr, "vertex", "vertex");

  int error = 0;
  char title[100];
  int exoid = 0;
  int num_elem_blk, num_node_sets, num_side_sets;
  error += im_ex_get_init(exoid, title, &dimension_,
			  &num_nodes_, &num_elem_, &num_elem_blk,
			  &num_node_sets, &num_side_sets);

  coords_ = new double [num_nodes_ * dimension_];

  error += im_ex_get_coord(exoid, coords_, coords_ + num_nodes_,
			   coords_ + 2 * num_nodes_);

  element_num_map_ = new int [num_elem_];
  error += im_ex_get_elem_num_map(exoid, element_num_map_);

  node_num_map_ = new int [num_nodes_];
  error += im_ex_get_node_num_map(exoid, node_num_map_);

  int *elem_blk_ids       = new int [num_elem_blk];
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
  Acoords_ = new double [num_elem_ * dimension_];
  int a = 0;
  std::vector<std::vector<int> > sur_elem;
  sur_elem.resize(num_nodes_);

  for(int b = 0; b < num_elem_blk; b++) {
    connect[b] = new int [num_nodes_per_elem[b]*num_elem_this_blk[b]];
    error += im_ex_get_elem_conn(exoid, elem_blk_ids[b], connect[b]);

    for(int i = 0; i < num_elem_this_blk[b]; i++) {
      Acoords_[a] = 0;
      Acoords_[num_elem_ + a] = 0;

      if (3 == dimension_) {
	Acoords_[2 * num_elem_ + a] = 0;
      }

      for(int j = 0; j < num_nodes_per_elem[b]; j++) {
	int node = connect[b][i * num_nodes_per_elem[b] + j] - 1;
	Acoords_[a] += coords_[node];
	Acoords_[num_elem_ + a] += coords_[num_nodes_ + node];

	if(3 == dimension_) {
	  Acoords_[2 * num_elem_ + a] += coords_[2 * num_nodes_ + node];
	}

	/*
	 * in the case of degenerate elements, where a node can be
	 * entered into the connect table twice, need to check to
	 * make sure that this element is not already listed as
	 * surrounding this node
	 */
	if (sur_elem[node].empty() ||
	    element_num_map_[a] != sur_elem[node][sur_elem[node].size()-1]) {
	  /* Add the element to the list */
	  sur_elem[node].push_back(element_num_map_[a]);
	}
      }

      Acoords_[a] /= num_nodes_per_elem[b];
      Acoords_[num_elem_ + a] /= num_nodes_per_elem[b];

      if(3 == dimension_) {
	Acoords_[2 * num_elem_ + a] /= num_nodes_per_elem[b];
      }

      a++;
    }

  }

  delete[] elem_blk_ids;
  elem_blk_ids = NULL;
  int nnodes_per_elem = num_nodes_per_elem[0];
  elemToNode_ = new int [num_elem_ * nnodes_per_elem];
  int telct = 0;
  elemOffsets_ = new int [num_elem_+1];
  tnoct_ = 0;
  int **reconnect = new int * [num_elem_];
  size_t max_nsur = 0;

  for (int b = 0; b < num_elem_blk; b++) {
    for (int i = 0; i < num_elem_this_blk[b]; i++) {
      elemOffsets_[telct] = tnoct_;
      reconnect[telct] = new int [num_nodes_per_elem[b]];

      for (int j = 0; j < num_nodes_per_elem[b]; j++) {
	elemToNode_[tnoct_]=
	  node_num_map_[connect[b][i*num_nodes_per_elem[b] + j]-1];
	reconnect[telct][j] = connect[b][i*num_nodes_per_elem[b] + j];
	++tnoct_;
      }

      ++telct;
    }
  }
  elemOffsets_[telct] = tnoct_;

  delete[] num_nodes_per_elem;
  num_nodes_per_elem = NULL;
  delete[] num_elem_this_blk;
  num_elem_this_blk = NULL;

  for(int b = 0; b < num_elem_blk; b++) {
    delete[] connect[b];
  }

  delete[] connect;
  connect = NULL;

  int max_side_nodes = nnodes_per_elem;
  int *side_nodes = new int [max_side_nodes];
  int *mirror_nodes = new int [max_side_nodes];

  /* Allocate memory necessary for the adjacency */
  start_ = new lno_t [num_elem_+1];
  std::vector<int> adj;

  for (int i=0; i < max_side_nodes; i++) {
    side_nodes[i]=-999;
    mirror_nodes[i]=-999;
  }

  /* Find the adjacency for a nodal based decomposition */
  nadj_ = 0;
  for(int ncnt=0; ncnt < num_nodes_; ncnt++) {
    if(sur_elem[ncnt].empty()) {
      printf("WARNING: Node = %d has no elements\n", ncnt+1);
    } else {
      size_t nsur = sur_elem[ncnt].size();
      if (nsur > max_nsur)
	max_nsur = nsur;
    }
  }

  int nprocs = comm.getSize();
  //if (nprocs > 1) {
    int neid=0,num_elem_blks_global,num_node_sets_global,num_side_sets_global;
    error += im_ne_get_init_global(neid,&num_nodes_global_,&num_elems_global_,
				   &num_elem_blks_global,&num_node_sets_global,
				   &num_side_sets_global);

    int num_internal_nodes, num_border_nodes, num_external_nodes;
    int num_internal_elems, num_border_elems, num_node_cmaps, num_elem_cmaps;
    int proc = 0;
    error += im_ne_get_loadbal_param(neid, &num_internal_nodes,
				     &num_border_nodes, &num_external_nodes,
				     &num_internal_elems, &num_border_elems,
				     &num_node_cmaps, &num_elem_cmaps, proc);

    int *node_cmap_ids = new int [num_node_cmaps];
    int *node_cmap_node_cnts = new int [num_node_cmaps];
    int *elem_cmap_ids = new int [num_elem_cmaps];
    int *elem_cmap_elem_cnts = new int [num_elem_cmaps];
    error += im_ne_get_cmap_params(neid, node_cmap_ids, node_cmap_node_cnts,
				   elem_cmap_ids, elem_cmap_elem_cnts, proc);
    delete[] elem_cmap_ids;
    elem_cmap_ids = NULL;
    delete[] elem_cmap_elem_cnts;
    elem_cmap_elem_cnts = NULL;

    int **node_ids = new int * [num_node_cmaps];
    int **node_proc_ids = new int * [num_node_cmaps];
    for(int j = 0; j < num_node_cmaps; j++) {
      node_ids[j] = new int [node_cmap_node_cnts[j]];
      node_proc_ids[j] = new int [node_cmap_node_cnts[j]];
      error += im_ne_get_node_cmap(neid, node_cmap_ids[j], node_ids[j],
				   node_proc_ids[j], proc);
    }
    delete[] node_cmap_ids;
    node_cmap_ids = NULL;
    int *sendCount = new int [nprocs];
    int *recvCount = new int [nprocs];

    // Post receives
    RCP<CommRequest<int> >*requests=new RCP<CommRequest<int> >[num_node_cmaps];
    for (int cnt = 0, i = 0; i < num_node_cmaps; i++) {
      try {
	requests[cnt++] =
	  Teuchos::ireceive<int,int>(comm,
				     rcp(&(recvCount[node_proc_ids[i][0]]),
					 false),
				     node_proc_ids[i][0]);
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    Teuchos::barrier<int>(comm);
    size_t totalsend = 0;

    for(int j = 0; j < num_node_cmaps; j++) {
      sendCount[node_proc_ids[j][0]] = 1;
      for(int i = 0; i < node_cmap_node_cnts[j]; i++) {
	sendCount[node_proc_ids[j][i]] += sur_elem[node_ids[j][i]-1].size()+2;
      }
      totalsend += sendCount[node_proc_ids[j][0]];
    }

    // Send data; can use readySend since receives are posted.
    for (int i = 0; i < num_node_cmaps; i++) {
      try {
	Teuchos::readySend<int,int>(comm, sendCount[node_proc_ids[i][0]],
				    node_proc_ids[i][0]);
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    // Wait for messages to return.
    try {
      Teuchos::waitAll<int>(comm, arrayView(requests, num_node_cmaps));
    }
    Z2_FORWARD_EXCEPTIONS;

    delete [] requests;

    // Allocate the receive buffer.
    size_t totalrecv = 0;
    int maxMsg = 0;
    int nrecvranks = 0;
    for(int i = 0; i < num_node_cmaps; i++) {
      if (recvCount[node_proc_ids[i][0]] > 0) {
	totalrecv += recvCount[node_proc_ids[i][0]];
	nrecvranks++;
	if (recvCount[node_proc_ids[i][0]] > maxMsg)
	  maxMsg = recvCount[node_proc_ids[i][0]];
      }
    }

    int *rbuf = NULL;
    if (totalrecv) rbuf = new int[totalrecv];

    requests = new RCP<CommRequest<int> > [nrecvranks];

    // Error checking for memory and message size.
    int OK[2] = {1,1};
    // OK[0] -- true/false indicating whether each message size fits in an int
    //          (for MPI).
    // OK[1] -- true/false indicating whether memory allocs are OK
    int gOK[2]; // For global reduce of OK.

    if (size_t(maxMsg) * sizeof(int) > INT_MAX) OK[0] = false;
    if (totalrecv && !rbuf) OK[1] = 0;
    if (!requests) OK[1] = 0;

    // Post receives

    size_t offset = 0;

    if (OK[0] && OK[1]) {
      int rcnt = 0;
      for (int i = 0; i < num_node_cmaps; i++) {
	if (recvCount[node_proc_ids[i][0]]) {
	  try {
	    requests[rcnt++] =
	      Teuchos::
	      ireceive<int,int>(comm,
				Teuchos::arcp(&rbuf[offset], 0,
					      recvCount[node_proc_ids[i][0]],
					      false),
				node_proc_ids[i][0]);
	  }
	  Z2_FORWARD_EXCEPTIONS;
	}
	offset += recvCount[node_proc_ids[i][0]];
      }
    }

    delete[] recvCount;

    // Use barrier for error checking
    Teuchos::reduceAll<int>(comm, Teuchos::REDUCE_MIN, 2, OK, gOK);
    if (!gOK[0] || !gOK[1]) {
      delete [] rbuf;
      delete [] requests;
      if (!gOK[0])
	throw std::runtime_error("Max single message length exceeded");
      else
	throw std::bad_alloc();
    }

    int *sbuf = NULL;
    if (totalsend) sbuf = new int[totalsend];
    a = 0;

    for(int j = 0; j < num_node_cmaps; j++) {
      sbuf[a++] = node_cmap_node_cnts[j];
      for(int i = 0; i < node_cmap_node_cnts[j]; i++) {
	sbuf[a++] = node_num_map_[node_ids[j][i]-1];
	sbuf[a++] = sur_elem[node_ids[j][i]-1].size();
	for(size_t ecnt=0; ecnt < sur_elem[node_ids[j][i]-1].size(); ecnt++) {
	  sbuf[a++] = sur_elem[node_ids[j][i]-1][ecnt];
	}
      }
    }

    delete[] node_cmap_node_cnts;
    node_cmap_node_cnts = NULL;

    for(int j = 0; j < num_node_cmaps; j++) {
      delete[] node_ids[j];
    }

    delete[] node_ids;
    node_ids = NULL;
    ArrayRCP<int> sendBuf;

    if (totalsend)
      sendBuf = ArrayRCP<int>(sbuf, 0, totalsend, true);
    else
      sendBuf = Teuchos::null;

    // Send data; can use readySend since receives are posted.
    offset = 0;
    for (int i = 0; i < num_node_cmaps; i++) {
      if (sendCount[node_proc_ids[i][0]]) {
	try{
	  Teuchos::readySend<int,
	    int>(comm,
		 Teuchos::arrayView(&sendBuf[offset],
				    sendCount[node_proc_ids[i][0]]),
		 node_proc_ids[i][0]);
	}
	Z2_FORWARD_EXCEPTIONS;
      }
      offset += sendCount[node_proc_ids[i][0]];
    }

    for(int j = 0; j < num_node_cmaps; j++) {
      delete[] node_proc_ids[j];
    }

    delete[] node_proc_ids;
    node_proc_ids = NULL;
    delete[] sendCount;

    // Wait for messages to return.
    try{
      Teuchos::waitAll<int>(comm, Teuchos::arrayView(requests, nrecvranks));
    }
    Z2_FORWARD_EXCEPTIONS;

    delete[] requests;
    a = 0;

    for (int i = 0; i < num_node_cmaps; i++) {
      int num_nodes_this_processor = rbuf[a++];

      for (int j = 0; j < num_nodes_this_processor; j++) {
	int this_node = rbuf[a++];
	int num_elem_this_node = rbuf[a++];

	for (int ncnt = 0; ncnt < num_nodes_; ncnt++) {
	  if (node_num_map_[ncnt] == this_node) {
	    for (int ecnt = 0; ecnt < num_elem_this_node; ecnt++) {
	      sur_elem[ncnt].push_back(rbuf[a++]);
	    }
	  }
	}
      }
    }

    delete[] rbuf;
    //}

  for(int ecnt=0; ecnt < num_elem_; ecnt++) {
    start_[ecnt] = nadj_;
    int nnodes = nnodes_per_elem;
    for(int ncnt=0; ncnt < nnodes; ncnt++) {
      int node = reconnect[ecnt][ncnt]-1;
      for(size_t i=0; i < sur_elem[node].size(); i++) {
	int entry = sur_elem[node][i];

	if(element_num_map_[ecnt] != entry &&
	   in_list(entry,
		   adj.size()-start_[ecnt],
		   &adj[start_[ecnt]]) < 0) {
	  adj.push_back(entry);
	  nadj_++;
	}
      }
    }
  }

  for(int b = 0; b < num_elem_; b++) {
    delete[] reconnect[b];
  }

  delete[] reconnect;
  reconnect = NULL;
  start_[num_elem_] = nadj_;

  adj_ = new zgid_t [nadj_];

  for (size_t i=0; i < nadj_; i++) {
    adj_[i] = adj[i];
  }

  delete[] side_nodes;
  side_nodes = NULL;
  delete[] mirror_nodes;
  mirror_nodes = NULL;
}

template <typename User>
void PamgenMeshAdapter<User>::print(int me)
{
  std::string fn(" PamgenMesh ");
  std::cout << me << fn
            << " dim = " << dimension_
            << " nnodes = " << num_nodes_
            << " nelems = " << num_elem_
            << std::endl;

  for (int i = 0; i < num_elem_; i++) {
    std::cout << me << fn << i 
              << " Elem " << element_num_map_[i]
              << " Coords: ";
    for (int j = 0; j < dimension_; j++)
      std::cout << Acoords_[i + j * num_elem_] << " ";
    std::cout << std::endl;
  }

  for (int i = 0; i < num_elem_; i++) {
    std::cout << me << fn << i 
              << " Elem " << element_num_map_[i]
              << " Graph: ";
    for (int j = start_[i]; j < start_[i+1]; j++)
      std::cout << adj_[j] << " ";
    std::cout << std::endl;
  }
}
  
}  //namespace Zoltan2
  
#endif
