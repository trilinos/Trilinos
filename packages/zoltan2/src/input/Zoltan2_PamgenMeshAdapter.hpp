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
// Questions? Contact Vitus Leung       (vjleung@sandia.gov)
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
#include <Zoltan2_StrideData.hpp>
#include <vector>

#include <im_exodusII_l.h

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
  typedef typename InputTraits<User>::gid_t    gid_t;
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

  PamgenMeshAdapter(string typestr = "region");

  ////////////////////////////////////////////////////////////////
  // The MeshAdapter interface.
  // This is the interface that would be called by a model or a problem .
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumOf(MeshEntityType etype) const
  {
    if (MESH_REGION == etype && 3 == dimension_ ||
	MESH_FACE == etype && 2 == dimension_) {
      return num_elem_;
    }

    if (MESH_VERTEX == etype) {
      return num_nodes_:
    }

    return 0;
  }
   
  size_t getIDsViewOf(MeshEntityType etype, const gid_t *&Ids) const
  {
    if (MESH_REGION == etype && 3 == dimension_ ||
	MESH_FACE == etype && 2 == dimension_) {
      Ids = element_num_map_;
      return num_elem_;
    }

    if (MESH_VERTEX == etype) {
      Ids = node_num_map_;
      return num_nodes_;
    }

    Ids = NULL;
    return 0;
  }

  void getWeigthsViewOf(MeshEntityType etype, const scalar_t *&weights,
			int &stride, int idx = 0) const
  {
    weights = NULL;
    stride = 0;
  }

  int getDimensionOf() const { return dimension_; }

  void getCoordinatesViewOf(MeshEntityType etype, const scalar_t *&coords,
			    int &stride, int dim} const {
    if (dim != dimension_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ";" <<__LINE__
	   << "  Invalid dimension " << dim << std::endl;
      throw std::runtime_error(emsg.str());
    } else if (MESH_REGION == etype && 3 == dimension_ ||
	       MESH_FACE == etype && 2 == dimension_) {
      coords = Acoords_;
      stride = 1;
    } else if (MESH_REGION == etype && 2 == dimension_) {
      coords = NULL;
      stride = 0;
    } else if (MESH_VERTEX == etype) {
      coords = coords_;
      stride = 1;
    } else {
      coords = NULL;
      stride = 0;
      Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
    }
  }

  bool availAdjs(MeshEntityType source, MeshEntityType target) {
    if (MESH_REGION == source && MESH_VERTEX == target && 3 == dimension_ ||
	MESH_FACE == source && MESH_VERTEX == target && 2 == dimension_) {
      return TRUE;
    }

    return FALSE;
  }

  size_t getLocalNumAdjs(MeshEntityType source, MeshEntityType target) const
  {
    if (availAdjs(source, target)) {
      return tnoct_;
    }

    return 0;
  }

  void getAdjsView(MeshEntityType source, MeshEntityType target,
		   const lno_t *&offsets, const gid_t *& adjacencyIds) const
  {
    if (MESH_REGION == source && MESH_VERTEX == target && 3 == dimension_ ||
	MESH_FACE == source && MESH_VERTEX == target && 2 == dimension) {
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

  bool avail2ndAdjs(MeshEntityType sourcetarget, MeshEntityType through)
  {
    if (MESH_REGION==sourcetarget && MESH_VERTEX==through && 3==dimension_ ||
	MESH_FACE==sourcetarget && MESH_VERTEX==through && 2==dimension_) {
      return TRUE;
    }

    return FALSE;
  }

  size_t getLocalNum2ndAdjs(MeshEntityType sourcetarget, 
			    MeshEntityType through) const
  {
    if (avail2ndAdjs(sourcetarget, through)) {
      return adj_.size();
    }

    return 0;
  }

  void get2ndAdjsView(MeshEntityType sourcetarget, MeshEntityType through, 
		      const lno_t *&offsets, const gid_t *& adjacencyIds) const
  {
    if (MESH_REGION==sourcetarget && MESH_VERTEX==through && 3==dimension_ ||
	MESH_FACE==sourcetarget && MESH_VERTEX==through && 2==dimension_) {
      offsets = start_;
      adjacencyIds = adj_;
    } else if (MESH_REGION == sourcetarget && 2 == dimension_) {
      offsets = NULL;
      adjacencyIds = NULL;
    } else {
      offsets = NULL;
      adjacencyIds = NULL;
      Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
    }
  }

private:
  long long dimension_, num_nodes_, num_elem_, *element_num_map_;
  long long *node_num_map_, *elemToNode_, tnoct_, *elemOffsets_;
  double *coords_, *Acoords_;
  std::vector<long long> start_, adj_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

  ssize_t in_list(const long long value, size_t count, long long *vector)
  {
    for(size_t i=0; i < count; i++) {
      if(vector[i] == value)
	return i;
    }
    return -1;
  }

template <typename User>
PamgenMeshAdapter<User>::PamgenMeshAdapter(string typestr = "region"):
  dimension_(0)
{
  setPrimaryEntityType(typestr);

  int error = 0;
  int exoid = 0;
  long long num_elem_blk, num_node_sets, num_side_sets;
  im_ex_get_init_l ( exoid, "PAMGEN Inline Mesh", &dimension_,
		     &num_nodes_, &num_elem_, &num_elem_blk,
		     &num_node_sets, &num_side_sets);

  coords_ = new double [num_nodes_ * dimension_];

  error += im_ex_get_coord_l(exoid, coords_, coords_ + num_nodes_,
			     coords_ + 2 * num_nodes_);

  *element_num_map_ = new long long [num_elem_];
  error += im_ex_get_elem_num_map_l(exoid, element_num_map_);

  *node_num_map_ = new long long [num_nodes_];
  error += im_ex_get_node_num_map_l(exoid, node_num_map_);

  long long *elem_blk_ids       = new long long [num_elem_blk];
  error += im_ex_get_elem_blk_ids_l(exoid, elem_blk_ids);

  long long *num_nodes_per_elem = new long long [num_elem_blk];
  long long *num_attr           = new long long [num_elem_blk];
  long long *num_elem_this_blk  = new long long [num_elem_blk];
  char **elem_type              = new char * [num_elem_blk];
  long long **connect           = new long long * [num_elem_blk];

  for(long long i = 0; i < num_elem_blk; i++){
    elem_type[i] = new char [MAX_STR_LENGTH + 1];
    error += im_ex_get_elem_block_l(exoid, elem_blk_id[i], elem_type[i],
				    (long long*)&(num_elem_this_blk[i]),
				    (long long*)&(num_nodes_per_elem[i]),
				    (long long*)&(num_attr[i]));
  }

  Acoords_ = new double [num_elem_ * dimension_];
  long long a = 0;
  std::vector<std::vector<long long> > sur_elem;
  sur_elem.resize(num_nodes_);

  for(long long b = 0; b < num_elem_blk; b++) {
    connect[b] = new long long [num_nodes_per_elem[b]*num_elem_this_blk[b]];
    error += im_ex_get_elem_conn_l(exoid, elem_blk_id[b], connect[b]);

    for(long long i = 0; i < num_elem_this_blk[b]; i++) {
      Acoords_[a] = 0;
      Acoords_[num_nodes_ + a] = 0;

      if (3 == dimension_) {
	Acoords_[2 * num_nodes_ + a] = 0;
      }

      for(long long j = 0; j < num_nodes_per_elem[b]; j++) {
	long long node = connect[b][i * num_nodes_per_elem[b] + j] - 1;
	Acoords_[a] += coords_[node];
	Acoords_[num_nodes_ + a] += coords_[num_nodes_ + node];

	if(3 == dimension_) {
	  Acoords_[2 * num_nodes_ + a] += coords_[2 * num_nodes_ + node];
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
      Acoords_[num_nodes_ + a] /= num_nodes_per_elem[b];

      if(3 == dimension_) {
	Acoords_[2 * num_nodes_ + a] /= num_nodes_per_elem[b];
      }

      a++;
    }
  }

  long long nnodes_per_elem = num_nodes_per_elem[0];
  elemToNode_ = new long long [num_elem_ * nnodes_per_elem];
  long long telct = 0;
  elemOffsets_ = new long long [num_elem_];
  tnoct_ = 0;
  long long **reconnect = new long long * [num_elem_];
  size_t max_nsur = 0;

  for (long long b = 0; b < num_elem_blk; b++) {
    for (long long i = 0; i < num_elem_this_blk[b]; i++) {
      elemOffsets_[telct] = tnoct_;
      reconnect[telct] = new long long [num_nodes_per_elem[b]];

      for (long long j = 0; j < num_nodes_per_elem[b]; j++) {
	elemToNode_[tnoct_] = connect[b][i*num_nodes_per_elem[b] + j]-1;
	reconnect[telct][j] = connect[b][i*num_nodes_per_elem[b] + j]-1;

	if(sur_elem[tnoct_].empty()) {
	  printf("WARNING: Node = "ST_ZU" has no elements\n", tnoct_+1);
	} else {
	  size_t nsur = sur_elem[tnoct_].size();
	  if (nsur > max_nsur)
	    max_nsur = nsur;
	}

	++tnoct_;
      }

      ++telct;
    }
  }

  long long max_side_nodes = nnodes_per_elem;
  long long side_nodes[max_side_nodes];
  long long mirror_nodes[max_side_nodes];

  /* Allocate memory necessary for the adjacency */
  start_.resize(num_nodes_);

  for (int i=0; i < max_side_nodes; i++) {
    side_nodes[i]=-999;
    mirror_nodes[i]=-999;
  }

  /* Find the adjacency for a nodal based decomposition */
  size_t nadj = 0;
  for(size_t ncnt=0; ncnt < num_nodes_; ncnt++) {
    start_[ncnt] = nadj;
    for(size_t ecnt=0; ecnt < sur_elem[ncnt].size(); ecnt++) {
      size_t elem = sur_elem[ncnt][ecnt];
      int nnodes = nnodes_per_elem;
      for(int i=0; i < nnodes; i++) {
	long long entry = reconnect[elem][i];

	if(ncnt != (size_t)entry &&
	   in_list(entry,
		   adj_.size()-start_[ncnt],
		   &adj_[start_[ncnt]]) < 0) {
	  adj_.push_back(entry);
	}
      }
    }
    nadj++;
  }

  delete[] elem_blk_ids;
  elem_blk_ids = NULL;
  delete[] num_nodes_per_elem;
  num_nodes_per_elem = NULL;
  delete[] num_attr;
  num_attr = NULL;
  delete[] num_elem_this_blk;
  num_elem_this_blk = NULL;

  for(long long i = 0; i < num_elem_blk; i++){
    delete[] elem_type[i];
  }

  delete[] elem_type;
  elem_type = NULL;

  for(long long b = 0; b < num_elem_blk; b++) {
    delete[] connect[b];
  }

  delete[] connect;
  connect = NULL;

  for(long long b = 0; b < num_elem_; b++) {
    delete[] reconnect[b];
  }

  delete[] reconnect;
  reconnect = NULL;
}

  
  
}  //namespace Zoltan2
  
#endif
