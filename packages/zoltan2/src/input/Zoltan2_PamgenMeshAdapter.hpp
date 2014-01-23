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

#include <Zoltan2_MeshInput.hpp>
#include <Zoltan2_StrideData.hpp>
#include <vector>

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
    if (MESH_REGION == etype) {
      return RnumIds_;
    }

    if (MESH_FACE == etype) {
      return FnumIds_;
    }

    if (MESH_EDGE == etype) {
      return EnumIds_;
    }

    if (MESH_VERTEX == etype) {
      return VnumIds_:
    }
  }
   
  size_t getIDsViewOf(MeshEntityType etype, const gid_t *&Ids) const
  {
    if (MESH_REGION == etype) {
      Ids = RidList_;
      return RnumIds_;
    }

    if (MESH_FACE == etype) {
      Ids = FidList_;
      return FnumIds_;
    }

    if (MESH_EDGE == etype) {
      Ids = EidList_;
      return EnumIds_;
    }

    if (MESH_VERTEX == etype) {
      Ids = VidList_;
      return VnumIds_;
    }
  }

  int getDimensionOf(MeshEntityType etype) const { return dimension_; }

  void getCoordinatesViewOf(MeshEntityType etype, const scalar_t *&coords,
			    int &stride, int dim} const {
    if (dim != dimension_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ";" <<__LINE__
	   << "  Invalid dimension " << dim << std::endl;
      throw std::runtime_error(emsg.str());
    }

    if (MESH_REGION == etype) {
      coords = Rcoords_;
      stride = 1;
    }

    if (MESH_FACE == etype) {
      coords = Fcoords_;
      stride = 1;
    }

    if (MESH_EDGE == etype) {
      coords = Ecoords_;
      stride = 1;
    }

    if (MESH_VERTEX == etype) {
      coords = Vcoords_;
      stride = 1;
    }
  }

private:
  lno_t RnumIds_, FnumIds_, EnumIds_, VnumIds_;
  const gid_t *RidList_, *FidList_, *EidList_, *VidList_;

  int dimension_;
  double * Rcoords_, Fcoords_, Ecoords_, Vcoords_, Acoords_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
PamgenMeshAdapter<User>::PamgenMeshAdapter(string typestr = "region"):
  dimension_(0)
{
  setPrimaryEntityType(typestr);

  int error = 0;
  int exoid = 0;
  int num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  error += im_ex_get_init ( exoid, "PAMGEN Inline Mesh", &dimension_,
			    &num_nodes, &num_elem, &num_elem_blk,
			    &num_node_sets, &num_side_sets);

  Vcoords_ = (double *)malloc(num_nodes * dimension_ * sizeof(double));

  error += im_ex_get_coord(exoid, Vcoords_, Vcoords_ + num_nodes,
			   Vcoords_ + 2 * num_nodes);

  if (3 == dimension_ && num_elem) {
    int * element_num_map = (int *)malloc(num_elem * sizeof(int));
    error += im_ex_get_elem_num_map(exoid, element_num_map);

    RnumIds_ = num_elem;
    RidList_ = element_num_map;
  } else {
    RnumIds_ = 0;
    RidList_ = NULL;
  }

  if (2 == dimension_ && num_elem) {
    int * element_num_map = (int *)malloc(num_elem * sizeof(int));
    error += im_ex_get_elem_num_map(exoid, element_num_map);

    FnumIds_ = num_elem;
    FidList_ = element_num_map;
  } else {
    FnumIds_ = 0;
    FidList_ = NULL;
  }

  EnumIds_ = 0;
  EidList_ = NULL;

  if (num_nodes) {
    int * node_num_map = (int *)malloc(num_nodes * sizeof(int));
    error += im_ex_get_node_num_map(exoid, node_num_map);

    VnumIds_ = num_nodes;
    VidList_ = node_num_map;
  } else {
    VnumIds_ = 0;
    VidList_ = NULL;
  }

  int * elem_blk_ids       = (int *)malloc(num_elem_blk * sizeof(int));
  int * num_nodes_per_elem = (int *)malloc(num_elem_blk * sizeof(int));
  int * num_attr           = (int *)malloc(num_elem_blk * sizeof(int));
  int * num_elem_this_blk  = (int *)malloc(num_elem_blk * sizeof(int));
  char ** elem_type        = (char **)malloc(num_elem_blk * sizeof(char *));
  int ** connect           = (int **)malloc(num_elem_blk * sizeof(int *));

  error += im_ex_get_elem_blk_ids(exoid, elem_blk_ids);

  for(int i = 0; i < num_elem_blk; i++){
    char * elem_type[i] = (char *)malloc((MAX_STR_LENGTH + 1) * sizeof(char));
    error += im_ex_get_elem_block(exoid, elem_blk_id[i], elem_type[i],
				  (int *)&(num_elem_this_blk[i]),
				  (int *)&(num_nodes_per_elem[i]),
				  (int *)&(num_attr[i]));
  }

  Acoords_ = (double *)malloc(num_elem * dimension_ * sizeof(double));
  int a = 0;

  for(int b = 0; b < num_elem_blk; b++){
    int * connect[b] = (int *)malloc(num_nodes_per_elem[b] *
				    num_elem_this_blk[b] * sizeof(int));
    error += im_ex_get_elem_conn(exoid, elem_blk_id[b], connect[b]);

    for(int i = 0; i < num_elem_this_blk[b]; i++){
      Acoords_[a] = 0;
      Acoords_[num_nodes + a] = 0;

      if (3 == dimension_) {
	Acoords_[2 * num_nodes + a] = 0;
      }

      for(int j = 0; j < num_nodes_per_elem[b]; j++){
	Acoords_[a] +=
	  Vcoords_[connect[b][i*num_elem_this_blk[b]+num_nodes_per_elem[b]] - 1];
	Acoords_[num_nodes + a] +=
	  Vcoords_[connect[b]
		 [num_nodes+i*num_elem_this_blk[b]+num_nodes_per_elem[b]] - 1];

	if(3 == dimension_) {
	  Acoords_[2 * num_nodes + a] +=
	    Vcoords_[connect[b]
		   [2*num_nodes+i*num_elem_this_blk[b]+num_nodes_per_elem[b]] -
		   1];
	}
      }

      Acoords_[a] /= num_nodes_per_elem[b];
      Acoords_[num_nodes + a] /= num_nodes_per_elem[b];

      if(3 == dimension_) {
	Acoords_[2 * num_nodes + a] /= num_nodes_per_elem[b];
      }

      a++;
    }
  }

  if (3 == dimension_) {
    Rcoords_ = Acoords_;
  } else {
    Rcoords_ = NULL;
  }

  if (2 == dimension_) {
    Fcoords_ = Acoords_;
  } else {
    Fcoords_ = NULL;
  }

  Ecoords_ = NULL;

}

  
  
}  //namespace Zoltan2
  
#endif
