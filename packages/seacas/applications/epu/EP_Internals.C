/*
 * Copyright(C) 2010 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#ifndef PRId64
#error "PRId64 not defined"
#endif

#include <EP_ExodusEntity.h>            // for Block, Mesh
#include <EP_Internals.h>               // for Internals, Redefine

#include <stdint.h>                 // for int64_t
#include <smart_assert.h>               // for SMART_ASSERT
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf
#include <stdlib.h>                     // for exit, EXIT_FAILURE
#include <algorithm>                    // for sort
#include <cstring>                      // for strncpy, strlen, memset

#include <iostream>                     // for operator<<, basic_ostream, etc
#include <string>                       // for string, basic_string
#include <vector>                       // for vector

extern "C" {
#define NO_NETCDF_2
#include <netcdf.h>
#include <exodusII.h>
#include <exodusII_int.h>
}

// Explicit Initialization of the functions used...
namespace Excn {
template int Internals<int>::write_meta_data(const Mesh &mesh,
					       const std::vector<Block>   &blocks,
					       const std::vector<NodeSet<int> > &nodesets,
					       const std::vector<SideSet<int> > &sidesets,
					       const CommunicationMetaData &comm);

template bool Internals<int>::check_meta_data(const Mesh &mesh,
					       const std::vector<Block>   &blocks,
					       const std::vector<NodeSet<int> > &nodesets,
					       const std::vector<SideSet<int> > &sidesets,
					       const CommunicationMetaData &comm);

template int Internals<int64_t>::write_meta_data(const Mesh &mesh,
					       const std::vector<Block>   &blocks,
					       const std::vector<NodeSet<int64_t> > &nodesets,
					       const std::vector<SideSet<int64_t> > &sidesets,
					       const CommunicationMetaData &comm);

template bool Internals<int64_t>::check_meta_data(const Mesh &mesh,
					       const std::vector<Block>   &blocks,
					       const std::vector<NodeSet<int64_t> > &nodesets,
					       const std::vector<SideSet<int64_t> > &sidesets,
					       const CommunicationMetaData &comm);
}

namespace {
  bool lessOffset(const Excn::Block& b1, const Excn::Block& b2) {
    return b1.offset_ < b2.offset_;
  }

  nc_type get_type(int exoid, unsigned int type)
  {
    if (ex_int64_status(exoid) & type)
      return NC_INT64;
    else
      return NC_INT;
  }

  int define_netcdf_vars(int exoid, const char *type, size_t count, const char *dim_num,
			 const char *stat_var, const char *id_var, const char *name_var);

  template <typename INT>
  int put_array(int exoid, const char *var_type, const std::vector<INT> &array);

  int put_id_array(int exoid, const char *var_type,  const std::vector<ex_entity_id> &array);

  int define_coordinate_vars(int exodusFilePtr, int64_t nodes, int node_dim,
			     int dimension, int dim_dim, int str_dim);
}
			 
Excn::Redefine::Redefine(int exoid)
  : exodusFilePtr(exoid)
{
  // Enter define mode...
  int status = nc_redef (exodusFilePtr);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg,
	    "Error: failed to put file id %d into define mode", exodusFilePtr);
    ex_err("Excn::Redefine::Redefine()",errmsg,status);
    exit(EXIT_FAILURE);
  }
}

Excn::Redefine::~Redefine()
{
  try {
    int status = nc_enddef (exodusFilePtr);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      char errmsg[MAX_ERR_LENGTH];
      sprintf(errmsg,
	      "Error: failed to complete variable definitions in file id %d",exodusFilePtr);
      ex_err("Excn::Redefine::~Redefine()",errmsg,status);
      exit(EXIT_FAILURE);
    }
  } catch (...) {
  }
}

template <typename INT>
int Excn::Internals<INT>::write_meta_data(const Mesh &mesh,
				     const std::vector<Block>   &blocks,
				     const std::vector<NodeSet<INT> > &nodesets,
				     const std::vector<SideSet<INT> > &sidesets,
				     const CommunicationMetaData &comm)
{
  SMART_ASSERT((int)blocks.size()   == mesh.blockCount);
  SMART_ASSERT((int)nodesets.size() == mesh.nodesetCount);
  SMART_ASSERT((int)sidesets.size() == mesh.sidesetCount);
  
  // May need to reorder the element blocks based on the 'offset_'
  // member. An element block contains the elements from 'offset_+1'
  // to 'offset_+elementCount'.  Typically, this is the order that
  // they appear in the 'blocks' array, but not always...
  //
  // This is a very hard error to track down, so we just copy the
  // 'blocks' array to a new array and sort it on the 'offset_'
  // member...
  //
  // All block-related calls after these meta-data calls are based on
  // the block id and not the block order, so we don't need to fix
  // them. 
  for (size_t i=0; i < blocks.size(); i++) {
    const_cast<Excn::Block&>(blocks[i]).position_ = i;
  }

  // Check whether current order is consistent.  The problem with just
  // sorting the data is that if a subset of the entire model is being
  // joined, there may be zero-length element blocks.  These blocks
  // have an offset of zero since there is no "lowest id" to get the
  // offset from. 
  //
  // If all that is wanted is a single subset of the entire model, it
  // doesn't matter if the blocks are reordered due to the zero
  // offsets; however, if the user then later joins multiple subsets
  // together, there will be issues of block ordering mismatches.
  //
  // So, to avoid sorting the data into incorrect order, we check
  // whether the data are consistent before sorting and if so, don't
  // do the sort.  Concsistent means that all blocks that have a
  // nonzero element count are in the correct order.
  int64_t last_offset = 0;
  bool order_ok = true;
  for (size_t i=0; i < blocks.size(); i++) {
    if (blocks[i].elementCount > 0) {
      if (blocks[i].offset_ < last_offset) {
	order_ok = false;
	break;
      }
      last_offset = blocks[i].offset_;
    }
  }
  
  std::vector<Block> sorted_blocks(blocks);
  if (!order_ok) {
    std::sort(sorted_blocks.begin(), sorted_blocks.end(), lessOffset);
    
    // Now, update the position_ field based on the sorted order.
    for (size_t i=0; i < blocks.size(); i++) {
      int orig_position = sorted_blocks[i].position_;
      const_cast<Excn::Block&>(blocks[orig_position]).position_ = i;
      SMART_ASSERT(blocks[orig_position].id == sorted_blocks[i].id);
    }
  }

  int ierr;
  {
    Excn::Redefine the_database(exodusFilePtr);

    ierr=put_metadata(mesh, comm);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(sorted_blocks);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(nodesets);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(sidesets);
    if (ierr != EX_NOERR) return(ierr);
  }

  // NON-Define mode output...
  ierr=put_non_define_data(mesh, comm);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(sorted_blocks);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(nodesets);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(sidesets);
  if (ierr != EX_NOERR) return(ierr);

  // For now, put entity names using the ExodusII api...
  {
    ssize_t max_entity = mesh.blockCount;
    if (mesh.nodesetCount > max_entity)
      max_entity = mesh.nodesetCount;
    if (mesh.sidesetCount > max_entity)
      max_entity = mesh.sidesetCount;
    if (mesh.blockCount > 0) {
      for (int i=0; i < mesh.blockCount; i++) {
	if (blocks[i].attributeCount > max_entity)
	  max_entity = blocks[i].attributeCount;
      }
    }
  	 

    int name_size = ex_inquire_int(exodusFilePtr, EX_INQ_MAX_READ_NAME_LENGTH);
    char **names = new char* [max_entity];
    for (ssize_t i=0; i < max_entity; i++) {
      names[i] = new char [name_size+1];
    }
  	 
    if (mesh.blockCount > 0) {
      for (int i=0; i < mesh.blockCount; i++) {
	std::strncpy(names[i], blocks[i].name_.c_str(), name_size);
	names[i][name_size] = 0;
      }
      ex_put_names(exodusFilePtr, EX_ELEM_BLOCK, names);

      for (int i=0; i < mesh.blockCount; i++) {
	if (blocks[i].attributeCount > 0) {
	  SMART_ASSERT((size_t)blocks[i].attributeCount == blocks[i].attributeNames.size());
	  for (int j=0; j < blocks[i].attributeCount; j++) {
	    std::memset(names[j], '\0', name_size+1);
	    if (blocks[i].attributeNames[j].size() > 0) {
	      std::strncpy(names[j], blocks[i].attributeNames[j].c_str(),
			   name_size);
	      names[j][name_size] = 0;
	    }
	  }
	  ierr=ex_put_attr_names(exodusFilePtr, EX_ELEM_BLOCK, blocks[i].id, names);
	  SMART_ASSERT(ierr == 0);
	}
      }
    }
  	 
    if (mesh.nodesetCount > 0) {
      for (int i=0; i < mesh.nodesetCount; i++) {
	std::strncpy(names[i], nodesets[i].name_.c_str(), name_size);
	names[i][name_size] = 0;
      }
      ex_put_names(exodusFilePtr, EX_NODE_SET, names);
    }

    if (mesh.sidesetCount > 0) {
      for (int i=0; i < mesh.sidesetCount; i++) {
	std::strncpy(names[i], sidesets[i].name_.c_str(), name_size);
	names[i][name_size] = 0;
      }
      ex_put_names(exodusFilePtr, EX_SIDE_SET, names);
    }

    for (ssize_t i=0; i < max_entity; i++) {
      delete [] names[i];
    }
    delete [] names;
  }

  ex_update(exodusFilePtr);
  return(EX_NOERR);
}

template <typename INT>
bool Excn::Internals<INT>::check_meta_data(const Mesh &mesh,
					   const std::vector<Block>   &blocks,
					   const std::vector<NodeSet<INT> > &nodesets,
					   const std::vector<SideSet<INT> > &sidesets,
					   const CommunicationMetaData &comm)
{
  ex_init_params init_data;
  ex_get_init_ext(exodusFilePtr, &init_data);

  bool matches = true;
  if (mesh.dimensionality != init_data.num_dim) {
    std::cerr << "ERROR: original mesh dimensionality ("
	      << mesh.dimensionality << ") does not match current dimensionality ("
	      << init_data.num_dim << ")\n";
    matches = false;
  }

  if (mesh.nodeCount != init_data.num_nodes) {
    std::cerr << "ERROR: original mesh node count ("
	      << mesh.nodeCount << ") does not match current node count ("
	      << init_data.num_nodes << ")\n";
    matches = false;
  }

  if (mesh.elementCount != init_data.num_elem) {
    std::cerr << "ERROR: original mesh element count ("
	      << mesh.elementCount << ") does not match current element count ("
	      << init_data.num_elem << ")\n";
    matches = false;
  }

  if (mesh.blockCount != init_data.num_elem_blk) {
    std::cerr << "ERROR: original mesh element block count ("
	      << mesh.blockCount << ") does not match current element block count ("
	      << init_data.num_elem_blk << ")\n";
    matches = false;
  }

  if (mesh.nodesetCount != init_data.num_node_sets) {
    std::cerr << "ERROR: original mesh nodeset count ("
	      << mesh.nodesetCount << ") does not match current nodeset count ("
	      << init_data.num_node_sets << ")\n";
    matches = false;
  }

  if (mesh.sidesetCount != init_data.num_side_sets) {
    std::cerr << "ERROR: original mesh sideset count ("
	      << mesh.sidesetCount << ") does not match current sideset count ("
	      << init_data.num_side_sets << ")\n";
    matches = false;
  }

  return matches;
}

template <typename INT>
int Excn::Internals<INT>::put_metadata(const Mesh &mesh,
				       const CommunicationMetaData&)
{
  int numdimdim  = 0;
  int numnoddim  = 0;
  int strdim     = 0;
  int namestrdim = 0;

  int map_type  = get_type(exodusFilePtr, EX_MAPS_INT64_DB);

  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Excn::Internals::put_metadata()";

  // define some attributes...
  int status = nc_put_att_text(exodusFilePtr, NC_GLOBAL, ATT_TITLE, 
			       strlen(mesh.title.c_str())+1, mesh.title.c_str());
  if (status != NC_NOERR) {
    sprintf(errmsg,
	    "Error: failed to define title attribute to file id %d", exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL); 
  }

  // For use later to help readers know how much memory to allocate
  // for name storage, we define an attribute containing the maximum
  // size of any name.
  {
    int current_len = 0;
    status=nc_put_att_int(exodusFilePtr, NC_GLOBAL, ATT_MAX_NAME_LENGTH, NC_INT, 1, &current_len);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
              "Error: failed to define ATT_MAX_NAME_LENGTH attribute to file id %d", exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }
  }

  // inquire previously defined dimensions
  status=nc_inq_dimid (exodusFilePtr, DIM_STR, &strdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get string length in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  // create name string length dimension
  if (maximumNameLength < 32)
    maximumNameLength = 32;
  status = nc_def_dim (exodusFilePtr, DIM_STR_NAME, maximumNameLength+1, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to define name string length in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  // ...and some dimensions..
  status=nc_def_dim(exodusFilePtr, DIM_NUM_DIM, mesh.dimensionality, &numdimdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to define number of dimensions in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  if (mesh.nodeCount > 0) {
    status=nc_def_dim(exodusFilePtr, DIM_NUM_NODES, mesh.nodeCount, &numnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of nodes in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // Define the node map here to avoid a later redefine call
    if (mesh.needNodeMap) {
      int dims[1];
      int varid = 0;
      dims[0] = numnoddim;
      status=nc_def_var(exodusFilePtr, VAR_NODE_NUM_MAP, map_type, 1, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	if (status == NC_ENAMEINUSE) {
	  sprintf(errmsg,
		  "Error: node numbering map already exists in file id %d",
		  exodusFilePtr);
	  ex_err(routine, errmsg, status);
	} else {
	  sprintf(errmsg,
		  "Error: failed to create node numbering map array in file id %d",
		  exodusFilePtr);
	  ex_err(routine, errmsg, status);
	}
	return (EX_FATAL);
      }
      ex_compress_variable(exodusFilePtr, varid, 1);
    }
  }

  if (mesh.elementCount > 0) {
    int numelemdim = 0;
    status=nc_def_dim(exodusFilePtr, DIM_NUM_ELEM, mesh.elementCount, &numelemdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of elements in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // Define the element map here to avoid a later redefine call
    if (mesh.needElementMap) {
      int dims[1];
      dims[0] = numelemdim;
      int varid = 0;
      status=nc_def_var(exodusFilePtr, VAR_ELEM_NUM_MAP, map_type, 1, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	if (status == NC_ENAMEINUSE) {
	  sprintf(errmsg,
		  "Error: element numbering map already exists in file id %d",
		  exodusFilePtr);
	  ex_err(routine, errmsg, status);
	} else {
	  sprintf(errmsg,
		  "Error: failed to create element numbering map in file id %d",
		  exodusFilePtr);
	  ex_err(routine, errmsg, status);
	}
	return (EX_FATAL);
      }
      ex_compress_variable(exodusFilePtr, varid, 1);
    }
  }

  if (mesh.blockCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "element block", mesh.blockCount,
				DIM_NUM_EL_BLK, VAR_STAT_EL_BLK, VAR_ID_EL_BLK, VAR_NAME_EL_BLK);
  }

  // node set id array:
  if (mesh.nodesetCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "node set", mesh.nodesetCount,
				DIM_NUM_NS, VAR_NS_STAT, VAR_NS_IDS, VAR_NAME_NS);
  }

  // side set id array:
  if (mesh.sidesetCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "side set", mesh.sidesetCount,
				DIM_NUM_SS, VAR_SS_STAT, VAR_SS_IDS, VAR_NAME_SS);
  }

  status = define_coordinate_vars(exodusFilePtr, mesh.nodeCount, numnoddim, mesh.dimensionality, numdimdim, namestrdim);
  if (status != EX_NOERR)
    return EX_FATAL;
  
  return (EX_NOERR);
}

template <typename INT>
int Excn::Internals<INT>::put_metadata(const std::vector<Block> &blocks)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::put_metadata(blocks)";
  int dims[2];

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status  = 0; // clear error code

  if (blocks.size() == 0)
    return (EX_NOERR);

  // Get number of element blocks defined for this file
  int dimid;
  size_t num_elem_blk = 0;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_EL_BLK, &dimid);
  if (status != NC_NOERR ) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: no element blocks defined in file id %d",
	    exodusFilePtr);
    ex_err(routine, errmsg, status);
    return (EX_FATAL);
  }

  int namestrdim;
  status=nc_inq_dimid (exodusFilePtr, DIM_STR_NAME, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get name string length in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  status=nc_inq_dimlen(exodusFilePtr,dimid,&num_elem_blk);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get number of element blocks in file id %d",
	    exodusFilePtr);
    ex_err(routine, errmsg, status);
    return (EX_FATAL);
  }

  SMART_ASSERT(blocks.size() == num_elem_blk);

  // Iterate over element blocks ...
  for (size_t iblk = 0; iblk < num_elem_blk; iblk++) {

    ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_ELEM_BLOCK));

    if (blocks[iblk].elementCount == 0)
      continue;

    // define some dimensions and variables
    int numelbdim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_EL_IN_BLK(iblk+1),
		       blocks[iblk].elementCount, &numelbdim);
    if (status != NC_NOERR) {
      if (status == NC_ENAMEINUSE) {	// duplicate entry
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: element block %"PRId64" already defined in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of elements/block for block %"PRId64" file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return (EX_FATAL);
    }

    int nelnoddim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_NOD_PER_EL(iblk+1),
		       blocks[iblk].nodesPerElement, &nelnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of nodes/element for block %"PRId64" in file id %d",
	      blocks[iblk].id, exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    // element attribute array
    if (blocks[iblk].attributeCount > 0) {
      int numattrdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_BLK(iblk+1),
			 blocks[iblk].attributeCount, &numattrdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of attributes in block %"PRId64" in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      dims[0] = numelbdim;
      dims[1] = numattrdim;
      int varid = 0;
      status=nc_def_var(exodusFilePtr, VAR_ATTRIB(iblk+1),
			nc_flt_code(exodusFilePtr), 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error:  failed to define attributes for element block %"PRId64" in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
      ex_compress_variable(exodusFilePtr, varid, 2);

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_ATTRIB(iblk+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for element block %"PRId64" in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // element connectivity array
    dims[0] = numelbdim;
    dims[1] = nelnoddim;

    int connid;
    status=nc_def_var(exodusFilePtr, VAR_CONN(iblk+1), bulk_type, 2, dims, &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to create connectivity array for block %"PRId64" in file id %d",
	      blocks[iblk].id, exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }
    ex_compress_variable(exodusFilePtr, connid, 1);

    // store element type as attribute of connectivity variable
    status=nc_put_att_text(exodusFilePtr, connid, ATT_NAME_ELB,
			   (int)std::strlen(blocks[iblk].elType)+1, blocks[iblk].elType);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to store element type name %s in file id %d",
	      blocks[iblk].elType,exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }
  }
  return (EX_NOERR);
}


template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const Mesh&,
				   const CommunicationMetaData&)
{
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<Block> &blocks)
{
  int num_elem_blk = blocks.size();  // Verified via SMART_ASSERT earlier...
  
  if (num_elem_blk > 0) {
    // first get id of element block ids array variable
    std::vector<ex_entity_id> elem_blk_id(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_id[iblk] = blocks[iblk].id;
    }

    if (put_id_array(exodusFilePtr, VAR_ID_EL_BLK, elem_blk_id) != NC_NOERR)
      return(EX_FATAL);

    // Now, write the element block status array
    std::vector<int> elem_blk_status(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_status[iblk] = blocks[iblk].elementCount > 0 ? 1 : 0;
    }

    if (put_array(exodusFilePtr, VAR_STAT_EL_BLK, elem_blk_status) != NC_NOERR)
      return(EX_FATAL);
  }
  return (EX_NOERR);
}

// ========================================================================
template <typename INT>
int Excn::Internals<INT>::put_metadata(const std::vector<NodeSet<INT> > &nodesets)
{
  const char *routine = "Excn::Internals::put_metadata(nodesets)";
  if (nodesets.empty())
    return EX_NOERR;

  char errmsg[MAX_ERR_LENGTH];
  int dims[2];

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status  = 0; // clear error code

  // Get number of node sets defined for this file
  int dimid;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_NS, &dimid);
  if (status  != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      sprintf(errmsg,
	      "Error: no node sets defined for file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    } else {
      sprintf(errmsg,
	      "Error: failed to locate node sets defined in file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    }
    return (EX_FATAL);
  }

  // inquire how many node sets are to be stored
  int num_node_sets = ex_inquire_int(exodusFilePtr, EX_INQ_NODE_SETS);
  SMART_ASSERT(static_cast<int>(nodesets.size()) == num_node_sets);

  for (int i=0; i<num_node_sets; i++) {

    //  NOTE: ex_inc_file_item is used to find the number of node sets
    // for a specific file and returns that value incremented.
    int cur_num_node_sets=(int)ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_NODE_SET));

    if (nodesets[i].nodeCount == 0)
      continue;

    status=nc_def_dim (exodusFilePtr, DIM_NUM_NOD_NS(cur_num_node_sets+1),
		       nodesets[i].nodeCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: node set %"PRId64" already defined in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to define number of nodes for set %"PRId64" in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }


    // define variable to store node set node list here instead of in expns
    dims[0] = dimid;
    int varid = 0;
    status=nc_def_var(exodusFilePtr,
		      VAR_NODE_NS(cur_num_node_sets+1),bulk_type,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: node set %"PRId64" node list already defined in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create node set %"PRId64" node list in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }
    ex_compress_variable(exodusFilePtr, varid, 1);

    // Create variable for distribution factors if required
    if (nodesets[i].dfCount > 0) {
      // num_dist_per_set should equal num_nodes_per_set
      if (nodesets[i].dfCount != nodesets[i].nodeCount) {
	status = EX_FATAL;
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: # dist fact (%"PRId64") not equal to # nodes (%"PRId64") "
		"in node set %"PRId64" file id %d",
		nodesets[i].dfCount, nodesets[i].nodeCount, nodesets[i].id,
		exodusFilePtr);
	ex_err(routine, errmsg, status);
	return(EX_FATAL);
      } else {
	// create variable for distribution factors
	status=nc_def_var(exodusFilePtr, VAR_FACT_NS(cur_num_node_sets+1),
			  nc_flt_code(exodusFilePtr), 1, dims, &varid);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  if (status == NC_ENAMEINUSE) {
	    sprintf(errmsg,
		    "Error: node set %"PRId64" dist factors already exist in file id %d",
		    nodesets[i].id,exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  } else {
	    sprintf(errmsg,
		    "Error: failed to create node set %"PRId64" dist factors in file id %d",
		    nodesets[i].id, exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  }
	  return(EX_FATAL);
	}
	ex_compress_variable(exodusFilePtr, varid, 2);
      }
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<NodeSet<INT> > &nodesets)
{
  if (nodesets.empty())
    return EX_NOERR;

  // Output nodeset ids...
  size_t num_nodesets = nodesets.size();
  std::vector<ex_entity_id> nodeset_id(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    nodeset_id[i] = nodesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_NS_IDS, nodeset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    status[i] = nodesets[i].nodeCount > 0 ? 1 : 0;
  }

  if (put_array(exodusFilePtr, VAR_NS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

// ========================================================================
template <typename INT>
int Excn::Internals<INT>::put_metadata(const std::vector<SideSet<INT> > &sidesets)
{
  const char *routine = "Excn::Internals::put_metadata(sidesets)";
  if (sidesets.empty())
    return EX_NOERR;

  char errmsg[MAX_ERR_LENGTH];
  int dims[2];

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status  = 0; // clear error code

  // Get number of side sets defined for this file
  int dimid;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_SS, &dimid);
  if (status  != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      sprintf(errmsg,
	      "Error: no side sets defined for file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    } else {
      sprintf(errmsg,
	      "Error: failed to locate side sets defined in file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    }
    return (EX_FATAL);
  }

  size_t num_side_sets = ex_inquire_int(exodusFilePtr, EX_INQ_SIDE_SETS);
  SMART_ASSERT(sidesets.size() == num_side_sets);

  for (size_t i=0; i<num_side_sets; i++) {

    //  NOTE: ex_inc_file_item is used to find the number of side sets
    // for a specific file and returns that value incremented.
    int cur_num_side_sets = (int)ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_SIDE_SET));

    if (sidesets[i].sideCount == 0)
      continue;

    status=nc_def_dim (exodusFilePtr, DIM_NUM_SIDE_SS(cur_num_side_sets+1),
		       sidesets[i].sideCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: side set %"PRId64" already defined in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to define number of sides for set %"PRId64" in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }


    dims[0] = dimid;
    int varid = 0;
    status=nc_def_var(exodusFilePtr, VAR_ELEM_SS(cur_num_side_sets+1),
		      bulk_type,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: side set %"PRId64" element list already defined in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create side set %"PRId64" element list in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }
    ex_compress_variable(exodusFilePtr, varid, 1);

    // create side list variable for side set
    status=nc_def_var(exodusFilePtr, VAR_SIDE_SS(cur_num_side_sets+1),
		      bulk_type,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: side list already exists for side set %"PRId64" in file id %d",
		sidesets[i].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create side list for side set %"PRId64" in file id %d",
		sidesets[i].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }
    ex_compress_variable(exodusFilePtr, varid, 1);

    // Create variable for distribution factors if required
    if (sidesets[i].dfCount > 0) {
      status=nc_def_dim(exodusFilePtr, DIM_NUM_DF_SS(cur_num_side_sets+1),
			sidesets[i].dfCount, &dimid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	if (status == NC_ENAMEINUSE) {
	  sprintf(errmsg,
		  "Error: side set df count %"PRId64" already defined in file id %d",
		  sidesets[i].id, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	} else {
	  sprintf(errmsg,
		  "Error: failed to define side set df count for set %"PRId64" in file id %d",
		  sidesets[i].id, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	}
	return(EX_FATAL);
      }

      // create distribution factor list variable for side set
      dims[0] = dimid;
      status=nc_def_var(exodusFilePtr, VAR_FACT_SS(cur_num_side_sets+1),
			nc_flt_code(exodusFilePtr), 1, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	if (status == NC_ENAMEINUSE) {
	  sprintf(errmsg,
		  "Error: dist factor list already exists for side set %"PRId64" in file id %d",
		  sidesets[i].id, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	} else {
	  sprintf(errmsg,
		  "Error: failed to create dist factor list for side set %"PRId64" in file id %d",
		  sidesets[i].id, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	}
	return(EX_FATAL);
      }
      ex_compress_variable(exodusFilePtr, varid, 2);
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<SideSet<INT> > &sidesets)
{
  if (sidesets.empty())
    return EX_NOERR;

  // Output sideset ids...
  int num_sidesets = (int)sidesets.size();
  std::vector<ex_entity_id> sideset_id(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    sideset_id[i] = sidesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_SS_IDS, sideset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    status[i] = sidesets[i].sideCount > 0 ? 1 : 0;
  }

  if (put_array(exodusFilePtr, VAR_SS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

namespace {
  template <typename INT>
  int put_array(int exoid, const char *var_type, const std::vector<INT> &array)
  {
    const char *routine = "Internals.C, put_int_array";
    char errmsg[MAX_ERR_LENGTH];
    int var_id;
    int status;

    status=nc_inq_varid(exoid, var_type, &var_id);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to locate %s in file id %d", var_type, exoid);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    if (sizeof(INT) == sizeof(int64_t))
      status=nc_put_var_longlong(exoid, var_id, (long long int*)&array[0]);
    else
      status=nc_put_var_int(exoid, var_id, &array[0]);
    
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to write %s array in file id %d", var_type, exoid);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }
    return (EX_NOERR);
  }

  int put_id_array(int exoid, const char *var_type, const std::vector<ex_entity_id> &ids)
  {
    const char *routine = "Internals.C, put_id_array";
    char errmsg[MAX_ERR_LENGTH];
    int var_id;

    int status=nc_inq_varid(exoid, var_type, &var_id);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to locate %s in file id %d", var_type, exoid);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    int id_type = get_type(exoid, EX_IDS_INT64_API);

    if (id_type == NC_INT64) {
      status=nc_put_var_longlong(exoid, var_id, (long long int*)&ids[0]);
    } else {
      // Have ex_entity_id (long long), need ints...
      std::vector<int> int_ids(ids.size());
      int_ids.assign(ids.begin(), ids.end());
      status=nc_put_var_int(exoid, var_id, &int_ids[0]);
    }

    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to write %s array in file id %d", var_type, exoid);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }
    return (EX_NOERR);
  }
  
  int define_coordinate_vars(int exodusFilePtr, int64_t nodes, int node_dim, int dimension, int dim_dim,
			     int str_dim)
  {
    const char *routine = "Ioex_Internals.C, define_coordinate_vars";
    char errmsg[MAX_ERR_LENGTH];
    int status;
    int dim[2];
    int varid;
    
    if (nodes > 0) {
      if (ex_large_model(exodusFilePtr) == 1) {
	// node coordinate arrays -- separate storage...
	dim[0] = node_dim;
	if (dimension > 0) {
	  status=nc_def_var(exodusFilePtr, VAR_COORD_X, nc_flt_code(exodusFilePtr), 1, dim, &varid);
	  if (status != NC_NOERR) {
	    ex_opts(EX_VERBOSE);
	    sprintf(errmsg,
		    "Error: failed to define node x coordinate array in file id %d",exodusFilePtr);
	    ex_err(routine,errmsg,status);
	    return(EX_FATAL);
	  }
	  ex_compress_variable(exodusFilePtr, varid, 2);
	}

	if (dimension > 1) {
	  status=nc_def_var(exodusFilePtr, VAR_COORD_Y, nc_flt_code(exodusFilePtr), 1, dim, &varid);
	  if (status != NC_NOERR) {
	    ex_opts(EX_VERBOSE);
	    sprintf(errmsg,
		    "Error: failed to define node y coordinate array in file id %d",exodusFilePtr);
	    ex_err(routine,errmsg,status);
	    return(EX_FATAL);
	  }
	  ex_compress_variable(exodusFilePtr, varid, 2);
	}

	if (dimension > 2) {
	  status=nc_def_var(exodusFilePtr, VAR_COORD_Z, nc_flt_code(exodusFilePtr), 1, dim, &varid);
	  if (status != NC_NOERR) {
	    ex_opts(EX_VERBOSE);
	    sprintf(errmsg,
		    "Error: failed to define node z coordinate array in file id %d",exodusFilePtr);
	    ex_err(routine,errmsg,status);
	    return(EX_FATAL);
	  }
	  ex_compress_variable(exodusFilePtr, varid, 2);
	}
      } else {
	// node coordinate arrays:  -- all stored together (old method)2
	dim[0] = dim_dim;
	dim[1] = node_dim;
	status=nc_def_var(exodusFilePtr, VAR_COORD, nc_flt_code(exodusFilePtr), 2, dim, &varid);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  sprintf(errmsg,
		  "Error: failed to define node coordinate array in file id %d",exodusFilePtr);
	  ex_err(routine,errmsg,status);
	  return(EX_FATAL);
	}
	ex_compress_variable(exodusFilePtr, varid, 2);
      }
    }

    // coordinate names array
    dim[0] = dim_dim;
    dim[1] = str_dim;

    status=nc_def_var(exodusFilePtr, VAR_NAME_COOR, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define coordinate name array in file id %d",
	      exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }
    return EX_NOERR;
  }

  int define_netcdf_vars(int exoid, const char *type, size_t count, const char *dim_num,
			 const char *stat_var, const char *id_var, const char *name_var)
  {
    int dimid = 0;
    int varid = 0;
    int dim[2];
    int namestrdim = 0;
    char errmsg[MAX_ERR_LENGTH];
    const char *routine = "Internals::define_netcdf_vars()";

    int status=nc_inq_dimid (exoid, DIM_STR_NAME, &namestrdim);

    status=nc_def_dim(exoid, dim_num, count, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of %ss in file id %d", type, exoid);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // id status array:
    dim[0] = dimid;
    status=nc_def_var(exoid, stat_var, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define side %s status in file id %d", type, exoid);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // id array:
    int ids_type  = get_type(exoid, EX_IDS_INT64_DB);
    status=nc_def_var(exoid, id_var, ids_type, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define %s property in file id %d", type, exoid);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // store property name as attribute of property array variable
    status=nc_put_att_text(exoid, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to store %s property name %s in file id %d",
	      type, "ID",exoid);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }

    dim[0] = dimid;
    dim[1] = namestrdim;

    status=nc_def_var(exoid, name_var, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define %s name array in file id %d",
	      type, exoid);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }
    return EX_NOERR;
  }

}

