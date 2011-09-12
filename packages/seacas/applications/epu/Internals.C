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
#include <Internals.h>

#include <ExodusEntity.h>

#include <algorithm>
#include <cstring>
#include <smart_assert.h>

extern "C" {
#define NO_NETCDF_2
#include <netcdf.h>
#include <exodusII.h>
#include <exodusII_int.h>
}

namespace {
  bool lessOffset(const Excn::Block& b1, const Excn::Block& b2) {
    return b1.offset_ < b2.offset_;
  }

  int put_int_array(int exoid, const char *var_type, const std::vector<int> &array);

  int define_coordinate_vars(int exodusFilePtr, int nodes, int node_dim,
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

Excn::Internals::Internals(int exoid, int maximum_name_length)
  : exodusFilePtr(exoid),
    nodeMapVarID(),
    elementMapVarID(),
    commIndexVar(0),
    elemCommIndexVar(0),
    maximumNameLength(maximum_name_length)
{}

bool Excn::Internals::check_meta_data(const Mesh &mesh,
				      const std::vector<Block>   &blocks,
				      const std::vector<NodeSet> &nodesets,
				      const std::vector<SideSet> &sidesets,
				      const CommunicationMetaData &comm)
{
  char* mytitle = new char[MAX_LINE_LENGTH + 1];
  int dimensionality;
  int node_count;
  int element_count;
  int block_count;
  int nodeset_count;
  int sideset_count;

  ex_get_init(exodusFilePtr, mytitle,  &dimensionality, &node_count, &element_count,
	      &block_count, &nodeset_count, &sideset_count);

  bool matches = true;
  if (mesh.dimensionality != dimensionality) {
    std::cerr << "ERROR: original mesh dimensionality ("
	      << mesh.dimensionality << ") does not match current dimensionality ("
	      << dimensionality << ")\n";
    matches = false;
  }

  if (mesh.nodeCount != node_count) {
    std::cerr << "ERROR: original mesh node count ("
	      << mesh.nodeCount << ") does not match current node count ("
	      << node_count << ")\n";
    matches = false;
  }

  if (mesh.elementCount != element_count) {
    std::cerr << "ERROR: original mesh element count ("
	      << mesh.elementCount << ") does not match current element count ("
	      << element_count << ")\n";
    matches = false;
  }

  if (mesh.blockCount != block_count) {
    std::cerr << "ERROR: original mesh block count ("
	      << mesh.blockCount << ") does not match current block count ("
	      << block_count << ")\n";
    matches = false;
  }

  if (mesh.nodesetCount != nodeset_count) {
    std::cerr << "ERROR: original mesh nodeset count ("
	      << mesh.nodesetCount << ") does not match current nodeset count ("
	      << nodeset_count << ")\n";
    matches = false;
  }

  if (mesh.sidesetCount != sideset_count) {
    std::cerr << "ERROR: original mesh sideset count ("
	      << mesh.sidesetCount << ") does not match current sideset count ("
	      << sideset_count << ")\n";
    matches = false;
  }

  return matches;
}

int Excn::Internals::write_meta_data(const Mesh &mesh,
				     const std::vector<Block>   &blocks,
				     const std::vector<NodeSet> &nodesets,
				     const std::vector<SideSet> &sidesets,
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
  int last_offset = 0;
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
    int max_entity = mesh.blockCount;
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
    for (int i=0; i < max_entity; i++) {
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

    for (int i=0; i < max_entity; i++) {
      delete [] names[i];
    }
    delete [] names;
  }

  ex_update(exodusFilePtr);
  return(EX_NOERR);
}

int Excn::Internals::put_metadata(const Mesh &mesh,
				  const CommunicationMetaData&)
{
  int numdimdim  = 0;
  int numnoddim  = 0;
  int numelemdim = 0;
  int elblkdim   = 0;
  int strdim     = 0;
  int namestrdim = 0;
  int dim[2];
  int dimid = 0;
  int varid = 0;

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
    int dims[1];
    dims[0] = numnoddim;
    status=nc_def_var(exodusFilePtr, VAR_NODE_NUM_MAP, NC_INT, 1, dims, &varid);
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
  }

  if (mesh.elementCount > 0) {
    status=nc_def_dim(exodusFilePtr, DIM_NUM_ELEM, mesh.elementCount, &numelemdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of elements in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // Define the element map here to avoid a later redefine call
    int dims[1];
    dims[0] = numelemdim;
    varid = 0;
    status=nc_def_var(exodusFilePtr, VAR_ELEM_NUM_MAP, NC_INT, 1, dims, &varid);
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
  }

  if (mesh.blockCount > 0) {
    status=nc_def_dim (exodusFilePtr, DIM_NUM_EL_BLK, mesh.blockCount, &elblkdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of element blocks in file id %d",
	      exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // ...and some variables
    // element block id status array...
    dim[0] = elblkdim;
    status=nc_def_var(exodusFilePtr, VAR_STAT_EL_BLK, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define element block status array in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // element block id array
    status=nc_def_var(exodusFilePtr, VAR_ID_EL_BLK, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define element block id array in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // store property name as attribute of property array variable
    status=nc_put_att_text(exodusFilePtr, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to store element block property name %s in file id %d",
	      "ID",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }

    dim[0] = elblkdim;
    dim[1] = namestrdim;

    status = nc_def_var (exodusFilePtr, VAR_NAME_EL_BLK, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define element block name array in file id %d",
	      exodusFilePtr);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }
  }

  // node set id array:
  if (mesh.nodesetCount > 0) {

    status=nc_def_dim (exodusFilePtr, DIM_NUM_NS, mesh.nodesetCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of node sets in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // node set id status array:
    dim[0] = dimid;
    status=nc_def_var(exodusFilePtr, VAR_NS_STAT, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to create node sets status array in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }


    // node set id array:
    dim[0] = dimid;
    status=nc_def_var(exodusFilePtr, VAR_NS_IDS, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to create node sets property array in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }


    // store property name as attribute of property array variable
    status=nc_put_att_text(exodusFilePtr, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to store node set property name %s in file id %d",
	      "ID",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }

    dim[0] = dimid;
    dim[1] = namestrdim;

    status=nc_def_var(exodusFilePtr, VAR_NAME_NS, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define nodeset name array in file id %d",
	      exodusFilePtr);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }
  }

  // side set id array:
  if (mesh.sidesetCount > 0) {

    status=nc_def_dim(exodusFilePtr, DIM_NUM_SS, mesh.sidesetCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of side sets in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // side set id status array:
    dim[0] = dimid;
    status=nc_def_var(exodusFilePtr, VAR_SS_STAT, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define side set status in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // side set id array:
    status=nc_def_var(exodusFilePtr, VAR_SS_IDS, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define side set property in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // store property name as attribute of property array variable
    status=nc_put_att_text(exodusFilePtr, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to store side set property name %s in file id %d",
	      "ID",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }

    dim[0] = dimid;
    dim[1] = namestrdim;

    status=nc_def_var(exodusFilePtr, VAR_NAME_SS, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define sideset name array in file id %d",
	      exodusFilePtr);
      ex_err(routine,errmsg,status);
      return (EX_FATAL);
    }
  }

  status = define_coordinate_vars(exodusFilePtr, mesh.nodeCount, numnoddim, mesh.dimensionality, numdimdim, namestrdim);
  if (status != EX_NOERR)
    return EX_FATAL;
  
  return (EX_NOERR);
}

int Excn::Internals::put_metadata(const std::vector<Block> &blocks)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::put_metadata(blocks)";
  int dims[2];

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
		"Error: element block %d already defined in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of elements/block for block %d file id %d",
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
	      "Error: failed to define number of nodes/element for block %d in file id %d",
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
		"Error: failed to define number of attributes in block %d in file id %d",
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
		"Error:  failed to define attributes for element block %d in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_ATTRIB(iblk+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for element block %d in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // element connectivity array
    dims[0] = numelbdim;
    dims[1] = nelnoddim;

    int connid;
    status=nc_def_var(exodusFilePtr, VAR_CONN(iblk+1), NC_INT, 2, dims, &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to create connectivity array for block %d in file id %d",
	      blocks[iblk].id, exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

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


int Excn::Internals::put_non_define_data(const Mesh&,
				   const CommunicationMetaData&)
{
  return EX_NOERR;
}

int Excn::Internals::put_non_define_data(const std::vector<Block> &blocks)
{
  int num_elem_blk = blocks.size();  // Verified via SMART_ASSERT earlier...
  
  if (num_elem_blk > 0) {
    // first get id of element block ids array variable
    std::vector<int> elem_blk_id(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_id[iblk] = blocks[iblk].id;
    }

    if (put_int_array(exodusFilePtr, VAR_ID_EL_BLK, elem_blk_id) != NC_NOERR)
      return(EX_FATAL);

    // Now, write the element block status array
    std::vector<int> elem_blk_status(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_status[iblk] = blocks[iblk].elementCount > 0 ? 1 : 0;
    }

    if (put_int_array(exodusFilePtr, VAR_STAT_EL_BLK, elem_blk_status) != NC_NOERR)
      return(EX_FATAL);
  }
  return (EX_NOERR);
}

// ========================================================================
int Excn::Internals::put_metadata(const std::vector<NodeSet> &nodesets)
{
  const char *routine = "Excn::Internals::put_metadata(nodesets)";
  if (nodesets.empty())
    return EX_NOERR;

  char errmsg[MAX_ERR_LENGTH];
  int dims[2];

  int status  = 0; // clear error code
  char *cdum = 0;
  float fdum;

  // Get number of node sets defined for this file
  int dimid;
  int num_node_sets = 0;
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
  if (ex_inquire(exodusFilePtr, EX_INQ_NODE_SETS, &num_node_sets, &fdum, cdum) != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get number of node sets defined for file id %d",
	    exodusFilePtr);
    ex_err(routine, errmsg, status);
    return (EX_FATAL);
  }

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
		"Error: node set %d already defined in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to define number of nodes for set %d in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }


    // define variable to store node set node list here instead of in expns
    dims[0] = dimid;
    int varid;
    status=nc_def_var(exodusFilePtr,
		      VAR_NODE_NS(cur_num_node_sets+1),NC_INT,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: node set %d node list already defined in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create node set %d node list in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }

    // Create variable for distribution factors if required
    if (nodesets[i].dfCount > 0) {
      // num_dist_per_set should equal num_nodes_per_set
      if (nodesets[i].dfCount != nodesets[i].nodeCount) {
	status = EX_FATAL;
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: # dist fact (%d) not equal to # nodes (%d) "
		"in node set %d file id %d",
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
		    "Error: node set %d dist factors already exist in file id %d",
		    nodesets[i].id,exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  } else {
	    sprintf(errmsg,
		    "Error: failed to create node set %d dist factors in file id %d",
		    nodesets[i].id, exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  }
	  return(EX_FATAL);
	}
      }
    }
  }
  return EX_NOERR;
}

int Excn::Internals::put_non_define_data(const std::vector<NodeSet> &nodesets)
{
  if (nodesets.empty())
    return EX_NOERR;

  // Output nodeset ids...
  int num_nodesets = nodesets.size();
  std::vector<int> nodeset_id(num_nodesets);
  for (int i = 0; i < num_nodesets; i++) {
    nodeset_id[i] = nodesets[i].id;
  }

  if (put_int_array(exodusFilePtr, VAR_NS_IDS, nodeset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_nodesets);
  for (int i = 0; i < num_nodesets; i++) {
    status[i] = nodesets[i].nodeCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_NS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

// ========================================================================
int Excn::Internals::put_metadata(const std::vector<SideSet> &sidesets)
{
  const char *routine = "Excn::Internals::put_metadata(sidesets)";
  if (sidesets.empty())
    return EX_NOERR;

  char errmsg[MAX_ERR_LENGTH];
  int dims[2];

  int status  = 0; // clear error code
  char *cdum = 0;
  float fdum;

  // Get number of side sets defined for this file
  int dimid;
  int num_side_sets = 0;
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

  // inquire how many side sets are to be stored
  if (ex_inquire(exodusFilePtr, EX_INQ_SIDE_SETS, &num_side_sets, &fdum, cdum) != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get number of side sets defined for file id %d",
	    exodusFilePtr);
    ex_err(routine, errmsg, status);
    return (EX_FATAL);
  }

  SMART_ASSERT(static_cast<int>(sidesets.size()) == num_side_sets);

  for (int i=0; i<num_side_sets; i++) {

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
		"Error: side set %d already defined in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to define number of sides for set %d in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }


    dims[0] = dimid;
    int varid = 0;
    status=nc_def_var(exodusFilePtr, VAR_ELEM_SS(cur_num_side_sets+1),
		      NC_INT,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: side set %d element list already defined in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create side set %d element list in file id %d",
		sidesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }

    // create side list variable for side set
    status=nc_def_var(exodusFilePtr, VAR_SIDE_SS(cur_num_side_sets+1),
		      NC_INT,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: side list already exists for side set %d in file id %d",
		sidesets[i].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create side list for side set %d in file id %d",
		sidesets[i].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }

    // Create variable for distribution factors if required
    if (sidesets[i].dfCount > 0) {
      status=nc_def_dim(exodusFilePtr, DIM_NUM_DF_SS(cur_num_side_sets+1),
			sidesets[i].dfCount, &dimid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	if (status == NC_ENAMEINUSE) {
	  sprintf(errmsg,
		  "Error: side set df count %d already defined in file id %d",
		  sidesets[i].id, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	} else {
	  sprintf(errmsg,
		  "Error: failed to define side set df count for set %d in file id %d",
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
		  "Error: dist factor list already exists for side set %d in file id %d",
		  sidesets[i].id, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	} else {
	  sprintf(errmsg,
		  "Error: failed to create dist factor list for side set %d in file id %d",
		  sidesets[i].id, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	}
	return(EX_FATAL);
      }
    }
  }
  return EX_NOERR;
}

int Excn::Internals::put_non_define_data(const std::vector<SideSet> &sidesets)
{
  if (sidesets.empty())
    return EX_NOERR;

  // Output sideset ids...
  int num_sidesets = sidesets.size();
  std::vector<int> sideset_id(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    sideset_id[i] = sidesets[i].id;
  }

  if (put_int_array(exodusFilePtr, VAR_SS_IDS, sideset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    status[i] = sidesets[i].sideCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_SS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

namespace {
  int put_int_array(int exoid, const char *var_type, const std::vector<int> &array)
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

  int define_coordinate_vars(int exodusFilePtr, int nodes, int node_dim, int dimension, int dim_dim, int str_dim)
  {
    const char *routine = "Ioex_Internals.C, define_coordinate_vars";
    char errmsg[MAX_ERR_LENGTH];
    int status;
    int dim[2];
    int varid;
    
    if (nodes > 0) {
      if (ex_large_model(exodusFilePtr) == 1) {
	// node coordinate arrays -- separate storage...

	/*
	 * Check that storage required for coordinates  is less
	 * than 2GB which is maximum size permitted by netcdf
	 * (in large file mode). 1<<29 == max number of integer items.
	 */
	int shift = nc_flt_code(exodusFilePtr) == NC_DOUBLE ? 28 : 29;
	if (nodes  > (1<<shift)) {
	  status  = EX_BADPARAM;
	  sprintf(errmsg,
		  "Error: Size to store nodal coordinates exceeds 2GB in file id %d",
		  exodusFilePtr);
	  ex_err(routine,errmsg,status);
	  return (EX_FATAL);
	}

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
}

