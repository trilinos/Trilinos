// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// Copyright 2001, 2008, 2009 Sandia Corporation, Albuquerque, NM.

#include <exodusII/Ioex_Internals.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeSet.h>
#include <Ioss_SideBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_Utils.h>
#include <algorithm>
extern "C" {
#include <exodusII_int.h>
#include <ne_nemesisI_int.h>
#include <ne_nemesisI.h>

#include <netcdf.h>
}
#include <assert.h>

#include <string>


using namespace Ioex;

namespace {
  int define_variable(int exodusFilePtr, int size, const char *dim, const char *var, nc_type type);
  int define_variables(int exodusFilePtr, int size, const char *dim, const char *var[], nc_type type);
  int conditional_define_variable(int exodusFilePtr, const char *var, int dimid, int *varid,  nc_type type);

  size_t max_string_length() {return MAX_STR_LENGTH;}
  int put_int_array(int exoid, const char *var_type,
		    const std::vector<int> &array);
  int define_coordinate_vars(int exodusFilePtr, int nodes, int node_dim,
			     int dimension, int dim_dim, int str_dim);
}

#define stringify(X) #X
#define version_string(X) stringify(X)

void Internals::update_last_time_attribute(double value)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::update_last_time_attribute()";
  int status = 0; // clear error code

  status=nc_put_att_double(exodusFilePtr, NC_GLOBAL, "last_written_time",
			   NC_DOUBLE, 1, &value);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to define 'last_written_time' attribute to file id %d",
	    exodusFilePtr);
    ex_err(routine,errmsg,status);
  }
}

bool Internals::read_last_time_attribute(double *value)
{
  // Check whether the "last_written_time" attribute exists.  If it does,
  // return the value of the attribute in 'value' and return 'true'.
  // If not, don't change 'value' and return 'false'.
  bool found = false;

  int status;
  nc_type att_type = NC_NAT;
  size_t att_len = 0;
  status = nc_inq_att(exodusFilePtr, NC_GLOBAL, "last_written_time", &att_type, &att_len);
  if (status == NC_NOERR && att_type == NC_DOUBLE) {
    // Attribute exists on this database, read it...
    double tmp = 0.0;
    status = nc_get_att_double(exodusFilePtr, NC_GLOBAL, "last_written_time", &tmp);
    if (status == NC_NOERR) {
      *value = tmp;
      found = true;
    } else {
      char errmsg[MAX_ERR_LENGTH];
      const char *routine = "Internals::read_last_time_attribute()";
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to read last_written_time attribute from file id %d", exodusFilePtr);
      ex_err(routine,errmsg,status);
      found = false;
    }
  }
  return found;
}

bool Internals::check_processor_info(int processor_count, int processor_id)
{
  // A restart file may contain an attribute which contains
  // information about the processor count and current processor id
  // when the file was written.  This code checks whether that
  // information matches the current processor count and id.  If it
  // exists, but doesn't match, a warning message is printed.
  // Eventually, this will be used to determine whether certain
  // decomposition-related data in the file is valid or has been
  // invalidated by a join/re-spread to a different number of
  // processors.
  bool matches = true;

  int status;
  nc_type att_type = NC_NAT;
  size_t att_len = 0;
  status = nc_inq_att(exodusFilePtr, NC_GLOBAL, "processor_info", &att_type, &att_len);
  if (status == NC_NOERR && att_type == NC_INT) {
    // Attribute exists on this database, read it and check that the information
    // matches the current processor count and procesor id.
    int proc_info[2];
    status = nc_get_att_int(exodusFilePtr, NC_GLOBAL, "processor_info", proc_info);
    if (status == NC_NOERR) {
      if (proc_info[0] != processor_count && proc_info[0] > 1) {
	IOSS_WARNING << "Processor decomposition count in file (" << proc_info[0]
		     << ") does not match current processor count (" << processor_count
		     << ").\n";
	matches = false;
    }
      if (proc_info[1] != processor_id) {
	IOSS_WARNING << "This file was originally written on processor " << proc_info[1]
		     << ", but is now being read on processor " << processor_id
		     << ". This may cause problems if there is any processor-dependent data on the file.\n";
	matches = false;
      }
    } else {
      char errmsg[MAX_ERR_LENGTH];
      const char *routine = "Internals::check_processor_info()";
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to read processor info attribute from file id %d", exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }
  }
  return matches;
}


Redefine::Redefine(int exoid)
  : exodusFilePtr(exoid)
{
  // Enter define mode...
  int status = nc_redef (exodusFilePtr);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    char errmsg[MAX_ERR_LENGTH];
    sprintf(errmsg,
	    "Error: failed to put file id %d into define mode", exodusFilePtr);
    ex_err("Ioex::Redefine::Redefine()",errmsg,status);
    exit(EXIT_FAILURE);
  }
}

Redefine::~Redefine()
{
  try {
    int status = nc_enddef (exodusFilePtr);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      char errmsg[MAX_ERR_LENGTH];
      sprintf(errmsg,
	      "Error: failed to complete variable definitions in file id %d",exodusFilePtr);
      ex_err("Ioex::Redefine::~Redefine()",errmsg,status);
      exit(EXIT_FAILURE);
    }
  } catch (...) {
  }
}

Block::Block(const Ioss::ElementBlock &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  elementCount = other.get_property("entity_count").get_int();
  nodesPerElement = other.get_property("topology_node_count").get_int();
  attributeCount = other.get_property("attribute_count").get_int();
  offset_ = other.get_offset();
  std::string el_type = other.get_property("topology_type").get_string();
  if (other.property_exists("original_element_type")) {
    el_type = other.get_property("original_element_type").get_string();
  }

  std::strncpy(elType, el_type.c_str(), max_string_length());
  elType[max_string_length()] = 0;

  // Fixup an exodusII kluge.  For triangular elements, the same
  // name is used for 2D elements and 3D shell elements.  Convert
  // to unambiguous names for the IO Subsystem.  The 2D name
  // stays the same, the 3D name becomes 'trishell#'
  // Here, we need to map back to the 'triangle' name...
  if (std::strncmp(elType, "trishell", 8) == 0)
    std::strncpy(elType, "triangle", 8);
}

Block& Block::operator=(const Block& other)
{
  name = other.name;
  id = other.id;
  elementCount = other.elementCount;
  nodesPerElement = other.nodesPerElement;
  attributeCount = other.attributeCount;
  offset_ = other.offset_;
  std::strcpy(elType, other.elType);
  return *this;
}

bool Block::operator==(const Block& other) const
{
  return name == other.name && id == other.id &&
    elementCount == other.elementCount &&
    nodesPerElement == other.nodesPerElement &&
    attributeCount == other.attributeCount &&
    std::strcmp(elType, other.elType) == 0;
}

NodeSet::NodeSet(const Ioss::NodeSet &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  nodeCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
}

bool NodeSet::operator==(const NodeSet& other) const
{
  return id == other.id &&
    nodeCount == other.nodeCount &&
    dfCount == other.dfCount &&
    name == other.name;
}

SideSet::SideSet(const Ioss::SideBlock &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  sideCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
  std::string io_name = other.name();

  // KLUGE: universal_sideset has side dfCount...
  if (io_name == "universal_sideset")
    dfCount = sideCount;
}

SideSet::SideSet(const Ioss::SideSet &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  sideCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
  std::string io_name = other.name();

  // KLUGE: universal_sideset has side dfCount...
  if (io_name == "universal_sideset")
    dfCount = sideCount;
}

bool SideSet::operator==(const SideSet& other) const
{
  return id == other.id &&
    sideCount == other.sideCount &&
    dfCount == other.dfCount &&
    name == other.name;
}

bool CommunicationMap::operator==(const CommunicationMap& other) const
{
  return id == other.id &&
    entityCount == other.entityCount &&
    type == other.type;
}

Internals::Internals(int exoid, int maximum_name_length)
  : exodusFilePtr(exoid),
    nodeMapVarID(),
    elementMapVarID(),
    commIndexVar(0),
    elemCommIndexVar(0),
    maximumNameLength(maximum_name_length)
{}

int Internals::write_meta_data(const Mesh &mesh,
			       const std::vector<Block>   &blocks,
			       const std::vector<NodeSet> &nodesets,
			       const std::vector<SideSet> &sidesets,
			       const CommunicationMetaData &comm)
{
  assert(static_cast<int>(blocks.size())   == mesh.blockCount);
  assert(static_cast<int>(nodesets.size()) == mesh.nodesetCount);
  assert(static_cast<int>(sidesets.size()) == mesh.sidesetCount);

  int ierr;
  {
    Redefine the_database(exodusFilePtr);

    // Set the database to NOFILL mode.  Only writes values we want written...
    int old_fill = 0;

    ierr=nc_set_fill(exodusFilePtr, NC_NOFILL, &old_fill);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh, comm);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(blocks);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(nodesets);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(sidesets);
    if (ierr != EX_NOERR) return(ierr);
  }

  // NON-Define mode output...
  ierr=put_non_define_data(mesh, comm);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(blocks);
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

    char **names = new char* [max_entity];

    for (int i=0; i < mesh.blockCount; i++) {
      names[i] = (char*)blocks[i].name.c_str();
    }
    ex_put_names(exodusFilePtr, EX_ELEM_BLOCK, names);


    if (mesh.nodesetCount > 0) {
      for (int i=0; i < mesh.nodesetCount; i++) {
	names[i] = (char*)nodesets[i].name.c_str();
      }
      ex_put_names(exodusFilePtr, EX_NODE_SET, names);
    }

    if (mesh.sidesetCount > 0) {
      for (int i=0; i < mesh.sidesetCount; i++) {
	names[i] = (char*)sidesets[i].name.c_str();
      }
      ex_put_names(exodusFilePtr, EX_SIDE_SET, names);
    }

    delete [] names;
  }

  return(EX_NOERR);
}

int Internals::put_metadata(const Mesh &mesh,
			    const CommunicationMetaData &comm)
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
  const char *routine = "Internals::put_metadata()";

  int status = nc_put_att_text(exodusFilePtr, NC_GLOBAL, ATT_TITLE,
			       (int)std::strlen(mesh.title)+1, mesh.title);

  // define some attributes...
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to define title attribute to file id %d", exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  // For use later as a consistency check, define the number of processors and
  // the current processor id as an attribute of the file...
  {
    int ltempsv[2];
    ltempsv[0] = comm.processorCount;
    ltempsv[1] = comm.processorId;
    status=nc_put_att_int(exodusFilePtr, NC_GLOBAL, "processor_info", NC_INT, 2, ltempsv);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define processor info attribute to file id %d", exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }
  }

  // For use later to determine whether a timestep is corrupt, we define an attribute
  // containing the last written time...
  {
    double fake_time = -1.0;
    status=nc_put_att_double(exodusFilePtr, NC_GLOBAL, "last_written_time",
			     NC_DOUBLE, 1, &fake_time);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define 'last_written_time' attribute to file id %d", exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }
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

  /* create name string length dimension */
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
  if (status != EX_NOERR) return EX_FATAL;

  // Define dimension for the number of processors
  if (comm.processorCount > 1) {
    int procdim;
    status=nc_inq_dimid(exodusFilePtr, DIM_NUM_PROCS, &procdim);
    if (status != NC_NOERR) {
      int ltempsv = comm.processorCount;
      status=nc_def_dim(exodusFilePtr, DIM_NUM_PROCS, ltempsv, &procdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to dimension \"%s\" in file ID %d",
		DIM_NUM_PROCS, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // If this is a parallel file then the status vectors are size 1
    int dimid_npf;
    status= nc_inq_dimid(exodusFilePtr, DIM_NUM_PROCS_F, &dimid_npf);
    if ((status) != NC_NOERR) {
      int ltempsv = 1; // 1 processor per file...
      status=nc_def_dim(exodusFilePtr, DIM_NUM_PROCS_F, ltempsv, &dimid_npf);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to dimension \"%s\" in file ID %d",
		DIM_NUM_PROCS_F, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // Define the file type variable...
    status=nc_inq_varid(exodusFilePtr, VAR_FILE_TYPE, &varid);
    if (status != NC_NOERR) {
      status=nc_def_var(exodusFilePtr, VAR_FILE_TYPE, NC_INT, 0, NULL, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define file type in file ID %d",
		exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // Output the file version
    int ierr = ne_put_version(exodusFilePtr);
    if (ierr < 0) return (ierr);

    if (comm.globalNodes > 0) {
      // Define dimension for number of global nodes
      int ltempsv = comm.globalNodes;
      int glonoddim = 0;
      status=nc_def_dim(exodusFilePtr, DIM_NUM_NODES_GLOBAL, ltempsv, &glonoddim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to dimension \"%s\" in file ID %d",
		DIM_NUM_NODES_GLOBAL, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    if (comm.globalElements > 0) {
      // Define dimension for number of global elements
      int ltempsv = comm.globalElements;
      int gloelemdim = 0;
      status=nc_def_dim(exodusFilePtr, DIM_NUM_ELEMS_GLOBAL, ltempsv, &gloelemdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to dimension \"%s\" in file ID %d",
		DIM_NUM_ELEMS_GLOBAL, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // Output the number of global element blocks. This is output as a
    // dimension since the vector of global element block IDs is sized
    // by this quantity.
    // For Sierra, global block count == local block count
    {
      const char *vars[] = {VAR_ELBLK_IDS_GLOBAL,
			    VAR_ELBLK_CNT_GLOBAL,
			    NULL};
      status = define_variables(exodusFilePtr, (int)mesh.blockCount, DIM_NUM_ELBLK_GLOBAL, vars, NC_INT);
      if (status != EX_NOERR) return EX_FATAL;
    }

    // Output the number of global node sets. This is output as a
    // dimension since the vector of global element block IDs is sized
    // by this quantity.
    {
      const char *vars[] = {VAR_NS_IDS_GLOBAL,
			    VAR_NS_NODE_CNT_GLOBAL,
			    VAR_NS_DF_CNT_GLOBAL,
			    NULL};
      status = define_variables(exodusFilePtr, (int)mesh.nodesetCount, DIM_NUM_NS_GLOBAL, vars, NC_INT);
      if (status != EX_NOERR) return EX_FATAL;
    }
	
    // Output the number of global side sets. This is output as a
    // dimension since the vector of global element block IDs is sized
    // by this quantity.
    {
      const char *vars[] = {VAR_SS_IDS_GLOBAL,
			    VAR_SS_SIDE_CNT_GLOBAL,
			    VAR_SS_DF_CNT_GLOBAL,
			    NULL};
      status = define_variables(exodusFilePtr, (int)mesh.sidesetCount, DIM_NUM_SS_GLOBAL, vars, NC_INT);
      if (status != EX_NOERR) return EX_FATAL;
    }

    // Internal Node status
    status = conditional_define_variable(exodusFilePtr, VAR_INT_N_STAT, dimid_npf, &nodeMapVarID[0], NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // Border node status
    status = conditional_define_variable(exodusFilePtr, VAR_BOR_N_STAT, dimid_npf, &nodeMapVarID[1], NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // External Node status
    status = conditional_define_variable(exodusFilePtr, VAR_EXT_N_STAT, dimid_npf, &nodeMapVarID[2], NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // Define the variable IDs for the elemental status vectors
    // Internal elements
    status = conditional_define_variable(exodusFilePtr, VAR_INT_E_STAT, dimid_npf, &elementMapVarID[0], NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // Border elements
    status = conditional_define_variable(exodusFilePtr, VAR_BOR_E_STAT, dimid_npf, &elementMapVarID[1], NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // Define variable for the internal element information
    status = define_variable(exodusFilePtr, comm.elementsInternal, DIM_NUM_INT_ELEMS, VAR_ELEM_MAP_INT, NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // Define variable for the border element information
    status = define_variable(exodusFilePtr, comm.elementsBorder, DIM_NUM_BOR_ELEMS, VAR_ELEM_MAP_BOR, NC_INT);
    if (status != EX_NOERR) return EX_FATAL;
    
    // Define variable for vector of internal FEM node IDs
    status = define_variable(exodusFilePtr, comm.nodesInternal, DIM_NUM_INT_NODES, VAR_NODE_MAP_INT, NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // Define variable for vector of border FEM node IDs
    status = define_variable(exodusFilePtr, comm.nodesBorder, DIM_NUM_BOR_NODES, VAR_NODE_MAP_BOR, NC_INT);
    if (status != EX_NOERR) return EX_FATAL;

    // Define dimension for vector of external FEM node IDs
    status = define_variable(exodusFilePtr, comm.nodesExternal, DIM_NUM_EXT_NODES, VAR_NODE_MAP_EXT, NC_INT);
    if (status != EX_NOERR) return EX_FATAL;
    
    // Add the nodal communication map count

    int ncnt_cmap = 0;
    for(int icm=0; icm < static_cast<int>(comm.nodeMap.size()); icm++) {
      ncnt_cmap += comm.nodeMap[icm].entityCount;
    }

    {
      const char *vars[] = {VAR_N_COMM_IDS,
			    VAR_N_COMM_STAT,
			    VAR_N_COMM_DATA_IDX,
			    NULL};
      
      status = define_variables(exodusFilePtr, (int)comm.nodeMap.size(), DIM_NUM_N_CMAPS, vars, NC_INT);
      if (status != EX_NOERR) return EX_FATAL;
    }
    {
      const char *vars[] = {VAR_N_COMM_NIDS,
			    VAR_N_COMM_PROC,
			    NULL};
      
      // Add dimensions for all of the nodal communication maps
      status = define_variables(exodusFilePtr, ncnt_cmap, DIM_NCNT_CMAP, vars, NC_INT);
      if (status != EX_NOERR) return EX_FATAL;
    }

    // Add the nodal communication map count
    int ecnt_cmap = 0;
    for (int icm=0; icm < static_cast<int>(comm.elementMap.size()); icm++)
      ecnt_cmap += comm.elementMap[icm].entityCount;

    {
      const char *vars[] = {VAR_E_COMM_IDS,
			    VAR_E_COMM_STAT,
			    VAR_E_COMM_DATA_IDX,
			    NULL};
      status = define_variables(exodusFilePtr, (int)comm.elementMap.size(), DIM_NUM_E_CMAPS, vars, NC_INT);
      if (status != EX_NOERR) return EX_FATAL;
    }
    {
      const char *vars[] = {VAR_E_COMM_EIDS,
			    VAR_E_COMM_PROC,
			    VAR_E_COMM_SIDS,
			    NULL};
      status = define_variables(exodusFilePtr, ecnt_cmap, DIM_ECNT_CMAP, vars, NC_INT);
      if (status != EX_NOERR) return EX_FATAL;
    }
  }
  return(EX_NOERR);
}

int Internals::put_metadata(const std::vector<Block> &blocks)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::put_metadata(blocks)";
  int dims[2];

  int status  = 0; // clear error code

  if (blocks.empty())
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
	    "Error: failed to get string length in file id %d",exodusFilePtr);
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

  assert(blocks.size() == num_elem_blk);

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


int Internals::put_non_define_data(const Mesh&,
				   const CommunicationMetaData &comm)
{
  const char *routine = "Internals::put_non_define_data(mesh)";
  char errmsg[MAX_ERR_LENGTH];
  int status = 0;
  // Metadata that must be written outside of define mode...

  // Output the file type
  if (comm.processorCount > 1) {
    int varid;
    status=nc_inq_varid(exodusFilePtr, VAR_FILE_TYPE, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to locate file type in file ID %d",
	      exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    int lftype = 0; // Parallel file...
    status=nc_put_var1_int(exodusFilePtr, varid, 0, &lftype);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: unable to output file type variable in file ID %d",
	      exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    int nmstat;
    size_t start[1];
    size_t count[1];

    nmstat = comm.nodesInternal == 0 ? 0 : 1;
    status=nc_put_var_int(exodusFilePtr, nodeMapVarID[0], &nmstat);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to output status for internal node map in file ID %d",
	      exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    nmstat = comm.nodesBorder == 0 ? 0 : 1;
    status=nc_put_var_int(exodusFilePtr, nodeMapVarID[1], &nmstat);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to output status for border node map in file ID %d",
	      exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    nmstat = comm.nodesExternal == 0 ? 0 : 1;
    status=nc_put_var_int(exodusFilePtr, nodeMapVarID[2], &nmstat);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to output status for external node map in file ID %d",
	      exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    nmstat = comm.elementsInternal == 0 ? 0 : 1;
    status=nc_put_var_int(exodusFilePtr, elementMapVarID[0], &nmstat);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to output status for internal elem map in file ID %d",
	      exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    nmstat = comm.elementsBorder == 0 ? 0 : 1;
    status=nc_put_var_int(exodusFilePtr, elementMapVarID[1], &nmstat);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to output status for border elem map in file ID %d",
	      exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    int ncnt_cmap = 0;
    for(int icm=0; icm < static_cast<int>(comm.nodeMap.size()); icm++) {
      ncnt_cmap += comm.nodeMap[icm].entityCount;
    }

    if (comm.nodeMap.size() > 0 && ncnt_cmap > 0) {
      int n_varid;
      status=nc_inq_varid(exodusFilePtr, VAR_N_COMM_STAT, &n_varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to find variable ID for \"%s\" in file ID %d",
		VAR_N_COMM_STAT, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      count[0] = 1;
      int nl_ncnt_cmap = 0;
      for (int icm=0; icm < static_cast<int>(comm.nodeMap.size()); icm++) {

	start[0] = icm;
	nmstat = comm.nodeMap[icm].entityCount > 0 ? 1 : 0;
	status=nc_put_vara_int(exodusFilePtr, n_varid, start, count, &nmstat);
	if(status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  sprintf(errmsg,
		  "Error: unable to output variable in file ID %d", exodusFilePtr);
	  ex_err(routine, errmsg, status);
	  return (EX_FATAL);
	}

	// increment to the next starting position
	nl_ncnt_cmap += comm.nodeMap[icm].entityCount;

	// fill the cmap data index
	nc_inq_varid(exodusFilePtr, VAR_N_COMM_DATA_IDX, &commIndexVar);
	status=nc_put_vara_int(exodusFilePtr, commIndexVar, start, count,
			       &nl_ncnt_cmap);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  sprintf(errmsg,
		  "Error: failed to output int elem map index in file ID %d",
		  exodusFilePtr);
	  ex_err(routine, errmsg, status);
	  return (EX_FATAL);
	}
      } // End "for(icm=0; icm < num_n_comm_maps; icm++)"

      // Put Communication set ids...
      std::vector<int> node_cmap_ids(comm.nodeMap.size());
      for (int i=0; i < static_cast<int>(comm.nodeMap.size()); i++) {
	node_cmap_ids[i] = comm.nodeMap[i].id;
      }
      if (put_int_array(exodusFilePtr, VAR_N_COMM_IDS, node_cmap_ids) != NC_NOERR)
	return(EX_FATAL);
    }
    // Set the status of the elemental communication maps
    int ecnt_cmap = 0;
    for (int icm=0; icm < static_cast<int>(comm.elementMap.size()); icm++)
      ecnt_cmap += comm.elementMap[icm].entityCount;

    if (comm.elementMap.size() > 0 && ecnt_cmap > 0) {

      // Get variable ID for elemental status vector
      int e_varid;
      status=nc_inq_varid(exodusFilePtr, VAR_E_COMM_STAT, &e_varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to find variable ID for \"%s\" in file ID %d",
		VAR_E_COMM_STAT, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      count[0] = 1;
      int nl_ecnt_cmap = 0; // reset this for index
      for (int icm=0; icm < static_cast<int>(comm.elementMap.size()); icm++) {

	start[0] = icm;
	nmstat = comm.elementMap[icm].entityCount > 0 ? 1 : 0;

	status=nc_put_vara_int(exodusFilePtr, e_varid, start, count, &nmstat);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  sprintf(errmsg,
		  "Error: unable to output variable in file ID %d", exodusFilePtr);
	  ex_err(routine, errmsg, status);
	  return (EX_FATAL);
	}

	// increment to the next starting position
	nl_ecnt_cmap += comm.elementMap[icm].entityCount;

	// fill the cmap data index
	nc_inq_varid(exodusFilePtr, VAR_E_COMM_DATA_IDX, &elemCommIndexVar);
	status=nc_put_vara_int(exodusFilePtr, elemCommIndexVar, start, count,
			       &nl_ecnt_cmap);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  sprintf(errmsg,
		  "Error: failed to output int elem map index in file ID %d",
		  exodusFilePtr);
	  ex_err(routine, errmsg, status);
	  return (EX_FATAL);
	}
      } // End "for(icm=0; icm < num_e_comm_maps; icm++)"

      // Get the variable ID for the elemental comm map IDs vector
      std::vector<int> elem_cmap_ids(comm.elementMap.size());
      for (int i=0; i < static_cast<int>(comm.elementMap.size()); i++) {
	elem_cmap_ids[i] = comm.elementMap[i].id;
      }
      if (put_int_array(exodusFilePtr, VAR_E_COMM_IDS, elem_cmap_ids) != NC_NOERR)
	return(EX_FATAL);
    }
  }
  return EX_NOERR;
}

int Internals::put_non_define_data(const std::vector<Block> &blocks)
{
  int num_elem_blk = (int)blocks.size();  // Verified via assert earlier...

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

    size_t  start[2];
    size_t  count[2];
    std::string text("");
    count[0] = 1;
    start[1] = 0;
    count[1] = text.size()+1;
  
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      if (blocks[iblk].attributeCount > 0 && blocks[iblk].elementCount > 0) {
	int varid;
	nc_inq_varid(exodusFilePtr, VAR_NAME_ATTRIB(iblk+1), &varid);
	
	for (int i = 0; i < blocks[iblk].attributeCount; i++) {
	  start[0] = i;
	  nc_put_vara_text(exodusFilePtr, varid, start, count, text.c_str());
	}
      }
    }
  }
  return (EX_NOERR);
}

// ========================================================================
int Internals::put_metadata(const std::vector<NodeSet> &nodesets)
{
  const char *routine = "Internals::put_metadata(nodesets)";
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

  assert(static_cast<int>(nodesets.size()) == num_node_sets);

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

int Internals::put_non_define_data(const std::vector<NodeSet> &nodesets)
{
  if (nodesets.empty())
    return EX_NOERR;

  // Output nodeset ids...
  size_t num_nodesets = nodesets.size();
  std::vector<int> nodeset_id(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    nodeset_id[i] = nodesets[i].id;
  }

  if (put_int_array(exodusFilePtr, VAR_NS_IDS, nodeset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    status[i] = nodesets[i].nodeCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_NS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

// ========================================================================
int Internals::put_metadata(const std::vector<SideSet> &sidesets)
{
  const char *routine = "Internals::put_metadata(sidesets)";
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

  assert(static_cast<int>(sidesets.size()) == num_side_sets);

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

int Internals::put_non_define_data(const std::vector<SideSet> &sidesets)
{
  if (sidesets.empty())
    return EX_NOERR;

  // Output sideset ids...
  int num_sidesets = (int)sidesets.size();
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
  int conditional_define_variable(int exodusFilePtr, const char *var, int dimid, int *varid,  nc_type type)
  {
    int status;
    
    char errmsg[MAX_ERR_LENGTH];
    const char *routine = "Internals::conditional_define_variable()";
    
    status=nc_inq_varid(exodusFilePtr, var, varid);
    if (status != NC_NOERR) {
      status=nc_def_var(exodusFilePtr, var, type, 1, &dimid, varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: Failed to define variable \"%s\" in file ID %d",
		var, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }
    return EX_NOERR;
  }

  int define_variable(int exodusFilePtr, int size, const char *dim, const char *var, nc_type type)
  {
    int dimid[1];
    int varid;
    int status;
    
    char errmsg[MAX_ERR_LENGTH];
    const char *routine = "Internals::define_variable()";

    if (size > 0) {
      status=nc_def_dim(exodusFilePtr, dim, size, &dimid[0]);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to dimension \"%s\" in file id %d",
		DIM_NUM_BOR_ELEMS, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      status=nc_def_var(exodusFilePtr, var, type, 1, dimid, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define variable \"%s\" in file ID %d",
		VAR_ELEM_MAP_BOR, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    } 
    return EX_NOERR;
  }

  int define_variables(int exodusFilePtr, int size, const char *dim, const char *var[], nc_type type)
  {
    int dimid[1];
    int varid;
    int status;
    
    char errmsg[MAX_ERR_LENGTH];
    const char *routine = "Internals::define_variables()";

    if (size > 0) {
      status=nc_def_dim(exodusFilePtr, dim, size, &dimid[0]);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to dimension \"%s\" in file id %d",
		DIM_NUM_BOR_ELEMS, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      int i = 0;
      while (var[i] != NULL) {
	status=nc_def_var(exodusFilePtr, var[i], type, 1, dimid, &varid);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  sprintf(errmsg,
		  "Error: failed to define variable \"%s\" in file ID %d",
		  VAR_ELEM_MAP_BOR, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	  return (EX_FATAL);
	}
	i++;
      }
    }
    return EX_NOERR;
  }

  int put_int_array(int exoid, const char *var_type, const std::vector<int> &array)
  {
    const char *routine = "Ioex_Internals.C, put_int_array";
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

