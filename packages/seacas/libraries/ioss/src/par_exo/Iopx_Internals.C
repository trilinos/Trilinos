// Copyright(C) 2012
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

// This must appear before exodusII_int.h include.
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#ifndef PRId64
#error "PRId64 not defined"
#endif

#include <Ioss_Utils.h>                 // for IOSS_WARNING
#include <assert.h>                     // for assert
#include <par_exo/Iopx_Internals.h>    // for Internals, ElemBlock, etc
extern "C" {
#include <exodusII_int.h>               // for EX_FATAL, EX_NOERR, etc
}
#include <netcdf.h>                     // for NC_NOERR, nc_def_var, etc
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf, NULL
#include <stdlib.h>                     // for exit, EXIT_FAILURE
#include <string.h>                     // for strlen, strncpy, strcpy, etc
#include <ostream>                      // for operator<<, etc
#include <string>                       // for string, operator==, etc
#include <vector>                       // for vector

#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_VariableType.h"
#include "exodusII.h"                   // for ex_err, ex_opts, etc

using namespace Iopx;

namespace {
  nc_type get_type(int exoid, unsigned int type)
  {
    if (ex_int64_status(exoid) & type)
      return NC_INT64;
    else
      return NC_INT;
  }
  int define_netcdf_vars(int exoid, const char *type, size_t count, const char *dim_num,
			 const char *stat_var, const char *id_var, const char *name_var);
  int define_variable(int exodusFilePtr, int64_t size, const char *dim, const char *var, nc_type type);
  int define_variables(int exodusFilePtr, int64_t size, const char *dim, const char *var[], const nc_type *types);
  int conditional_define_variable(int exodusFilePtr, const char *var, int dimid, int *varid,  nc_type type);

  size_t max_string_length() {return MAX_STR_LENGTH;}
  int put_int_array(int exoid, const char *var_type,  const std::vector<int> &array);
  int put_id_array(int exoid, const char *var_type,  const std::vector<entity_id> &array);
  int define_coordinate_vars(int exodusFilePtr, int64_t nodes, int node_dim,
			     int dimension, int dim_dim, int str_dim);
  template <typename T>
  int output_names(const std::vector<T> &entities, int exoid, ex_entity_type nc_type);
}

#define stringify(X) #X
#define version_string(X) stringify(X)

void Internals::update_last_time_attribute(double value)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::update_last_time_attribute()";

  int status=nc_put_att_double(exodusFilePtr, NC_GLOBAL, "last_written_time",
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

  nc_type att_type = NC_NAT;
  size_t att_len = 0;
  int status = nc_inq_att(exodusFilePtr, NC_GLOBAL, "last_written_time", &att_type, &att_len);
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
    ex_err("Iopx::Redefine::Redefine()",errmsg,status);
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
      ex_err("Iopx::Redefine::~Redefine()",errmsg,status);
      exit(EXIT_FAILURE);
    }
  } catch (...) {
  }
}

NodeBlock::NodeBlock(const Ioss::NodeBlock &other)
{
  name = other.name();
  if (other.property_exists("id")) {
    id = other.get_property("id").get_int();
  } else {
    id = 1;
  }
  if (other.property_exists("locally_owned_count")) {
    entityCount = other.get_property("locally_owned_count").get_int();
  } else {
    entityCount = other.get_property("entity_count").get_int();
  }
  attributeCount = other.get_property("attribute_count").get_int();
}

NodeBlock& NodeBlock::operator=(const NodeBlock& other)
{
  name = other.name;
  id = other.id;
  entityCount = other.entityCount;
  attributeCount = other.attributeCount;
  return *this;
}

EdgeBlock::EdgeBlock(const Ioss::EdgeBlock &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  nodesPerEntity = other.get_property("topology_node_count").get_int();
  attributeCount = other.get_property("attribute_count").get_int();

  std::string el_type = other.get_property("topology_type").get_string();
  if (other.property_exists("original_topology_type")) {
    el_type = other.get_property("original_topology_type").get_string();
  }

  std::strncpy(elType, el_type.c_str(), max_string_length());
  elType[max_string_length()] = 0;
}

EdgeBlock& EdgeBlock::operator=(const EdgeBlock& other)
{
  name = other.name;
  id = other.id;
  entityCount = other.entityCount;
  nodesPerEntity = other.nodesPerEntity;
  attributeCount = other.attributeCount;
  std::strcpy(elType, other.elType);
  return *this;
}

bool EdgeBlock::operator==(const EdgeBlock& other) const
{
  return name == other.name &&
    id == other.id &&
    entityCount == other.entityCount &&
    nodesPerEntity == other.nodesPerEntity &&
    attributeCount == other.attributeCount;
}

FaceBlock::FaceBlock(const Ioss::FaceBlock &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  nodesPerEntity = other.get_property("topology_node_count").get_int();
  if (other.field_exists("connectivty_edge")) {
    edgesPerEntity = other.get_field("connectivity_edge").raw_storage()->component_count();
  } else {
    edgesPerEntity = 0;
  }
  attributeCount = other.get_property("attribute_count").get_int();

  std::string el_type = other.get_property("topology_type").get_string();
  if (other.property_exists("original_topology_type")) {
    el_type = other.get_property("original_topology_type").get_string();
  }

  std::strncpy(elType, el_type.c_str(), max_string_length());
  elType[max_string_length()] = 0;
}

FaceBlock& FaceBlock::operator=(const FaceBlock& other)
{
  name = other.name;
  id = other.id;
  entityCount = other.entityCount;
  nodesPerEntity = other.nodesPerEntity;
  edgesPerEntity = other.edgesPerEntity;
  attributeCount = other.attributeCount;
  std::strcpy(elType, other.elType);
  return *this;
}

bool FaceBlock::operator==(const FaceBlock& other) const
{
  return name == other.name &&
    id == other.id &&
    entityCount == other.entityCount &&
    nodesPerEntity == other.nodesPerEntity &&
    edgesPerEntity == other.edgesPerEntity &&
    attributeCount == other.attributeCount;
}

ElemBlock::ElemBlock(const Ioss::ElementBlock &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  nodesPerEntity = other.get_property("topology_node_count").get_int();

  if (other.field_exists("connectivity_edge")) {
    edgesPerEntity = other.get_field("connectivity_edge").raw_storage()->component_count();
  } else {
    edgesPerEntity = 0;
  }

  if (other.field_exists("connectivity_face")) {
    facesPerEntity = other.get_field("connectivity_face").raw_storage()->component_count();
  } else {
    facesPerEntity = 0;
  }

  attributeCount = other.get_property("attribute_count").get_int();
  std::string el_type = other.get_property("topology_type").get_string();
  if (other.property_exists("original_topology_type")) {
    el_type = other.get_property("original_topology_type").get_string();
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

ElemBlock& ElemBlock::operator=(const ElemBlock& other)
{
  name = other.name;
  id = other.id;
  entityCount = other.entityCount;
  nodesPerEntity = other.nodesPerEntity;
  edgesPerEntity = other.edgesPerEntity;
  facesPerEntity = other.facesPerEntity;
  attributeCount = other.attributeCount;
  std::strcpy(elType, other.elType);
  return *this;
}

bool ElemBlock::operator==(const ElemBlock& other) const
{
  return name == other.name && id == other.id &&
    entityCount == other.entityCount &&
    nodesPerEntity == other.nodesPerEntity &&
    edgesPerEntity == other.edgesPerEntity &&
    facesPerEntity == other.facesPerEntity &&
    attributeCount == other.attributeCount &&
    std::strcmp(elType, other.elType) == 0;
}

NodeSet::NodeSet(const Ioss::NodeSet &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  if (other.property_exists("locally_owned_count")) {
    entityCount = other.get_property("locally_owned_count").get_int();
  } else {
    entityCount = other.get_property("entity_count").get_int();
  }
  attributeCount = other.get_property("attribute_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
  if (dfCount > 0 && dfCount != entityCount) {
    dfCount = entityCount;
  }
}

bool NodeSet::operator==(const NodeSet& other) const
{
  return id == other.id &&
    entityCount == other.entityCount &&
    dfCount == other.dfCount &&
    name == other.name;
}

EdgeSet::EdgeSet(const Ioss::EdgeSet &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  attributeCount = other.get_property("attribute_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
}

bool EdgeSet::operator==(const EdgeSet& other) const
{
  return id == other.id &&
    entityCount == other.entityCount &&
    dfCount == other.dfCount &&
    name == other.name;
}

FaceSet::FaceSet(const Ioss::FaceSet &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  attributeCount = other.get_property("attribute_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
}

bool FaceSet::operator==(const FaceSet& other) const
{
  return id == other.id &&
    entityCount == other.entityCount &&
    dfCount == other.dfCount &&
    name == other.name;
}

ElemSet::ElemSet(const Ioss::ElementSet &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  attributeCount = other.get_property("attribute_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
}

bool ElemSet::operator==(const ElemSet& other) const
{
  return id == other.id &&
    entityCount == other.entityCount &&
    dfCount == other.dfCount &&
    name == other.name;
}

SideSet::SideSet(const Ioss::SideBlock &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
  std::string io_name = other.name();

  // KLUGE: universal_sideset has side dfCount...
  if (io_name == "universal_sideset")
    dfCount = entityCount;
}

SideSet::SideSet(const Ioss::SideSet &other)
{
  name = other.name();
  id = other.get_property("id").get_int();
  entityCount = other.get_property("entity_count").get_int();
  dfCount = other.get_property("distribution_factor_count").get_int();
  std::string io_name = other.name();

  // KLUGE: universal_sideset has side dfCount...
  if (io_name == "universal_sideset")
    dfCount = entityCount;
}

bool SideSet::operator==(const SideSet& other) const
{
  return id == other.id &&
    entityCount == other.entityCount &&
    dfCount == other.dfCount &&
    name == other.name;
}

Internals::Internals(int exoid, int maximum_name_length, const Ioss::ParallelUtils &util)
  : exodusFilePtr(exoid),
    nodeMapVarID(),
    elementMapVarID(),
    commIndexVar(0),
    elemCommIndexVar(0),
    maximumNameLength(maximum_name_length),
    parallelUtil(util)
{}

int Internals::write_meta_data(Mesh &mesh)
{
  int ierr;
  {
    // Determine global counts...
    get_global_counts(mesh);

    Redefine the_database(exodusFilePtr);

    // Set the database to NOFILL mode.  Only writes values we want written...
    int old_fill = 0;

    ierr=nc_set_fill(exodusFilePtr, NC_NOFILL, &old_fill);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.edgeblocks);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.faceblocks);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.elemblocks);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.nodesets);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.edgesets);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.facesets);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.elemsets);
    if (ierr != EX_NOERR) return(ierr);

    ierr=put_metadata(mesh.sidesets);
    if (ierr != EX_NOERR) return(ierr);
  }


  // NON-Define mode output...
  size_t my_proc = parallelUtil.parallel_rank();
  ierr=put_non_define_data(mesh);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(mesh.edgeblocks);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(mesh.faceblocks);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(mesh.elemblocks);
  if (ierr != EX_NOERR) return(ierr);

  ierr=put_non_define_data(mesh.nodesets);
  if (ierr != EX_NOERR) return(ierr);
  
  ierr=put_non_define_data(mesh.edgesets);
  if (ierr != EX_NOERR) return(ierr);
  
  ierr=put_non_define_data(mesh.facesets);
  if (ierr != EX_NOERR) return(ierr);
  
  ierr=put_non_define_data(mesh.elemsets);
  if (ierr != EX_NOERR) return(ierr);
  
  ierr=put_non_define_data(mesh.sidesets);
  if (ierr != EX_NOERR) return(ierr);

  // For now, put entity names using the ExodusII api...
  output_names(mesh.edgeblocks, exodusFilePtr, EX_EDGE_BLOCK);
  output_names(mesh.faceblocks, exodusFilePtr, EX_FACE_BLOCK);
  output_names(mesh.elemblocks, exodusFilePtr, EX_ELEM_BLOCK);
  output_names(mesh.nodesets,   exodusFilePtr, EX_NODE_SET);
  output_names(mesh.edgesets,   exodusFilePtr, EX_EDGE_SET);
  output_names(mesh.facesets,   exodusFilePtr, EX_FACE_SET);
  output_names(mesh.elemsets,   exodusFilePtr, EX_ELEM_SET);
  output_names(mesh.sidesets,   exodusFilePtr, EX_SIDE_SET);

  return(EX_NOERR);
}

void Internals::get_global_counts(Mesh &mesh)
{
  std::vector<int64_t> counts;
  std::vector<int64_t> global_counts;
  
  counts.push_back(mesh.nodeblocks[0].entityCount);
  for (size_t i=0; i < mesh.edgeblocks.size(); i++) {
    counts.push_back(mesh.edgeblocks[i].entityCount);
  }
  for (size_t i=0; i < mesh.faceblocks.size(); i++) {
    counts.push_back(mesh.faceblocks[i].entityCount);
  }
  for (size_t i=0; i < mesh.elemblocks.size(); i++) {
    counts.push_back(mesh.elemblocks[i].entityCount);
  }
  for (size_t i=0; i < mesh.nodesets.size(); i++) {
    counts.push_back(mesh.nodesets[i].entityCount);
    counts.push_back(mesh.nodesets[i].dfCount);
  }
  for (size_t i=0; i < mesh.edgesets.size(); i++) {
    counts.push_back(mesh.edgesets[i].entityCount);
    counts.push_back(mesh.edgesets[i].dfCount);
  }
  for (size_t i=0; i < mesh.facesets.size(); i++) {
    counts.push_back(mesh.facesets[i].entityCount);
    counts.push_back(mesh.facesets[i].dfCount);
  }
  for (size_t i=0; i < mesh.elemsets.size(); i++) {
    counts.push_back(mesh.elemsets[i].entityCount);
    counts.push_back(mesh.elemsets[i].dfCount);
  }
  for (size_t i=0; i < mesh.sidesets.size(); i++) {
    counts.push_back(mesh.sidesets[i].entityCount);
    counts.push_back(mesh.sidesets[i].dfCount);
  }

  // Now gather this information on each processor so
  // they can determine the offsets and totals...
  global_counts.resize(counts.size() * parallelUtil.parallel_size());
  
  MPI_Allgather(&counts[0],        counts.size(), MPI_LONG_LONG_INT,
		&global_counts[0], counts.size(), MPI_LONG_LONG_INT,
		parallelUtil.communicator());
  
  std::vector<int64_t> offsets(counts.size());
  
  size_t my_proc = parallelUtil.parallel_rank();
  size_t proc_count = parallelUtil.parallel_size();
  
  // Calculate offsets for each entity on each processor
  for (size_t j=0; j < offsets.size(); j++) {
    for (size_t i=0; i < my_proc; i++) {
      offsets[j] += global_counts[i*offsets.size() + j];
    }
  }
  
  // Now calculate the total count of entities over all processors
  for (size_t j=0; j < offsets.size(); j++) {
    for (size_t i=1; i < proc_count; i++) {
      global_counts[j] += global_counts[i*offsets.size() + j];
    }
  }

  size_t j = 0;
  mesh.nodeblocks[0].procOffset = offsets[j];
  mesh.nodeblocks[0].entityCount = global_counts[j++];

  for (size_t i=0; i < mesh.edgeblocks.size(); i++) {
    mesh.edgeblocks[i].procOffset = offsets[j];
    mesh.edgeblocks[i].entityCount = global_counts[j++];
  }
  for (size_t i=0; i < mesh.faceblocks.size(); i++) {
    mesh.faceblocks[i].procOffset = offsets[j];
    mesh.faceblocks[i].entityCount = global_counts[j++];
  }
  for (size_t i=0; i < mesh.elemblocks.size(); i++) {
    mesh.elemblocks[i].procOffset = offsets[j];
    mesh.elemblocks[i].entityCount = global_counts[j++];
  }
  for (size_t i=0; i < mesh.nodesets.size(); i++) {
    mesh.nodesets[i].procOffset = offsets[j];
    mesh.nodesets[i].entityCount = global_counts[j++];
    mesh.nodesets[i].dfCount = global_counts[j++];
  }
  for (size_t i=0; i < mesh.edgesets.size(); i++) {
    mesh.edgesets[i].procOffset = offsets[j];
    mesh.edgesets[i].entityCount = global_counts[j++];
    mesh.edgesets[i].dfCount = global_counts[j++];
  }
  for (size_t i=0; i < mesh.facesets.size(); i++) {
    mesh.facesets[i].procOffset = offsets[j];
    mesh.facesets[i].entityCount = global_counts[j++];
    mesh.facesets[i].dfCount = global_counts[j++];
  }
  for (size_t i=0; i < mesh.elemsets.size(); i++) {
    mesh.elemsets[i].procOffset = offsets[j];
    mesh.elemsets[i].entityCount = global_counts[j++];
    mesh.elemsets[i].dfCount = global_counts[j++];
  }
  for (size_t i=0; i < mesh.sidesets.size(); i++) {
    mesh.sidesets[i].procOffset = offsets[j];
    mesh.sidesets[i].entityCount = global_counts[j++];
    mesh.sidesets[i].dfProcOffset = offsets[j];
    mesh.sidesets[i].dfCount = global_counts[j++];
  }
}

int Internals::put_metadata(const Mesh &mesh)
{
  int numdimdim  = 0;
  int numnoddim  = 0;
  int strdim     = 0;
  int namestrdim = 0;
  int varid = 0;

  int map_type  = get_type(exodusFilePtr, EX_MAPS_INT64_DB);
  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);
  int ids_type  = get_type(exodusFilePtr, EX_IDS_INT64_DB);

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

  if (mesh.nodeblocks[0].entityCount > 0) {
    status=nc_def_dim(exodusFilePtr, DIM_NUM_NODES, mesh.nodeblocks[0].entityCount, &numnoddim);
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

  if (mesh.nodeblocks[0].attributeCount > 0) {
    int numattrdim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_NBLK,
		       mesh.nodeblocks[0].attributeCount, &numattrdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of attributes in node block %"PRId64" in file id %d",
	      (entity_id)mesh.nodeblocks[0].id,exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    int dims[2];
    dims[0] = numnoddim;
    dims[1] = numattrdim;
    status=nc_def_var(exodusFilePtr, VAR_NATTRIB,
		      nc_flt_code(exodusFilePtr), 2, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error:  failed to define attributes for node block %"PRId64" in file id %d",
	      mesh.nodeblocks[0].id,exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    // Attribute name array...
    dims[0] = numattrdim;
    dims[1] = namestrdim;

    status = nc_def_var(exodusFilePtr, VAR_NAME_NATTRIB,
			NC_CHAR, 2, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define attribute name array for node block %"PRId64" in file id %d",
	      mesh.nodeblocks[0].id,exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }
  }

  size_t elem_count = 0;
  for (size_t i=0; i < mesh.elemblocks.size(); i++) {
    elem_count += mesh.elemblocks[i].entityCount;
  }

  if (elem_count > 0) {
    int numelemdim;
    status=nc_def_dim(exodusFilePtr, DIM_NUM_ELEM, elem_count, &numelemdim);
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

  size_t face_count = 0;
  for (size_t i=0; i < mesh.faceblocks.size(); i++) {
    face_count += mesh.faceblocks[i].entityCount;
  }

  if (face_count > 0) {
    int numfacedim;
    status=nc_def_dim(exodusFilePtr, DIM_NUM_FACE, face_count, &numfacedim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of faces in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // Define the face map here to avoid a later redefine call
    int dims[1];
    dims[0] = numfacedim;
    varid = 0;
    status=nc_def_var(exodusFilePtr, VAR_FACE_NUM_MAP, map_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: face numbering map already exists in file id %d",
		exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create face numbering map in file id %d",
		exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return (EX_FATAL);
    }
  }

  size_t edge_count = 0;
  for (size_t i=0; i < mesh.edgeblocks.size(); i++) {
    edge_count += mesh.edgeblocks[i].entityCount;
  }

  if (edge_count > 0) {
    int numedgedim;
    status=nc_def_dim(exodusFilePtr, DIM_NUM_EDGE, edge_count, &numedgedim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of edges in file id %d",exodusFilePtr);
      ex_err(routine,errmsg,status);
      return(EX_FATAL);
    }

    // Define the edge map here to avoid a later redefine call
    int dims[1];
    dims[0] = numedgedim;
    varid = 0;
    status=nc_def_var(exodusFilePtr, VAR_EDGE_NUM_MAP, map_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: edge numbering map already exists in file id %d",
		exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create edge numbering map in file id %d",
		exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return (EX_FATAL);
    }
  }

  // ========================================================================
  // Blocks...
  if (!mesh.edgeblocks.empty()) {
    status = define_netcdf_vars(exodusFilePtr, "edge block", mesh.edgeblocks.size(),
				DIM_NUM_ED_BLK, VAR_STAT_ED_BLK, VAR_ID_ED_BLK, VAR_NAME_ED_BLK);
  }

  if (!mesh.faceblocks.empty()) {
    status = define_netcdf_vars(exodusFilePtr, "face block", mesh.faceblocks.size(),
				DIM_NUM_FA_BLK, VAR_STAT_FA_BLK, VAR_ID_FA_BLK, VAR_NAME_FA_BLK);
  }

  if (!mesh.elemblocks.empty()) {
    status = define_netcdf_vars(exodusFilePtr, "element block", mesh.elemblocks.size(),
				DIM_NUM_EL_BLK, VAR_STAT_EL_BLK, VAR_ID_EL_BLK, VAR_NAME_EL_BLK);
  }

  // ========================================================================
  // Sets...
  if (!mesh.nodesets.empty()) {
    status = define_netcdf_vars(exodusFilePtr, "node set", mesh.nodesets.size(),
				DIM_NUM_NS, VAR_NS_STAT, VAR_NS_IDS, VAR_NAME_NS);
  }

  if (!mesh.edgesets.empty()) {
    status = define_netcdf_vars(exodusFilePtr, "edge set", mesh.edgesets.size(),
				DIM_NUM_ES, VAR_ES_STAT, VAR_ES_IDS, VAR_NAME_ES);
  }

  if (!mesh.facesets.empty()) {
    status = define_netcdf_vars(exodusFilePtr, "face set", mesh.facesets.size(),
				DIM_NUM_FS, VAR_FS_STAT, VAR_FS_IDS, VAR_NAME_FS);
  }

  if (!mesh.elemsets.empty()) {
    status = define_netcdf_vars(exodusFilePtr, "element set", mesh.elemsets.size(),
				DIM_NUM_ELS, VAR_ELS_STAT, VAR_ELS_IDS, VAR_NAME_ELS);
  }

  // ========================================================================
  // side sets...
  if (!mesh.sidesets.empty() > 0) {
    status = define_netcdf_vars(exodusFilePtr, "side set", mesh.sidesets.size(),
				DIM_NUM_SS, VAR_SS_STAT, VAR_SS_IDS, VAR_NAME_SS);
  }

  // ========================================================================
  status = define_coordinate_vars(exodusFilePtr, mesh.nodeblocks[0].entityCount, numnoddim,
				  mesh.dimensionality, numdimdim, namestrdim);
  if (status != EX_NOERR) return EX_FATAL;

  return(EX_NOERR);
}

int Internals::put_metadata(const std::vector<ElemBlock> &blocks)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::put_metadata(blocks)";
  int dims[2];

  int status  = 0; // clear error code

  if (blocks.empty())
    return (EX_NOERR);

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

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

  // Iterate over blocks ...
  for (size_t iblk = 0; iblk < num_elem_blk; iblk++) {
    ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_ELEM_BLOCK));

    if (blocks[iblk].entityCount == 0)
      continue;

    // define some dimensions and variables
    int numelbdim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_EL_IN_BLK(iblk+1),
		       blocks[iblk].entityCount, &numelbdim);
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
		       blocks[iblk].nodesPerEntity, &nelnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of nodes/element for block %"PRId64" in file id %d",
	      blocks[iblk].id, exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    // element connectivity array
    {
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

    if (blocks[iblk].edgesPerEntity > 0) {
      int neledgdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_EDG_PER_EL(iblk+1),
			 blocks[iblk].edgesPerEntity, &neledgdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of edges/element for block %"PRId64" in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // element->edge connectivity array
      dims[0] = numelbdim;
      dims[1] = neledgdim;

      int connid;
      status=nc_def_var(exodusFilePtr, VAR_ECONN(iblk+1), bulk_type, 2, dims, &connid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to create element->edge connectivity array for block %"PRId64" in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    if (blocks[iblk].facesPerEntity > 0) {
      int nelfacdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_FAC_PER_EL(iblk+1),
			 blocks[iblk].facesPerEntity, &nelfacdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of faces/element for block %"PRId64" in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // element->face connectivity array
      dims[0] = numelbdim;
      dims[1] = nelfacdim;

      int connid;
      status=nc_def_var(exodusFilePtr, VAR_FCONN(iblk+1), bulk_type, 2, dims, &connid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to create element->edge connectivity array for block %"PRId64" in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
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
  }
  return (EX_NOERR);
}

int Internals::put_metadata(const std::vector<FaceBlock> &blocks)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::put_metadata(blocks)";
  int dims[2];

  int status  = 0; // clear error code

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  if (blocks.empty())
    return (EX_NOERR);

  // Get number of face blocks defined for this file
  int dimid;
  size_t num_face_blk = 0;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_FA_BLK, &dimid);
  if (status != NC_NOERR ) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: no face blocks defined in file id %d",
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

  status=nc_inq_dimlen(exodusFilePtr,dimid,&num_face_blk);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get number of face blocks in file id %d",
	    exodusFilePtr);
    ex_err(routine, errmsg, status);
    return (EX_FATAL);
  }

  assert(blocks.size() == num_face_blk);

  // Iterate over blocks ...
  for (size_t iblk = 0; iblk < num_face_blk; iblk++) {
    ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_FACE_BLOCK));

    if (blocks[iblk].entityCount == 0)
      continue;

    // define some dimensions and variables
    int numelbdim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_FA_IN_FBLK(iblk+1),
		       blocks[iblk].entityCount, &numelbdim);
    if (status != NC_NOERR) {
      if (status == NC_ENAMEINUSE) {	// duplicate entry
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: face block %"PRId64" already defined in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of faces/block for block %"PRId64" file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return (EX_FATAL);
    }

    int nelnoddim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_NOD_PER_FA(iblk+1),
		       blocks[iblk].nodesPerEntity, &nelnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of nodes/face for block %"PRId64" in file id %d",
	      blocks[iblk].id, exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    // face attribute array
    if (blocks[iblk].attributeCount > 0) {
      int numattrdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_FBLK(iblk+1),
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
      status=nc_def_var(exodusFilePtr, VAR_FATTRIB(iblk+1),
			nc_flt_code(exodusFilePtr), 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error:  failed to define attributes for face block %"PRId64" in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_FATTRIB(iblk+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for face block %"PRId64" in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // face connectivity array
    dims[0] = numelbdim;
    dims[1] = nelnoddim;

    int connid;
    status=nc_def_var(exodusFilePtr, VAR_FBCONN(iblk+1), bulk_type, 2, dims, &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to create connectivity array for block %"PRId64" in file id %d",
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


int Internals::put_metadata(const std::vector<EdgeBlock> &blocks)
{
  char errmsg[MAX_ERR_LENGTH];
  const char *routine = "Internals::put_metadata(blocks)";
  int dims[2];

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status  = 0; // clear error code

  if (blocks.empty())
    return (EX_NOERR);

  // Get number of edge blocks defined for this file
  int dimid;
  size_t num_edge_blk = 0;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_ED_BLK, &dimid);
  if (status != NC_NOERR ) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: no edge blocks defined in file id %d",
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

  status=nc_inq_dimlen(exodusFilePtr,dimid,&num_edge_blk);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get number of edge blocks in file id %d",
	    exodusFilePtr);
    ex_err(routine, errmsg, status);
    return (EX_FATAL);
  }

  assert(blocks.size() == num_edge_blk);

  // Iterate over blocks ...
  for (size_t iblk = 0; iblk < num_edge_blk; iblk++) {
    ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_EDGE_BLOCK));

    if (blocks[iblk].entityCount == 0)
      continue;

    // define some dimensions and variables
    int numelbdim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_ED_IN_EBLK(iblk+1),
		       blocks[iblk].entityCount, &numelbdim);
    if (status != NC_NOERR) {
      if (status == NC_ENAMEINUSE) {	// duplicate entry
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: edge block %"PRId64" already defined in file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of edges/block for block %"PRId64" file id %d",
		blocks[iblk].id, exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return (EX_FATAL);
    }

    int nelnoddim;
    status=nc_def_dim (exodusFilePtr, DIM_NUM_NOD_PER_ED(iblk+1),
		       blocks[iblk].nodesPerEntity, &nelnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to define number of nodes/edge (%"PRId64") for block %"PRId64" in file id %d",
	      blocks[iblk].nodesPerEntity, blocks[iblk].id, exodusFilePtr);
      ex_err(routine, errmsg, status);
      return (EX_FATAL);
    }

    // edge attribute array
    if (blocks[iblk].attributeCount > 0) {
      int numattrdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_EBLK(iblk+1),
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
      status=nc_def_var(exodusFilePtr, VAR_EATTRIB(iblk+1),
			nc_flt_code(exodusFilePtr), 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error:  failed to define attributes for edge block %"PRId64" in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_EATTRIB(iblk+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for edge block %"PRId64" in file id %d",
		blocks[iblk].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }

    // edge connectivity array
    dims[0] = numelbdim;
    dims[1] = nelnoddim;

    int connid;
    status=nc_def_var(exodusFilePtr, VAR_EBCONN(iblk+1), bulk_type, 2, dims, &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      sprintf(errmsg,
	      "Error: failed to create connectivity array for block %"PRId64" in file id %d",
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


int Internals::put_non_define_data(const Mesh&)
{
  return EX_NOERR;
}

int Internals::put_non_define_data(const std::vector<ElemBlock> &blocks)
{
  int num_elem_blk = (int)blocks.size();  // Verified via assert earlier...

  if (num_elem_blk > 0) {
    // first get id of element block ids array variable
    std::vector<entity_id> elem_blk_id(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_id[iblk] = blocks[iblk].id;
    }

    if (put_id_array(exodusFilePtr, VAR_ID_EL_BLK, elem_blk_id) != NC_NOERR)
      return(EX_FATAL);

    // Now, write the element block status array
    std::vector<int> elem_blk_status(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_status[iblk] = blocks[iblk].entityCount > 0 ? 1 : 0;
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
      if (blocks[iblk].attributeCount > 0 && blocks[iblk].entityCount > 0) {
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

int Internals::put_non_define_data(const std::vector<FaceBlock> &blocks)
{
  int num_face_blk = (int)blocks.size();  // Verified via assert earlier...

  if (num_face_blk > 0) {
    // first get id of face block ids array variable
    std::vector<entity_id> face_blk_id(num_face_blk);
    for (int iblk = 0; iblk < num_face_blk; iblk++) {
      face_blk_id[iblk] = blocks[iblk].id;
    }

    if (put_id_array(exodusFilePtr, VAR_ID_FA_BLK, face_blk_id) != NC_NOERR)
      return(EX_FATAL);

    // Now, write the face block status array
    std::vector<int> face_blk_status(num_face_blk);
    for (int iblk = 0; iblk < num_face_blk; iblk++) {
      face_blk_status[iblk] = blocks[iblk].entityCount > 0 ? 1 : 0;
    }

    if (put_int_array(exodusFilePtr, VAR_STAT_FA_BLK, face_blk_status) != NC_NOERR)
      return(EX_FATAL);

    size_t  start[2];
    size_t  count[2];
    std::string text("");
    count[0] = 1;
    start[1] = 0;
    count[1] = text.size()+1;
  
    for (int iblk = 0; iblk < num_face_blk; iblk++) {
      if (blocks[iblk].attributeCount > 0 && blocks[iblk].entityCount > 0) {
	int varid;
	nc_inq_varid(exodusFilePtr, VAR_NAME_FATTRIB(iblk+1), &varid);
	
	for (int i = 0; i < blocks[iblk].attributeCount; i++) {
	  start[0] = i;
	  nc_put_vara_text(exodusFilePtr, varid, start, count, text.c_str());
	}
      }
    }
  }
  return (EX_NOERR);
}

int Internals::put_non_define_data(const std::vector<EdgeBlock> &blocks)
{
  int num_edge_blk = (int)blocks.size();  // Verified via assert earlier...

  if (num_edge_blk > 0) {
    // first get id of edge block ids array variable
    std::vector<entity_id> edge_blk_id(num_edge_blk);
    for (int iblk = 0; iblk < num_edge_blk; iblk++) {
      edge_blk_id[iblk] = blocks[iblk].id;
    }

    if (put_id_array(exodusFilePtr, VAR_ID_ED_BLK, edge_blk_id) != NC_NOERR)
      return(EX_FATAL);

    // Now, write the edge block status array
    std::vector<int> edge_blk_status(num_edge_blk);
    for (int iblk = 0; iblk < num_edge_blk; iblk++) {
      edge_blk_status[iblk] = blocks[iblk].entityCount > 0 ? 1 : 0;
    }

    if (put_int_array(exodusFilePtr, VAR_STAT_ED_BLK, edge_blk_status) != NC_NOERR)
      return(EX_FATAL);

    size_t  start[2];
    size_t  count[2];
    std::string text("");
    count[0] = 1;
    start[1] = 0;
    count[1] = text.size()+1;
  
    for (int iblk = 0; iblk < num_edge_blk; iblk++) {
      if (blocks[iblk].attributeCount > 0 && blocks[iblk].entityCount > 0) {
	int varid;
	nc_inq_varid(exodusFilePtr, VAR_NAME_EATTRIB(iblk+1), &varid);
	
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

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

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
  num_node_sets = ex_inquire_int(exodusFilePtr, EX_INQ_NODE_SETS);

  int namestrdim;
  status=nc_inq_dimid (exodusFilePtr, DIM_STR_NAME, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get string length in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  assert(static_cast<int>(nodesets.size()) == num_node_sets);

  for (int i=0; i<num_node_sets; i++) {

    //  NOTE: ex_inc_file_item is used to find the number of node sets
    // for a specific file and returns that value incremented.
    int cur_num_node_sets=(int)ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_NODE_SET));

    if (nodesets[i].entityCount == 0)
      continue;

    status=nc_def_dim (exodusFilePtr, DIM_NUM_NOD_NS(cur_num_node_sets+1),
		       nodesets[i].entityCount, &dimid);
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
    int varid;
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
      if (nodesets[i].dfCount != nodesets[i].entityCount) {
	status = EX_FATAL;
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: # dist fact (%"PRId64") not equal to # nodes (%"PRId64") "
		"in node set %"PRId64" file id %d",
		nodesets[i].dfCount, nodesets[i].entityCount, nodesets[i].id,
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
      }
      ex_compress_variable(exodusFilePtr, varid, 2);
    }
    if (nodesets[i].attributeCount > 0) {
      int numattrdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_NS(cur_num_node_sets+1),
			 nodesets[i].attributeCount, &numattrdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of attributes in nodeset %"PRId64" in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      dims[0] = dimid;
      dims[1] = numattrdim;
      status=nc_def_var(exodusFilePtr, VAR_NSATTRIB(cur_num_node_sets+1),
			nc_flt_code(exodusFilePtr), 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error:  failed to define attributes for element nodeset %"PRId64" in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_NSATTRIB(cur_num_node_sets+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for nodeset %"PRId64" in file id %d",
		nodesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }
  }
  return EX_NOERR;
}

// ========================================================================
int Internals::put_metadata(const std::vector<EdgeSet> &edgesets)
{
  const char *routine = "Internals::put_metadata(edgesets)";
  if (edgesets.empty())
    return EX_NOERR;

  char errmsg[MAX_ERR_LENGTH];
  int dims[2];

  int status  = 0; // clear error code

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  // Get number of edge sets defined for this file
  int dimid;
  int num_edge_sets = 0;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_ES, &dimid);
  if (status  != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      sprintf(errmsg,
	      "Error: no edge sets defined for file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    } else {
      sprintf(errmsg,
	      "Error: failed to locate edge sets defined in file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    }
    return (EX_FATAL);
  }

  // inquire how many edge sets are to be stored
  num_edge_sets = ex_inquire_int(exodusFilePtr, EX_INQ_EDGE_SETS);

  assert(static_cast<int>(edgesets.size()) == num_edge_sets);

  int namestrdim;
  status=nc_inq_dimid (exodusFilePtr, DIM_STR_NAME, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get string length in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  for (int i=0; i<num_edge_sets; i++) {

    //  NOTE: ex_inc_file_item is used to find the number of edge sets
    // for a specific file and returns that value incremented.
    int cur_num_edge_sets=(int)ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_EDGE_SET));

    if (edgesets[i].entityCount == 0)
      continue;

    status=nc_def_dim (exodusFilePtr, DIM_NUM_EDGE_ES(cur_num_edge_sets+1),
		       edgesets[i].entityCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: edge set %"PRId64" already defined in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to define number of edges for set %"PRId64" in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }


    // define variable to store edge set edge list here instead of in expns
    dims[0] = dimid;
    int varid;
    status=nc_def_var(exodusFilePtr,
		      VAR_EDGE_ES(cur_num_edge_sets+1),bulk_type,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: edge set %"PRId64" edge list already defined in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create edge set %"PRId64" edge list in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }

    // Orientation variable
    status = nc_def_var(exodusFilePtr,
			VAR_ORNT_ES(cur_num_edge_sets+1), bulk_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: extra list already exists for edge set %"PRId64" in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine,errmsg,exerrval);
      } else {
	sprintf(errmsg,
		"Error: failed to create extra list for edge set %"PRId64" in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine,errmsg,exerrval);
      }
      return(EX_FATAL);
    }

    // Create variable for distribution factors if required
    if (edgesets[i].dfCount > 0) {
      // num_dist_per_set should equal num_edges_per_set
      if (edgesets[i].dfCount != edgesets[i].entityCount) {
	status = EX_FATAL;
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: # dist fact (%"PRId64") not equal to # edges (%"PRId64") "
		"in edge set %"PRId64" file id %d",
		edgesets[i].dfCount, edgesets[i].entityCount, edgesets[i].id,
		exodusFilePtr);
	ex_err(routine, errmsg, status);
	return(EX_FATAL);
      } else {
	// create variable for distribution factors
	status=nc_def_var(exodusFilePtr, VAR_FACT_ES(cur_num_edge_sets+1),
			  nc_flt_code(exodusFilePtr), 1, dims, &varid);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  if (status == NC_ENAMEINUSE) {
	    sprintf(errmsg,
		    "Error: edge set %"PRId64" dist factors already exist in file id %d",
		    edgesets[i].id,exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  } else {
	    sprintf(errmsg,
		    "Error: failed to create edge set %"PRId64" dist factors in file id %d",
		    edgesets[i].id, exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  }
	  return(EX_FATAL);
	}
      }
    }
    if (edgesets[i].attributeCount > 0) {
      int numattrdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_ES(cur_num_edge_sets+1),
			 edgesets[i].attributeCount, &numattrdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of attributes in edgeset %"PRId64" in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      dims[0] = dimid;
      dims[1] = numattrdim;

      status=nc_def_var(exodusFilePtr, VAR_ESATTRIB(cur_num_edge_sets+1),
			nc_flt_code(exodusFilePtr), 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error:  failed to define attributes for element edgeset %"PRId64" in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_ESATTRIB(cur_num_edge_sets+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for edgeset %"PRId64" in file id %d",
		edgesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }
  }
  return EX_NOERR;
}

// ========================================================================
int Internals::put_metadata(const std::vector<FaceSet> &facesets)
{
  const char *routine = "Internals::put_metadata(facesets)";
  if (facesets.empty())
    return EX_NOERR;

  char errmsg[MAX_ERR_LENGTH];
  int dims[2];

  int status  = 0; // clear error code

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  // Get number of face sets defined for this file
  int dimid;
  int num_face_sets = 0;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_FS, &dimid);
  if (status  != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      sprintf(errmsg,
	      "Error: no face sets defined for file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    } else {
      sprintf(errmsg,
	      "Error: failed to locate face sets defined in file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    }
    return (EX_FATAL);
  }

  // inquire how many face sets are to be stored
  num_face_sets = ex_inquire_int(exodusFilePtr, EX_INQ_FACE_SETS);

  assert(static_cast<int>(facesets.size()) == num_face_sets);

  int namestrdim;
  status=nc_inq_dimid (exodusFilePtr, DIM_STR_NAME, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get string length in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  for (int i=0; i<num_face_sets; i++) {

    //  NOTE: ex_inc_file_item is used to find the number of face sets
    // for a specific file and returns that value incremented.
    int cur_num_face_sets=(int)ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_FACE_SET));

    if (facesets[i].entityCount == 0)
      continue;

    status=nc_def_dim (exodusFilePtr, DIM_NUM_FACE_FS(cur_num_face_sets+1),
		       facesets[i].entityCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: face set %"PRId64" already defined in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to define number of faces for set %"PRId64" in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }


    // define variable to store face set face list here instead of in expns
    dims[0] = dimid;
    int varid;
    status=nc_def_var(exodusFilePtr,
		      VAR_FACE_FS(cur_num_face_sets+1),bulk_type,1,dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: face set %"PRId64" face list already defined in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create face set %"PRId64" face list in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }

    // Orientation variable
    status = nc_def_var(exodusFilePtr,
			VAR_ORNT_FS(cur_num_face_sets+1), bulk_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: extra list already exists for face set %"PRId64" in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine,errmsg,exerrval);
      } else {
	sprintf(errmsg,
		"Error: failed to create extra list for face set %"PRId64" in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine,errmsg,exerrval);
      }
      return(EX_FATAL);
    }

    // Create variable for distribution factors if required
    if (facesets[i].dfCount > 0) {
      // num_dist_per_set should equal num_faces_per_set
      if (facesets[i].dfCount != facesets[i].entityCount) {
	status = EX_FATAL;
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: # dist fact (%"PRId64") not equal to # faces (%"PRId64") "
		"in face set %"PRId64" file id %d",
		facesets[i].dfCount, facesets[i].entityCount, facesets[i].id,
		exodusFilePtr);
	ex_err(routine, errmsg, status);
	return(EX_FATAL);
      } else {
	// create variable for distribution factors
	status=nc_def_var(exodusFilePtr, VAR_FACT_FS(cur_num_face_sets+1),
			  nc_flt_code(exodusFilePtr), 1, dims, &varid);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  if (status == NC_ENAMEINUSE) {
	    sprintf(errmsg,
		    "Error: face set %"PRId64" dist factors already exist in file id %d",
		    facesets[i].id,exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  } else {
	    sprintf(errmsg,
		    "Error: failed to create face set %"PRId64" dist factors in file id %d",
		    facesets[i].id, exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  }
	  return(EX_FATAL);
	}
      }
    }
    if (facesets[i].attributeCount > 0) {
      int numattrdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_FS(cur_num_face_sets+1),
			 facesets[i].attributeCount, &numattrdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of attributes in faceset %"PRId64" in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      dims[0] = dimid;
      dims[1] = numattrdim;
      status=nc_def_var(exodusFilePtr, VAR_FSATTRIB(cur_num_face_sets+1),
			nc_flt_code(exodusFilePtr), 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error:  failed to define attributes for element faceset %"PRId64" in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_FSATTRIB(cur_num_face_sets+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for faceset %"PRId64" in file id %d",
		facesets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }
    }
  }
  return EX_NOERR;
}

// ========================================================================
int Internals::put_metadata(const std::vector<ElemSet> &elemsets)
{
  const char *routine = "Internals::put_metadata(elemsets)";
  if (elemsets.empty())
    return EX_NOERR;

  char errmsg[MAX_ERR_LENGTH];
  int dims[2];

  int status  = 0; // clear error code

  // Get number of element sets defined for this file
  int dimid;
  int num_elem_sets = 0;
  status=nc_inq_dimid (exodusFilePtr, DIM_NUM_ELS, &dimid);
  if (status  != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      sprintf(errmsg,
	      "Error: no element sets defined for file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    } else {
      sprintf(errmsg,
	      "Error: failed to locate element sets defined in file id %d", exodusFilePtr);
      ex_err(routine, errmsg, status);
    }
    return (EX_FATAL);
  }

  // inquire how many element sets are to be stored
  num_elem_sets = ex_inquire_int(exodusFilePtr, EX_INQ_ELEM_SETS);

  assert(static_cast<int>(elemsets.size()) == num_elem_sets);

  int namestrdim;
  status=nc_inq_dimid (exodusFilePtr, DIM_STR_NAME, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    sprintf(errmsg,
	    "Error: failed to get string length in file id %d",exodusFilePtr);
    ex_err(routine,errmsg,status);
    return(EX_FATAL);
  }

  for (int i=0; i<num_elem_sets; i++) {

    //  NOTE: ex_inc_file_item is used to find the number of elem sets
    // for a specific file and returns that value incremented.
    int cur_num_elem_sets=(int)ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_ELEM_SET));

    if (elemsets[i].entityCount == 0)
      continue;

    status=nc_def_dim (exodusFilePtr, DIM_NUM_ELE_ELS(cur_num_elem_sets+1),
		       elemsets[i].entityCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: elem set %"PRId64" already defined in file id %d",
		elemsets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to define number of elems for set %"PRId64" in file id %d",
		elemsets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }


    // define variable to store element set element list here instead of in expns
    dims[0] = dimid;
    int varid;
    status=nc_def_var(exodusFilePtr,
		      VAR_ELEM_ELS(cur_num_elem_sets+1),
		      get_type(exodusFilePtr, EX_BULK_INT64_DB),1,dims, &varid);

    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
	sprintf(errmsg,
		"Error: element set %"PRId64" element list already defined in file id %d",
		elemsets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      } else {
	sprintf(errmsg,
		"Error: failed to create element set %"PRId64" element list in file id %d",
		elemsets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
      }
      return(EX_FATAL);
    }

    // Create variable for distribution factors if required
    if (elemsets[i].dfCount > 0) {
      // num_dist_per_set should equal num_elems_per_set
      if (elemsets[i].dfCount != elemsets[i].entityCount) {
	status = EX_FATAL;
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: # dist fact (%"PRId64") not equal to # elements (%"PRId64") "
		"in element set %"PRId64" file id %d",
		elemsets[i].dfCount, elemsets[i].entityCount, elemsets[i].id,
		exodusFilePtr);
	ex_err(routine, errmsg, status);
	return(EX_FATAL);
      } else {
	// create variable for distribution factors
	status=nc_def_var(exodusFilePtr, VAR_FACT_ELS(cur_num_elem_sets+1),
			  nc_flt_code(exodusFilePtr), 1, dims, &varid);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  if (status == NC_ENAMEINUSE) {
	    sprintf(errmsg,
		    "Error: element set %"PRId64" dist factors already exist in file id %d",
		    elemsets[i].id,exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  } else {
	    sprintf(errmsg,
		    "Error: failed to create element set %"PRId64" dist factors in file id %d",
		    elemsets[i].id, exodusFilePtr);
	    ex_err(routine, errmsg, status);
	  }
	  return(EX_FATAL);
	}
      }
    }
    if (elemsets[i].attributeCount > 0) {
      int numattrdim;
      status=nc_def_dim (exodusFilePtr, DIM_NUM_ATT_IN_NS(cur_num_elem_sets+1),
			 elemsets[i].attributeCount, &numattrdim);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define number of attributes in elemset %"PRId64" in file id %d",
		elemsets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      dims[0] = dimid;
      dims[1] = numattrdim;
      status=nc_def_var(exodusFilePtr, VAR_NSATTRIB(cur_num_elem_sets+1),
			nc_flt_code(exodusFilePtr), 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error:  failed to define attributes for element elemset %"PRId64" in file id %d",
		elemsets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
      }

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_NSATTRIB(cur_num_elem_sets+1),
			  NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define attribute name array for elemset %"PRId64" in file id %d",
		elemsets[i].id,exodusFilePtr);
	ex_err(routine, errmsg, status);
	return (EX_FATAL);
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
  std::vector<entity_id> nodeset_id(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    nodeset_id[i] = nodesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_NS_IDS, nodeset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    status[i] = nodesets[i].entityCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_NS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

int Internals::put_non_define_data(const std::vector<EdgeSet> &edgesets)
{
  if (edgesets.empty())
    return EX_NOERR;

  // Output edgeset ids...
  size_t num_edgesets = edgesets.size();
  std::vector<entity_id> edgeset_id(num_edgesets);
  for (size_t i = 0; i < num_edgesets; i++) {
    edgeset_id[i] = edgesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_ES_IDS, edgeset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_edgesets);
  for (size_t i = 0; i < num_edgesets; i++) {
    status[i] = edgesets[i].entityCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_ES_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

int Internals::put_non_define_data(const std::vector<FaceSet> &facesets)
{
  if (facesets.empty())
    return EX_NOERR;

  // Output faceset ids...
  size_t num_facesets = facesets.size();
  std::vector<entity_id> faceset_id(num_facesets);
  for (size_t i = 0; i < num_facesets; i++) {
    faceset_id[i] = facesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_FS_IDS, faceset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_facesets);
  for (size_t i = 0; i < num_facesets; i++) {
    status[i] = facesets[i].entityCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_FS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

int Internals::put_non_define_data(const std::vector<ElemSet> &elemsets)
{
  if (elemsets.empty())
    return EX_NOERR;

  // Output elemset ids...
  size_t num_elemsets = elemsets.size();
  std::vector<entity_id> elemset_id(num_elemsets);
  for (size_t i = 0; i < num_elemsets; i++) {
    elemset_id[i] = elemsets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_ELS_IDS, elemset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_elemsets);
  for (size_t i = 0; i < num_elemsets; i++) {
    status[i] = elemsets[i].entityCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_ELS_STAT, status) != NC_NOERR)
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

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);
  int status  = 0; // clear error code

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
  num_side_sets = ex_inquire_int(exodusFilePtr, EX_INQ_SIDE_SETS);

  assert(static_cast<int>(sidesets.size()) == num_side_sets);

  for (int i=0; i<num_side_sets; i++) {

    //  NOTE: ex_inc_file_item is used to find the number of side sets
    // for a specific file and returns that value incremented.
    int cur_num_side_sets = (int)ex_inc_file_item(exodusFilePtr, ex_get_counter_list(EX_SIDE_SET));

    if (sidesets[i].entityCount == 0)
      continue;

    status=nc_def_dim (exodusFilePtr, DIM_NUM_SIDE_SS(cur_num_side_sets+1),
		       sidesets[i].entityCount, &dimid);
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
    }
    ex_compress_variable(exodusFilePtr, varid, 2);
  }
  return EX_NOERR;
}

int Internals::put_non_define_data(const std::vector<SideSet> &sidesets)
{
  if (sidesets.empty())
    return EX_NOERR;

  // Output sideset ids...
  int num_sidesets = (int)sidesets.size();
  std::vector<entity_id> sideset_id(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    sideset_id[i] = sidesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_SS_IDS, sideset_id) != NC_NOERR)
    return(EX_FATAL);

  // Now, write the status array
  std::vector<int> status(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    status[i] = sidesets[i].entityCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_SS_STAT, status) != NC_NOERR)
    return(EX_FATAL);

  return (EX_NOERR);
}

namespace {
  template <typename T>
  int output_names(const std::vector<T> &entities, int exoid, ex_entity_type ent_type)
  {
    if (!entities.empty()) {
      std::vector<char*> names;
      
      names.resize(entities.size());
      for (size_t i=0; i < entities.size(); i++) {
	names[i] = (char*)entities[i].name.c_str();
      }
      return ex_put_names(exoid, ent_type, &names[0]);
    }
    return EX_NOERR;
  }

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
    ex_compress_variable(exodusFilePtr, *varid, 1);
    return EX_NOERR;
  }

  int define_variable(int exodusFilePtr, int64_t size, const char *dim, const char *var, nc_type type)
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
      ex_compress_variable(exodusFilePtr, varid, 1);
    } 
    return EX_NOERR;
  }

  int define_variables(int exodusFilePtr, int64_t size, const char *dim, const char *var[], const nc_type *types)
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
	status=nc_def_var(exodusFilePtr, var[i], types[i], 1, dimid, &varid);
	if (status != NC_NOERR) {
	  ex_opts(EX_VERBOSE);
	  sprintf(errmsg,
		  "Error: failed to define variable \"%s\" in file ID %d",
		  VAR_ELEM_MAP_BOR, exodusFilePtr);
	  ex_err(routine, errmsg, status);
	  return (EX_FATAL);
	}
	ex_compress_variable(exodusFilePtr, varid, 1);
	i++;
      }
    }
    return EX_NOERR;
  }

  int put_int_array(int exoid, const char *var_type, const std::vector<int> &array)
  {
    const char *routine = "Iopx_Internals.C, put_int_array";
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

  int put_id_array(int exoid, const char *var_type, const std::vector<entity_id> &ids)
  {
    const char *routine = "Iopx_Internals.C, put_id_array";
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
      // Have entity_id (long long), need ints...
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
  
  int define_coordinate_vars(int exodusFilePtr, int64_t nodes, int node_dim, int dimension, int dim_dim, int str_dim)
  {
    const char *routine = "Iopx_Internals.C, define_coordinate_vars";
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
	  ex_compress_variable(exodusFilePtr, varid, 1);
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
