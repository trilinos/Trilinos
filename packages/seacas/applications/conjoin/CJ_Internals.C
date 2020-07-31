/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <CJ_Internals.h>

#define NO_NETCDF_2
#include <exodusII.h>
#include <exodusII_int.h>
#include <fmt/format.h>
#include <netcdf.h>

#include <CJ_ExodusEntity.h>

#include <algorithm>
#include <copy_string_cpp.h>
#include <cstring>
#include <smart_assert.h>

using entity_id = int64_t;

namespace {
  nc_type get_type(int exoid, unsigned int type)
  {
    if ((ex_int64_status(exoid) & type) != 0U) {
      return NC_INT64;
    }

    return NC_INT;
  }
  bool lessOffset(const Excn::Block &b1, const Excn::Block &b2) { return b1.offset_ < b2.offset_; }

  int put_int_array(int exoid, const char *var_type, const std::vector<int> &array);
  int put_id_array(int exoid, const char *var_type, const std::vector<entity_id> &ids);

  int define_coordinate_vars(int exodusFilePtr, size_t nodes, int node_dim, int dimension,
                             int dim_dim, int str_dim);
} // namespace

Excn::Redefine::Redefine(int exoid) : exodusFilePtr(exoid)
{
  // Enter define mode...
  int status = nc_redef(exodusFilePtr);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    std::string errmsg;
    errmsg = fmt::format("Error: failed to put file id %d into define mode", exodusFilePtr);
    ex_err_fn(exoid, __func__, errmsg.c_str(), status);
    exit(EXIT_FAILURE);
  }
}

Excn::Redefine::~Redefine()
{
  try {
    int status = nc_enddef(exodusFilePtr);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      std::string errmsg;
      errmsg = fmt::format("Error: failed to complete variable definitions in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      exit(EXIT_FAILURE);
    }
  }
  catch (...) {
  }
}

Excn::Internals::Internals(int exoid, int maximum_name_length)
    : exodusFilePtr(exoid), maximumNameLength(maximum_name_length)
{
}

template <typename INT>
int Excn::Internals::write_meta_data(const Mesh<INT> &mesh, const std::vector<Block> &blocks,
                                     const std::vector<NodeSet<INT>> &nodesets,
                                     const std::vector<SideSet<INT>> &sidesets,
                                     const CommunicationMetaData &    comm)
{
  SMART_ASSERT(blocks.size() == mesh.blockCount);
  SMART_ASSERT(nodesets.size() == mesh.nodesetCount);
  SMART_ASSERT(sidesets.size() == mesh.sidesetCount);

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
  for (size_t i = 0; i < blocks.size(); i++) {
    const_cast<Excn::Block &>(blocks[i]).position_ = i;
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
  size_t last_offset = 0;
  bool   order_ok    = true;
  for (auto &block : blocks) {
    if (block.elementCount > 0) {
      if (block.offset_ < last_offset) {
        order_ok = false;
        break;
      }
      last_offset = block.offset_;
    }
  }

  std::vector<Block> sorted_blocks(blocks);
  if (!order_ok) {
    std::sort(sorted_blocks.begin(), sorted_blocks.end(), lessOffset);

    // Now, update the position_ field based on the sorted order.
    for (size_t i = 0; i < blocks.size(); i++) {
      size_t orig_position                                       = sorted_blocks[i].position_;
      const_cast<Excn::Block &>(blocks[orig_position]).position_ = i;
      SMART_ASSERT(blocks[orig_position].id == sorted_blocks[i].id);
    }
  }

  int ierr;
  {
    Excn::Redefine the_database(exodusFilePtr);

    ierr = put_metadata(mesh, comm);
    if (ierr != EX_NOERR) {
      return (ierr);
    }

    ierr = put_metadata(sorted_blocks);
    if (ierr != EX_NOERR) {
      return (ierr);
    }

    ierr = put_metadata(nodesets);
    if (ierr != EX_NOERR) {
      return (ierr);
    }

    ierr = put_metadata(sidesets);
    if (ierr != EX_NOERR) {
      return (ierr);
    }
  }

  // NON-Define mode output...
  ierr = put_non_define_data(mesh, comm);
  if (ierr != EX_NOERR) {
    return (ierr);
  }

  ierr = put_non_define_data(sorted_blocks);
  if (ierr != EX_NOERR) {
    return (ierr);
  }

  ierr = put_non_define_data(nodesets);
  if (ierr != EX_NOERR) {
    return (ierr);
  }

  ierr = put_non_define_data(sidesets);
  if (ierr != EX_NOERR) {
    return (ierr);
  }

  // For now, put entity names using the ExodusII api...
  {
    size_t max_entity = mesh.blockCount;
    if (mesh.nodesetCount > max_entity) {
      max_entity = mesh.nodesetCount;
    }
    if (mesh.sidesetCount > max_entity) {
      max_entity = mesh.sidesetCount;
    }
    if (mesh.blockCount > 0) {
      for (size_t i = 0; i < mesh.blockCount; i++) {
        if (blocks[i].attributeCount > max_entity) {
          max_entity = blocks[i].attributeCount;
        }
      }
    }

    size_t name_size = ex_inquire_int(exodusFilePtr, EX_INQ_MAX_READ_NAME_LENGTH);
    auto   names     = new char *[max_entity];
    for (size_t i = 0; i < max_entity; i++) {
      names[i] = new char[name_size + 1];
    }

    if (mesh.blockCount > 0) {
      for (size_t i = 0; i < mesh.blockCount; i++) {
        copy_string(names[i], sorted_blocks[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_ELEM_BLOCK, names);

      for (size_t i = 0; i < mesh.blockCount; i++) {
        if (blocks[i].attributeCount > 0) {
          SMART_ASSERT(blocks[i].attributeCount == blocks[i].attributeNames.size());
          for (size_t j = 0; j < blocks[i].attributeCount; j++) {
            std::memset(names[j], '\0', name_size + 1);
            if (!blocks[i].attributeNames[j].empty()) {
              copy_string(names[j], blocks[i].attributeNames[j], name_size + 1);
            }
          }
          ierr = ex_put_attr_names(exodusFilePtr, EX_ELEM_BLOCK, blocks[i].id, names);
          SMART_ASSERT(ierr == 0);
        }
      }
    }

    if (mesh.nodesetCount > 0) {
      for (size_t i = 0; i < mesh.nodesetCount; i++) {
        copy_string(names[i], nodesets[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_NODE_SET, names);
    }

    if (mesh.sidesetCount > 0) {
      for (size_t i = 0; i < mesh.sidesetCount; i++) {
        copy_string(names[i], sidesets[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_SIDE_SET, names);
    }

    for (size_t i = 0; i < max_entity; i++) {
      delete[] names[i];
    }
    delete[] names;
  }

  ex_update(exodusFilePtr);
  return (EX_NOERR);
}

template <typename INT>
int Excn::Internals::put_metadata(const Mesh<INT> &mesh, const CommunicationMetaData & /*unused*/)
{
  int numdimdim  = 0;
  int numnoddim  = 0;
  int numelemdim = 0;
  int timedim    = 0;
  int elblkdim   = 0;
  int strdim     = 0;
  int namestrdim = 0;
  int dim[2];
  int dimid = 0;
  int varid = 0;

  int map_type = get_type(exodusFilePtr, EX_MAPS_INT64_DB);
  int ids_type = get_type(exodusFilePtr, EX_IDS_INT64_DB);

  std::string errmsg;

  // define some attributes...
  int status = nc_put_att_text(exodusFilePtr, NC_GLOBAL, ATT_TITLE, mesh.title.length() + 1,
                               mesh.title.c_str());
  if (status != NC_NOERR) {
    errmsg = fmt::format("Error: failed to define title attribute to file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  // For use later to help readers know how much memory to allocate
  // for name storage, we define an attribute containing the maximum
  // size of any name.
  {
    int current_len = 0;
    status = nc_put_att_int(exodusFilePtr, NC_GLOBAL, ATT_MAX_NAME_LENGTH, NC_INT, 1, &current_len);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define ATT_MAX_NAME_LENGTH attribute to file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
  }

  // inquire previously defined dimensions
  status = nc_inq_dimid(exodusFilePtr, DIM_STR, &strdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg = fmt::format("Error: failed to get string length in file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  // create name string length dimension
  if (maximumNameLength < 32) {
    maximumNameLength = 32;
  }
  status = nc_def_dim(exodusFilePtr, DIM_STR_NAME, maximumNameLength + 1, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg = fmt::format("Error: failed to define name string length in file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  // ...and some dimensions..
  status = nc_def_dim(exodusFilePtr, DIM_NUM_DIM, mesh.dimensionality, &numdimdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg =
        fmt::format("Error: failed to define number of dimensions in file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  if ((status = nc_def_dim(exodusFilePtr, DIM_TIME, NC_UNLIMITED, &timedim)) != NC_NOERR) {
    errmsg = fmt::format("Error: failed to define time dimension in file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  {
    dim[0] = timedim;
    if ((status = nc_def_var(exodusFilePtr, VAR_WHOLE_TIME, nc_flt_code(exodusFilePtr), 1, dim,
                             &varid)) != NC_NOERR) {
      errmsg = fmt::format("Error: failed to define whole time step variable in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
  }

  if (mesh.nodeCount > 0) {
    status = nc_def_dim(exodusFilePtr, DIM_NUM_NODES, mesh.nodeCount, &numnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define number of nodes in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // Define the node map here to avoid a later redefine call
    int dims[1];
    dims[0] = numnoddim;
    status  = nc_def_var(exodusFilePtr, VAR_NODE_NUM_MAP, map_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg =
            fmt::format("Error: node numbering map already exists in file id %d", exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create node numbering map array in file id %d",
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }
    ex__compress_variable(exodusFilePtr, varid, 1);
  }

  if (mesh.elementCount > 0) {
    status = nc_def_dim(exodusFilePtr, DIM_NUM_ELEM, mesh.elementCount, &numelemdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of elements in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // Define the element map here to avoid a later redefine call
    int dims[1];
    dims[0] = numelemdim;
    varid   = 0;
    status  = nc_def_var(exodusFilePtr, VAR_ELEM_NUM_MAP, map_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg =
            fmt::format("Error: element numbering map already exists in file id %d", exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create element numbering map in file id %d",
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }
    ex__compress_variable(exodusFilePtr, varid, 1);
  }

  if (mesh.blockCount > 0) {
    status = nc_def_dim(exodusFilePtr, DIM_NUM_EL_BLK, mesh.blockCount, &elblkdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define number of element blocks in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // ...and some variables
    // element block id status array...
    dim[0] = elblkdim;
    status = nc_def_var(exodusFilePtr, VAR_STAT_EL_BLK, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define element block status array in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // element block id array
    status = nc_def_var(exodusFilePtr, VAR_ID_EL_BLK, ids_type, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define element block id array in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // store property name as attribute of property array variable
    status = nc_put_att_text(exodusFilePtr, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store element block property name %s in file id %d",
                           "ID", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    dim[0] = elblkdim;
    dim[1] = namestrdim;

    status = nc_def_var(exodusFilePtr, VAR_NAME_EL_BLK, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define element block name array in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
  }

  // node set id array:
  if (mesh.nodesetCount > 0) {

    status = nc_def_dim(exodusFilePtr, DIM_NUM_NS, mesh.nodesetCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of node sets in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // node set id status array:
    dim[0] = dimid;
    status = nc_def_var(exodusFilePtr, VAR_NS_STAT, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to create node sets status array in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // node set id array:
    dim[0] = dimid;
    status = nc_def_var(exodusFilePtr, VAR_NS_IDS, ids_type, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to create node sets property array in file id %d",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // store property name as attribute of property array variable
    status = nc_put_att_text(exodusFilePtr, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store node set property name %s in file id %d", "ID",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    dim[0] = dimid;
    dim[1] = namestrdim;

    status = nc_def_var(exodusFilePtr, VAR_NAME_NS, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define nodeset name array in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
  }

  // side set id array:
  if (mesh.sidesetCount > 0) {

    status = nc_def_dim(exodusFilePtr, DIM_NUM_SS, mesh.sidesetCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of side sets in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // side set id status array:
    dim[0] = dimid;
    status = nc_def_var(exodusFilePtr, VAR_SS_STAT, NC_INT, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define side set status in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // side set id array:
    status = nc_def_var(exodusFilePtr, VAR_SS_IDS, ids_type, 1, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define side set property in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // store property name as attribute of property array variable
    status = nc_put_att_text(exodusFilePtr, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store side set property name %s in file id %d", "ID",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    dim[0] = dimid;
    dim[1] = namestrdim;

    status = nc_def_var(exodusFilePtr, VAR_NAME_SS, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define sideset name array in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
  }

  status = define_coordinate_vars(exodusFilePtr, mesh.nodeCount, numnoddim, mesh.dimensionality,
                                  numdimdim, namestrdim);
  if (status != EX_NOERR) {
    return EX_FATAL;
  }

  return (EX_NOERR);
}

int Excn::Internals::put_metadata(const std::vector<Block> &blocks)
{
  std::string errmsg;
  int         dims[2];

  int status = 0; // clear error code

  if (blocks.empty()) {
    return (EX_NOERR);
  }

  // Get number of element blocks defined for this file
  int    dimid;
  size_t num_elem_blk = 0;
  status              = nc_inq_dimid(exodusFilePtr, DIM_NUM_EL_BLK, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg = fmt::format("Error: no element blocks defined in file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  int namestrdim;
  status = nc_inq_dimid(exodusFilePtr, DIM_STR_NAME, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg = fmt::format("Error: failed to get name string length in file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  status = nc_inq_dimlen(exodusFilePtr, dimid, &num_elem_blk);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg =
        fmt::format("Error: failed to get number of element blocks in file id %d", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return (EX_FATAL);
  }

  SMART_ASSERT(blocks.size() == num_elem_blk);

  // Iterate over element blocks ...
  for (size_t iblk = 0; iblk < num_elem_blk; iblk++) {

    ex__inc_file_item(exodusFilePtr, ex__get_counter_list(EX_ELEM_BLOCK));

    if (blocks[iblk].elementCount == 0) {
      continue;
    }

    // define some dimensions and variables
    int numelbdim;
    status = nc_def_dim(exodusFilePtr, DIM_NUM_EL_IN_BLK(iblk + 1), blocks[iblk].elementCount,
                        &numelbdim);
    if (status != NC_NOERR) {
      if (status == NC_ENAMEINUSE) { // duplicate entry
        ex_opts(EX_VERBOSE);
        errmsg = fmt::format("Error: element block {} already defined in file id %d",
                             blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        ex_opts(EX_VERBOSE);
        errmsg =
            fmt::format("Error: failed to define number of elements/block for block {} file id %d",
                        blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }

    int nelnoddim;
    status = nc_def_dim(exodusFilePtr, DIM_NUM_NOD_PER_EL(iblk + 1), blocks[iblk].nodesPerElement,
                        &nelnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of nodes/element for block {} in file id %d",
                      blocks[iblk].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    // element attribute array
    if (blocks[iblk].attributeCount > 0) {
      int numattrdim;
      status = nc_def_dim(exodusFilePtr, DIM_NUM_ATT_IN_BLK(iblk + 1), blocks[iblk].attributeCount,
                          &numattrdim);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        errmsg =
            fmt::format("Error: failed to define number of attributes in block {} in file id %d",
                        blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return (EX_FATAL);
      }

      dims[0]   = numelbdim;
      dims[1]   = numattrdim;
      int varid = 0;
      status = nc_def_var(exodusFilePtr, VAR_ATTRIB(iblk + 1), nc_flt_code(exodusFilePtr), 2, dims,
                          &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        errmsg =
            fmt::format("Error:  failed to define attributes for element block {} in file id %d",
                        blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return (EX_FATAL);
      }
      ex__compress_variable(exodusFilePtr, varid, 2);

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_ATTRIB(iblk + 1), NC_CHAR, 2, dims, &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        errmsg = fmt::format("Error: failed to define attribute name array for element block {}"
                             " in file id %d",
                             blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return (EX_FATAL);
      }
    }

    // element connectivity array
    dims[0] = numelbdim;
    dims[1] = nelnoddim;

    int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);
    int connid;
    status = nc_def_var(exodusFilePtr, VAR_CONN(iblk + 1), bulk_type, 2, dims, &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to create connectivity array for block {} in file id %d",
                           blocks[iblk].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
    ex__compress_variable(exodusFilePtr, connid, 1);

    // store element type as attribute of connectivity variable
    status = nc_put_att_text(exodusFilePtr, connid, ATT_NAME_ELB,
                             static_cast<int>(std::strlen(blocks[iblk].elType)) + 1,
                             blocks[iblk].elType);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store element type name %s in file id %d",
                           blocks[iblk].elType, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
  }
  return (EX_NOERR);
}

template <typename INT>
int Excn::Internals::put_non_define_data(const Mesh<INT> & /*unused*/,
                                         const CommunicationMetaData & /*unused*/)
{
  return EX_NOERR;
}

int Excn::Internals::put_non_define_data(const std::vector<Block> &blocks)
{
  int num_elem_blk = blocks.size(); // Verified via SMART_ASSERT earlier...

  if (num_elem_blk > 0) {
    // first get id of element block ids array variable
    std::vector<entity_id> elem_blk_id(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_id[iblk] = blocks[iblk].id;
    }

    if (put_id_array(exodusFilePtr, VAR_ID_EL_BLK, elem_blk_id) != NC_NOERR) {
      return (EX_FATAL);
    }

    // Now, write the element block status array
    std::vector<int> elem_blk_status(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_status[iblk] = blocks[iblk].elementCount > 0 ? 1 : 0;
    }

    if (put_int_array(exodusFilePtr, VAR_STAT_EL_BLK, elem_blk_status) != NC_NOERR) {
      return (EX_FATAL);
    }
  }
  return (EX_NOERR);
}

// ========================================================================
template <typename INT> int Excn::Internals::put_metadata(const std::vector<NodeSet<INT>> &nodesets)
{
  if (nodesets.empty()) {
    return EX_NOERR;
  }

  std::string errmsg;
  int         dims[2];

  int status = 0; // clear error code

  // Get number of node sets defined for this file
  int dimid;
  status = nc_inq_dimid(exodusFilePtr, DIM_NUM_NS, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      errmsg = fmt::format("Error: no node sets defined for file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    else {
      errmsg =
          fmt::format("Error: failed to locate node sets defined in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    return (EX_FATAL);
  }

  // inquire how many node sets are to be stored
  int num_node_sets = ex_inquire_int(exodusFilePtr, EX_INQ_NODE_SETS);
  SMART_ASSERT(static_cast<int>(nodesets.size()) == num_node_sets);

  for (int i = 0; i < num_node_sets; i++) {

    //  NOTE: ex__inc_file_item is used to find the number of node sets
    // for a specific file and returns that value incremented.
    int cur_num_node_sets = ex__inc_file_item(exodusFilePtr, ex__get_counter_list(EX_NODE_SET));

    if (nodesets[i].nodeCount == 0) {
      continue;
    }

    status = nc_def_dim(exodusFilePtr, DIM_NUM_NOD_NS(cur_num_node_sets + 1), nodesets[i].nodeCount,
                        &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: node set {} already defined in file id %d", nodesets[i].id,
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to define number of nodes for set {} in file id %d",
                             nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }

    // define variable to store node set node list here instead of in expns
    int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);
    dims[0]       = dimid;
    int varid;
    status =
        nc_def_var(exodusFilePtr, VAR_NODE_NS(cur_num_node_sets + 1), bulk_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: node set {} node list already defined in file id %d",
                             nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create node set {} node list in file id %d",
                             nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }
    ex__compress_variable(exodusFilePtr, varid, 1);

    // Create variable for distribution factors if required
    if (nodesets[i].dfCount > 0) {
      // num_dist_per_set should equal num_nodes_per_set
      if (nodesets[i].dfCount != nodesets[i].nodeCount) {
        status = EX_FATAL;
        ex_opts(EX_VERBOSE);
        errmsg = fmt::format("Error: # dist fact ({}) not equal to # nodes ({}) "
                             "in node set {} file id %d",
                             (int64_t)nodesets[i].dfCount, (int64_t)nodesets[i].nodeCount,
                             nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return (EX_FATAL);
      }

      // create variable for distribution factors
      status = nc_def_var(exodusFilePtr, VAR_FACT_NS(cur_num_node_sets + 1),
                          nc_flt_code(exodusFilePtr), 1, dims, &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg = fmt::format("Error: node set {} dist factors already exist in file id %d",
                               nodesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg = fmt::format("Error: failed to create node set {} dist factors in file id %d",
                               nodesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return (EX_FATAL);
      }
      ex__compress_variable(exodusFilePtr, varid, 2);
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals::put_non_define_data(const std::vector<NodeSet<INT>> &nodesets)
{
  if (nodesets.empty()) {
    return EX_NOERR;
  }

  // Output nodeset ids...
  int                    num_nodesets = nodesets.size();
  std::vector<entity_id> nodeset_id(num_nodesets);
  for (int i = 0; i < num_nodesets; i++) {
    nodeset_id[i] = nodesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_NS_IDS, nodeset_id) != NC_NOERR) {
    return (EX_FATAL);
  }

  // Now, write the status array
  std::vector<int> status(num_nodesets);
  for (int i = 0; i < num_nodesets; i++) {
    status[i] = nodesets[i].nodeCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_NS_STAT, status) != NC_NOERR) {
    return (EX_FATAL);
  }

  return (EX_NOERR);
}

// ========================================================================
template <typename INT> int Excn::Internals::put_metadata(const std::vector<SideSet<INT>> &sidesets)
{
  if (sidesets.empty()) {
    return EX_NOERR;
  }

  std::string errmsg;
  int         dims[2];

  int status = 0; // clear error code

  // Get number of side sets defined for this file
  int dimid;
  status = nc_inq_dimid(exodusFilePtr, DIM_NUM_SS, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      errmsg = fmt::format("Error: no side sets defined for file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    else {
      errmsg =
          fmt::format("Error: failed to locate side sets defined in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    return (EX_FATAL);
  }

  // inquire how many side sets are to be stored
  int num_side_sets = ex_inquire_int(exodusFilePtr, EX_INQ_SIDE_SETS);
  SMART_ASSERT(static_cast<int>(sidesets.size()) == num_side_sets);

  for (int i = 0; i < num_side_sets; i++) {

    //  NOTE: ex__inc_file_item is used to find the number of side sets
    // for a specific file and returns that value incremented.
    int cur_num_side_sets = ex__inc_file_item(exodusFilePtr, ex__get_counter_list(EX_SIDE_SET));

    if (sidesets[i].sideCount == 0) {
      continue;
    }

    status = nc_def_dim(exodusFilePtr, DIM_NUM_SIDE_SS(cur_num_side_sets + 1),
                        sidesets[i].sideCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: side set {} already defined in file id %d", sidesets[i].id,
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to define number of sides for set {} in file id %d",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }

    int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);
    dims[0]       = dimid;
    int varid     = 0;
    status =
        nc_def_var(exodusFilePtr, VAR_ELEM_SS(cur_num_side_sets + 1), bulk_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: side set {} element list already defined in file id %d",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create side set {} element list in file id %d",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }
    ex__compress_variable(exodusFilePtr, varid, 1);

    // create side list variable for side set
    status =
        nc_def_var(exodusFilePtr, VAR_SIDE_SS(cur_num_side_sets + 1), bulk_type, 1, dims, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: side list already exists for side set {} in file id %d",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create side list for side set {} in file id %d",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return (EX_FATAL);
    }
    ex__compress_variable(exodusFilePtr, varid, 1);

    // Create variable for distribution factors if required
    if (sidesets[i].dfCount > 0) {
      status = nc_def_dim(exodusFilePtr, DIM_NUM_DF_SS(cur_num_side_sets + 1), sidesets[i].dfCount,
                          &dimid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg = fmt::format("Error: side set df count {} already defined in file id %d",
                               sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg = fmt::format("Error: failed to define side set df count for set {} in file id %d",
                               sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return (EX_FATAL);
      }

      // create distribution factor list variable for side set
      dims[0] = dimid;
      status  = nc_def_var(exodusFilePtr, VAR_FACT_SS(cur_num_side_sets + 1),
                          nc_flt_code(exodusFilePtr), 1, dims, &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg =
              fmt::format("Error: dist factor list already exists for side set {} in file id %d",
                          sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg =
              fmt::format("Error: failed to create dist factor list for side set {} in file id %d",
                          sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return (EX_FATAL);
      }
      ex__compress_variable(exodusFilePtr, varid, 2);
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals::put_non_define_data(const std::vector<SideSet<INT>> &sidesets)
{
  if (sidesets.empty()) {
    return EX_NOERR;
  }

  // Output sideset ids...
  int                    num_sidesets = sidesets.size();
  std::vector<entity_id> sideset_id(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    sideset_id[i] = sidesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_SS_IDS, sideset_id) != NC_NOERR) {
    return (EX_FATAL);
  }

  // Now, write the status array
  std::vector<int> status(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    status[i] = sidesets[i].sideCount > 0 ? 1 : 0;
  }

  if (put_int_array(exodusFilePtr, VAR_SS_STAT, status) != NC_NOERR) {
    return (EX_FATAL);
  }

  return (EX_NOERR);
}

namespace {
  int put_int_array(int exoid, const char *var_type, const std::vector<int> &array)
  {
    std::string errmsg;
    int         var_id;
    int         status;

    status = nc_inq_varid(exoid, var_type, &var_id);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to locate %s in file id %d", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    status = nc_put_var_int(exoid, var_id, array.data());
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to write %s array in file id %d", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
    return (EX_NOERR);
  }

  int put_id_array(int exoid, const char *var_type, const std::vector<entity_id> &ids)
  {
    std::string errmsg;
    int         var_id;

    int status = nc_inq_varid(exoid, var_type, &var_id);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to locate %s in file id %d", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }

    int id_type = get_type(exoid, EX_IDS_INT64_API);

    if (id_type == NC_INT64) {
      status = nc_put_var_longlong(exoid, var_id, (long long int *)ids.data());
    }
    else {
      // Have entity_id (long long), need ints...
      std::vector<int> int_ids(ids.size());
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4244)
#endif
      int_ids.assign(ids.begin(), ids.end());
#ifdef _MSC_VER
#pragma warning(pop)
#endif
      status = nc_put_var_int(exoid, var_id, int_ids.data());
    }

    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to write %s array in file id %d", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
    return (EX_NOERR);
  }

  int define_coordinate_vars(int exodusFilePtr, size_t nodes, int node_dim, int dimension,
                             int dim_dim, int str_dim)
  {
    std::string errmsg;
    int         status;
    int         dim[2];
    int         varid;

    if (nodes > 0) {
      if (ex_large_model(exodusFilePtr) == 1) {
        // node coordinate arrays -- separate storage...

        dim[0] = node_dim;
        if (dimension > 0) {
          status =
              nc_def_var(exodusFilePtr, VAR_COORD_X, nc_flt_code(exodusFilePtr), 1, dim, &varid);
          if (status != NC_NOERR) {
            ex_opts(EX_VERBOSE);
            errmsg = fmt::format("Error: failed to define node x coordinate array in file id %d",
                                 exodusFilePtr);
            ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
            return (EX_FATAL);
          }
          ex__compress_variable(exodusFilePtr, varid, 1);
        }

        if (dimension > 1) {
          status =
              nc_def_var(exodusFilePtr, VAR_COORD_Y, nc_flt_code(exodusFilePtr), 1, dim, &varid);
          if (status != NC_NOERR) {
            ex_opts(EX_VERBOSE);
            errmsg = fmt::format("Error: failed to define node y coordinate array in file id %d",
                                 exodusFilePtr);
            ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
            return (EX_FATAL);
          }
          ex__compress_variable(exodusFilePtr, varid, 1);
        }

        if (dimension > 2) {
          status =
              nc_def_var(exodusFilePtr, VAR_COORD_Z, nc_flt_code(exodusFilePtr), 1, dim, &varid);
          if (status != NC_NOERR) {
            ex_opts(EX_VERBOSE);
            errmsg = fmt::format("Error: failed to define node z coordinate array in file id %d",
                                 exodusFilePtr);
            ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
            return (EX_FATAL);
          }
          ex__compress_variable(exodusFilePtr, varid, 1);
        }
      }
      else {
        // node coordinate arrays:  -- all stored together (old method)2
        dim[0] = dim_dim;
        dim[1] = node_dim;
        status = nc_def_var(exodusFilePtr, VAR_COORD, nc_flt_code(exodusFilePtr), 2, dim, &varid);
        if (status != NC_NOERR) {
          ex_opts(EX_VERBOSE);
          errmsg = fmt::format("Error: failed to define node coordinate array in file id %d",
                               exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
          return (EX_FATAL);
        }
      }
    }

    // coordinate names array
    dim[0] = dim_dim;
    dim[1] = str_dim;

    status = nc_def_var(exodusFilePtr, VAR_NAME_COOR, NC_CHAR, 2, dim, &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define coordinate name array in file id %d", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return (EX_FATAL);
    }
    return EX_NOERR;
  }
} // namespace

template int Excn::Internals::write_meta_data(const Excn::Mesh<int> &                mesh,
                                              const std::vector<Excn::Block> &       blocks,
                                              const std::vector<Excn::NodeSet<int>> &nodesets,
                                              const std::vector<Excn::SideSet<int>> &sidesets,
                                              const Excn::CommunicationMetaData &    comm);
template int Excn::Internals::write_meta_data(const Excn::Mesh<int64_t> &                mesh,
                                              const std::vector<Excn::Block> &           blocks,
                                              const std::vector<Excn::NodeSet<int64_t>> &nodesets,
                                              const std::vector<Excn::SideSet<int64_t>> &sidesets,
                                              const Excn::CommunicationMetaData &        comm);
