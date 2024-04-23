/*
 * Copyright(C) 1999-2024,  National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <EP_ExodusEntity.h> // for Block, Mesh
#include <EP_Internals.h>    // for Internals, Redefine

#include "vector_data.h"
#include <algorithm> // for sort
#include <array>
#include <copy_string_cpp.h>
#include <cstddef> // for size_t
#include <cstdint> // for int64_t
#include <cstdlib> // for exit, EXIT_FAILURE
#include <cstring> // for strlen, memset
#include <fmt/ostream.h>
#include <smart_assert.h> // for SMART_ASSERT
#include <string>         // for string, basic_string
#include <vector>         // for vector

extern "C" {
#define NO_NETCDF_2
#include <exodusII.h>
#include <exodusII_int.h>
#include <netcdf.h>
}

// Explicit Initialization of the functions used...
namespace Excn {
  template int Internals<int>::write_meta_data(const Mesh &mesh, const std::vector<Block> &blocks,
                                               const std::vector<NodeSet<int>>   &nodesets,
                                               const std::vector<SideSet<int>>   &sidesets,
                                               const std::vector<EdgeBlock<int>> &edgeblock,
                                               const std::vector<FaceBlock<int>> &faceblock,
                                               const CommunicationMetaData       &comm);

  template bool Internals<int>::check_meta_data(const Mesh &mesh, const std::vector<Block> &blocks,
                                                const std::vector<NodeSet<int>>   &nodesets,
                                                const std::vector<SideSet<int>>   &sidesets,
                                                const std::vector<EdgeBlock<int>> &edgeblock,
                                                const std::vector<FaceBlock<int>> &faceblock,
                                                const CommunicationMetaData       &comm);

  template int Internals<int64_t>::write_meta_data(const Mesh                            &mesh,
                                                   const std::vector<Block>              &blocks,
                                                   const std::vector<NodeSet<int64_t>>   &nodesets,
                                                   const std::vector<SideSet<int64_t>>   &sidesets,
                                                   const std::vector<EdgeBlock<int64_t>> &edgeblock,
                                                   const std::vector<FaceBlock<int64_t>> &faceblock,
                                                   const CommunicationMetaData           &comm);

  template bool Internals<int64_t>::check_meta_data(
      const Mesh &mesh, const std::vector<Block> &blocks,
      const std::vector<NodeSet<int64_t>> &nodesets, const std::vector<SideSet<int64_t>> &sidesets,
      const std::vector<EdgeBlock<int64_t>> &edgeblock,
      const std::vector<FaceBlock<int64_t>> &faceblock, const CommunicationMetaData &comm);
} // namespace Excn

namespace {
  bool lessOffset(const Excn::Block &b1, const Excn::Block &b2) { return b1.offset_ < b2.offset_; }

  nc_type get_type(int exoid, unsigned int type)
  {
    if ((ex_int64_status(exoid) & type) != 0U) {
      return NC_INT64;
    }

    return NC_INT;
  }

  int define_netcdf_vars(int exoid, const char *type, size_t count, const char *dim_num,
                         const char *stat_var, const char *id_var, const char *name_var);

  template <typename INT>
  int put_array(int exoid, const char *var_type, const std::vector<INT> &array);

  int put_id_array(int exoid, const char *var_type, const std::vector<ex_entity_id> &ids);

  int define_coordinate_vars(int exodusFilePtr, int64_t nodes, int node_dim, int dimension,
                             int dim_dim, int str_dim);
} // namespace

bool Excn::is_path_absolute(const std::string &path)
{
#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
  return path[0] == '\\' || path[1] == ':';
#else
  return path[0] == '/';
#endif
}

Excn::Redefine::Redefine(int exoid) : exodusFilePtr(exoid)
{
  // Enter define mode...
  int status = nc_redef(exodusFilePtr);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    std::string errmsg;
    errmsg = fmt::format("Error: failed to put file id {} into define mode", exodusFilePtr);
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
      errmsg = fmt::format("Error: failed to complete variable definitions in file id {}",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      exit(EXIT_FAILURE);
    }
  }
  catch (...) {
  }
}

template <typename INT>
int Excn::Internals<INT>::write_meta_data(const Mesh &mesh, const std::vector<Block> &blocks,
                                          const std::vector<NodeSet<INT>>   &nodesets,
                                          const std::vector<SideSet<INT>>   &sidesets,
                                          const std::vector<EdgeBlock<INT>> &edgeblocks,
                                          const std::vector<FaceBlock<INT>> &faceblocks,
                                          const CommunicationMetaData       &comm)
{
  SMART_ASSERT((int)blocks.size() == mesh.blockCount);
  SMART_ASSERT((int)nodesets.size() == mesh.nodesetCount);
  SMART_ASSERT((int)sidesets.size() == mesh.sidesetCount);
  SMART_ASSERT((int)edgeblocks.size() == mesh.edgeBlockCount);
  SMART_ASSERT((int)faceblocks.size() == mesh.faceBlockCount);

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
  int64_t last_offset = 0;
  bool    order_ok    = true;
  for (const auto &block : blocks) {
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
      int orig_position                                          = sorted_blocks[i].position_;
      const_cast<Excn::Block &>(blocks[orig_position]).position_ = i;
      SMART_ASSERT(blocks[orig_position].id == sorted_blocks[i].id);
    }
  }

  int ierr;
  {
    Excn::Redefine the_database(exodusFilePtr);

    ierr = put_metadata(mesh, comm);
    if (ierr != EX_NOERR) {
      return ierr;
    }

    ierr = put_metadata(sorted_blocks);
    if (ierr != EX_NOERR) {
      return ierr;
    }

    ierr = put_metadata(nodesets);
    if (ierr != EX_NOERR) {
      return ierr;
    }

    ierr = put_metadata(sidesets);
    if (ierr != EX_NOERR) {
      return ierr;
    }

    ierr = put_metadata(edgeblocks);
    if (ierr != EX_NOERR) {
      return ierr;
    }

    ierr = put_metadata(faceblocks);
    if (ierr != EX_NOERR) {
      return ierr;
    }
  }

  // NON-Define mode output...
  ierr = put_non_define_data(mesh, comm);
  if (ierr != EX_NOERR) {
    return ierr;
  }

  ierr = put_non_define_data(sorted_blocks);
  if (ierr != EX_NOERR) {
    return ierr;
  }

  ierr = put_non_define_data(nodesets);
  if (ierr != EX_NOERR) {
    return ierr;
  }

  ierr = put_non_define_data(sidesets);
  if (ierr != EX_NOERR) {
    return ierr;
  }

  ierr = put_non_define_data(edgeblocks);
  if (ierr != EX_NOERR) {
    return ierr;
  }

  ierr = put_non_define_data(faceblocks);
  if (ierr != EX_NOERR) {
    return ierr;
  }

  // For now, put entity names using the ExodusII api...
  {
    int max_entity = std::max({mesh.blockCount, mesh.nodesetCount, mesh.sidesetCount,
                               mesh.edgeBlockCount, mesh.faceBlockCount});
    for (const auto &block : blocks) {
      max_entity = std::max(max_entity, block.attributeCount);
    }

    size_t name_size = ex_inquire_int(exodusFilePtr, EX_INQ_MAX_READ_NAME_LENGTH);
    auto  *names     = new char *[max_entity];
    for (int i = 0; i < max_entity; i++) {
      names[i] = new char[name_size + 1];
    }

    if (mesh.blockCount > 0) {
      for (int i = 0; i < mesh.blockCount; i++) {
        copy_string(names[i], sorted_blocks[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_ELEM_BLOCK, names);

      for (int i = 0; i < mesh.blockCount; i++) {
        if (blocks[i].attributeCount > 0) {
          SMART_ASSERT((size_t)blocks[i].attributeCount == blocks[i].attributeNames.size());
          for (int j = 0; j < blocks[i].attributeCount; j++) {
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
      for (int i = 0; i < mesh.nodesetCount; i++) {
        copy_string(names[i], nodesets[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_NODE_SET, names);
    }

    if (mesh.sidesetCount > 0) {
      for (int i = 0; i < mesh.sidesetCount; i++) {
        copy_string(names[i], sidesets[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_SIDE_SET, names);
    }

    if (mesh.edgeBlockCount > 0) {
      for (int i = 0; i < mesh.edgeBlockCount; i++) {
        copy_string(names[i], edgeblocks[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_EDGE_BLOCK, names);
    }

    if (mesh.faceBlockCount > 0) {
      for (int i = 0; i < mesh.faceBlockCount; i++) {
        copy_string(names[i], faceblocks[i].name_, name_size + 1);
      }
      ex_put_names(exodusFilePtr, EX_FACE_BLOCK, names);
    }

    for (int i = 0; i < max_entity; i++) {
      delete[] names[i];
    }
    delete[] names;
  }

  ex_update(exodusFilePtr);
  return EX_NOERR;
}

template <typename INT>
bool Excn::Internals<INT>::check_meta_data(const Mesh &mesh, const std::vector<Block> & /*unused*/,
                                           const std::vector<NodeSet<INT>> & /*unused*/,
                                           const std::vector<SideSet<INT>> & /*unused*/,
                                           const std::vector<EdgeBlock<INT>> & /*unused*/,
                                           const std::vector<FaceBlock<INT>> & /*unused*/,
                                           const CommunicationMetaData & /*unused*/)
{
  ex_init_params init_data{};
  ex_get_init_ext(exodusFilePtr, &init_data);

  bool matches = true;
  if (mesh.dimensionality != init_data.num_dim) {
    fmt::print(stderr,
               "ERROR: (EPU) original mesh dimensionality ({}) does not match current "
               "dimensionality ({})\n",
               mesh.dimensionality, init_data.num_dim);
    matches = false;
  }

  if (mesh.nodeCount != init_data.num_nodes) {
    fmt::print(
        stderr,
        "ERROR: (EPU) original mesh node count ({}) does not match current node count ({})\n",
        mesh.nodeCount, init_data.num_nodes);
    matches = false;
  }

  if (mesh.elementCount != init_data.num_elem) {
    fmt::print(
        stderr,
        "ERROR: (EPU) original mesh element count ({}) does not match current element count ({})\n",
        mesh.elementCount, init_data.num_elem);
    matches = false;
  }

  if (mesh.blockCount != init_data.num_elem_blk) {
    fmt::print(stderr,
               "ERROR: (EPU) original mesh element block count ({}) does not match current element "
               "block count ({})\n",
               mesh.blockCount, init_data.num_elem_blk);
    matches = false;
  }

  if (mesh.nodesetCount != init_data.num_node_sets) {
    fmt::print(
        stderr,
        "ERROR: (EPU) original mesh nodeset count ({}) does not match current nodeset count ({})\n",
        mesh.nodesetCount, init_data.num_node_sets);
    matches = false;
  }

  if (mesh.sidesetCount != init_data.num_side_sets) {
    fmt::print(
        stderr,
        "ERROR: (EPU) original mesh sideset count ({}) does not match current sideset count ({})\n",
        mesh.sidesetCount, init_data.num_side_sets);
    matches = false;
  }

  if (mesh.edgeBlockCount != init_data.num_edge_blk) {
    fmt::print(stderr,
               "ERROR: (EPU) original mesh edgeblock count ({}) does not match current edgeblock "
               "count ({})\n",
               mesh.edgeBlockCount, init_data.num_edge_blk);
    matches = false;
  }

  if (mesh.faceBlockCount != init_data.num_face_blk) {
    fmt::print(stderr,
               "ERROR: (EPU) original mesh faceblock count ({}) does not match current faceblock "
               "count ({})\n",
               mesh.faceBlockCount, init_data.num_face_blk);
    matches = false;
  }

  return matches;
}

template <typename INT>
int Excn::Internals<INT>::put_metadata(const Mesh &mesh, const CommunicationMetaData & /*unused*/)
{
  int numdimdim  = 0;
  int timedim    = 0;
  int numnoddim  = 0;
  int namestrdim = 0;

  int map_type = get_type(exodusFilePtr, EX_MAPS_INT64_DB);

  std::string errmsg;

  // define some attributes...
  int status = nc_put_att_text(exodusFilePtr, NC_GLOBAL, ATT_TITLE, mesh.title.length() + 1,
                               mesh.title.c_str());
  if (status != NC_NOERR) {
    errmsg = fmt::format("Error: failed to define title attribute to file id {}", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return EX_FATAL;
  }

  // For use later to help readers know how much memory to allocate
  // for name storage, we define an attribute containing the maximum
  // size of any name.
  {
    int current_len = 0;
    status = nc_put_att_int(exodusFilePtr, NC_GLOBAL, ATT_MAX_NAME_LENGTH, NC_INT, 1, &current_len);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define ATT_MAX_NAME_LENGTH attribute to file id {}",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
  }

  // create name string length dimension
  if (maximumNameLength < 32) {
    maximumNameLength = 32;
  }
  status = nc_def_dim(exodusFilePtr, DIM_STR_NAME, maximumNameLength + 1, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg = fmt::format("Error: failed to define name string length in file id {}", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return EX_FATAL;
  }

  // ...and some dimensions..
  status = nc_def_dim(exodusFilePtr, DIM_NUM_DIM, mesh.dimensionality, &numdimdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg =
        fmt::format("Error: failed to define number of dimensions in file id {}", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return EX_FATAL;
  }

  if ((status = nc_def_dim(exodusFilePtr, DIM_TIME, NC_UNLIMITED, &timedim)) != NC_NOERR) {
    errmsg = fmt::format("Error: failed to define time dimension in file id {}", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return EX_FATAL;
  }

  {
    std::array dim{timedim};
    int        varid = 0;
    if ((status = nc_def_var(exodusFilePtr, VAR_WHOLE_TIME, nc_flt_code(exodusFilePtr), 1,
                             Data(dim), &varid)) != NC_NOERR) {
      errmsg = fmt::format("Error: failed to define whole time step variable in file id {}",
                           exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    struct exi_file_item *file = exi_find_file_item(exodusFilePtr);
    if (file) {
      file->time_varid = varid;
    }

    exi_compress_variable(exodusFilePtr, varid, -2); /* don't compress, but do set collective io */
  }

  if (mesh.nodeCount > 0) {
    status = nc_def_dim(exodusFilePtr, DIM_NUM_NODES, mesh.nodeCount, &numnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define number of nodes in file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // Define the node map here to avoid a later redefine call
    if (mesh.needNodeMap) {
      std::array dims{numnoddim};
      int        varid = 0;
      status = nc_def_var(exodusFilePtr, VAR_NODE_NUM_MAP, map_type, 1, Data(dims), &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg =
              fmt::format("Error: node numbering map already exists in file id {}", exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg = fmt::format("Error: failed to create node numbering map array in file id {}",
                               exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return EX_FATAL;
      }
      exi_compress_variable(exodusFilePtr, varid, 1);
    }
  }

  if (mesh.elementCount > 0) {
    int numelemdim = 0;
    status         = nc_def_dim(exodusFilePtr, DIM_NUM_ELEM, mesh.elementCount, &numelemdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of elements in file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // Define the element map here to avoid a later redefine call
    if (mesh.needElementMap) {
      std::array dims{numelemdim};
      int        varid = 0;
      status = nc_def_var(exodusFilePtr, VAR_ELEM_NUM_MAP, map_type, 1, Data(dims), &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg = fmt::format("Error: element numbering map already exists in file id {}",
                               exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg = fmt::format("Error: failed to create element numbering map in file id {}",
                               exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return EX_FATAL;
      }
      exi_compress_variable(exodusFilePtr, varid, 1);
    }
  }

  if (mesh.blockCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "element block", mesh.blockCount, DIM_NUM_EL_BLK,
                                VAR_STAT_EL_BLK, VAR_ID_EL_BLK, VAR_NAME_EL_BLK);
    if (status != EX_NOERR) {
      return EX_FATAL;
    }
  }

  // node set id array:
  if (mesh.nodesetCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "node set", mesh.nodesetCount, DIM_NUM_NS,
                                VAR_NS_STAT, VAR_NS_IDS, VAR_NAME_NS);
    if (status != EX_NOERR) {
      return EX_FATAL;
    }
  }

  // side set id array:
  if (mesh.sidesetCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "side set", mesh.sidesetCount, DIM_NUM_SS,
                                VAR_SS_STAT, VAR_SS_IDS, VAR_NAME_SS);
    if (status != EX_NOERR) {
      return EX_FATAL;
    }
  }

  // edge block id array:
  if (mesh.edgeBlockCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "edge block", mesh.edgeBlockCount, DIM_NUM_ED_BLK,
                                VAR_STAT_ED_BLK, VAR_ID_ED_BLK, VAR_NAME_ED_BLK);
    if (status != EX_NOERR) {
      return EX_FATAL;
    }
  }

  // face block id array:
  if (mesh.faceBlockCount > 0) {
    status = define_netcdf_vars(exodusFilePtr, "face block", mesh.faceBlockCount, DIM_NUM_FA_BLK,
                                VAR_STAT_FA_BLK, VAR_ID_FA_BLK, VAR_NAME_FA_BLK);
    if (status != EX_NOERR) {
      return EX_FATAL;
    }
  }

  status = define_coordinate_vars(exodusFilePtr, mesh.nodeCount, numnoddim, mesh.dimensionality,
                                  numdimdim, namestrdim);
  if (status != EX_NOERR) {
    return EX_FATAL;
  }

  return EX_NOERR;
}

template <typename INT> int Excn::Internals<INT>::put_metadata(const std::vector<Block> &blocks)
{
  std::string errmsg;

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status = 0; // clear error code

  if (blocks.empty()) {
    return EX_NOERR;
  }

  // Get number of element blocks defined for this file
  int    dimid;
  size_t num_elem_blk = 0;
  status              = nc_inq_dimid(exodusFilePtr, DIM_NUM_EL_BLK, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg = fmt::format("Error: no element blocks defined in file id {}", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return EX_FATAL;
  }

  int namestrdim;
  status = nc_inq_dimid(exodusFilePtr, DIM_STR_NAME, &namestrdim);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg = fmt::format("Error: failed to get name string length in file id {}", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return EX_FATAL;
  }

  status = nc_inq_dimlen(exodusFilePtr, dimid, &num_elem_blk);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    errmsg =
        fmt::format("Error: failed to get number of element blocks in file id {}", exodusFilePtr);
    ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    return EX_FATAL;
  }

  SMART_ASSERT(blocks.size() == num_elem_blk);

  // Iterate over element blocks ...
  for (size_t iblk = 0; iblk < num_elem_blk; iblk++) {

    exi_inc_file_item(exodusFilePtr, exi_get_counter_list(EX_ELEM_BLOCK));

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
        errmsg = fmt::format("Error: element block {} already defined in file id {}",
                             blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        ex_opts(EX_VERBOSE);
        errmsg =
            fmt::format("Error: failed to define number of elements/block for block {} file id {}",
                        blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }

    int nelnoddim;
    status = nc_def_dim(exodusFilePtr, DIM_NUM_NOD_PER_EL(iblk + 1), blocks[iblk].nodesPerElement,
                        &nelnoddim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of nodes/element for block {} in file id {}",
                      blocks[iblk].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // element attribute array
    if (blocks[iblk].attributeCount > 0) {
      int numattrdim;
      status = nc_def_dim(exodusFilePtr, DIM_NUM_ATT_IN_BLK(iblk + 1), blocks[iblk].attributeCount,
                          &numattrdim);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        errmsg =
            fmt::format("Error: failed to define number of attributes in block {} in file id {}",
                        blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return EX_FATAL;
      }

      std::array dims{numelbdim, numattrdim};
      int        varid = 0;
      status = nc_def_var(exodusFilePtr, VAR_ATTRIB(iblk + 1), nc_flt_code(exodusFilePtr), 2,
                          Data(dims), &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        errmsg =
            fmt::format("Error:  failed to define attributes for element block {} in file id {}",
                        blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return EX_FATAL;
      }
      exi_compress_variable(exodusFilePtr, varid, 2);

      // Attribute name array...
      dims[0] = numattrdim;
      dims[1] = namestrdim;

      status = nc_def_var(exodusFilePtr, VAR_NAME_ATTRIB(iblk + 1), NC_CHAR, 2, Data(dims), &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        errmsg = fmt::format(
            "Error: failed to define attribute name array for element block {} in file id {}",
            blocks[iblk].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return EX_FATAL;
      }
    }

    // element connectivity array
    std::array dims{numelbdim, nelnoddim};

    int connid;
    status = nc_def_var(exodusFilePtr, VAR_CONN(iblk + 1), bulk_type, 2, Data(dims), &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to create connectivity array for block {} in file id {}",
                           blocks[iblk].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
    exi_compress_variable(exodusFilePtr, connid, 1);

    // store element type as attribute of connectivity variable
    status = nc_put_att_text(exodusFilePtr, connid, ATT_NAME_ELB,
                             static_cast<int>(std::strlen(blocks[iblk].elType)) + 1,
                             blocks[iblk].elType);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store element type name {} in file id {}",
                           blocks[iblk].elType, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const Mesh & /*unused*/,
                                              const CommunicationMetaData & /*unused*/)
{
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<Block> &blocks)
{
  int num_elem_blk = blocks.size(); // Verified via SMART_ASSERT earlier...

  if (num_elem_blk > 0) {
    // first get id of element block ids array variable
    std::vector<ex_entity_id> elem_blk_id(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_id[iblk] = blocks[iblk].id;
    }

    if (put_id_array(exodusFilePtr, VAR_ID_EL_BLK, elem_blk_id) != NC_NOERR) {
      return EX_FATAL;
    }

    // Now, write the element block status array
    std::vector<int> elem_blk_status(num_elem_blk);
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      elem_blk_status[iblk] = blocks[iblk].elementCount > 0 ? 1 : 0;
    }

    if (put_array(exodusFilePtr, VAR_STAT_EL_BLK, elem_blk_status) != NC_NOERR) {
      return EX_FATAL;
    }
  }
  return EX_NOERR;
}

// ========================================================================
template <typename INT>
int Excn::Internals<INT>::put_metadata(const std::vector<NodeSet<INT>> &nodesets)
{
  if (nodesets.empty()) {
    return EX_NOERR;
  }

  std::string errmsg;

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status = 0; // clear error code

  // Get number of node sets defined for this file
  int dimid;
  status = nc_inq_dimid(exodusFilePtr, DIM_NUM_NS, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      errmsg = fmt::format("Error: no node sets defined for file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    else {
      errmsg =
          fmt::format("Error: failed to locate node sets defined in file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    return EX_FATAL;
  }

  // inquire how many node sets are to be stored
  int num_node_sets = ex_inquire_int(exodusFilePtr, EX_INQ_NODE_SETS);
  SMART_ASSERT(static_cast<int>(nodesets.size()) == num_node_sets);

  for (int i = 0; i < num_node_sets; i++) {

    //  NOTE: exi_inc_file_item is used to find the number of node sets
    // for a specific file and returns that value incremented.
    int cur_num_node_sets =
        (int)exi_inc_file_item(exodusFilePtr, exi_get_counter_list(EX_NODE_SET));

    if (nodesets[i].nodeCount == 0) {
      continue;
    }

    status = nc_def_dim(exodusFilePtr, DIM_NUM_NOD_NS(cur_num_node_sets + 1), nodesets[i].nodeCount,
                        &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: node set {} already defined in file id {}", nodesets[i].id,
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to define number of nodes for set {} in file id {}",
                             nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }

    // define variable to store node set node list here instead of in expns
    std::array dims{dimid};
    int        varid = 0;
    status = nc_def_var(exodusFilePtr, VAR_NODE_NS(cur_num_node_sets + 1), bulk_type, 1, Data(dims),
                        &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: node set {} node list already defined in file id {}",
                             nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create node set {} node list in file id {}",
                             nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }
    exi_compress_variable(exodusFilePtr, varid, 1);

    // Create variable for distribution factors if required
    if (nodesets[i].dfCount > 0) {
      // num_dist_per_set should equal num_nodes_per_set
      if (nodesets[i].dfCount != nodesets[i].nodeCount) {
        status = EX_FATAL;
        ex_opts(EX_VERBOSE);
        errmsg =
            fmt::format("Error: # dist fact ({}) not equal to # nodes ({}) "
                        "in node set {} file id {}",
                        nodesets[i].dfCount, nodesets[i].nodeCount, nodesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        return EX_FATAL;
      }

      // create variable for distribution factors
      status = nc_def_var(exodusFilePtr, VAR_FACT_NS(cur_num_node_sets + 1),
                          nc_flt_code(exodusFilePtr), 1, Data(dims), &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg = fmt::format("Error: node set {} dist factors already exist in file id {}",
                               nodesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg = fmt::format("Error: failed to create node set {} dist factors in file id {}",
                               nodesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return EX_FATAL;
      }
      exi_compress_variable(exodusFilePtr, varid, 2);
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<NodeSet<INT>> &nodesets)
{
  if (nodesets.empty()) {
    return EX_NOERR;
  }

  // Output nodeset ids...
  size_t                    num_nodesets = nodesets.size();
  std::vector<ex_entity_id> nodeset_id(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    nodeset_id[i] = nodesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_NS_IDS, nodeset_id) != NC_NOERR) {
    return EX_FATAL;
  }

  // Now, write the status array
  std::vector<int> status(num_nodesets);
  for (size_t i = 0; i < num_nodesets; i++) {
    status[i] = nodesets[i].nodeCount > 0 ? 1 : 0;
  }

  if (put_array(exodusFilePtr, VAR_NS_STAT, status) != NC_NOERR) {
    return EX_FATAL;
  }

  return EX_NOERR;
}

// ========================================================================
template <typename INT>
int Excn::Internals<INT>::put_metadata(const std::vector<SideSet<INT>> &sidesets)
{
  if (sidesets.empty()) {
    return EX_NOERR;
  }

  std::string errmsg;

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status = 0; // clear error code

  // Get number of side sets defined for this file
  int dimid;
  status = nc_inq_dimid(exodusFilePtr, DIM_NUM_SS, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      errmsg = fmt::format("Error: no side sets defined for file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    else {
      errmsg =
          fmt::format("Error: failed to locate side sets defined in file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    return EX_FATAL;
  }

  size_t num_side_sets = ex_inquire_int(exodusFilePtr, EX_INQ_SIDE_SETS);
  SMART_ASSERT(sidesets.size() == num_side_sets);

  for (size_t i = 0; i < num_side_sets; i++) {

    //  NOTE: exi_inc_file_item is used to find the number of side sets
    // for a specific file and returns that value incremented.
    int cur_num_side_sets =
        (int)exi_inc_file_item(exodusFilePtr, exi_get_counter_list(EX_SIDE_SET));

    if (sidesets[i].sideCount == 0) {
      continue;
    }

    status = nc_def_dim(exodusFilePtr, DIM_NUM_SIDE_SS(cur_num_side_sets + 1),
                        sidesets[i].sideCount, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: side set {} already defined in file id {}", sidesets[i].id,
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to define number of sides for set {} in file id {}",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }

    std::array dims{dimid};
    int        varid = 0;
    status = nc_def_var(exodusFilePtr, VAR_ELEM_SS(cur_num_side_sets + 1), bulk_type, 1, Data(dims),
                        &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: side set {} element list already defined in file id {}",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create side set {} element list in file id {}",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }
    exi_compress_variable(exodusFilePtr, varid, 1);

    // create side list variable for side set
    status = nc_def_var(exodusFilePtr, VAR_SIDE_SS(cur_num_side_sets + 1), bulk_type, 1, Data(dims),
                        &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: side list already exists for side set {} in file id {}",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg = fmt::format("Error: failed to create side list for side set {} in file id {}",
                             sidesets[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }
    exi_compress_variable(exodusFilePtr, varid, 1);

    // Create variable for distribution factors if required
    if (sidesets[i].dfCount > 0) {
      status = nc_def_dim(exodusFilePtr, DIM_NUM_DF_SS(cur_num_side_sets + 1), sidesets[i].dfCount,
                          &dimid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg = fmt::format("Error: side set df count {} already defined in file id {}",
                               sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg = fmt::format("Error: failed to define side set df count for set {} in file id {}",
                               sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return EX_FATAL;
      }

      // create distribution factor list variable for side set
      dims[0] = dimid;
      status  = nc_def_var(exodusFilePtr, VAR_FACT_SS(cur_num_side_sets + 1),
                           nc_flt_code(exodusFilePtr), 1, Data(dims), &varid);
      if (status != NC_NOERR) {
        ex_opts(EX_VERBOSE);
        if (status == NC_ENAMEINUSE) {
          errmsg =
              fmt::format("Error: dist factor list already exists for side set {} in file id {}",
                          sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        else {
          errmsg =
              fmt::format("Error: failed to create dist factor list for side set {} in file id {}",
                          sidesets[i].id, exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
        }
        return EX_FATAL;
      }
      exi_compress_variable(exodusFilePtr, varid, 2);
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<SideSet<INT>> &sidesets)
{
  if (sidesets.empty()) {
    return EX_NOERR;
  }

  // Output sideset ids...
  int                       num_sidesets = (int)sidesets.size();
  std::vector<ex_entity_id> sideset_id(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    sideset_id[i] = sidesets[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_SS_IDS, sideset_id) != NC_NOERR) {
    return EX_FATAL;
  }

  // Now, write the status array
  std::vector<int> status(num_sidesets);
  for (int i = 0; i < num_sidesets; i++) {
    status[i] = sidesets[i].sideCount > 0 ? 1 : 0;
  }

  if (put_array(exodusFilePtr, VAR_SS_STAT, status) != NC_NOERR) {
    return EX_FATAL;
  }

  return EX_NOERR;
}

// ========================================================================
template <typename INT>
int Excn::Internals<INT>::put_metadata(const std::vector<EdgeBlock<INT>> &edgeblocks)
{
  if (edgeblocks.empty()) {
    return EX_NOERR;
  }

  std::string errmsg;

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status = 0; // clear error code

  // Get number of edge blocks defined for this file
  int dimid;
  status = nc_inq_dimid(exodusFilePtr, DIM_NUM_ED_BLK, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      errmsg = fmt::format("Error: no edge blocks defined for file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    else {
      errmsg =
          fmt::format("Error: failed to locate edge blocks defined in file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    return EX_FATAL;
  }

  size_t num_edge_blocks = ex_inquire_int(exodusFilePtr, EX_INQ_EDGE_BLK);
  SMART_ASSERT(edgeblocks.size() == num_edge_blocks);

  for (size_t i = 0; i < num_edge_blocks; i++) {

    //  NOTE: exi_inc_file_item is used to find the number of edge blocks
    // for a specific file and returns that value incremented.
    int cur_num_edge_blocks =
        (int)exi_inc_file_item(exodusFilePtr, exi_get_counter_list(EX_EDGE_BLOCK));

    if (edgeblocks[i].edgeCount == 0) {
      continue;
    }

    int num_edges_in_edgeblock_dim;
    status = nc_def_dim(exodusFilePtr, DIM_NUM_ED_IN_EBLK(cur_num_edge_blocks + 1),
                        edgeblocks[i].edgeCount, &num_edges_in_edgeblock_dim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: edge block {} already defined in file id {}", edgeblocks[i].id,
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg =
            fmt::format("Error: failed to define number of edges for edge block {} in file id {}",
                        edgeblocks[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }

    int num_nodes_per_edge_dim;
    status = nc_def_dim(exodusFilePtr, DIM_NUM_NOD_PER_ED(cur_num_edge_blocks + 1),
                        edgeblocks[i].nodesPerEdge, &num_nodes_per_edge_dim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of nodes/edge for block {} in file id {}",
                      edgeblocks[i].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // element connectivity array
    std::array dims{num_edges_in_edgeblock_dim, num_nodes_per_edge_dim};
    int        connid;
    status = nc_def_var(exodusFilePtr, VAR_EBCONN(cur_num_edge_blocks + 1), bulk_type, 2,
                        Data(dims), &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to create connectivity array for edge block {} in file id {}",
                      edgeblocks[i].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
    exi_compress_variable(exodusFilePtr, connid, 1);

    // store edge type as attribute of connectivity variable
    status = nc_put_att_text(exodusFilePtr, connid, ATT_NAME_ELB,
                             static_cast<int>(std::strlen(edgeblocks[i].elType)) + 1,
                             edgeblocks[i].elType);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store edge type name {} in file id {}",
                           edgeblocks[i].elType, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<EdgeBlock<INT>> &edgeblocks)
{
  if (edgeblocks.empty()) {
    return EX_NOERR;
  }

  // Output edgeblock ids...
  int                       num_edgeblocks = (int)edgeblocks.size();
  std::vector<ex_entity_id> edgeblock_id(num_edgeblocks);
  for (int i = 0; i < num_edgeblocks; i++) {
    edgeblock_id[i] = edgeblocks[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_ID_ED_BLK, edgeblock_id) != NC_NOERR) {
    return EX_FATAL;
  }

  // Now, write the status array
  std::vector<int> status(num_edgeblocks);
  for (int i = 0; i < num_edgeblocks; i++) {
    status[i] = edgeblocks[i].edgeCount > 0 ? 1 : 0;
  }

  if (put_array(exodusFilePtr, VAR_STAT_ED_BLK, status) != NC_NOERR) {
    return EX_FATAL;
  }

  return EX_NOERR;
}

// ========================================================================
template <typename INT>
int Excn::Internals<INT>::put_metadata(const std::vector<FaceBlock<INT>> &faceblocks)
{
  if (faceblocks.empty()) {
    return EX_NOERR;
  }

  std::string errmsg;

  int bulk_type = get_type(exodusFilePtr, EX_BULK_INT64_DB);

  int status = 0; // clear error code

  // Get number of face blocks defined for this file
  int dimid;
  status = nc_inq_dimid(exodusFilePtr, DIM_NUM_FA_BLK, &dimid);
  if (status != NC_NOERR) {
    ex_opts(EX_VERBOSE);
    if (status == NC_EBADDIM) {
      errmsg = fmt::format("Error: no face blocks defined for file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    else {
      errmsg =
          fmt::format("Error: failed to locate face blocks defined in file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
    }
    return EX_FATAL;
  }

  size_t num_face_blocks = ex_inquire_int(exodusFilePtr, EX_INQ_FACE_BLK);
  SMART_ASSERT(faceblocks.size() == num_face_blocks);

  for (size_t i = 0; i < num_face_blocks; i++) {

    //  NOTE: exi_inc_file_item is used to find the number of face blocks
    // for a specific file and returns that value incremented.
    int cur_num_face_blocks =
        (int)exi_inc_file_item(exodusFilePtr, exi_get_counter_list(EX_FACE_BLOCK));

    if (faceblocks[i].faceCount == 0) {
      continue;
    }

    int num_faces_in_faceblock_dim;
    status = nc_def_dim(exodusFilePtr, DIM_NUM_FA_IN_FBLK(cur_num_face_blocks + 1),
                        faceblocks[i].faceCount, &num_faces_in_faceblock_dim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      if (status == NC_ENAMEINUSE) {
        errmsg = fmt::format("Error: face block {} already defined in file id {}", faceblocks[i].id,
                             exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      else {
        errmsg =
            fmt::format("Error: failed to define number of faces for face block {} in file id {}",
                        faceblocks[i].id, exodusFilePtr);
        ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      }
      return EX_FATAL;
    }

    int num_nodes_per_face_dim;
    status = nc_def_dim(exodusFilePtr, DIM_NUM_NOD_PER_FA(cur_num_face_blocks + 1),
                        faceblocks[i].nodesPerFace, &num_nodes_per_face_dim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define number of nodes/face for block {} in file id {}",
                      faceblocks[i].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // element connectivity array
    std::array dims{num_faces_in_faceblock_dim, num_nodes_per_face_dim};

    int connid;
    status = nc_def_var(exodusFilePtr, VAR_FBCONN(cur_num_face_blocks + 1), bulk_type, 2,
                        Data(dims), &connid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to create connectivity array for face block {} in file id {}",
                      faceblocks[i].id, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
    exi_compress_variable(exodusFilePtr, connid, 1);

    // store face type as attribute of connectivity variable
    status = nc_put_att_text(exodusFilePtr, connid, ATT_NAME_ELB,
                             static_cast<int>(std::strlen(faceblocks[i].elType)) + 1,
                             faceblocks[i].elType);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store face type name {} in file id {}",
                           faceblocks[i].elType, exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
  }
  return EX_NOERR;
}

template <typename INT>
int Excn::Internals<INT>::put_non_define_data(const std::vector<FaceBlock<INT>> &faceblocks)
{
  if (faceblocks.empty()) {
    return EX_NOERR;
  }

  // Output faceblock ids...
  int                       num_faceblocks = (int)faceblocks.size();
  std::vector<ex_entity_id> faceblock_id(num_faceblocks);
  for (int i = 0; i < num_faceblocks; i++) {
    faceblock_id[i] = faceblocks[i].id;
  }

  if (put_id_array(exodusFilePtr, VAR_ID_FA_BLK, faceblock_id) != NC_NOERR) {
    return EX_FATAL;
  }

  // Now, write the status array
  std::vector<int> status(num_faceblocks);
  for (int i = 0; i < num_faceblocks; i++) {
    status[i] = faceblocks[i].faceCount > 0 ? 1 : 0;
  }

  if (put_array(exodusFilePtr, VAR_STAT_FA_BLK, status) != NC_NOERR) {
    return EX_FATAL;
  }

  return EX_NOERR;
}

namespace {
  template <typename INT>
  int put_array(int exoid, const char *var_type, const std::vector<INT> &array)
  {
    std::string errmsg;
    int         var_id;
    int         status;

    status = nc_inq_varid(exoid, var_type, &var_id);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to locate {} in file id {}", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    if (sizeof(INT) == sizeof(int64_t)) {
      status = nc_put_var_longlong(exoid, var_id, (long long int *)Data(array));
    }
    else {
      status = nc_put_var_int(exoid, var_id, Data(array));
    }

    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to write {} array in file id {}", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
    return EX_NOERR;
  }

  int put_id_array(int exoid, const char *var_type, const std::vector<ex_entity_id> &ids)
  {
    std::string errmsg;
    int         var_id;

    int status = nc_inq_varid(exoid, var_type, &var_id);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to locate {} in file id {}", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    int id_type = get_type(exoid, EX_IDS_INT64_API);

    if (id_type == NC_INT64) {
      status = nc_put_var_longlong(exoid, var_id, (long long int *)Data(ids));
    }
    else {
      // Have ex_entity_id (long long), need ints...
      std::vector<int> int_ids(ids.size());
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4244)
#endif
      int_ids.assign(ids.begin(), ids.end());
#ifdef _MSC_VER
#pragma warning(pop)
#endif
      status = nc_put_var_int(exoid, var_id, Data(int_ids));
    }

    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to write {} array in file id {}", var_type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
    return EX_NOERR;
  }

  int define_coordinate_vars(int exodusFilePtr, int64_t nodes, int node_dim, int dimension,
                             int dim_dim, int str_dim)
  {
    std::string errmsg;
    int         status;
    int         varid;

    if (nodes > 0) {
      // node coordinate arrays -- separate storage...
      std::array dim{node_dim};
      if (dimension > 0) {
        status = nc_def_var(exodusFilePtr, VAR_COORD_X, nc_flt_code(exodusFilePtr), 1, Data(dim),
                            &varid);
        if (status != NC_NOERR) {
          ex_opts(EX_VERBOSE);
          errmsg = fmt::format("Error: failed to define node x coordinate array in file id {}",
                               exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
          return EX_FATAL;
        }
        exi_compress_variable(exodusFilePtr, varid, 2);
      }

      if (dimension > 1) {
        status = nc_def_var(exodusFilePtr, VAR_COORD_Y, nc_flt_code(exodusFilePtr), 1, Data(dim),
                            &varid);
        if (status != NC_NOERR) {
          ex_opts(EX_VERBOSE);
          errmsg = fmt::format("Error: failed to define node y coordinate array in file id {}",
                               exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
          return EX_FATAL;
        }
        exi_compress_variable(exodusFilePtr, varid, 2);
      }

      if (dimension > 2) {
        status = nc_def_var(exodusFilePtr, VAR_COORD_Z, nc_flt_code(exodusFilePtr), 1, Data(dim),
                            &varid);
        if (status != NC_NOERR) {
          ex_opts(EX_VERBOSE);
          errmsg = fmt::format("Error: failed to define node z coordinate array in file id {}",
                               exodusFilePtr);
          ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
          return EX_FATAL;
        }
        exi_compress_variable(exodusFilePtr, varid, 2);
      }
    }

    // coordinate names array
    std::array dims{dim_dim, str_dim};
    status = nc_def_var(exodusFilePtr, VAR_NAME_COOR, NC_CHAR, 2, Data(dims), &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg =
          fmt::format("Error: failed to define coordinate name array in file id {}", exodusFilePtr);
      ex_err_fn(exodusFilePtr, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
    return EX_NOERR;
  }

  int define_netcdf_vars(int exoid, const char *type, size_t count, const char *dim_num,
                         const char *stat_var, const char *id_var, const char *name_var)
  {
    int         dimid      = 0;
    int         varid      = 0;
    int         namestrdim = 0;
    std::string errmsg;

    int status = nc_inq_dimid(exoid, DIM_STR_NAME, &namestrdim);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to get string name dimension in file id {}", exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    status = nc_def_dim(exoid, dim_num, count, &dimid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define number of {}s in file id {}", type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // id status array:
    std::array dim{dimid};
    status = nc_def_var(exoid, stat_var, NC_INT, 1, Data(dim), &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define side {} status in file id {}", type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // id array:
    int ids_type = get_type(exoid, EX_IDS_INT64_DB);
    status       = nc_def_var(exoid, id_var, ids_type, 1, Data(dim), &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define {} property in file id {}", type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    // store property name as attribute of property array variable
    status = nc_put_att_text(exoid, varid, ATT_PROP_NAME, 3, "ID");
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to store {} property name {} in file id {}", type, "ID",
                           exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }

    std::array dims{dimid, namestrdim};
    status = nc_def_var(exoid, name_var, NC_CHAR, 2, Data(dims), &varid);
    if (status != NC_NOERR) {
      ex_opts(EX_VERBOSE);
      errmsg = fmt::format("Error: failed to define {} name array in file id {}", type, exoid);
      ex_err_fn(exoid, __func__, errmsg.c_str(), status);
      return EX_FATAL;
    }
    return EX_NOERR;
  }
} // namespace
