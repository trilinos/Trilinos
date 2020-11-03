/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef SEACAS_Internals_h
#define SEACAS_Internals_h

#include <string> // for string
#include <vector> // for vector
namespace Excn {
  class Block;
} // namespace Excn
namespace Excn {
  class CommunicationMetaData;
} // namespace Excn
namespace Excn {
  class Mesh;
} // namespace Excn
namespace Excn {
  template <typename INT> class NodeSet;
} // namespace Excn
namespace Excn {
  template <typename INT> class SideSet;
} // namespace Excn
namespace Excn {
} // namespace Excn
namespace Excn {
} // namespace Excn
namespace Excn {
} // namespace Excn
namespace Excn {
} // namespace Excn
namespace Excn {
} // namespace Excn

/*!
 * This set of classes provides a thin wrapper around the exodusII
 * internals.  It supplants several of the exodusII API calls in
 * order to avoid ncredef calls which totally rewrite the existing
 * database and can be very expensive.  These routines provide all
 * required variable, dimension, and attribute definitions to the
 * underlying netcdf file with only a single ncredef call.
 *
 * To use the application must create an Excn::Internals instance
 * and call the Excn::Internals::write_meta_data() function.  This
 * function requires several classes as arguments including:
 * <ul>
 * <li> Mesh -- defines mesh global metadata
 * <li> Block -- defines metadata for each block
 * <li> NodeSet -- defines metadata for each nodeset
 * <li> SideSet -- defines metadata for each sideset
 * <li> CommunicationMetaData -- global metadata relating to
 * parallel info.
 * </ul>
 *
 * Calling Excn::Internals::write_meta_data(), replaces the
 * following exodusII and nemesis API calls:
 * <ul>
 * <li> ex_put_init(),
 * <li> ex_put_elem_block(),
 * <li> ex_put_node_set_param(),
 * <li> ex_put_side_set_param(),
 * <li> ne_put_init_info(),
 * <li> ne_put_loadbal_param(),
 * <li> ne_put_cmap_params(),
 * </ul>
 */

namespace Excn {

  class Redefine
  {
  public:
    explicit Redefine(int exoid);
    ~Redefine();

  private:
    int exodusFilePtr;
  };

  bool is_path_absolute(const std::string &path);

  template <typename INT> class Internals
  {
  public:
    Internals(int exoid, int maximum_name_length)
        : exodusFilePtr(exoid), nodeMapVarID(), elementMapVarID(), commIndexVar(0),
          elemCommIndexVar(0), maximumNameLength(maximum_name_length)
    {
    }

    int write_meta_data(const Mesh &mesh, const std::vector<Block> &blocks,
                        const std::vector<NodeSet<INT>> &nodesets,
                        const std::vector<SideSet<INT>> &sidesets,
                        const CommunicationMetaData &    comm);

    bool check_meta_data(const Mesh &mesh, const std::vector<Block> &blocks,
                         const std::vector<NodeSet<INT>> &nodesets,
                         const std::vector<SideSet<INT>> &sidesets,
                         const CommunicationMetaData &    comm);

  private:
    int put_metadata(const Mesh &mesh, const CommunicationMetaData &comm);
    int put_metadata(const std::vector<Block> &blocks);
    int put_metadata(const std::vector<NodeSet<INT>> &nodesets);
    int put_metadata(const std::vector<SideSet<INT>> &sidesets);

    int put_non_define_data(const Mesh &mesh, const CommunicationMetaData &comm);
    int put_non_define_data(const std::vector<Block> &blocks);
    int put_non_define_data(const std::vector<NodeSet<INT>> &nodesets);
    int put_non_define_data(const std::vector<SideSet<INT>> &sidesets);

    int exodusFilePtr;
    int nodeMapVarID[3];
    int elementMapVarID[2];
    int commIndexVar;
    int elemCommIndexVar;
    int maximumNameLength;
  };
} // namespace Excn
#endif /* SEACAS_Internals_h */
