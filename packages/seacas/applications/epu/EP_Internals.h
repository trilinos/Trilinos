/*
 * Copyright(C) 2010-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
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
 *     * Neither the name of NTESS nor the names of its
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
