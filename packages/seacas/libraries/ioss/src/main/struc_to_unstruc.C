// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_FileInfo.h>
#include <Ioss_Hex8.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SubSystem.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_TerminalColor.h>
#include <Ioss_Utils.h>
#include <cgns/Iocgns_StructuredZoneData.h>
#include <cgns/Iocgns_Utils.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <unistd.h>
#include <vector>

#define OUTPUT                                                                                     \
  if (rank == 0)                                                                                   \
  std::cerr

// ========================================================================

namespace {

  int rank = 0;

  void show_step(int istep, double time);

  void transfer_nodal(const Ioss::Region &region, Ioss::Region &output_region);
  void transfer_connectivity(Ioss::Region &region, Ioss::Region &output_region);
  void output_sidesets(Ioss::Region &region, Ioss::Region &output_region);

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region);
  void transfer_elementblocks(Ioss::Region &region, Ioss::Region &output_region);
  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region);
  void create_unstructured(const std::string &inpfile, const std::string &outfile);

  void transfer_sb_fields(const Ioss::Region &region, Ioss::Region &output_region,
                          Ioss::Field::RoleType role);

  template <typename T>
  void transfer_fields(const std::vector<T *> &entities, Ioss::Region &output_region,
                       Ioss::Field::RoleType role);

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role);

  void transfer_sb_field_data(const Ioss::Region &region, Ioss::Region &output_region,
                              Ioss::Field::RoleType role);

  void transfer_coord(std::vector<double> &to, std::vector<double> &from,
                      std::vector<size_t> &node_id_list)
  {
    assert(from.empty() || !node_id_list.empty());
    for (size_t i = 0; i < from.size(); i++) {
      size_t idx = node_id_list[i];
      to[idx]    = from[i];
    }
  }
} // namespace
// ========================================================================

namespace {
  std::string codename;
  std::string version = "4.7";
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  Ioss::Init::Initializer io;
  std::string             in_file = argv[1];
  std::string             out_file;
  if (argc > 2) {
    out_file = argv[2];
  }
  else {
    OUTPUT << "ERROR: Syntax is " << argv[0] << " {structured_input} {unstructured_output}\n";
    return EXIT_FAILURE;
  }

  OUTPUT << "Structured Input:    '" << in_file << "'\n";
  OUTPUT << "Unstructured Output: '" << out_file << "'\n";
  OUTPUT << '\n';

  double begin = Ioss::Utils::timer();
  create_unstructured(in_file, out_file);
  double end = Ioss::Utils::timer();

  OUTPUT << "\n\tElapsed time = " << end - begin << " seconds.\n";

  OUTPUT << "\n" << codename << " execution successful.\n";
#ifdef SEACAS_HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

namespace {
  void create_unstructured(const std::string &inpfile, const std::string &outfile)
  {
    Ioss::PropertyManager properties;
    Ioss::DatabaseIO *    dbi =
        Ioss::IOFactory::create("cgns", inpfile, Ioss::READ_MODEL, MPI_COMM_WORLD, properties);
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    if (region.mesh_type() != Ioss::MeshType::STRUCTURED) {
      int myProcessor = region.get_database()->util().parallel_rank();
      if (myProcessor == 0) {
        std::cerr << "\nERROR: The input mesh is not of type STRUCTURED.\n";
      }
      return;
    }

      //========================================================================
      // OUTPUT ...
      //========================================================================
#if 0
    if (dbi->util().parallel_size() > 1) {
      properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
    }
#endif

    Ioss::DatabaseIO *dbo =
        Ioss::IOFactory::create("exodus", outfile, Ioss::WRITE_RESTART, MPI_COMM_WORLD, properties);
    if (dbo == nullptr || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "region_2");
    // Set the qa information...
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
      OUTPUT << "ERROR: Could not put output region into define model state\n";
      std::exit(EXIT_FAILURE);
    }

    transfer_nodeblock(region, output_region);
    transfer_elementblocks(region, output_region);
    transfer_sidesets(region, output_region);

    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    // Model defined, now fill in the model data...
    output_region.begin_mode(Ioss::STATE_MODEL);

    transfer_connectivity(region, output_region);
    transfer_nodal(region, output_region);
    output_sidesets(region, output_region);
    output_region.end_mode(Ioss::STATE_MODEL);

    if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
      OUTPUT << "\n Number of time steps on database     =" << std::setw(12)
             << region.get_property("state_count").get_int() << "\n\n";

      output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

      // For each 'TRANSIENT' field in the node blocks and element
      // blocks, transfer to the output node and element blocks.
      transfer_fields(&region, &output_region, Ioss::Field::TRANSIENT);
      transfer_sb_fields(region, output_region, Ioss::Field::TRANSIENT);

      //      transfer_fields(region.get_structured_blocks(), output_region,
      //      Ioss::Field::TRANSIENT);

      output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

      output_region.begin_mode(Ioss::STATE_TRANSIENT);
      // Get the timesteps from the input database.  Step through them
      // and transfer fields to output database...

      int step_count = region.get_property("state_count").get_int();

      for (int istep = 1; istep <= step_count; istep++) {
        double time = region.get_state_time(istep);

        int ostep = output_region.add_state(time);
        show_step(istep, time);

        output_region.begin_state(ostep);
        region.begin_state(istep);

        transfer_sb_field_data(region, output_region, Ioss::Field::TRANSIENT);

        region.end_state(istep);
        output_region.end_state(ostep);
      }
      output_region.end_mode(Ioss::STATE_TRANSIENT);
    }
  }

  void transfer_nodal(const Ioss::Region &region, Ioss::Region &output_region)
  {
    size_t num_nodes = region.get_node_blocks()[0]->entity_count();
    auto   nb        = output_region.get_node_blocks()[0];

    if (!output_region.get_database()->needs_shared_node_information()) {
      std::vector<int> ids(num_nodes); // To hold the global node id map.
      auto &           blocks = region.get_structured_blocks();
      for (auto &block : blocks) {
        std::vector<int> cell_id;
        block->get_field_data("cell_node_ids", cell_id);

        for (size_t i = 0; i < cell_id.size(); i++) {
          size_t idx = block->m_blockLocalNodeIndex[i];
          assert(idx < num_nodes);
          if (ids[idx] == 0) {
            ids[idx] = cell_id[i];
          }
        }
      }
      assert(nb != nullptr);
      nb->put_field_data("ids", ids);
    }

    std::vector<double> coordinate_x(num_nodes);
    std::vector<double> coordinate_y(num_nodes);
    std::vector<double> coordinate_z(num_nodes);

    auto &blocks = region.get_structured_blocks();
    for (auto &block : blocks) {
      std::vector<double> coord_tmp;
      block->get_field_data("mesh_model_coordinates_x", coord_tmp);
      transfer_coord(coordinate_x, coord_tmp, block->m_blockLocalNodeIndex);

      block->get_field_data("mesh_model_coordinates_y", coord_tmp);
      transfer_coord(coordinate_y, coord_tmp, block->m_blockLocalNodeIndex);

      block->get_field_data("mesh_model_coordinates_z", coord_tmp);
      transfer_coord(coordinate_z, coord_tmp, block->m_blockLocalNodeIndex);
    }

    nb->put_field_data("mesh_model_coordinates_x", coordinate_x);
    nb->put_field_data("mesh_model_coordinates_y", coordinate_y);
    nb->put_field_data("mesh_model_coordinates_z", coordinate_z);
  }

  void transfer_connectivity(Ioss::Region &region, Ioss::Region &output_region)
  {
    auto &blocks = region.get_structured_blocks();
    for (auto &block : blocks) {
      // We have a structured block of size ni x nj x nk.
      // Need to convert that to element connectivity
      // Node numbers are zero-based offset into this structured block
      // After generated, then map zero-based block-local into one-based global.
      // Since we are outputting "connectivity_raw", there is no difference
      // in parallel or serial.

      size_t ni = block->get_property("ni").get_int();
      size_t nj = block->get_property("nj").get_int();
      size_t nk = block->get_property("nk").get_int();

      size_t xp1yp1 = (ni + 1) * (nj + 1);

      // Find matching element block in output region...
      const auto &name   = block->name();
      auto *      output = output_region.get_element_block(name);
      assert(output != nullptr);

      {
        std::vector<int> connect;
        connect.reserve(ni * nj * nk * 8);
        for (size_t k = 0; k < nk; k++) {
          for (size_t j = 0, m = 0; j < nj; j++) {
            for (size_t i = 0; i < ni; i++, m++) {
              size_t base = (k * xp1yp1) + m + j;

              connect.push_back(base);
              connect.push_back(base + 1);
              connect.push_back(base + ni + 2);
              connect.push_back(base + ni + 1);

              connect.push_back(xp1yp1 + base);
              connect.push_back(xp1yp1 + base + 1);
              connect.push_back(xp1yp1 + base + ni + 2);
              connect.push_back(xp1yp1 + base + ni + 1);
            }
          }
        }
        // 'connect' contains 0-based block-local node ids at this point
        // Now, map them to processor-global values...
        // NOTE: "processor-global" is 1..num_node_on_processor
        if (!connect.empty()) {
          const auto &gnil = block->m_blockLocalNodeIndex;
          assert(!gnil.empty());
          for (int &i : connect) {
            i = gnil[i] + 1;
          }
        }

        output->put_field_data("connectivity_raw", connect);
      }

      {
        std::vector<int> ids;
        block->get_field_data("cell_ids", ids);
        output->put_field_data("ids", ids);
      }
    }
  }

  void output_sidesets(Ioss::Region &region, Ioss::Region &output_region)
  {
    // Maps the 'parent_face'+3 returned from which_parent_face()
    // to the local 1-based face of the hex elements in that block.
    static int face_map[] = {5, 1, 4, 0, 2, 3, 6};

    const auto &ssets = region.get_sidesets();
    for (auto ss : ssets) {
      // Get corresponding sidset on output region...
      auto ofs = output_region.get_sideset(ss->name());
      assert(ofs != nullptr);

      const auto &fbs      = ss->get_side_blocks();
      size_t      fb_index = 0;
      for (auto fb : fbs) {
        // Get corresponding sideblock on output sideset 'ofs'
        // Assumes sideblocks are ordered the same on input and output.
        auto ofb = ofs->get_block(fb_index++);
        assert(ofb != nullptr);

        // Get parent structured block for this side block...
        auto parent = fb->parent_block();
        assert(parent != nullptr);
        assert(parent->type() == Ioss::STRUCTUREDBLOCK);
        auto sb_parent = dynamic_cast<const Ioss::StructuredBlock *>(parent);
        Ioss::Utils::check_dynamic_cast(sb_parent);

        // Find this sideblock on the parent block...
        auto &bc_name = fb->name();
        for (auto &bc : sb_parent->m_boundaryConditions) {
          auto bc_compose = bc.m_bcName + "/" + sb_parent->name();
          if (bc_name == bc_compose) {
            std::vector<int> elem_side;
            if (bc.get_face_count() > 0) {
              Ioss::IJK_t range_beg      = bc.m_rangeBeg;
              Ioss::IJK_t cell_range_end = bc.m_rangeEnd;

              // The range_beg/end are current points and not cells.
              // Need to convert cell_range_end to cells which is typically just point-1
              // except if the ordinal is the plane of this surface. In this case,
              // range_beg[ord] == cell_range_end[ord].  If this is the case,
              // and we are at the top end of the range, then the beg and end must
              // both be reduced by 1.  If at the bottom end of the range, then both
              // are equal to 1.
              //
              for (int i = 0; i < 3; i++) {
                if (cell_range_end[i] == range_beg[i]) {
                  if (cell_range_end[i] != 1) {
                    cell_range_end[i]--;
                    range_beg[i]--;
                  }
                }
                else {
                  cell_range_end[i]--;
                }
              }

#if IOSS_DEBUG_OUTPUT
              std::cerr << bc << "\n";
#endif
              auto parent_face = face_map[bc.which_parent_face() + 3];
              elem_side.reserve(bc.get_face_count() * 2);
              for (auto k = range_beg[2]; k <= cell_range_end[2]; k++) {
                for (auto j = range_beg[1]; j <= cell_range_end[1]; j++) {
                  for (auto i = range_beg[0]; i <= cell_range_end[0]; i++) {
                    auto cell_id = sb_parent->get_global_cell_id(i, j, k);
                    assert(cell_id > 0);
                    elem_side.push_back(cell_id);
                    elem_side.push_back(parent_face);
                  }
                }
              }
            }
            ofb->put_field_data("element_side", elem_side);
            break;
          }
        }
      }
    }
  }

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region)
  {
    const auto &nbs = region.get_node_blocks();
    assert(nbs.size() == 1);
    size_t degree    = nbs[0]->get_property("component_degree").get_int();
    size_t num_nodes = nbs[0]->entity_count();
    auto nb = new Ioss::NodeBlock(output_region.get_database(), nbs[0]->name(), num_nodes, degree);
    output_region.add(nb);

    if (output_region.get_database()->needs_shared_node_information()) {
      std::vector<int> ids(num_nodes); // To hold the global node id map.
      auto &           blocks = region.get_structured_blocks();
      for (auto &block : blocks) {
        std::vector<int> cell_id;
        block->get_field_data("cell_node_ids", cell_id);

        for (size_t i = 0; i < cell_id.size(); i++) {
          size_t idx = block->m_blockLocalNodeIndex[i];
          assert(idx < num_nodes);
          if (ids[idx] == 0) {
            ids[idx] = cell_id[i];
          }
        }
      }
      assert(nb != nullptr);
      nb->put_field_data("ids", ids);

      // Each structured block on the incoming mesh has a list of the nodes it shares with
      // other blocks.  Use this to construct the "node owning
      // processor" information.  Assume that if a node is shared with
      // a lower-numbered processor, then that processor owns the
      // node...

      auto shared_nodes = Iocgns::Utils::resolve_processor_shared_nodes(
          region, region.get_database()->util().parallel_rank());

      int              myProcessor = output_region.get_database()->util().parallel_rank();
      std::vector<int> owning_processor(num_nodes, myProcessor);
      for (auto &block : blocks) {
        int zone = block->get_property("zone").get_int();
        for (const auto &shared : shared_nodes[zone]) {
          size_t idx = block->m_blockLocalNodeIndex[shared.first];
          if (owning_processor[idx] > (int)shared.second) {
            owning_processor[idx] = shared.second;
          }
        }
      }
      nb->put_field_data("owning_processor", owning_processor);
    }

    std::cout << "P[" << rank << "] Number of coordinates per node =" << std::setw(12) << degree
              << "\n";
    std::cout << "P[" << rank << "] Number of nodes                =" << std::setw(12) << num_nodes
              << "\n";
  }

  void transfer_elementblocks(Ioss::Region &region, Ioss::Region &output_region)
  {
    const auto &blocks         = region.get_structured_blocks();
    size_t      total_entities = 0;
    for (auto iblock : blocks) {
      const std::string &name  = iblock->name();
      std::string        type  = Ioss::Hex8::name;
      size_t             count = iblock->get_property("cell_count").get_int();
      auto block = new Ioss::ElementBlock(output_region.get_database(), name, type, count);
      output_region.add(block);
#if IOSS_DEBUG_OUTPUT
      std::cerr << "P[" << rank << "] Created Element Block '" << name << "' with " << count
                << " elements.\n";
#endif
      total_entities += count;
    }
    std::cout << "P[" << rank << "] Number of Element Blocks       =" << std::setw(12)
              << blocks.size() << ", Number of elements (cells) =" << std::setw(12)
              << total_entities << "\n";
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region)
  {
    size_t      total_sides = 0;
    const auto &ssets       = region.get_sidesets();
    for (auto ss : ssets) {
      const std::string &name = ss->name();

      int                      ss_sides = 0;
      auto                     surf     = new Ioss::SideSet(output_region.get_database(), name);
      Ioss::SideBlockContainer fbs      = ss->get_side_blocks();
      for (auto fb : fbs) {
        const std::string &fbname   = fb->name();
        std::string        fbtype   = fb->get_property("topology_type").get_string();
        std::string        partype  = fb->get_property("parent_topology_type").get_string();
        size_t             num_side = fb->entity_count();
        total_sides += num_side;

        auto block =
            new Ioss::SideBlock(output_region.get_database(), fbname, fbtype, partype, num_side);
        surf->add(block);
        block->property_add(Ioss::Property("set_offset", ss_sides));

        ss_sides += num_side;
      }
      output_region.add(surf);
    }
    std::cout << "P[" << rank << "] Number of SideSets             =" << std::setw(12)
              << ssets.size() << ", Number of cell faces       =" << std::setw(12) << total_sides
              << "\n";
  }

  void show_step(int istep, double time)
  {
    OUTPUT.setf(std::ios::scientific);
    OUTPUT.setf(std::ios::showpoint);
    OUTPUT << "     Time step " << std::setw(5) << istep << " at time " << std::setprecision(5)
           << time << '\n';
  }

  void transfer_sb_fields(const Ioss::Region &region, Ioss::Region &output_region,
                          Ioss::Field::RoleType role)
  {
    auto   nb        = output_region.get_node_blocks()[0];
    size_t num_nodes = region.get_node_blocks()[0]->entity_count();
    auto & blocks    = region.get_structured_blocks();
    for (auto &block : blocks) {
      Ioss::NameList fields;
      block->field_describe(role, &fields);

      const auto &name   = block->name();
      auto *      eblock = output_region.get_element_block(name);
      assert(eblock != nullptr);

      for (const auto &field_name : fields) {
        Ioss::Field field      = block->get_field(field_name);
        bool        cell_field = Iocgns::Utils::is_cell_field(field);
        if (cell_field) {
          eblock->field_add(field);
        }
        else {
          if (!nb->field_exists(field_name)) {
            field.reset_count(num_nodes);
            nb->field_add(field);
          }
        }
      }
    }
  }

  void transfer_sb_field_data(const Ioss::Region &region, Ioss::Region &output_region,
                              Ioss::Field::RoleType role)
  {
    {
      auto                nb        = output_region.get_node_blocks()[0];
      size_t              num_nodes = region.get_node_blocks()[0]->entity_count();
      std::vector<double> node_data(num_nodes);
      std::vector<double> data;

      // Handle nodal fields first...
      Ioss::NameList fields;
      nb->field_describe(role, &fields);

      for (const auto &field_name : fields) {
        assert(nb->field_exists(field_name));
        Ioss::Field               field      = nb->get_field(field_name);
        const Ioss::VariableType *var_type   = field.raw_storage();
        size_t                    comp_count = var_type->component_count();

        auto &blocks = region.get_structured_blocks();
        for (auto &block : blocks) {
          auto node_id_list = block->m_blockLocalNodeIndex;
          assert(!node_id_list.empty());

          for (size_t i = 0; i < node_id_list.size(); i++) {
            size_t node = node_id_list[i];
            for (size_t j = 0; j < comp_count; j++) {
              node_data[comp_count * node + j] = data[comp_count * i + j];
            }
          }
        }
        nb->put_field_data(field_name, node_data);
      }
    }
    // Now handle the cell-center fields on the structured blocks...
    auto &blocks = region.get_structured_blocks();
    for (auto &block : blocks) {
      Ioss::NameList fields;
      block->field_describe(role, &fields);

      for (auto field_name : fields) {
        Ioss::Field field = block->get_field(field_name);
        if (Iocgns::Utils::is_cell_field(field)) {
          std::vector<double> data;
          block->get_field_data(field_name, data);

          const auto &name   = block->name();
          auto *      eblock = output_region.get_element_block(name);
          assert(eblock != nullptr);
          eblock->put_field_data(field_name, data);
        }
      }
    }
  }

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role)
  {
    // Check for transient fields...
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    // Iterate through results fields and transfer to output
    // database...
    for (const auto &field_name : fields) {
      Ioss::Field field = ige->get_field(field_name);
      if (field_name != "ids" && !oge->field_exists(field_name)) {
        oge->field_add(field);
      }
    }
  }

} // namespace
