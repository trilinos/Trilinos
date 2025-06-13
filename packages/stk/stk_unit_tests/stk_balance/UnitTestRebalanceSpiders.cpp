// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
//
#include "MeshFixtureRebalance.hpp"
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_balance/rebalance.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include "UnitTestSpiderMeshSetup.hpp"
#include <vector>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {
using stk::unit_test_util::build_mesh;

using MeshSetupFunction = void (*)(stk::mesh::BulkData & bulk);

class RebalanceSpiders : public MeshFixtureRebalance
{
public:
  RebalanceSpiders()
    : MeshFixtureRebalance()
  {
    m_balanceSettings.setShouldFixSpiders(true);
    m_balanceSettings.set_is_rebalancing(true);
    m_balanceSettings.set_output_filename(get_output_file_name());
    m_balanceSettings.set_num_input_processors(stk::parallel_machine_size(get_comm()));
  }

  void initialize_mesh(MeshSetupFunction meshSetupFunction)
  {
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, get_comm());
    meshSetupFunction(*bulk);
    stk::io::write_mesh(get_input_file_name(), *bulk);

    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    get_meta().set_coordinate_field_name(m_balanceSettings.getCoordinateFieldName());
    stk::balance::internal::register_internal_fields_and_parts(get_bulk(), m_balanceSettings);

    stk::io::fill_mesh_preexisting(m_ioBroker, get_input_file_name(), get_bulk());
  }

  void rebalance_mesh(int numFinalProcs, const std::string & decompMethod = "rcb")
  {
    m_balanceSettings.set_num_output_processors(numFinalProcs);
    m_balanceSettings.setDecompMethod(decompMethod);

    stk::set_outputP0(&stk::outputNull());
    stk::balance::rebalance(m_ioBroker, m_balanceSettings);
    stk::reset_default_output_streams();
  }

  void check_spider_body_ownership(const stk::mesh::EntityKeyProcVec & spiderBody)
  {
    check_entity_key_ownership(spiderBody);
  }

  void check_spider_legs_ownership(const stk::mesh::EntityKeyProcVec & spiderLegs)
  {
    check_entity_key_ownership(spiderLegs);
  }

private:
  void check_entity_key_ownership(const stk::mesh::EntityKeyProcVec & entityKeyProcs)
  {
    // Make sure all procs have written their files before p0 tries to read them
    MPI_Barrier(get_comm());

    const int numOutputProcs = m_balanceSettings.get_num_output_processors();
    if (get_parallel_rank() == 0) {
      for (int proc = 0; proc < numOutputProcs; ++proc) {
        std::shared_ptr<stk::mesh::BulkData> tempBulk = build_mesh(MPI_COMM_SELF);
        std::string subdomainFileName = get_subdomain_filename(numOutputProcs, proc);
        stk::io::fill_mesh(subdomainFileName, *tempBulk);

        for (const auto & entityKeyProc : entityKeyProcs) {
          if (entityKeyProc.second == proc) {
            const stk::mesh::Entity entity = tempBulk->get_entity(entityKeyProc.first);
            EXPECT_EQ(tempBulk->is_valid(entity), true) << "Entity " << entityKeyProc.first
                                                       << " not found on proc " << proc;
          }
        }
      }
    }
  }

};

//
// Initial:                      7                           Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                   (2:2)//   | |   \  \(6:3)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(5:3)
//                     //(1:2) | |     \   \                                                                    .
//                    //      /  |(4:2) |   \                                                                   .
//                   //       |  |      |    \                                                                  .
//                  2/        |  4      |     6
//                  /    (3:2)|          \                                                                      .
//                 /          |           \                                                                     .
//               1            3            5
//
//
// Final:                        7                           Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                   (2:0)//   | |   \  \(6:3)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(5:3)
//                     //(1:0) | |     \   \                                                                    .
//                    //      /  |(4:1) |   \                                                                   .
//                   //       |  |      |    \                                                                  .
//                  2/        |  4      |     6
//                  /    (3:2)|          \                                                                      .
//                 /          |           \                                                                     .
//               1            3            5
//
TEST_F(RebalanceSpiders, notASpider_EnoughLegsNoVolumeElements_PartitionedNormally)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  initialize_mesh(make_mesh_non_spider_no_volume_elements);

  rebalance_mesh(4);

  check_spider_body_ownership({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 7), 0}});

  // Beams not part of a spider, so partitioned with RCB
  check_spider_legs_ownership({{stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 3), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 4), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 5), 3},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 6), 3}});

  clean_up_temporary_files();
  clean_up_temporary_input_files(); 
}


//
// Initial:                      37                          Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            / || |\                                                                           .
//                           /  || | \                                                                          .
//                          /   || |  \                                                                         .
//                         /   / |  \  \                                                                        .
//                   (9:2)/    | |   \  \(13:3)
//                       /     | |    |  \                                                                      .
//                      /      | |    |(12:3)
//                     /       | |     \   \                                                                    .
//                    /       /  |(11:2)|   \                                                                   .
//                   /        |  |      |    \                                                                  .
//     4------------8---------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           /|   (10:2)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          / |         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
//
// Final                         37                          Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            / || |\                                                                           .
//                           /  || | \                                                                          .
//                          /   || |  \                                                                         .
//                         /   / |  \  \                                                                        .
//                   (9:0)/    | |   \  \(13:1)
//                       /     | |    |  \                                                                      .
//                      /      | |    |(12:1)
//                     /       | |     \   \                                                                    .
//                    /       /  |(11:1)|   \                                                                   .
//                   /        |  |      |    \                                                                  .
//     4------------8---------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           /|   (10:1)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          / |         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:2)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
TEST_F(RebalanceSpiders, notASpider_NotEnoughLegs_PartitionedNormally)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  initialize_mesh(make_mesh_non_spider_not_enough_legs);

  rebalance_mesh(4);

  check_spider_body_ownership({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0}});

  // Beams not part of a spider, so partitioned with RCB along with the other elements
  check_spider_legs_ownership({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1}});

  clean_up_temporary_files();
  clean_up_temporary_input_files(); 
}


//
// Initial:                      37                          Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                  (10:2)//   | |   \  \(14:3)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(13:3)
//                     //(9:2) | |     \   \                                                                    .
//                    //      /  |(12:2)|   \                                                                   .
//                   //       |  |      |    \                                                                  .
//     4------------8/--------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           //   (11:2)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          //|         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
//
// Final:                        37                          Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                  (10:0)//   | |   \  \(14:1)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(13:1)
//                     //(9:0) | |     \   \                                                                    .
//                    //      /  |(12:0)|   \                                                                   .
//                   //       |  |      |    \                                                                  .
//     4------------8/--------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           //   (11:0)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          //|         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
//
TEST_F(RebalanceSpiders, oneSpider_NoBodyElement)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  initialize_mesh(make_mesh_one_spider_no_body_element);

  rebalance_mesh(4);

  check_spider_body_ownership({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0}});

  check_spider_legs_ownership({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});

  clean_up_temporary_files();
  clean_up_temporary_input_files(); 
}


//
// Initial:                      37 (15:3)                   Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                  (10:2)//   | |   \  \(14:3)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(13:3)
//                     //(9:2) | |     \   \                                                                    .
//                    //      /  |(12:2)|   \                                                                   .
//                   //       |  |      |    \                                                                  .
//     4------------8/--------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           //   (11:2)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          //|         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
//
// Final:                        37 (15:0)                   Elem: (id:proc)
//                              /|\                          Node: id
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                  (10:0)//   | |   \  \(14:1)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(13:1)
//                     //(9:0) | |     \   \                                                                    .
//                    //      /  |(12:0)|   \                                                                   .
//                   //       |  |      |    \                                                                  .
//     4------------8/--------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           //   (11:0)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          //|         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
TEST_F(RebalanceSpiders, oneSpider_ParticleBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  initialize_mesh(make_mesh_one_spider_particle_body);

  rebalance_mesh(4);

  check_spider_body_ownership({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0}});

  check_spider_legs_ownership({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});

  clean_up_temporary_files();
  clean_up_temporary_input_files(); 
}


// Initial:                           (15:3)                       Elem: (id:proc)
//                               37------------38                  Node: id
//                              /|\                                                                             .
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                  (10:2)//   | |   \  \(14:3)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(13:3)
//                     //(9:2) | |     \   \                                                                    .
//                    //      /  |(12:2)|   \                                                                   .
//                   //       |  |      |    \                                                                  .
//     4------------8/--------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           //   (11:2)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          //|         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
//
// Final:                             (15:0)                       Elem: (id:proc)
//                               37------------38                  Node: id
//                              /|\                                                                             .
//                             /||\\                                                                            .
//                            //|| |\                                                                           .
//                           // || | \                                                                          .
//                          //  || |  \                                                                         .
//                         //  / |  \  \                                                                        .
//                  (10:0)//   | |   \  \(14:1)
//                       //    | |    |  \                                                                      .
//                      //     | |    |(13:1)
//                     //(9:0) | |     \   \                                                                    .
//                    //      /  |(12:0)|   \                                                                   .
//                   //       |  |      |    \                                                                  .
//     4------------8/--------|--12-----|-----16-----------20-----------24-----------28-----------32-----------36
//    /|           //   (11:0)| /|       \   /|           /|           /|           /|           /|           /|
//   / |          //|         |/ |        \ / |          / |          / |          / |          / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
TEST_F(RebalanceSpiders, oneSpider_BeamBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  initialize_mesh(make_mesh_one_spider_beam_body);

  rebalance_mesh(4);

  check_spider_body_ownership({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0}});

  check_spider_legs_ownership({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});

  clean_up_temporary_files();
  clean_up_temporary_input_files(); 
}


//
// Initial:                      37--------------------------------------------------38             Elem: (id:proc)
//                              /|\                  (21:3)                         /|\             Node: id
//                             /||\\                                               /||\\                        .
//                            //|| |\                                             //|| |\                       .
//                           // || | \                                           // || | \                      .
//                          //  || |  \                                         //  || |  \                     .
//                         //  / |  \  \                                       //  / |  \  \                    .
//                  (10:2)//   | |   \  \(14:3)                         (16:0)//   | |   \  \(20:1)
//                       //    | |    |  \                                   //    | |    |  \                  .
//                      //     | |    |(13:3)                               //     | |    |(19:1)
//                     //(9:2) | |     \   \                               //(15:0)| |     \   \                .
//                    //      /  |(12:2)|   \                             //      /  |(18:0)|   \               .
//                   //       |  |      |    \                           //       |  |      |    \              .
//     4------------8/--------|--12-----|-----16-----------20-----------24--------|--28-----|-----32-----------36
//    /|           //   (11:2)| /|       \   /|           /|           //   (17:0)| /|       \   /|           /|
//   / |          //|         |/ |        \ / |          / |          //|         |/ |        \ / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
//
// Final:                        37--------------------------------------------------38             Elem: (id:proc)
//                              /|\                  (21:0)                         /|\             Node: id
//                             /||\\                                               /||\\                        .
//                            //|| |\                                             //|| |\                       .
//                           // || | \                                           // || | \                      .
//                          //  || |  \                                         //  || |  \                     .
//                         //  / |  \  \                                       //  / |  \  \                    .
//                  (10:0)//   | |   \  \(14:1)                         (16:2)//   | |   \  \(20:3)
//                       //    | |    |  \                                   //    | |    |  \                  .
//                      //     | |    |(13:1)                               //     | |    |(19:3)
//                     //(9:0) | |     \   \                               //(15:2)| |     \   \                .
//                    //      /  |(12:0)|   \                             //      /  |(18:2)|   \               .
//                   //       |  |      |    \                           //       |  |      |    \              .
//     4------------8/--------|--12-----|-----16-----------20-----------24--------|--28-----|-----32-----------36
//    /|           //   (11:0)| /|       \   /|           /|           //   (17:2)| /|       \   /|           /|
//   / |          //|         |/ |        \ / |          / |          //|         |/ |        \ / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
TEST_F(RebalanceSpiders, compoundSpider_BeamBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  initialize_mesh(make_mesh_compound_spider_beam_body);

  rebalance_mesh(4);

  check_spider_body_ownership({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                               {stk::mesh::EntityKey(stk::topology::NODE_RANK, 38), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 21), 0}});

  check_spider_legs_ownership({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1},

                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 16), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 17), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 18), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 19), 3},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 20), 3}});

  clean_up_temporary_files();
  clean_up_temporary_input_files(); 
}


//
// Initial:                      37 (15:3)                                           38 (22:1)      Elem: (id:proc)
//                              /|\                                                 /|\             Node: id
//                             /||\\                                               /||\\                        .
//                            //|| |\                                             //|| |\                       .
//                           // || | \                                           // || | \                      .
//                          //  || |  \                                         //  || |  \                     .
//                         //  / |  \  \                                       //  / |  \  \                    .
//                  (10:2)//   | |   \  \(14:3)                         (17:0)//   | |   \  \(21:1)
//                       //    | |    |  \                                   //    | |    |  \                  .
//                      //     | |    |(13:3)                               //     | |    |(20:1)
//                     //(9:2) | |     \   \                               //(16:0)| |     \   \                .
//                    //      /  |(12:2)|   \                             //      /  |(19:0)|   \               .
//                   //       |  |      |    \                           //       |  |      |    \              .
//     4------------8/--------|--12-----|-----16-----------20-----------24--------|--28-----|-----32-----------36
//    /|           //   (11:2)| /|       \   /|           /|           //   (18:0)| /|       \   /|           /|
//   / |          //|         |/ |        \ / |          / |          //|         |/ |        \ / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
//
// Final:                        37 (15:0)                                           38 (22:2)      Elem: (id:proc)
//                              /|\                                                 /|\             Node: id
//                             /||\\                                               /||\\                        .
//                            //|| |\                                             //|| |\                       .
//                           // || | \                                           // || | \                      .
//                          //  || |  \                                         //  || |  \                     .
//                         //  / |  \  \                                       //  / |  \  \                    .
//                  (10:0)//   | |   \  \(14:1)                         (17:2)//   | |   \  \(21:3)
//                       //    | |    |  \                                   //    | |    |  \                  .
//                      //     | |    |(13:1)                               //     | |    |(20:3)
//                     //(9:0) | |     \   \                               //(16:2)| |     \   \                .
//                    //      /  |(12:0)|   \                             //      /  |(19:2)|   \               .
//                   //       |  |      |    \                           //       |  |      |    \              .
//     4------------8/--------|--12-----|-----16-----------20-----------24--------|--28-----|-----32-----------36
//    /|           //   (11:0)| /|       \   /|           /|           //   (18:2)| /|       \   /|           /|
//   / |          //|         |/ |        \ / |          / |          //|         |/ |        \ / |          / |
//  1------------5------------9------------13-----------17-----------21-----------25-----------29-----------33 |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  |  (1:0)  |  |  (2:0)  |  |  (3:1)  |  |  (4:1)  |  |  (5:2)  |  |  (6:2)  |  |  (7:3)  |  |  (8:3)  |  |
//  |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |         |  |
//  |  3---------|--7---------|--11--------|--15--------|--19--------|--23--------|--27--------|--31--------|--35
//  | /          | /          | /          | /          | /          | /          | /          | /          | /
//  |/           |/           |/           |/           |/           |/           |/           |/           |/
//  2------------6------------10-----------14-----------18-----------22-----------26-----------30-----------34
//
TEST_F(RebalanceSpiders, twoSpiders_ParticleBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  initialize_mesh(make_mesh_two_spiders_particle_body);

  rebalance_mesh(4);

  check_spider_body_ownership({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0},

                               {stk::mesh::EntityKey(stk::topology::NODE_RANK, 38), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 22), 2}});

  check_spider_legs_ownership({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1},

                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 16), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 17), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 18), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 19), 2},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 20), 3},
                               {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 21), 3}});

  clean_up_temporary_files();
  clean_up_temporary_input_files(); 
}

}
