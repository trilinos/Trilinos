#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include "UnitTestSpiderMeshSetup.hpp"

namespace {

class SpiderElement : public stk::unit_test_util::simple_fields::MeshFixture
{
protected:
  SpiderElement()
  {
    m_balanceSettings.setShouldFixSpiders(true);

    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::balance::internal::register_internal_fields(get_bulk(), m_balanceSettings);
  }

  void fill_decomp_list_from_current_ownership(stk::mesh::EntityProcVec & decomp) {
    stk::mesh::EntityVector localElems;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), localElems);
    for (const stk::mesh::Entity & elem : localElems) {
      decomp.push_back(std::make_pair(elem, get_bulk().parallel_owner_rank(elem)));
    }
  }

  void fix_spider_elements()
  {
    stk::balance::internal::fill_spider_connectivity_count_fields(get_bulk(), m_balanceSettings);
    stk::mesh::EntityProcVec decomp;
    fill_decomp_list_from_current_ownership(decomp);
    stk::balance::internal::fill_output_subdomain_field(get_bulk(), m_balanceSettings, decomp);
    m_changeList = std::make_unique<stk::balance::DecompositionChangeList>(get_bulk(), decomp);

    stk::balance::internal::fix_spider_elements(m_balanceSettings, get_bulk(), *m_changeList);
  }

  void check_spider_body(const stk::mesh::EntityKeyProcVec & spiderBody)
  {
    check_entity_key_ownership(spiderBody);
  }

  void check_spider_legs(const stk::mesh::EntityKeyProcVec & spiderLegs)
  {
    check_entity_key_ownership(spiderLegs);
  }

private:
  void check_entity_key_ownership(const stk::mesh::EntityKeyProcVec & entityKeyProcs)
  {
    for (const stk::mesh::EntityKeyProc & entityKeyProc : entityKeyProcs) {
      const stk::mesh::Entity entity = get_bulk().get_entity(entityKeyProc.first);
      if (get_bulk().is_valid(entity) && get_bulk().bucket(entity).owned()) {
        const int expectedNewOwner = entityKeyProc.second;
        const int changedNewOwner = m_changeList->get_entity_destination(entity);
        const int actualNewOwner = (changedNewOwner >= 0) ? changedNewOwner : get_bulk().parallel_owner_rank(entity);

        EXPECT_EQ(actualNewOwner, expectedNewOwner)
            << "New owner for " << entityKeyProc.first << " is " << actualNewOwner << " and not " << expectedNewOwner;
      }
    }
  }

  stk::balance::StkBalanceSettings m_balanceSettings;
  std::unique_ptr<stk::balance::DecompositionChangeList> m_changeList;
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
// Final: (Same as initial)
//

TEST_F(SpiderElement, notASpider_EnoughLegsNoVolumeElements_NoMovement)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  make_mesh_non_spider_no_volume_elements(get_bulk());

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 7), 2}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 3), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 4), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 5), 3},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 6), 3}});
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
// Final: (Same as initial)
//

TEST_F(SpiderElement, notASpider_NotEnoughLegs_NoMovement)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  make_mesh_non_spider_not_enough_legs(get_bulk());

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 2}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 3},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 3}});
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

TEST_F(SpiderElement, oneSpider_NoBodyElement)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  make_mesh_one_spider_no_body_element(get_bulk());

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});
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
//

TEST_F(SpiderElement, oneSpider_ParticleBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  make_mesh_one_spider_particle_body(get_bulk());

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});
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
//

TEST_F(SpiderElement, oneSpider_BeamBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  make_mesh_one_spider_beam_body(get_bulk());

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});
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
//
//
//
TEST_F(SpiderElement, compoundSpider_BeamBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  make_mesh_compound_spider_beam_body(get_bulk());

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 38), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 21), 0}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
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
//
//
//
TEST_F(SpiderElement, twoSpiders_ParticleBody)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  make_mesh_two_spiders_particle_body(get_bulk());

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0},

                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 38), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 22), 2}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
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
}

}
