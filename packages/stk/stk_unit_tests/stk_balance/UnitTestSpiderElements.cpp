#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

namespace {

class SpiderElement : public stk::unit_test_util::MeshFixture
{
protected:
  SpiderElement()
  {
    m_balanceSettings.setShouldFixSpiders(true);
    const int initValue = 0;

    stk::mesh::Field<int> & beamConnectivityField =
        get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::NODE_RANK,
                                                        m_balanceSettings.getSpiderBeamConnectivityCountFieldName());
    stk::mesh::put_field_on_mesh(beamConnectivityField, get_meta().universal_part(), &initValue);

    stk::mesh::Field<int> & volumeConnectivityField =
        get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK,
                                                        m_balanceSettings.getSpiderVolumeConnectivityCountFieldName());
    stk::mesh::put_field_on_mesh(volumeConnectivityField, get_meta().universal_part(), &initValue);

    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }

  void make_mesh_non_spider_no_volume_elements()
  {
    std::string meshDesc = "2,1,BEAM_2,1,7\n"
                           "2,2,BEAM_2,2,7\n"
                           "2,3,BEAM_2,3,7\n"
                           "2,4,BEAM_2,4,7\n"
                           "3,5,BEAM_2,5,7\n"
                           "3,6,BEAM_2,6,7\n";

    std::vector<double> coordinates {
      0,0,1, 0,0,0, 1,0,1, 1,0,0, 2,0,1, 2,0,0,
      1,1,0.5
    };

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
  }

  void make_mesh_non_spider_not_enough_legs()
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                           "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                           "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                           "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                           "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                           "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                           "2,9,BEAM_2,8,37\n"
                           "2,10,BEAM_2,9,37\n"
                           "2,11,BEAM_2,12,37\n"
                           "3,12,BEAM_2,13,37\n"
                           "3,13,BEAM_2,16,37\n";

    std::vector<double> coordinates {
      0,1,1, 0,0,1, 0,0,0, 0,1,0,
      1,1,1, 1,0,1, 1,0,0, 1,1,0,
      2,1,1, 2,0,1, 2,0,0, 2,1,0,
      3,1,1, 3,0,1, 3,0,0, 3,1,0,
      4,1,1, 4,0,1, 4,0,0, 4,1,0,
      5,1,1, 5,0,1, 5,0,0, 5,1,0,
      6,1,1, 6,0,1, 6,0,0, 6,1,0,
      7,1,1, 7,0,1, 7,0,0, 7,1,0,
      8,1,1, 8,0,1, 8,0,0, 8,1,0,
      2,2,0.5
    };

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
  }

  void make_mesh_one_spider_no_body_element()
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                           "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                           "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                           "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                           "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                           "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                           "2,9,BEAM_2,5,37\n"
                           "2,10,BEAM_2,8,37\n"
                           "2,11,BEAM_2,9,37\n"
                           "2,12,BEAM_2,12,37\n"
                           "3,13,BEAM_2,13,37\n"
                           "3,14,BEAM_2,16,37\n";

    std::vector<double> coordinates {
      0,1,1, 0,0,1, 0,0,0, 0,1,0,
      1,1,1, 1,0,1, 1,0,0, 1,1,0,
      2,1,1, 2,0,1, 2,0,0, 2,1,0,
      3,1,1, 3,0,1, 3,0,0, 3,1,0,
      4,1,1, 4,0,1, 4,0,0, 4,1,0,
      5,1,1, 5,0,1, 5,0,0, 5,1,0,
      6,1,1, 6,0,1, 6,0,0, 6,1,0,
      7,1,1, 7,0,1, 7,0,0, 7,1,0,
      8,1,1, 8,0,1, 8,0,0, 8,1,0,
      2,2,0.5
    };

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
  }

  void make_mesh_one_spider_particle_body()
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                           "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                           "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                           "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                           "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                           "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                           "2,9,BEAM_2,5,37\n"
                           "2,10,BEAM_2,8,37\n"
                           "2,11,BEAM_2,9,37\n"
                           "2,12,BEAM_2,12,37\n"
                           "3,13,BEAM_2,13,37\n"
                           "3,14,BEAM_2,16,37\n"
                           "3,15,PARTICLE,37\n";

    std::vector<double> coordinates {
      0,1,1, 0,0,1, 0,0,0, 0,1,0,
      1,1,1, 1,0,1, 1,0,0, 1,1,0,
      2,1,1, 2,0,1, 2,0,0, 2,1,0,
      3,1,1, 3,0,1, 3,0,0, 3,1,0,
      4,1,1, 4,0,1, 4,0,0, 4,1,0,
      5,1,1, 5,0,1, 5,0,0, 5,1,0,
      6,1,1, 6,0,1, 6,0,0, 6,1,0,
      7,1,1, 7,0,1, 7,0,0, 7,1,0,
      8,1,1, 8,0,1, 8,0,0, 8,1,0,
      2,2,0.5
    };

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
  }

  void make_mesh_one_spider_beam_body()
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                           "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                           "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                           "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                           "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                           "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                           "2,9,BEAM_2,5,37\n"
                           "2,10,BEAM_2,8,37\n"
                           "2,11,BEAM_2,9,37\n"
                           "2,12,BEAM_2,12,37\n"
                           "3,13,BEAM_2,13,37\n"
                           "3,14,BEAM_2,16,37\n"
                           "3,15,BEAM_2,37,38\n";

    std::vector<double> coordinates {
      0,1,1, 0,0,1, 0,0,0, 0,1,0,
      1,1,1, 1,0,1, 1,0,0, 1,1,0,
      2,1,1, 2,0,1, 2,0,0, 2,1,0,
      3,1,1, 3,0,1, 3,0,0, 3,1,0,
      4,1,1, 4,0,1, 4,0,0, 4,1,0,
      5,1,1, 5,0,1, 5,0,0, 5,1,0,
      6,1,1, 6,0,1, 6,0,0, 6,1,0,
      7,1,1, 7,0,1, 7,0,0, 7,1,0,
      8,1,1, 8,0,1, 8,0,0, 8,1,0,
      2,2,0.5,
      2,3,0.5
    };

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
  }

  void make_mesh_compound_spider_beam_body()
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                           "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                           "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                           "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                           "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                           "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                           "2,9,BEAM_2,5,37\n"
                           "2,10,BEAM_2,8,37\n"
                           "2,11,BEAM_2,9,37\n"
                           "2,12,BEAM_2,12,37\n"
                           "3,13,BEAM_2,13,37\n"
                           "3,14,BEAM_2,16,37\n"
                           "0,15,BEAM_2,21,38\n"
                           "0,16,BEAM_2,24,38\n"
                           "0,17,BEAM_2,25,38\n"
                           "0,18,BEAM_2,28,38\n"
                           "1,19,BEAM_2,29,38\n"
                           "1,20,BEAM_2,32,38\n"
                           "3,21,BEAM_2,37,38\n";

    std::vector<double> coordinates {
      0,1,1, 0,0,1, 0,0,0, 0,1,0,
      1,1,1, 1,0,1, 1,0,0, 1,1,0,
      2,1,1, 2,0,1, 2,0,0, 2,1,0,
      3,1,1, 3,0,1, 3,0,0, 3,1,0,
      4,1,1, 4,0,1, 4,0,0, 4,1,0,
      5,1,1, 5,0,1, 5,0,0, 5,1,0,
      6,1,1, 6,0,1, 6,0,0, 6,1,0,
      7,1,1, 7,0,1, 7,0,0, 7,1,0,
      8,1,1, 8,0,1, 8,0,0, 8,1,0,
      2,2,0.5,
      6,2,0.5
    };

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
  }

  void make_mesh_two_spiders_particle_body()
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                           "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                           "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                           "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                           "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                           "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                           "2,9,BEAM_2,5,37\n"
                           "2,10,BEAM_2,8,37\n"
                           "2,11,BEAM_2,9,37\n"
                           "2,12,BEAM_2,12,37\n"
                           "3,13,BEAM_2,13,37\n"
                           "3,14,BEAM_2,16,37\n"
                           "3,15,PARTICLE,37\n"
                           "0,16,BEAM_2,21,38\n"
                           "0,17,BEAM_2,24,38\n"
                           "0,18,BEAM_2,25,38\n"
                           "0,19,BEAM_2,28,38\n"
                           "1,20,BEAM_2,29,38\n"
                           "1,21,BEAM_2,32,38\n"
                           "1,22,PARTICLE,38\n";

    std::vector<double> coordinates {
      0,1,1, 0,0,1, 0,0,0, 0,1,0,
      1,1,1, 1,0,1, 1,0,0, 1,1,0,
      2,1,1, 2,0,1, 2,0,0, 2,1,0,
      3,1,1, 3,0,1, 3,0,0, 3,1,0,
      4,1,1, 4,0,1, 4,0,0, 4,1,0,
      5,1,1, 5,0,1, 5,0,0, 5,1,0,
      6,1,1, 6,0,1, 6,0,0, 6,1,0,
      7,1,1, 7,0,1, 7,0,0, 7,1,0,
      8,1,1, 8,0,1, 8,0,0, 8,1,0,
      2,2,0.5,
      6,2,0.5
    };

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);
  }

  void fix_spider_elements()
  {
    stk::balance::internal::fill_spider_connectivity_count_fields(get_bulk(), m_balanceSettings);
    stk::balance::internal::fix_spider_elements(m_balanceSettings, get_bulk());
  }

  void check_spider_body(const stk::mesh::EntityKeyProcVec & spiderBody)
  {
    check_entity_key_ownership(spiderBody);
  }

  void check_spider_legs(const stk::mesh::EntityKeyProcVec & spiderLegs)
  {
    check_entity_key_ownership(spiderLegs);
  }

  void check_spider_feet(const stk::mesh::EntityKeyProcVec & spiderFeet)
  {
    check_entity_key_ownership(spiderFeet);
  }

private:
  void check_entity_key_ownership(const stk::mesh::EntityKeyProcVec & entityKeyProcs)
  {
    for (const stk::mesh::EntityKeyProc & entityKeyProc : entityKeyProcs) {
      if (entityKeyProc.second == get_parallel_rank()) {
        const stk::mesh::Entity entity = get_bulk().get_entity(entityKeyProc.first);
        EXPECT_TRUE(get_bulk().is_valid(entity) && get_bulk().bucket(entity).owned())
                    << entityKeyProc.first << " not owned on proc " << entityKeyProc.second;
      }
    }
  }

  stk::balance::StkBalanceSettings m_balanceSettings;
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

  make_mesh_non_spider_no_volume_elements();

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 7), 2}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 2), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 3), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 4), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 5), 3},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 6), 3}});

  check_spider_feet({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 2), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 3), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 4), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 5), 3},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 6), 3}});
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

  make_mesh_non_spider_not_enough_legs();

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 2}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 2},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 3},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 3}});

  check_spider_feet({{stk::mesh::EntityKey(stk::topology::NODE_RANK,  5), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  8), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 16), 1}});
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

  make_mesh_one_spider_no_body_element();

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});

  check_spider_feet({{stk::mesh::EntityKey(stk::topology::NODE_RANK,  5), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  8), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 16), 1}});
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

  make_mesh_one_spider_particle_body();

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});

  check_spider_feet({{stk::mesh::EntityKey(stk::topology::NODE_RANK,  5), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  8), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 16), 1}});
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

  make_mesh_one_spider_beam_body();

  fix_spider_elements();

  check_spider_body({{stk::mesh::EntityKey(stk::topology::NODE_RANK, 37), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 15), 0}});

  check_spider_legs({{stk::mesh::EntityKey(stk::topology::ELEM_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 10), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 11), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::ELEM_RANK, 14), 1}});

  check_spider_feet({{stk::mesh::EntityKey(stk::topology::NODE_RANK,  5), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  8), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 16), 1}});
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

  make_mesh_compound_spider_beam_body();

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

  check_spider_feet({{stk::mesh::EntityKey(stk::topology::NODE_RANK,  5), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  8), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 16), 1},

                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 21), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 24), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 25), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 28), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 29), 3},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 32), 3}});
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

  make_mesh_two_spiders_particle_body();

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

  check_spider_feet({{stk::mesh::EntityKey(stk::topology::NODE_RANK,  5), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  8), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK,  9), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 12), 0},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 13), 1},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 16), 1},

                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 21), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 24), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 25), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 28), 2},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 29), 3},
                     {stk::mesh::EntityKey(stk::topology::NODE_RANK, 32), 3}});
}

}
