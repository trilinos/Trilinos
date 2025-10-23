#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/OutputStreams.hpp>
#include "stk_mesh/base/CreateFaces.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_unit_test_utils/timer.hpp"
#include "stk_tools/mesh_tools/EntityDisconnectTool.hpp"


inline constexpr size_t MAX_NEIGHBORS = 40;
class ElementDisconnectPerfFixture : public stk::unit_test_util::MeshFixture
{
 public:
  ElementDisconnectPerfFixture() : stk::unit_test_util::MeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }

  void create_mesh_with_faces(unsigned nx_, unsigned ny_, unsigned nz_)
  {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    stk::io::fill_mesh(get_generated_mesh_description(), get_bulk());
    stk::mesh::create_faces(get_bulk());
  }

  std::string get_generated_mesh_description() const
  {
    return "generated:" + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + "|sideset:xXyYzZ";
  }

  template <typename FacePred>
  auto get_faces(FacePred predicate)
  {
    const auto& bulk = get_bulk();
    auto predBound = [&bulk, predicate](const stk::mesh::Entity& face) { return predicate(bulk, face); };
    auto facesFiltered = stk::mesh::get_entities(bulk, stk::topology::FACE_RANK);
    auto endFiltered = std::stable_partition(facesFiltered.begin(), facesFiltered.end(), predBound);
    facesFiltered.erase(endFiltered, facesFiltered.end());
    stk::util::sort_and_unique(facesFiltered);
    return facesFiltered;
  }

  template <typename Func>
  void run_repeated(const unsigned numRepeats, Func func)
  {
    auto comm = stk::parallel_machine_world();
    stk::unit_test_util::BatchTimer batchTimer(comm);
    batchTimer.initialize_batch_timer();

    for (unsigned iRepeat = 0; iRepeat < numRepeats; ++iRepeat) {
      func(batchTimer);
    }
    batchTimer.print_batch_timing(numRepeats);
  }

  int nx;
  int ny;
  int nz;
};

TEST_F(ElementDisconnectPerfFixture, disconnect_all_faces_batched_destroy_declare) {
  static constexpr auto numRepeats = 3u;
  run_repeated(numRepeats, [this](stk::unit_test_util::BatchTimer& batchTimer) {
    reset_mesh();
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    create_mesh_with_faces(100, 100, 100);
    auto allFaces = get_faces([](const stk::mesh::BulkData& bulk, const stk::mesh::Entity& face) { return true; });

    batchTimer.start_batch_timer();
    stk::experimental::EntityDisconnectTool disconnecter(get_bulk(), allFaces);
    disconnecter.determine_new_nodes();
    disconnecter.modify_mesh();
    batchTimer.stop_batch_timer();
  });
}