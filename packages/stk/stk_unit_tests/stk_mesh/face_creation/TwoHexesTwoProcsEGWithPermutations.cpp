#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "FaceCreatorFixture.hpp"

namespace
{

class FaceCreatorVaryingPerumutations : public FaceCreatorFixture
{
protected:

  void run_test(const std::vector<std::vector<int>> input_permutations)
  {
    set_permutations_for_both_procs(input_permutations);
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
  }

  void set_permutations_for_both_procs(const std::vector<std::vector<int>> input_permutations)
  {
    m_permutations_per_proc = input_permutations;
  }

  virtual unsigned get_permuted_index(unsigned i)
  {
    return m_permutations_per_proc[get_bulk().parallel_rank()][i];
  }

private:
  std::vector<std::vector<int> > m_permutations_per_proc;
};


TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation0)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {0, 1, 2, 3},
      {0, 3, 2, 1}
    };
    run_test(permutation);
  }
}

TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation1)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {1, 2, 3, 0},
      {3, 2, 1, 0}
    };
    run_test(permutation);
  }
}

TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation2)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {2, 3, 0, 1},
      {2, 1, 0, 3}
    };
    run_test(permutation);
  }
}

TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation3)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {3, 0, 1, 2},
      {1, 0, 3, 2}
    };
    run_test(permutation);
  }
}

TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation4)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {0, 3, 2, 1},
      {0, 1, 2, 3}
    };
    run_test(permutation);
  }
}

TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation5)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {3, 2, 1, 0},
      {1, 2, 3, 0}
    };
    run_test(permutation);
  }
}

TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation6)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {2, 1, 0, 3},
      {2, 3, 0, 1}
    };
    run_test(permutation);
  }
}

TEST_F(FaceCreatorVaryingPerumutations, twoHexesTwoProcsCreateTwoFacesPermutation7)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    std::vector<std::vector<int>> permutation = {
      {1, 0, 3, 2},
      {3, 0, 1, 2}
    };
    run_test(permutation);
  }
}

}
