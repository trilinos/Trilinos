
#include <use_cases/UseCase_1.hpp>
#include <use_cases/UseCase_2.hpp>
#include <use_cases/UseCase_3.hpp>
#include <use_cases/UseCase_4.hpp>
#include <use_cases/UseCase_ElementDeath_1.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

void printStatus(bool status)
{
  if (status) {
    std::cout << "passed" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
  }
}

int main ( int argc, char * argv[] )
{
  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  stk::ParallelMachine parallel_machine = use_case_environment.m_comm;

  bool status = true;
  {
    std::cout << "Use Case 1 ... ";
    stk::mesh::use_cases::UseCase_1_Mesh mesh(parallel_machine);
    stk::mesh::use_cases::populate(mesh,1,3);
    bool local_status = stk::mesh::use_cases::verifyMesh(mesh,1,3);
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case 2 ... ";
    stk::mesh::use_cases::UseCase_2_Mesh mesh(parallel_machine);
    stk::mesh::use_cases::populate(mesh,1,3);
    bool local_status = stk::mesh::use_cases::verifyMesh(mesh,1,3);
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case 3 ... ";
    stk::mesh::use_cases::UseCase_3_Mesh mesh(parallel_machine);
    stk::mesh::use_cases::populate(mesh);
    bool local_status = stk::mesh::use_cases::verifyMesh(mesh);
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case 4 ... ";
    stk::mesh::use_cases::UseCase_4_Mesh mesh(parallel_machine);
    stk::mesh::use_cases::populate(mesh);
    stk::mesh::use_cases::runAlgorithms(mesh);
    bool local_status = stk::mesh::use_cases::verifyMesh(mesh);
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case Element Death 1 ... ";
    bool local_status = element_death_use_case(parallel_machine);
    printStatus(local_status);
    status = status && local_status;
  }

  bool result = -1;
  if (status) {
    result = 0;
  }
  std::cout << "Overall Use Case testing results: ";
  printStatus(status);
  std::cout << std::endl;

  return result;
}
