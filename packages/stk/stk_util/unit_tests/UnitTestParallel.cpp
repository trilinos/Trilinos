
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

class UnitTestParallel : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestParallel );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION( UnitTestParallel);

void UnitTestParallel::testUnit()
{
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  std::string s;
  std::ostringstream strout;

  std::cout << "all_write_string " << std::flush;
  
//  for (size_t i = 0; i < 250000; ++i) {
  for (size_t i = 0; i < 100; ++i) {
    if (mpi_rank == 0 && i%1000 == 0)
      std::cout << "." << std::flush;
    
    stk::all_write_string(MPI_COMM_WORLD, strout, s);
  }
  
  CPPUNIT_ASSERT(mpi_rank < mpi_size);
}
