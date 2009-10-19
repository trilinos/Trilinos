#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <mpi.h>

class UnitTestTemplate : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestTemplate );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION( UnitTestTemplate);

void UnitTestTemplate::testUnit()
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  CPPUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  CPPUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
  
  CPPUNIT_ASSERT(mpi_rank < mpi_size);
}
