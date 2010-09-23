#include <iostream>
#include <string>

#include "src_file.hpp"

#ifdef MYAPP_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif

#ifdef MYAPP_EPETRA
#include "Epetra_SerialDenseVector.h"
#endif



int main(int argc, char *argv[]) {

  int status = 0;

  // Initialize MPI and timer
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  // Use of flag set in CMakeLists.txt by how Trilinos was configured
#ifdef MYAPP_MPI
    Teuchos::MpiComm<int> comm =
      Teuchos::MpiComm<int>(Teuchos::opaqueWrapper((MPI_Comm)MPI_COMM_WORLD));
#else
    Teuchos::SerialComm<int> comm = Teuchos::SerialComm<int>();
#endif


  try {
    std::string infile("input.xml");

    // Function from another file
    status = buildDemo::src_file(infile, comm);

   // Flag set in CMakeLists.txt that detects if Epetra was enabled in Trilinos
#ifdef MYAPP_EPETRA
    const int len=10;
    Epetra_SerialDenseVector vec(len);
    if (vec.Length() != len) status += 1000;
    std::cout << "\nEpetra called for vec of length " << len << std::endl;
#endif
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    status = 10;
  }
  catch (std::string& s) {
    std::cout << s << std::endl;
    status = 20;
  }
  catch (char *s) {
    std::cout << s << std::endl;
    status = 30;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    status = 40;
  }

  // Status=0 signals to ctest that the test passed.
  return status;
}
