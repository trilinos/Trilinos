#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_FileInputSource.hpp"

using namespace Teuchos;

int main(int argc, char* argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

 try {
 
  FileInputSource InputFile("xml-data");

  XMLObject Object = InputFile.getObject();
  
  cout << Object.toString();
 
 } catch(std::exception& e) {
    cerr << "caught exception " << e.what() << endl;
 }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
