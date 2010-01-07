#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_HybridPlatform.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FileInputSource.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListReader.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>

std::string fnMatrix("bcsstk17.rsa");

template <class Node>
class runTest {
  public:
  static void run(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) {
    using std::cout; 
    using std::endl;
    cout << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << endl;
    //
    // Get the data from the HB file and build the Map,Matrix
    //
    Teuchos::RCP< Tpetra::CrsMatrix<float,int,int,Node> > A;
    Tpetra::Utils::readHBMatrix(fnMatrix,comm,node,A);
  }
};

int main(int argc, char **argv) {
  using std::string;
  using std::cout;
  using std::endl;
  using Teuchos::FileInputSource;
  using Teuchos::XMLObject;
  using Teuchos::ParameterList;
  using Teuchos::XMLParameterListReader;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;

  Teuchos::GlobalMPISession mpisess(&argc,&argv,&cout);
  RCP<const Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));

  //
  // Get test parameters from command-line processor
  //  
  Teuchos::CommandLineProcessor cmdp(false,true);
  string fnMachine("machine.xml");
  cmdp.setOption("matrix-file",&fnMatrix,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("machine-file",&fnMachine,"Filename for XML machine description file.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // read machine file and initialize platform
  // 
  FileInputSource fileSrc(fnMachine);
  XMLObject machXML = fileSrc.getObject();
  XMLParameterListReader pl2xml;
  ParameterList machPL = pl2xml.toParameterList(machXML);
  Tpetra::HybridPlatform platform(comm,machPL);
  platform.runUserCode<runTest>();

  cout << "End Result: TEST PASSED" << endl;
  return 0;
}
