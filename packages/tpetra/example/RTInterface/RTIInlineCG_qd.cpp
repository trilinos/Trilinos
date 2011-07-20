#include "RTIInlineCG.hpp"
#include <qd/qd_real.h>

int main(int argc, char *argv[])
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::ParameterList;

  // 
  // Get the communicator
  //
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  auto comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myImageID = comm->getRank();

  //
  // Get example parameters from command-line processor
  //  
  bool verbose = (myImageID==0);
  std::string xmlfile;
  std::string machineFile;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("param-file", &xmlfile,"XML file for solver parameters");
  cmdp.setOption("machine-file",&machineFile,"Filename for XML machine description file.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // read machine file and initialize platform
  // 
  RCP<Teuchos::ParameterList> machinePL = Teuchos::parameterList();
  std::string defaultMachine(
    " <ParameterList>                                                               "
    "   <ParameterList name='%1=0'>                                                 "
    "     <Parameter name='NodeType'     type='string' value='Kokkos::SerialNode'/> "
    "   </ParameterList>                                                            "
    " </ParameterList>                                                              "
  );
  Teuchos::updateParametersFromXmlString(defaultMachine,machinePL.getRawPtr());
  if (machineFile != "") Teuchos::updateParametersFromXmlFile(machineFile,machinePL.getRawPtr());

  // 
  // create the platform object
  // 
  Tpetra::HybridPlatform platform(comm,*machinePL);

  //
  // instantiate a driver on the scalar stack
  //
  CGDriver<qd_real> driver;
  // hand output stream to driver
  if (verbose) driver.out = Teuchos::getFancyOStream(Teuchos::rcp(&std::cout,false));
  else         driver.out = Teuchos::getFancyOStream(Teuchos::rcp(new Teuchos::oblackholestream()));

  //
  // get the solver parameters
  // 
  RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
  // default solver stack parameters
  std::string xmlString(
    " <ParameterList>                                                  \n"
    "   <Parameter name='tolerance' value='5e-15' type='double'/>      \n"
    "   <Parameter name='verbose' value='1' type='int'/>               \n"
    "   <Parameter name='kappa' value='100' type='double'/>            \n"
    "   <Parameter name='size'  value='100' type='int'/>               \n"
    "   <Parameter name='numIters' value='100' type='int'/>            \n"
    " </ParameterList>                                                 \n"
  );
  Teuchos::updateParametersFromXmlString(xmlString,params.getRawPtr());
  if (xmlfile != "") Teuchos::updateParametersFromXmlFile(xmlfile,params.getRawPtr());
  // hand solver parameters to driver
  driver.params = params;

  // 
  // run the driver
  // 
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  platform.runUserCode(driver);
  fpu_fix_end(&old_cw);

  //
  // Print result
  if (driver.testPassed) {
    *driver.out << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}
