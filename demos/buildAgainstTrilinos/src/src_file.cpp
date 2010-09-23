#include "src_file.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

int buildDemo::src_file(std::string& infile, Teuchos::Comm<int>& comm) {

  Teuchos::ParameterList p;
  Teuchos::updateParametersFromXmlFileAndBroadcast(infile, &p, comm);

  std::cout << "\nProcessor " << comm.getRank()
            << "  has param list : \n" << p << std::endl;

  return 0;
}
