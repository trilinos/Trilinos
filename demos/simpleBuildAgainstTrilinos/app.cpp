#include "Teuchos_Comm.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include <string>

namespace buildDemo {

int getLocalVectorSize(std::string& infile, const Teuchos::Comm<int> &comm) {
  using Teuchos::inOutArg;
  Teuchos::ParameterList p;
  Teuchos::updateParametersFromXmlFileAndBroadcast(infile, inOutArg(p), comm);
  return Teuchos::get<int>(p.sublist("Application"), "my local vector size");
}

} // namespace buildDemo

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard mpiSession(&argc, &argv);
  auto comm = Tpetra::getDefaultComm();
  auto out = Teuchos::VerboseObjectBase::getDefaultOStream();

  std::string inputFileName("input.xml");
  const int localSize = buildDemo::getLocalVectorSize(inputFileName, *comm);
  const int globalSize = localSize * comm->getSize();

  typedef Tpetra::Vector<>::global_ordinal_type GO_t;
  auto map = Tpetra::createUniformContigMap<int, GO_t>(globalSize, comm);
  auto vec = Tpetra::createVector<double>(map);
  vec->putScalar(1.0);

  *out << "vec.norm1() = " << vec->norm1() << "\n";

  return 0;
}

// NOTE: To understand the Teuchos and Tpetra code, see the Teuchos and Tpetra
// Doxygen documentation online.
