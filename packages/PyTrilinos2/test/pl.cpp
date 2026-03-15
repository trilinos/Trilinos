#include <exception>
#include <iostream>
#include <cxxabi.h>
#include <source_location>
#include <Tpetra_Core.hpp>
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>
#include "Teuchos_stacktrace.hpp"
#include "MueLu_MasterList.hpp"

#include "pybind11/embed.h"


void handle_exception(std::exception &e,
                      const std::source_location location = std::source_location::current()) {

  std::string error_message = e.what();
  std::string stacktrace;
  std::string mangled_exception_type = typeid(e).name();

  int status                     = 0;
  char* demangled_exception_type     = 0;
  demangled_exception_type           = abi::__cxa_demangle(mangled_exception_type.c_str(), 0, 0, &status);
  std::string exception_type     = demangled_exception_type;

#ifdef HAVE_TEUCHOS_STACKTRACE
  stacktrace = Teuchos::get_stacktrace();
#endif

  // init python interpreter
  pybind11::scoped_interpreter guard{};

  auto data = pybind11::dict(pybind11::arg("filename") = location.file_name(),
                             pybind11::arg("line") = location.line(),
                             pybind11::arg("column") = location.column(),
                             pybind11::arg("function_name") = location.function_name(),
                             pybind11::arg("error_message") = error_message,
                             pybind11::arg("stacktrace") = stacktrace,
                             pybind11::arg("exception_type") = exception_type);

  pybind11::exec(R"(
      print("filename:", filename)
      print("line:", line)
      print("column:", column)
      print("function_name:", function_name)
      print("exception_type", exception_type)
      print("error_message", error_message)
      print("stacktrace", stacktrace)
                  )", pybind11::globals(), data);

}

int main(int argc, char **argv) {
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    std::string xmlFileName = "test.xml";
    auto comm = Tpetra::getDefaultComm();
    auto paramList = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName,
                                                     paramList.ptr(), *comm);
    const auto &validList = *MueLu::MasterList::List();
    try {
      paramList->validateParameters(validList);
    } catch (Teuchos::Exceptions::InvalidParameterName &e) {
      handle_exception(e);
      // throw;
    }
  }
}
