#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "MueLu_ParameterListAcceptor.hpp"

namespace MueLu {

// Another version of ParameterListAcceptorImpl with conditional parameters: parameters "T" and "K" are available only according to the value of the parameter "Solver".
//
// Note: There is something in Teuchos/Optika to deal with such cases. I haven't figured out yet how to use it.
// See for instance: http://trilinos.sandia.gov/packages/docs/r10.2/packages/optika/browser/doc/html/test_2dependencyandexec_2main_8cpp-source.html

using Teuchos::ParameterList;
using Teuchos::RCP;

class ParameterListAcceptorAdvImpl : public ParameterListAcceptorImpl {
 public:
  ParameterListAcceptorAdvImpl() {}

  virtual ~ParameterListAcceptorAdvImpl() {}

  // This functions add *all* the extra parameters recursively using GetValidParameterListSimple
  RCP<const ParameterList> GetValidParameterList() const {
    RCP<const ParameterList> validParamList = GetValidParameterListSimple();

    int numParams;
    do {
      numParams      = validParamList->numParams();
      validParamList = GetValidParameterListSimple();
    } while (validParamList->numParams() != numParams);

    return validParamList;
  }

  // GetValidParameterListSimple only add one extra level of default parameters. Ex: if "Solver" is not set in the input list "pL",
  // extra parameters "T" or "K" are not added to the validParamList.
  virtual RCP<const ParameterList> GetValidParameterListSimple() const = 0;

  void GetDocumentation(std::ostream& os) const {
    GetAdvancedDocumentation(os);
  }

 private:
  virtual void GetAdvancedDocumentation(std::ostream& os) const = 0;
};

class MyFactory : public ParameterListAcceptorAdvImpl {
 public:
  MyFactory() {}

  virtual ~MyFactory() {}

  RCP<const ParameterList> GetValidParameterListSimple() const {
    RCP<ParameterList> validParamList = Teuchos::rcp(new ParameterList());  // output list

    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validator_type;
    validParamList->set("Solver", "ILUT", "The type of solver to use.", Teuchos::rcp(new validator_type(Teuchos::tuple<std::string>("ILUT", "ILUK"), "Solver")));

    AddILUTParameters(*validParamList);
    AddILUKParameters(*validParamList);

    return validParamList;
  }

  // Main algorithm
  void Build() {
    const ParameterList& pL = GetParameterList();
    std::string type        = pL.get<std::string>("Solver");
    if (type == "ILUT") {
      pL.get<double>("T");

    } else if (type == "ILUK") {
      pL.get<int>("K");
    }
  }

 private:
  // Separates functions to be used by both GetValidParameterListSimple and GetDocumentation.
  static void AddILUTParameters(ParameterList& paramList) {
    paramList.set("T", 0.1, "ILUT threshold");
  }

  static void AddILUKParameters(ParameterList& paramList) {
    paramList.set("K", 1, "ILUK level of fill");
  }

  void GetAdvancedDocumentation(std::ostream& os) const {
    os << "## Parameters:" << std::endl;
    printParameterListOptions(os, *GetValidParameterListSimple());

    os << "# ILUT specific parameters:" << std::endl;
    {
      ParameterList p;
      AddILUKParameters(p);
      printParameterListOptions(os, p);
    }

    os << "# ILUK specific parameters:" << std::endl;
    {
      ParameterList p;
      AddILUTParameters(p);
      printParameterListOptions(os, p);
    }

    os << "## Fully described default method:" << std::endl;
    GetValidParameterList()->print(os, 2, true, false);
    os << std::endl;
  }
};

}  // namespace MueLu

int main(int argc, char* argv[]) {
  using MueLu::MyFactory;
  using Teuchos::ParameterList;

  bool success = false;
  try {
    //
    // Documentation
    //
    std::cout << "\n#\n# Documentation\n#\n"
              << std::endl;
    MyFactory dummy;
    dummy.GetDocumentation(std::cout);

    //

    std::cout << "#\n# main()\n#\n"
              << std::endl;

    //
    // User parameter list
    //

    ParameterList paramList;
    // paramList.set("Solver", "ILUK");

    std::cout << "# Input parameter list:" << std::endl;
    std::cout << paramList << std::endl
              << std::endl;

    //
    // Validation of the user parameter list
    //

    MyFactory f;
    f.SetParameterList(paramList);

    std::cout << "# Parameter list after validation:" << std::endl;
    std::cout << paramList << std::endl
              << std::endl;

    //
    // Algorithm
    //

    f.Build();

    std::cout << "# Parameter list after algorithm (flag used/unused):" << std::endl;

    std::cout << f.GetParameterList() << std::endl
              << std::endl;

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
