#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// I'm trying to figure out what is the best way to use parameter list in MueLu.
// ParameterList provides a lot of capabilities but it's not always obvious how to use this class properly as
// there is a lot of important decisions to be made.
//
// We want to:
// - provide a list of valid parameters for documentation
// - delegate parameter validation directly to the class that are using the parameters.
// - output the used/unused parameters + added default parameters to the list
//
// See also:
// - http://trilinos.sandia.gov/packages/docs/dev/packages/teuchos/doc/html/classTeuchos_1_1ParameterList.html
//   (this example is weird because the same list is used for both input parameters and validation)
// - the Teuchos::ParameterListAcceptor/ParameterListAcceptorDefaultBase class
// - Other package usage of parameter list (Ifpack2, ...)

#include "MueLu_ParamListTest.hpp"

class MyFactory : public ParameterListAcceptorImpl {

public:

  MyFactory() { }

  virtual ~MyFactory() { }

  RCP<const ParameterList> getValidParameters(const ParameterList& pL = ParameterList()) const {
    //std::cout << "MyFactory::getValidParameters()" << std::endl;
    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validator_type;

    Teuchos::ParameterList paramList(pL); // make a copy to avoid setting [use]/[unused] flags here. Even if the input list is const, these flags are modified!
    RCP<ParameterList> validParamList = rcp(new ParameterList()); // output list

    validParamList->set("ParamA", 0.1, "Documentation string for param A");
    validParamList->set("ParamB", 0.2, "Documentation string for param B");
    validParamList->set("ParamC", 0.3, "Documentation string for param C");
    validParamList->set("ParamD", 0.4, "Documentation string for param D");

    return validParamList;
  }

  // Main algorithm
  //
  // - We do not want to set default parameters inside of the algorithm.
  //  => Otherwise, it's difficult to track the defaults. Cf. ML.
  //  => use ParameterList::get() without the default value input parameter.
  void Build() {

    if (paramList_.get<double>("ParamA") == 0.5) { } // change "[used]"/["unused"] flag
    if (paramList_.get<double>("ParamC") == 0.5) { }

    // statsParamList_.set(...);

  }

  // - Do we want to store output stats on the same parameter list?
  //   => better to distinguish as stats parameters are not valid input.
  // - RCP? => a view seems enough.
  Teuchos::ParameterList statsParamList_;

};

int main(int argc, char* argv[]) {
  //
  // Documentation
  //
  std::cout << "\n#\n# Documentation\n#\n" << std::endl;
  MyFactory dummy; dummy.getDocumentation(std::cout);

  //

  std::cout << "#\n# main()\n#\n" << std::endl;

  //
  // User parameter list
  //

  // Note: users can also directly modify the validParamList list but this is not really useful.

  ParameterList paramList;
  paramList.set("ParamA", 0.001);
  paramList.set("ParamB", 0.002);

  std::cout << "# Input parameter list:" << std::endl;
  std::cout << paramList << std::endl << std::endl;

  //
  // Validation of the user parameter list
  //

  MyFactory f;
  f.setParameterList(paramList);

  if (0) {// if users want to keep their list untouched:
    ParameterList tmp(paramList);
    f.setParameterList(tmp);
  }

  std::cout << "# Parameter list after validation:" << std::endl;
  std::cout << paramList << std::endl << std::endl;

  //
  // Algorithm
  //

  f.Build();

  std::cout << "# Parameter list after algorithm (flag used/unused):" << std::endl;
  if (0) // do not work with my design: flags used/unused are not set for the initial parameter list
    std::cout << paramList << std::endl << std::endl;

  std::cout << f.getParameterList() << std::endl << std::endl;

  // See also ~MyFactory()

  return 0;
}
