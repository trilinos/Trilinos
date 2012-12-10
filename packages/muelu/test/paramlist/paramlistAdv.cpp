#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "MueLu_ParameterListAcceptor.hpp"

namespace MueLu{

  // Another version of ParameterListAcceptorImpl with conditional parameters: parameters "T" and "K" are available only according to the value of the parameter "Solver".
  //
  // Note: There is something in Teuchos/Optika to deal with such cases. I haven't figured out yet how to use it.
  // See for instance: http://trilinos.sandia.gov/packages/docs/r10.2/packages/optika/browser/doc/html/test_2dependencyandexec_2main_8cpp-source.html

  class ParameterListAcceptorAdvImpl: public ParameterListAcceptorImpl {

  public:

    ParameterListAcceptorAdvImpl() { }

    virtual ~ParameterListAcceptorAdvImpl() { }

    // This functions add *all* the extra parameters recursively using getValidParametersSimple
    RCP<const ParameterList> getValidParameters(const ParameterList& pL = ParameterList()) const {
      RCP<const ParameterList> validParamList = getValidParametersSimple(pL);

      int numParams;
      do {
        numParams = validParamList->numParams();
        validParamList = getValidParametersSimple(*validParamList);
      } while (validParamList->numParams() != numParams);

      return validParamList;
    }

    // getValidParametersSimple only add one extra level of default parameters. Ex: if "Solver" is not set in the input list "pL",
    // extra parameters "T" or "K" are not added to the validParamList.
    virtual RCP<const ParameterList> getValidParametersSimple(const ParameterList& pL = ParameterList()) const = 0;

    void setParameterList(ParameterList & paramList) {
      // Validate and add defaults parameters.
      paramList.validateParametersAndSetDefaults(*getValidParameters(paramList));

      // Do we need to validate range of int/double parameters here? Can it be done automatically as for valid std::string using getValidParameters()?
      paramList_ = paramList; // copy
    }

    const Teuchos::ParameterList & getParameterList() const {
      return paramList_;
    }

    void getDocumentation(std::ostream &os) const {
      getAdvancedDocumentation(os);
    }

  private:

    virtual void getAdvancedDocumentation(std::ostream &os) const = 0;

  };



  class MyFactory : public ParameterListAcceptorAdvImpl {

  public:

    MyFactory() { }

    virtual ~MyFactory() { }

    RCP<const ParameterList> getValidParametersSimple(const ParameterList& pL = ParameterList()) const {
      //std::cout << "MyFactory::getValidParameters()" << std::endl;
      typedef Teuchos::StringToIntegralParameterEntryValidator<int> validator_type;

      Teuchos::ParameterList paramList(pL); // make a copy to avoid setting [use]/[unused] flags here. Even if the input list is const, these flags are modified!
      RCP<ParameterList> validParamList = rcp(new ParameterList()); // output list

      validParamList->set("Solver", "ILUT", "The type of solver to use.", rcp(new validator_type(Teuchos::tuple<std::string>("ILUT", "ILUK"), "Solver")));

      if (paramList.isParameter("Solver")) {
        // conditional parameters

        std::string type = paramList.get<std::string>("Solver");
        // std::cout << "getValidParameters::pL         =>" << std::cout << paramList << std::endl; // previous get() should not set [used] flag.
        // std::cout << "getValidParameters::paramList: =>" << std::cout << pL << std::endl;

        if (type == "ILUT") {
          addILUTParameters(*validParamList);
        } else if (type == "ILUK") {
          addILUKParameters(*validParamList);
        } else {
          // not a valid parameter value. What to do? Ignore. We are not validating at this point.
        }
      }

      return validParamList;
    }

    // Main algorithm
    void Build() {

      std::string type = paramList_.get<std::string>("Solver");
      if (type == "ILUT") {
        paramList_.get<double>("T");

      } else if (type == "ILUK") {
        paramList_.get<int>("K");
      }

    }

  private:
    // Separates functions to be used by both getValidParametersSimple and getDocumentation.
    static void addILUTParameters(ParameterList& paramList) {
      paramList.set("T", 0.1, "ILUT threshold");
    }

    static void addILUKParameters(ParameterList& paramList) {
      paramList.set("K", 1, "ILUK level of fill");
    }

    void getAdvancedDocumentation(std::ostream &os) const {

      os << "## Parameters:" << std::endl;
      printParameterListOptions(os, *getValidParametersSimple());

      os << "# ILUT specific parameters:" << std::endl;
      { ParameterList p; addILUKParameters(p); printParameterListOptions(os, p); }


      os << "# ILUK specific parameters:" << std::endl;
      { ParameterList p; addILUTParameters(p); printParameterListOptions(os, p); }

      os << "## Fully described default method:" << std::endl;
      getValidParameters()->print(os, 2, true, false);
      os << std::endl;
    }

  };

}

int main(int argc, char* argv[]) {
  using Teuchos::ParameterList;
  using MueLu::MyFactory;

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

  ParameterList paramList;
  // paramList.set("Solver", "ILUK");

  std::cout << "# Input parameter list:" << std::endl;
  std::cout << paramList << std::endl << std::endl;

  //
  // Validation of the user parameter list
  //

  MyFactory f;
  f.setParameterList(paramList);

  std::cout << "# Parameter list after validation:" << std::endl;
  std::cout << paramList << std::endl << std::endl;

  //
  // Algorithm
  //

  f.Build();

  std::cout << "# Parameter list after algorithm (flag used/unused):" << std::endl;

  std::cout << f.getParameterList() << std::endl << std::endl;

  return 0;
}
