#ifndef MUELU_PARAMLISTTEST_HPP
#define MUELU_PARAMLISTTEST_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::tuple;

// Printing the valid parameter list only showing documentation fields
inline void printParameterListOptions(std::ostream &os, const Teuchos::ParameterList & p) {
  p.print(os, Teuchos::ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));
  os << std::endl;
}

// Abstract interface of a class accepting parameter lists.
class ParameterListAcceptor {

public:

  ParameterListAcceptor() { }

  virtual ~ParameterListAcceptor() { }

  // Return a const parameter list of valid parameters that setParameterList() will accept.
  // Also define the default values of parameters according to the input parameter list.
  //
  // Questions:
  // - Do we want this function to be static?
  // => yes, but we also want it virtual. So static is not possible. But this method does not access any instance specific data.
  // - Do we want to pass an optional param list as input parameter?
  //   yes => if some parameters depends of other parameters, this would allow to return different default / valid parameters according to already set parameters.
  //   Ex: fact: ILUT => extra param threshold; fact: ILUK => extra param level-of-fill
  //
  // If a parameter is unused AND default, it is printed as [default] by std::cout << paramList but printed as unused by paramList.unused(std::cout).
  // So if we set default parameters in getValidParameters() that are not used, user will get a warning message. We don't want that for [default].
  // One solution is to never set any unused parameter in getValidParameters().
  // If some parameters are available only conditionnaly, do not set them by default in setValidParameters when the conditions are not meet.
  //
  virtual RCP<const ParameterList> getValidParameters(const ParameterList& pL = ParameterList()) const = 0;

  // Set parameters from a parameter list and return with default values.
  //
  // Questions:
  //  - Input: const ParameterList or not const? IE: do we want to modify user paramlist directly?
  //    => not const: user can make a copy outside if he do not want us to modify his initial paramlist.
  //  - Input: RCP or just reference? RCP avoid the copy of the parameter list but
  //           I'm afraid that a user modify the list outside of the class after doing a setParameterList
  //           (this by-pass the validation and can create some confusion).
  //           In another hand, if we use a copy of the list, we do not set the "[used]"/["unused"] flag
  //           on the original list during Build(). getParameterList has to be used to retrieve the list with the correct flags.
  //
  // What we really want is for the user to have a const version outside and MueLu having a non-const version inside.
  // Is there a C++ pattern to do that?
  virtual void setParameterList(ParameterList & paramList) = 0;

  // Add a parameter to the current parameter list
  // virtual void setParameter(ParameterEntry)

  // Get the parameter list.
  // Teuchos::ParameterListAcceptor also provides a non-const version of this but I don't see why.
  // We don't want a user to mess with the internal parameter list.
  // To change a parameter, one can make a copy and call setParameterList again.
  // No need for RCP, it's just a view
  virtual const Teuchos::ParameterList & getParameterList() const = 0;

  virtual void getDocumentation(std::ostream &os) const = 0;

};


// Partial implementation of ParameterListAcceptor that stores the object parameters in an internal parameterList
//
// Another possibility is to copy the parameters to internal variables (as Ifpack) ?
// => I don't think it's a good solution because:
//      - we need to rebuild a parameter list for getParameterList()
//      - we will read all the parameters even if they are unused (wrong "[used]"/["unused"] flag)
// double paramA_;
// double paramB_;
class ParameterListAcceptorImpl: public ParameterListAcceptor {

public:

  ParameterListAcceptorImpl() { }

  virtual ~ParameterListAcceptorImpl() {
    bool warnings = true;
    if (warnings) {
      paramList_.unused(std::cout);
    }
  }

  virtual void setParameterList(ParameterList & paramList) {
    // Validate and add defaults parameters.
    paramList.validateParametersAndSetDefaults(*getValidParameters(paramList));

    // Do we need to validate range of int/double parameters here? Can it be done automatically as for valid std::string in getValidParameters()?
    paramList_ = paramList; // copy
  }


  virtual const Teuchos::ParameterList & getParameterList() const {
    return paramList_;
  }

  virtual void getDocumentation(std::ostream &os) const {
    // default implementation

    os << "## Parameters:" << std::endl;
    printParameterListOptions(os, *getValidParameters());

    os << "## Fully described default method:" << std::endl;
    getValidParameters()->print(os, 2, true, false);
    os << std::endl;
  }

protected:

  // Note: Teuchos::ParameterListAcceptor provides function to access this field on sub-classes (getMyParamList()). We might want to do that as well.
  Teuchos::ParameterList paramList_;

};

#endif // MUELU_PARAMLISTTEST_HPP
