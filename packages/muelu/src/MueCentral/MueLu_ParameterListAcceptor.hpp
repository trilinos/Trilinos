#ifndef MUELU_PARAMETERLISTACCEPTOR_HPP
#define MUELU_PARAMETERLISTACCEPTOR_HPP

#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// TODO See also: Teuchos::ParameterListAcceptor, Teko::Clonable

namespace MueLu {

using Teuchos::ParameterEntry;

// Printing the valid parameter list only showing documentation fields
void printParameterListOptions(std::ostream &os, const Teuchos::ParameterList &p);

/*! @class ParameterListAcceptor
    @brief Abstract interface of a class accepting parameter lists.

*/
class ParameterListAcceptor {
 public:
  //! @name Constructors/Destructors
  //! @{
  ParameterListAcceptor();

  virtual ~ParameterListAcceptor() = default;

  //! @}
  //! @name Get/Set methods.
  //! @{

  /*! @brief Return a const parameter list of valid parameters that setParameterList() will accept.

     Also define the default values of parameters according to the input parameter list.
  */
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
  // If some parameters are available only conditionnaly, do not set them by default in setValidParameters when the conditions are not met.
  //
  virtual Teuchos::RCP<const Teuchos::ParameterList> GetValidParameterList() const = 0;

  //! @brief Set parameters from a parameter list and return with default values.
  //
  // Questions:
  //  - Input: const ParameterList or not const? IE: do we want to modify user paramlist directly?
  //    => not const: user can make a copy outside if he do not want us to modify his initial paramlist.
  //  - Input: RCP or just reference? RCP avoid the copy of the parameter list but
  //           I'm afraid that a user modify the list outside of the class after doing a SetParameterList
  //           (this by-pass the validation and can create some confusion).
  //           In another hand, if we use a copy of the list, we do not set the "[used]"/["unused"] flag
  //           on the original list during Build(). GetParameterList has to be called to retrieve the list with the correct flags.
  //
  // What we really want is for the user to have a const version outside and MueLu having a non-const version inside.
  // Is there a C++ pattern to do that?
  virtual void SetParameterList(const Teuchos::ParameterList &paramList) = 0;

  // Get the parameter list.
  // Teuchos::ParameterListAcceptor also provides a non-const version of this but I don't see why.
  // We don't want a user to mess with the internal parameter list.
  // To change a parameter, one can make a copy and call SetParameterList again.
  // No need for RCP, it's just a view
  virtual const Teuchos::ParameterList &GetParameterList() const = 0;

  //! @brief Set a parameter directly as a ParameterEntry.
  virtual void SetParameter(const std::string &name, const ParameterEntry &entry) = 0;

  //! @brief Retrieves a const entry with the name name.
  virtual const ParameterEntry &GetParameter(const std::string &name) const = 0;

  virtual void GetDocumentation(std::ostream &os) const = 0;

  //! @}

};  // class ParameterListAcceptor

// Partial implementation of ParameterListAcceptor that stores the object parameters in an internal parameterList
class ParameterListAcceptorImpl : public ParameterListAcceptor {
 public:
  ParameterListAcceptorImpl();

  virtual ~ParameterListAcceptorImpl() = default;

  virtual void SetParameterList(const Teuchos::ParameterList &paramList);

  // The returned list always has an entry for each valid parameter.
  // Therefore, there is not need to test if a parameter is present before getting it.
  virtual const Teuchos::ParameterList &GetParameterList() const;

  void SetParameter(const std::string &name, const ParameterEntry &entry);

  const ParameterEntry &GetParameter(const std::string &name) const;

  virtual void GetDocumentation(std::ostream &os) const;

 private:
  mutable                                 // mutable can be avoid by changing return type of GetParameterList() to RCP but conceptually, I like the fact that GetParameterList returns ref to param list.
      Teuchos::ParameterList paramList_;  // Private: Use GetParameterList() to access this list
                                          // This list might be empty before calling GetParameterList()
                                          // but it is populate with automatically if needed in SetParameterList/GetParameterList/...
                                          // Therefore, there is no need to test if a parameter exist before accessing it with pL.get("paramName")
                                          // in sub classes.

  // Note: Teuchos::ParameterListAcceptor has a getMyParamList() function to access paramList_ without going through the logic of getParameterList()
  // Do we need that too?
};

}  // namespace MueLu

#endif  // MUELU_PARAMETERLISTACCEPTOR_HPP
