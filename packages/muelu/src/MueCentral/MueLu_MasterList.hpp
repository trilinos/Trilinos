// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MASTERLIST_HPP
#define MUELU_MASTERLIST_HPP

#include <sstream>
#include <map>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_any.hpp>

namespace MueLu {

/*!
  @class MasterList class.
  @brief Static class that holds the complete list of valid MueLu parameters.

  This class creates a ParameterList that is used to validate user-supplied ParameterLists.  This same list
  is the source of default values if a parameter does not appear in the user ParameterList.
  This class also provides ParameterLists for certain common problem types, such as 2D and 3D isotropic Poisson.
  @todo Add method to register user-defined problem type.  This will need both string identifier
   and ParameterList (or string version of parameter list).

*/

template <typename T, typename U>
class DefaultProblemStrings;

class MasterList {
 public:
  //! @brief Return a "master" list of all valid parameters and their default values.
  static Teuchos::RCP<const Teuchos::ParameterList> List();
  //! @brief Return default parameter settings for the specified problem type.
  static Teuchos::RCP<Teuchos::ParameterList> GetProblemSpecificList(std::string const& problemType);

  //! @brief Returns default value on the "master" list for a parameter with the specified name and type.
  template <typename T>
  static const T& getDefault(const std::string& name) {
    return List()->get<T>(name);
  }

  //! @brief Returns default entry from the "master" list corresponding to the specified name.
  static const Teuchos::ParameterEntry& getEntry(const std::string& name) {
    return List()->getEntry(name);
  }

  //! @brief Create xml string for given MueLu parameter (easy xml input format)
  //!
  //! @note: We should check whether template type T is the same as the expected parameter type in the parameter list
  template <typename T>
  static std::string generateXMLParameterString(const std::string& name, T data) {
    Teuchos::ParameterEntry entry = getEntry(name);
    std::stringstream ss;
    // TODO: check whether getAny().typeName() matches the typeid(T)
    //       One could use a technique as described here: http://stackoverflow.com/questions/4484982/how-to-convert-typename-t-to-string-in-c
    ss << "<Parameter name=\"" << name << "\" type=\"" << entry.getAny().typeName() << "\" value=\"" << data << "\"/>";
    return ss.str();
  }

  //! @brief Translate ML parameter to corresponding MueLu parameter
  static std::string ML2MueLu(const std::string& name) {
    std::map<std::string, std::string>::iterator it;
    it = ML2MueLuLists_.find(name);
    if (it != ML2MueLuLists_.end())
      return it->second;
    else
      return "";
  }

  static std::string interpretParameterName(const std::string& name, const std::string& value);

 private:
  MasterList();
  MasterList(const MasterList&);
  MasterList& operator=(const MasterList&);

  //! @brief A ParameterList that holds all valid parameters and their default values.
  static Teuchos::RCP<Teuchos::ParameterList> masterList_;
  //! @brief String equivalent of the masterList_.
  static const std::string stringList_;
  //! @brief A ParameterList that holds all valid parameters and their default values for a particular problem type.
  static Teuchos::RCP<Teuchos::ParameterList> problemSpecificList_;
  //! @brief The problem type associated with the current problem-specific ParameterList.
  static std::string problemType_;
  //! @brief Map of string equivalents of the problemSpecificList_.  The first entry is the problem type, the second is the string equivalent.
  static std::map<std::string, std::string> DefaultProblemTypeLists_;
  //! @brief Map of ML parameter strings to corresponding MueLu parametes
  static std::map<std::string, std::string> ML2MueLuLists_;
};

/*!
  @class DefaultProblemStrings class.
  @brief Helper class to initialize DefaultProblemTypeLists_ in class MasterList.
*/
template <typename T, typename U>
class DefaultProblemStrings {
 public:
  DefaultProblemStrings(const T& key, const U& val) {
    map_[key] = val;
  }

  DefaultProblemStrings<T, U>& operator()(const T& key, const U& val) {
    map_[key] = val;
    return *this;
  }

  operator std::map<T, U>() const {
    return map_;
  }

 private:
  std::map<T, U> map_;
};

}  // namespace MueLu

#endif
