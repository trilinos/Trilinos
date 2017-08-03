#include <map>
#include <string>
#include <list>

#ifndef _TPETRA_DETAILS_ENVIRONMENT_HPP_
#define _TPETRA_DETAILS_ENVIRONMENT_HPP_

namespace Tpetra {
namespace Details {

class DefaultEnvironmentVariables {
  /*
   * This class merely provides a list of variables that should be cached on
   * instantiation of the singleton Environment().  This can be moved and made
   * more general to include default variables for other packages.
   *
   */
 public:
  static std::list<std::string> getDefaults();
};

/// \class Environment
/// \brief Access environment variables, and cache their values.
///
/// "Environment" here refers to the operating system or shell.  For
/// example, PATH is a commonly used environment variable.  Tpetra can
/// read any environment variable, but it would be wise to restrict
/// yourself to variables that start with the prefix "TPETRA_".
///
/// First, get the singleton instance of this class:
/// \code
/// Environment& env = Environment::getInstance ();
/// \endcode
/// Then, use it to read the environment:
/// \code
/// const bool exists = env.variableExists ("TPETRA_DEBUG");
/// const bool debug = exists && env.getBooleanValue ("TPETRA_DEBUG");
/// const std::string val =
///   env.getValues ("TPETRA_NAME_OF_THING", "DEFAULT_NAME");
/// \endcode
class Environment {
public:
  //! Get singleton instance.
  static Environment& getInstance();

private:
  Environment ();

  std::map<std::string, const char *> environCache_;
  void cacheVariable(const std::string& variableName,
                     const char* variableValue);

 public:
  /* Don't allow any copies of the singleton to appear... (requires CXX11)
   * Note: Scott Meyers mentions in his Effective Modern
   *       C++ book, that deleted functions should generally
   *       be public as it results in better error messages
   *       due to the compilers behavior to check accessibility
           before deleted status
  */
  Environment(Environment const&)    = delete;
  void operator=(Environment const&) = delete;

  /** @brief Checks if the environment variable variableName is cached
   *
   *  @param variableName the name of the environment variable
   *  @return whether variableName is cached
   */
  bool variableIsCached(const std::string& variableName);

  /** @brief Checks for existence of the environment variable variableName
   *
   *  This is nonconst because the implementation reserves the right to cache
   *  results.
   *
   *  @param variableName the name of the environment variable
   *  @return whether variableName exists
  */
  bool variableExists(const std::string& variableName);

  /** @brief Gets value of the given environment variable.
   *
   *  This is nonconst because the implementation reserves the right to
   *  cache the results.
   *
   *  @param variableName the name of the environment variable
   *  @param defaultValue [optional, ""] the default value to return in the case
   *         that the environment variable variableName does not exist
   *
   *  @return false if the variable does not exist, is one of [0, NO, FALSE],
   *          else true
   *
   */
  bool getBooleanValue(const std::string& variableName,
                       const std::string& defaultValue="");

  /** @brief Gets value of the given environment variable.
   *
   *  This is nonconst because the implementation reserves the right to
   *  cache the results.
   *
   *  @param variableName the name of the environment variable
   *  @param defaultValue [optional, ""] the default value to return in the case
   *         that the environment variable variableName does not exist
   *
   *  @return the value of the environment variable (or default if it does not
   *          exist)
   *
   */
  std::string getValue(const std::string& variableName,
                       const std::string& defaultValue="");

  // /* @brief Sets value of the given environment variable.
  //  *
  //  *  @param variableName the name of the environment variable
  //  *  @param variableValue the value of the environment variable
  //  *  @param overwrite [optional, 1] overwrite the variable if it already exists
  //  *
  //  *  @return void
  //  *
  //  */
  // void setValue(const std::string& variableName,
  //               const std::string& variableValue,
  //               const int overwrite=1);


  void clearCache();
};
}  // namespace Details
}  // namespace Tpetra
#endif  // _TPETRA_DETAILS_ENVIRONMENT_HPP_
