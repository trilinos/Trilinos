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
///
/// This class is only allowed to get environment variable values; it
/// is not allowed to set them.  This is in part because many
/// third-party run-time libraries only read environment variables
/// once, at some time outside of Tpetra's control.  Thus, we want to
/// discourage Tpetra developers from assuming that they can set
/// environment variables in order to control other libraries.  We
/// also want to discourage use of environment variables as a means of
/// communication between software components.
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

  /** @brief Get Boolean value of the given environment variable.
   *
   *  This is nonconst because the implementation reserves the right
   *  to cache the results.
   *
   *  @param variableName the name of the environment variable
   *
   *  @param defaultValue [optional, false] the default value to
   *    return in the case that the environment variable does not
   *    exist.
   *
   *  @return If the variable does not exist, or if the upper-case
   *    version of it is 0, "NO", or "FALSE", return false.
   *    Otherwise, return true.
   */
  bool
  getBooleanValue (const std::string& variableName,
                   const bool defaultValue = false);

  /** @brief Get value of the given environment variable as a string.
   *
   *  This is nonconst because the implementation reserves the right
   *  to cache the results.
   *
   *  @param variableName [in] Name of the environment variable.  This
   *    is case sensitive.
   *
   *  @param defaultValue [in, optional, ""] Default value to return
   *    in the case that the environment variable does not exist.
   *
   *  @return If the environment variable exists, its value as a
   *    string; otherwise, defaultValue.
   *
   */
  std::string getValue(const std::string& variableName,
                       const std::string& defaultValue="");

  void clearCache();
};
}  // namespace Details
}  // namespace Tpetra
#endif  // _TPETRA_DETAILS_ENVIRONMENT_HPP_
