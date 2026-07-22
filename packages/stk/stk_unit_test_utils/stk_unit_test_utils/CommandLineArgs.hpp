#pragma once

#include <cstddef>
#include <vector>
#include <string>

namespace stk
{
namespace unit_test_util
{
class GlobalCommandLineArguments
{
 public:
  static GlobalCommandLineArguments& self()
  {
    static GlobalCommandLineArguments s;
    return s;
  }
  void set_values(int argc, char** argv)
  {
    m_argc = argc;
    m_argv = argv;
  }
  int& get_argc() { return m_argc; }
  char**& get_argv() { return m_argv; }

 private:
  int m_argc = 0;
  char** m_argv = nullptr;
};

class Args
{
public:
  Args(const std::vector<std::string> & strArgs)
    : m_stringArgs(strArgs),
      m_argc(m_stringArgs.size()),
      m_argv(strArgs.empty() ? nullptr : new const char*[m_argc])
  {
    for (int i = 0; i < m_argc; ++i) {
      m_argv[i] = const_cast<char*>(m_stringArgs[i].c_str());
    }
  }

  ~Args()
  {
    delete [] m_argv;
  }

  int argc() { return m_argc; }
  const char** argv() { return m_argv; }
  char** nonconst_argv() { return const_cast<char**>(argv()); }

private:
  const std::vector<std::string> m_stringArgs;
  int m_argc;
  const char** m_argv;
};

}  // namespace unit_test_util
}  // namespace stk
