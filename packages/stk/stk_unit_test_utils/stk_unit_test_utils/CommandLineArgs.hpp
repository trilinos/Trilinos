#pragma once

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

}  // namespace unit_test_util
}  // namespace stk
