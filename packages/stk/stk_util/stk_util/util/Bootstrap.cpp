#include <stk_util/util/Bootstrap.hpp>

namespace stk {

Bootstrap *
Bootstrap::s_front = 0;

bool
Bootstrap::s_bootstrapped = false;


void
Bootstrap:: bootstrap()
{
  s_bootstrapped = true;
  for (Bootstrap *f = s_front; f; f = f->m_next)
    (*f->m_f)();
}


Bootstrap::Bootstrap(
  FunctionPtr           f)
  : m_next(s_front),
    m_f(f)
{
  s_front = this;

  // IF already bootstrapped, execute immediately
  if (s_bootstrapped)
    (*f)();
}

} // namespace stk
