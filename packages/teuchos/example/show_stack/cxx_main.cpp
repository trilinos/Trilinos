#include "Teuchos_stacktrace.hpp"


void g()
{
#ifdef HAVE_TEUCHOS_STACKTRACE
    Teuchos::show_stacktrace();
#endif
}


void f()
{
    g();
}


int main()
{
#ifdef HAVE_TEUCHOS_STACKTRACE
    Teuchos::print_stack_on_segfault();
#endif
    f();

    // This will segfault:
    char *p = NULL; *p = 0;

    return 0;
}
