#include "Teuchos_stacktrace.hpp"

void g()
{
    Teuchos::show_backtrace();
}

void f()
{
    g();
}

int main()
{
    Teuchos::print_stack_on_segfault();
    f();

    // This will segfault:
    char *p = NULL; *p = 0;

    return 0;
}
