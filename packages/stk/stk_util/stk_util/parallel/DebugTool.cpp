#include "DebugTool.hpp"
#include <iostream>
#include <unistd.h>

int wait_for_debugger()
{
  volatile int i = 0;
  int j = 0;
  int pid = getpid();
  std::cout << "Waiting for debugger to attach to process " << pid << ".  Once the debugger is attached, do \"set var i = 1\" in gdb (or similar syntax in your debugger)" << std::endl;
  while (i == 0)
  {
    j++;
  }

  return j;
}