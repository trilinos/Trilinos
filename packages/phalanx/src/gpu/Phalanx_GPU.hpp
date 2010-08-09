#ifndef PHALANX_GPU_HPP
#define PHALANX_GPU_HPP

#include <iostream>

namespace PHX {
  void ListDevices(std::ostream& os);
  void EvaluateGPU();
}

#endif
