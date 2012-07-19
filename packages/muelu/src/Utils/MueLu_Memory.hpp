#ifndef MUELU_MEMORY_HPP
#define MUELU_MEMORY_HPP

#include <string>

#include <Teuchos_Comm.hpp>
#include <Teuchos_Time.hpp>
#include "MueLu_config.hpp"

namespace MueLu {
  
  namespace MemUtils {
    
    std::string PrintMemoryUsage();
    std::string PrintMemoryInfo();
    void ReportTimeAndMemory(Teuchos::Time const &timer, Teuchos::Comm<int> const &Comm);
    
  } //namespace MemUtils
  
} //namespace MueLu

#endif //ifndef MUELU_MEMORY_HPP
