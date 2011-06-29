#include <sstream>
#include <fstream>
#include "MueLu_Memory.hpp"

#include <iostream> // TODO: remove
#include <unistd.h>
#include <time.h>
#include <malloc.h>

//#define MUELU_USE_MALLINFO

namespace MueLu {
  
  namespace MemUtils {
    
    std::string PrintMemoryUsage() {


#ifdef MUELU_USE_MALLINFO
      struct mallinfo mem_stats = mallinfo();
      double memory = mem_stats.hblkhd + mem_stats.usmblks + mem_stats.uordblks;

      char memchar[128];
      sprintf(memchar,"%12.1f MB",memory/1048576.0);
      std::string mem(memchar);

      return mem;
#else
      std::ostringstream mem;
      std::ifstream proc("/proc/self/status");
      std::string s;

      mem << PrintMemoryInfo() << " ";
      while(getline(proc, s), !proc.fail()) {
        if(s.substr(0, 6) == "VmSize") {
          mem << s;
          return mem.str();
        }
      }      
      return mem.str();
#endif

    }

    std::string PrintMemoryInfo() {

#ifdef MUELU_USE_MALLINFO
      struct mallinfo mem_stats = mallinfo();
      double memory = mem_stats.hblkhd + mem_stats.usmblks + mem_stats.uordblks;

      char memchar[128];
      sprintf(memchar,"%12.1f MB",memory/1048576.0);
      std::string mem(memchar);

      return mem;
#else
      std::ostringstream mem;
      std::ifstream proc("/proc/meminfo");
      std::string s;
      while(getline(proc, s), !proc.fail()) {
        if(s.substr(0, 7) == "MemFree") {
          mem << s;
          return mem.str();
        }
    
      }
      return mem.str();
#endif
    }

    void ReportTimeAndMemory(Teuchos::Time const &timer, Teuchos::Comm<int> const &Comm)
    {
      double maxTime=0,minTime=0,avgTime=0;
      double localTime = timer.totalElapsedTime();
      int ntimers=1, root=0;
      MPI_Reduce(&localTime,&maxTime,ntimers,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD);
      MPI_Reduce(&localTime,&minTime,ntimers,MPI_DOUBLE,MPI_MIN,root,MPI_COMM_WORLD);
      MPI_Reduce(&localTime,&avgTime,ntimers,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
      avgTime /= Comm.getSize();
      //std::cout << "(" << Comm.getRank() << ") " << localTime << std::endl; 
      if (Comm.getRank()==0) {
        std::cout << "&&& " << timer.name() << " time    &&& "
                 << "max=" << maxTime << "   min=" << minTime << "  avg=" << avgTime << std::endl;
        std::cout << "&&& " << timer.name() << " memory  &&& " << MemUtils::PrintMemoryUsage() << std::endl;
      }
    } //ReportTimeAndMemory
    
  } //namespace MemUtils
  
} //namespace MueLu
