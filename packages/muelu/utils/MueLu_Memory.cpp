#include <sstream>
#include <fstream>
#include "MueLu_Memory.hpp"

#include <iostream> // TODO: remove

namespace MueLu {
  
  namespace MemUtils {
    
    std::string PrintMemoryUsage() {
      std::ostringstream mem;
      std::ifstream proc("/proc/self/status");
      std::string s;

      //      std::cout << "Free memory of the node:" << PrintMemoryInfo() << std::endl; //TODO: remove
      mem << PrintMemoryInfo() << " ";

      while(getline(proc, s), !proc.fail()) {
	if(s.substr(0, 6) == "VmSize") {
	  mem << s;
	  return mem.str();
	}
      }
      
      return mem.str();
    }

    std::string PrintMemoryInfo() {
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
