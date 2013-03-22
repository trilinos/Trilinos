#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <locale>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include <stk_percept/Util.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <sys/resource.h>

// FIXME
    double s_timers[10] = {0,0,0,0,0,0,0,0,0,0};




namespace shards {

    std::ostream& operator<<(std::ostream& os, const shards::Array<double, shards::NaturalOrder>& container) 
    {
      // Save the format state of the original ostream os.
      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(os);  

      os.setf(std::ios_base::scientific, std::ios_base::floatfield);
      os.setf(std::ios_base::right);
      int myprec = os.precision();
  
      int size = container.size();
      int rank = container.rank();
      Teuchos::Array<int> multiIndex(rank);
      //Teuchos::Array<int> dimensions(rank);
      std::vector<int> dimensions(rank);
      //container.dimensions(dimensions);
      for (int irank = 0; irank < rank; irank++)
        {
          dimensions[irank] = container.dimension(irank);
        }
  
      os<< "===============================================================================\n"\
        << "\t Container size = " << size << "\n"
        << "\t Container rank = " << rank << "\n" ;
  
      if( (rank == 0 ) && (size == 0) ) {
        os<< "====================================================================================\n"\
          << "|                        *** This is an empty container ****                       |\n";
      }
      else {
        os<< "\t Dimensions     = ";
    
        for(int r = 0; r < rank; r++){
          os << " (" << dimensions[r] <<") ";
        }
        os << "\n";
    
        os<< "====================================================================================\n"\
          << "|              Multi-index          Enumeration             Value                  |\n"\
          << "====================================================================================\n";
      }
  
      int address=0;
      
      std::vector<int> dims = dimensions;
      dims.resize(8,1);
      std::vector<int> idim(8);
      for (idim[0] = 0; idim[0] < dims[0]; idim[0]++)
        {
          for (idim[1] = 0; idim[1] < dims[1]; idim[1]++)
            {
              for (idim[2] = 0; idim[2] < dims[2]; idim[2]++)
                {
                  for (idim[3] = 0; idim[3] < dims[3]; idim[3]++)
                    {
                      for (idim[4] = 0; idim[4] < dims[4]; idim[4]++)
                        {
                          for (idim[5] = 0; idim[5] < dims[5]; idim[5]++)
                            {
                              for (idim[6] = 0; idim[6] < dims[6]; idim[6]++)
                                {
                                  std::ostringstream mistring;
                                  for (int jr = 0; jr < rank; jr++)
                                    {
                                      mistring << idim[jr] << std::dec << " " ;
                                    }
                                  os.setf(std::ios::right, std::ios::adjustfield);
                                  os << std::setw(27) << mistring.str(); 
                                  os << std::setw(20) << address;
                                  os << "             ";
                                  os.setf(std::ios::left, std::ios::adjustfield);
                                  os << std::setw(myprec+8) 
                                     << container[address]
                                     << "\n";
                                  ++address;
                                }
                            }
                        }
                    }
                }
            }
        }
#if 0
      for(int address = 0; address < size; address++){
        container.getMultiIndex(multiIndex,address);
        std::ostringstream mistring;
        for(int r = 0; r < rank; r++){
          mistring <<  multiIndex[r] << std::dec << " "; 
        }
        os.setf(std::ios::right, std::ios::adjustfield);
        os << std::setw(27) << mistring.str(); 
        os << std::setw(20) << address;
        os << "             ";
        os.setf(std::ios::left, std::ios::adjustfield);
        os << std::setw(myprec+8) << container[address] << "\n";
      }
#endif
  
      os<< "====================================================================================\n\n";

      // reset format state of os
      os.copyfmt(oldFormatState);

      return os;
    }

}

namespace stk { 

  namespace mesh { 

    std::ostream &operator<<(std::ostream& out, const stk::mesh::Entity& entity)
    {
      if (entity.entity_rank() != stk::mesh::fem::FEMMetaData::NODE_RANK)
        {
          out << "Elem: " << entity.identifier() << " rank= " << entity.entity_rank() << " nodes: ";

          const mesh::PairIterRelation elem_nodes = entity.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity & node = * elem_nodes[ inode ].entity();

              out << node.identifier() << " ";
            }
          //out << std::endl;

        }

      else if (entity.entity_rank() == stk::mesh::fem::FEMMetaData::NODE_RANK)
        {
          out << "Node: " << entity.identifier();

        }
      else 
        {
          out << "rank unknown: " << entity.entity_rank();
        }

      return out;
    }    

  }

  namespace percept { 

#define ENABLE_PAUSE 1
    //#define ENABLE_PAUSE (1 && NDEBUG)

    double Util::s_timers[10] = {0,0,0,0,0,0,0,0,0,0};

#if ENABLE_PAUSE
    static bool s_doPause = true;
    static bool s_doPauseNeverAgain = false;  // set to false to enable pausing
#else
    static bool s_doPause = false;
    static bool s_doPauseNeverAgain = true;  // set to false to enable pausing
#endif
    static std::map<int, bool> s_flag;
    static int s_p_rank=0;
    void Util::setFlag(int which_flag, bool val) { s_flag[which_flag]=val; }
    bool Util::getFlag(int which_flag) { return s_flag[which_flag]; }
    void Util::setRank(int rank){ s_p_rank=rank; }
    int Util::get_rank() { return s_p_rank; }

    void Util::setDoPause(bool enableOrDisablePause)
    {
      s_doPause = enableOrDisablePause;
    }
    void Util::pause(bool doPause, const char *msg)
    {
      if (s_doPause && doPause && !s_doPauseNeverAgain)
        {
          if (msg)
            {
              std::cout << msg << std::endl;
            }

          std::cout << "pause...:: c or <any-key> to continue, q to quit, C to continue without further pauses: " << std::endl;
          char buf[100];
          std::cin >> buf;
          std::string buf_s(buf);
          if (buf_s == "q")
            {
              std::exit(1);
            }
          //             else if (buf_s == "c")
          //               {
          //               }
          else if (buf_s == "C")
            {
              s_doPauseNeverAgain = true;
            }

        }

    }
    void Util::debug_stop()
    {
      std::cout << "in debug_stop" << std::endl;
    }

    static unsigned s_trace_mem_0[] = {0u,0u,0u,0u,0u,0u,0u,0u,0u,0u};
    static unsigned s_trace_mem_1[] = {0u,0u,0u,0u,0u,0u,0u,0u,0u,0u};
    static double s_trace_time_0[] = {0,0,0,0,0,0,0,0,0,0};
    static double s_trace_time_1[] = {0,0,0,0,0,0,0,0,0,0};

    void 
    Util::
    trace_cpu_time_and_mem_0(unsigned index)
    {
        {
          size_t heap_in_bytes = 0;
          size_t memory_in_bytes = Util::memory(heap_in_bytes);

          double cpu_in_min = Util::cpu_time()/60.0;
        
          s_trace_mem_0[index] = memory_in_bytes;
          s_trace_time_0[index] = cpu_in_min;
        }
    }

    void 
    Util::
    trace_cpu_time_and_mem_1(unsigned index)
    {
        {
          size_t heap_in_bytes = 0;
          size_t memory_in_bytes = Util::memory(heap_in_bytes);

          double cpu_in_min = Util::cpu_time()/60.0;
        
          s_trace_mem_1[index] += memory_in_bytes - s_trace_mem_0[index];
          s_trace_time_1[index] += cpu_in_min - s_trace_time_0[index];
        }
    }


    void Util::trace_cpu_time_and_mem_print(int index, std::string msg)
    {
      std::cout << "tmp trace_cpu_time_and_mem_print "              
                << msg << " = " << ((double)s_trace_mem_1[index])/(1024.0*1024.0) << " [Mb] " 
                << s_trace_time_1[index] << " [min] " << std::endl ;
    }

    void Util::replace(std::string &str, const std::string &find_what, const std::string &replace_with)
    {
      std::string::size_type pos = 0;
      while((pos = str.find(find_what, pos)) != std::string::npos)
        {
          str.erase(pos, find_what.length());
          str.insert(pos, replace_with);
          pos += replace_with.length();
        }
    }



static void
get_memory_info(
  size_t &		memory_usage,
  size_t &		faults)
{
  memory_usage = 0;
  faults = 0;

  std::ifstream proc("/proc/self/stat", std::ios_base::in|std::ios_base::binary);
  if (proc) {

    std::string s;
    int i;
    for (i = 0; i < 11; ++i)
      proc >> s;

    proc >> faults;
    ++i;

    for (; i < 22; ++i)
      proc >> s;

    proc >> memory_usage;
    ++i;
  }
}

static void
get_heap_info(
  size_t &		heap_size,
  size_t &		largest_free)
{
  heap_size = 0;
  largest_free = 0;

#if defined(SIERRA_HEAP_INFO)
# if 0 // if defined(REDS) // Redstorm now links in gnu's malloc
  static size_t reds_fragments;
  static unsigned long reds_total_free;
  static unsigned long reds_heap_size;
  static unsigned long reds_largest_free;

  ::heap_info(&reds_fragments, &reds_total_free, &reds_largest_free, &reds_heap_size);

  heap_size = reds_heap_size;
  largest_free = reds_largest_free;

  slibout.m(Slib::LOG_MEMORY) <<"reds_fragments " << reds_fragments
			      << ", reds_total_free " << reds_total_free
			      << ", reds_largest_free " << reds_largest_free
			      << ", reds_heap_size " << reds_heap_size << Diag::dendl;

// # elif defined(__linux__)
# elif defined(__linux__) || defined(REDS) && ! defined(__IBMCPP__)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = (unsigned int) minfo.uordblks + (unsigned int) minfo.hblkhd;
  largest_free = (unsigned int) minfo.fordblks;


# elif defined(__sun)
  pstatus_t proc_status;

  std::ifstream proc("/proc/self/status", std::ios_base::in|std::ios_base::binary);
  if (proc) {
    proc.read((char *)&proc_status, sizeof(proc_status));
    heap_size = proc_status.pr_brksize;
    slibout.m(Slib::LOG_MEMORY) <<"pr_brksize " << proc_status.pr_brksize
				<< ", pr_stksize " << proc_status.pr_stksize << Diag::dendl;
  }
# endif
#endif // defined(SIERRA_HEAP_INFO)
}



    double
    Util::cpu_time()
    {
#if defined(REDS)
      struct rusage my_rusage;

      ::getrusage(RUSAGE_SELF, &my_rusage);

      double seconds = my_rusage.ru_utime.tv_sec;
      double micro_seconds = my_rusage.ru_utime.tv_usec;
  
      return seconds + micro_seconds*1.0e-6;

#else
      static struct rusage my_rusage;

      ::getrusage(RUSAGE_SELF, &my_rusage);

      double seconds = (double)(my_rusage.ru_utime.tv_sec + my_rusage.ru_stime.tv_sec);
      double micro_seconds = (double)(my_rusage.ru_utime.tv_usec + my_rusage.ru_stime.tv_usec);
      return seconds + micro_seconds*1.0e-6;

#endif
    }

    size_t
    Util::memory(size_t& heap)
    {
#if defined(REDS)
      struct rusage my_rusage;

      ::getrusage(RUSAGE_SELF, &my_rusage);

      double seconds = my_rusage.ru_utime.tv_sec;
      double micro_seconds = my_rusage.ru_utime.tv_usec;
  
      return seconds + micro_seconds*1.0e-6;

#else

      /* Maximum resident set size (in kilobytes).  */
      size_t largest_free = 0;
      size_t faults = 0;
      size_t memory_in_bytes = 0;
      get_heap_info(heap, largest_free);
      get_memory_info(memory_in_bytes, faults);
      return memory_in_bytes;

#endif
    }

    //========================================================================================================================

    // FIXME
    bool Util::isLinearElement(shards::CellTopology& cell_topo)
    {
      if (cell_topo.getVertexCount() == cell_topo.getNodeCount())
        return true;
      else
        return false;
    }



    //========================================================================================================================
    stk::diag::TimerSet& perceptTimerSet()
    {
      static stk::diag::TimerSet s_perceptTimerSet(PERCEPT_TIMER_ROOT);

      return s_perceptTimerSet;
    }

    stk::diag::Timer& perceptTimer() {
      const std::string name("PerceptTimer");
      static stk::diag::Timer s_perceptTimer (stk::diag::createRootTimer(name, perceptTimerSet()));

      return s_perceptTimer;
    }

    LapTimeType getLapTime(stk::diag::Timer& lap_timer) { return lap_timer.getMetric<stk::diag::WallTime>().getLap(); }
    LapCountType getAccumulatedLap(stk::diag::Timer& timer, bool option) { return timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(option); }



  }
}
