#include <iostream>
#include <string>
#include <cstdlib>

#include <stk_percept/Util.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <sys/resource.h>

    double s_timers[10] = {0,0,0,0,0,0,0,0,0,0};


double
srk_cpu_time()
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

  double seconds = my_rusage.ru_utime.tv_sec + my_rusage.ru_stime.tv_sec;
  double micro_seconds = my_rusage.ru_utime.tv_usec + my_rusage.ru_stime.tv_usec;
  
  return seconds + micro_seconds*1.0e-6;

#endif
}


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
      if (entity.entity_rank() == mesh::Element || entity.entity_rank() == mesh::Face || entity.entity_rank() == mesh::Edge  )
        {
          out << "Elem: " << entity.identifier() << " rank= " << entity.entity_rank() << " nodes: ";

          const mesh::PairIterRelation elem_nodes = entity.relations( mesh::Node );
          unsigned num_node = elem_nodes.size();
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity & node = * elem_nodes[ inode ].entity();

              out << node.identifier() << " ";
            }
          //out << std::endl;

        }

      else if (entity.entity_rank() == mesh::Node)
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
    int Util::getRank() { return s_p_rank; }

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

    // FIXME
    bool Util::isLinearElement(shards::CellTopology& cell_topo)
    {
      if (cell_topo.getVertexCount() == cell_topo.getNodeCount())
        return true;
      else
        return false;
    }

    //========================================================================================================================

    void Util::printEntity(std::ostream& out, const stk::mesh::Entity& entity, stk::mesh::FieldBase* field)
    {
      if (entity.entity_rank() == mesh::Element || entity.entity_rank() == mesh::Face || entity.entity_rank() == mesh::Edge  )
        {
          int fieldStride = 3;
          {
            unsigned nfr = field->restrictions().size();
            //if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
            for (unsigned ifr = 0; ifr < nfr; ifr++)
              {
                const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                //mesh::Part& frpart = eMesh.getMetaData()->get_part(fr.ordinal());
                fieldStride = fr.stride[0] ;
              }
          }

          out << "Elem: " << entity.identifier() << " rank= " << entity.entity_rank() << " nodes: \n";

          const mesh::PairIterRelation elem_nodes = entity.relations( mesh::Node );
          unsigned num_node = elem_nodes.size();
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity & node = * elem_nodes[ inode ].entity();

              out << " id= " << node.identifier() << " ";
              double *f_data = PerceptMesh::field_data(field, node);
              out << " data = " ;
              for (int ifd=0; ifd < fieldStride; ifd++)
                {
                  out << f_data[ifd] << " ";
                }
              out << std::endl;
            }

        }
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
