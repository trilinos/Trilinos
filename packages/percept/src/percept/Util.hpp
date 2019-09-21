// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Util_hpp
#define percept_Util_hpp

#include <iostream>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <math.h>


#include "Shards_Array.hpp"
#include "Shards_CellTopology.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"
// #include "Shards_Array.hpp"
#include "Teuchos_RCP.hpp"
// #include "Teuchos_BLAS.hpp"
#include "Teuchos_oblackholestream.hpp"
//#include "Teuchos_Assert.hpp"

#include <percept/ExceptionWatch.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/tokenize.hpp>

#include <boost/lexical_cast.hpp>
#include <percept/stk_mesh.hpp>


// output operators
// bcarnes: needed to put these in the namespace matching the second
//   arg to operator << in order to get compiler (gcc 4.4.4) to work

namespace std {

    template<class T, class L, class A>
    std::ostream& operator<<(std::ostream& out,  const std::set<T,L,A>& val)
    {
      copy(val.begin(), val.end(), std::ostream_iterator<T>(out, " "));
      return out;
    }

    template<class T, class L>
    std::ostream& operator<<(std::ostream& out,  const std::vector<T,L>& val)
    {
      copy(val.begin(), val.end(), std::ostream_iterator<T>(out, " "));
      return out;
    }

    template<class Key, class Val, class Comp, class Alloc>
    std::ostream& operator<<(std::ostream& out, std::map<Key, Val, Comp, Alloc>& val)
    {
      typename std::map<Key, Val, Comp, Alloc >::iterator it;
      for (it = val.begin();
           it != val.end();
           it++)
        {
          out << "map[ " << (*it).first << "]= " << (*it).second << " \n";
        }

      return out;
    }

    template<typename Scalar>
    std::ostream& operator<<(std::ostream& os, const std::vector<Scalar>& vec)
    {
      os << "vector size()= " << vec.size() << " entries= [ \n";
      for (unsigned i = 0; i < vec.size(); i++)
        os << " " << vec[i] << " , "; // << "\n";
      os << " ]\n";
      return os;
    }

}

namespace shards {
  std::ostream& operator<<(std::ostream& os, const shards::Array<double, shards::NaturalOrder>& container) ;
}

  namespace percept {

      enum TraceTypes {
        CONNECT_LOCAL,
        CONNECT_LOCAL_createNewNeededNodes,
        CONNECT_LOCAL_createNewElements,
        CONNECT_LOCAL_URP_createOrGetNode,
        CONNECT_LOCAL_URP_declare_relation
      };

    //========================================================================================================================
    enum TurboOption {
      TURBO_ELEMENT,
      TURBO_BUCKET,
      TURBO_NONE
    };

    class PerceptMesh;

    //========================================================================================================================
    class Util
    {
    public:
      static double s_timers[10];

      static bool isLinearElement(shards::CellTopology& cell_topo);
      static double cpu_time();
      static size_t memory(size_t& heap);
      static void pause(bool doPause=true, const char *msg=0);
      static void setDoPause(bool enableOrDisablePause);
      static void setInfo(int val);
      static int  getInfo();
      static void setFlag(int which_flag, bool val);
      static bool getFlag(int which_flag);
      static void setRank(int rank);
      static int get_rank();
      static void debug_stop();
      static void replace(std::string &str, const std::string &find_what, const std::string &replace_with);

      template<class T>
      static void makeUnique(std::vector<T>& vec)
      {
        std::sort(vec.begin(), vec.end());
        typename std::vector<T>::iterator del = std::unique(vec.begin(), vec.end());
        vec.resize(std::distance(vec.begin(), del));
      }

      // return true if input starts with the fragment
      static bool startsWith(const std::string& input, const std::string& fragment) { return input.substr(0, fragment.length()) == fragment; }
      static bool startsWith(const char* input, const char* fragment) 
      {
        for (int i = 0; fragment[i] != '\0'; ++i)
        {
          if (input[i] != fragment[i]) // e.g. input[i] == '\0'
            return false;
        }
        return true; 
      }
      //static std::string convert_to_mm(double d, int precision=15, std::string mm_prec="`15");
      static std::string convert_to_mm(double d, int precision=6, std::string mm_prec="");
      static std::string convertToMathematica(double d, int precision=6, std::string mm_prec="") {
        return convert_to_mm(d,precision,mm_prec);
      }
      static std::string convert_to_mathematica(double d, int precision=6, std::string mm_prec="") {
        return convert_to_mm(d,precision,mm_prec);
      }

      static void trace_cpu_time_and_mem_0(unsigned index);
      static void trace_cpu_time_and_mem_1(unsigned index);
      static void trace_cpu_time_and_mem_print(int index, std::string msg);

      static std::vector<std::string> split(std::string str, std::string separators, bool remove_spaces=true)
      {
        if (remove_spaces)
          {
            replace(str, " ", "");
          }
        std::vector<std::string> tokens;
        stk::util::tokenize(str, separators, tokens);
        return tokens;
      }
      static std::string join(std::vector<std::string> tokens, std::string separator, bool add_spaces=true)
      {
        if (add_spaces)
          {
            separator += " ";
          }
        std::string ret;
        for (unsigned i = 0; i < tokens.size(); i++)
          {
            ret += tokens[i];
            if (i != tokens.size() - 1)
              ret += separator;
          }
        return ret;
      }

      static bool approx_equal_relative(double a, double b, double tol)
      {
        if (fabs(a-b) <= tol*(fabs(a)+fabs(b))*0.5)
          return true;
        else
          return false;
      }
      static bool approx_equal_absolute(double a, double b, double tol)
      {
        if (fabs(a-b) < tol)
          return true;
        else
          return false;
      }

      static bool file_exists(std::string filename)
      {
        std::fstream fin;
        fin.open(filename.c_str(), std::ios::in);
        if( fin.is_open() )
          {
            fin.close();
            return true;
          }
        fin.close();
        return false;
      }

    };

    //========================================================================================================================
    template<class VecType>
    class GenericVectorOfObjectPointers : public std::vector<VecType *>
    {
    public:
      GenericVectorOfObjectPointers(int n) : std::vector<VecType *>(n) {}
      GenericVectorOfObjectPointers(VecType *vt1=0,
                    VecType *vt2=0,
                    VecType *vt3=0,
                    VecType *vt4=0,
                    VecType *vt5=0,
                    VecType *vt6=0,
                    VecType *vt7=0,
                    VecType *vt8=0)
      {
        if (vt1) push_back(vt1);
        if (vt2) push_back(vt2);
        if (vt3) push_back(vt3);
        if (vt4) push_back(vt4);
        if (vt5) push_back(vt5);
        if (vt6) push_back(vt6);
        if (vt7) push_back(vt7);
        if (vt8) push_back(vt8);
      }
    };


    //========================================================================================================================
    enum {
      PERCEPT_TIMER_ROOT = 0x00001000		///< Enable root timer
    };
    stk::diag::TimerSet& perceptTimerSet();
    stk::diag::Timer& perceptTimer() ;
    typedef stk::diag::MetricTraits<stk::diag::WallTime>::Type LapTimeType;
    LapTimeType getLapTime(stk::diag::Timer& lap_timer);
    typedef stk::diag::MetricTraits<stk::diag::LapCount>::Type LapCountType;
    LapCountType getAccumulatedLap(stk::diag::Timer& timer, bool option=false);

    //========================================================================================================================

    template<class T>
    std::string toString(T t) { return boost::lexical_cast<std::string>(t); }

    inline std::string toString(stk::topology::rank_t t) { return toString<unsigned int>(static_cast<unsigned int>(t)); }

    template<class T>
    inline
    T square(T t) { return t*t; }

    template<class T>
    inline
    T SQR(T t) { return t*t; }

    inline int toInt(std::string t) { return boost::lexical_cast<int>(t); }

    //========================================================================================================================

#define QUOTE(A) #A
#define EXPAND_AND_QUOTE(A) QUOTE(A)

#define PERCEPT_OUT(a) " " << QUOTE(a) << " = " << a

#define TOKENPASTE2(x,y) x ## y
#define TOKENPASTE(x,y) TOKENPASTE2(x,y)

#define VERIFY_OP_ON( val1, op, val2, message) do { if (1) { if (!(val1 op val2)) { std::ostringstream msg_loc; \
                                                               msg_loc << "\n\n  ERROR: " << __FILE__ << " line: " << __LINE__; \
                                                               msg_loc << "\n          " << message << "\n" \
                                                                       << "   { " #val1 << " = " << val1 << " } \n" \
                                                                       << "   { " #val2 << " = " << val2 << " } \n"; \
                                                               msg_loc << " ERROR condition is: \n    " << #val1 " " #op " " #val2 << "\n"; \
                                                               std::cout << msg_loc.str() << std::endl; \
                                                               percept::Util::debug_stop(); \
                                                               throw std::runtime_error(msg_loc.str()); } } } while (0)

#define VERIFY_OP_ON_BOOL( val1, op, val2, message, result) do { if (1) { if (!(val1 op val2)) { result = false; std::ostringstream msg_loc; \
                                                               msg_loc << "\n\n  ERROR: " << __FILE__ << " line: " << __LINE__; \
                                                               msg_loc << "\n          " << message << "\n" \
                                                                       << "   { " #val1 << " = " << val1 << " } \n" \
                                                                       << "   { " #val2 << " = " << val2 << " } \n"; \
                                                               msg_loc << " ERROR condition is: \n    " << #val1 " " #op " " #val2 << "\n"; \
                                                               std::cout << msg_loc.str() << std::endl; \
                                                               percept::Util::debug_stop(); \
                                                                } } } while (0)
#define VERIFY_OP_ON_BOOL_NO_FAIL( val1, op, val2, message, result) do { if (1) { if (!(val1 op val2)) { result = false; } } } while(0)

#define VERIFY_MSG( message) do { if (1) { { std::ostringstream msg_loc; \
                                                               msg_loc << "\n\n  ERROR: " << __FILE__ << " line: " << __LINE__; \
                                                               msg_loc << "\n          " << message << "\n"; \
                                                               std::cout << msg_loc.str() << std::endl; \
                                                               percept::Util::debug_stop(); \
                                                               throw std::runtime_error(msg_loc.str()); } } } while (0)


#if defined(NDEBUG) || defined(__CUDA_ARCH__)

#define VERIFY_1(message)  do {} while(0)
#define VERIFY(expr, val1, val2, message)  do {} while(0)
#define VERIFY_OP( val1, op, val2, message)  do {} while(0)
#define VERIFY_EQ1( val1, val2, message)  do {} while(0)
#define VERIFY_NE( val1, val2, message)  do {} while(0)

#else

#define VERIFY_ON 1

#define VERIFY_1(message) do { if(VERIFY_ON) {  std::ostringstream msg_loc; msg_loc << message; \
                                              std::cout << msg_loc.str() << std::endl;          \
                                              throw std::runtime_error(msg_loc.str());  } } while (0)


#define VERIFY(expr, val1, val2, message) do { if (VERIFY_ON) { if (!(expr)) { std::ostringstream msg_loc; msg_loc << message << " " << val1 << " " << val2; \
                                              std::cout << msg_loc.str() << std::endl; \
                                              percept::Util::debug_stop();       \
                                              throw std::runtime_error(msg_loc.str()); } } } while (0)

#define VERIFY_OP( val1, op, val2, message) do { if (VERIFY_ON) { if (!(val1 op val2)) { std::ostringstream msg_loc; \
                                                msg_loc << "\n\n  ERROR: " << __FILE__ << " line: " << __LINE__; \
                                                msg_loc << "\n          " << message << "\n" \
                                                        << "   { " #val1 << " = " << val1 << " } \n" \
                                                        << "   { " #val2 << " = " << val2 << " } \n"; \
                                                msg_loc << " ERROR condition is: \n    " << #val1 " " #op " " #val2 << "\n"; \
                                                std::cout << msg_loc.str() << std::endl; \
                                                percept::Util::debug_stop();     \
                                                throw std::runtime_error(msg_loc.str()); } } } while (0)



#define VERIFY_EQ1( val1, val2, message) do { if (VERIFY_ON) { if (!(val1==val2)) { std::ostringstream msg_loc; msg_loc << message << " " << val1 << " " << val2; \
                                             std::cout << msg_loc.str() << std::endl; \
                                             percept::Util::debug_stop();        \
                                             throw std::runtime_error(msg_loc.str()); } } } while (0)

#define VERIFY_NE( val1, val2, message) do { if (VERIFY_ON) { if (!(val1 != val2)) { std::ostringstream msg_loc; msg_loc << message << " " << val1 << " " << val2; \
                                            std::cout << msg_loc.str() << std::endl; \
                                            percept::Util::debug_stop();         \
                                            throw std::runtime_error(msg_loc.str()); } } } while (0)

#endif

#define VERIFY_EQ( val1, val2, message) VERIFY_OP(val1, ==, val2, message)
#define VERIFY_EQ_ON( val1, val2, message) VERIFY_OP_ON(val1, ==, val2, message)

#define VERIFY_TRUE( val1) VERIFY_OP( (val1), ==, true, EXPAND_AND_QUOTE(val1))
#define VERIFY_TRUE_ON( val1) VERIFY_OP_ON( (val1), ==, true, EXPAND_AND_QUOTE(val1))


  }

#endif
