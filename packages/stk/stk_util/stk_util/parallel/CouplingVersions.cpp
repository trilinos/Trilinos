#include "stk_util/parallel/CouplingVersions.hpp"
#include "stk_util/parallel/CouplingVersions_impl.hpp"

#include "stk_util/stk_config.h"
#include <stdexcept>
#include <map>
#include <string>
#include <iostream>
#include <cassert>
#include <array>

namespace stk {
namespace util {


#ifdef STK_HAS_MPI

void MPI_Op_MaxMinReduction(void* invec, void* inoutvec, int* len, MPI_Datatype* datatype)
{
  int* invec_int    = reinterpret_cast<int*>(invec);
  int* inoutvec_int = reinterpret_cast<int*>(inoutvec);

  inoutvec_int[1] = std::min(invec_int[1], inoutvec_int[1]);
  inoutvec_int[2] = std::max(invec_int[2], inoutvec_int[2]);
}

std::pair<int, int> allreduce_minmax(MPI_Comm comm, int localVersion)
{
  // for compatibility with the ParallelReduce code, the buffer has to be
  // large enough for 3 ints (empty struct + padding + 2 ints), even though 
  // only the ints are used
  constexpr int bufSize = 3;
  std::array<int, bufSize> inbuf{-1, localVersion, localVersion}, outbuf;

  MPI_Op mpiOp = MPI_OP_NULL ;
  MPI_Op_create( MPI_Op_MaxMinReduction , false , &mpiOp );

  const int root = 0;
  int result = MPI_Reduce(inbuf.data(), outbuf.data(), bufSize*sizeof(int), MPI_BYTE, mpiOp, root, comm);
  if (result != MPI_SUCCESS) {
    std::cerr << "MPI_Reduce failed" << std::endl;
    MPI_Abort(comm, 1);
  }

  result = MPI_Bcast(outbuf.data(), bufSize*sizeof(int), MPI_BYTE, root, comm);
  if (result != MPI_SUCCESS) {
    std::cerr << "MPI_Bcast failed" << std::endl;
    MPI_Abort(comm, 1);
  }

  MPI_Op_free(&mpiOp);

  return {outbuf[1], outbuf[2]};
}


class StkCompatibleVersion
{
  public:

    int get_version() const { return m_version; }

    int get_global_max_version() const { return m_globalMaxVersion; }

    void set_version_for_testing(MPI_Comm comm, unsigned int version)
    {
      set_version_impl(comm, version);
    }

    void set_version(MPI_Comm comm)
    {
      set_version_impl(comm, std::min(impl::SHORT_TERM_STK_MAX_COUPLING_VERSION, m_version) /*m_version*/);
    }

    void set_error_on_reset(bool val)
    {
      m_errorOnResetVersion = val;
    }

    void reset_global_max_coupling_version()
    {
      m_globalMaxVersion = impl::SHORT_TERM_STK_MAX_COUPLING_VERSION; // STK_MAX_COUPLING_VERSION;
    }

  private:

    void set_version_impl(MPI_Comm comm, int maxVersion)
    {
      int oldVersion = m_version;
      auto p = allreduce_minmax(comm, maxVersion);
      int globalMinVersion = p.first;
      m_globalMaxVersion = std::max(m_globalMaxVersion, p.second);
      m_version = globalMinVersion;

      abort_if_incompatible(comm, m_version, maxVersion, oldVersion);
      m_isVersionSet = true;
    }

    bool is_compatible(int requestedVersion)
    {
      return requestedVersion >= STK_MIN_COUPLING_VERSION && requestedVersion <= STK_MAX_COUPLING_VERSION;
    }

    void abort_if_incompatible(MPI_Comm comm, int requestedVersion, int maxVersion, int oldVersion)
    {
      print_unsupported_version_warning(1, __LINE__, __FILE__);

      if (m_isVersionSet && (m_errorOnResetVersion && requestedVersion != oldVersion)) {
        int myRank;
        MPI_Comm_rank(comm, &myRank);
        if (myRank) { 
          std::cerr << "STK version compatibility has already been set, and cannot be changed" << std::endl;
        }

        MPI_Abort(comm, 1);
      }

      if (!is_compatible(requestedVersion)) {
        bool doOutput = m_version >= 2 ? is_output_proc(comm, maxVersion == STK_MAX_COUPLING_VERSION) : true;

        if (doOutput) {
          std::cerr << get_version_error_string(requestedVersion) << std::endl;
        }

        MPI_Abort(comm, 1);
      }
    }

    std::string get_version_error_string(int requestedVersion)
    {
      return std::string("Requested version is not compatible with this version of STK.  Requested version: ") + 
             std::to_string(requestedVersion) + ", minimum supported version " + std::to_string(STK_MIN_COUPLING_VERSION) +
             ", maximum supported version " + std::to_string(STK_MAX_COUPLING_VERSION);
    }


    int is_output_proc(MPI_Comm comm, bool condition)
    {
      print_unsupported_version_warning(1, __LINE__, __FILE__);
      if (m_version < 2) {
        throw std::runtime_error("This function cannot be used with STK versions prior to 2");
      }

      int myRank;
      MPI_Comm_rank(comm, &myRank);
      std::array<int, 2> localData{condition, myRank}, globalData;

      MPI_Allreduce(localData.data(), globalData.data(), 1, MPI_2INT, MPI_MAXLOC, comm);

      return globalData[1] == myRank && condition;
    }

    int m_version = STK_MAX_COUPLING_VERSION;
    int m_globalMaxVersion = impl::SHORT_TERM_STK_MAX_COUPLING_VERSION; // STK_MAX_COUPLING_VERSION;
    bool m_isVersionSet = false;
    bool m_errorOnResetVersion = true;
};



StkCompatibleVersion& get_stk_coupling_version()
{
  static StkCompatibleVersion version;
  return version;
}

#endif

int get_common_coupling_version()
{
#ifdef STK_HAS_MPI
  return get_stk_coupling_version().get_version();
#else
  return STK_impl::SHORT_TERM_MAX_COUPLING_VERSION; //STK_MAX_COUPLING_VERSION;
#endif
}

int get_local_max_coupling_version()
{
  return impl::SHORT_TERM_STK_MAX_COUPLING_VERSION; //STK_MAX_COUPLING_VERSION;
}

int get_local_min_coupling_version()
{
  return STK_MIN_COUPLING_VERSION;
}


int get_global_max_coupling_version()
{
#ifdef STK_HAS_MPI
  return get_stk_coupling_version().get_global_max_version();
#else
  return impl::SHORT_TERM_STK_MAX_COUPLING_VERSION; // STK_MAX_COUPLING_VERSION
#endif
}

std::string get_deprecation_date(int version)
{
  static std::map<int, std::string> deprecationDates{
                                                      std::make_pair(0, "5/22/2022"),
                                                      std::make_pair(1, "7/6/2022"),
                                                      std::make_pair(2, "7/26/2022"),
                                                      std::make_pair(3, "7/26/2022"),
                                                      std::make_pair(4, "7/27/2022"),
                                                      std::make_pair(5, "9/13/2022"),
                                                      std::make_pair(6, "9/18/2022"),
                                                      std::make_pair(7, "10/16/2022"),
                                                      std::make_pair(8, "11/19/2022"),
                                                      std::make_pair(9, "1/31/2023"),
                                                      std::make_pair(10, "4/12/2023"),
                                                      std::make_pair(11, "4/19/2023"),
                                                      std::make_pair(12, "3/11/2024"),
                                                      std::make_pair(13, "3/28/2024"),
                                                      std::make_pair(14, "")
                                                    };

  return deprecationDates.at(version);
}


void set_coupling_version(MPI_Comm comm)
{
#ifdef STK_HAS_MPI
  get_stk_coupling_version().set_version(comm);
#endif
}

bool is_local_stk_coupling_deprecated()
{
  return STK_MAX_COUPLING_VERSION < get_global_max_coupling_version();
}


void print_unsupported_version_warning(int version, int line, const char* file)
{                                                                                      
  if ( STK_MIN_COUPLING_VERSION > version ) {
    std::cerr  << "The function at line " << line << " of file " << file
               << " can be simplified now that STK_MIN_COUPLING_VERSION is greater than "
               << (version) << std::endl;
  }
}


namespace impl {

void set_coupling_version(MPI_Comm comm, int version)
{
#ifdef STK_HAS_MPI
  get_stk_coupling_version().set_version_for_testing(comm, static_cast<unsigned int>(version));
#endif
}

void set_error_on_reset(bool val)
{
#ifdef STK_HAS_MPI
  get_stk_coupling_version().set_error_on_reset(val);
#endif
}

void reset_global_max_coupling_version()
{
  get_stk_coupling_version().reset_global_max_coupling_version();
}

}

}
}

