//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef KOKKOS_TBBNODE_HPP_
#define KOKKOS_TBBNODE_HPP_

#include "Kokkos_StandardNodeMemoryModel.hpp"
#include "Kokkos_NodeHelpers.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>

namespace Teuchos {
  // forward declarations
  class ParameterList;
}

#include <stdlib.h>

namespace KokkosClassic {

  template <class WDPin>
  struct BlockedRangeWDP {
    mutable WDPin *wd;
    inline BlockedRangeWDP (WDPin &in) : wd (&in) {}
    inline void operator () (tbb::blocked_range<int>& rng) const {
      for (int i = rng.begin (); i != rng.end (); ++i) {
        wd->execute (i);
      }
    }
  };

  template <class WDPin>
  struct BlockedRangeWDPReducer {
    WDPin &wd;
    typename WDPin::ReductionType result;
    BlockedRangeWDPReducer (WDPin &in) :
      wd (in),
      result (in.identity ())
    {}

    BlockedRangeWDPReducer (BlockedRangeWDPReducer &in, tbb::split) :
      wd (in.wd)
    {
      result = wd.identity ();
    }

    void operator () (tbb::blocked_range<int> &rng) {
      typename WDPin::ReductionType tmpi;
      const int end = rng.end ();
      for (int i = rng.begin (); i != end; ++i) {
        tmpi = wd.generate (i);
        result = wd.reduce (result, tmpi);
      }
    }
    inline void join (const BlockedRangeWDPReducer<WDPin> &other) {
      result = wd.reduce (result, other.result);
    }
  };

  /** \brief %Kokkos node interface to the Intel Threading Building Blocks threading library.
      \ingroup kokkos_node_api
   */
  class TBBNode : public StandardNodeMemoryModel {
  public:
    //! Constructor that sets default parameters.
    TBBNode ();

    /*! \brief Constructor that takes a list of parameters.

      This constructor accepts the following parameters:
      - "Num Threads" [int] Specifies the number of threads; calls
        TBBNode::init() if non-negative. Otherwise, late
        initialization. Default: -1.
    */
    TBBNode (Teuchos::ParameterList &pl);

    //! Default destructor, calls tbb::task_scheduler_init::terminate().
    ~TBBNode ();

    /*! \brief Get default parameters for this node */
    static Teuchos::ParameterList getDefaultParameters ();

    /// \brief Initialize TBB with the given number of threads.
    ///
    /// This method calls TBB's tbb::task_scheduler_init::initialize()
    /// method, with \c numThreads as the argument if it is greater
    /// than 0, and tbb::task_scheduler_init::automatic otherwise.  If
    /// init has already been called, this calls
    /// tbb:task_scheduler_init::terminate() first.
    void init (const int numThreads);

    /// \brief Parallel for "skeleton"; a wrapper around tbb::parallel_for.
    ///
    /// See \ref kokkos_node_api "Kokkos Node API"
    template <class WDP>
    static void parallel_for (const int begin, const int end, WDP wd) {
      BlockedRangeWDP<WDP> tbb_wd(wd);
      tbb::parallel_for (tbb::blocked_range<int> (begin, end), tbb_wd,
                         tbb::auto_partitioner ());
    }

    /// \brief Parallel reduction "skeleton"; a wrapper around tbb::parallel_reduce.
    ///
    /// See \ref kokkos_node_api "Kokkos Node API"
    template <class WDP>
    static typename WDP::ReductionType
    parallel_reduce (const int begin, const int end, WDP wd) {
      BlockedRangeWDPReducer<WDP> tbb_wd (wd);
      tbb::parallel_reduce (tbb::blocked_range<int> (begin, end), tbb_wd,
                            tbb::auto_partitioner ());
      return tbb_wd.result;
    }

    //! No-op for TBBNode.
    inline void sync() const {}

    /// \brief Return the name of the node type
    /// See \ref kokkos_node_api "Kokkos Node API"
    static std::string name ();

  private:
    tbb::task_scheduler_init tsi_;
    bool alreadyInit_;
  };

#ifdef _MSC_VER
#pragma warning(push)
// destructor could not be generated because a base class destructor is inaccessible
#pragma warning(disable : 4624)
#endif

  template <>
  class ArrayOfViewsHelper<TBBNode> :
    public ArrayOfViewsHelperTrivialImpl<TBBNode>
  {};

#ifdef _MSC_VER
#pragma warning(pop)
#endif

} // namespace KokkosClassic

#endif
