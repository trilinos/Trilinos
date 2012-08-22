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


#ifndef KOKKOS_OPENMPNODE_HPP_
#define KOKKOS_OPENMPNODE_HPP_

#include "Kokkos_StandardNodeMemoryModel.hpp"
#include "Kokkos_NodeHelpers.hpp"

#include <omp.h>

namespace Teuchos {
  //forward declarations
  class ParameterList;
}

namespace Kokkos {

  /// \class OpenMPNode
  /// \brief Kokkos Node using OpenMP for parallelization.
  class OpenMPNode : public StandardNodeMemoryModel {
  public:

    /// \brief Constructor.
    ///
    /// The constructor currently accepts the following (optional)
    /// parameters:
    /// - "Num Threads" (int): The number of threads that OpenMP
    ///   should use.  If not provided, or if -1 or 0, OpenMP will
    ///   pick the number of threads in the usual way.
    /// - "Verbose" (int): If nonzero, the Kokkos Node will print
    ///    status output to std::cout.
    ///
    /// \warning If you set the "Num Threads" parameter to a positive
    ///   value, this will set the number of threads that _all_
    ///   OpenMPNode instances will use, not just this one.
    OpenMPNode(Teuchos::ParameterList &pl);

    /// \brief Default constructor (sets default parameter values).
    ///
    /// For parameters and their default values, see the documentation
    /// of the constructor that accepts a ParameterList input / output
    /// argument.
    OpenMPNode();

    //! Destructor.
    ~OpenMPNode();

    /// \brief Set the number of threads that OpenMP should use.
    ///
    /// It's not necessary to call this method unless you want to
    /// change the number of threads that OpenMP should use, after
    /// constructing the Kokkos Node instance.
    ///
    /// \warning This method will affect the number of threads used by
    ///   all OpenMPNode instances.
    ///
    /// \param numThreads [in] The number of threads that OpenMP
    ///   should use.  Ignored if -1 or 0.
    void init(int numThreads);

    //! Perform a parallel for loop on the given half-exclusive index range.
    template <class WDP>
    static void parallel_for(int beg, int end, WDP wd) {
#pragma omp parallel for schedule(guided) default(shared)
      for (int i = beg; i < end; ++i) {
        wd.execute(i);
      }
    }

    //! Perform a parallel reduction on the given half-exclusive index range.
    template <class WDP>
    static typename WDP::ReductionType
    parallel_reduce(int beg, int end, WDP wd) {
      typedef typename WDP::ReductionType ReductionType;
      ReductionType threadResult = WDP::identity();
      ReductionType globalResult = WDP::identity();
#pragma omp parallel default(shared) firstprivate(threadResult)
      {
#pragma omp for
        for (int i = beg; i < end; ++i) {
          threadResult = wd.reduce(threadResult, wd.generate(i));
        }
#pragma omp critical
        globalResult = wd.reduce(globalResult, threadResult);
      }
      return globalResult;
    }

    //! Synchronize threads; this is a no-op for OpenMPNode.
    inline void sync() const {};

  private:
    /// \brief "Num Threads" parameter value.
    ///
    /// If <= 0, OpenMPNode ignores this and lets OpenMP pick its own
    /// number of threads.
    int curNumThreads_;

    //! Whether to print verbose status output to std::cout.
    bool verbose_;
  };

  template <> class ArrayOfViewsHelper<OpenMPNode> : public ArrayOfViewsHelperTrivialImpl<OpenMPNode> {};

}

#endif
