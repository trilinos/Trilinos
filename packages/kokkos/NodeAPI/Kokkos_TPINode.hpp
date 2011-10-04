// @HEADER
// ***********************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_TPINODE_HPP_
#define KOKKOS_TPINODE_HPP_

#include "Kokkos_StandardNodeMemoryModel.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include <TPI.h>

namespace Teuchos {
  // forward declarations
  class ParameterList;
}

namespace Kokkos {

  template <class WDP>
  struct WDPPlusRange {
    WDPPlusRange(int Beg, int End, WDP Wdp) : wdp(Wdp), beg(Beg), end(End){}
    WDP wdp;
    int beg, end;
  };

  inline
  void tpi_work_span(TPI_Work* work, int beg, int end, int& ibeg, int& iend)
  {
    const int chunk = ( end - beg + work->count - 1 ) / work->count ;
  
    iend = chunk * ( work->rank + 1 ) + beg;
    ibeg = chunk * ( work->rank ) + beg;
  
    if ( end < iend ) { iend = end; }
  }

  template<class WDP>
  void tpi_execute(TPI_Work * work)
  {
    // get work/data pair
    const WDPPlusRange<WDP>* const_wdp_wrapper = static_cast<const WDPPlusRange<WDP>*>(work->info);
    WDPPlusRange<WDP>* wdp_wrapper = const_cast<WDPPlusRange<WDP>*>(const_wdp_wrapper);
    WDP wdp = wdp_wrapper->wdp;
    int beg = wdp_wrapper->beg, end = wdp_wrapper->end;
    int ibeg, iend;
    // determine my share of the work
    tpi_work_span(work, beg, end, ibeg, iend);
    // do my share of the work
    for (int i=ibeg; i<iend; ++i) {
      wdp.execute(i);
    }
  }

  template<class WDP>
  void tpi_reduction_work(TPI_Work * work)
  {
    const WDPPlusRange<WDP>* wdp_wrapper = static_cast<const WDPPlusRange<WDP>*>(work->info);
    int beg = wdp_wrapper->beg, end = wdp_wrapper->end;
    WDP wdp = wdp_wrapper->wdp;
    int ibeg, iend;
    tpi_work_span(work, beg, end, ibeg, iend);
  
    typedef typename WDP::ReductionType ReductionType;
    ReductionType tmpi;
    ReductionType &res = *(static_cast<ReductionType*>(work->reduce));
  
    for (int i=ibeg; i<iend; ++i) {
      tmpi = wdp.generate(i);
      res = wdp.reduce(res, tmpi);
    }
  }

  template<class WDP>
  void tpi_reduction_join(TPI_Work * work, const void* src)
  {
    typedef typename WDP::ReductionType ReductionType;
  
    const WDPPlusRange<WDP>* wdp_wrapper = static_cast<const WDPPlusRange<WDP>*>(work->info);
    WDP wdp = wdp_wrapper->wdp;
  
    ReductionType& work_reduce = *(static_cast<ReductionType*>(work->reduce));
    const ReductionType& src_reduce  = *(static_cast<const ReductionType*>(src));
  
    work_reduce = wdp.reduce(work_reduce, src_reduce);
  }

  template<class WDP>
  void tpi_reduction_init(TPI_Work * work)
  {
    typedef typename WDP::ReductionType ReductionType;
    *(static_cast<ReductionType*>(work->reduce)) = WDP::identity();
  }

  /** \brief %Kokkos node interface to the ThreadPool threading library.
      \ingroup kokkos_node_api
   */
  class TPINode : public StandardNodeMemoryModel {
    public:

      /*! \brief Constructor acceptings a list of parameters
          
          This constructor accepts the parameters:
          \param "Num Threads" [int] Specifies the number of threads, calls TPINode::init(). Throws std::runtime_error if less than zero. Default: 0.
          \param "Verbose"     [int] Non-zero parameter specifies that the constructor is verbose, printing information about the number of threads. Default: 0.
          
       */
      TPINode(Teuchos::ParameterList &plist);

      /*! \brief Thread initialization method.
          If \c numThreads is greater than zero, this calls TPI_Init(). If the threads have already been initialized by this node, it first calls TPI_Finalize().
       */
      void init(int numThreads);

      /*! \brief Default destructor calls TPI_Finalize().
          TPI_Finalize() is called if the number of initialized threads is greater than zero; otherwise, the destructor has no effect.
      */
      ~TPINode();

      //! \begin parallel for skeleton, a wrapper around TPI_Run_threads. See \ref kokkos_node_api "Kokkos Node API"
      template <class WDP>
      static void parallel_for(int beg, int end, WDP wd) {
        WDPPlusRange<WDP> wdp_plus(beg,end,wd);
        TPI_Run_threads(tpi_execute<WDP>, &wdp_plus, 0 );
      }

      //! \begin parallel reduction skeleton, a wrapper around TPI_Run_threads_reduce. See \ref kokkos_node_api "Kokkos Node API"
      template <class WDP>
      static typename WDP::ReductionType 
      parallel_reduce(int beg, int end, WDP wd) {
        typedef typename WDP::ReductionType ReductionType;
        ReductionType result = WDP::identity();
        WDPPlusRange<WDP> wdp_plus(beg,end,wd);
        TPI_Run_threads_reduce(tpi_reduction_work<WDP>, &wdp_plus,
                               tpi_reduction_join<WDP>,
                               tpi_reduction_init<WDP>, sizeof(result), &result);
        return result;
      }

      //! \begin No-op for TPINode.
      inline void sync() const {};

    private:
      int curNumThreads_;
  };

  template <> class ArrayOfViewsHelper<TPINode> : public ArrayOfViewsHelperTrivialImpl<TPINode> {};

} // end namespace Kokkos

#endif
