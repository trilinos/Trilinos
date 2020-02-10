/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSSPARSE_IMPL_SPILUK_SYMBOLIC_HPP_
#define KOKKOSSPARSE_IMPL_SPILUK_SYMBOLIC_HPP_

/// \file KokkosSparse_spiluk_symbolic_impl.hpp
/// \brief Implementation of the symbolic phase of sparse ILU(k).

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_spiluk_handle.hpp>

//#define SYMBOLIC_OUTPUT_INFO

namespace KokkosSparse {
namespace Impl {
namespace Experimental {


template <class IlukHandle,
          class RowMapType, 
          class EntriesType,
          class LevelType1,
          class LevelType2,
          class size_type>
void level_sched ( IlukHandle& thandle,
                   const RowMapType row_map, const EntriesType entries, const size_type nrows, 
                   LevelType1& level_list, LevelType2& level_ptr, LevelType2& level_idx, size_type &nlevels ) {
  // Scheduling currently compute on host

  typedef typename IlukHandle::nnz_lno_t nnz_lno_t;

  nlevels      = 0;
  level_ptr(0) = 0;

  for ( size_type i = 0; i < nrows; ++i ) {
    size_type l = 0;
    size_type rowstart= row_map(i);
    size_type rowend  = row_map(i+1);
    for ( size_type j = rowstart; j < rowend; ++j ) {
      nnz_lno_t col = entries(j);
      l = std::max(l, level_list(col));
    }
    level_list(i)   = l+1;
    level_ptr(l+1) += 1;
    nlevels         = std::max(nlevels, l+1);
  }

  for ( size_type i = 1; i <= nlevels; ++i ) {
    level_ptr(i) += level_ptr(i-1);
  }

  for ( size_type i = 0; i < nrows; i++ ) {
    level_idx(level_ptr(level_list(i)-1)) = i;
    level_ptr(level_list(i)-1) += 1;
  }

  for ( size_type i = nlevels-1; i > 0; --i ) {
    level_ptr(i) = level_ptr(i-1);
  }

  level_ptr(0) = 0;

  //Find the maximum number of rows of levels
  size_type maxrows = 0;
  for ( size_type i = 0; i < nlevels; ++i ) {
    size_type lnrows = level_ptr(i+1) - level_ptr(i);
    if( maxrows < lnrows ) {
      maxrows = lnrows;
    }
  }

  thandle.set_num_levels(nlevels);
  thandle.set_level_maxrows(maxrows);
 
}

 //Linear Search for the smallest row index
template <class size_type, class nnz_lno_t, class ViewType>
size_type search_col_index ( nnz_lno_t j, size_type lenl, ViewType h_iL, ViewType h_llev, ViewType h_iw ) {

  nnz_lno_t irow = h_iL(j);
  nnz_lno_t ipos = j;

  //Find the smallest col index
  for(size_type k = j+1; k < lenl; ++k) {
    if( h_iL(k) < irow ) {
      irow = h_iL(k);
      ipos = k;
    }
  }

  if( ipos != j ) {//Swap entries
    nnz_lno_t row = h_iL(j);
    h_iL(j)       = h_iL(ipos);
    h_iL(ipos)    = row;

    nnz_lno_t t   = h_llev(j);
    h_llev(j)     = h_llev(ipos);
    h_llev(ipos)  = t;

    h_iw(irow)    = j;
    h_iw(row)     = ipos;
  }
  return ((size_type)irow);
}

template <class IlukHandle,
          class ARowMapType,
          class AEntriesType,
          class LRowMapType,
          class LEntriesType,
          class URowMapType,
          class UEntriesType>
void iluk_symbolic ( IlukHandle& thandle,
                     const typename IlukHandle::const_nnz_lno_t &fill_lev,
                     const ARowMapType&  A_row_map_d,
                     const AEntriesType& A_entries_d,
                           LRowMapType&  L_row_map_d,
                           LEntriesType& L_entries_d,
                           URowMapType&  U_row_map_d,
                           UEntriesType& U_entries_d ) {

 if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_RP ||
      thandle.get_algorithm() == KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1 )
/*   || thandle.get_algorithm() == KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHED_TP2 )*/
 {
  // Scheduling and symbolic phase currently compute on host - need host copy of all views

  typedef typename ARowMapType::HostMirror  AHostRowMapType;
  typedef typename AEntriesType::HostMirror AHostEntriesType;
  typedef typename LRowMapType::HostMirror  LHostRowMapType;
  typedef typename LEntriesType::HostMirror LHostEntriesType;
  typedef typename URowMapType::HostMirror  UHostRowMapType;
  typedef typename UEntriesType::HostMirror UHostEntriesType;

  typedef typename IlukHandle::size_type size_type;
  typedef typename IlukHandle::nnz_lno_t nnz_lno_t;

  typedef typename IlukHandle::nnz_lno_view_t             HandleDeviceEntriesType;
  typedef typename IlukHandle::nnz_lno_view_t::HostMirror HandleHostEntriesType;

  typedef typename IlukHandle::nnz_row_view_t             HandleDeviceRowMapType;
  typedef typename IlukHandle::nnz_row_view_t::HostMirror HandleHostRowMapType;

  //typedef typename IlukHandle::signed_integral_t signed_integral_t;

  size_type nrows = thandle.get_nrows();

  AHostRowMapType A_row_map = Kokkos::create_mirror_view(A_row_map_d);
  Kokkos::deep_copy(A_row_map, A_row_map_d);

  AHostEntriesType A_entries = Kokkos::create_mirror_view(A_entries_d);
  Kokkos::deep_copy(A_entries, A_entries_d);

  LHostRowMapType  L_row_map = Kokkos::create_mirror_view(L_row_map_d);
  LHostEntriesType L_entries = Kokkos::create_mirror_view(L_entries_d);
  UHostRowMapType  U_row_map = Kokkos::create_mirror_view(U_row_map_d);
  UHostEntriesType U_entries = Kokkos::create_mirror_view(U_entries_d);
  
  HandleDeviceRowMapType dlevel_list = thandle.get_level_list();
  HandleHostRowMapType level_list    = Kokkos::create_mirror_view(dlevel_list);
  Kokkos::deep_copy(level_list, dlevel_list);
  
  HandleDeviceEntriesType dlevel_ptr = thandle.get_level_ptr();
  HandleHostEntriesType level_ptr    = Kokkos::create_mirror_view(dlevel_ptr);
  Kokkos::deep_copy(level_ptr, dlevel_ptr);

  HandleDeviceEntriesType dlevel_idx = thandle.get_level_idx();
  HandleHostEntriesType level_idx    = Kokkos::create_mirror_view(dlevel_idx);
  Kokkos::deep_copy(level_idx, dlevel_idx);
 
  size_type nlev = 0;

  //Level scheduling on A???
  //level_sched<IlukHandle, AHostRowMapType, AHostEntriesType, HandleHostRowMapType, HandleHostEntriesType, size_type > 
  //                                    (thandle, A_row_map, A_entries, nrows, level_list, level_ptr, level_idx, nlev);
  //level_sched (thandle, A_row_map, A_entries, nrows, level_list, level_ptr, level_idx, nlev);

  //Symbolic phase
  //Kokkos::resize(L_row_map_d, nrows-3);// error: static assertion failed: Can only resize managed views
  //Kokkos::resize(L_entries_d, L_entries_d.extent(0)-3);
  //thandle.set_nnzL(L_entries_d.extent(0)+5);

  typedef Kokkos::View<nnz_lno_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> HostTmpViewType;
    
  HostTmpViewType h_lev ( "h_lev",  thandle.get_nnzL() );
  HostTmpViewType h_iw  ( "h_iw",   nrows              );
  HostTmpViewType h_iL  ( "h_iL",   nrows              );
  HostTmpViewType h_llev( "h_llev", nrows              );

  size_type cntL = 0;
  size_type cntU = 0;
  size_type iU, ulev, lenu, lenl;
  
  L_row_map(0) = 0;
  U_row_map(0) = 0;

  Kokkos::deep_copy( h_iw, nnz_lno_t(-1) );

  //Main loop
  for (size_type i = 0; i < nrows; ++i) {
    iU = i;
    ulev = i;
    lenl = lenu = 0;

    //Unpack the ith row
    size_type k1 = A_row_map(i);
    size_type k2 = A_row_map(i+1);

    for ( size_type k = k1; k < k2; ++k ) {
      size_type col = static_cast<size_type>(A_entries(k));
      if (col > i) {//U part
        h_iw(col)         = lenu;
        h_iL(iU+lenu)     = col;
        h_llev(ulev+lenu) = 0;
        lenu++;
      } 
      else if (col < i) {//L part
        h_iw(col)    = lenl;
        h_iL(lenl)   = col;
        h_llev(lenl) = 0;
        lenl++;
      }
    }

    //Eliminate rows
    nnz_lno_t j = -1;
    while (static_cast<size_type>(++j) < lenl) {
      size_type row  = search_col_index(j, lenl, h_iL, h_llev, h_iw);
      nnz_lno_t jlev = h_llev(j);
      k1 = U_row_map(row)  +1;
      k2 = U_row_map(row+1);
      for (size_type k = k1; k < k2; ++k) {
        size_type col  = static_cast<size_type>(U_entries(k));
        nnz_lno_t lev1 = jlev + h_lev(k) + 1;
        if (lev1 > fill_lev) continue;
        nnz_lno_t ipos = h_iw(col);
        if (ipos == -1) {//Fill-in
          if (col > i) {//U part
            h_iw(col)         = lenu;
            h_iL(iU+lenu)     = col;
            h_llev(ulev+lenu) = lev1;
            lenu++;
          }
          else if (col < i) {//L part
            h_iw(col)    = lenl;
            h_iL(lenl)   = col;
            h_llev(lenl) = lev1;
            lenl++;
          }
        }
        else {//Not a fill-in
          if (col > i) 
            h_llev(ulev+ipos) = std::min(h_llev(ulev+ipos), lev1);
          else if (col < i)
            h_llev(ipos) = std::min(h_llev(ipos), lev1);
        }
      }
    }

    //Reset iw
    for (size_type k = 0; k < lenl; ++k) h_iw(h_iL(k))    = -1;
    for (size_type k = 0; k < lenu; ++k) h_iw(h_iL(iU+k)) = -1;

    //Copy U part+diag and levels
    //if (cntU+lenu+1 > U_entries_d.extent(0)) {
    //  size_type newsize = (size_type)(U_entries_d.extent(0)*EXPAND_FACT);
    //  //thandle.set_nnzU(newsize);
    //  Kokkos::resize(h_lev, newsize);
    //}
    //U diag entry
    U_entries(cntU) = i;
    cntU++;
    //U part
    for (size_type k = 0; k < lenu; ++k) {
      U_entries(cntU) = h_iL(iU+k);
      h_lev(cntU)     = h_llev(ulev+k);
      cntU++;
    }
    U_row_map(i+1) = cntU;

    //Copy L part
    //if (cntL+lenl > L_entries_d.extent(0)) {
    //  size_type newsize = (size_type) (L_entries_d.extent(0)*EXPAND_FACT);
    //  //thandle.set_nnzL(newsize);
    //}
    for (size_type k = 0; k < lenl; ++k) {
      L_entries(cntL) = h_iL(k);
      cntL++;
    }
#ifdef KEEP_DIAG
    //L diag entry
    L_entries(cntL) = i;
    cntL++;
#endif
    L_row_map(i+1) = cntL;
  }//End main loop i

  thandle.set_nnzL(cntL);
  thandle.set_nnzU(cntU);

  //Level scheduling on L
  level_sched (thandle, L_row_map, L_entries, nrows, level_list, level_ptr, level_idx, nlev);  
  
  thandle.set_symbolic_complete();

  // Output check
#ifdef SYMBOLIC_OUTPUT_INFO
  std::cout << "  ILU(k) fill_level: " << fill_lev << std::endl;
  std::cout << "  symbolic complete: " << thandle.is_symbolic_complete() << std::endl;
  std::cout << "  num levels: " << thandle.get_num_levels() << std::endl;
  std::cout << "  max num rows levels: " << thandle.get_level_maxrows() << std::endl;

  std::cout << "  iluk_symbolic result: " << std::endl;

  std::cout << "  level_list = ";
  for ( size_type i = 0; i < nrows; ++i )
  { std::cout << level_list(i) << " "; }
  std::cout << std::endl;

  std::cout << "  level_ptr = ";
  for ( size_type i = 0; i < nlev+1; ++i )
  { std::cout << level_ptr(i) << " "; }
  std::cout << std::endl;

  std::cout << "  level_idx = ";
  for ( size_type i = 0; i < nrows; ++i )
  { std::cout << level_idx(i) << " "; }
  std::cout << std::endl;

  std::cout << "  nnzL: " << thandle.get_nnzL() << std::endl;
  std::cout << "  L_row_map = ";
  for ( size_type i = 0; i < nrows+1; ++i )
  { std::cout << L_row_map(i) << " "; }
  std::cout << std::endl;

  std::cout << "  L_entries = ";
  for ( size_type i = 0; i < thandle.get_nnzL(); ++i )
  { std::cout << L_entries(i) << " "; }
  std::cout << std::endl;

  std::cout << "  nnzU: " << thandle.get_nnzU() << std::endl;
  std::cout << "  U_row_map = ";
  for ( size_type i = 0; i < nrows+1; ++i )
  { std::cout << U_row_map(i) << " "; }
  std::cout << std::endl;

  std::cout << "  U_entries = ";
  for ( size_type i = 0; i < thandle.get_nnzU(); ++i )
  { std::cout << U_entries(i) << " "; }
  std::cout << std::endl;
#endif

  Kokkos::deep_copy(dlevel_ptr, level_ptr);
  Kokkos::deep_copy(dlevel_idx, level_idx);
  Kokkos::deep_copy(dlevel_list, level_list);
  
  Kokkos::deep_copy(L_row_map_d, L_row_map);
  Kokkos::deep_copy(L_entries_d, L_entries);
  Kokkos::deep_copy(U_row_map_d, U_row_map);
  Kokkos::deep_copy(U_entries_d, U_entries);
 }
} // end iluk_symbolic

} // namespace Experimental
} // namespace Impl
} // namespace KokkosSparse

#endif
