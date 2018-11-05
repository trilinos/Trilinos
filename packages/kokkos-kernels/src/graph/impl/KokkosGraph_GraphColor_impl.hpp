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

#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <vector>
#include "KokkosGraph_GraphColorHandle.hpp"

#include <bitset>

#ifndef _KOKKOSCOLORINGIMP_HPP
#define _KOKKOSCOLORINGIMP_HPP


#define EBCOLORING_HIGHER_QUALITY        //suggested
namespace KokkosGraph{

namespace Impl{

#define VB_COLORING_FORBIDDEN_SIZE 64
#define VBBIT_COLORING_FORBIDDEN_SIZE 64
/*! \brief Base class for graph coloring purposes.
 *  Each color represents the set of the vertices that are independent,
 *  e.g. no vertex having same color shares an edge.
 *  General aim is to find the minimum number of colors, minimum number of independent sets.
 */
template <typename HandleType, typename lno_row_view_t_, typename lno_nnz_view_t_>
class GraphColor {
public:

  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;

  typedef typename HandleType::color_t color_t;
  typedef typename HandleType::color_view_t color_view_t;

  typedef typename HandleType::size_type size_type;
  typedef typename HandleType::nnz_lno_t nnz_lno_t;


  typedef typename in_lno_row_view_t::HostMirror row_lno_host_view_t; //Host view type


  typedef typename in_lno_nnz_view_t::HostMirror nnz_lno_host_view_t; //Host view type

  typedef typename HandleType::color_host_view_t color_host_view_t; //Host view type

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename HandleType::const_size_type const_size_type;
  typedef typename lno_row_view_t_::const_type const_lno_row_view_t;

  typedef typename lno_nnz_view_t_::const_type const_lno_nnz_view_t;
  typedef typename lno_nnz_view_t_::non_const_type non_const_lno_nnz_view_t;

protected:
  nnz_lno_t nv; //# vertices
  size_type ne; //# edges
  const_lno_row_view_t xadj; //rowmap
  const_lno_nnz_view_t adj; // entries
  const_lno_nnz_view_t kok_src, kok_dst; //Edge list storage of the graph
  HandleType *cp;

public:
  /**
   * \brief GraphColor constructor.
   * \param nv_: number of vertices in the graph
   * \param ne_: number of edges in the graph
   * \param row_map: the xadj array of the graph. Its size is nv_ +1
   * \param entries: adjacency array of the graph. Its size is ne_
   * \param coloring_handle: GraphColoringHandle object that holds the specification about the graph coloring,
   *    including parameters.
   */
  GraphColor(
      nnz_lno_t nv_,
      size_type ne_,
      const_lno_row_view_t row_map,
      const_lno_nnz_view_t entries,
      HandleType *coloring_handle):
        nv (nv_), ne(ne_),xadj(row_map), adj (entries),
        kok_src(), kok_dst(), cp(coloring_handle){}

  /** \brief GraphColor destructor.
   */
  virtual ~GraphColor(){}


  /** \brief Function to color the vertices of the graphs. This is the base class,
   * therefore, it only performs sequential coloring on the host device, ignoring the execution space.
   * \param colors is the output array corresponding the color of each vertex.Size is this->nv.
   *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
   *   it should be all initialized with zeros. Any positive value in the given array, will make the
   *   algorithm to assume that the color is fixed for the corresponding vertex.
   * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
   */
  virtual void color_graph(
      color_view_t d_colors,
      int &num_phases){

    num_phases = 1;


    color_host_view_t colors = Kokkos::create_mirror_view (d_colors);
    typename const_lno_row_view_t::HostMirror h_xadj = Kokkos::create_mirror_view (this->xadj);
    typename const_lno_nnz_view_t::HostMirror h_adj = Kokkos::create_mirror_view (this->adj);

    //typename nnz_lno_host_view_t::HostMirror::HostMirror::HostMirror h_adj = tmp;

    Kokkos::deep_copy (h_xadj, this->xadj);
    Kokkos::deep_copy (h_adj, this->adj);

    MyExecSpace::fence();




    //create a ban color array to keep track of
    //which colors have been taken by the neighbor vertices.
    nnz_lno_t *banned_colors = new nnz_lno_t[this->nv];


    for (nnz_lno_t i = 0; i < this->nv; ++i) banned_colors[i] = 0;

    color_t max_color = 0;
    //traverse vertices greedily
    for (nnz_lno_t i = 0; i < this->nv; ++i){

      size_type nbegin = h_xadj(i);
      size_type nend = h_xadj(i + 1);
      //std::cout << "nb:" << nbegin << " ne:" << nend << std::endl;
      //check the colors of neighbors
      for (size_type j = nbegin; j < nend; ++j){
        nnz_lno_t n = h_adj(j);
        if (n >= nv) continue;
        //set the banned_color of the color of the neighbor vertex to my vertex index.
        //the entries in the banned_color array that has my vertex index will be the set of prohibeted colors.
        banned_colors[colors(n)] = i;
      }
      //check the prohibeted colors, and pick the first available one.
      for (color_t j = 1; j <= max_color; ++j) {
        if(banned_colors[j] != i){
          colors(i) = j;
          break;
        }
      }
      //if no color is available, pick a new color.
      if (colors(i) == 0) colors(i) = ++max_color;
    }
    delete [] banned_colors;

    Kokkos::deep_copy (d_colors, colors); // Copy from host to device.
  }


  /** \brief Function to distance-2 color the vertices of the graphs. This is the base class,
   * therefore, it only performs sequential coloring on the host device, ignoring the execution space.
   * \param colors is the output array corresponding the color of each vertex. Size is this->nv.
   *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
   *   it should be all initialized with zeros. Any positive value in the given array, will make the
   *   algorithm to assume that the color is fixed for the corresponding vertex.
   * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
   */
  template <typename col_view_t_, typename lno_colnnz_view_t_>
  void d2_color_graph(
      color_view_t d_colors,
      int &num_phases,
      nnz_lno_t num_cols,
      col_view_t_ col_map, lno_colnnz_view_t_ col_entries){

    num_phases = 1;
    color_host_view_t colors = Kokkos::create_mirror_view (d_colors);

    row_lno_host_view_t h_xadj = Kokkos::create_mirror_view (this->xadj);
    typename const_lno_nnz_view_t::HostMirror h_adj = Kokkos::create_mirror_view (this->adj);
    Kokkos::deep_copy (h_xadj, this->xadj);
    Kokkos::deep_copy (h_adj, this->adj);
    MyExecSpace::fence();


    typedef typename col_view_t_::HostMirror col_host_view_t; //Host view type
    typedef typename lno_colnnz_view_t_::HostMirror col_nnz_view_t; //Host view type
    col_host_view_t h_c_xadj = Kokkos::create_mirror_view (col_map);
    col_nnz_view_t h_c_adj = Kokkos::create_mirror_view (col_entries);
    Kokkos::deep_copy (h_c_xadj, col_map);
    Kokkos::deep_copy (h_c_adj, col_entries);
    MyExecSpace::fence();


    Kokkos::Impl::Timer timer;
    //create a ban color array to keep track of
    //which colors have been taking by the neighbor vertices.
    nnz_lno_t *banned_colors = new nnz_lno_t[this->nv];

    //size_t howmany = 0;

    for (nnz_lno_t i = 0; i < this->nv; ++i) {
      colors(i) = banned_colors[i] = 0;
    }

    color_t max_color = 0;
    //traverse vertices greedily
    for (nnz_lno_t i = 0; i < this->nv; ++i){

      size_type col_begin = h_xadj(i);
      size_type col_end = h_xadj(i + 1);
      //check the colors of neighbors
      for (size_type j = col_begin; j < col_end; ++j){
        nnz_lno_t col = h_adj(j);
        if (col >= num_cols) continue;
        //set the banned_color of the color of the neighbor vertex to my vertex index.
        //the entries in the banned_color array that has my vertex index will be the set of prohibeted colors.
        //banned_colors[colors(col)] = i;


        size_type d2_rowbegin = h_c_xadj(col);
        size_type d2_rowend = h_c_xadj(col + 1);

        for (size_type j2 = d2_rowbegin; j2 < d2_rowend; ++j2){

          //howmany++;
          nnz_lno_t d2row = h_c_adj(j2);

          if (d2row >= nv || d2row == i) continue;
          banned_colors[colors(d2row)] = i;
        }
      }
      //check the prohibeted colors, and pick the first available one.
      for (color_t j = 1; j <= max_color; ++j) {
        if(banned_colors[j] != i){
          colors(i) = j;
          break;
        }
      }
      //if no color is available, pick a new color.
      if (colors(i) == 0) colors(i) = ++max_color;
    }
    //std::cout << "howmany:" << howmany << std::endl;
    delete [] banned_colors;
    //std::cout << "D2 Seq Color time:" << timer.seconds()  << std::endl;
    Kokkos::deep_copy (d_colors, colors); // Copy from host to device.
  }


  /** \brief Function to distance-2 color the vertices of the graphs. This is the base class,
   * therefore, it only performs sequential coloring on the host device, ignoring the execution space.
   * \param colors is the output array corresponding the color of each vertex. Size is this->nv.
   *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
   *   it should be all initialized with zeros. Any positive value in the given array, will make the
   *   algorithm to assume that the color is fixed for the corresponding vertex.
   * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
   */
  virtual void d2_color_graph2(
      color_view_t d_colors,
      int &num_phases){

    num_phases = 1;
    color_host_view_t colors = Kokkos::create_mirror_view (d_colors);
    typename const_lno_row_view_t::HostMirror h_xadj = Kokkos::create_mirror_view (this->xadj);
    typename const_lno_nnz_view_t::HostMirror h_adj = Kokkos::create_mirror_view (this->adj);
    Kokkos::deep_copy (h_xadj, this->xadj);
    Kokkos::deep_copy (h_adj, this->adj);
    MyExecSpace::fence();

    //create a ban color array to keep track of
    //which colors have been taking by the neighbor vertices.
    //row_index_type *banned_colors = new row_index_type[this->nv];
    //std::vector <nonzero_index_type> banned_colors(this->nv);

    std::vector <nnz_lno_t> vertex_cs(this->nv);

    std::vector <nnz_lno_t> col_d1_max_cs(this->nv);
    std::vector <color_t> col_d1_forbid_colors(this->nv);
    std::vector <nnz_lno_t> col_d1forbid_cs(this->nv);






    for (nnz_lno_t i = 0; i < this->nv; ++i) {
      col_d1_max_cs[i] = col_d1_forbid_colors[i] =  col_d1forbid_cs[i] = 0;
      colors(i) = 0;
      vertex_cs[i] = 0;
    }

    //std::cout << "coloring now" << std::endl;
    //color_t max_color = 0;
    //traverse vertices greedily
    for (nnz_lno_t i = 0; i < this->nv; ++i){

      size_type nbegin = h_xadj(i);
      nnz_lno_t row_size = h_xadj(i + 1) - nbegin;

      nnz_lno_t my_cs = vertex_cs[i];
      nnz_lno_t my_max_ncs = col_d1_max_cs[i];


      //now find maximum d1 cs
      for (nnz_lno_t jj = 0; jj < row_size; ++jj){
        size_type j = jj +  nbegin;
        nnz_lno_t n = h_adj(j);
        if (n >= nv || n == i) continue;
        if (col_d1forbid_cs[n] > my_cs) my_cs = col_d1forbid_cs[n];
        if (vertex_cs[n] > my_max_ncs) my_max_ncs = vertex_cs[n];
      }
      col_d1_max_cs[i] = my_max_ncs;



      while (true){

        color_t my_d2_forbidden = 0;
        //color_type my_d1_forbidden = 0;

        for (nnz_lno_t jj = 0; jj < row_size; ++jj){
          size_type j = jj +  nbegin;
          nnz_lno_t n = h_adj(j);
          if (n >= nv || n == i) continue;

          if (my_cs == vertex_cs[n]){
            my_d2_forbidden = my_d2_forbidden | colors(n);
          }
          if (my_cs == col_d1forbid_cs[n]){
            my_d2_forbidden = my_d2_forbidden | col_d1_forbid_colors[n];
          }
          else if (my_cs > col_d1forbid_cs[n] && my_cs <= col_d1_max_cs[n]){
            size_type nbegin2 = h_xadj(n);
            nnz_lno_t row_size2 = h_xadj(n + 1) - nbegin2;
            for (nnz_lno_t jj2 = 0; jj2 < row_size2; ++jj2){
              size_type j2 = jj2 +  nbegin2;
              nnz_lno_t n2 = h_adj(j2);
              if (n2 >= nv || n2 == i) continue;
              if (my_cs == vertex_cs[n2]){
                my_d2_forbidden = my_d2_forbidden | colors(n2);
              }
            }
          }
          if (~my_d2_forbidden == 0) {
            break;
          }
        }

        color_t avaliable_colors = ~my_d2_forbidden;
        if (avaliable_colors == 0) {
          ++my_cs;
        }
        else {
          colors(i) = avaliable_colors & (-avaliable_colors);
          vertex_cs[i] = my_cs;


          for (nnz_lno_t jj = 0; jj < row_size; ++jj){
            size_type j = jj +  nbegin;
            nnz_lno_t n = h_adj(j);
            if (n >= nv || n == i) continue;
            if (my_cs == col_d1forbid_cs[n]){
              col_d1_forbid_colors[n] = col_d1_forbid_colors[n] | colors(i);
              if (~col_d1_forbid_colors[n] == 0){
                nnz_lno_t my_d1_cs = my_cs + 1;
                color_t my_d1_forbidden = 0;

                size_type nbegin2 = h_xadj(n);
                nnz_lno_t row_size2 = h_xadj(n + 1) - nbegin2;

                while (true){
                  for (nnz_lno_t kk = 0; kk < row_size2; ++kk){
                    size_type ll = kk +  nbegin2;
                    nnz_lno_t n2 = h_adj(ll);
                    if (n2 >= nv || n2 == n) continue;
                    if (vertex_cs[n2] == my_d1_cs) my_d1_forbidden = my_d1_forbidden | colors(n2);
                    if (~my_d1_forbidden == 0){
                      break;
                    }
                  }

                  if (~my_d1_forbidden == 0){
                    my_d1_cs++;
                  }
                  else {
                    col_d1forbid_cs[n] = my_d1_cs;
                    col_d1_forbid_colors[n] = my_d1_forbidden;
                    break;
                  }
                }
              }
            }
            if(my_cs > col_d1_max_cs[n]) col_d1_max_cs[n] = my_cs;
          }
          break;
        }
      }
    }

    int color_size = sizeof (color_t) * 8;

    for (nnz_lno_t ii = 0; ii < this->nv; ++ii){
      color_t val = colors(ii);
      if (val){
        //find the position in the bit.
        int i = 1;
        while ((val & 1) == 0) {
          ++i;
          val = val >> 1;
        }
        colors(ii) = i + vertex_cs[ii] * color_size;
      }
    }

    Kokkos::deep_copy (d_colors, colors); // Copy from host to device.
  }

};


template <typename HandleType, typename lno_row_view_t_, typename lno_nnz_view_t_>
class GraphColor2: public GraphColor <HandleType,lno_row_view_t_,lno_nnz_view_t_> {
public:

  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef typename HandleType::color_view_t color_view_t;


  typedef typename HandleType::size_type size_type;
  typedef typename in_lno_row_view_t::HostMirror row_lno_host_view_t; //Host view type


  typedef typename in_lno_nnz_view_t::HostMirror nnz_lno_host_view_t; //Host view type



  typedef typename HandleType::nnz_lno_t nnz_lno_t;
  //typedef typename idx_edge_memory_traits::MemorySpace MyEdgeMemorySpace;



  typedef typename HandleType::color_t color_t;
  typedef typename HandleType::color_host_view_t color_host_view_t; //Host view type


  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;






  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;
  //typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;


  typedef typename in_lno_nnz_view_t::const_type const_lno_nnz_view_t;
  //typedef typename in_lno_nnz_view_t::non_const_type non_const_lno_nnz_view_t;


public:
  /**
   * \brief GraphColor constructor.
   * \param nv_: number of vertices in the graph
   * \param ne_: number of edges in the graph
   * \param row_map: the xadj array of the graph. Its size is nv_ +1
   * \param entries: adjacency array of the graph. Its size is ne_
   * \param coloring_handle: GraphColoringHandle object that holds the specification about the graph coloring,
   *    including parameters.
   */
  GraphColor2(
      nnz_lno_t nv_,
      size_type ne_,
      const_lno_row_view_t row_map,
      const_lno_nnz_view_t entries,
      HandleType *coloring_handle):
        GraphColor<HandleType,lno_row_view_t_,lno_nnz_view_t_>(nv_, ne_, row_map, entries, coloring_handle)
        {}

  /** \brief GraphColor destructor.
   */
  virtual ~GraphColor2(){}

  /** \brief Function to distance-2 color the vertices of the graphs. This is the base class,
   * therefore, it only performs sequential coloring on the host device, ignoring the execution space.
   * \param colors is the output array corresponding the color of each vertex. Size is this->nv.
   *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
   *   it should be all initialized with zeros. Any positive value in the given array, will make the
   *   algorithm to assume that the color is fixed for the corresponding vertex.
   * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
   */
  template <typename col_view_t_, typename lno_colnnz_view_t_>
  void d2_color_graph(
      color_view_t d_colors,
      int &num_phases,
      nnz_lno_t num_cols,
      col_view_t_ col_map, lno_colnnz_view_t_ col_entries){

    num_phases = 1;
    color_host_view_t colors = Kokkos::create_mirror_view (d_colors);
    row_lno_host_view_t h_xadj = Kokkos::create_mirror_view (this->xadj);
    typename const_lno_nnz_view_t::HostMirror h_adj = Kokkos::create_mirror_view (this->adj);
    Kokkos::deep_copy (h_xadj, this->xadj);
    Kokkos::deep_copy (h_adj, this->adj);
    MyExecSpace::fence();

    typedef typename col_view_t_::HostMirror col_host_view_t; //Host view type
    typedef typename lno_colnnz_view_t_::HostMirror col_nnz_view_t; //Host view type
    col_host_view_t h_c_xadj = Kokkos::create_mirror_view (col_map);
    col_nnz_view_t h_c_adj = Kokkos::create_mirror_view (col_entries);
    Kokkos::deep_copy (h_c_xadj, col_map);
    Kokkos::deep_copy (h_c_adj, col_entries);
    MyExecSpace::fence();


    //create a ban color array to keep track of
    //which colors have been taking by the neighbor vertices.
    //row_index_type *banned_colors = new row_index_type[this->nv];
    //std::vector <nonzero_index_type> banned_colors(this->nv);

    Kokkos::Impl::Timer timer;

    std::vector <nnz_lno_t> row_cs(this->nv);
    std::vector <nnz_lno_t> col_d1_max_cs(num_cols);
    std::vector <color_t> col_d1_forbid_colors(num_cols);
    std::vector <nnz_lno_t> col_d1forbid_cs(num_cols);

    for (nnz_lno_t i = 0; i < this->nv; ++i) {
      colors(i) = row_cs[i] = 0;
    }

    for (nnz_lno_t i = 0; i < num_cols; ++i) {
      col_d1_max_cs[i] = col_d1_forbid_colors[i] =  col_d1forbid_cs[i] = 0;
    }

    //color_t max_color = 0;
    //traverse vertices greedily

    //size_t howmany = 0;
    for (nnz_lno_t i = 0; i < this->nv; ++i){

      //std::cout << "i:" << i << std::endl;
      size_type colbegin = h_xadj(i);
      nnz_lno_t row_size = h_xadj(i + 1) - colbegin;

      nnz_lno_t my_row_cs = row_cs[i];
      //nonzero_index_type my_max_ncs = col_d1_max_cs[i];

      //now find maximum d1 cs
      for (nnz_lno_t jj = 0; jj < row_size; ++jj){
        //howmany++;
        size_type j = jj +  colbegin;
        nnz_lno_t col = h_adj(j);
        if (col >= num_cols) continue;
        if (col_d1forbid_cs[col] > my_row_cs) my_row_cs = col_d1forbid_cs[col];
      }

      while (true){
        color_t my_d2_forbidden = 0;

        for (nnz_lno_t jj = 0; jj < row_size; ++jj){


          //howmany++;
          size_type j = jj +  colbegin;

          nnz_lno_t col = h_adj(j);

          //std::cout << "\ti:" << i << " neighbor:" << jj << " col:" << col << std::endl;
          if (col >= num_cols) continue;

          //TODO: do we forbid distance-1 neighbors too?
          /*
          if (my_row_cs == row_cs[col]){
            my_d2_forbidden = my_d2_forbidden | colors(col);
          }
          */
          if (my_row_cs == col_d1forbid_cs[col]){
            my_d2_forbidden = my_d2_forbidden | col_d1_forbid_colors[col];
          }
          else if (my_row_cs > col_d1forbid_cs[col] && my_row_cs <= col_d1_max_cs[col]){
            //howmany++;
            size_type row_begin = h_c_xadj(col);
            nnz_lno_t col_size = h_c_xadj(col + 1) - row_begin;
            for (nnz_lno_t jj2 = 0; jj2 < col_size; ++jj2){
              size_type j2 = jj2 +  row_begin;
              nnz_lno_t d2_row = h_c_adj(j2);
              if (d2_row >= this->nv || d2_row == i) continue;
              if (my_row_cs == row_cs[d2_row]){
                my_d2_forbidden = my_d2_forbidden | colors(d2_row);
              }
            }
          }
          if (~my_d2_forbidden == 0) {
            break;
          }
        }

        color_t avaliable_colors = ~my_d2_forbidden;
        if (avaliable_colors == 0) {
          ++my_row_cs;
        }
        else {
          colors(i) = avaliable_colors & (-avaliable_colors);
          row_cs[i] = my_row_cs;


          //std::cout << "\ti:" << i << " found color cs:" << my_row_cs << " colors:" << colors(i) << std::endl;


          for (nnz_lno_t jj = 0; jj < row_size; ++jj){
            //howmany++;
            size_type j = jj +  colbegin;
            nnz_lno_t col = h_adj(j);
            if (col >= num_cols) continue;
            //std::cout << "\ti:" << i << " sending color to neighbor:" << jj << " col:" << col << std::endl;


            if (my_row_cs == col_d1forbid_cs[col]){
              col_d1_forbid_colors[col] = col_d1_forbid_colors[col] | colors(i);
              if (~col_d1_forbid_colors[col] == 0){
                nnz_lno_t col_d1_cs = my_row_cs + 1;
                color_t col_d1_forbidden = 0;

                if (col_d1_cs <= col_d1_max_cs[col]){
                  size_type row_begin2 = h_c_xadj(col);
                  nnz_lno_t col_size2 = h_c_xadj(col + 1) - row_begin2;

                  while (true){
                    col_d1_forbidden = 0;
                    //std::cout << "\ti:" << i << " creating forbid for col:" << col << " in cs:" << col_d1_cs << std::endl;
                    for (nnz_lno_t kk = 0; kk < col_size2; ++kk){
                      //howmany++;
                      size_type ll = kk +  row_begin2;
                      nnz_lno_t d2_row = h_c_adj(ll);
                      if (d2_row >= this->nv) continue;
                      if (row_cs[d2_row] == col_d1_cs) col_d1_forbidden = col_d1_forbidden | colors(d2_row);
                      if (~col_d1_forbidden == 0){
                        break;
                      }
                    }

                    if (~col_d1_forbidden == 0){
                      col_d1_cs++;
                    }
                    else {
                      col_d1forbid_cs[col] = col_d1_cs;
                      col_d1_forbid_colors[col] = col_d1_forbidden;
                      break;
                    }
                  }
                }
                else {
                  col_d1forbid_cs[col] = col_d1_cs;
                  col_d1_forbid_colors[col] = col_d1_forbidden;
                }
              }
            }
            if(my_row_cs > col_d1_max_cs[col]) col_d1_max_cs[col] = my_row_cs;
          }
          break;
        }
      }
    }

    //std::cout << "howmany:" << howmany << std::endl;

    int color_size = sizeof (color_t) * 8;

    for (nnz_lno_t ii = 0; ii < this->nv; ++ii){
      color_t val = colors(ii);
      if (val){
        //find the position in the bit.
        int i = 1;
        while ((val & 1) == 0) {
          ++i;
          val = val >> 1;
        }
        colors(ii) = i + row_cs[ii] * color_size;
      }
    }
    Kokkos::deep_copy (d_colors, colors); // Copy from host to device.
  }

};

/*! \brief Class for the vertex based graph coloring algorithms.
 *  They work better on CPUs and Xeon Phis, but edge-based ones are better on GPUs.
 *  Includes 3 algorithms:
 *  VB: Speculative parallel vertex based algorithm using a forbidden array of size 64 per thread.
 *  Best on Xeon Phi among the vertex based algorithms.
 *  VBBIT: Speculative parallel vertex based using a long integer for forbidden colors per vertex.
 *  Best on GPUs among the vertex based algorithms.
 *  VBCS: Speculative parallel vertex based using color set implementation.
 */
template <typename HandleType, typename lno_row_view_t_, typename lno_nnz_view_t_>
class GraphColor_VB:public GraphColor <HandleType,lno_row_view_t_,lno_nnz_view_t_>{
public:

  typedef long long int ban_type;

  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef typename HandleType::color_view_t color_view_type;


  typedef typename HandleType::size_type size_type;
  typedef typename lno_row_view_t_::device_type row_lno_view_device_t;

  typedef typename HandleType::nnz_lno_t nnz_lno_t;

  typedef typename HandleType::color_t color_t;
  typedef typename HandleType::color_host_view_t color_host_view_t; //Host view type

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename Kokkos::View<nnz_lno_t, row_lno_view_device_t> single_dim_index_view_type;
  //typedef typename Kokkos::View<row_index_type, Kokkos::MemoryUnmanaged> um_array_type;
  typedef typename single_dim_index_view_type::HostMirror single_dim_index_host_view_type; //Host view type

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  typedef typename HandleType::size_type_temp_work_view_t size_type_temp_work_view_t;
  typedef typename HandleType::size_type_persistent_work_view_t size_type_persistent_work_view_t;


  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;



  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;


  typedef typename lno_nnz_view_t_::const_type const_lno_nnz_view_t;
  typedef typename lno_nnz_view_t_::non_const_type non_const_lno_nnz_view_t;


protected:



  //typedef Kokkos::View<idx, /*Kokkos::Serial::array_layout,*//* Kokkos::Serial,*/ Kokkos::MemoryUnmanaged> um_array_type;


  bool _serialConflictResolution; //if true use serial conflict resolution
  bool _ticToc; //if true print info in each step
  char _conflictlist; //0 for no conflictlist, 1 for atomic, 2 for pps

  double _pps_ratio; //the minimum number of reduction on the size of the conflictlist to create a new conflictlist
  nnz_lno_t _min_vertex_cut_off; //minimum number of vertices to reduce the conflictlist further.
  bool _edge_filtering; //if true, edge-filtering is applied by swaps on adjacency array.
  int _chunkSize; //the size of the minimum work unit assigned to threads. Changes the convergence on GPUs
  char _use_color_set; //the VB algorithm type.
                        // 0 for VB:
                        // 1: for VBCS
                        // 2: for VBBIT

  int _max_num_iterations;

public:
  /**
   * \brief GraphColor_VB constructor.
   * \param nv_: number of vertices in the graph
   * \param ne_: number of edges in the graph
   * \param row_map: the xadj array of the graph. Its size is nv_ +1
   * \param entries: adjacency array of the graph. Its size is ne_
   * \param coloring_handle: GraphColoringHandle object that holds the specification about the graph coloring,
   *    including parameters.
   */
  GraphColor_VB(
      nnz_lno_t nv_, size_type ne_,
      const_lno_row_view_t row_map, const_lno_nnz_view_t entries,
      HandleType *coloring_handle):
    GraphColor<HandleType,lno_row_view_t_,lno_nnz_view_t_>(nv_, ne_, row_map, entries, coloring_handle),
    _serialConflictResolution(coloring_handle->get_serial_conflict_resolution()),
    _ticToc(coloring_handle->get_tictoc()),
    _conflictlist(),
    _pps_ratio(coloring_handle->get_min_reduction_for_conflictlist()),
    _min_vertex_cut_off(coloring_handle->get_min_elements_for_conflictlist()),
    _edge_filtering(coloring_handle->get_vb_edge_filtering()),
    _chunkSize(coloring_handle->get_vb_chunk_size()),
    _use_color_set(),
    _max_num_iterations(coloring_handle->get_max_number_of_iterations())
    {
      switch (coloring_handle->get_coloring_algo_type()){
      case COLORING_VB:
        this->_use_color_set = 0;
        break;
      case COLORING_VBBIT:
        this->_use_color_set = 2;
        break;
      case COLORING_VBCS:
        this->_use_color_set = 1;
        break;
      default: //cannnot get in here.
        this->_use_color_set = 0;
        break;

      }

      switch (coloring_handle->get_conflict_list_type()){
      case COLORING_NOCONFLICT:
        this->_conflictlist = 0;
        break;
      case COLORING_ATOMIC:
        this->_conflictlist = 1;
        break;
      case COLORING_PPS:
        this->_conflictlist = 2;
        break;
      }
    }

  /** \brief GraphColor_VB destructor.
    */
  virtual ~GraphColor_VB(){}

  /** \brief Function to color the vertices of the graphs. Performs a vertex-based coloring.
   * \param colors is the output array corresponding the color of each vertex. Size is this->nv.
   *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
   *   it should be all initialized with zeros. Any positive value in the given array, will make the
   *   algorithm to assume that the color is fixed for the corresponding vertex.
   * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
   */
  virtual void color_graph(color_view_type colors,int &num_loops){

    if (this->_ticToc) {
      std::cout
          << "\tVB params:" << std::endl
          << "\tuseConflictList:" << int (this->_conflictlist) << std::endl
          << "\talgorithm:" << (int)this->_use_color_set << std::endl
          << "\tserialConflictResolution:"  << (int) this->_serialConflictResolution << std::endl
          << "\tticToc:" << (int) this->_ticToc << std::endl
          << "\tuse_color_set:" << (int) this->_use_color_set << std::endl
          << "\tpps_ratio:" << this->_pps_ratio << std::endl
          << "\tmin_vertex_cut_off:" << this->_min_vertex_cut_off << std::endl
          << "\tedge_filtering:" << (int)this->_edge_filtering << std::endl
          << "\tmax_num_iterations:" << this->_max_num_iterations << std::endl
          << "\tchunkSize:" << this->_chunkSize << std::endl;
    }

    //if the edge filtering is selected, then we do swaps on the adj array.
    //to not to touch the given one, we copy the adj array.
    nnz_lno_temp_work_view_t adj_copy;

    //if we use edge-filtering, we perform swaps.
    //We need to copy the adjacency array so that we dont harm the original one.
    if (this->_edge_filtering){
      adj_copy = nnz_lno_temp_work_view_t(Kokkos::ViewAllocateWithoutInitializing("adj copy"), this->ne);
      Kokkos::deep_copy(adj_copy, this->adj);
    }

    //if color set algorithm is used, we need one more array to represent the range.
    nnz_lno_temp_work_view_t vertex_color_set;
    if (this->_use_color_set == 1){
      vertex_color_set = nnz_lno_temp_work_view_t("colorset", this->nv);
    }

    //the conflictlist
    nnz_lno_temp_work_view_t current_vertexList =
        nnz_lno_temp_work_view_t(Kokkos::ViewAllocateWithoutInitializing("vertexList"), this->nv);

    //init vertexList sequentially.
    Kokkos::parallel_for("KokkosGraph::GraphColoring::InitList",
        my_exec_space(0, this->nv), functorInitList<nnz_lno_temp_work_view_t> (current_vertexList));


    // the next iteration's conflict list
    nnz_lno_temp_work_view_t next_iteration_recolorList;
    // the size of the current conflictlist

    //the size of the next iteration's conflictlist
    single_dim_index_view_type next_iteration_recolorListLength;
    //if parallel prefix sum is selected instead of atomic operations,
    //we need one more work array to do the prefix sum.
    nnz_lno_temp_work_view_t pps_work_view;

    // if a conflictlist is used
    if (this->_conflictlist > 0){
      // Vertices to recolor. Will swap with vertexList.
      next_iteration_recolorList = nnz_lno_temp_work_view_t(Kokkos::ViewAllocateWithoutInitializing("recolorList"), this->nv);
      next_iteration_recolorListLength = single_dim_index_view_type("recolorListLength");
      if (this->_conflictlist == 2) {
        pps_work_view = nnz_lno_temp_work_view_t("pps_view", this->nv);
      }
    }

    nnz_lno_t numUncolored = this->nv;
    nnz_lno_t current_vertexListLength = this->nv;


    double t, total=0.0;
    Kokkos::Impl::Timer timer;


    int iter=0;
    for (; (iter < this->_max_num_iterations) && (numUncolored>0); iter++){

      if (this->_edge_filtering)
      {
        // First color greedy speculatively,
        //some conflicts expected
        this->colorGreedyEF(
            this->xadj,
            adj_copy,
            colors,
            vertex_color_set,
            current_vertexList,
            current_vertexListLength);
      }
      else {
        // First color greedy speculatively,
        //some conflicts expected
        this->colorGreedy(
            this->xadj,
            this->adj,
            colors,
            vertex_color_set,
            current_vertexList,
            current_vertexListLength);
      }

      MyExecSpace::fence();

      if (this->_ticToc){
        t = timer.seconds();
        total += t;
        std::cout << "\tTime speculative greedy phase " << iter << " : " << t << std::endl;
        timer.reset();
      }

      bool swap_work_arrays = true;
      if (this->_edge_filtering)
      {
      numUncolored = this->findConflicts(
          swap_work_arrays,
          this->xadj, adj_copy,
          colors, vertex_color_set,
          current_vertexList, current_vertexListLength,
          next_iteration_recolorList, next_iteration_recolorListLength,
          pps_work_view);
      }
      else {
        numUncolored = this->findConflicts(
            swap_work_arrays,
            this->xadj, this->adj,
            colors, vertex_color_set,
            current_vertexList, current_vertexListLength,
            next_iteration_recolorList, next_iteration_recolorListLength,
            pps_work_view);
      }

      MyExecSpace::fence();

      if (_ticToc){
        t = timer.seconds();
        total += t;
        std::cout << "\tTime conflict detection " << iter << " : " << t << std::endl;
        timer.reset();
      }

      if (this->_serialConflictResolution) break; // Break after first iteration.
      if (this->_conflictlist && swap_work_arrays && (iter + 1< this->_max_num_iterations)){
        // Swap recolorList and vertexList
        nnz_lno_temp_work_view_t temp = current_vertexList;
        current_vertexList = next_iteration_recolorList;
        next_iteration_recolorList = temp;
        current_vertexListLength = numUncolored;
        next_iteration_recolorListLength = single_dim_index_view_type("recolorListLength");
      }

    }

    //if VBCS algorithm is used, the colors are converted back to original form.
    if (this->_use_color_set == 1){
      Kokkos::parallel_for("KokkosGraph::GraphColoring::SetFinalColors",
          my_exec_space(0, this->nv), set_final_colors (colors, vertex_color_set));
    }
    if (numUncolored > 0){

    	if (this->_edge_filtering)
    	{
    		// Resolve conflicts by recoloring in serial
    		this->resolveConflicts(
    				this->nv,
					this->xadj, adj_copy,
					colors,
					current_vertexList, current_vertexListLength
    		);
    	}
    	else {
    		// Resolve conflicts by recoloring in serial
    		this->resolveConflicts(
    				this->nv,
					this->xadj, this->adj,
					colors,
					current_vertexList, current_vertexListLength
    		);


    	}
    	MyExecSpace::fence();
    	if (_ticToc){
    		t = timer.seconds();
    		total += t;
    		std::cout << "\tTime serial conflict resolution: " << t << std::endl;
    	}
    }
    num_loops = iter;
  }    // color_graph (end)


private:
  /** \brief Performs speculative coloring based on the given colorings.
   *  \param xadj_: row map of the graph
   *  \param adj_: entries, columns of the graph
   *  \param vertex_colors_: colors corresponding to each vertex
   *  \param vertex_color_set: if VBCS is used, color set of each vertex
   *  \param current_vertexList_: current conflictlist
   *  \param current_vertexListLength_: size of current conflictlist
   */
  void colorGreedy(
      const_lno_row_view_t xadj_,
      const_lno_nnz_view_t adj_,
      color_view_type vertex_colors_,
      nnz_lno_temp_work_view_t vertex_color_set,
      nnz_lno_temp_work_view_t current_vertexList_,
      nnz_lno_t current_vertexListLength_) {

    nnz_lno_t chunkSize_ = this->_chunkSize; // Process chunkSize vertices in one chunk

    if (current_vertexListLength_ < 100*chunkSize_) chunkSize_ = 1;

    //if the algorithm VBBIT
    if (this->_use_color_set == 2) {

      functorGreedyColor_IMPLOG gc(
          this->nv,
          xadj_, adj_,
          vertex_colors_,  current_vertexList_,
          current_vertexListLength_, chunkSize_);
      Kokkos::parallel_for("KokkosGraph::GraphColoring::GreedyColor_IMPLOG",
          my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);

    }
    // VBCS algorithm
    else if (this->_use_color_set == 1){
      functorGreedyColor_IMP gc(
          this->nv,
          xadj_, adj_,
          vertex_colors_, vertex_color_set, current_vertexList_,
          current_vertexListLength_, chunkSize_);
      Kokkos::parallel_for("KokkosGraph::GraphColoring::GreedyColor_IMP",
          my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);

    }
    //VB algorithm
    else if (this->_use_color_set == 0)
    {

      functorGreedyColor  gc(
          this->nv,
          xadj_, adj_,
          vertex_colors_,
          current_vertexList_, current_vertexListLength_, chunkSize_);
      Kokkos::parallel_for("KokkosGraph::GraphColoring::GreedyColor",
          my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);

    }
  }      // colorGreedy (end)



  /** \brief Performs speculative coloring based on the given colorings.
   *  \param xadj_: row map of the graph
   *  \param adj_: entries, columns of the graph
   *  \param vertex_colors_: colors corresponding to each vertex
   *  \param vertex_color_set: if VBCS is used, color set of each vertex
   *  \param current_vertexList_: current conflictlist
   *  \param current_vertexListLength_: size of current conflictlist
   */
  void  colorGreedyEF(
      const_lno_row_view_t xadj_,
      nnz_lno_temp_work_view_t adj_,
      color_view_type vertex_colors_,
      nnz_lno_temp_work_view_t vertex_color_set,
      nnz_lno_temp_work_view_t current_vertexList_,
      nnz_lno_t current_vertexListLength_) {

//    std::cout << ">>> WCMCLEN GraphColor_VB::colorGreedyEF (KokkosGraph_GraphColor_impl.hpp)" << std::endl;

    nnz_lno_t chunkSize_ = this->_chunkSize; // Process chunkSize vertices in one chunk

    if (current_vertexListLength_ < 100*chunkSize_) chunkSize_ = 1;

    //if the algorithm VBBIT
    if (this->_use_color_set == 2) {

      //If edge filtering is applied

      functorGreedyColor_IMPLOG_EF gc(
          this->nv,
          xadj_, adj_,
          vertex_colors_, current_vertexList_,
          current_vertexListLength_, chunkSize_);
      Kokkos::parallel_for("KokkosGraph::GraphColoring::GreedyColor_IMPLOG_EF",
          my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);


    }
    // VBCS algorithm
    else if (this->_use_color_set == 1){
      functorGreedyColor_IMP_EF gc(
          this->nv,
          xadj_, adj_,
          vertex_colors_, vertex_color_set, current_vertexList_,
          current_vertexListLength_, chunkSize_);
      Kokkos::parallel_for("KokkosGraph::GraphColoring::GreedyColor_IMP_EF",
          my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);
    }
    //VB algorithm
    else if (this->_use_color_set == 0)
    {
//      std::cout << ">>> WCMCLEN GraphColor_VB::colorGreedyEF --- (call the functor)" << std::endl;
      functorGreedyColor_EF  gc(
          this->nv,
          xadj_, adj_,
          vertex_colors_,
          current_vertexList_, current_vertexListLength_, chunkSize_);
      Kokkos::parallel_for("KokkosGraph::GraphColoring::GreedyColor_EF",
          my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);

    }
  }

  /** \brief Performs conflict resolution
   *  \param swap_work_arrays: An output parameter whether indicating whether the work arrays should be swapped or not.
   *  \param xadj_: row map of the graph
   *  \param adj_: entries, columns of the graph
   *  \param vertex_colors_: colors corresponding to each vertex
   *  \param vertex_color_set: if VBCS is used, color set of each vertex
   *  \param current_vertexList_: current conflictlist
   *  \param current_vertexListLength_: size of current conflictlist
   *  \param next_iteration_recolorList_: current conflictlist
   *  \param next_iteration_recolorListLength_: size of current conflictlist
   *  \param pps_work_view: size of current conflictlist
   */
  template <typename adj_view_t>
  nnz_lno_t findConflicts(
      bool &swap_work_arrays,
      const_lno_row_view_t xadj_,
      adj_view_t adj_,
      color_view_type vertex_colors_,
      nnz_lno_temp_work_view_t vertex_color_set_,
      nnz_lno_temp_work_view_t current_vertexList_,
      nnz_lno_t current_vertexListLength_,
      nnz_lno_temp_work_view_t next_iteration_recolorList_,
      single_dim_index_view_type next_iteration_recolorListLength_,
      nnz_lno_temp_work_view_t pps_work_view) {

    swap_work_arrays = true;
    nnz_lno_t numUncolored = 0;
    if (this->_conflictlist == 0){
      if (this->_use_color_set == 0 || this->_use_color_set == 2){
        functorFindConflicts_No_Conflist<adj_view_t> conf( this->nv, xadj_, adj_, vertex_colors_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
      else {
        functorFindConflicts_No_Conflist_IMP<adj_view_t> conf(this->nv, xadj_, adj_,vertex_colors_, vertex_color_set_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
    }
    else if (this->_conflictlist == 2){ //IF PPS
      if (this->_use_color_set == 0 || this->_use_color_set == 2){
        // Check for conflicts. Compute numUncolored == numConflicts.
        functorFindConflicts_PPS<adj_view_t> conf(this->nv, xadj_, adj_,vertex_colors_,current_vertexList_,next_iteration_recolorList_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
      else {
        functorFindConflicts_PPS_IMP<adj_view_t> conf(this->nv,
            xadj_, adj_,vertex_colors_, vertex_color_set_,
            current_vertexList_,next_iteration_recolorList_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }


      if( numUncolored && (current_vertexListLength_ >= this->_min_vertex_cut_off) &&
          (double (numUncolored) / current_vertexListLength_  <  (1.0 - this->_pps_ratio))){
        if (this->_ticToc){
          std::cout << "\tcreating work array with pps current_vertexListLength_:" <<
              current_vertexListLength_ << " params->min_vertex_cut_off:" << this->_min_vertex_cut_off << std::endl;
        }
        single_dim_index_host_view_type h_numUncolored(&numUncolored);
        Kokkos::deep_copy (next_iteration_recolorListLength_, h_numUncolored);

        MyExecSpace::fence();

        Kokkos::parallel_scan ("KokkosGraph::GraphColoring::PrefixSum",
            my_exec_space(0, current_vertexListLength_),
            parallel_prefix_sum<nnz_lno_temp_work_view_t>(current_vertexList_, next_iteration_recolorList_, pps_work_view));

        MyExecSpace::fence();
        Kokkos::parallel_for ("KokkosGraph::GraphColoring::CreateNewWorkArray",
            my_exec_space(0, current_vertexListLength_),
            create_new_work_array<nnz_lno_temp_work_view_t>(current_vertexList_, next_iteration_recolorList_, pps_work_view));
      }
      else {
        swap_work_arrays = false;
      }
    }
    else { //IF ATOMIC
      if (this->_use_color_set == 0 || this->_use_color_set == 2){
        // Check for conflicts. Compute numUncolored == numConflicts.
        functorFindConflicts_Atomic<adj_view_t> conf(this->nv,
            xadj_, adj_,vertex_colors_,current_vertexList_,
            next_iteration_recolorList_, next_iteration_recolorListLength_);
        Kokkos::parallel_reduce("KokkosGraph::GraphColoring::FindConflictsAtomic",
            my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
      else {
        functorFindConflicts_Atomic_IMP<adj_view_t> conf(this->nv,
            xadj_, adj_,vertex_colors_, vertex_color_set_,
            current_vertexList_,next_iteration_recolorList_, next_iteration_recolorListLength_);
        Kokkos::parallel_reduce("KokkosGraph::GraphColoring::FindConflictsAtomic_IMP",
            my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
    }
    if (this->_ticToc){
      std::cout << "\tnumUncolored:" << numUncolored << std::endl;
    }
    return numUncolored;
  }

  /** \brief Performs serial conflict resolution
   *  \param nv: Number of vertices
   *  \param xadj_: row map of the graph
   *  \param adj_: entries, columns of the graph
   *  \param vertex_colors_: colors corresponding to each vertex
   *  \param current_vertexList_: current conflictlist
   *  \param current_vertexListLength_: size of current conflictlist
   */
  template <typename adj_view_t>
  void  resolveConflicts(
      nnz_lno_t _nv,
      const_lno_row_view_t xadj_,
      adj_view_t adj_,
      color_view_type vertex_colors_,
      nnz_lno_temp_work_view_t current_vertexList_,
      size_type current_vertexListLength_) {

    color_t *forbidden = new color_t[_nv];
    nnz_lno_t i=0;
    nnz_lno_t end = _nv;
    typename nnz_lno_temp_work_view_t::HostMirror h_recolor_list;

    if (this->_conflictlist){
      end = current_vertexListLength_;
      h_recolor_list = Kokkos::create_mirror_view (current_vertexList_);
      Kokkos::deep_copy (h_recolor_list, current_vertexList_);
    }
    color_host_view_t h_colors = Kokkos::create_mirror_view (vertex_colors_);
    typename const_lno_row_view_t::HostMirror h_idx = Kokkos::create_mirror_view (xadj_);
    typename adj_view_t::HostMirror h_adj = Kokkos::create_mirror_view (adj_);


    Kokkos::deep_copy (h_colors, vertex_colors_);
    Kokkos::deep_copy (h_idx, xadj_);
    Kokkos::deep_copy (h_adj, adj_);

    for (nnz_lno_t k=0; k <end; k++){
      if (this->_conflictlist){
        i = h_recolor_list(k);
      }
      else {
        // Check for uncolored vertices
        i = k;
      }
      if (h_colors(i) > 0) continue;
      for (size_type j=h_idx(i); j<h_idx(i+1); j++){
        if (h_adj(j) == i) continue; // Skip self-loops
        forbidden[h_colors(h_adj(j))] = i;
      }
      // color vertex i with smallest available color
      int c=1;
      while (forbidden[c]==i) c++;
      h_colors(i) = c;
    }
    Kokkos::deep_copy (vertex_colors_, h_colors);
    delete [] forbidden;
  }

public:
  //Speculative Coloring Functors

  /**
   * Functor for VBBIT algorithms speculative coloring with edge filtering.
   */
  struct functorGreedyColor_IMPLOG_EF {
    nnz_lno_t nv;
    const_lno_row_view_t _idx; //rowmap
    nnz_lno_temp_work_view_t _adj; //entries
    color_view_type _colors; // vertex colors
    nnz_lno_temp_work_view_t _vertexList; // conflictlist
    nnz_lno_t _vertexListLength;
    nnz_lno_t _chunkSize;

    functorGreedyColor_IMPLOG_EF(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
		nnz_lno_temp_work_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_t vertexListLength,
        nnz_lno_t chunkSize) :nv(nv_), _idx(xadj_), _adj(adj_), _colors(colors),
      _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii) const {
      nnz_lno_t i = 0;

      //outer loop is on chunks, a thread is assigned as many vertex as the chunksize
      for (nnz_lno_t ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength){
          i = _vertexList(ii*_chunkSize +ichunk);
        }
        else{
          continue;
        }

        if (_colors(i) > 0) continue; // Already colored this vertex

        size_type my_xadj_end = _idx(i+1);
        size_type xadjbegin = _idx(i);

        // Do multiple passes if array is too small.
        color_t degree = my_xadj_end-xadjbegin; // My degree
        color_t offset = 0;

        //we parse the neigborlist multiple times,
        //each time we look for a certain range of colors.
        for (; (offset <= degree + VBBIT_COLORING_FORBIDDEN_SIZE); offset += VBBIT_COLORING_FORBIDDEN_SIZE){

          // Forbidden colors
          //we use a single (long) int for forbidden colors
          ban_type forbidden = 0;

          // Check nbors, fill forbidden array.
          for (size_type j=xadjbegin; j<my_xadj_end; ++j){
            nnz_lno_t n  = _adj(j);
            if (n == i || n >= nv) continue; // Skip self-loops
            color_t c = _colors(n);

            color_t color_offset = c-offset;
            //if color is within the current range, or if its color is in a previously traversed range
            if (c && color_offset <= VBBIT_COLORING_FORBIDDEN_SIZE){
              //apply edge filtering, place it to front of the adjacency list,
              //so that we wont see that anymore.
              if (j > xadjbegin){
                _adj(j) = _adj(xadjbegin);
                _adj(xadjbegin) = n;
              }
              ++xadjbegin;

              //if it is in the current range, then add the color to banned colors
              if (c > offset){
                //convert color to bit representation.
                ban_type ban_color_bit = 1;
                ban_color_bit = ban_color_bit << (color_offset - 1);
                //add it to forbidden colors
                forbidden = forbidden | ban_color_bit;
                //if there are no available colors in this range,
                //early exit, no need to traverse the rest.
                if (~forbidden == 0) {
                  break;
                }
              }
            }
          }

          forbidden = ~(forbidden);
          //check if an available color exits.
          if (forbidden){
            //if there is an available color, choose the first color,
            //using 2s complement.
            ban_type my_new_color = forbidden & (-forbidden);
            color_t val = 1;
            //convert it back to decimal color.
            while ((my_new_color & 1) == 0) {
              ++val;
              my_new_color = my_new_color >> 1;
            }
            _colors(i) = val + offset;
            break;
          }
        }
      }
    }
  };

  /**
   * Functor for VBBIT algorithms speculative coloring without edge filtering.
   */
  struct functorGreedyColor_IMPLOG {

    nnz_lno_t nv;
    const_lno_row_view_t _idx;
    const_lno_nnz_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_t _vertexListLength;
    nnz_lno_t _chunkSize;

    functorGreedyColor_IMPLOG(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
        const_lno_nnz_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_t vertexListLength,
        nnz_lno_t chunkSize) : nv(nv_),
          _idx(xadj_), _adj(adj_), _colors(colors),
          _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii) const {
      nnz_lno_t i = 0;
      for (nnz_lno_t ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex

        size_type my_xadj_end = _idx(i+1);
        size_type xadjbegin = _idx(i);
        // Do multiple passes if array is too small.
        color_t degree = my_xadj_end - xadjbegin; // My degree
        color_t offset = 0;

        for (; (offset <= degree + VBBIT_COLORING_FORBIDDEN_SIZE); offset += VBBIT_COLORING_FORBIDDEN_SIZE){

          ban_type forbidden = 0; // Forbidden colors

          // Check nbors, fill forbidden array.
          for (size_type j=xadjbegin; j<my_xadj_end; ++j){
            nnz_lno_t n  = _adj(j);
            if (n == i || n >= nv) continue; // Skip self-loops
            color_t c = _colors(n);
            color_t color_offset = c-offset;
            //if color is in the current range
            //convert it to binary and add it to forbidden
            if (color_offset <= VBBIT_COLORING_FORBIDDEN_SIZE && c > offset){

              ban_type ban_color_bit = 1;
              ban_color_bit = ban_color_bit << (color_offset - 1);

              forbidden = forbidden | ban_color_bit;
              if (~forbidden == 0) {
                break;
              }
            }
          }
          forbidden = (~forbidden);

          if (forbidden){
            ban_type my_new_color = forbidden & (-forbidden);

            color_t val = 1;

            while ((my_new_color & 1) == 0) {
              ++val;
              my_new_color = my_new_color >> 1;
            }
            _colors(i) = val + offset;
            break;
          }
        }
      }
    }

  };

  /**
   * Functor for VBCS algorithms speculative coloring with edge filtering.
   */
  struct functorGreedyColor_IMP_EF {

    nnz_lno_t nv;
    const_lno_row_view_t _xadj;
    nnz_lno_temp_work_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _color_set ;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_t _vertexListLength;
    nnz_lno_t _chunkSize;

    functorGreedyColor_IMP_EF(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
        nnz_lno_temp_work_view_t adj_,
        color_view_type colors, nnz_lno_temp_work_view_t color_set,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_t vertexListLength,
        nnz_lno_t chunkSize): nv(nv_),
        _xadj(xadj_), _adj(adj_),
        _colors(colors), _color_set(color_set),
        _vertexList(vertexList), _vertexListLength(vertexListLength),
        _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t &ii) const {
      nnz_lno_t i = 0;
      for (nnz_lno_t ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex
        size_type xadj_end = _xadj(i+1);
        size_type xadj_begin = _xadj(i);

        //my color set starts from zero, but if we are leaving vertices
        //that cannot be colored in this iteration, we retrieve it from their previous color_sets.
        nnz_lno_t my_color_set = 0;
        while (1){
          color_t ban_colors = 0;

          for (size_type j = xadj_begin; j < xadj_end && ~ban_colors; ++j){
            nnz_lno_t n = _adj(j);
            if (n == i|| n >= nv) continue; // Skip self-loops

            nnz_lno_t neighbor_color_set = _color_set(n);
            //only if the neigbor has the same color set
            if (neighbor_color_set <= my_color_set ){
              color_t ncolor = _colors(n);
              if (ncolor){
                if(j > xadj_begin){
                  _adj(j) = _adj(xadj_begin);
                  _adj(xadj_begin) = n;
                }
                ++xadj_begin;
                if (neighbor_color_set == my_color_set){
                  ban_colors = ban_colors | ncolor;
                }
              }
            }
          }

          ban_colors = ~(ban_colors);
          if (ban_colors){
            color_t my_color = ban_colors & (-ban_colors);
            _color_set(i) = my_color_set;
            _colors(i) = my_color;
            break;
          }
          else{
            my_color_set += 1;
          }
        }
      }
    }
  };

  /**
   * Functor for VBCS algorithms speculative coloring without edge filtering.
   */
  struct functorGreedyColor_IMP {

    nnz_lno_t nv;
    const_lno_row_view_t _xadj;
    const_lno_nnz_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _color_set ;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_t _vertexListLength;
    nnz_lno_t _chunkSize;


    functorGreedyColor_IMP(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
        const_lno_nnz_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t color_set,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_t vertexListLength,
        nnz_lno_t chunkSize
    ) : nv (nv_),
      _xadj(xadj_), _adj(adj_),
      _colors(colors), _color_set(color_set),
      _vertexList(vertexList), _vertexListLength(vertexListLength),
      _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t &ii) const {
      nnz_lno_t i = 0;
      for (nnz_lno_t ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex
        size_type xadj_end = _xadj(i+1);
        size_type xadj_begin = _xadj(i);

        //my color set starts from zero, but if we are leaving vertices
        //that cannot be colored in this iteration, we retrieve it from their previous color_sets.
        nnz_lno_t my_color_set =  0;
        //idx degree = xadj_end - xadj_begin;
        for (; ;){
          color_t ban_colors = 0;
          for (size_type j = xadj_begin; j < xadj_end ;++j){
            nnz_lno_t n = _adj(j);
            if (n == i|| n >= nv) continue; // Skip self-loops
            if ( my_color_set == _color_set(n)){
              ban_colors = ban_colors | _colors(n);
              if (~ban_colors == 0){
                break;
              }
            }
          }

          ban_colors = ~(ban_colors);

          if (ban_colors){
            color_t my_color = ban_colors & (-ban_colors);
            _color_set(i) = my_color_set;
            _colors(i) = my_color;
            break;
          }
          else{
            my_color_set += 1;
          }
        }
      }
    }
  };

  /**
   * Functor for VB algorithm speculative coloring with edge filtering.
   */
  struct functorGreedyColor_EF {
    nnz_lno_t nv;
    const_lno_row_view_t _idx;
    nnz_lno_temp_work_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_t _vertexListLength;
    nnz_lno_t _chunkSize;

    functorGreedyColor_EF(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
		nnz_lno_temp_work_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_t vertexListLength,
        nnz_lno_t chunkSize
    ) : nv (nv_),
      _idx(xadj_), _adj(adj_), _colors(colors),
      _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii) const {
      // Color vertex i with smallest available color.
      //
      // Each thread colors a chunk of vertices to prevent all
      // vertices getting the same color.
      //
      // This version uses a bool array of size FORBIDDEN_SIZE.
      // TODO: With chunks, the forbidden array should be char/int
      //       and reused for all vertices in the chunk.
      //
//      std::cout << ">>> WCMCLEN functorGreedyColor_EF (KokkosGraph_GraphColor_impl.hpp) <<<" << std::endl;
      nnz_lno_t i = 0;
      for (nnz_lno_t ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex

        bool foundColor = false; // Have we found a valid color?

        // Use forbidden array to find available color.
        // This array should be small enough to fit in fast memory (use Kokkos memoryspace?)
        bool forbidden[VB_COLORING_FORBIDDEN_SIZE]; // Forbidden colors

        // Do multiple passes if array is too small.
        color_t degree = _idx(i+1)-_idx(i); // My degree
        size_type my_xadj_end = _idx(i+1);
        color_t offset = 0;
        size_type xadjbegin = _idx(i);

        for (; (offset <= degree + VB_COLORING_FORBIDDEN_SIZE) && (!foundColor); offset += VB_COLORING_FORBIDDEN_SIZE){
          // initialize
          for (int j=0; j< VB_COLORING_FORBIDDEN_SIZE; j++){
            forbidden[j] = false;
          }
          if (offset == 0) forbidden[0] = true; // by convention, start at 1

          // Check nbors, fill forbidden array.
          for (size_type j=xadjbegin; j<my_xadj_end; ++j){
            nnz_lno_t n  = _adj(j);
            if (n == i|| n >= nv) {
              continue; // Skip self-loops
            }
            color_t c= _colors(n);
            // Removed option to leave potentially conflicted vertices uncolored.
            //if (c== -1){ // Nbor is being colored at same time
            //  _colors[i] = 0; // Neutral color, skip and recolor later
            //  foundColor = true;
            //  return;
            //}
            if ((c>= offset) && (c-offset < VB_COLORING_FORBIDDEN_SIZE)){
              forbidden[c-offset] = true;
            }
            if (c && c-offset < VB_COLORING_FORBIDDEN_SIZE){
              if (j > xadjbegin){
                _adj(j) = _adj(xadjbegin);
                _adj(xadjbegin) = n;
              }
              ++xadjbegin;
            }

          }

          // color vertex i with smallest available color (FirstFit)
          // TODO: Add options for other color choices (Random, LeastUsed)
          for (int c=0; c< VB_COLORING_FORBIDDEN_SIZE; c++){
            if (!forbidden[c]){
              _colors(i) = offset+c;
              //_colors[i] += (i&1); // RandX strategy to reduce conflicts
              foundColor = true;
              break;
            }
          }
        }
      }
    }

  };


  /**
   * Functor for VB algorithm speculative coloring without edge filtering.
   */
  struct functorGreedyColor {
    nnz_lno_t nv;
    const_lno_row_view_t _idx;
    const_lno_nnz_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_t _vertexListLength;
    nnz_lno_t _chunkSize;

    functorGreedyColor(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
        const_lno_nnz_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_t vertexListLength,
        nnz_lno_t chunkSize
    ) : nv (nv_),
      _idx(xadj_), _adj(adj_), _colors(colors),
      _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii) const {
      // Color vertex i with smallest available color.
      //
      // Each thread colors a chunk of vertices to prevent all
      // vertices getting the same color.
      //
      // This version uses a bool array of size FORBIDDEN_SIZE.
      // TODO: With chunks, the forbidden array should be char/int
      //       and reused for all vertices in the chunk.
      //
//      std::cout << ">>> WCMCLEN functorGreedyColor (KokkosGraph_GraphColor_impl.hpp)" << std::endl;

      nnz_lno_t i = 0;
      for (nnz_lno_t ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex

        bool foundColor = false; // Have we found a valid color?

        // Use forbidden array to find available color.
        // This array should be small enough to fit in fast memory (use Kokkos memoryspace?)
        bool forbidden[VB_COLORING_FORBIDDEN_SIZE]; // Forbidden colors

        // Do multiple passes if array is too small.
        color_t degree = _idx(i+1)-_idx(i); // My degree
        color_t offset = 0;
        for (; (offset <= degree + VB_COLORING_FORBIDDEN_SIZE) && (!foundColor); offset += VB_COLORING_FORBIDDEN_SIZE){
          // initialize
          for (int j=0; j< VB_COLORING_FORBIDDEN_SIZE; j++){
            forbidden[j] = false;
          }
          if (offset == 0) forbidden[0] = true; // by convention, start at 1

          // Check nbors, fill forbidden array.
          for (size_type j=_idx(i); j<_idx(i+1); j++){
            if (_adj(j) == i|| _adj(j)  >= nv) continue; // Skip self-loops
            color_t c= _colors(_adj(j));
            // Removed option to leave potentially conflicted vertices uncolored.
            //if (c== -1){ // Nbor is being colored at same time
            //  _colors[i] = 0; // Neutral color, skip and recolor later
            //  foundColor = true;
            //  return;
            //}
            if ((c>= offset) && (c-offset < VB_COLORING_FORBIDDEN_SIZE))
              forbidden[c-offset] = true;
          }

          // color vertex i with smallest available color (FirstFit)
          // TODO: Add options for other color choices (Random, LeastUsed)
          for (int c=0; c< VB_COLORING_FORBIDDEN_SIZE; c++){
            if (!forbidden[c]){
              _colors(i) = offset+c;
              //_colors[i] += (i&1); // RandX strategy to reduce conflicts
              foundColor = true;
              break;
            }
          }
        }
      }
    }
  };      // functorGreedyColor  (end)


  //Conflict find and worklist creation functors.

  /**
   * Finds conflicts without creating a new worklist
   */
  template <typename adj_view_t>
  struct functorFindConflicts_No_Conflist {

    nnz_lno_t nv;
    const_lno_row_view_t _idx;
    adj_view_t _adj;
    color_view_type _colors;

    functorFindConflicts_No_Conflist(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
		adj_view_t adj_,
        color_view_type colors) : nv (nv_),
          _idx(xadj_), _adj(adj_),_colors(colors)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii, nnz_lno_t &numConflicts) const {

      color_t my_color = _colors(ii);
      size_type xadjend = _idx(ii+1);
      size_type j=_idx(ii);
#ifdef DEGREECOMP
      idx myDegree = xadjend - j;
#endif
      for (; j< xadjend; j++){
        nnz_lno_t neighbor = _adj(j);

        if (
#ifndef DEGREECOMP
            ii < neighbor && neighbor < nv &&
#endif
            _colors(neighbor) == my_color
#ifdef DEGREECOMP
            && (myDegree < _idx(neighbor + 1) - _idx(neighbor) ||
                (myDegree == _idx(neighbor + 1) - _idx(neighbor) && ii < neighbor))
#endif
        ) {
          //std::cout << "me:" << ii << " n:" << neighbor << " color:" << my_color << std::endl;
          _colors(ii) = 0; // Uncolor vertex i
          numConflicts += 1;
          break; // Once i is uncolored and marked conflict
        }
      }
    }
  };


  /**
   * Finds conflicts by marking the work vertices to be used later for creation of new worklist with PPS
   */
  template <typename adj_view_t>
  struct functorFindConflicts_PPS {
    nnz_lno_t nv;
    const_lno_row_view_t _idx;
    adj_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_temp_work_view_t _recolorList;



    functorFindConflicts_PPS(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
		adj_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_temp_work_view_t recolorList) :
          nv (nv_),
          _idx(xadj_), _adj(adj_), _colors(colors),
          _vertexList(vertexList),
          _recolorList(recolorList){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii, nnz_lno_t &numConflicts) const {
      nnz_lno_t i = _vertexList(ii);
      color_t my_color = _colors(i);
      _recolorList(i) = 0;
      // check vertex i conflicts

      size_type xadjend = _idx(i+1);
      size_type j=_idx(i);
#ifdef DEGREECOMP
      idx myDegree = xadjend - j;
#endif
      for (; j<xadjend; j++){
        nnz_lno_t neighbor = _adj(j);
        if (
#ifndef DEGREECOMP
            i < neighbor && neighbor < nv &&
#endif
            _colors(neighbor) == my_color
#ifdef DEGREECOMP
            && (myDegree < _idx(neighbor + 1) - _idx(neighbor)
                                                     || (myDegree == _idx(neighbor + 1) - _idx(neighbor) && i < neighbor))
#endif
        ) {
          _colors(i) = 0; // Uncolor vertex i
          _recolorList(i) = 1;
          numConflicts += 1;
          break; // Once i is uncolored and marked conflict
        }
      }
    }
  };


  /**
   * Finds conflicts and creates new worklist using atomic operations.
   */
  template <typename adj_view_t>
  struct functorFindConflicts_Atomic {
    nnz_lno_t nv;
    const_lno_row_view_t _idx;
    adj_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_temp_work_view_t _recolorList;
    single_dim_index_view_type _recolorListLength;


    functorFindConflicts_Atomic(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
        adj_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_temp_work_view_t recolorList,
        single_dim_index_view_type recolorListLength
    ) : nv (nv_),
      _idx(xadj_), _adj(adj_), _colors(colors),
      _vertexList(vertexList),
      _recolorList(recolorList),
      _recolorListLength(recolorListLength){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii, nnz_lno_t &numConflicts) const {
      typedef typename std::remove_reference< decltype( _recolorListLength() ) >::type atomic_incr_type;

      nnz_lno_t i = _vertexList(ii);
      color_t my_color = _colors(i);

      size_type xadjend = _idx(i+1);
      size_type j=_idx(i);
#ifdef DEGREECOMP
      idx myDegree = xadjend - j;
#endif

      for (; j < xadjend; j++){
        nnz_lno_t neighbor = _adj(j);
        if (
#ifndef DEGREECOMP
            i < neighbor && neighbor < nv &&
#endif
            _colors(neighbor) == my_color
#ifdef DEGREECOMP
            && (myDegree < _idx(neighbor + 1) - _idx(neighbor) ||
                (myDegree == _idx(neighbor + 1) - _idx(neighbor) && i < neighbor))
#endif
        ) {
          _colors(i) = 0; // Uncolor vertex i
          // Atomically add vertex i to recolorList
          const nnz_lno_t k = Kokkos::atomic_fetch_add( &_recolorListLength(), atomic_incr_type(1));
          _recolorList(k) = i;
          numConflicts += 1;
          break; // Once i is uncolored and marked conflict
        }
      }
    }
  };    // struct functorFindConflicts_Atomic (end)


  /**
   * VBCS:  Finds conflicts without creating a new worklist
   */
  template <typename adj_view_t>

  struct functorFindConflicts_No_Conflist_IMP {

    nnz_lno_t nv;
    const_lno_row_view_t _xadj;
    adj_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _color_sets;


    functorFindConflicts_No_Conflist_IMP(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
		adj_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t color_sets
    ) : nv (nv_),
      _xadj(xadj_), _adj(adj_), _colors(colors), _color_sets(color_sets){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii, nnz_lno_t &numConflicts) const {
      color_t my_color = _colors(ii);
      if (my_color == 0){
        // this should only happen when one_color_set_per_iteration is set to true.
        numConflicts++;
      }
      else {
        nnz_lno_t my_color_set = _color_sets(ii);
        size_type my_xadj_end = _xadj(ii+1);
        // check vertex i conflicts

        size_type j=_xadj(ii);
#ifdef DEGREECOMP
        idx myDegree = my_xadj_end - j;
#endif

        for (; j<my_xadj_end; j++){
          nnz_lno_t neighbor = _adj(j);
          if (
#ifndef DEGREECOMP
              ii < neighbor && neighbor < nv &&
#endif
              _colors(neighbor) == my_color && my_color_set == _color_sets(neighbor)
#ifdef DEGREECOMP
              && (myDegree < _xadj(neighbor + 1) - _xadj(neighbor)||
                  (myDegree == _xadj(neighbor + 1) - _xadj(neighbor) && ii < neighbor))
#endif
          ) {
            _colors(ii) = 0; // Uncolor vertex i
            _color_sets(ii) = 0;
            numConflicts++;
            break; // Once i is uncolored and marked conflict
          }
        }

      }
    }
  };    // functorFindConflicts_No_Conflist_IMP (end)


  /**
   * VBCS: Finds conflicts by marking the work vertices to be used later for creation of new worklist with PPS
   */
  template <typename adj_view_t>
  struct functorFindConflicts_PPS_IMP {
    nnz_lno_t nv;
    const_lno_row_view_t _xadj;
    adj_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _color_sets;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_temp_work_view_t _recolorList;

    functorFindConflicts_PPS_IMP(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
		adj_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t color_sets,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_temp_work_view_t recolorList
    ) : nv (nv_),
      _xadj(xadj_), _adj(adj_), _colors(colors), _color_sets(color_sets),
      _vertexList(vertexList),
      _recolorList(recolorList){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii, nnz_lno_t &numConflicts) const {
      nnz_lno_t i = _vertexList(ii);
      _recolorList(i) = 0;
      color_t my_color = _colors(i);
      if (my_color == 0){
        _recolorList(i) = 1;
        numConflicts++;
      }
      else {
        nnz_lno_t my_color_set = _color_sets(i);
        size_type my_xadj_end = _xadj(i+1);
        // check vertex i conflicts

        size_type j=_xadj(i);
#ifdef DEGREECOMP
        idx myDegree = my_xadj_end - j;
#endif
        for (; j<my_xadj_end; j++){
          nnz_lno_t neighbor = _adj(j);
          if (
#ifndef DEGREECOMP
              i < neighbor && neighbor < nv &&
#endif
              _colors(neighbor) == my_color && my_color_set == _color_sets(neighbor)
#ifdef DEGREECOMP
              && (myDegree < _xadj(neighbor + 1) - _xadj(neighbor)||
                  (myDegree == _xadj(neighbor + 1) - _xadj(neighbor) && i < neighbor))
#endif
          ) {
            _colors(i) = 0; // Uncolor vertex i
            _color_sets(i) = 0;
            _recolorList(i) = 1;
            numConflicts++;
            break; // Once i is uncolored and marked conflict
          }
        }
      }
    }
  };    // functorFindConflicts_PPS_IMP (end)


  /**
   * VBCS:Finds conflicts and creates new worklist using atomic operations.
   */
  template <typename adj_view_t>
  struct functorFindConflicts_Atomic_IMP {

    nnz_lno_t nv;
    const_lno_row_view_t _xadj;
    adj_view_t _adj;
    color_view_type _colors;
    nnz_lno_temp_work_view_t _color_sets;
    nnz_lno_temp_work_view_t _vertexList;
    nnz_lno_temp_work_view_t _recolorList;
    single_dim_index_view_type _recolorListLength;


    functorFindConflicts_Atomic_IMP(
        nnz_lno_t nv_,
        const_lno_row_view_t xadj_,
		adj_view_t adj_,
        color_view_type colors,
        nnz_lno_temp_work_view_t color_sets,
        nnz_lno_temp_work_view_t vertexList,
        nnz_lno_temp_work_view_t recolorList,
        single_dim_index_view_type recolorListLength
    ) : nv (nv_),
      _xadj(xadj_), _adj(adj_), _colors(colors), _color_sets(color_sets),
      _vertexList(vertexList),
      _recolorList(recolorList),
      _recolorListLength(recolorListLength){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii, nnz_lno_t &numConflicts) const {
      typedef typename std::remove_reference< decltype( _recolorListLength() ) >::type atomic_incr_type;
      nnz_lno_t i = _vertexList(ii);
      color_t my_color = _colors(i);
      if (my_color == 0){
        // this should only happen when one_color_set_per_iteration is set to true.
        const nnz_lno_t k = Kokkos::atomic_fetch_add( &_recolorListLength(), atomic_incr_type(1));
        _recolorList(k) = i;
        numConflicts++;
      }
      else {
        nnz_lno_t my_color_set = _color_sets(i);
        size_type my_xadj_end = _xadj(i+1);
        // check vertex i conflicts

        size_type j=_xadj(i);
#ifdef DEGREECOMP
        idx myDegree = my_xadj_end - j;
#endif
        for (; j< my_xadj_end; j++){
          nnz_lno_t neighbor = _adj(j);
          if (
#ifndef DEGREECOMP
              i < neighbor && neighbor < nv &&
#endif
              _colors(neighbor) == my_color && my_color_set == _color_sets(neighbor)
#ifdef DEGREECOMP
              && (myDegree < _xadj(neighbor + 1) - _xadj(neighbor)||
                  (myDegree == _xadj(neighbor + 1) - _xadj(neighbor) && i < neighbor))
#endif
          ) {
            _colors(i) = 0; // Uncolor vertex i
            _color_sets(i) = 0;
            // Atomically add vertex i to recolorList
            const nnz_lno_t k = Kokkos::atomic_fetch_add( &_recolorListLength(), atomic_incr_type(1));
            _recolorList(k) = i;
            numConflicts++;
            break; // Once i is uncolored and marked conflict
          }
        }
      }
    }
  };      // functorFindConflicts_Atomic_IMP (end)

  //Helper Functors
  /**
   * Functor to init a list sequentialy, that is list[i] = i
   */
  template <typename view_type>
  struct functorInitList{
    view_type _vertexList;
    functorInitList (view_type vertexList) : _vertexList(vertexList) { }
    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t i) const {
      // Natural order
      _vertexList(i) = i;
    }
  };


  /**
   * Functor for parallel prefix sum
   */
  template <typename view_type>
  struct parallel_prefix_sum{
    view_type _vertexList;
    view_type _recolorList;
    view_type _pps_view;

    parallel_prefix_sum(
        view_type vertexList,
        view_type recolorList,
        view_type pps_view):
          _vertexList(vertexList),_recolorList(recolorList),_pps_view(pps_view){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const typename view_type::non_const_value_type ii, size_t& update, const bool final) const {
      typename view_type::non_const_value_type w = _vertexList(ii);
      update += _recolorList(w);
      if (final) {
        _pps_view(w) = (update);
      }
    }
  };


  /**
   * Functor for creating new worklist using pps
   */
  template <typename view_type>
  struct create_new_work_array{
    view_type _vertexList;
    view_type _recolorList;
    view_type _pps_view;

    create_new_work_array(
        view_type vertexList,
        view_type recolorList,
        view_type pps_view):
          _vertexList(vertexList),_recolorList(recolorList),_pps_view(pps_view){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const typename view_type::non_const_value_type ii) const {
      typename view_type::non_const_value_type w = _vertexList(ii);
      typename view_type::non_const_value_type left_work = 0;
      if (ii > 0){
        left_work = _pps_view(_vertexList(ii - 1));
      }
      typename view_type::non_const_value_type pps_current = _pps_view(w);
      if(pps_current != left_work){
        typename view_type::non_const_value_type future_index = pps_current;
        _recolorList(future_index - 1) = w;
      }
    }
  };


  /**
   * Converting VBCS colors to final colors.
   */
  struct set_final_colors{
    color_view_type kokcol;
    nnz_lno_temp_work_view_t kokcolset; //the colors that are represented with bits, and the colors set that the color is in.
    color_t color_size;

    /** \brief functor constructor.
     * \param kokcol_  the colors of the vertices. Represented with bits.
     * \param kokcolset_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     */
    set_final_colors(color_view_type kokcol_, nnz_lno_temp_work_view_t kokcolset_
    ): kokcol(kokcol_),kokcolset(kokcolset_), color_size ( sizeof(color_t) * 8){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t &ii) const {

      color_t val = kokcol(ii);
      if (val){
        //find the position in the bit.
        nnz_lno_t i = 1;
        while ((val & 1) == 0) {
          ++i;
          val = val >> 1;
        }

        //idx i = log2(val) + 1;
        //set the final color.
        kokcol(ii) = i + kokcolset(ii) * color_size;
      }
    }
  };
};  // class GraphColor_VB




/*! \brief Class for the deterministic vertex based graph coloring algorithms.
 */
template <typename HandleType, typename lno_row_view_t_, typename lno_nnz_view_t_>
class GraphColor_VBD:public GraphColor <HandleType,lno_row_view_t_,lno_nnz_view_t_>{
public:

  typedef long long int ban_type;

  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef typename HandleType::color_view_t color_view_type;


  typedef typename HandleType::size_type size_type;
  typedef typename lno_row_view_t_::device_type row_lno_view_device_t;

  typedef typename HandleType::nnz_lno_t nnz_lno_t;

  typedef typename HandleType::color_t color_t;
  typedef typename HandleType::color_host_view_t color_host_view_t; //Host view type

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename Kokkos::View<nnz_lno_t, row_lno_view_device_t> single_dim_index_view_type;
  typedef typename single_dim_index_view_type::HostMirror single_dim_index_host_view_type; //Host view type

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  typedef typename HandleType::size_type_temp_work_view_t size_type_temp_work_view_t;
  typedef typename HandleType::size_type_persistent_work_view_t size_type_persistent_work_view_t;


  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;



  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;


  typedef typename lno_nnz_view_t_::const_type const_lno_nnz_view_t;
  typedef typename lno_nnz_view_t_::non_const_type non_const_lno_nnz_view_t;


protected:



  bool _ticToc; //if true print info in each step
  int _chunkSize; //the size of the minimum work unit assigned to threads. Changes the convergence on GPUs
  char _use_color_set; //the VBD algorithm type.
                        // 0 for VBD:

public:
  /**
   * \brief GraphColor_VBD constructor.
   * \param nv_: number of vertices in the graph
   * \param ne_: number of edges in the graph
   * \param row_map: the xadj array of the graph. Its size is nv_ +1
   * \param entries: adjacency array of the graph. Its size is ne_
   * \param coloring_handle: GraphColoringHandle object that holds the specification about the graph coloring,
   *    including parameters.
   */
  GraphColor_VBD(
      nnz_lno_t nv_, size_type ne_,
      const_lno_row_view_t row_map, const_lno_nnz_view_t entries,
      HandleType *coloring_handle):
    GraphColor<HandleType,lno_row_view_t_,lno_nnz_view_t_>(nv_, ne_, row_map, entries, coloring_handle),
    _ticToc(coloring_handle->get_tictoc()),
    _chunkSize(coloring_handle->get_vb_chunk_size()),
    _use_color_set()
    {
      switch (coloring_handle->get_coloring_algo_type()){
      case COLORING_VBD:
        this->_use_color_set = 0;
        break;
      case COLORING_VBDBIT:
        this->_use_color_set = 1;
        break;
      default: //cannnot get in here.
        this->_use_color_set = 0;
        break;

      }

    }

  /** \brief GraphColor_VBD destructor.
    */
  virtual ~GraphColor_VBD(){}

  /** \brief Function to color the vertices of the graphs. Performs a vertex-based coloring.
   * \param colors is the output array recording the color of each vertex. Size is this->nv.
   *   Attn: Color array must contain only positive numbers. If there are no initial colors,
   *   it should be all initialized with zeros. Any strictly positive value in the color array,
   *   will make the algorithm assume that the corresponding vertex is already .
   * \param num_loops: The number of loops in the while statement required to color the graph.
   */
  virtual void color_graph(color_view_type colors, int &num_loops){

    if (this->_ticToc) {
      std::cout
          << "\tVBD params:" << std::endl
          << "\talgorithm:" << (int)this->_use_color_set << std::endl
          << "\tticToc:" << (int) this->_ticToc << std::endl
          << "\tchunkSize:" << this->_chunkSize << std::endl;
    }

    nnz_lno_t numVertices = this->nv;

    size_type maxColors = 0;
    nnz_lno_persistent_work_view_t score
      = nnz_lno_persistent_work_view_t(Kokkos::ViewAllocateWithoutInitializing("score"), this->nv);
    functorScoreCalculation<size_type, MyExecSpace> scoreCalculation(score, this->xadj);
    
    Kokkos::parallel_reduce("Deterministic Coloring: compute initial scores", my_exec_space(0, this->nv),
                            scoreCalculation, Kokkos::Max<size_type>(maxColors));

   if (this->_ticToc) {
     std::cout << "maxColors: " << maxColors << std::endl;
   }

    // Create the dependency list of the graph
    nnz_lno_persistent_work_view_t dependency("dependency", numVertices);
    Kokkos::View<size_type, MyTempMemorySpace> frontierSize("frontierSize");
    typename Kokkos::View<size_type, MyTempMemorySpace>::HostMirror host_frontierSize =
      Kokkos::create_mirror_view(frontierSize);
    Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize("newFrontierSize");
    typename Kokkos::View<size_type, MyTempMemorySpace>::HostMirror host_newFrontierSize =
      Kokkos::create_mirror_view(newFrontierSize);
    nnz_lno_temp_work_view_t frontier("frontier", numVertices);
    nnz_lno_temp_work_view_t newFrontier("newFrontier", numVertices);
    functorInitialDependency myInitialDependency(this->xadj, this->adj, score, dependency,
                                                 newFrontier, newFrontierSize);
    Kokkos::parallel_for("Deterministic Coloring: compute dependency list",
                         my_exec_space(0, numVertices),
                         myInitialDependency);

    Kokkos::deep_copy(host_newFrontierSize, newFrontierSize);
    while(host_newFrontierSize() > 0) {
      ++num_loops;
      // First swap fontier with newFrontier and fontierSize with newFrontierSize
      // reset newFrontierSize
      functorSwapOnDevice mySwapOnDevice(frontierSize, newFrontierSize);
      Kokkos::parallel_for("Swap frontier sizes", my_exec_space(0, 1), mySwapOnDevice);
      Kokkos::deep_copy(host_frontierSize, frontierSize);
      {
	auto swap_tmp = frontier;
	frontier = newFrontier;
	newFrontier = swap_tmp;
      }

      // Loop over nodes in the frontier
      // First variant without bit array, easier to understand/program
      if(this->_use_color_set == 0) {
        functorDeterministicColoring myDeterministicColoring(this->xadj, this->adj,
                                                             dependency, frontier, frontierSize,
                                                             newFrontier, newFrontierSize,
                                                             maxColors, colors);
        Kokkos::parallel_for("Deterministic Coloring: color nodes in frontier",
                             my_exec_space(0, host_frontierSize()),
                             myDeterministicColoring);

      } else if(this->_use_color_set == 1) {
        // Second variant with bit array for efficiency on GPU
        // The bit array is of size 64 so if maxColors > 64,
        // we need to use successive color ranges of width 64
        // to represent all the possible colors on the graph.
        functorDeterministicColoringBitArray myDeterministicColoringBitArray(this->xadj, this->adj,
                                                                             dependency, frontier,
                                                                             frontierSize,
                                                                             newFrontier,
                                                                             newFrontierSize,
                                                                             maxColors, colors);
        Kokkos::parallel_for("Deterministic Coloring: color nodes in frontier",
                             my_exec_space(0, host_frontierSize()),
                             myDeterministicColoringBitArray); // Loop over current frontier
      }
      Kokkos::deep_copy(host_newFrontierSize, newFrontierSize);
    } // while newFrontierSize

  } // color_graph()


  template <class max_type, class execution_space>
  struct functorScoreCalculation {
    nnz_lno_persistent_work_view_t score_;
    const_lno_row_view_t numNeighbors_;

    functorScoreCalculation(nnz_lno_persistent_work_view_t& score, const_lno_row_view_t& numNeighbors)
      : score_(score), numNeighbors_(numNeighbors) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int i,  size_type &update) const {
      score_(i) = numNeighbors_(i + 1) - numNeighbors_(i);
      update = ( (size_type) score_(i) < update ? update : (size_type) score_(i) );
    }
  }; // functorScoreCalculation()

  struct functorSwapOnDevice {
    Kokkos::View<size_type, MyTempMemorySpace> frontierSize_;
    Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize_;

    functorSwapOnDevice(Kokkos::View<size_type, MyTempMemorySpace> frontierSize,
                        Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize)
      : frontierSize_(frontierSize), newFrontierSize_(newFrontierSize) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int dummy) const {
      frontierSize_() = newFrontierSize_();
      newFrontierSize_() = 0;
    }

  }; // functorSwapOnDevice

  struct functorInitialDependency {
    const_lno_row_view_t xadj_;
    const_lno_nnz_view_t adj_;
    nnz_lno_persistent_work_view_t score_;
    nnz_lno_persistent_work_view_t dependency_;
    nnz_lno_temp_work_view_t newFrontier_;
    Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize_;

    functorInitialDependency(const_lno_row_view_t rowPtr,
                             const_lno_nnz_view_t colInd,
                             nnz_lno_persistent_work_view_t score,
                             nnz_lno_persistent_work_view_t dependency,
                             nnz_lno_temp_work_view_t newFrontier,
                             Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize)
      : xadj_(rowPtr), adj_(colInd), score_(score), dependency_(dependency),
        newFrontier_(newFrontier), newFrontierSize_(newFrontierSize) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int node) const {
      int myScore = score_(node);
      int numNeighs = xadj_(node + 1) - xadj_(node);
      for(int neigh = 0; neigh < numNeighs; ++neigh) {
        if(myScore < score_(adj_(xadj_(node) + neigh))) {
          dependency_(node) = dependency_(node) + 1;
        }
        if(( myScore == score_(adj_(xadj_(node) + neigh)) ) &&
           ( node < adj_(xadj_(node) + neigh) )) {
          dependency_(node) = dependency_(node) + 1;
        }
      }
      if(dependency_(node) == 0) {
        const size_type newFrontierIdx
          = Kokkos::atomic_fetch_add(&newFrontierSize_(), 1);
        newFrontier_(newFrontierIdx) = node;
      }
    }

  }; // functorInitialDependency

  struct functorDeterministicColoring {
    const_lno_row_view_t xadj_;
    const_lno_nnz_view_t adj_;
    nnz_lno_persistent_work_view_t dependency_;
    nnz_lno_temp_work_view_t frontier_;
    Kokkos::View<size_type, MyTempMemorySpace> frontierSize_;
    nnz_lno_temp_work_view_t newFrontier_;
    Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize_;
    size_type maxColors_;
    color_view_type colors_;

    functorDeterministicColoring(const_lno_row_view_t rowPtr,
                                 const_lno_nnz_view_t colInd,
                                 nnz_lno_persistent_work_view_t dependency,
                                 nnz_lno_temp_work_view_t frontier,
                                 Kokkos::View<size_type, MyTempMemorySpace> frontierSize,
                                 nnz_lno_temp_work_view_t newFrontier,
                                 Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize,
                                 size_type maxColors,
                                 color_view_type colors)
      : xadj_(rowPtr), adj_(colInd), dependency_(dependency), frontier_(frontier),
      frontierSize_(frontierSize), newFrontier_(newFrontier), newFrontierSize_(newFrontierSize),
      maxColors_(maxColors), colors_(colors) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type frontierIdx) const {
      size_type frontierNode = frontier_(frontierIdx);
      int* bannedColors = new int[maxColors_];
      for(size_type colorIdx= 0; colorIdx < maxColors_; ++colorIdx) {
        bannedColors[colorIdx] = 0;
      }

      // Loop over neighbors, find banned colors, decrement dependency and update newFrontier
      for(size_type neigh = xadj_(frontierNode); neigh < xadj_(frontierNode + 1); ++neigh) {
        bannedColors[colors_(adj_(neigh))] = 1;

        // We want to avoid the cost of atomic operations when not needed
        // so let's check that the node is not already colored, i.e.
        // its dependency is not -1.
        if(dependency_(adj_(neigh)) >= 0) {
          nnz_lno_t myDependency = Kokkos::atomic_fetch_add(&dependency_(adj_(neigh)), -1);
          // dependency(myAdj(neigh)) = dependency(myAdj(neigh)) - 1;
          if(myDependency - 1 == 0) {
            const size_type newFrontierIdx
              = Kokkos::atomic_fetch_add(&newFrontierSize_(), 1);
            newFrontier_(newFrontierIdx) = adj_(neigh);
          }
        }
      } // Loop over neighbors

      for(size_type color = 1; color < maxColors_; ++color) {
        if(bannedColors[color] == 0) {
          colors_(frontierNode) = color;
          break;
        }
      } // Loop over banned colors
      delete [] bannedColors;
    }
  }; // functorDeterministicColoring

  struct functorDeterministicColoringBitArray {
    const_lno_row_view_t xadj_;
    const_lno_nnz_view_t adj_;
    nnz_lno_persistent_work_view_t dependency_;
    nnz_lno_temp_work_view_t frontier_;
    Kokkos::View<size_type, MyTempMemorySpace> frontierSize_;
    nnz_lno_temp_work_view_t newFrontier_;
    Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize_;
    size_type maxColors_;
    color_view_type colors_;

    functorDeterministicColoringBitArray(const_lno_row_view_t rowPtr,
                                 const_lno_nnz_view_t colInd,
                                 nnz_lno_persistent_work_view_t dependency,
                                 nnz_lno_temp_work_view_t frontier,
                                 Kokkos::View<size_type, MyTempMemorySpace> frontierSize,
                                 nnz_lno_temp_work_view_t newFrontier,
                                 Kokkos::View<size_type, MyTempMemorySpace> newFrontierSize,
                                 size_type maxColors,
                                 color_view_type colors)
    : xadj_(rowPtr), adj_(colInd), dependency_(dependency), frontier_(frontier),
      frontierSize_(frontierSize), newFrontier_(newFrontier), newFrontierSize_(newFrontierSize),
      maxColors_(maxColors), colors_(colors) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type frontierIdx) const {

      size_type frontierNode = frontier_(frontierIdx);
      // Initialize bit array to all bits = 0
      unsigned long long bannedColors = 0;
      color_t myColor = 0, colorOffset = 0;

      while(myColor == 0) {
        // Loop over neighbors, find banned colors in the range:
        // [colorOffset + 1, colorOffset + 64]
        for(size_type neigh = xadj_(frontierNode); neigh < xadj_(frontierNode + 1); ++neigh) {
          color_t neighColor = colors_(adj_(neigh));
          // Check that the color is in the current range
          if(neighColor > colorOffset && neighColor < colorOffset + 65) {
            // Set bannedColors' bit in location colors(adj_(neigh)) to 1.
            bannedColors |= (1ULL << (neighColor - 1));
          }

          // We want to avoid the cost of atomic operations when not needed
          // so let's check that the node is not already colored, i.e.
          // its dependency is not -1.
          if(colorOffset == 0 && dependency_(adj_(neigh)) >= 0) {
            nnz_lno_t myDependency =
              Kokkos::atomic_fetch_add(&dependency_(adj_(neigh)), -1);
            if(myDependency - 1 == 0) {
              const size_type newFrontierIdx
                = Kokkos::atomic_fetch_add(&newFrontierSize_(), 1);
              newFrontier_(newFrontierIdx) = adj_(neigh);
            }
          }
        } // Loop over neighbors

        if(~bannedColors == 0ULL) {
          colorOffset += 64;
          // Reset bannedColors to all 0 bits
          bannedColors |= ~bannedColors;
        } else {
          color_t colorIdx = 1;
          // Check if index colordIdx - 1, is set to one in bannedColors
          while(bannedColors & (1ULL << (colorIdx - 1))) {++colorIdx;}
          myColor = colorOffset + colorIdx;
        }
      }
      colors_(frontierNode) = myColor;
    }
  }; // functorDeterministicColoringBitArray


};  // class GraphColor_VBD


/*! \brief Class for modular parallel graph coloring using Kokkos.
 *  Performs a edge_base coloring, with the hope of better load balance
 *  as well as better memory accesses on GPUs.
 */
template <typename HandleType, typename in_row_index_view_type_, typename in_nonzero_index_view_type_>
class GraphColor_EB:public GraphColor <HandleType,in_row_index_view_type_,in_nonzero_index_view_type_>{
public:

  typedef long long int ban_type;

  typedef in_row_index_view_type_ in_row_index_view_type;
  typedef in_nonzero_index_view_type_ in_nonzero_index_view_type;
  typedef typename HandleType::color_view_t color_view_type;


  typedef typename HandleType::size_type size_type;
  typedef typename in_row_index_view_type_::device_type row_lno_view_device_t;


  typedef typename HandleType::nnz_lno_t nnz_lno_t;


  typedef typename HandleType::color_t color_t;
  typedef typename HandleType::color_host_view_t color_host_view_t; //Host view type




  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename Kokkos::View<nnz_lno_t, row_lno_view_device_t> single_dim_index_view_type;

  typedef typename single_dim_index_view_type::HostMirror single_dim_index_host_view_type; //Host view type
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  typedef typename HandleType::size_type_temp_work_view_t size_type_temp_work_view_t;
  typedef typename HandleType::size_type_persistent_work_view_t size_type_persistent_work_view_t;

  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;

  typedef typename Kokkos::View<color_t *, MyTempMemorySpace> color_temp_work_view_type;

  typedef Kokkos::View<char *, MyTempMemorySpace> char_temp_work_view_type;
  typedef typename char_temp_work_view_type::HostMirror char_temp_work_host_view_type; //Host view type



  typedef typename in_row_index_view_type::const_type const_lno_row_view_t;
  typedef typename in_nonzero_index_view_type::const_type const_nonzero_index_view_type;

public:

  /**
   * \brief GraphColor_EB constructor.
   * \param nv_ number of vertices in the graph
   * \param ne_ number of edges in the graph
   * \param xadj_ the xadj array of the graph. Its size is nv_ +1
   * \param adj_ adjacency array of the graph. Its size is ne_
   */
  GraphColor_EB(nnz_lno_t nv_,
                size_type ne_,
                const_lno_row_view_t row_map,
                const_nonzero_index_view_type entries,
                HandleType *coloring_handle):
    GraphColor<HandleType,in_row_index_view_type_,in_nonzero_index_view_type_>(nv_, ne_, row_map, entries, coloring_handle)
    {}

  /**
   * \brief Class Destructor.
   */
  virtual ~GraphColor_EB(){}



  /** \brief function to color the vertices of the graphs. Performs an edge based graph coloring.
   *  the algorithm uses kokkos, so it is modular.
   * \param colors is the output array corresponding the color of each vertex.Size is this->nv.
   * \param num_loops is the output for the number of phases that the algorithm took to converge.
   */
  virtual void color_graph(color_view_type kok_colors, int &num_loops ){


//    std::cout << ">>> WCMCLEN GraphColor_EB::color_graph()" << std::endl;

    //get EB parameters
    color_t numInitialColors = this->cp->get_eb_num_initial_colors();
    double pps_cutoff = this->cp->get_min_reduction_for_conflictlist();
    size_type ps_min = this->cp->get_min_elements_for_conflictlist();
    bool use_pps = (this->cp->get_conflict_list_type() == COLORING_PPS);

    bool tictoc = this->cp->get_tictoc();

    Kokkos::Impl::Timer *timer = NULL;

    if (tictoc){
      timer = new Kokkos::Impl::Timer();
      std::cout << "\tRewriting EB params. num_initial_colors:" << numInitialColors
          << " prefix_sum_shrink_min:"  << ps_min
          << " ps_cutoff:" << pps_cutoff
          << std::endl;
    }

    size_type numEdges = 0;
    nnz_lno_persistent_work_view_t _kok_src, _kok_dst;


    this->cp->get_lower_diagonal_edge_list (this->nv, this->ne, this->xadj, this->adj, numEdges, _kok_src, _kok_dst);
    size_type num_work_edges = numEdges;

    //allocate memory for vertex ban colors, and tentative bans
    color_temp_work_view_type color_ban (Kokkos::ViewAllocateWithoutInitializing("color_ban"), this->nv);
    color_temp_work_view_type tentative_color_ban(Kokkos::ViewAllocateWithoutInitializing("tentative_color_ban"), this->nv);//views are initialized with zero
    //allocate memory for vertex color set shifts.
    nnz_lno_temp_work_view_t color_set ("color_set", this->nv); //initialized with zero.
    //initialize colors, color bans
    Kokkos::parallel_for ("KokkosGraph::GraphColoring::initColors",
        my_exec_space (0, this->nv) , init_colors (kok_colors, color_ban, numInitialColors));
    //std::cout << "nv:" << this->nv << " init_colors" << std::endl;

    //worklist
    size_type_temp_work_view_t edge_conflict_indices
    (Kokkos::ViewAllocateWithoutInitializing("edge_conflict_indices"), num_work_edges);
    //next iterations conflict list
    size_type_temp_work_view_t new_edge_conflict_indices
    (Kokkos::ViewAllocateWithoutInitializing("new_edge_conflict_indices"), num_work_edges);

    size_type_temp_work_view_t
    pps(Kokkos::ViewAllocateWithoutInitializing("prefix_sum"), num_work_edges);

    char_temp_work_view_type edge_conflict_marker
    (Kokkos::ViewAllocateWithoutInitializing("edge_conflict_marker"), num_work_edges);


    //initialize the worklist sequentiall, and markers as 1.
    Kokkos::parallel_for ("KokkosGraph::GraphColoring::InitWorkArrays",
        my_exec_space (0, num_work_edges),
        init_work_arrays(edge_conflict_indices, edge_conflict_marker)
    );
    MyExecSpace::fence();
    //std::cout << "nv:" << this->nv << " init_work_arrays" << std::endl;

    double inittime = 0;
    if (tictoc){
      inittime = timer->seconds();
      timer->reset();
    }
    double mc_time = 0, cnt_time = 0, ban_time = 0, expand_ban_time = 0, color_time = 0, pps_time = 0;

    size_type i = 0;


    if (tictoc){
    std::cout << "\t"; KokkosKernels::Impl::print_1Dview(_kok_src);
    std::cout << "\t"; KokkosKernels::Impl::print_1Dview(_kok_dst);
    std::cout << "\t"; KokkosKernels::Impl::print_1Dview(kok_colors);
    std::cout << "\t"; KokkosKernels::Impl::print_1Dview(color_set);
    }

    while(1){

      ++i;
      //std::cout << "nv:" << this->nv << " i:" << i  << " num_work_edges:" << num_work_edges<< std::endl;
      //conflict detection mark conflicts as color 0.
      //update their bans
      Kokkos::parallel_for(
          "KokkosGraph::GraphColoring::HalfEdgeMarkConflicts",
          my_exec_space(0,num_work_edges),
          halfedge_mark_conflicts (
              _kok_src, _kok_dst,
              kok_colors, color_set,
              color_ban, tentative_color_ban
              ,edge_conflict_indices
          )
      );

      MyExecSpace::fence();
      //std::cout << "nv:" << this->nv << " i:" << i <<  " halfedge_mark_conflicts" << std::endl;

      if (tictoc){
        mc_time += timer->seconds();
        timer->reset();
      }


      size_type num_conflict_reduction = 0;

      //count conflicts, and mark the edges that does not need to be processed.


      if (num_work_edges > 0)
      Kokkos::parallel_reduce(
          "KokkosGraph::GraphColoring::HalfEdgeConflictsCount",
          my_exec_space(0, num_work_edges),
          halfedge_conflict_count(
              _kok_src, _kok_dst,
              kok_colors, color_set,
              edge_conflict_indices,
              edge_conflict_marker
          ),num_conflict_reduction);

      MyExecSpace::fence();

      /*
      std::cout << "nv:" << this->nv
          << " i:" << i
          << " num_work_edges:" << num_work_edges
          << " num_conflict_reduction:" << num_conflict_reduction
          << " kok_src:" << kok_src.extent(0)
          << " kok_dst:" << kok_dst.extent(0)
          << " kok_colors:" << kok_colors.extent(0)
          << " color_set:" << color_set.extent(0)
          << " edge_conflict_indices:" << edge_conflict_indices.extent(0)
          << " edge_conflict_marker:" << edge_conflict_marker.extent(0)
          << std::endl;
      */

      if (tictoc){
        cnt_time += timer->seconds();
        timer->reset();
      }
      //std::cout << "nv:" << this->nv << " i:" << i <<  " before break:" << num_work_edges - num_conflict_reduction << std::endl;

      //if (num_work_edges <= num_conflict_reduction) break;
      if (num_work_edges - num_conflict_reduction == 0) break;

      //if the reduction is good enough w.r.t. parameters, create new worklist.
      if (num_work_edges > ps_min && num_conflict_reduction / double (num_work_edges) > pps_cutoff)
      {
        //use_pps = false;
        if (use_pps){
          //calculate new positions of the edges in new worklist
          Kokkos::parallel_scan (
              "KokkosGraph::GraphColoring::CalcEdgePositions",
              my_exec_space(0, num_work_edges),
              parallel_prefix_sum(edge_conflict_indices, edge_conflict_marker, pps)
          );
          MyExecSpace::fence();

          //write the edge indices to new worklist.
          Kokkos::parallel_for (
              "KokkosGraph::GraphColoring::CreateNewWorkArray",
              my_exec_space(0, num_work_edges),
              create_new_work_array(edge_conflict_indices, edge_conflict_marker, pps, new_edge_conflict_indices));
          MyExecSpace::fence();
        }
        else {
          //create new worklist
          single_dim_index_view_type new_index = single_dim_index_view_type("recolorListLength");;
          Kokkos::parallel_for (
              "KokkosGraph::GraphColoring::CreateNewWorkArrayAtomic",
              my_exec_space(0, num_work_edges),
              atomic_create_new_work_array(new_index, edge_conflict_indices, edge_conflict_marker, new_edge_conflict_indices));
          MyExecSpace::fence();
        }

        //swap old and new worklist
        size_type_temp_work_view_t tmp = new_edge_conflict_indices;
        new_edge_conflict_indices =  edge_conflict_indices;
        edge_conflict_indices = tmp;
        num_work_edges -= num_conflict_reduction;
        num_conflict_reduction = 0;

        if (tictoc){
          pps_time += timer->seconds();
          timer->reset();
        }
      }

      //create ban colors using the colored neighbors
      Kokkos::parallel_for (
          "KokkosGraph::GraphColoring::HalfEdgeBancColors",
          my_exec_space(0,num_work_edges),
          halfedge_ban_colors(
              _kok_src, _kok_dst,
              kok_colors, color_set,
              color_ban,edge_conflict_indices, edge_conflict_marker
          )
      );

      MyExecSpace::fence();

      if (tictoc){
        ban_time += timer->seconds();
        timer->reset();
      }


      //create tentative ban using the uncolored neighbors.
      Kokkos::parallel_for (
          "KokkosGraph::GraphColoring::HalfEdgeExpandBanForUnmatchedNeighbors",
          my_exec_space(0,num_work_edges),
          halfedge_expand_ban_for_unmatched_neighbors(
              _kok_src, _kok_dst,
              kok_colors, color_set,
              color_ban,
              tentative_color_ban,edge_conflict_indices)
      );

      if (tictoc){
        expand_ban_time += timer->seconds();
        timer->reset();
      }



      //chose a color based on the ban arrays.
      //if all colors in current set are taken, increase the color set, try again in the next iteration.
      Kokkos::parallel_for("KokkosGraph::GraphColoring::ChooseColors",
          my_exec_space(0,this->nv), choose_colors(kok_colors, color_set, color_ban, tentative_color_ban));
      if (tictoc){
        color_time += timer->seconds();
        timer->reset();
      }
    }
    if (tictoc){
      std::cout <<
          "\tinit_time:" << inittime <<
          " mc:" << mc_time <<
          " cnt_time:" << cnt_time <<
          " ban_time:" << ban_time <<
          " expand ban time:" << expand_ban_time <<
          " pps time:" << pps_time <<
          " color time:" << color_time << std::endl<< std::endl;
    }

    //set the final colors.
    Kokkos::parallel_for("KokkosGraph::GraphColoring::SetFinalColors",
        my_exec_space(0,this->nv), set_final_colors (kok_colors, color_set));

    num_loops = i;

    if (tictoc){
      delete timer;
    }
  }


  /*! \brief Functor to initialize the colors of the vertices randomly,
   *  with the hope that it will reduce the conflict in parallel execution.
   *  It also initializes the color bans.
   */
  struct init_colors{

    color_view_type kokcolors;
    color_temp_work_view_type color_ban; //colors
    color_t hash; //the number of colors to be assigned initially.

    //the value to initialize the color_ban_. We avoid using the first bit representing the sign.
    //Therefore if idx is int, it can represent 32-1 colors. Use color_set to represent more.
    color_t color_ban_init_val;


    init_colors (color_view_type colors,color_temp_work_view_type color_ban_,color_t hash_):
      kokcolors(colors), color_ban(color_ban_), hash(hash_){
      color_t tmp = 1;
      color_ban_init_val = tmp <<( sizeof(color_t) * 8 -1);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type &ii) const {
      //set colors based on their indices.
      color_t tmp1 = 1;
      kokcolors(ii) = tmp1 << (ii % hash);
      color_ban(ii) = color_ban_init_val;
    }
  };


  /*! \brief Functor to initialize the worklist
     */
  struct init_work_arrays{
    size_type_temp_work_view_t _edge_conflict_indices;
    char_temp_work_view_type _edge_conflict_marker;

    init_work_arrays (
        size_type_temp_work_view_t edge_conflict_indices,
        char_temp_work_view_type edge_conflict_marker):
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker){};

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type &ii) const {
      _edge_conflict_indices(ii)= ii; //every edge needs to be worked on initially.
      _edge_conflict_marker(ii) = 1; //every edge is a conflict initially.
    }
  };

  /**\brief Functor to mark the conflicts.
   * It goes to all the edges, and checks if two ends of an edge have the same color.
   * If they do, then it marks the one with the larger index as conflict.
   */
  struct halfedge_mark_conflicts {
    //edge list, source and destinations of the edge list.
    nnz_lno_persistent_work_view_t srcs, dsts;
    color_view_type kokcolors;

    nnz_lno_temp_work_view_t color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_view_type color_ban, tentative_color_ban; //color ban for each vertex represented with bit, as well as tentative color ban.
    size_type_temp_work_view_t edge_conf_indices;


    //      idx color_ban_init_val;  //the value to initialize the color_ban_. We avoid using the first bit representing the sign.
    //                                //Therefore if idx is int, it can represent 32-1 colors. Use color_set to represent more.

    /** \brief functor constructor.
     * \param srcs_ sources of the edgelist
     * \param dsts_ destinations of the edgelist
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     * \param tentative_color_ban_ :bit reprensentation of the neighbor vertex colors, that have been colored in the current iteration.
     *                              it is tentative, because coloring might have conflicts.
     * \param edge_conf_indices_ : The worklist for the edges.
     */
    halfedge_mark_conflicts (
        nnz_lno_persistent_work_view_t srcs_,
        nnz_lno_persistent_work_view_t dsts_,
        color_view_type kokcolors_,
        nnz_lno_temp_work_view_t color_set_,
        color_temp_work_view_type color_ban_,
        color_temp_work_view_type tentative_color_ban_,
        size_type_temp_work_view_t edge_conf_indices_):
      srcs(srcs_), dsts (dsts_),
      kokcolors(kokcolors_), color_set(color_set_),
      color_ban(color_ban_), tentative_color_ban(tentative_color_ban_),
      edge_conf_indices(edge_conf_indices_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type &ii) const {
      size_type work_index = edge_conf_indices(ii);
      //traverse edges,
      nnz_lno_t src_id = srcs(work_index);
      nnz_lno_t dst_id = dsts(work_index);


      color_t source_color = kokcolors(src_id);
      color_t dst_color = kokcolors(dst_id);

      //if the source and destionation have the same color, e.g. same color and same color_set.
      //then we have a conflict.
      char is_conflicted = (source_color != 0 && (source_color == dst_color) && (color_set(src_id) == color_set(dst_id)));
      if (is_conflicted){
        //this functor works both sides, although there is a reverse edge that will be encountered.
        //we want to mark the conflicts as soon as possible, so that those conflicted vertices neighbors wont have unnecessary conflicts.
        //idx conflict_ver = (src_id < dst_id) ? src_id : dst_id;
        //
        //TODO: dst_id seems to reduce the num colors, without increaisng runtime
        kokcolors(dst_id) = 0;
        tentative_color_ban(dst_id) = 0;
      }
    }
  };

  /**\brief Functor to count the number of conflicts
   * Also, as a side effect, it marks edge_conflict_marker and
   * remove those that are not needed to be looked further
   */
  struct halfedge_conflict_count{
    nnz_lno_persistent_work_view_t _kok_src;
    nnz_lno_persistent_work_view_t _kok_dst;
    color_view_type _kok_colors;
    nnz_lno_temp_work_view_t _color_set; //the colors that are represented with bits, and the colors set that the color is in.
    size_type_temp_work_view_t _edge_conflict_indices;
    char_temp_work_view_type _edge_conflict_marker;

    halfedge_conflict_count(
        nnz_lno_persistent_work_view_t kok_src_,
        nnz_lno_persistent_work_view_t kok_dst_,
        color_view_type kok_colors,
        nnz_lno_temp_work_view_t color_set,
        size_type_temp_work_view_t edge_conflict_indices,
        char_temp_work_view_type edge_conflict_marker):
          _kok_src( kok_src_),
          _kok_dst( kok_dst_),
          _kok_colors( kok_colors),
          _color_set( color_set),
          _edge_conflict_indices( edge_conflict_indices),
          _edge_conflict_marker( edge_conflict_marker){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type &ii, size_type &sum) const {

      size_type w = _edge_conflict_indices(ii);

      if (_edge_conflict_marker(w) == 0){
        sum += 1;
      }
      else{
        nnz_lno_t d = this->_kok_dst(w);
        nnz_lno_t s = this->_kok_src(w);

        color_t dc = _kok_colors(d);
        color_t sc = _kok_colors(s);

        if ( (dc && sc) || //if both colored
            (sc && (_color_set(d) > _color_set(s))) || //if source is colored, and destination color set is larger than source
            (dc && (_color_set(s) > _color_set(d))) //or if destionation is colored, and the source color set is larger
        ){
          //then no need to look at this edge anymore.
          _edge_conflict_marker(w) = 0;
          sum += 1;
        }
      }
    }
  };


  /**
   * \brief Functor to perform parallel prefix sum for edges so that the position
   * on the next conflictlist is calculated.
   */
  struct parallel_prefix_sum{
    size_type_temp_work_view_t _edge_conflict_indices;
    char_temp_work_view_type _edge_conflict_marker;
    size_type_temp_work_view_t _pps_view;

    parallel_prefix_sum(
        size_type_temp_work_view_t edge_conflict_indices,
        char_temp_work_view_type edge_conflict_marker,
        size_type_temp_work_view_t pps_view):
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker),
          _pps_view(pps_view){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type ii, size_t& update, const bool final) const {
      size_type w = _edge_conflict_indices(ii);
      if (final) {
        _pps_view(w) =  size_type(update);
      }
      update += _edge_conflict_marker(w);
    }
  };

  /**
   * \brief Functor to create the new work array.
   */
  struct create_new_work_array{
    size_type_temp_work_view_t _edge_conflict_indices;
    char_temp_work_view_type _edge_conflict_marker;
    size_type_temp_work_view_t _pps_view;
    size_type_temp_work_view_t _new_edge_conflict_indices;

    create_new_work_array(
        size_type_temp_work_view_t edge_conflict_indices,
        char_temp_work_view_type edge_conflict_marker,
        size_type_temp_work_view_t pps_view,
        size_type_temp_work_view_t new_edge_conflict_indices):
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker),
          _pps_view(pps_view),
          _new_edge_conflict_indices(new_edge_conflict_indices){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type ii) const {
      size_type w = _edge_conflict_indices(ii);
      if(_edge_conflict_marker(w)){
        size_type future_index = _pps_view(w);
        _new_edge_conflict_indices(future_index) = w;
      }
    }
  };

  /**
   * \brief Functor to create the new work array with atomic operations.
   */
  struct atomic_create_new_work_array{
    single_dim_index_view_type _new_index;
    size_type_temp_work_view_t _edge_conflict_indices;
    char_temp_work_view_type _edge_conflict_marker;
    size_type_temp_work_view_t _new_edge_conflict_indices;

    atomic_create_new_work_array(
        single_dim_index_view_type new_index,
        size_type_temp_work_view_t edge_conflict_indices,
        char_temp_work_view_type edge_conflict_marker,
        size_type_temp_work_view_t new_edge_conflict_indices):
          _new_index(new_index),
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker),
          _new_edge_conflict_indices(new_edge_conflict_indices){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type ii) const {
      typedef typename std::remove_reference< decltype( _new_index() ) >::type atomic_incr_type;
      size_type w = _edge_conflict_indices(ii);
      if(_edge_conflict_marker(w)){
        const size_type future_index = Kokkos::atomic_fetch_add( &_new_index(), atomic_incr_type(1));
        _new_edge_conflict_indices(future_index) = w;
      }
    }

  };


  /** \brief Functor for creating the ban colors for uncolored vertices.
   * It only creates ban_color, based on the certain information.
   * That is it works on the edges where one end is colored (for good)
   * and the other part is not colored.
   */
  struct halfedge_ban_colors {
    nnz_lno_persistent_work_view_t srcs, dsts;  //edge list, source and destinations of the edge list.
    color_view_type kokcolors;
    nnz_lno_temp_work_view_t color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_view_type color_ban; //color ban for each vertex represented with bit
    size_type_temp_work_view_t conflict_indices;
    char_temp_work_view_type edge_conflict_marker;

    /** \brief Functor constructor.
     * \param srcs_ sources of the edgelist
     * \param dsts_ destinations of the edgelist
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     */
    halfedge_ban_colors (
        nnz_lno_persistent_work_view_t srcs_, nnz_lno_persistent_work_view_t dsts_,
        color_view_type kokcolors_, nnz_lno_temp_work_view_t  color_set_,
        color_temp_work_view_type color_ban_,
        size_type_temp_work_view_t conflict_indices_, char_temp_work_view_type edge_conflict_marker_):
          srcs(srcs_), dsts (dsts_),
          kokcolors(kokcolors_), color_set(color_set_),
          color_ban(color_ban_),
          conflict_indices(conflict_indices_), edge_conflict_marker(edge_conflict_marker_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type &ii) const {
      size_type work_index = conflict_indices(ii);
      nnz_lno_t dst_id = dsts(work_index);
      color_t dst_col = kokcolors(dst_id);
      nnz_lno_t src_id = srcs(work_index);
      color_t src_col = kokcolors(src_id);

      //check destionation color.
      //continue only if it is not colored
      if ((!dst_col  && src_col) || (!src_col  && dst_col)){
        //check src color, send its color to ban colors only if it is colored.
        nnz_lno_t dest_col_set = color_set (dst_id);
        nnz_lno_t src_col_set = color_set (src_id);
        //check if they are in the same color set.
        //if they are not, we do not ban the color, as it represents a different color.
        if (src_col_set == dest_col_set){
          //atomic or, as no threads owns 'dst' (neither src)
          nnz_lno_t uncolored_vertex = dst_col? src_id: dst_id;
          Kokkos::atomic_fetch_or<color_t>(&(color_ban(uncolored_vertex)), src_col | dst_col);
          edge_conflict_marker(work_index) = 0;
        }
      }
    }
  };


  /**
   * \brief Functor to tentatively color vertices. It propogates the color information
   * to other end.
   */
  struct halfedge_expand_ban_for_unmatched_neighbors{
    nnz_lno_persistent_work_view_t srcs, dsts; //edge list, source and destinations of the edge list.
    color_view_type kokcolors;
    nnz_lno_temp_work_view_t color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_view_type color_ban, tentative_color_ban; //color ban for each vertex represented with bit, as well as tentative color ban.
    color_t first_digit;
    size_type_temp_work_view_t conflict_indices;


    /** \brief functor constructor.
     * \param srcs_ sources of the edgelist
     * \param dsts_ destinations of the edgelist
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     * \param tentative_color_ban_ :bit reprensentation of the neighbor vertex colors, that have been colored in the current iteration.
     *                              it is tentative, because coloring might have conflicts.
     */
    halfedge_expand_ban_for_unmatched_neighbors (
        nnz_lno_persistent_work_view_t srcs_, nnz_lno_persistent_work_view_t dsts_,
        color_view_type kokcolors_, nnz_lno_temp_work_view_t  color_set_ ,
        color_temp_work_view_type color_ban_, color_temp_work_view_type tentative_color_ban_,
        size_type_temp_work_view_t conflict_indices_):
          srcs(srcs_), dsts (dsts_),
          kokcolors(kokcolors_), color_set(color_set_),
          color_ban(color_ban_), tentative_color_ban(tentative_color_ban_),
          conflict_indices(conflict_indices_){
      color_t tmp = 1;
      first_digit = tmp <<( sizeof(color_t) * 8 -1);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type &ii) const {
      size_type work_index = conflict_indices(ii);
      nnz_lno_t dst_id = dsts(work_index);
      color_t dst_col = kokcolors(dst_id);

      //if the destionation is colored already, we have nothing to do.
      //otherwise, if destionation is uncolored, or if its color < 0 (it has been tentatively colored)
      //then we need to check the source.
      if (dst_col == 0 || (dst_col & first_digit) ){
        nnz_lno_t src_id = srcs(work_index);
        color_t src_col = kokcolors(src_id);
        //if source is colored, again we have nothing to do.
        //if it is tentatively colored or uncolored, then we have work to do.
        if (src_col == 0 || (src_col & first_digit)){
          //check their colors sets, if they are on different color sets,
          //we dont need to care about the prohibted colors on each other -- at least in this iteration.
          nnz_lno_t dest_col_set = color_set (dst_id);
          nnz_lno_t src_col_set = color_set (src_id);

          if (src_col_set == dest_col_set){
            if ((dst_col & first_digit)  && (src_col & first_digit)){
              //if both ends are tentatively colored, we can check for the conflict here,
              //and fix it as below. However doing so increased the number of colors,
              //so, it has been turned of for now.

              if (src_col == dst_col && dest_col_set == src_col_set){
                //idx smaller_index = (src_id > dst_id) ? src_id : dst_id;
                nnz_lno_t smaller_index = dst_id; //TODO which one is better? this seems to be not much changing
                //idx smaller_index = src_id;
                //then both have been colored tentavitely. propoagate the color of src to dst.
                Kokkos::atomic_fetch_or<color_t>(&(tentative_color_ban(smaller_index)), -src_col);
                nnz_lno_t banned_colors = ~(color_ban(smaller_index) | tentative_color_ban(smaller_index));
                nnz_lno_t larger_col = banned_colors & (-banned_colors);
                kokcolors(smaller_index) = -(larger_col);
              }
            }
            else if(src_col != 0) {
              //if src is tentavily colored, and dst is not colored,
              //then we send the color information to dst's tentative_ban.

              //Kokkos::atomic_fetch_or<color_type>(&(color_ban(dst_id)), -src_col);
              Kokkos::atomic_fetch_or<color_t>(&(tentative_color_ban(dst_id)), -src_col);
            }
            else if (dst_col != 0){
              //if it is dst tentatively colors, but src is not colored,
              //then we send the dst color info to src's tentative_ban

              //Kokkos::atomic_fetch_or<color_type>(&(color_ban(src_id)), -dst_col);
              Kokkos::atomic_fetch_or<color_t>(&(tentative_color_ban(src_id)), -dst_col);
            }
            else {
              //idx smaller_index = src_id < dst_id > 0 ? src_id: dst_id;
              //idx larger_index = src_id < dst_id > 0 ? dst_id : src_id;
#ifndef TOOHIGHQUALITY
              nnz_lno_t smaller_index = src_id;
              nnz_lno_t larger_index = dst_id;
#endif
#ifdef TOOHIGHQUALITY
              row_index_type smaller_index = dst_id;
              row_index_type larger_index = src_id;
#endif

              //idx smaller_col =  src_id < dst_id > 0 ? src_col: dst_col;
              //if both ends are uncolored, tentatively color the the source if its index is smaller than dst.
              //make an 'bitwise or' of color_ban and tentative_color_ban to get the all prohibited colors.
              //we need to find the right most zero here. it is easier to find right most 1, so we do a not of the result color ban.
              color_t banned_colors = ~(color_ban(smaller_index) | tentative_color_ban(smaller_index));
              //the 'bitwise and' of banned_colors with two's complement result in only the rightmost 1 to be set, which is our color.
              src_col = banned_colors & (-banned_colors);
              //set it to minus of the color, as it is tentative coloring.
              kokcolors(smaller_index) = -(src_col);
              //send the color information to dst's tentative color ban.
              Kokkos::atomic_fetch_or<color_t>(&(tentative_color_ban(larger_index)), src_col);
              //Kokkos::atomic_fetch_or<color_type>(&(color_ban(dst_id)), src_col);
            }
          }
        }
      }
    }
  };





  /** \brief Functor responsible for choosing a color for each uncolored vertex,
   * given the color_ban and tentative_color_ban
   */
  struct choose_colors {
    color_view_type kokcolors;
    nnz_lno_temp_work_view_t color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_view_type color_ban, tentative_color_ban;  //color ban for each vertex represented with bit, as well as tentative color ban.
    color_t color_ban_init_val;  //the value to initialize the color_ban_. We avoid using the first bit representing the sign.
    //Therefore if idx is int, it can represent 32-1 colors. Use color_set to represent more

    /** \brief functor constructor.
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     * \param tentative_color_ban_ :bit reprensentation of the neighbor vertex colors, that have been colored in the current iteration.
     *                              it is tentative, because coloring might have conflicts.
     */
    choose_colors ( color_view_type kokcolors_, nnz_lno_temp_work_view_t  color_set_,
        color_temp_work_view_type color_ban_,  color_temp_work_view_type tentative_color_ban_):
          kokcolors(kokcolors_), color_set(color_set_),
          color_ban(color_ban_), tentative_color_ban(tentative_color_ban_){
      //color_ban should always have 1 at the first bit, so that that color is not allowed.
      color_t tmp = 1;
      color_ban_init_val = tmp <<( sizeof(color_t) * 8 -1);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t &ii) const {
      //if the vertex is uncolored, we will choose a new color for the vertex.
      if (kokcolors[ii] == 0){
        color_t certain_info = color_ban(ii);
        //get the banned_color_set by taking 'bitwise or' or color ban and tentative_color_ban
        color_t banned_colors = ~(certain_info | tentative_color_ban(ii));
        //my color is the first non set bit in the banned_colors. We perform a not operation,
        //and make a 'bitwise and' with its 2's complement to find the first zero bit.
        color_t my_color = banned_colors & (-banned_colors);
        if (my_color == 0){
#ifdef EBCOLORING_HIGHER_QUALITY
          //if my color is zero, that is all the available colors in this set has been taken by the neighbors
          //then I might need to change my color set. But we need to be sure about this, so we check if
          //color_ban is full as well, since tentative_color_ban might be too pessimist.
          banned_colors = ~(certain_info);
          my_color = banned_colors & (-banned_colors);
          if (my_color == 0){
#endif
            //if all colors are taken by the neighbors (certainly excluding tentative colors), then
            //I need to change my color set.
            //if there are still available colors w.r.t. color_ban, then try one more time without increasing the color_set.
            color_set(ii) += 1; // increase color set.
            color_ban(ii) = color_ban_init_val; //set the color ban to its initial value.
#ifdef EBCOLORING_HIGHER_QUALITY
          }
#endif
          //in each case we cannot color this vertex. set the tentative_color_ban to 0
          //try to color it at the next iteration.
          tentative_color_ban(ii) = 0;
          //color_ban(ii) = color_ban_init_val; //set the color ban to its initial value.
        }
        else {
          kokcolors(ii) = my_color;
        }
      }
      else if (kokcolors(ii) & color_ban_init_val){
        kokcolors(ii) = -kokcolors(ii);
      }
    }
  };

  /** Functor responsible for setting the final color for each vertex,
   *  The color of the vertex is found ast color_set * (sizeof(color_type) * 8 -1) + log2(color)
   */
  struct set_final_colors{
    color_view_type kokcol;
    nnz_lno_temp_work_view_t kokcolset; //the colors that are represented with bits, and the colors set that the color is in.
    color_t color_size;

    /** \brief functor constructor.
     * \param kokcol_  the colors of the vertices. Represented with bits.
     * \param kokcolset_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     */
    set_final_colors(color_view_type kokcol_, nnz_lno_temp_work_view_t  kokcolset_): kokcol(kokcol_),kokcolset(kokcolset_), color_size ( sizeof(color_t) * 8 -1){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t &ii) const {
      nnz_lno_t i = 0;
      color_t val = kokcol(ii);
      //if check below is necessary.
      // this happens when a vertices all neighbors are colored,
      //so the information from all neighbors are taken, no edge to be processed by this vertex.
      //the algorithm works on the number of edges, if the edges are all consumed, the loop
      //might terminate with an early exit without coloring this vertex.
      //this happens when all neighbors consumes all the colors in the current vertex set,
      //and the vertex left to be colored in the next iteration.
      //but the vertex couldnt be colored, because there is no more edge left to be worked on.

      if(val == 0) val = 1;

      //find the position in the bit.
      while (val != 0) {
        ++i;
        val = val >> 1;
      }
      //set the final color.
      kokcol(ii) = i + kokcolset(ii) * color_size;
    }
  };
};

}
}

#endif //_KOKKOSCOLORINGIMP_HPP
