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
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdexcept>
#ifndef _KOKKOSKERNELSIOUTILS_HPP
#define _KOKKOSKERNELSIOUTILS_HPP

#include "Kokkos_ArithTraits.hpp"
#include <Kokkos_Core.hpp>
#include "KokkosKernels_SimpleUtils.hpp"
#include <sys/stat.h>

namespace KokkosKernels{


namespace Impl{


//MD: Bases on Christian's sparseMatrix_generate function in test_crsmatrix.cpp file.
template< typename ScalarType , typename OrdinalType, typename SizeType>
void kk_sparseMatrix_generate(
    OrdinalType nrows,
    OrdinalType ncols,
    SizeType &nnz,
    OrdinalType row_size_variance,
    OrdinalType bandwidth,
    ScalarType* &values,
    SizeType* &rowPtr,
    OrdinalType* &colInd)
{
  rowPtr = new SizeType[nrows+1];

  OrdinalType elements_per_row = nnz/nrows;
  srand(13721);
  rowPtr[0] = 0;
  for(int row=0;row<nrows;row++)
  {
    int varianz = (1.0*rand()/INT_MAX-0.5)*row_size_variance;
    rowPtr[row+1] = rowPtr[row] + elements_per_row+varianz;
  }
  nnz = rowPtr[nrows];
  values = new ScalarType[nnz];
  colInd = new OrdinalType[nnz];
  for(OrdinalType row=0;row<nrows;row++)
  {

    for(SizeType k=rowPtr[row] ;k<rowPtr[row+1];k++)
    {

      while (true){
        OrdinalType pos = (1.0*rand()/INT_MAX-0.5)*bandwidth+row;
        if(pos<0) pos+=ncols;
        if(pos>=ncols) pos-=ncols;

        bool is_already_in_the_row = false;
        for(SizeType j = rowPtr[row] ; j<k ;j++){
          if (colInd[j] == pos){
            is_already_in_the_row = true;
            break;
          }
        }
        if (!is_already_in_the_row) {

          colInd[k]= pos;
          values[k] = 100.0*rand()/INT_MAX-50.0;
          break;
        }
      }


    }
  }
}

template< typename ScalarType , typename OrdinalType, typename SizeType>
void kk_sparseMatrix_generate_lower_upper_triangle(
    char uplo,
    OrdinalType nrows,
    OrdinalType ncols,
    SizeType &nnz,
    OrdinalType row_size_variance,
    OrdinalType bandwidth,
    ScalarType* &values,
    SizeType* &rowPtr,
    OrdinalType* &colInd)
{
  rowPtr = new SizeType[nrows+1];

  //OrdinalType elements_per_row = nnz/nrows;
  srand(13721);
  rowPtr[0] = 0;
  for(int row=0;row<nrows;row++)
  {
    //int varianz = (1.0*rand()/INT_MAX-0.5)*row_size_variance;
    if (uplo =='L')
      rowPtr[row+1] = rowPtr[row] + row + 1;
    else
      rowPtr[row+1] = rowPtr[row] + ncols - (row) ;
  }
  nnz = rowPtr[nrows];
  values = new ScalarType[nnz];
  colInd = new OrdinalType[nnz];
  for(OrdinalType row=0;row<nrows;row++)
  {

    for(SizeType k=rowPtr[row]; k<rowPtr[row+1]; k++) {
      if (uplo =='L')
        colInd[k]= k - rowPtr[row];
      else
        colInd[k]= row + (k - rowPtr[row]);
      values[k] = 1.0;
    }
  }
}

template< typename ScalarType , typename OrdinalType, typename SizeType>
void kk_diagonally_dominant_sparseMatrix_generate(
    OrdinalType nrows,
    OrdinalType ncols,
    SizeType &nnz,
    OrdinalType row_size_variance,
    OrdinalType bandwidth,
    ScalarType* &values,
    SizeType* &rowPtr,
    OrdinalType* &colInd)
{
  rowPtr = new SizeType[nrows+1];

  OrdinalType elements_per_row = nnz/nrows;
  srand(13721);
  rowPtr[0] = 0;
  for(int row=0;row<nrows;row++)
  {
    int varianz = (1.0*rand()/INT_MAX-0.5)*row_size_variance;
    rowPtr[row+1] = rowPtr[row] + elements_per_row+varianz;
  }
  nnz = rowPtr[nrows];
  values = new ScalarType[nnz];
  colInd = new OrdinalType[nnz];
  const ScalarType temp = 10;
  for(OrdinalType row=0;row<nrows;row++)
  {
    ScalarType total_values = 0;
    for(SizeType k=rowPtr[row] ;k<rowPtr[row+1] - 1;k++)
    {
      OrdinalType pos = row;
      while (pos == row){
        pos = ((1.0*rand())/INT_MAX-0.5)*bandwidth+row;
      }
      if(pos<0) pos+=ncols;

      if(pos>=ncols) pos-=ncols;
      colInd[k]= pos;
      values[k] = 100.0*rand()/INT_MAX-50.0;
      total_values += Kokkos::Details::ArithTraits<ScalarType>::abs(values[k]);
    }

    colInd[rowPtr[row+1] - 1]= row;
    values[rowPtr[row+1] - 1] = total_values * temp;
  }
}

template <typename crsMat_t>
crsMat_t kk_generate_diagonally_dominant_sparse_matrix(
    typename crsMat_t::const_ordinal_type nrows,
    typename crsMat_t::const_ordinal_type ncols,
    typename crsMat_t::non_const_size_type &nnz,
    typename crsMat_t::const_ordinal_type row_size_variance,
    typename crsMat_t::const_ordinal_type bandwidth){

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
  typedef typename graph_t::entries_type::non_const_type   cols_view_t;
  typedef typename crsMat_t::values_type::non_const_type values_view_t;


  typedef typename row_map_view_t::non_const_value_type size_type;
  typedef typename cols_view_t::non_const_value_type lno_t;
  typedef typename values_view_t::non_const_value_type scalar_t;
  lno_t *adj;
  size_type *xadj;//, nnzA;
  scalar_t *values;

  kk_diagonally_dominant_sparseMatrix_generate<scalar_t, lno_t, size_type>(
      nrows, ncols, nnz, row_size_variance,  bandwidth,
      values, xadj, adj);

  row_map_view_t rowmap_view("rowmap_view", nrows+1);
  cols_view_t columns_view("colsmap_view", nnz);
  values_view_t values_view("values_view", nnz);

  {
    typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
    typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
    typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

    for (lno_t i = 0; i <= nrows; ++i){
      hr(i) = xadj[i];
    }

    for (size_type i = 0; i < nnz; ++i){
      hc(i) = adj[i];
      hv(i) = values[i];
    }
    Kokkos::deep_copy (rowmap_view , hr);
    Kokkos::deep_copy (columns_view , hc);
    Kokkos::deep_copy (values_view , hv);
  }

  graph_t static_graph (columns_view, rowmap_view);
  crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
  delete [] xadj; delete [] adj; delete [] values;
  return crsmat;
}


template <typename crsMat_t>
crsMat_t kk_generate_triangular_sparse_matrix(
    char uplo,
    typename crsMat_t::const_ordinal_type nrows,
    typename crsMat_t::const_ordinal_type ncols,
    typename crsMat_t::non_const_size_type &nnz,
    typename crsMat_t::const_ordinal_type row_size_variance,
    typename crsMat_t::const_ordinal_type bandwidth){

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
  typedef typename graph_t::entries_type::non_const_type   cols_view_t;
  typedef typename crsMat_t::values_type::non_const_type values_view_t;


  typedef typename row_map_view_t::non_const_value_type size_type;
  typedef typename cols_view_t::non_const_value_type lno_t;
  typedef typename values_view_t::non_const_value_type scalar_t;
  lno_t *adj;
  size_type *xadj;//, nnzA;
  scalar_t *values;

  kk_sparseMatrix_generate_lower_upper_triangle<scalar_t, lno_t, size_type>(
      uplo,
      nrows, ncols, nnz, row_size_variance,  bandwidth,
      values, xadj, adj);

  row_map_view_t rowmap_view("rowmap_view", nrows+1);
  cols_view_t columns_view("colsmap_view", nnz);
  values_view_t values_view("values_view", nnz);

  {
    typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
    typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
    typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

    for (lno_t i = 0; i <= nrows; ++i){
      hr(i) = xadj[i];
    }

    for (size_type i = 0; i < nnz; ++i){
      hc(i) = adj[i];
      hv(i) = values[i];
    }
    Kokkos::deep_copy (rowmap_view , hr);
    Kokkos::deep_copy (columns_view , hc);
    Kokkos::deep_copy (values_view , hv);
    Kokkos::fence();
  }

  graph_t static_graph (columns_view, rowmap_view);
  crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
  delete [] xadj; delete [] adj; delete [] values;
  return crsmat;
}

template <typename crsMat_t>
crsMat_t kk_generate_sparse_matrix(
    typename crsMat_t::const_ordinal_type nrows,
    typename crsMat_t::const_ordinal_type ncols,
    typename crsMat_t::non_const_size_type &nnz,
    typename crsMat_t::const_ordinal_type row_size_variance,
    typename crsMat_t::const_ordinal_type bandwidth){

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
  typedef typename graph_t::entries_type::non_const_type   cols_view_t;
  typedef typename crsMat_t::values_type::non_const_type values_view_t;


  typedef typename row_map_view_t::non_const_value_type size_type;
  typedef typename cols_view_t::non_const_value_type lno_t;
  typedef typename values_view_t::non_const_value_type scalar_t;
  lno_t *adj;
  size_type *xadj;//, nnzA;
  scalar_t *values;

  kk_sparseMatrix_generate<scalar_t, lno_t, size_type>(
      nrows, ncols, nnz, row_size_variance,  bandwidth,
      values, xadj, adj);

  row_map_view_t rowmap_view("rowmap_view", nrows+1);
  cols_view_t columns_view("colsmap_view", nnz);
  values_view_t values_view("values_view", nnz);

  {
    typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
    typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
    typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

    for (lno_t i = 0; i <= nrows; ++i){
      hr(i) = xadj[i];
    }

    for (size_type i = 0; i < nnz; ++i){
      hc(i) = adj[i];
      hv(i) = values[i];
    }
    Kokkos::deep_copy (rowmap_view , hr);
    Kokkos::deep_copy (columns_view , hc);
    Kokkos::deep_copy (values_view , hv);
    Kokkos::fence();
  }

  graph_t static_graph (columns_view, rowmap_view);
  crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
  delete [] xadj; delete [] adj; delete [] values;
  return crsmat;
}




//TODO: need to fix the size_type. All over the reading inputs are lno_t.

template <typename stype>
void md_malloc(stype **arr, size_t n, std::string alloc_str = ""){
  *arr = new stype[n];
  if (*arr == NULL){
    throw std::runtime_error ("Memory Allocation Problem\n");
  }
}

template <typename idx, typename wt>
struct Edge{
  idx src;
  idx dst;
  wt ew;
  bool operator<(const Edge <idx,wt> & a) const
  {
    //return !((this->src < a.src) || (this->src == a.src && this->dst < a.dst));
    return (this->src < a.src) || (this->src == a.src && this->dst < a.dst);
  }
};


////////////////////////////////////////////////////////////////////////////////
// From MTGL
////////////////////////////////////////////////////////////////////////////////
inline size_t kk_get_file_size(const char* file)
{
  struct stat stat_buf;

#ifdef _WIN32
  int retval = _stat(file, &stat_buf);
#else
  int retval = stat(file, &stat_buf);
#endif

  return retval == 0 ? size_t(stat_buf.st_size) : size_t(0);
}

template <typename lno_t>
void buildEdgeListFromBinSrcTarg_undirected(
    const char *fnameSrc, const char*fnameTarg,
    size_t &numEdges,
    lno_t **srcs, lno_t **dst){
   size_t srcFileSize = kk_get_file_size(fnameSrc);
   size_t trgFileSize = kk_get_file_size(fnameTarg);
   // test these values

   size_t srcSize = srcFileSize / sizeof(lno_t);
   size_t trgSize = trgFileSize / sizeof(lno_t);
   if (srcSize != trgSize)
   {
     throw std::runtime_error ("Src and Target file needs to be the same size");
   }
   //Assumption that each edge is listed once
   numEdges = srcSize;

   md_malloc<lno_t>(srcs, numEdges);
   md_malloc<lno_t>(dst, numEdges);

   ////////////////////////////////////////////////////////
   // Read source data into buffer
   ////////////////////////////////////////////////////////

   std::ifstream myFile (fnameSrc, std::ios::in | std::ios::binary);

   myFile.read((char *) *srcs, sizeof(lno_t) * (numEdges));

   myFile.close();

   std::ifstream myFile2 (fnameTarg, std::ios::in | std::ios::binary);

   myFile2.read((char *) *dst, sizeof(lno_t) * (numEdges));

   myFile2.close();
   //

}



template <typename idx_array_type>
inline void kk_write_1Dview_to_file(idx_array_type view, const char *filename){

  typedef typename idx_array_type::HostMirror host_type;
  //typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view (view);
  Kokkos::deep_copy (host_view , view);
  Kokkos::fence();
  std::ofstream myFile (filename, std::ios::out );
  for (size_t i = 0; i < view.dimension_0(); ++i){
	  myFile << host_view(i) << std::endl;
  }
  myFile.close();
}

template <typename idx_array_type>
inline void kk_read_1Dview_from_file(idx_array_type &view, const char *filename){

  typedef typename idx_array_type::HostMirror host_type;
  //typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view (view);
  std::ifstream myFile (filename, std::ios::in );

  for (size_t i = 0; i < view.dimension_0(); ++i){
	  myFile >> host_view(i);
  }
  myFile.close();
  Kokkos::deep_copy (view, host_view);
  Kokkos::fence();
}




template <typename idx>
void convert_crs_to_lower_triangle_edge_list(idx nv, idx *xadj, idx *adj, idx *lower_triangle_srcs, idx *lower_triangle_dests){
  idx ind = 0;
  for (idx i = 0; i < nv; ++i){
    idx xb = xadj[i];
    idx xe = xadj[i+1];
    for (idx j = xb; j < xe; ++j){
      idx dst = adj[j];
      if (i < dst){
        lower_triangle_srcs[ind] = i;
        lower_triangle_dests[ind++] = dst;
      }
    }
  }
}

template <typename idx>
void convert_crs_to_edge_list(idx nv, idx *xadj, idx *srcs){
  for (idx i = 0; i < nv; ++i){
    idx xb = xadj[i];
    idx xe = xadj[i+1];
    for (idx j = xb; j < xe; ++j){
      srcs[j] = i;
    }
  }
}

template <typename size_type, typename lno_t, typename wt>
void convert_edge_list_to_csr (lno_t nv, size_type ne,
    lno_t *srcs, lno_t *dests, wt *ew,
    size_type *xadj, lno_t *adj, wt *crs_ew){

  std::vector <struct Edge<lno_t, wt> > edges (ne);
  for(size_type i = 0; i < ne; ++i){
    edges[i].src = srcs[i];
    edges[i].dst = dests[i];
    edges[i].ew = ew[i];
  }
  std::sort (edges.begin(), edges.begin() + ne);

  size_type eind = 0;
  for (lno_t i = 0; i < nv; ++i){
    (xadj)[i] = eind;
    while (edges[eind].src == i){
      (adj)[eind] = edges[eind].dst;
      (*crs_ew)[eind] = edges[eind].ew;
      ++eind;
    }
  }
  xadj[nv] = eind;

}

template <typename in_lno_t, typename size_type, typename lno_t>
void convert_undirected_edge_list_to_csr (
    lno_t nv, size_type ne,
    in_lno_t *srcs, in_lno_t *dests,
    size_type *xadj, lno_t *adj){

  std::vector <struct Edge<lno_t, double> > edges (ne * 2);
  for(size_type i = 0; i < ne; ++i){
    edges[i * 2].src = srcs[i];
    edges[i * 2].dst = dests[i];

    edges[i * 2 + 1].src = dests[i];
    edges[i * 2 + 1].dst = srcs[i];
  }
#ifdef KOKKOSKERNELS_HAVE_OUTER
#include<parallel/multiseq_selection.h>
#include<parallel/multiway_merge.h>
#include<parallel/merge.h>
#include<parallel/multiway_mergesort.h>
  __gnu_parallel::parallel_sort_mwms<false,true, struct Edge<lno_t, double> *>
  (&(edges[0]), &(edges[0])+ne*2,
      std::less<struct Edge<lno_t, double> >(), 64);
#else
  std::sort (edges.begin(), edges.begin() + ne * 2);
#endif


  size_type eind = 0;
  for (lno_t i = 0; i < nv; ++i){
    (xadj)[i] = eind;
    while (edges[eind].src == i){
      (adj)[eind] = edges[eind].dst;
      //(*crs_ew)[eind] = edges[eind].ew;
      ++eind;
    }
  }
  xadj[nv] = eind;
}
/*

template <typename lno_t, typename size_type, typename scalar_t>
void read_graph_src_dst_bin(
    lno_t *nv, size_type *ne
    ,size_type **xadj, lno_t **adj, scalar_t **ew,
    const char *fnameSrc, const char *fnameTarg){

  size_t numEdges = 0;
  size_t *srcs, *dst; //this type is hard coded
  buildEdgeListFromBinSrcTarg_undirected(
      fnameSrc, fnameTarg,
      &numEdges,
      &srcs, &dst);

  lno_t num_vertex = 0;
  for (size_t i = 0; i < numEdges; ++i){
    if (num_vertex < srcs[i]) num_vertex = srcs[i];
    if (num_vertex < dst[i]) num_vertex = dst[i];
  }
  num_vertex += 1;

  *nv = num_vertex;
  *ne = numEdges * 2;

  md_malloc<size_type>(xadj, num_vertex + 1);
  md_malloc<lno_t>(adj, numEdges * 2);
  convert_undirected_edge_list_to_csr (
      num_vertex, numEdges,
      srcs, dst,
      *xadj, *adj);

  delete [] srcs;
  delete [] dst;
}
*/


template <typename idx, typename wt>
void write_edgelist_bin(
    size_t ne,
    const idx *edge_begins,
    const  idx *edge_ends,
    const  wt *ew,
    const  char *filename){
  std::ofstream myFile (filename, std::ios::out | std::ios::binary);
  myFile.write((char *) &ne, sizeof(idx));
  myFile.write((char *) edge_begins, sizeof(idx) * (ne));
  myFile.write((char *) edge_ends, sizeof(idx) * (ne));
  myFile.write((char *) ew, sizeof(wt) * (ne));
  myFile.close();
}

template <typename idx, typename wt>
void read_edgelist_bin(
    idx *ne,
    idx **edge_begins,
    idx **edge_ends,
    wt **ew,
    const  char *filename){

  std::ifstream myFile (filename, std::ios::in | std::ios::binary);


  myFile.read((char *) ne, sizeof(idx));
  md_malloc<idx>(edge_begins, *ne);
  md_malloc<idx>(edge_ends, *ne);
  md_malloc<wt> (ew, *ne);
  myFile.read((char *) *edge_begins, sizeof(idx) * (*ne));
  myFile.read((char *) *edge_ends, sizeof(idx) * (*ne));
  myFile.read((char *) *ew, sizeof(wt) * (*ne));
  myFile.close();
}





template <typename lno_t, typename size_type, typename scalar_t>
void write_graph_bin(lno_t nv, size_type ne,const size_type *xadj,const  lno_t *adj,const  scalar_t *ew,const  char *filename){
  std::ofstream myFile (filename, std::ios::out | std::ios::binary);
  myFile.write((char *) &nv, sizeof(lno_t));
  myFile.write((char *) &ne, sizeof(size_type));
  myFile.write((char *) xadj, sizeof(size_type) * (nv + 1));

  myFile.write((char *) adj, sizeof(lno_t) * (ne));

  myFile.write((char *) ew, sizeof(scalar_t) * (ne));

  myFile.close();
}

template <typename lno_t, typename size_type, typename scalar_t>
void write_graph_crs(lno_t nv, size_type ne,const size_type *xadj,const  lno_t *adj,const  scalar_t *ew,const  char *filename){
  std::ofstream myFile (filename, std::ios::out );
  myFile << nv << " " << ne << std::endl;

  for (lno_t i = 0; i <= nv; ++i){
    myFile  << xadj[i] << " ";
  }
  myFile  << std::endl;

  for (lno_t i = 0; i < nv; ++i){
    size_type b = xadj[i];
    size_type e = xadj[i + 1];
    for (size_type j = b; j < e; ++j){
      myFile  << adj[j] << " ";
    }
    myFile  << std::endl;
  }
  for (lno_t i = 0; i < ne; ++i){
    myFile  << ew[i] << " ";
  }
  myFile  << std::endl;

  myFile.close();
}

template <typename lno_t, typename size_type, typename scalar_t>
void write_graph_ligra(lno_t nv, size_type ne,const size_type *xadj,const  lno_t *adj,const  scalar_t *ew,const  char *filename){

  std::ofstream ff (filename);
  ff << "AdjacencyGraph" << std::endl;
  ff << nv << std::endl << ne << std::endl;
  for (lno_t i = 0; i < nv; ++i){
    ff << xadj[i] << std::endl;
  }
  for (size_type i = 0; i < ne; ++i){
    ff << adj[i] << std::endl;
  }
  ff.close();
}


template <typename lno_t, typename size_type, typename scalar_t>
void write_graph_mtx(lno_t nv, size_type ne,const size_type *xadj,const  lno_t *adj,const  scalar_t *ew,const  char *filename){

  std::ofstream myFile (filename);
  myFile  << "%%MatrixMarket matrix coordinate real general" << std::endl;
  myFile  << nv << " " << nv << " " << ne << std::endl;
  for (lno_t i = 0; i < nv; ++i){
    size_type b = xadj[i];
    size_type e = xadj[i + 1];
    for (size_type j = b; j < e; ++j){
      myFile  << i + 1 << " " << (adj)[j] + 1 << " " << ew [j] << std::endl;
    }
  }

  myFile.close();
}



template <typename lno_t, typename size_type, typename scalar_t>
void read_graph_bin(lno_t *nv, size_type *ne,size_type **xadj, lno_t **adj, scalar_t **ew, const char *filename){

  std::ifstream myFile (filename, std::ios::in | std::ios::binary);

  myFile.read((char *) nv, sizeof(lno_t));
  myFile.read((char *) ne, sizeof(size_type));
  md_malloc<size_type>(xadj, *nv+1);
  md_malloc<lno_t>(adj, *ne);
  md_malloc<scalar_t> (ew, *ne);
  myFile.read((char *) *xadj, sizeof(size_type) * (*nv + 1));
  myFile.read((char *) *adj, sizeof(lno_t) * (*ne));
  myFile.read((char *) *ew, sizeof(scalar_t) * (*ne));
  myFile.close();
}

template <typename lno_t, typename size_type, typename scalar_t>
void read_graph_crs(lno_t *nv, size_type *ne,size_type **xadj, lno_t **adj, scalar_t **ew, const char *filename){

  std::ifstream myFile (filename, std::ios::in );
  myFile >> *nv >> *ne;

  md_malloc<size_type>(xadj, *nv+1);
  md_malloc<lno_t>(adj, *ne);
  md_malloc<scalar_t> (ew, *ne);

  for (lno_t i = 0; i <= *nv; ++i){
    myFile  >> (*xadj)[i];
  }

  for (size_type i = 0; i < *ne; ++i){
    myFile  >> (*adj)[i];
  }
  for (size_type i = 0; i < *ne; ++i){
    myFile  >> (*ew)[i];
  }
  myFile.close();
}




inline bool endswith (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}


template <typename crs_matrix_t>
void write_kokkos_crst_matrix(crs_matrix_t a_crsmat,const  char *filename){


  typedef typename crs_matrix_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
  typedef typename graph_t::entries_type::non_const_type   cols_view_t;
  typedef typename crs_matrix_t::values_type::non_const_type   values_view_t;

  typedef typename cols_view_t::value_type lno_t;
  typedef typename row_map_view_t::value_type size_type;
  typedef typename values_view_t::value_type scalar_t;

  std::string strfilename(filename);
  if (endswith(strfilename, ".mtx")){
    write_graph_mtx<lno_t, size_type, scalar_t>(a_crsmat.numRows(),
        a_crsmat.graph.entries.dimension_0(),
        a_crsmat.graph.row_map.data(),
        a_crsmat.graph.entries.data(),
        a_crsmat.values.data(),filename);
  }

  else if (endswith(strfilename, ".bin")){
    write_graph_bin<lno_t, size_type, scalar_t>(a_crsmat.numRows(),
        a_crsmat.graph.entries.dimension_0(),
        a_crsmat.graph.row_map.data(),
        a_crsmat.graph.entries.data(),
        a_crsmat.values.data(),filename);
  }
  else if (endswith(strfilename, ".ligra")){
    write_graph_ligra<lno_t, size_type, scalar_t>(a_crsmat.numRows(),
        a_crsmat.graph.entries.dimension_0(),
        a_crsmat.graph.row_map.data(),
        a_crsmat.graph.entries.data(),
        a_crsmat.values.data(),filename);
  }
  else if (endswith(strfilename, ".crs")){
    write_graph_crs<lno_t, size_type, scalar_t>(a_crsmat.numRows(),
        a_crsmat.graph.entries.dimension_0(),
        a_crsmat.graph.row_map.data(),
        a_crsmat.graph.entries.data(),
        a_crsmat.values.data(),filename);

  }
  else {
    throw std::runtime_error ("Writer is not available\n");
  }


}




template <typename lno_t, typename size_type, typename scalar_t>
int read_mtx (
    const char *fileName,
    lno_t *nv, size_type *ne,
    size_type **xadj, lno_t **adj, scalar_t **ew,
    bool symmetrize = false, bool remove_diagonal = true,
    bool transpose = false){

  std::ifstream mmf (fileName, std::ifstream::in);
  if (!mmf.is_open()) {
    throw std::runtime_error ("File cannot be opened\n");
  }

  std::string fline = "";
  getline(mmf, fline);

  if (fline.size() < 2 || fline[0] != '%' || fline[1] != '%'){
    throw std::runtime_error ("Invalid MM file. Line-1\n");
  }


  int mtx_object = 0; // 0- matrix, 1-vector
  int mtx_format = 0; // 0- coordinate; 1- array
  //int mtx_field = 0; //0-real, 1-double, 2-complex, 3- integer, 4-pattern
  int mtx_sym = 0; //0-general, 1-symmetric, 2-skew-symmetric, 3-hermitian

  if (fline.find("matrix") != std::string::npos){
    mtx_object = 0;
  } else if (fline.find("vector") != std::string::npos){
    mtx_object = 1;
  }

  if (fline.find("coordinate") != std::string::npos){
    mtx_format = 0;
  }
  else if (fline.find("array") == std::string::npos){
    mtx_format = 1;
  }

  if (fline.find("real") != std::string::npos){
    //mtx_field = 0;
  }
  else if (fline.find("double") != std::string::npos){
    //mtx_field = 1;
  }
  else if (fline.find("complex") != std::string::npos){
    //mtx_field = 2;
  }
  else if (fline.find("integer") != std::string::npos){
    //mtx_field = 3;
  }
  else if (fline.find("pattern") != std::string::npos){
    //mtx_field = 4;
  }

  if (fline.find("skew-symmetric") != std::string::npos){
    mtx_sym = 2;
  }
  else if (fline.find("symmetric") != std::string::npos){
    mtx_sym = 1;
  }
  else if (fline.find("hermitian") != std::string::npos){
    mtx_sym = 3;
  }
  else if (fline.find("general") != std::string::npos){
    mtx_sym = 0;
  }

  if (mtx_object == 1) {
    throw std::runtime_error ("VECTOR TYPE NOT HANDLED YET\n");
  }
  if (mtx_format == 1) {
    throw std::runtime_error ("ARRAY TYPE NOT HANDLED YET\n");
  }
  if (!symmetrize && (mtx_sym == 2 || mtx_sym == 3)) {
    throw std::runtime_error ("SKEW-SYMMETRIC and HERMITIAN TYPE NOT HANDLED YET\n");
  }


  while(1){
    getline(mmf, fline);
    if(fline[0] != '%') break;
  }
  std::stringstream ss (fline);
  lno_t nr = 0, nc = 0;
  size_type nnz = 0;

  ss >> nr >> nc >> nnz;


  //if (nr != nc) {std::cerr << "NON-SQUARE MATRIX TYPE NOT HANDLED YET"<< std::endl; return (1); }
  size_type noEdges = nnz;
  if (mtx_sym == 1 || symmetrize) noEdges = 2 * nnz;

  std::vector <struct Edge<lno_t, scalar_t> > edges (noEdges);
  size_type nE = 0;
  lno_t noDiagonal = 0;
  for (size_type i = 0; i < nnz; ++i){
    getline(mmf, fline);
    std::stringstream ss2 (fline);
    struct Edge<lno_t, scalar_t> tmp;
    lno_t s,d;
    scalar_t w;
    ss2 >> s >> d >> w;
    if (!transpose){
      tmp.src = s - 1;
      tmp.dst = d - 1;
      tmp.ew = w;
    }
    else {
      tmp.src = d - 1;
      tmp.dst = s - 1;
      tmp.ew = w;
    }

    if (tmp.src == tmp.dst){
      noDiagonal++;
      if (!remove_diagonal){
        edges[nE++] = tmp;
      }
      continue;
    }
    edges[nE++] = tmp;
    if (mtx_sym == 1 || symmetrize){
      struct Edge<lno_t, scalar_t> tmp2;
      tmp2.src = tmp.dst;
      tmp2.dst = tmp.src;
      tmp2.ew = tmp.ew;
      edges[nE++] = tmp2;
    }
  }

  mmf.close();

  std::sort (edges.begin(), edges.begin() + nE);

  if (transpose){
    lno_t tmp = nr;
    nr = nc;
    nc = tmp;
  }

  //idx *nv, idx *ne, idx **xadj, idx **adj, wt **wt

  *nv = nr;

  *ne = nE;
  //*xadj = new idx[nr + 1];
  md_malloc<size_type>(xadj, nr+1);
  //*adj = new idx[nE];
  md_malloc<lno_t>(adj, nE);
  //*ew = new wt[nE];
  md_malloc<scalar_t>(ew, nE);

  size_type eind = 0;
  size_type actual = 0;
  for (lno_t i = 0; i < nr; ++i){
    (*xadj)[i] = actual;
    bool is_first = true;
    while (edges[eind].src == i){
      if (is_first || !symmetrize || eind == 0 || (eind > 0 && edges[eind - 1].dst != edges[eind].dst)){

        (*adj)[actual] = edges[eind].dst;
        (*ew)[actual] = edges[eind].ew;
        ++actual;
      }
      is_first = false;
      ++eind;
    }
  }
  (*xadj)[nr] = actual;
  *ne = actual;
  return 0;
}

template <typename lno_t, typename size_type, typename scalar_t>
void read_matrix(lno_t *nv, size_type *ne,size_type **xadj, lno_t **adj, scalar_t **ew, const char *filename){

  std::string strfilename(filename);
  if (endswith(strfilename, ".mtx")){
    read_mtx (
        filename,
        nv, ne,
        xadj, adj, ew,false,false,false);
  }

  else if (endswith(strfilename, ".bin")){
    read_graph_bin(nv, ne,xadj, adj, ew, filename);
  }

  else if (endswith(strfilename, ".crs")){

    read_graph_crs(nv, ne,xadj, adj, ew, filename);
  }

  else {
    throw std::runtime_error ("Reader is not available\n");
  }
}

template <typename crsMat_t>
crsMat_t read_kokkos_crst_matrix(const char * filename_){

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
  typedef typename graph_t::entries_type::non_const_type   cols_view_t;
  typedef typename crsMat_t::values_type::non_const_type values_view_t;

  typedef typename row_map_view_t::value_type size_type;
  typedef typename cols_view_t::value_type   lno_t;
  typedef typename values_view_t::value_type scalar_t;


  lno_t nv, *adj;
  size_type *xadj, nnzA;
  scalar_t *values;
  read_matrix<lno_t, size_type, scalar_t>(
      &nv, &nnzA, &xadj, &adj, &values, filename_);

  row_map_view_t rowmap_view("rowmap_view", nv+1);
  cols_view_t columns_view("colsmap_view", nnzA);
  values_view_t values_view("values_view", nnzA);


  {
    typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
    typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
    typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

    for (lno_t i = 0; i <= nv; ++i){
      hr(i) = xadj[i];
    }

    for (size_type i = 0; i < nnzA; ++i){
      hc(i) = adj[i];
      hv(i) = values[i];
    }
    Kokkos::deep_copy (rowmap_view , hr);
    Kokkos::deep_copy (columns_view , hc);
    Kokkos::deep_copy (values_view , hv);
  }

  lno_t ncols = 0;
  KokkosKernels::Impl::kk_view_reduce_max
      <cols_view_t, typename crsMat_t::execution_space>(nnzA, columns_view, ncols);
  ncols += 1;
  
  graph_t static_graph (columns_view, rowmap_view);
  crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
  delete [] xadj; delete [] adj; delete [] values;
  return crsmat;
}


template <typename crsGraph_t>
crsGraph_t read_kokkos_crst_graph(const char * filename_){

  typedef typename crsGraph_t::row_map_type::non_const_type row_map_view_t;
  typedef typename crsGraph_t::entries_type::non_const_type   cols_view_t;

  typedef typename row_map_view_t::value_type size_type;
  typedef typename cols_view_t::value_type   lno_t;
  typedef double scalar_t ;

  lno_t nv, *adj;
  size_type *xadj, nnzA;
  scalar_t *values;
  read_matrix<lno_t, size_type, scalar_t>(
      &nv, &nnzA, &xadj, &adj, &values, filename_);

  row_map_view_t rowmap_view("rowmap_view", nv+1);
  cols_view_t columns_view("colsmap_view", nnzA);



  {
    typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
    typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);

    for (lno_t i = 0; i <= nv; ++i){
      hr(i) = xadj[i];
    }

    for (size_type i = 0; i < nnzA; ++i){
      hc(i) = adj[i];
    }
    Kokkos::deep_copy (rowmap_view , hr);
    Kokkos::deep_copy (columns_view , hc);
  }

  lno_t ncols = 0;
  KokkosKernels::Impl::kk_view_reduce_max
      <cols_view_t, typename crsGraph_t::execution_space>(nnzA, columns_view, ncols);
  ncols += 1;

  crsGraph_t static_graph (columns_view, rowmap_view, ncols);
  delete [] xadj; delete [] adj; delete [] values;
  return static_graph;
}



template <typename size_type, typename nnz_lno_t>
inline void kk_sequential_create_incidence_matrix(
    nnz_lno_t num_rows,
    const size_type *xadj,
    const nnz_lno_t *adj,
    size_type *i_adj //output. preallocated
  ){

  std::vector<size_type> c_xadj(num_rows);
  for (nnz_lno_t i = 0; i < num_rows; i++){
    c_xadj[i] = xadj[i];
  }
  int eCnt=0;
  for (nnz_lno_t i = 0; i < num_rows; i++){
    size_type begin = xadj[i];
    size_type end = xadj[i + 1];
    nnz_lno_t adjsize = end - begin;

    for (nnz_lno_t j = 0; j < adjsize; j++){
      size_type aind = j + begin;
      nnz_lno_t col = adj[aind];
      if (i < col){
        i_adj[c_xadj[i]++] = eCnt;
        i_adj[c_xadj[col]++] = eCnt++;
      }
    }
  }

  for (nnz_lno_t i = 0; i < num_rows; i++){
    if (c_xadj[i] != xadj[i+1]){
      std::cout << "i:" << i << " c_xadj[i]:" << c_xadj[i] << " xadj[i+1]:" << xadj[i+1] << std::endl;
    }
  }
}

template <typename size_type, typename nnz_lno_t>
inline void kk_sequential_create_incidence_matrix_transpose(
    const nnz_lno_t num_rows,
    const size_type num_edges,
    const size_type *xadj,
    const nnz_lno_t *adj,
    size_type *i_xadj, //output. preallocated
    nnz_lno_t *i_adj //output. preallocated
  ){

  for (nnz_lno_t i = 0; i < num_edges/2 + 1; i++){
    i_xadj[i] = i * 2;
  }
  int eCnt=0;
  for (nnz_lno_t i = 0; i < num_rows; i++){
    size_type begin = xadj[i];
    size_type end = xadj[i + 1];
    nnz_lno_t adjsize = end - begin;

    for (nnz_lno_t j = 0; j < adjsize; j++){
      size_type aind = j + begin;
      nnz_lno_t col = adj[aind];
      if (i < col){
        i_adj[eCnt++] = i;
        i_adj[eCnt++] = col;
      }
    }
  }
}


}
}

#endif
