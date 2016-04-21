
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <vector>


#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Timer.hpp>

#ifndef _KOKKOSKERNELSGRAPHUTILS_HPP
#define _KOKKOSKERNELSGRAPHUTILS_HPP



namespace KokkosKernels{
namespace Experimental{



namespace Graph{

namespace Utils{


template <typename stype>
void md_malloc(stype **arr, size_t n, std::string alloc_str = ""){
  *arr = new stype[n];
  if (*arr == NULL){
    std::cerr << "Memory Allocation Problem " << alloc_str << std::endl;
    exit(1);
  }
}

template <typename idx, typename wt>
struct Edge{
  idx src;
  idx dst;
  wt ew;
  bool operator<(const Edge <idx,wt> & a) const
  {
    return (this->src < a.src) || (this->src == a.src && this->dst < a.dst);
  }
};


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

template <typename idx, typename wt>
void convert_edge_list_to_csr (idx nv, idx ne, idx *srcs, idx *dests, wt *ew, idx *xadj, idx *adj, wt *crs_ew){

  std::vector <struct Edge<idx, wt> > edges (ne);
  for(idx i = 0; i < ne; ++i){
    edges[i].src = srcs[i];
    edges[i].dst = dests[i];
    edges[i].ew = ew[i];
  }
  std::sort (edges.begin(), edges.begin() + ne);

  idx eind = 0;
  for (idx i = 0; i < nv; ++i){
    (xadj)[i] = eind;
    while (edges[eind].src == i){
      (adj)[eind] = edges[eind].dst;
      (*crs_ew)[eind] = edges[eind].ew;
      ++eind;
    }
  }
  xadj[nv] = eind;

}

template <typename idx, typename wt>
void write_graph_bin(idx nv, idx ne,idx *xadj, idx *adj, wt *ew, char *filename){
  std::ofstream myFile (filename, std::ios::out | std::ios::binary);
  myFile.write((char *) &nv, sizeof(idx));
  myFile.write((char *) &ne, sizeof(idx));
  myFile.write((char *) xadj, sizeof(idx) * (nv + 1));
  myFile.write((char *) adj, sizeof(idx) * (ne));
  myFile.write((char *) ew, sizeof(wt) * (ne));
  myFile.close();
}

template <typename idx, typename wt>
void read_graph_bin(idx *nv, idx *ne,idx **xadj, idx **adj, wt **ew, char *filename){

  std::cout << "filename:" << filename << std::endl;
  std::ifstream myFile (filename, std::ios::in | std::ios::binary);

  myFile.read((char *) nv, sizeof(idx));
  myFile.read((char *) ne, sizeof(idx));
  md_malloc<idx>(xadj, *nv+1);
  md_malloc<idx>(adj, *ne);
  md_malloc<wt> (ew, *ne);
  myFile.read((char *) *xadj, sizeof(idx) * (*nv + 1));
  myFile.read((char *) *adj, sizeof(idx) * (*ne));
  myFile.read((char *) *ew, sizeof(wt) * (*ne));
  myFile.close();
  std::cout << "nr:" << *nv << " nnz:" << *ne << std::endl;
}


template <typename idx, typename wt>
int read_mtx (
    char *fileName,
    idx *nv, idx *ne,
    idx **xadj, idx **adj, wt **ew,
    bool symmetrize = false, bool remove_diagonal = true,
    bool transpose = false){

  std::ifstream mmf (fileName, std::ifstream::in);
  if (!mmf.is_open()) {std::cerr << "File cannot be opened" << std::endl; return 1;}
  std::cout << "Reading MTX file with name:" << fileName << std::endl;
  std::string fline = "";
  getline(mmf, fline);

  if (fline.size() < 2 || fline[0] != '%' || fline[1] != '%'){
    std::cerr << "Invalid MM file. Line-1" << std::endl;
    std::cerr << fline << std::endl;
    return 1;
  }


  int mtx_object = 0; // 0- matrix, 1-vector
  int mtx_format = 0; // 0- coordinate; 1- array
  int mtx_field = 0; //0-real, 1-double, 2-complex, 3- integer, 4-pattern
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
    mtx_field = 0;
  }
  else if (fline.find("double") != std::string::npos){
    mtx_field = 1;
  }
  else if (fline.find("complex") != std::string::npos){
    mtx_field = 2;
  }
  else if (fline.find("integer") != std::string::npos){
    mtx_field = 3;
  }
  else if (fline.find("pattern") != std::string::npos){
    mtx_field = 4;
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

  if (mtx_object == 1) {std::cerr << "VECTOR TYPE NOT HANDLED YET"<< std::endl; return (1); }
  if (mtx_format == 1) {std::cerr << "ARRAY TYPE NOT HANDLED YET"<< std::endl; return (1); }
  if (!symmetrize && (mtx_sym == 2 || mtx_sym == 3)) {std::cerr << "SKEW-SYMMETRIC and HERMITIAN TYPE NOT HANDLED YET"<< std::endl; return (1);}


  while(1){
    getline(mmf, fline);
    if(fline[0] != '%') break;
  }
  std::stringstream ss (fline);
  idx nr = 0, nc = 0, nnz = 0;

  ss >> nr >> nc >> nnz;


  //if (nr != nc) {std::cerr << "NON-SQUARE MATRIX TYPE NOT HANDLED YET"<< std::endl; return (1); }
  idx noEdges = nnz;
  if (mtx_sym == 1 || symmetrize) noEdges = 2 * nnz;
  std::cout << "nr:" << nr << " nc:" << nc << " nnz:" << nnz << " noEdges:" << noEdges << std::endl;

  std::vector <struct Edge<idx, wt> > edges (noEdges);
  idx nE = 0, noDiagonal = 0;
  for (idx i = 0; i < nnz; ++i){
    getline(mmf, fline);
    std::stringstream ss (fline);
    struct Edge<idx, wt> tmp;
    idx s,d;
    wt w;
    ss >> s >> d >> w;
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
      struct Edge<idx, wt> tmp2;
      tmp2.src = tmp.dst;
      tmp2.dst = tmp.src;
      tmp2.ew = tmp.ew;
      edges[nE++] = tmp2;
    }
  }

  mmf.close();

  std::sort (edges.begin(), edges.begin() + nE);

  if (transpose){
    idx tmp = nr;
    nr = nc;
    nc = tmp;
  }
  std::cout   << "MTX READ" << std::endl
    << "NV:" << nr << " NNZ:" << nnz << std::endl
    << "No Diagonals:" << noDiagonal << " no edges:" << nE << std::endl;

  //idx *nv, idx *ne, idx **xadj, idx **adj, wt **wt

  *nv = nr;

  *ne = nE;
  //*xadj = new idx[nr + 1];
  md_malloc<idx>(xadj, nr+1);
  //*adj = new idx[nE];
  md_malloc<idx>(adj, nE);
  //*ew = new wt[nE];
  md_malloc<wt>(ew, nE);

  idx eind = 0;
  idx actual = 0;
  for (idx i = 0; i < nr; ++i){
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

template <typename idx, typename color_type>
void count_colors(idx nv, color_type *colors, idx *histogram, color_type *numColors){
  color_type nc = 0;

  for (idx i = 0; i < nv + 1; ++i ){
    histogram[i] = 0;
  }
  for (idx i = 0; i < nv; ++i ){
    if (histogram[colors[i]]++ == 0) nc++;
  }
  *numColors = nc;
}



template <typename idx, typename color_type, typename MyExecSpace>
void count_colors(
    idx nv,
    typename Kokkos::View<color_type * , MyExecSpace> colors,
    typename Kokkos::View<idx *, MyExecSpace> histogram,
    color_type *numColors){
  typedef Kokkos::View<color_type * , MyExecSpace> color_array_type;
  typedef Kokkos::View<idx *, MyExecSpace> idx_array_type;
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  color_type nc = 0;

  struct count_colors{
    color_array_type colors_;
    idx_array_type histogram_;
    idx nvertex_;
    idx increment;
    count_colors(color_array_type cs, idx_array_type hs_, idx nv_):
      colors_(cs), histogram_(hs_), nvertex_(nv_), increment(1){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, color_type &numCols) const {

      idx val = Kokkos::atomic_fetch_add<idx>(&(histogram_[colors_[ii]]), increment);
      if (val == 0) numCols++;
    }
  };

  Kokkos::parallel_reduce(
      my_exec_space(0, nv),
      count_colors(colors, histogram, nv),
      nc);
  Kokkos::fence();
  *numColors = nc;
}

template <typename idx, typename color_type>
int is_coloring_valid(idx nv, idx *xadj, idx *adj, color_type* colors){
  for (idx i = 0; i < nv; ++i) {
    idx nb = xadj[i];
    idx ne = xadj[i+1];
    color_type scol = colors[i];
    for(idx j = nb; j < ne; ++j) {
      idx d = adj[j];
      color_type dcol = colors[d];
      if (scol == dcol){
      std::cout << "#############################ERROR##########################################"<< std::endl;
      std::cout << "s:" << i << " d:" << d << " scol:" << scol << " dcol:" << dcol << std::endl;
      std::cout << "############################################################################"<< std::endl;
        return 0;
      }
    }
  }
  return 1;
}


template <typename idx>
double get_std_degree(idx nv, idx ne,idx *xadj, idx *adj){
  double avg_degree = ne / (double) (nv);

  double std = 0;
  for (idx i = 0; i < nv; ++i){
    double degree = xadj[i+1] - xadj[i];
    degree = degree - avg_degree;
    std += degree * degree;
  }
  return std / nv;
}

template <typename idx>
void get_warp_properties(
    long & total_warp_work, long & total_wasted_cycles, long &actual_work,
    idx nv, idx ne, idx *xadj, idx *adj, idx chunksize = 8, idx warp_size = 32){
  total_warp_work = 0;
  total_wasted_cycles = 0;
  actual_work = 0;
  for (idx i = 0; i < nv; i += warp_size * chunksize){

    idx cutoff = i + warp_size * chunksize;
    if (nv < cutoff) cutoff = nv;

    for (idx k = 0; k < chunksize; ++k){
      idx max_thread_work = 0;

      for (idx j = i + k; j < cutoff ; j+= chunksize){
        idx thread_work = xadj[j+1 ] - xadj[j];
        actual_work += thread_work;
        if (thread_work > max_thread_work) max_thread_work = thread_work;
      }

      idx wasted_cycles = 0;
      for (idx j = i + k; j < cutoff ; j+= chunksize){
        idx thread_work = xadj[j+1] - xadj[j];
        wasted_cycles += max_thread_work - thread_work;
      }

      total_warp_work += max_thread_work * warp_size;
      total_wasted_cycles += wasted_cycles;
    }
  }
}

template <typename idx>
double get_thread_imbalance(
    idx nv, idx ne, idx *xadj, idx *adj, idx num_threads = 228){

  idx block_size = ceil (nv / double (num_threads));

  idx max_work = 0;
  for (idx i = 0; i < num_threads; ++i){
    idx work_begin = i * block_size;
    idx work_end = (i + 1) * block_size;
    if (work_end > nv) work_end = nv;

    idx thread_work = xadj [work_end] - xadj[work_begin];
    if (max_work < thread_work){
      max_work = thread_work;
    }
  }

  double expected_work = ne / double(num_threads) + 0.5;
  return max_work / expected_work;

}



template <typename idx, typename wt>
int is_symmetric_graph(idx nv, idx *xadj, idx *adj, wt *ew = NULL){
  for (idx i = 0; i < nv; ++i){
    idx nb = xadj[i];
    idx ne = xadj[i + 1];
    for (idx j = nb; j < ne; ++j){
      idx n = adj[j];
      bool am_i_ns_neighbor = false;

      wt inew = 1;
      if (ew){
        inew = ew[j];
      }

      idx nnb = xadj[n];
      idx nne = xadj[n + 1];


      for (idx nj = nnb; nj < nne; ++nj){
        idx nn = adj[nj];

        if (nn == i){
          if (ew){
            wt inewr = ew[nj];
            if (inew != inewr) {
              am_i_ns_neighbor = false;
            }

          }
          else {
            am_i_ns_neighbor = true;
          }

          break;
        }


      }
      if (!am_i_ns_neighbor){
        std::cout << "not symetric i:" << i << " n:" << n << " bot not reverse." << std::endl;
        return 0;
      }
    }
  }
  return 1;
}
}
}
}
}
#endif
