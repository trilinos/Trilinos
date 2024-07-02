// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//#include <Teuchos_ConfigDefs.hpp>
//#include <Ifpack2_ConfigDefs.hpp>

//#include <Ifpack2_Version.hpp>
#include <iostream>
#include <fstream>

#include <Ifpack2_Relaxation.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <vector>

using Teuchos::RCP;
using Teuchos::rcp;
using namespace Teuchos;

template <typename stype>
void md_malloc(stype **arr, size_t n, std::string alloc_str = ""){
  *arr = new stype[n];
  if (*arr == NULL){
    std::cerr << "Memory Allocation Problem " << alloc_str << std::endl;
    exit(1);
  }
}

/**
 * Reads binary graph, the graphs are generated from mtx file using TpetraKernels_conMTX2BIN.exe.
 * See the convert TpetraKernels_conMTX2BIN.exe executable in KokkosKernels package.
 *  /perf_test/graph/TpetraKernels_conMTX2BIN.exe
 *
 */
template <typename idx, typename wt>
void read_graph_bin(idx *nv, idx *ne,idx **xadj, idx **adj, wt **ew, char *filename){


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
}


typedef int zlno_t;
typedef int zgno_t;
typedef double zscalar_t;

int main( int argc, char* argv[] )
{
  if (argc < 2){
    std::cerr << "Usage:" << argv[0] << " input_bin_file [max-num-iter=500]" << std::endl;
    exit(1);
  }
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
  int rank = tcomm->getRank();
  int max_num_iter = 500;
  if (argc >= 3){
    max_num_iter = atoi(argv[2]);
  }
  double min_norm = 0.0000001;
  if (argc == 4){
    min_norm = atof(argv[3]);
  }

  std::cout << "max_num_iter:" << max_num_iter << " min_norm:" << min_norm << std::endl;

  int nr = 0, ne = 0;
  int *xadj = NULL, *adj= NULL;
  double *vals= NULL;

  if (rank == 0){
    read_graph_bin<int, double> (
        &nr, &ne, &xadj, &adj, &vals, argv[1]);
  }

  
  tcomm->broadcast (0, sizeof(int), (char *)&nr);
  tcomm->broadcast (0, sizeof(int), (char *)&ne);
  std::cout << "nr:" << nr << " ne:" << ne <<std::endl;
  if (rank != 0){
    xadj = new int [nr+1];
    adj = new int [ne];
    vals = new double [ne];
  }
  tcomm->broadcast (0, sizeof(int) * (nr + 1), (char *)xadj);
  tcomm->broadcast (0, sizeof(int) * (ne), (char *)adj);
  tcomm->broadcast (0, sizeof(double) * (ne), (char *)vals);


  int max_num_elements = 0;
  for (zlno_t lclRow = 0; lclRow < nr; ++lclRow) {
    int begin = xadj[lclRow];
    int end = xadj[lclRow + 1];
    if (end - begin > max_num_elements) max_num_elements = end - begin;
  }

  typedef Tpetra::Map<>::node_type znode_t;
  typedef Tpetra::Map<zlno_t, zgno_t, znode_t> map_t;
  size_t numGlobalCoords = nr;
  RCP<const map_t> map = rcp (new map_t (numGlobalCoords, 0, tcomm));
  RCP<const map_t> domainMap = map;
  RCP<const map_t> rangeMap = map;

  typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> tcrsMatrix_t;
  typedef Tpetra::RowMatrix<zscalar_t, zlno_t, zgno_t, znode_t> row_matrix_type;
  RCP<tcrsMatrix_t> TpetraCrsMatrix(new tcrsMatrix_t (map, 0));
  const zlno_t numMyElements = map->getLocalNumElements ();
  const zgno_t myBegin = map->getGlobalElement (0);

  //std::vector<zgno_t> tmp_indices(max_num_elements);
  //std::vector<zscalar_t> tmp_vals(max_num_elements);
  zgno_t  *tmp_indices = new zgno_t[max_num_elements];
  zscalar_t *tmp_vals = new zscalar_t [max_num_elements];
  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    zgno_t begin = xadj[gblRow];
    zgno_t end = xadj[gblRow + 1];


    for (int i = begin; i < end; ++i) {
      tmp_indices[i - begin] = adj[i];
      tmp_vals[i - begin] = vals[i];

    }
    const Teuchos::ArrayView< const zgno_t > indices(tmp_indices, end-begin);
    const Teuchos::ArrayView< const zscalar_t > values(tmp_vals, end-begin);
    TpetraCrsMatrix->insertGlobalValues (gblRow, indices, values);
  }
  TpetraCrsMatrix->fillComplete ();


  
  delete []xadj; delete [] adj; delete []vals;
  delete []tmp_indices; delete [] tmp_vals;

  typedef Tpetra::Vector<zscalar_t, zlno_t, zgno_t, znode_t> vec_type;

  vec_type X_wanted (domainMap);
  vec_type Y_result (rangeMap);
  X_wanted.randomize ();
  TpetraCrsMatrix->apply(X_wanted, Y_result);
  vec_type X_seek (domainMap);
  X_seek.putScalar(0);

  ParameterList params, params_mt;
  params.set ("relaxation: type", "Symmetric Gauss-Seidel");
  params.set ("relaxation: sweeps", 1);
  params.set ("relaxation: zero starting solution", false);

  params_mt.set ("relaxation: type", "MT Symmetric Gauss-Seidel");
  params_mt.set ("relaxation: sweeps", 1);
  params_mt.set ("relaxation: zero starting solution", false);


  Ifpack2::Relaxation<row_matrix_type> prec_mt (TpetraCrsMatrix);
  Ifpack2::Relaxation<row_matrix_type> prec (TpetraCrsMatrix);


  prec.setParameters (params);

  prec.initialize () ;
  prec.compute () ;

  prec_mt.setParameters (params_mt);

  prec_mt.initialize () ;

  prec_mt.compute () ;

  X_seek.putScalar(0);
  for (int i = 0; i < max_num_iter; ++i){
    prec_mt.apply (Y_result, X_seek);
    vec_type X_diff(X_seek, Teuchos::Copy);
    X_diff.update (1, X_wanted, -1);
    double normInf = X_diff.normInf ();
    if (i % 10 == 0)
    std::cout << "i:" << i << " norm:" << normInf << std::endl;
    if (normInf < min_norm) break;
  }
  std::cout << "MT Flops:" << prec_mt.getApplyFlops ()
                    << " App Time:" <<  prec_mt.getApplyTime ()
                    << " Comp Time:" <<  prec_mt.getComputeTime ()
                    << " Init Time:" <<  prec_mt.getInitializeTime ()
                      << std::endl;


  X_seek.putScalar(0);
  for (int i = 0; i < max_num_iter; ++i){
    prec.apply (Y_result, X_seek);
    vec_type X_diff(X_seek, Teuchos::Copy);
    X_diff.update (1, X_wanted, -1);
    double normInf = X_diff.normInf ();
    if (i % 10 == 0)
    std::cout << "i:" << i << " norm:" << normInf << std::endl;
    if (normInf < min_norm) break;
  }
  std::cout << "Sequential Flops:" << prec.getApplyFlops ()
                << " App Time:" <<  prec.getApplyTime ()
                << " Comp Time:" <<  prec.getComputeTime ()
                << " Init Time:" <<  prec.getInitializeTime ()
                  << std::endl;
  return 0;
}
