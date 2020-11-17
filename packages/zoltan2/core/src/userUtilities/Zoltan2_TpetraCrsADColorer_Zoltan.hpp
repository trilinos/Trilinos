#pragma once

#include "Zoltan2_TpetraCrsADColorer.hpp"
#include "Zoltan2_TpetraCrsADColorerUtils.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include "zoltan_cpp.h"

namespace Zoltan
{

// Implementation of CrsColorer<> using Zoltan partial distance-2 coloring.
// This is distributed-parallel, but not shared.
template <typename CrsMatrixType>
class ZoltanCrsColorer
{
public:

  typedef Zoltan2_TpetraCrsADColorer<CrsMatrixType> base_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::graph_type graph_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::local_ordinal_type local_ordinal_type;
  typedef typename base_type::global_ordinal_type global_ordinal_type;
  typedef typename base_type::node_type node_type;
  typedef typename base_type::execution_space execution_space;
  typedef typename base_type::list_of_colors_type list_of_colors_type;
  typedef typename base_type::list_of_colors_host_type list_of_colors_host_type;

  // Constructor
  ZoltanCrsColorer(const Teuchos::RCP<matrix_type> &matrix_)
    : matrix(matrix_), graph(matrix_->getCrsGraph()), transpose_graph()
  {}

  // Destructor
  ~ZoltanCrsColorer() {}

  // Compute coloring data
  void
  computeColoring(
    const Teuchos::ParameterList &coloring_params, 
    const std::string &method, 
    int &num_colors,
    list_of_colors_host_t &list_of_colors_host,
    list_of_colors_t &list_of_colors);

private:

  Teuchos::RCP<const matrix_type> graph;
  Teuchos::RCP<const graph_type> graph;
  Teuchos::RCP<const graph_type> transpose_graph;

  //
  // Call-back functions for Zoltan interface
  //

  // Data for call-back functions
  struct ZoltanData
  {
    Teuchos::RCP<const graph_type> graph, trans_graph;
    Teuchos::Array<int> col_procs, trans_col_procs;

    // Set matrix graphs (trans_graph_ may be null for symmetric/symmetrized
    // problems).  Computes the remote processor ID's for each entry in the
    // graph's/trans_graph's column maps (involves communication, and so can't
    // be done inside the Zoltan call-back functions).
    void
    setGraphs(
        const Teuchos::RCP<const graph_type> &graph_,
        const Teuchos::RCP<const graph_type> &trans_graph_ = Teuchos::null)
    {
      graph = graph_;
      trans_graph = trans_graph_;
      col_procs.resize(graph->getColMap()->getNodeNumElements());
      auto gids = graph->getColMap()->getNodeElementList();

      Tpetra::LookupStatus ret = 
              graph->getRowMap()->getRemoteIndexList(gids, col_procs());
      TEUCHOS_TEST_FOR_EXCEPTION(ret != Tpetra::AllIDsPresent, std::logic_error,
                                 "Zoltan2::CrsADColorer: getRemoteIndexList() "
                                 "failed!");

      if (trans_graph != Teuchos::null)
      {
        trans_col_procs.resize(trans_graph->getColMap()->getNodeNumElements());
        gids = trans_graph->getColMap()->getNodeElementList();
        ret = trans_graph->getRowMap()->getRemoteIndexList(gids,
                                                           trans_col_procs());
        TEUCHOS_TEST_FOR_EXCEPTION(ret != Tpetra::AllIDsPresent,
                                   std::logic_error,
                                   "Zoltan2::CrsADColorer getRemoteIndexList() "
                                   "failed!");
      }
    }
  };

  // Number of vertices
  static int
  get_number_of_vertices(void *data, int *ierr);

  // Vertex IDs on this processor
  static void
  get_vertex_list(
      void *data,
      int sizeGID,
      int sizeLID,
      ZOLTAN_ID_PTR global_ids,
      ZOLTAN_ID_PTR local_ids,
      int wgt_dim,
      float *obj_wgts,
      int *ierr);

  // Get number of edges on this processor
  static int
  get_number_of_edges(
      void *data,
       int sizeGID,
       int sizeLID,
       ZOLTAN_ID_PTR global_id,
       ZOLTAN_ID_PTR local_id,
       int *ierr);

  // Get edge ids on this processor
  static void
  get_edge_list(
      void *data,
      int sizeGID,
      int sizeLID,
      ZOLTAN_ID_PTR global_id,
      ZOLTAN_ID_PTR local_id,
      ZOLTAN_ID_PTR nbor_global_ids,
      int *nbor_procs,
      int wgt_dim,
      float *ewgts,
      int *ierr);

  //
  // Call-back functions for Zoltan interface with symmetric graph
  //

  // Number of vertices
  static int
  sym_get_number_of_vertices(void *data, int *ierr);

  // Vertex IDs on this processor
  static void
  sym_get_vertex_list(
      void *data,
      int sizeGID,
      int sizeLID,
      ZOLTAN_ID_PTR global_ids,
      ZOLTAN_ID_PTR local_ids,
      int wgt_dim,
      float *obj_wgts,
      int *ierr);

  // Get number of edges on this processor
  static int
  sym_get_number_of_edges(
      void *data,
      int sizeGID,
      int sizeLID,
      ZOLTAN_ID_PTR global_id,
      ZOLTAN_ID_PTR local_id,
      int *ierr);

  // Get edge ids on this processor
  static void
  sym_get_edge_list(
      void *data,
      int sizeGID,
      int sizeLID,
      ZOLTAN_ID_PTR global_id,
      ZOLTAN_ID_PTR local_id,
      ZOLTAN_ID_PTR nbor_global_ids,
      int *nbor_procs,
      int wgt_dim,
      float *ewgts,
      int *ierr);
};

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
void
ZoltanCrsColorer<CrsMatrixType>::computeColoring(
  const Teuchos::RCP<const matrix_t> &matrix, 
  const Teuchos::ParameterList &coloring_params, 
  const std::string &method, 
  int &num_colors,
  list_of_colors_host_t &list_of_colors_host,
  list_of_colors_t &list_of_colors
{
  const bool symmetric  = coloring_params.get("Symmetric", false);
  const bool symmetrize = coloring_params.get("Symmetrize", false);

  // Get MPI communicator, and throw an exception if our comm isn't MPI
  Teuchos::RCP<const Teuchos::Comm<int>> comm = 
           this->graph->getRowMap()->getComm();
  Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm, true);

  Teuchos::RCP<const graph_type> transpose_graph;

  // Compute transpose graph
  if (!symmetric && !symmetrize)
  {
    transpose_graph = Impl::compute_transpose_graph(*(this->matrix));
  }

  // Create Zoltan for distance-2 coloring
  Zoltan *zz = new Zoltan(*(mpi_comm->getRawMpiComm()));
  zz->Set_Param("COLORING_PROBLEM", method.c_str());
  zz->Set_Param("DEBUG_LEVEL", "0");
  if (!symmetric && symmetrize)
    zz->Set_Param("GRAPH_SYMMETRIZE", "TRANSPOSE");

  // Set Zoltan params
  Teuchos::ParameterList &zoltan_params = coloring_params.sublist("Zoltan");
  for (auto p : zoltan_params)
    zz->Set_Param(p.first, Teuchos::getValue<std::string>(p.second));

  // Setup interface functions
  ZoltanData zd;
  if (symmetric || symmetrize)
  {
    zd.setGraphs(this->graph);
    zz->Set_Num_Obj_Fn(sym_get_number_of_vertices, &zd);
    zz->Set_Obj_List_Fn(sym_get_vertex_list, &zd);
    zz->Set_Num_Edges_Fn(sym_get_number_of_edges, &zd);
    zz->Set_Edge_List_Fn(sym_get_edge_list, &zd);
  }
  else
  {
    zd.setGraphs(this->graph, transpose_graph);
    zz->Set_Num_Obj_Fn(get_number_of_vertices, &zd);
    zz->Set_Obj_List_Fn(get_vertex_list, &zd);
    zz->Set_Num_Edges_Fn(get_number_of_edges, &zd);
    zz->Set_Edge_List_Fn(get_edge_list, &zd);
  }

  // Do coloring of columns with Zoltan -- we can request colors for
  // columns we don't own
  const size_t num_local_cols  = this->graph->getNodeNumCols();
  const size_t num_global_rows = this->graph->getGlobalNumRows();

  Teuchos::Array<ZOLTAN_ID_TYPE> col_gids(num_local_cols);
  auto gids = this->graph->getColMap()->getNodeElementList();

  if (symmetric || symmetrize)
    for (size_t i = 0; i < num_local_cols; ++i)
      col_gids[i] = gids[i];
  else
    for (size_t i = 0; i < num_local_cols; ++i)
      col_gids[i] = gids[i] + num_global_rows;

  list_of_colors_type my_list_of_colors("ZoltanCrsColorer::list_of_colors",
                                        num_local_cols);
  list_of_colors_host = Kokkos::create_mirror_view(my_list_of_colors);

  int num_gid_entries = 1;
  int ret = zz->Color(num_gid_entries, num_local_cols, col_gids.getRawPtr(),
                      list_of_colors_host.data());

  TEUCHOS_TEST_FOR_EXCEPTION(ret != ZOLTAN_OK, std::logic_error,
                             "Zoltan::Color returned " << ret << std::endl);

  Kokkos::deep_copy(my_list_of_colors, list_of_colors_host);
  list_of_colors = my_list_of_colors;

  const bool dump_zoltan = coloring_params.get("Dump Zoltan Data", false);
  if (dump_zoltan)
  {
    std::string zoltan_dump_file =
         coloring_params.get("Zoltan Dump File Name", "zoltan_graph.txt");
    zz->Generate_Files(zoltan_dump_file, 0, 0, 1, 0);
  }

  delete zz;

  // Compute global number of colors
  int local_num_colors = 0;
  Kokkos::parallel_reduce("ZoltanCrsColorer::find_num_colors",
      Kokkos::RangePolicy<execution_space>(0, num_local_cols),
      KOKKOS_LAMBDA(const size_t i, int &lcl_max) {
        if (my_list_of_colors[i] > lcl_max)
          lcl_max = my_list_of_colors[i];
      },
      Kokkos::Max<int>(local_num_colors));

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &local_num_colors,
                     &num_colors);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
int
ZoltanCrsColorer<CrsMatrixType>::get_number_of_vertices(void *data, int *ierr)
{
  ZoltanData *zoltan_data = static_cast<ZoltanData *>(data);
  *ierr = ZOLTAN_OK;
  return zoltan_data->graph->getNodeNumRows() +
         zoltan_data->trans_graph->getNodeNumRows();
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
void
ZoltanCrsColorer<CrsMatrixType>::get_vertex_list(
    void *data,
    int sizeGID,
    int sizeLID,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int wgt_dim,
    float *obj_wgts,
    int *ierr)
{
  assert(sizeGID == 1 && sizeLID == 1);

  ZoltanData *zoltan_data      = static_cast<ZoltanData *>(data);
  *ierr                        = ZOLTAN_OK;

  const size_t num_local_rows  = zoltan_data->graph->getNodeNumRows();
  const size_t num_local_cols  = zoltan_data->trans_graph->getNodeNumRows();
  const size_t num_global_rows = zoltan_data->graph->getGlobalNumRows();
  auto row_gids = zoltan_data->graph->getRowMap()->getNodeElementList();
  auto col_gids = zoltan_data->trans_graph->getRowMap()->getNodeElementList();

  for (size_t i = 0; i < num_local_rows; ++i)
  {
    local_ids[i] = i;
    global_ids[i] = row_gids[i];
  }
  for (size_t i = 0; i < num_local_cols; ++i)
  {
    local_ids[num_local_rows + i]  = num_local_rows + i;
    global_ids[num_local_rows + i] = num_global_rows + col_gids[i];
  }
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
int
ZoltanCrsColorer<CrsMatrixType>::get_number_of_edges(
    void *data,
    int sizeGID,
    int sizeLID,
    ZOLTAN_ID_PTR global_id,
    ZOLTAN_ID_PTR local_id,
    int *ierr)
{
  assert(sizeGID == 1 && sizeLID == 1);

  ZoltanData *zoltan_data = static_cast<ZoltanData *>(data);
  *ierr = ZOLTAN_OK;

  const size_t num_local_rows = zoltan_data->graph->getNodeNumRows();
  const ZOLTAN_ID_TYPE lid = *local_id;
  int num_edges = 0;

  if (lid < num_local_rows)
    num_edges = zoltan_data->graph->getNumEntriesInLocalRow(lid);
  else
    num_edges = 
        zoltan_data->trans_graph->getNumEntriesInLocalRow(lid - num_local_rows);

  return num_edges;
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
void
ZoltanCrsColorer<CrsMatrixType>::get_edge_list(
    void *data,
    int sizeGID,
    int sizeLID,
    ZOLTAN_ID_PTR global_id,
    ZOLTAN_ID_PTR local_id,
    ZOLTAN_ID_PTR nbor_global_ids,
    int *nbor_procs,
    int wgt_dim,
    float *ewgts,
    int *ierr)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;

  ZoltanData *zoltan_data = static_cast<ZoltanData *>(data);
  *ierr = ZOLTAN_OK;

  const size_t num_local_rows = zoltan_data->graph->getNodeNumRows();
  const size_t num_global_rows = zoltan_data->graph->getGlobalNumRows();
  const ZOLTAN_ID_TYPE lid = *local_id;

  if (lid < num_local_rows)
  {
    const int num_nbr = zoltan_data->graph->getNumEntriesInLocalRow(lid);
    const auto colMap = zoltan_data->graph->getColMap();

    ArrayView<const local_ordinal_type> lcl_ids;
    zoltan_data->graph->getLocalRowView(lid, lcl_ids);

    for (int j = 0; j < num_nbr; ++j) {
      nbor_global_ids[j] = num_global_rows 
                         + colMap->getGlobalElement(lcl_ids[j]);
      nbor_procs[j] = zoltan_data->col_procs[lcl_ids[j]];
    }
  }
  else
  {
    const int num_nbr = 
          zoltan_data->trans_graph->getNumEntriesInLocalRow(lid-num_local_rows);
    const auto colMap = zoltan_data->trans_graph->getColMap();

    ArrayView<const local_ordinal_type> lcl_ids;
    zoltan_data->trans_graph->getLocalRowView(lid - num_local_rows, lcl_ids);
    for (int j = 0; j < num_nbr; ++j)
    {
      nbor_global_ids[j] = colMap->getGlobalElement(lcl_ids[j]);
      nbor_procs[j] = zoltan_data->trans_col_procs[lcl_ids[j]];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
int
ZoltanCrsColorer<CrsMatrixType>::sym_get_number_of_vertices(
  void *data,
  int *ierr)
{
  ZoltanData *zoltan_data = static_cast<ZoltanData *>(data);
  *ierr = ZOLTAN_OK;
  return zoltan_data->graph->getNodeNumRows();
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
void
ZoltanCrsColorer<CrsMatrixType>::sym_get_vertex_list(
    void *data,
    int sizeGID,
    int sizeLID,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int wgt_dim,
    float *obj_wgts,
    int *ierr)
{
  ZoltanData *zoltan_data = static_cast<ZoltanData *>(data);
  *ierr = ZOLTAN_OK;

  const size_t num_local_rows = zoltan_data->graph->getNodeNumRows();
  auto row_gids = zoltan_data->graph->getRowMap()->getNodeElementList();
  for (size_t i = 0; i < num_local_rows; ++i)
  {
    local_ids[i]  = i;
    global_ids[i] = row_gids[i];
  }
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
int
ZoltanCrsColorer<CrsMatrixType>::sym_get_number_of_edges(
    void *data,
    int sizeGID,
    int sizeLID,
    ZOLTAN_ID_PTR global_id,
    ZOLTAN_ID_PTR local_id,
    int *ierr)
{
  ZoltanData *zoltan_data  = static_cast<ZoltanData *>(data);
  *ierr = ZOLTAN_OK;

  const ZOLTAN_ID_TYPE lid = *local_id;
  int num_edges = zoltan_data->graph->getNumEntriesInLocalRow(lid);
  return num_edges;
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
void
ZoltanCrsColorer<CrsMatrixType>::sym_get_edge_list(
    void *data,
    int sizeGID,
    int sizeLID,
    ZOLTAN_ID_PTR global_id,
    ZOLTAN_ID_PTR local_id,
    ZOLTAN_ID_PTR nbor_global_ids,
    int *nbor_procs,
    int wgt_dim,
    float *ewgts,
    int *ierr)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;

  ZoltanData *zoltan_data  = static_cast<ZoltanData *>(data);
  *ierr = ZOLTAN_OK;

  const ZOLTAN_ID_TYPE lid = *local_id;
  const int num_nbr = zoltan_data->graph->getNumEntriesInLocalRow(lid);

  ArrayView<const local_ordinal_type> lcl_ids;
  zoltan_data->graph->getLocalRowView(lid, lcl_ids);
  const auto colMap = zoltan_data->graph->getColMap();

  for (int j = 0; j < num_nbr; ++j)
  {
    nbor_global_ids[j] = colMap->getGlobalElement(lcl_ids[j]);
    nbor_procs[j]      = zoltan_data->col_procs[lcl_ids[j]];
  }
}
} // namespace Zoltan2
