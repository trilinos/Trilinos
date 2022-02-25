#pragma once

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "TpetraExt_MatrixMatrix.hpp"

#include "Zoltan2_XpetraCrsMatrixAdapter.hpp"
#include "Zoltan2_ColoringProblem.hpp"
#include "Zoltan2_ColoringSolution.hpp"

namespace Zoltan2
{

// Implementation of CrsColorer<> using Zoltan2.  Currently this is a distance-1
// algorithm, which requires forming A^T*A, and only works in serial
template <typename CrsMatrixType>
class Zoltan2CrsColorer 
{
public:
  typedef CrsMatrixType matrix_t;
  typedef typename matrix_t::crs_graph_type graph_t;
  typedef typename matrix_t::node_type node_t;
  typedef typename node_t::device_type device_t;
  typedef Kokkos::View<int *, device_t> list_of_colors_t;
  typedef typename list_of_colors_t::HostMirror list_of_colors_host_t;

  // Constructor
  Zoltan2CrsColorer(const Teuchos::RCP<matrix_t> &matrix_) 
    : matrix(matrix_), graph(matrix_->getCrsGraph()) 
  {}

  // Destructor
  ~Zoltan2CrsColorer() {}

  // Compute coloring data
  void
  computeColoring(
    Teuchos::ParameterList &coloring_params,
    int &num_colors,
    list_of_colors_host_t &list_of_colors_host,
    list_of_colors_t &list_of_colors)
  {
    // TODO:  logic for symmetric/Jacobian/Hessian a la ZoltanCrsColorer
    // User can tell us that the matrix is symmetric; 
    // otherwise, guess based on the matrix type
    // const std::string matrixType = coloring_params.get("matrixType","Jacobian");
    // const bool symmetric = coloring_params.get("symmetric",
    //                                           (matrixType=="Jacobian" ? false
    //                                                                   : true));

    // TODO:  Until Ian's code is ready...
    // TODO:  Check the logic here:  doing Partial Distance-2 via local
    // TODO:  distance-1 on J^T*J

    // Compute A^T*A where A = jac
    // We have to do this because Zoltan2 can only do distance-1 coloring
    this->matrix->setAllToScalar(1.0);

    if (!this->matrix->isFillComplete())
      this->matrix->fillComplete();

    const size_t nzpr = this->matrix->getGlobalMaxNumRowEntries();
    matrix_t C(this->matrix->getRowMap(), nzpr * nzpr); // TODO: Check this

    Tpetra::MatrixMatrix::Multiply(*(this->matrix), true,
                                   *(this->matrix), false, C);

    // Create Zoltan2 coloring problem and solve
    typedef Zoltan2::XpetraCrsMatrixAdapter<matrix_t> Z2Adapter_t;
    Z2Adapter_t z2_adapter(rcp(&C, false));

    Teuchos::ParameterList z2_params = coloring_params.sublist("Zoltan2");
    Zoltan2::ColoringProblem<Z2Adapter_t> z2_problem(&z2_adapter, &z2_params);
    z2_problem.solve();

    // Extract colors
    Zoltan2::ColoringSolution<Z2Adapter_t> *z2_solution = 
                                            z2_problem.getSolution();
    int local_num_colors = z2_solution->getNumColors();

    Teuchos::ArrayRCP<int> local_list_of_colors = z2_solution->getColorsRCP();
    const size_t len = local_list_of_colors.size();

    TEUCHOS_TEST_FOR_EXCEPTION(
        len != this->graph->getColMap()->getLocalNumElements(), std::logic_error,
        "Incorrect length of color list!");

    list_of_colors_host_t list_of_colors_tmp(
                                     local_list_of_colors.getRawPtr(), len);

    list_of_colors = list_of_colors_t("list_of_colors", len);
    list_of_colors_host = Kokkos::create_mirror_view(list_of_colors);

    Kokkos::deep_copy(list_of_colors_host, list_of_colors_tmp);
    Kokkos::deep_copy(list_of_colors, list_of_colors_host);

    // Compute global number of colors
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
             this->graph->getRowMap()->getComm();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1,
                       &local_num_colors, &num_colors);
  }

private:

  const Teuchos::RCP<matrix_t> matrix;
  const Teuchos::RCP<const graph_t> graph;
};

///////////////////////////////////////////////////////////////////////////////
// Specialization of Tpetra::Zoltan2CrsColorer for BlockCrsMatrix
// Zoltan2 does not support BlockCrs, so this just throws an error
template <typename SC, typename LO, typename GO, typename NO>
class Zoltan2CrsColorer<Tpetra::BlockCrsMatrix<SC, LO, GO, NO> > 
{
public:
  typedef Tpetra::BlockCrsMatrix<SC, LO, GO, NO> matrix_t;
  typedef typename matrix_t::crs_graph_type graph_t;

  // Constructor
  Zoltan2CrsColorer(const Teuchos::RCP<matrix_t> &matrix_) 
    : matrix(matrix_), graph(matrix_->getCrsGraph())
  {}

  // Destructor
  ~Zoltan2CrsColorer() {}

  // Compute coloring data
  void
  computeColoring(Teuchos::ParameterList &coloring_params) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Zoltan2 colorer does not support "
                               "Tpetra::BlockCrsMatrix!");
  }

private:

  const Teuchos::RCP<matrix_t> matrix;
  const Teuchos::RCP<const graph_t> graph;
};
} // namespace Zoltan2
