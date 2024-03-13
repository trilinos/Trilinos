// Tpetra provides distributed sparse linear algebra
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_Export_decl.hpp>

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include <Ifpack2_AdditiveSchwarz_decl.hpp>

namespace MueLu {

// Define default types
typedef double scalar_type;
typedef int local_ordinal_type;
typedef int global_ordinal_type;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType node_type;

typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> operator_type;

typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> crs_matrix_type;
typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;
typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> multivector_type;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> driver_map_type;

typedef MueLu::TpetraOperator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> muelu_tpetra_operator_type;
typedef MueLu::Utilities<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MueLuUtilities;

typedef Xpetra::Matrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> xpetra_matrix;

typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> precond_type;

struct DomainPartitioning {
  global_ordinal_type nx = 0;
  global_ordinal_type ny = 0;
  global_ordinal_type nz = 0;

  global_ordinal_type bricksize_x = 0;
  global_ordinal_type bricksize_y = 0;
  global_ordinal_type bricksize_z = 0;
};

class AdditiveVariant : public operator_type {
 public:
  AdditiveVariant(RCP<crs_matrix_type>, RCP<multivector_type>, DomainPartitioning);

  void apply(const multivector_type &, multivector_type &, Teuchos::ETransp, scalar_type, scalar_type) const;

  bool hasTransposeApply() const { return false; }

  RCP<const driver_map_type> getDomainMap() const { return DomainMap_; };

  RCP<const driver_map_type> getRangeMap() const { return RangeMap_; };

 private:
  DomainPartitioning domain_;

  RCP<const Teuchos::Comm<int> > GlobalComm_;

  RCP<multivector_type> coords_;

  RCP<const driver_map_type> DomainMap_;

  RCP<const driver_map_type> RangeMap_;

  // MueLu Preconditioner to store the smoother at the fine level
  RCP<muelu_tpetra_operator_type> B_fine_ = Teuchos::null;

  // RCP<Ifpack2::AdditiveSchwarz< row_matrix_type, precond_type > > B_DD_;

  // MueLu Preconditioner to store
  RCP<muelu_tpetra_operator_type> B_coarse_ = Teuchos::null;

  void AdditiveFineSmoother(RCP<crs_matrix_type>);

  void AdditiveCoarseSolver(RCP<crs_matrix_type>);
};

}  // namespace MueLu
