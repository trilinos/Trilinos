#include "shylu_directsolver_interface_decl.hpp"
#include "shylu_directsolver_interface_def.hpp"
#include "Tpetra_ConfigDefs.hpp"

namespace ShyLU {

  template class DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>;

#if defined(HAVE_SHYLUCORE_TPETRA)

#if defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)
  template class DirectSolverInterface<
    Tpetra::CrsMatrix<float, int, int>,
    Tpetra::MultiVector<float, int, int> >;
#endif // defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
  template class DirectSolverInterface<
    Tpetra::CrsMatrix<double, int, int> ,
    Tpetra::MultiVector<double, int, int> >;
#endif // defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)

#endif // defined(HAVE_SHYLUCORE_TPETRA)

} // namespace ShyLU
