#include "shylu_directsolver_interface_decl.hpp"
#include "shylu_directsolver_interface_def.hpp"


namespace ShyLU{

  template class DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>;

#ifdef HAVE_SHYLUCORE_TPETRA
  //Tpetra
  typedef KokkosClassic::DefaultNode::DefaultNodeType node;
 
  template class DirectSolverInterface<
    Tpetra::CrsMatrix<float, int, int, node>,
    Tpetra::MultiVector<float,int,int, node> >;

  template class DirectSolverInterface<
    Tpetra::CrsMatrix<double, int, int, node> ,
    Tpetra::MultiVector<double, int, int, node> >;
#endif


} //end namespace
