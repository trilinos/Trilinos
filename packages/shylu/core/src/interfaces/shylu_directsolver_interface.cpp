#include "shylu_directsolver_interface_decl.hpp"
#include "shylu_directsolver_interface_def.hpp"


namespace ShyLU{

  template class DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>;



#if defined(HAVE_SHYLUCORE_TPETRA)


   
  typedef Tpetra::Details::DefaultTypes::node_type node_type;

  template class DirectSolverInterface<
    Tpetra::CrsMatrix<float, int, int, node_type>,
    Tpetra::MultiVector<float,int,int, node_type> >;

  template class DirectSolverInterface<
    Tpetra::CrsMatrix< double, int, int, node_type> ,
    Tpetra::MultiVector< double, int, int, node_type> >;


#endif


} //end namespace
