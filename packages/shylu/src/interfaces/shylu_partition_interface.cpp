#include "shylu_partition_interface_decl.hpp"
#include "shylu_partition_interface_def.hpp"


namespace ShyLU{

  //Epetra
  template class PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>;

#ifdef HAVE_SHYLU_TPETRA
  //Tpetra
  typedef KokkosClassic::DefaultNode::DefaultNodeType node;
 
  template class PartitionInterface<
    Tpetra::CrsMatrix<float, int, int, node>,
    Tpetra::MultiVector<float,int,int, node> >;

  template class PartitionInterface<
    Tpetra::CrsMatrix<double, int, int, node> ,
    Tpetra::MultiVector<double, int, int, node> >;
#endif
  /*
  template class PartitionInterface<
    Tpetra::CrsMatrix<float, int, long, node> ,
       Tpetra::MultiVector<double, int, long, node> >;
  */

}//end namespace
