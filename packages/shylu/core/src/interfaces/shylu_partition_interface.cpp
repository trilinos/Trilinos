#include "shylu_partition_interface_decl.hpp"
#include "shylu_partition_interface_def.hpp"


namespace ShyLU{

  //Epetra
  template class PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>;

#ifdef HAVE_SHYLUCORE_TPETRA
  //Tpetra
  //typedef KokkosClassic::DefaultNode::DefaultNodeType node;
 
  template class PartitionInterface<
    Tpetra::CrsMatrix<float, int, int>,
    Tpetra::MultiVector<float,int,int> >;

  template class PartitionInterface<
    Tpetra::CrsMatrix<double, int, int> ,
    Tpetra::MultiVector<double, int, int> >;
#endif
  /*
  template class PartitionInterface<
    Tpetra::CrsMatrix<float, int, long, node> ,
       Tpetra::MultiVector<double, int, long, node> >;
  */

}//end namespace
