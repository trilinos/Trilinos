#include "shylu_test_interface.hpp"


int shylu_interface::partitioning_interface(Epetra_CrsMatrix *A, Epetra_CrsMatrix *AHat,
					    Epetra_MultiVector *b, Epetra_MultiVector *bHat, 
					    ParameterList shyLUList)
{

  string partitioningPackage = Teuchos::getParameter<string>(shyLUList,
							       "Partitioning Package");
  ParameterList pList;
  pList = shyLUList.sublist(partitioningPackage + " Input");
  if(partitioningPackage.compare("Zoltan2")==0)
    {
      cout << "**ERROR** Zoltan2 not supported for Epetra Matrix" << endl;
      exit(1);
    }
  else if(partitioningPackage.compare("Isorropia")==0)
  {
    Isorropia::Epetra::Partitioner *partitioner = new
    Isorropia::Epetra::Partitioner(A, pList, false);
    partitioner->partition();
    Isorropia::Epetra::Redistributor rd(partitioner);
    rd.redistribute(*A, AHat);
    rd.redistribute(*b, bHat); 
  }
  else
  {
     cout << "----ERROR---" << endl;
     exit(1);
   }
}

#ifdef HAVE_TPETRAA
template <class ST, class LT, class GT>
int
shylu_interface::partitioning_interface(Tpetra::CrsMatrix<ST, LT, GT> *A, Tpetra::CrsMatrix<ST,LT,GT> *AHat, Tpetra::MultiVector<ST,LT, GT> *b, Tpetra::MultiVector<ST,LT,GT> *bHat, ParameterList shyLUList)
{
  cout << "Fill in zoltan2 example";
}
#endif
