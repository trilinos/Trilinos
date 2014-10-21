#include "shylu_partition_interface.hpp"

#include "Zoltan2_config.h"
#include <Teuchos_XMLParameterListHelpers.hpp>


#if defined(HAVE_ZOLTAN2_PARMETIS) || defined(HAVE_ZOLTAN2_SCOTCH) 
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#ifdef HAVE_TPETRA_EPETRA
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#endif

#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>


//using namespace ShyLU;


template <class Matrix, class Vector>
PartitionInterface<Matrix, Vector>::PartitionInterface(Matrix* inA, Teuchos::ParameterList *inpList)
{ 
  A = inA;
  pList = inpList;  
#ifdef HAVE_ISORROPIA
  ipart = NULL;
  ird = NULL;
#endif
  zadapter = NULL;
  zproblem = NULL;
}
template <class Matrix, class Vector>
PartitionInterface<Matrix,Vector>::~PartitionInterface()
{
#ifdef HAVE_ISORROPIA
  if(ipart!=NULL) delete ipart;
  if(ird!=NULL) delete ird;
#endif
  if(zadapter!=NULL) delete zadapter;
  if(zproblem!=NULL) delete zproblem;
}
template <class Matrix, class Vector>
int PartitionInterface<Matrix, Vector>::partitionIsorropia()
{
  cout << " Not Supported \n";
  return 1;
}

template <>
int PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::partitionIsorropia()
{
  ParameterList subList = pList->sublist("Isorropia Input");
  ipart = new Isorropia::Epetra::Partitioner(A,subList,false);
  ipart->partition();
  ird = new Isorropia::Epetra::Redistributor(ipart);
  return 0;
}

template <class Matrix, class Vector>
int PartitionInterface<Matrix, Vector>::partitionZoltan2()
{
 
  ParameterList subList = pList->sublist("Zoltan2 Input");
  Teuchos::RCP<Matrix> rA(A);
  zadapter = new Zoltan2::XpetraCrsMatrixAdapter<Matrix, Vector>(rA);
  zproblem = new Zoltan2::PartitioningProblem<Zoltan2::XpetraCrsMatrixAdapter <Matrix, Vector> >(zadapter, &subList);
  zproblem->solve();
  return 0;
}


template <class Matrix, class Vector>
int PartitionInterface<Matrix,Vector>::partition()
{
  return partitionZoltan2();
}
template <>
int PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::partition()
{
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
  if(partitioningPackage.compare("Isorropia") == 0)
    {
      return partitionIsorropia();
    }
  else if (partitioningPackage.compare("Zoltan2") == 0)
    {
      return partitionZoltan2();
    }
  else
    {
      cout << "**Error**: Paritioning package selected is not supported\n";
    }
  return 1;
}

template <class Matrix, class Vector>
Matrix* PartitionInterface<Matrix,Vector>::reorderMatrix()
{
  Matrix *B;
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
 
  if (partitioningPackage.compare("Zoltan2") == 0)
    {
      zadapter->applyPartitioningSolution(*A, B, zproblem->getSolution());     
    }
  return B;
}
template < >
Epetra_CrsMatrix* PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::reorderMatrix()
{
  Epetra_CrsMatrix *B;
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
  if(partitioningPackage.compare("Isorropia") == 0)
    {
      ird->redistribute(*A, B);  
    }
  else if (partitioningPackage.compare("Zoltan2") == 0)
    {
      zadapter->applyPartitioningSolution(*A, B, zproblem->getSolution());     
    }
  return B;

}
template <class Matrix, class Vector>
Vector* PartitionInterface<Matrix, Vector>::reorderVector(Vector* x)
{
  Vector *b;
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
 
      Teuchos::RCP<Vector> rx(x);
      Zoltan2::XpetraMultiVectorAdapter<Vector> tempVecAdapter(rx);
      tempVecAdapter.applyPartitioningSolution(*x, b, zproblem->getSolution());
    
  return b;
}
template < >
Epetra_MultiVector* PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::reorderVector(Epetra_MultiVector* x )
{
  Epetra_MultiVector *b;
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
  if(partitioningPackage.compare("Isorropia") == 0)
    {
      ird->redistribute(*x, b);  
    }
  else if (partitioningPackage.compare("Zoltan2") == 0)
    {
      Teuchos::RCP<Epetra_MultiVector> rx(x);
      Zoltan2::XpetraMultiVectorAdapter<Epetra_MultiVector> tempVecAdapter(rx);
      tempVecAdapter.applyPartitioningSolution(*x, b, zproblem->getSolution());
    }
  return b;
}
template class PartitionInterface<Epetra_CrsMatrix,Epetra_MultiVector>;


typedef double scalar_type;
  typedef int local_o_type;
  typedef int global_o_type;
  typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_o_type, global_o_type, node_type> Matrix_t;
typedef Tpetra::MultiVector<scalar_type, local_o_type, global_o_type, node_type> Vector_t;
template class PartitionInterface<Matrix_t, Vector_t>;

 
