#ifndef SHYLU_PARTITION_INTERFACE_DEF_HPP
#define SHYLU_PARTITION_INTERFACE_DEF_HPP

#include "ShyLUCore_config.h"
#include "shylu_partition_interface_decl.hpp"

//#include "Zoltan2_config.h"
#include <Teuchos_XMLParameterListHelpers.hpp>


#ifdef HAVE_SHYLU_ZOLTAN2  
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#ifdef HAVE_SHYLU_TPETRA
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#endif

#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>


namespace ShyLU {


template <class Matrix, class Vector>
PartitionInterface<Matrix, Vector>::PartitionInterface(Matrix* inA, Teuchos::ParameterList *inpList)
{ 
  A = inA;
  pList = inpList;  
  ipart = NULL;
  ird = NULL;
#ifdef HAVE_SHLY_ZOLTAN2
  zadapter = NULL;
  zproblem = NULL;
#endif

}
template <class Matrix, class Vector>
PartitionInterface<Matrix,Vector>::~PartitionInterface()
{

  if(ipart!=NULL) delete ipart;
  if(ird!=NULL) delete ird;
  //if(zadapter!=NULL) delete zadapter;
  //if(zproblem!=NULL) delete zproblem;
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
  Teuchos::ParameterList subList = pList->sublist("Isorropia Input");
  ipart = new Isorropia::Epetra::Partitioner(A,subList,false);
  ipart->partition();
  ird = new Isorropia::Epetra::Redistributor(ipart);
  return 0;
}


#ifdef HAVE_SHYLU_ZOLTAN2
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
#endif


template <class Matrix, class Vector>
int PartitionInterface<Matrix,Vector>::partition()
{
#ifdef HAVE_SHYLU_ZOLTAN2
  return partitionZoltan2();
#else
  return 1;
#endif
}
template <>
int PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::partition()
{
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
  if(partitioningPackage.compare("Isorropia") == 0)
    {
      return partitionIsorropia();
    }
#ifdef HAVE_SHYLU_ZOLTAN2
  else if (partitioningPackage.compare("Zoltan2") == 0)
    {
      return partitionZoltan2();
    }
#endif
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

#if defined(HAVE_ZOLTAN2_PARMETIS) || defined(HAVE_ZOLTAN2_SCOTCH)
  if(partitioningPackage.compare("Zoltan2") == 0)
    {
      zadapter->applyPartitioningSolution(*A, B, zproblem->getSolution());     
    }
#else
B = NULL;

#endif
  return B;
}
template < >
Epetra_CrsMatrix* PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::reorderMatrix()
{
  Epetra_CrsMatrix *B = NULL;
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
  if(partitioningPackage.compare("Isorropia") == 0)
    {
      ird->redistribute(*A, B);  
    }
#if defined(HAVE_ZOLTAN2_PARMETIS) || defined(HAVE_ZOLTAN2_SCOTCH)
  else if (partitioningPackage.compare("Zoltan2") == 0)
    {
      zadapter->applyPartitioningSolution(*A, B, zproblem->getSolution());     
    }
#endif
  return B;

}
template <class Matrix, class Vector>
Vector* PartitionInterface<Matrix, Vector>::reorderVector(Vector* x)
{
  Vector *b = NULL;
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
#if defined(HAVE_ZOLTAN2_PARMETIS) || defined(HAVE_ZOLTAN2_SCOTCH) 
      Teuchos::RCP<Vector> rx(x);
      Zoltan2::XpetraMultiVectorAdapter<Vector> tempVecAdapter(rx);
      tempVecAdapter.applyPartitioningSolution(*x, b, zproblem->getSolution());
#endif
  return b;
}
template < >
Epetra_MultiVector* PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::reorderVector(Epetra_MultiVector* x )
{
  Epetra_MultiVector *b = NULL;
  string partitioningPackage = Teuchos::getParameter<string>(*pList, "Partitioning Package");
  if(partitioningPackage.compare("Isorropia") == 0)
    {
      ird->redistribute(*x, b);  
    }
#if defined(HAVE_ZOLTAN2_PARMETIS) || defined(HAVE_ZOLTAN2_SCOTCH)
  else if (partitioningPackage.compare("Zoltan2") == 0)
    {
      Teuchos::RCP<Epetra_MultiVector> rx(x);
      Zoltan2::XpetraMultiVectorAdapter<Epetra_MultiVector> tempVecAdapter(rx);
      tempVecAdapter.applyPartitioningSolution(*x, b, zproblem->getSolution());
    }
#endif
  return b;
}

} 
#endif
