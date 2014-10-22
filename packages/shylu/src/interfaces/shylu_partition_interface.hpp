#ifndef SHYLU_PARTITION_INTERFACE_HPP
#define SHYLU_PARTITION_INTERFACE_HPP

//shylu
#include <shylu.h>
#include <shylu_util.h>

//Epetra
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

//Isorropia
#include <Isorropia_config.h>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraProber.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

//Tperta
#include <Tpetra_config.h>
#ifdef HAVE_TPETRA_EPETRA
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#endif

#include <Zoltan2_config.h>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>


#include <Teuchos_XMLParameterListHelpers.hpp>


//namespace ShyLU{

template <class Matrix, class Vector>
class PartitionInterface
{

public:
  ~PartitionInterface();
  PartitionInterface(Matrix* inA, Teuchos::ParameterList* pList);
  int partition();
  Matrix* reorderMatrix();
  Vector* reorderVector(Vector* x);

private:
  int partitionIsorropia();
  int partitionZoltan2();
  Teuchos::ParameterList* pList;
  Matrix* A;
  ///other handlers needed by zoltan2 and isorropia
  Zoltan2::XpetraCrsMatrixAdapter<Matrix,Vector> *zadapter;
  Zoltan2::PartitioningProblem<Zoltan2::XpetraCrsMatrixAdapter<Matrix,Vector> > *zproblem;
  Isorropia::Epetra::Partitioner *ipart;
  Isorropia::Epetra::Redistributor *ird;


};


#endif
/*
template <> 
class PartitionInterface <Epetra_CrsMatrix, Epetra_MultiVector> 
{
public:
  int partition();
  Epetra_CrsMatrix* reorderMatrix();
  Epetra_MultiVector* reorderVector(Epetra_MultiVector* x);

  };

*/
//}




