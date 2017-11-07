//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

/** \file shylu_partition_interface_def.hpp

    \brief Epetra/Tpetra templated interface for calls to Zoltan(Isorropia)/Zoltans

    \author Joshua Dennis Booth
*/

#ifndef SHYLU_PARTITION_INTERFACE_DEF_HPP
#define SHYLU_PARTITION_INTERFACE_DEF_HPP

#include "ShyLU_DDCore_config.h"
#include "shylu_partition_interface_decl.hpp"

//#include "Zoltan2_config.h"
#include <Teuchos_XMLParameterListHelpers.hpp>


#ifdef HAVE_SHYLU_DDCORE_ZOLTAN2
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#ifdef HAVE_SHYLU_DDCORE_TPETRA
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
  std::cout << " Not Supported " << std::endl;
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


#ifdef HAVE_SHYLU_DDCORE_ZOLTAN2
template <class Matrix, class Vector>
int PartitionInterface<Matrix, Vector>::partitionZoltan2()
{

  Teuchos::ParameterList subList = pList->sublist("Zoltan2 Input");
  Teuchos::RCP<Matrix> rA(A, false);
  zadapter = new Zoltan2::XpetraCrsMatrixAdapter<Matrix, Vector>(rA);
  zproblem = new Zoltan2::PartitioningProblem<Zoltan2::XpetraCrsMatrixAdapter <Matrix, Vector> >(zadapter, &subList);
  zproblem->solve();
  return 0;
}
#endif


template <class Matrix, class Vector>
int PartitionInterface<Matrix,Vector>::partition()
{
#ifdef HAVE_SHYLU_DDCORE_ZOLTAN2
  return partitionZoltan2();
#else
  return 1;
#endif
}
template <>
int PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::partition()
{
  std::string partitioningPackage = Teuchos::getParameter<std::string>(*pList, "Partitioning Package");
  if(partitioningPackage.compare("Isorropia") == 0)
    {
      return partitionIsorropia();
    }
#ifdef HAVE_SHYLU_DDCORE_ZOLTAN2
  else if (partitioningPackage.compare("Zoltan2") == 0)
    {
      return partitionZoltan2();
    }
#endif
  else
    {
      std::cout << "**Error**: Paritioning package selected is not supported" << std::endl;
    }
  return 1;
}

template <class Matrix, class Vector>
Matrix* PartitionInterface<Matrix,Vector>::reorderMatrix()
{
  Matrix *B;
  std::string partitioningPackage = Teuchos::getParameter<std::string>(*pList, "Partitioning Package");

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
  std::string partitioningPackage = Teuchos::getParameter<std::string>(*pList, "Partitioning Package");
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
  std::string partitioningPackage = Teuchos::getParameter<std::string>(*pList, "Partitioning Package");
#if defined(HAVE_ZOLTAN2_PARMETIS) || defined(HAVE_ZOLTAN2_SCOTCH)
      Teuchos::RCP<Vector> rx(x, false);
      Zoltan2::XpetraMultiVectorAdapter<Vector> tempVecAdapter(rx);
      tempVecAdapter.applyPartitioningSolution(*x, b, zproblem->getSolution());
#endif
  return b;
}
template < >
Epetra_MultiVector* PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>::reorderVector(Epetra_MultiVector* x )
{
  Epetra_MultiVector *b = NULL;
  std::string partitioningPackage = Teuchos::getParameter<std::string>(*pList, "Partitioning Package");
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
