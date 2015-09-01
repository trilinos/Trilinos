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

/** \file shylu_partition_interface_decl.hpp
    
    \brief Epetra/Tpetra templated interface for calls to Zoltan(Isorropia)/Zoltans

    \author Joshua Dennis Booth
*/
#ifndef SHYLU_PARTITION_INTERFACE_DECL_HPP
#define SHYLU_PARTITION_INTERFACE_DECL_HPP

#include "ShyLUCore_config.h"

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
#ifdef HAVE_SHYLUCORE_TPETRA
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#endif

//#include <Zoltan2_config.h>
#ifdef HAVE_SHYLUCORE_ZOLTAN2
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#include <Teuchos_XMLParameterListHelpers.hpp>

namespace ShyLU{

  /** \brief PartitionInterface class templated on Epetra/Tpetra Matrix and Vector
   *
   * This class acts as an interface that will allow Shylu to call either Zoltan(Isorropia)/Zoltan2
   * without having to address if the matrix/submatrix is either Epetra or Tpetra form.
   * Currently:  Only limited support for different partitioners has been added
   */
template <class Matrix, class Vector>
class PartitionInterface
{

public:
  ~PartitionInterface();

  /** \brief Main constructor of class
   *
   * This constructor requires a Teuchos ParameterList that provide information on the partitioning method.
   * It assumes that the correct sublist matches the packaged called based on matrix type (Epetra/Tpetra)
   */
  PartitionInterface(Matrix* inA, Teuchos::ParameterList* pList);
  int partition();
  Matrix* reorderMatrix();
  Vector* reorderVector(Vector* x);

private:
  int partitionIsorropia();
  Teuchos::ParameterList* pList;
  Matrix* A;
  ///other handlers needed by zoltan2 and isorropia
#ifdef HAVE_SHYLUCORE_ZOLTAN2
  int partitionZoltan2();
  Zoltan2::XpetraCrsMatrixAdapter<Matrix,Vector> *zadapter;
  Zoltan2::PartitioningProblem<Zoltan2::XpetraCrsMatrixAdapter<Matrix,Vector> > *zproblem;
#endif
  Isorropia::Epetra::Partitioner *ipart;
  Isorropia::Epetra::Redistributor *ird;

};

}// end namespace ShyLU

#endif




