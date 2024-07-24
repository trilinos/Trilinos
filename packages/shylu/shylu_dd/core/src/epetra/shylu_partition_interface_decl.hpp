// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_partition_interface_decl.hpp
    
    \brief Epetra/Tpetra templated interface for calls to Zoltan(Isorropia)/Zoltans

    \author Joshua Dennis Booth
*/
#ifndef SHYLU_PARTITION_INTERFACE_DECL_HPP
#define SHYLU_PARTITION_INTERFACE_DECL_HPP

#include "ShyLU_DDCore_config.h"

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
#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#endif

//#include <Zoltan2_config.h>
#if defined(HAVE_SHYLU_DDCORE_ZOLTAN2CORE)
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
#if defined(HAVE_SHYLU_DDCORE_ZOLTAN2CORE) && (!defined(HAVE_SHYLU_DDCORE_TPETRA) || defined(HAVE_TPETRA_INST_INT_INT))
  int partitionZoltan2();
  Zoltan2::XpetraCrsMatrixAdapter<Matrix,Vector> *zadapter;
  Zoltan2::PartitioningProblem<Zoltan2::XpetraCrsMatrixAdapter<Matrix,Vector> > *zproblem;
#endif
  Isorropia::Epetra::Partitioner *ipart;
  Isorropia::Epetra::Redistributor *ird;

};

}// end namespace ShyLU

#endif





#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

