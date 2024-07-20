// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __test_support_hpp__
#define __test_support_hpp__

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Epetra_LocalMap.h"

#include "Stokhos_Epetra.hpp"
  
Teuchos::RCP<Teuchos::ParameterList> buildAppParams(int num_KL,bool full_expansion);

Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > buildBasis(int num_KL,int porder);

Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> >
buildBasis(int num_KL,const std::vector<int> & order);

Teuchos::RCP<Stokhos::ParallelData> buildParallelData(bool full_expansion,int num_KL,
                                                      const Teuchos::RCP<const Epetra_Comm> & globalComm,
                                                      const Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > & basis);

Teuchos::RCP<Stokhos::ParallelData> buildParallelData(bool full_expansion,int num_KL,
                                                      const Teuchos::RCP<const Epetra_Comm> & globalComm,
                                                      const Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > & basis);

#endif
