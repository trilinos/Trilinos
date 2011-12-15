/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include <iostream>

#include "Ifpack2_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond(Teuchos::ParameterList& test_params,
              const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> >& A)
{
  Teuchos::Time timer("precond");

  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
  Teuchos::RCP<Tprec> prec;
  Ifpack2::Factory factory;

  std::string prec_name("not specified");
  Ifpack2::getParameter(test_params, "Ifpack2::Preconditioner", prec_name);
  prec = factory.create(prec_name, A);

  Teuchos::ParameterList tif_params;
  if (test_params.isSublist("Ifpack2")) {
    tif_params = test_params.sublist("Ifpack2");
  }

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "Configuring/Initializing/Computing Ifpack2 preconditioner..."
       << std::endl;
  }

  prec->setParameters(tif_params);
  prec->initialize();
  timer.start();
  prec->compute();
  timer.stop();

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "... Finished Computing Ifpack2 preconditioner (time: "<<timer.totalElapsedTime() << "s)"
       << std::endl;
  }

  // typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  // magnitudeType condest = prec->computeCondEst(Ifpack2::Cheap);
  // if (A->getRowMap()->getComm()->getRank() == 0) {
  //   std::cout << "Condition estimate(cheap) for preconditioner on proc 0: "
  //             << condest << std::endl;
  // }
  return prec;
}

#endif

