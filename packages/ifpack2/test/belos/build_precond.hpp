/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include <iostream>

#include "Ifpack2_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond (Teuchos::ParameterList& test_params,
               const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A)
{
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using std::cout;
  using std::endl;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  Teuchos::Time timer_init("init");
  Teuchos::Time timer("precond");
  Teuchos::Time timer2("precond_reuse");
  const int myRank = A->getRowMap ()->getComm ()->getRank ();

  RCP<FancyOStream> out = getFancyOStream (rcpFromRef (cout));

  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
  Teuchos::RCP<Tprec> prec;
  Ifpack2::Factory factory;

  std::string prec_name("not specified");
  Ifpack2::getParameter(test_params, "Ifpack2::Preconditioner", prec_name);
  prec = factory.create<row_matrix_type> (prec_name, A);

  Teuchos::ParameterList tif_params;
  if (test_params.isSublist("Ifpack2")) {
    tif_params = test_params.sublist("Ifpack2");
  }
  bool reuse_pattern = false;
  if (test_params.isParameter("Reuse Pattern"))
    {
      reuse_pattern = test_params.get<bool>("Reuse Pattern");
    }

  if (myRank == 0) {
    *out << "Configuring, initializing, and computing Ifpack2 preconditioner" << endl;
  }
  {
    OSTab tab (*out);
    prec->setParameters (tif_params);
    {
      Teuchos::TimeMonitor timeMon (timer_init);
      prec->initialize ();
    }
     if (myRank == 0) {
      *out << "Time Init: " << timer_init.totalElapsedTime() << endl;
     }
    {
      Teuchos::TimeMonitor timeMon (timer);
      prec->compute ();
    }
    if (myRank == 0) {
      *out << "Finished computing Ifpack2 preconditioner" << endl;
      OSTab tab2 (*out);
      *out << "Time (s): " << timer.totalElapsedTime () << endl;
    }
    if (reuse_pattern == true)
      {
	{
	  Teuchos::TimeMonitor timeMon (timer2);
	  prec->compute ();
	}
	if (myRank == 0) {
	  *out << "Finished recomputing Ifpack2 preconditioner" 
	       << endl;
	  OSTab tab2 (*out);
	  *out << "Time (s): " << timer2.totalElapsedTime () 
	       << endl;
	}
      }
  }

  return prec;
}

#endif

