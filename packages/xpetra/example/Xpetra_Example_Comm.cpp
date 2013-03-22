// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_Comm.hpp>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include <Xpetra_EpetraComm.hpp>

#include "Tpetra_DefaultPlatform.hpp"

// This driver simply tests Teuchos2Epetra_Comm

int main(int argc, char** argv)
{
  typedef int                                                       Ordinal;
  typedef double                                                    Scalar;

  typedef Tpetra::MpiPlatform<Kokkos::SerialNode>                   MpiPlatform;
  typedef Tpetra::SerialPlatform<Kokkos::SerialNode>                SerialPlatform;

  typedef MpiPlatform::NodeType                                     MpiNodeType;
  typedef SerialPlatform::NodeType                                  SerialNodeType;

  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  ParameterList pl;

  Ordinal numThreads=1;
  pl.set("Num Threads",numThreads);
  RCP<MpiNodeType> mpiNode = rcp(new MpiNodeType(pl));
  RCP<SerialNodeType> serialNode = rcp(new SerialNodeType(pl));

  MpiPlatform    myMpiPlat(mpiNode);
  SerialPlatform mySerialPlat(serialNode);

  {
    RCP<const Comm<int> > teuchosComm = mySerialPlat.getComm();
    RCP<const Epetra_Comm> epetraComm = Teuchos2Epetra_Comm(teuchosComm);

    assert(epetraComm != Teuchos::null);
  }

  {
    RCP<const Comm<int> > teuchosComm = myMpiPlat.getComm();
    RCP<const Epetra_Comm> epetraComm = Teuchos2Epetra_Comm(teuchosComm);

    assert(epetraComm != Teuchos::null);
  }

  return(0);
} //main
