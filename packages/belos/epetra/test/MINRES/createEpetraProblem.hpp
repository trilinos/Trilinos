//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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

#include "Epetra_config.h"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"

class Epetra_Comm;
class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_MultiVector;

namespace Belos {

  namespace Test {
    class MPISession {
    public:
      MPISession (Teuchos::Ptr<int> argc, Teuchos::Ptr<char**> argv);
      ~MPISession ();
      Teuchos::RCP<const Epetra_Comm> getComm ();

    private:
      // Lazily initialized Epetra communicator wrapper.
      Teuchos::RCP<Epetra_Comm> comm_;
    };
  } // namespace Test


void
createEpetraProblem (const Teuchos::RCP<const Epetra_Comm>& epetraComm,
		     const std::string& filename,
		     Teuchos::RCP<Epetra_Map>& rowMap,
		     Teuchos::RCP<Epetra_CrsMatrix>& A,
		     Teuchos::RCP<Epetra_MultiVector>& B,
		     Teuchos::RCP<Epetra_MultiVector>& X,
		     int &numRHS); // in/out

} // namespace Belos
