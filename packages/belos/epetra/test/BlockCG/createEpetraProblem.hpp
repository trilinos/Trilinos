// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifdef EPETRA_MPI
#include "mpi.h"
#endif

class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_MultiVector;

#include "Teuchos_RCP.hpp"
using Teuchos::RCP;
using Teuchos::rcp;

namespace Belos {

class MPIFinalize {
public:
  ~MPIFinalize() {
#ifdef EPETRA_MPI 
    MPI_Finalize();
#endif
  }
};

int createEpetraProblem(
			std::string                      &filename
			,RCP<Epetra_Map>         *rowMap
			,RCP<Epetra_CrsMatrix>   *A
			,RCP<Epetra_MultiVector> *B
			,RCP<Epetra_MultiVector> *X
			,int                             *MyPID
			);

} // namespace Belos
