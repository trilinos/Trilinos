//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "EpetraExt_readEpetraLinearSystem.h"
#include "Trilinos_Util.h"

void EpetraExt::readEpetraLinearSystem(
  const std::string                               &fileName
  ,const Epetra_Comm                              &comm
  ,Teuchos::RefCountPtr<Epetra_CrsMatrix>         *A
  ,Teuchos::RefCountPtr<Epetra_Map>               *map
  ,Teuchos::RefCountPtr<Epetra_Vector>            *x
  ,Teuchos::RefCountPtr<Epetra_Vector>            *b
  ,Teuchos::RefCountPtr<Epetra_Vector>            *xExact
  )
{

  Epetra_Map       *readMap;
  Epetra_CrsMatrix *readA; 
  Epetra_Vector    *readx; 
  Epetra_Vector    *readb;
  Epetra_Vector    *readxexact;

  const std::string::size_type ext_dot = fileName.rfind(".");
  TEST_FOR_EXCEPT( ext_dot == std::string::npos );
  std::string ext = fileName.substr(ext_dot+1);
  //std::cout << "\nfileName = " << fileName << "\next = " << ext << std::endl;

  char *hacked_file_str = const_cast<char*>(fileName.c_str());
  
  if ( ext == "triU" ) { 
    const bool NonContiguousMap = true; 
    TEST_FOR_EXCEPT(
      0!=Trilinos_Util_ReadTriples2Epetra(
        hacked_file_str, false, comm, readMap, readA, readx, 
        readb, readxexact, NonContiguousMap
        )
      );
  }
  else if ( ext == "triS" ) { 
    const bool NonContiguousMap = true; 
    TEST_FOR_EXCEPT(
      0!=Trilinos_Util_ReadTriples2Epetra(
        hacked_file_str, true, comm, readMap, readA, readx, 
        readb, readxexact, NonContiguousMap
        )
      );
  }
  else if( ext == "mtx" ) { 
    TEST_FOR_EXCEPT(
      0!=Trilinos_Util_ReadMatrixMarket2Epetra(
        hacked_file_str, comm, readMap, 
        readA, readx, readb, readxexact
        )
      );
  }
  else if ( ext == "hb" ) {
    Trilinos_Util_ReadHb2Epetra(
      hacked_file_str, comm, readMap, readA, readx, 
      readb, readxexact
      ); // No error return???
  }
  else {
    TEST_FOR_EXCEPTION(
      true, std::logic_error
      ,"Error, the file = \'"<<hacked_file_str<<"\' has the extension "
      "\'*."<<ext<<"\' is not \'*.triU\', \'*.triS\', \'*.mtx\', or \'*.hb\'!"
      );
  }

  Teuchos::RefCountPtr<Epetra_CrsMatrix>    loc_A         = Teuchos::rcp(readA);
  Teuchos::RefCountPtr<Epetra_Map>          loc_map       = Teuchos::rcp(readMap);
  Teuchos::RefCountPtr<Epetra_Vector>       loc_x         = Teuchos::rcp(readx);
  Teuchos::RefCountPtr<Epetra_Vector>       loc_b         = Teuchos::rcp(readb);
  Teuchos::RefCountPtr<Epetra_Vector>       loc_xExact    = Teuchos::rcp(readxexact);

  if(A)       *A       = loc_A;
  if(map)     *map     = loc_map;
  if(x)       *x       = loc_x;
  if(b)       *b       = loc_b;
  if(xExact)  *xExact  = loc_xExact;
  
}
