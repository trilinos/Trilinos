//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
