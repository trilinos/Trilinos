//@HEADER
// ***********************************************************************
// 
//                       MATLAB Engine Package
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
//@HEADER

#ifndef MATLAB_ENGINE_H
#define MATLAB_ENGINE_H
#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"

// the following deal with matlab provided headers:
#include "engine.h"
#include "mex.h"
#undef printf
//! MatlabEngine: 

class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_CrsGraph;
class Epetra_SerialDenseMatrix;
class Epetra_BlockMap;
class Epetra_IntSerialDenseMatrix;

/*! The EpetraExt_MatlabEngine class allows Epetra data objects to be exported to Matlab. 

<b>Error Codes</b>
<ul>
  <li> -1 engPutVariable returned a nonzero result
  <li> -2 internal gather of epetra object and copy to matlab object returned a nonzero result
  <li> -3 engEvalString returned a nonzero result
  <li> -4 engOutputBuffer returned a nonzero result
</ul>


*/

//using namespace EpetraExt;
namespace EpetraExt {
//=========================================================================
class MatlabEngine {

  public:

  //@{ \name Constructors/destructors.
  //! Basic MatlabEngine constuctor.
  /*! Creates a MATLAB_Engine object and fills with default values.  

    \param Comm In
           An Epetra Communicator 

    \return  MATLAB_Engine object

  */
  MatlabEngine(const Epetra_Comm& Comm);
  ~MatlabEngine();
  //! MatlabEngine copy constructor.  There is no copy...
  
  // MatlabEngine(const MATLAB_Engine& Source);
  //@}
  
  //@{ \name LotsHere methods

  //! EvalString method
  int EvalString (char* command);
  
  /*! err codes
   -1 engOutputBuffer matlab call returned an error
   -2 engEvalString matlab call returned an error */
  int EvalString (char* command, char* outputBuffer, int outputBufferSize);

  //@}

  int PutMultiVector(const Epetra_MultiVector& multiVector, const char* variableName);
  int PutRowMatrix(const Epetra_RowMatrix& A, const char* variableName, bool transA);
  int PutCrsGraph(const Epetra_CrsGraph& A, const char* variableName, bool transA);
  int PutSerialDenseMatrix(const Epetra_SerialDenseMatrix& A, const char* variableName);
  int PutIntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& A, const char* variableName);
  int PutBlockMap(const Epetra_BlockMap& blockMap, const char* variableName);
  
 private:

    Engine* Engine_ ;
    int MyPID_ ;

    const Epetra_Comm& Comm_ ;

};
} // namespace EpetraExt

#endif /* MATLAB_ENGINE_H */
