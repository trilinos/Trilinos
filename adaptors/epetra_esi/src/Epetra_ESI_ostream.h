#ifndef _Epetra_ESI_ostream_h_
#define _Epetra_ESI_ostream_h_
//
//A convenience: an ostream << operator for epetra_esi::CrsMatrix.
//
ostream& operator<<(ostream& os, epetra_esi::CrsMatrix<double,int>& mat);

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_ostream.cpp"
#endif

#endif

