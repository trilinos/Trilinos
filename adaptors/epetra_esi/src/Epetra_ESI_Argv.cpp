#ifndef _Epetra_ESI_Argv_cpp_
#define _Epetra_ESI_Argv_cpp_

//include Epetra_Object.h, which is where iostream, etc., get included.
#include "Epetra_Object.h"

#ifndef EPETRA_ESI_INCLUDE_IMPLEMENTATION
//we're going to turn on the flag to include implementations, but if it isn't
//already on, we'll need to turn it back off afterwards... is that clear?
#define Epet_ESI_Arg_UNDEF
#endif

#define EPETRA_ESI_INCLUDE_IMPLEMENTATION

#include "Epetra_ESI_Argv.h"

#ifdef Epet_ESI_Arg_UNDEF
#undef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#undef Epet_ESI_Arg_UNDEF
#endif

//------------------------------------------------------------------------------
epetra_esi::Argv::Argv()
 : argv_(0, 1)
{
}

#endif

