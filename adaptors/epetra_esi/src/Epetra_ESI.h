#ifndef _Epetra_ESI_h_
#define _Epetra_ESI_h_

#include "Epetra_Object.h"

/* 
This is a package-header, it includes all of the separate Epetra_ESI
headers. This is just for convenience. Each of the Epetra_ESI_* headers
can be included separately, they each in turn include the headers they
depend on, and they all have include-guards to prevent them from being
included more than once in a compilation unit.
*/

#include "Epetra_ESI_platforms.h"

#include "Epetra_ESI_CHK_ERR.h"

#include "Epetra_Array.h"
#include "Epetra_ESI_utils.h"

#include "Epetra_ESI_Argv.h"
#include "Epetra_ESI_Object.h"
#include "Epetra_ESI_IndexSpace.h"
#include "Epetra_ESI_Vector.h"
#include "Epetra_ESI_CrsMatrix.h"
#include "Epetra_ESI_Operator.h"

#include "Aztec_ESI_Solver.h"

#endif

