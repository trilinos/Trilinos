/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_file_id = "$Id$";
#endif

#include "spice.h"
#include <stdio.h>
#include "spmatrix.h"
#include "smpdefs.h"
#include "spdefs.h"

void
SMPlinkRows( Matrix )
SMPmatrix *Matrix;
{
    void spcLinkRowsandCreateInternalVectors(char *);

    spcLinkRowsandCreateInternalVectors( (char *)Matrix );
}

void spcLinkRowsandCreateInternalVectors( eMatrix )
char *eMatrix;
{
    MatrixPtr  Matrix = (MatrixPtr)eMatrix;
    
    if (NOT Matrix->RowsLinked)
        spcLinkRows( Matrix );
    if (NOT Matrix->InternalVectorsAllocated)
        spcCreateInternalVectors( Matrix );
}
