/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <ctype.h>
#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */

#include "zz_util_const.h"
#include "zoltan_mem.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Remove leading & trailing white space and convert to upper case. */

int Zoltan_Clean_String(
char *string1,			/* original string */
char **pstring2) 		/* cleaned string to return */
{

    char     *string2;		/* cleaned up string */
    int       start, end;	/* indices bounding true string */
    int       length1;		/* length of string 1 */
    int       i;		/* loop counter */

    length1 = strlen(string1);
    start = 0;
    end = length1;
    while (start < length1 && isspace((int)(string1[start])))
	start++;
    while (end > start && isspace((int)(string1[end])))
	end--;

    string2 = (char *) ZOLTAN_MALLOC((end - start + 1) * sizeof(char));
    *pstring2 = string2;

    if (string2 == NULL)
	return (ZOLTAN_MEMERR);

    for (i = start; i < end; i++) {
	*string2++ = toupper(string1[i]);
    }
    *string2 = '\0';

    return (ZOLTAN_OK);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
