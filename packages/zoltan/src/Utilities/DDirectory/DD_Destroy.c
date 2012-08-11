

#include <stdio.h>
#include <stdlib.h>

#include "DD.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif



/**********************  Zoltan_DD_Destroy()  *********************/

void Zoltan_DD_Destroy (
 Zoltan_DD_Directory **dd)     /* contains directory state information */
   {
   int i;
   DD_Node *ptr;
   DD_Node *next;

   char *yo = "ZOLTAN_DD_Destroy";

   /* input sanity check */
   if (dd == NULL || *dd == NULL) {
      ZOLTAN_PRINT_ERROR (0, yo, "Input argument dd is NULL");
      return;
   }
   if ((*dd)->debug_level > 4)
      ZOLTAN_TRACE_IN ((*dd)->my_proc, yo, NULL);

   /* for each linked list head, walk its list freeing memory */
   for (i = 0; i < (*dd)->table_length; i++)
      for (ptr = (*dd)->table[i]; ptr != NULL; ptr = next)  {
         next = ptr->next;            /* save before deletion         */
         ZOLTAN_FREE (&ptr);          /* destroy node                 */
      }

   /* execute user registered cleanup function, if needed */
   if ((*dd)->cleanup != NULL)
       (*dd)->cleanup((*dd)->hashdata);

   MPI_Comm_free (&((*dd)->comm));    /* free MPI Comm, ignore errors */

   if ((*dd)->debug_level > 4)
      ZOLTAN_TRACE_OUT ((*dd)->my_proc, yo, NULL);

   ZOLTAN_FREE (dd);                  /* free directory structure     */
   return;
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
