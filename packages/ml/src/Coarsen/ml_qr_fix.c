#include "ml_include.h"
#include "ml_qr_fix.h"

int ML_qr_fix_Create(ML_qr_fix** ptr)
{
  *ptr = (ML_qr_fix*) ML_allocate(sizeof(ML_qr_fix));
  (*ptr)->level = 0;

  return(0);
}

int ML_qr_fix_Destroy(ML_qr_fix** ptr)
{
  ML_qr_fix* fix = *ptr;
  if (fix == NULL)
    return(0);

  /* loop over all fields and release them */

  /* free the structure itself */

  ML_free(*ptr);
  *ptr = NULL;

  return(0);
}

int ML_qr_fix_Print(ML_qr_fix* ptr)
{
  if (ptr == NULL)
    return(-1);

  printf("level = %d\n", ptr->level);

  return(0);
}
