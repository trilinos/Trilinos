/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 */

#include <ctype.h>

#if defined(ADDC_)
void exupcs_(char *string, long int StrLength)
#else
void exupcs(char *string, long int StrLength)
#endif
{
  int i;

  for (i = 0; i < StrLength; i++, string++) {
    if (!isprint(*string))
      *string = ' ';
    else if (isalpha(*string) && islower(*string))
      *string = toupper(*string);
  }
  return;
}
