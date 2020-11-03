/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int count_digits(int value)
{
  char tmp[32];
  sprintf(tmp, "%d", value);
  return strlen(tmp);
}

/* Output a zero-padded sequence of digits from [0..limit).  For
   example zpad 10 will output 00, 01, 02, 03, ..., 09 each on a
   single line.
 */
int main(int argc, char **argv)
{
  int start = 0;
  if (argc < 2) {
    fprintf(stderr, "Usage: %s limit [start=0]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (argc == 3) {
    start = strtol(argv[2], NULL, 10);
  }

  {
    int  i;
    char format[] = "%.0Xd\n";
    char digits[2];
    int  limit = strtol(argv[1], NULL, 10);

    /* Count number of digits needed to represent 'limit' */
    int width = count_digits(limit);
    sprintf(digits, "%d", width);

    /* Create an output format that will zero-pad to that width */
    format[3] = digits[0];
    for (i = start; i < limit; i++) {
      printf(format, i);
    }
  }
  return EXIT_SUCCESS;
}
