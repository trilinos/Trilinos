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
 * $Name$
 *====================================================================*/

#include <stdio.h>
 
main(argc, argv)
int argc;
char *argv[];
 
{
 
  char first,second,third;
  int count = 0;
  int request;
  int first_int;
  char str[80];
  int  max = 8;

  sprintf(str,"#!/bin/sh");
  count = 0;
  while ( (first_int = getchar()) != EOF ) {
    first = (char)first_int;
    if ( first == str[count]) count++;
    else count = 0;

    if (count == max) break;
  }
  if (count == max) {
     for (count = 0 ; count < max ; count++) putchar(str[count]);
     while ( (first_int = getchar()) != EOF ) {
        putchar(first_int);
     }
  }
}
