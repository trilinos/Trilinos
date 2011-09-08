#include "fortranc.h"
#include <stdio.h>
#include <string.h>

#if defined ADDC_
void wildcmp_(char *fwild, char *fstring, FTNINT *ret, FTNINT wlen, FTNINT slen)
#else
void wildcmp(char *fwild, char *fstring, FTNINT *ret, FTNINT wlen, FTNINT slen)
#endif
{
  /* Written by Jack Handy - jakkhandy@hotmail.com */

  char *cp = 0, *mp = 0;
  char cwild[2048], cstring[2048];
  char *wild = cwild;
  char *string = cstring;

  strncpy(wild, fwild, wlen); wild[wlen] = '\0';
  strncpy(string, fstring, slen); string[slen] = '\0';
  
  fprintf(stderr, "%s %s %d %d %d %d\n", wild, string, wlen, slen, strlen(wild), strlen(string));
  while ((*string) && (*wild != '*')) {
    if ((*wild != *string) && (*wild != '?')) {
      fprintf(stderr, "false\n");
      *ret= 0;
      return;
    }
    wild++;
    string++;
  }

  while (*string) {
    if (*wild == '*') {
      if (!*++wild) {
	fprintf(stderr, "true\n");
	*ret= 1;
        return;
      }
      mp = wild;
      cp = string+1;
    } else if ((*wild == *string) || (*wild == '?')) {
      wild++;
      string++;
    } else {
      wild = mp;
      string = cp++;
    }
  }

  while (*wild == '*') {
    wild++;
  }
  *ret = !*wild;
  if (*ret) 
    fprintf(stderr, "true\n");
  else
    fprintf(stderr, "false\n");
    
  return;
}

#if 0
#include <stdio.h>
#include <assert.h>
void test(char *pat, char *stng, int expect)
{
  if (wildcmp(pat, stng) != expect) {
    fprintf(stderr, "Test failed\n");
  }
    
}
int main()
{
  test( "", "", 1 );
  test( "*", "", 1 );
  test( "*", "A", 1 );
  test( "", "A", 0 );
  test( "A*", "", 0 );
  test( "A*", "AAB", 1 );
  test( "A*", "BAA", 0 );
  test( "A*", "A", 1 );
  test( "A*B", "", 0 );
  test( "A*B", "AAB", 1 );
  test( "A*B", "AB", 1 );
  test( "A*B", "AABA", 0 );
  test( "A*B", "ABAB", 1 );
  test( "A*B", "ABBBB", 1 );
  test( "A*B*C", "", 0 );
  test( "A*B*C", "ABC", 1 );
  test( "A*B*C", "ABCC", 1 );
  test( "A*B*C", "ABBBC", 1 );
  test( "A*B*C", "ABBBBCCCC", 1 );
  test( "A*B*C", "ABCBBBCBCCCBCBCCCC", 1 );
  test( "A*B*", "AB", 1 );
  test( "A*B*", "AABA", 1 );
  test( "A*B*", "ABAB", 1 );
  test( "A*B*", "ABBBB", 1 );
  test( "A*B*C*", "", 0 );
  test( "A*B*C*", "ABC", 1 );
  test( "A*B*C*", "ABCC", 1 );
  test( "A*B*C*", "ABBBC", 1 );
  test( "A*B*C*", "ABBBBCCCC", 1 );
  test( "A*B*C*", "ABCBBBCBCCCBCBCCCC", 1 );
  test( "A?", "AAB", 0 );
  test( "A?B", "AAB", 1 );
  test( "A?*", "A", 0 );
  test( "A?*", "ABBCC", 1 );
  test( "A?*", "BAA", 0 );
  test("a*bc",   "abbc", 1);
  test("a*c",    "abbc", 1);
  test("a*b",    "a", 0);
  test("a*?b",   "axb", 1);
  test("a**b",   "axb", 1);
  test("bl?h.*", "blah.jpg", 1);
  test("bl?h.*", "blaah.kk", 0);
  test("bl?h.*", "blah.", 1);
}
#endif
