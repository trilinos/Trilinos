#include <stdio.h>
#include <malloc.h>

/**************************************************************************/
/* Look for any of the patterns given by argv[] in the standard input.    */
/* If any patterns match with the standard input we print a '1',          */
/* otherwise we print a '0'.                                              */
/**************************************************************************/
int main(int argc, char *argv[]){

  int ch, i, j, Nstrings, *Nmatched;
  int start_candidate, str_length;

  Nstrings = argc - 1;
  Nmatched = (int *) malloc(sizeof(int)*(Nstrings+1));
  for (i = 0; i <= Nstrings; i++) Nmatched[i] = 0;

  while ( (ch = getchar()) != EOF) {
    for (i = 1; i <= Nstrings; i++) {
      if ( argv[i][Nmatched[i]] == ch ) {
	Nmatched[i]++;
        /* Have we matched argv[i]? */
        if (Nmatched[i] == strlen(argv[i])) {
	  /* woops it could match and be a substring of a */
	  /* longer word ... and I don't want that!!!!    */
	  ch = getchar();
	  if ( ((ch >= 'a') && (ch <= 'z')) ||
	       ((ch >= '0') && (ch <= '9')) ||
	       ((ch >= 'A') && (ch <= 'Z'))) { 
	    Nmatched[i] = 0; ungetc(ch,stdin);
	  }
	  else {
	    printf("1\n");
	    /*	    printf("matched %s  %d\n",argv[i],i); */
	    return 0;
	  }
        }
      }
      else {
        /* The character 'ch' does not fit the string that we were */
        /* trying to match. However, there is a chance that a      */
        /* substring could match. For example:                     */
        /*    argv[i] = 'ababc'                                    */
        /*    stdin   = 'abababc'                                  */
        /* When the 3rd 'a' is encountered in stdin, we need to    */
        /* recognize that even though 'ababa' does not match,      */
        /* the tail end ('aba') is the start of a new string which */
        /* could potentially match. This is done by checking all   */
        /* the substring possibilities.                            */

	start_candidate = 1;
	while (start_candidate <= Nmatched[i]) {
	  str_length = Nmatched[i] - start_candidate;
	  for (j = 0; j < str_length; j++) {
	    if  (argv[i][start_candidate+j] != argv[i][j]) break;
	  }
	  if ((j==str_length) && (argv[i][str_length] == ch)) {
	    Nmatched[i] = str_length + 1; 
            start_candidate = Nmatched[i]+1;
	  }
          else {
	    start_candidate++;
	    if (start_candidate > Nmatched[i]) Nmatched[i] = 0;
	  }
	}
      }
    }
  }
  printf("0\n");
  return 0;
}
