#include <stdio.h>
#include <string.h>

main ()
   {
   FILE *fp, *fout;
   char input[200];
   char *a, *b, *c;
   char current[26];
   
   fp = fopen ("./reportfile", "r");
   fout = fopen ("./reportfile.tmp", "w");
   current[0] = 0;

   while (fgets (input, 200, fp) != NULL)
      {
      if (strncmp (input, current, 25) != 0)
         {
         strncpy (current, input, 25);
                  
         a = input + 8;
         input[13] = 0;
         input[22] = 0;
         b = input + 14;
         c = input + 28;

         fprintf (fout, "%s %s %s", a, b, c);
         }
      }
   }
