
#include <stdio.h>
#include <string.h>
#include <math.h>

main (int argc, char **argv)
   {
   char testname[20], parseline[200] ;
   FILE *input, *output ;
   int cuts[25], cut, sum=0, count, i, n=0, k=0 ;
   double mean, var=0.0;
   char current[25]="ibm01";
   int baseline[]={7750, 12751, 31607, 13444, 47633, 25065, 20030, 33830, 20509,
    45227, 26861, 62568, 36791, 56294, 93675, 77278, 81782, 63324};
   int sdev_base[]= {67, 376, 453, 85, 403, 407, 154, 251, 477, 481, 500,
    858, 559, 676, 962, 1087, 999, 1030};
   int ave_base[] ={243, 272, 781, 440, 1718, 376, 762, 1157, 523, 778, 701,
    2006, 884, 1636, 1809, 1723, 2397, 1539};
   double ave_ss = 0.0;
   int best_sum=0, best_sdev=0;
   double secs=0.0, total_secs = 0.0;
   int parts;

   input = fopen ("./reportfile.sorted", "r") ;
   if (input == NULL)
      printf ("unable to open reportfile.sorted file");

   output = fopen ("./reportfile.eval", "w") ;
   if (output == NULL)
      printf ("unable to open output file\n") ;

   while (fgets (parseline, 199, input) != NULL)
      {
      count = sscanf (parseline, "%s %*s p=%d, %*s %d",
       &testname, &parts, &cut);
      if (count != 3)
         {
         fprintf (stderr, "Error reading line\n");
         continue;
         }

      if (strcmp(testname, current) != 0)
         {
         sum = 0;
         var = 0.0;
         for (i = 0; i < n; i++)
            sum += cuts[i];
         mean = sum / (double)n;
         ave_ss += ( ((mean-ave_base[k]) * (mean-ave_base[k])) / (ave_base[k] * ave_base[k]));
         for (i = 0; i < n; i++)
            var += (((double) cuts[i] - mean) * ((double) cuts[i] - mean));

         fprintf (output, "       Sum %d, Mean %.1f, StdDev %.1f\n",
          sum, mean, sqrt(var / (double) (n-1)));
         if (parts == 2)
            {
            fprintf (output, "       Sum is %s, SDEV is %s\n\n",
             sum < baseline[k] ? "+" : "-",
            (int) sqrt (var / (double) (n-1)) < sdev_base[k] ? "+" : "-");
            if (sum < baseline[k])
               best_sum++;
            if ((int) sqrt (var / (double) (n-1)) < sdev_base[k])
               best_sdev++;
            }
         else
            fprintf (output, "\n");

         n = 0;
         k++;
         strcpy (current, testname);
         }

      cuts[n++] = cut;
      fprintf (output, "%s", parseline);
      }
   sum = 0;
   var = 0.0;
   for (i = 0; i < n; i++)
      sum += cuts[i];
   mean = sum / (double)n;
   ave_ss += (((mean-ave_base[k]) * (mean-ave_base[k])) / (ave_base[k] * ave_base[k]));
   for (i = 0; i < n; i++)
      var += (((double) cuts[i] - mean) * ((double) cuts[i] - mean));

   fprintf (output, "       Sum %d, Mean %.1f, StdDev %.1f\n", sum, mean,
    sqrt(var / (double) (n-1)));

    if (parts == 2)
       {
       fprintf (output, "       Sum is %s, SDEV is %s\n\n",
        sum < baseline[k] ? "+" : "-",
        (int) sqrt (var / (double) (n-1)) < sdev_base[k] ? "+" : "-");
       if (sum < baseline[k])
          best_sum++;
       if ((int) sqrt (var / (double) (n-1)) < sdev_base[k])
          best_sdev++;

       fprintf (output, "FINAL:  best sum %d, best sdev %d\nrms of per dev from target %.3f\n",
        best_sum, best_sdev, sqrt(ave_ss/18.0));
       }
    else
       fprintf (output, "\n\n");

   fclose (input) ;
   fclose (output) ;
   }

