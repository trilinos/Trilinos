
#include <stdio.h>

main (int argc, char **argv)
   {
   char testname[20], method[20], inputname[100], parseline[200] ;
   FILE *input, *output ;
   int hours, minutes, i ;
   float seconds, total_time=-1.0, cuts1=-1, cuts2=-1 ;
   float balance=0.0, wbalance=-1.0 ;
   int init_mem_total, init_mem_max, final_mem_total, final_mem_max ;
   int V, E, P, parts=0 ;
   int count, maxcount=0, temp, temp_vsmall, vsmall=0, redl=-1 ;
   int small_cutt, small_cutl, final ;
   float temp_small_cutt, temp_small_cutl ;
   char local[80], global[80] ;


   sprintf (inputname, "./answers/%s", argv[1]) ;
   sscanf (argv[1], "%[^.].%s", testname, method) ;
   for (i=0 ; i < 80 ; i++)
      if (testname[i] == '-')
          testname[i] = 0 ;

   input = fopen (inputname, "r") ;
   if (input == NULL)
      printf ("unable to open input file %s\n", inputname) ;

   output = fopen ("./reportfile", "a") ;
   if (output == NULL)
      printf ("unable to open output file\n") ;

   while (fgets (parseline, 199, input) != NULL)
      {
      temp = sscanf (parseline, " TIME : %d:%d:%f", &hours, &minutes, &seconds) ;
      if (temp == 3)
         total_time = 3600 * hours + 60 * minutes + seconds ;
      temp = sscanf (parseline, " START %d |V|= %d", &count, &temp_vsmall) ;
      if (temp == 2)
         if (count > maxcount)
            {
            maxcount = count ;
            vsmall = temp_vsmall ;
            }
      temp = sscanf (parseline, " FINAL %d %*s %*s %*s %*s %*s p=%d cutt=%f cutl=%f",
       &final, &parts, &temp_small_cutt, &temp_small_cutl) ;
      if (temp == 4 && final == maxcount)
         {
         small_cutt = (int) temp_small_cutt ;
         small_cutl = (int) temp_small_cutl ;
         }
      temp = sscanf (parseline, " Cuts(total/links) : %f %f", &cuts1, &cuts2) ;
      temp = sscanf (parseline, " Balance : %f", &balance) ;
      temp = sscanf (parseline, " Balance w. : %f", &wbalance) ;
      temp = sscanf (parseline, " info:%*[^|]|V|=%d |E|=%d |P|=%d", &V,&E,&P) ;
      temp = sscanf (parseline, " Initial Memory: %d %d", &init_mem_total, &init_mem_max) ;
      temp = sscanf (parseline, " Final Memory: %d %d", &final_mem_total, &final_mem_max) ;
      temp = sscanf (parseline, " local %[^,], global %[^,], redl %d", local, global, &redl) ;
      }
/*
   if (total_time == -1.0)
      fprintf (output, "%9s %3s    ????\n", testname, method) ;
   else */
        if (wbalance < 0.0)
      fprintf (output, "%9s %2d %3s %3s %3d %3s %4d %4d %3.1f %3.1f %6.3f %5d %5d %5.1f %6.2f\n",
       testname, parts, method, global, redl, local, maxcount, vsmall,
       /*cuts1 == 0 ? 0.0 : (float)small_cutt/(float)cuts1, */ 0.0,
       /*cuts2 == 0 ? 0.0 : (float)small_cutl/(float)cuts2, */ 0.0,
       balance, (int) cuts1, (int) cuts2,
       (float)final_mem_max/(float)init_mem_max, total_time) ;
   else
      fprintf (output, "%9s %2d %3s %3s %3d %3s %4d %4d %3.1f %3.1f %6.3f %5d %5d %5.1f %6.2f \n",
       testname, parts, method, global, redl, local, maxcount, vsmall,
       /* cuts1 == 0 ? 0.0 : (float)small_cutt/(float)cuts1, */ 0.0,
       /* cuts2 == 0 ? 0.0 : (float)small_cutl/(float)cuts2, */ 0.0,
       wbalance, (int) cuts1, (int) cuts2,
       (float)final_mem_max/(float)init_mem_max, total_time) ;


   fclose (input) ;
   fclose (output) ;
   }

