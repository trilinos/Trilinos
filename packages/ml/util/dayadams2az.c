#include <stdio.h>
int main(int argc, char *argv[]) {

   int Nrows, Nnzs, current_row, row, col, i;
   int flag, first_time = 1, j;
   char string[80];
   double val;
   extern int mygetline(char string[]);

#ifdef adams
#define both
#endif
#ifdef day
#define both
#endif
#ifndef both
    /* read of the header 'Start of MSR matrix' */


    flag = mygetline(string);

    /* read of the number of rows and nonzeros  */

    if ((flag != 0) && (string[0] >= '0') && (string[0] <= '9'))
       flag = sscanf(string,"%d%d",&Nrows, &Nnzs);
    else 
       flag = scanf("%d%d",&Nrows, &Nnzs);

    if (flag != 2) {
       printf("Error: Couldn't read number of rows and nonzeros\n");
       exit(1);
    }
#else
   if (argc != 2) {
       printf("Usage: az_capt2read Nrows\n");
       exit(1);
   }
   if (sscanf(argv[1],"%d",&Nrows) != 1) {
      printf("Couldn't parse command line\n");
      exit(1);
   }
#ifdef adams
#ifdef old
flag = scanf("%d%d",&Nrows, &Nnzs);
#endif
#endif
#endif
    printf("%d\n",Nrows);
    flag = mygetline(string);
    if (flag == 0) flag = mygetline(string);

    current_row = 1;
    while ( flag != 0) {
#ifdef day
/*
string[0] = ' ';
string[1] = ' ';
string[2] = ' ';
string[3] = ' ';
j = 4; while ( string[j] != ',') j++;
string[j] = ' ';
j = 4; while ( string[j] != ')') j++;
string[j] = ' ';
while ( string[j] != '=') j++;
string[j] = ' ';
while ( string[j] != ';') j++;
string[j] = ' ';
*/
#endif
       sscanf(string,"%d%d%lf",&row,&col,&val);
       if (row > Nrows) {
          fprintf(stderr,"Error: row (%d) exceeds total number of rows (%d).\n",
		 	  row, Nrows);
          exit(1);
       }
       if ((col > Nrows) && first_time) {
          fprintf(stderr,"Warning: col (%d) exceeds total number of rows (%d).\n",
		 	  col, Nrows);
          first_time = 0;
       }
       if (row == current_row+1) printf("-1\n");
       else if (row > current_row+1) {
          if (first_time) 
             fprintf(stderr,"Warning: Empty rows (e.g. %d)?\n",current_row+1);
          printf("-1\n"); 
          for (i = current_row+1; i < row; i++) printf("%d 0.0\n-1\n",i-1);
          first_time = 0;
       }
       else if (row < current_row) {
          fprintf(stderr,"Error: rows are not in order\n");
          exit(1);
       }
       if ((val != 0.0) || (col == row)) printf("%d  %20.13e\n",col-1,val);
       current_row = row;
       flag = mygetline(string);
   }
   printf("-1\n");
}

int mygetline(char string[])
{
   int ch, i = 0;

   while ( (ch=getchar()) != EOF) {
        string[i++] = ch;
        if (ch == '\n') break;
   }
   string[i] = '\0';
   if (i < 2) return(0);
   else return(1);
}
 
