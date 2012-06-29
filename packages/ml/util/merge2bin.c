#include <stdio.h>
int main(int argc, char *argv[]) {

   int Nrows, Nnzs, current_row, row1, row2, col1, col2, i;
   int which_one;
   int low_row, low_col;
   int flag1, flag2, first_time = 1, j, minusone = -1, itemp;
   char string1[80], string2[80];
   double val1, val2, low_val, dtemp;
   extern int mygetline(FILE *fp, char string[]);
   FILE *fp1, *fp2, *fp3;

fp1 = fopen("upper3","r");
fp2 = fopen("lower3","r");
if ( (fp1 == NULL) || (fp2 == NULL) ) {
  fprintf(stderr,"could not open files\n");
  exit(1);
}
#ifdef binary
if ( (fp3 = fopen("final","wb")) == NULL) {
     fprintf(stderr,"Could not open file\n");
     exit(1);
}
#endif
#ifdef adams
#define both
#endif
#ifdef day
#define both
#endif
#ifndef both
    /* read of the header 'Start of MSR matrix' */


    flag1 = mygetline(fp1,string1);

    /* read of the number of rows and nonzeros  */

    if ((flag1 != 0) && (string1[0] >= '0') && (string1[0] <= '9'))
       flag1 = sscanf(string1,"%d%d",&Nrows, &Nnzs);
    else 
       flag1 = scanf("%d%d",&Nrows, &Nnzs);

    if (flag1 != 2) {
       fprintf(stderr,"Error: Couldn't read number of rows and nonzeros\n");
       exit(1);
    }
#else
   if (argc != 2) {
       fprintf(stderr,"Usage: az_capt2read Nrows\n");
       exit(1);
   }
   if (sscanf(argv[1],"%d",&Nrows) != 1) {
      fprintf(stderr,"Couldn't parse command line\n");
      exit(1);
   }
#ifdef adams
flag1 = scanf("%d%d",&Nrows, &Nnzs);
#endif
#endif
#ifdef binary
    fwrite(&Nrows, sizeof(int), 1, fp3);
#elif
    printf("%d\n",Nrows);
#endif
    flag1 = mygetline(fp1,string1);
    flag2 = mygetline(fp2,string2);
    if (flag1 == 0) flag1 = mygetline(fp1,string1);
    if (flag2 == 0) flag2 = mygetline(fp2,string2);

    current_row = 1;
    while (( flag1 != 0) || (flag2 != 0) ) {
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
       sscanf(string1,"%d%d%lf",&row1,&col1,&val1);
       sscanf(string2,"%d%d%lf",&row2,&col2,&val2);
       if ( (row1 <= row2) || (flag2 == 0) ) {
           low_row = row1;
           low_col = col1;
           low_val = val1;
           which_one = 1;
       }
       else {
           low_row = row2;
           low_col = col2;
           low_val = val2;
           which_one = 2;
       }
       if (flag1 == 0) {
           low_row = row2;
           low_col = col2;
           low_val = val2;
           which_one = 2;
       }
       if (low_row > Nrows) {
          fprintf(stderr,"Error: row (%d) exceeds total number of rows (%d).\n",
		 	  low_row, Nrows);
          exit(1);
       }
       if ((low_col > Nrows) && first_time) {
          fprintf(stderr,"Warning: col (%d) exceeds total number of rows (%d).\n",
		 	  low_col, Nrows);
          first_time = 0;
       }
       if (low_row == current_row+1) 
#ifdef binary
          fwrite(&minusone, sizeof(int), 1, fp3);
#elif
          printf("-1\n");
#endif
       else if (low_row > current_row+1) {
          if (first_time) 
             fprintf(stderr,"Warning: Empty rows (e.g. %d)?\n",current_row+1);
#ifdef binary
          fwrite(&minusone, sizeof(int), 1, fp3);
#elif
          printf("-1\n");
#endif
          for (i = current_row+1; i < low_row; i++) {
#ifdef binary
             itemp = i-1;
             fwrite(&itemp, sizeof(int), 1, fp3);
             dtemp = 0.0;
             fwrite(&dtemp, sizeof(double), 1, fp3);
             fwrite(&minusone, sizeof(int), 1, fp3);
#elif
             printf("%d 0.0\n-1\n",i-1);
#endif
          }
          first_time = 0;
       }
       else if (low_row < current_row) {
          fprintf(stderr,"Error: rows are not in order\n");
          exit(1);
       }
       if ((low_val != 0.0) || (low_col == low_row)) {
#ifdef binary
          itemp = low_col-1;
          fwrite(&itemp, sizeof(int), 1, fp3);
          fwrite(&low_val, sizeof(double), 1, fp3);
#elif
          printf("%d  %20.13e\n",low_col-1,low_val);
#endif
       }
       current_row = low_row;
       if ( which_one == 1) {
          flag1 = mygetline(fp1,string1);
       }
       else flag2 = mygetline(fp2,string2);
   }
#ifdef binary
   fwrite(&minusone, sizeof(int), 1, fp3);
#elif
   printf("-1\n");
#endif
}

int mygetline(FILE *fp, char string[])
{
   int ch, i = 0;

   while ( (ch=getc(fp)) != EOF) {
        string[i++] = ch;
        if (ch == '\n') break;
   }
   string[i] = '\0';
   if (i < 2) return(0);
   else return(1);
}
 
