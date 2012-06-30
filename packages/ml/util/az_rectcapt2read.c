#include <stdio.h>
/*
 *  Convert matlab-type sparse matrix data to AZ_read_msr()
 *  format. This routine is mostly used for matrices (though
 *  it can be used for vectors as well). Each matrix entry
 *  should appear on a line of the form
 *      row_index col_index value
 *
 *  Note: 
 *     1) This routine should work with AZ_capture_matrix.dat and
 *        AZ_capture_rhs.dat files.
 *     2) row and column indices start from 1 (not 0).
 *     3) rows should be sorted. That is all entries for
 *        row i must appear before row i+k (k > 0).
 *     4) The number of rows must be given. This can either
 *        be passed on the command line (i.e. argv[1])
 *        or by preceeding the matrix data with a single
 *        line containing 2 numbers (the first is the number
 *        of matrix rows).
 *     5) String data proceeding the matrix/vector numbers is ignored.
 *
 *
 *  This routine can be used for vectors instead of matrices. In this
 *  case, each line should contain only 1 entry correspond to the vector 
 *  element value.
 *
 */
int main(int argc, char *argv[]) {

   int Nrows, Nnzs, current_row, row, col, i;
   int flag, first_time = 1, j, ch, empty;
   char string[80];
   double val;
   extern int mygetline(char string[]);

   /* We need to determine if there is some kind of header in the  */
   /* file. First strip off any strings.                           */ 

   
   ch = getchar();
   while ( ((ch >= 'a') && (ch <= 'z')) || ((ch >= 'A') && (ch <= 'Z'))) {
      while (ch != '\n') ch = getchar();
      ch = getchar();
   }

   /* Now strip off any values indicating the size and number of nonzeros. */
   i = 0;
   string[i++] = ch;
   while (ch != '\n') {
      ch = getchar();
      string[i++] = (char ) ch;
   }
   string[i] = '\0';

   flag = 0;
   if ( (i = sscanf(string,"%d%d%lf",&Nrows,&j,&val)) != 2) {
             /* No header found (line with just 2 numbers) */
      flag = 1;
      if (argc == 2) {
         if (sscanf(argv[1],"%d",&Nrows) != 1) {
            printf("Couldn't parse command line\n");
            exit(1);
         }
      }
      else {
         printf("No header found nor was the number of rows \n");
         printf("specified on the command line. You must either \n");
         printf("do\n\naz_capt2read Nrows\n\nor the file must have a line ");
         printf("preceeding the\nmatrix data with just two numbers. The ");
         printf("first\nindicating the number of matrix rows.\n");
         exit(1);
      }
   }
   printf("%d\n",Nrows);

   if (i == 1) {
      if (sscanf(string,"%lf",&val) != 1) {
         printf("can not parse first number?\n");
         exit(1);
      }
      printf("%d  %20.13e\n-1\n",0,val);
      for (i = 1; i < Nrows; i++) {
         scanf("%lf",&val);
         printf("%d  %20.13e\n-1\n",i,val);
      }
      exit(1);
   }


   if (flag == 0) flag = mygetline(string);

   current_row = 1;
   empty = 1;
   while ( flag != 0) {
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
	  if (empty) printf("0 0.0\n");
          printf("-1\n"); 
	  /*
          for (i = current_row+1; i < row; i++) printf("%d 0.0\n-1\n",i-1);
	  */
          for (i = current_row+1; i < row; i++) printf("%d 0.0\n-1\n",0);
          first_time = 0;
       }
       else if (row < current_row) {
          fprintf(stderr,"Error: rows are not in order\n");
          exit(1);
       }
       if ((val != 0.0) || (col == row)) {
          printf("%d  %20.13e\n",col-1,val);
          empty = 0;
       }
       current_row = row;
       flag = mygetline(string);
   }
   printf("-1\n");
   for (i = current_row; i < Nrows; i++) 
      printf("0 0.0\n-1\n");
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
 
