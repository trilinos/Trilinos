#include <stdio.h>
main() {

   extern int mygetline(char string[]);
int  flag, row, col, j;
char string[80];
double val;


   flag = mygetline(string);
/*   printf("%s",string); */

   flag = mygetline(string);
   while (flag != 0) {
/*      printf("%s",string); */
string[0] = ' ';
string[1] = ' ';
string[2] = ' ';
string[3] = ' ';
j = 3; while ( string[j] != ',') j++;
string[j] = ' ';
j = 4; while ( string[j] != ')') j++;
string[j] = ' ';
while ( string[j] != '=') j++;
string[j] = ' ';
while ( string[j] != ';') j++;
string[j] = ' ';
       sscanf(string,"%d%d%lf",&row,&col,&val);
/*
       printf("Kaa(%6d,%6d)=%20.13e;\n",row,col,val);
*/
       printf("%6d %6d %20.13e\n",row,col,val);
if (row != col)
       printf("%6d %6d %20.13e\n",col,row,val);
/*
       printf("Kaa(%6d,%6d)=%20.13e;\n",col,row,val);
*/
       flag = mygetline(string);
   }

}




mygetline(char string[])
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

