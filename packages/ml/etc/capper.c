#include <stdio.h>
int main(int argc, char *argv[]){

  int ch;

  while ( (ch = getchar()) != EOF) {
    if ( (ch >= 'a') && (ch <= 'z')) {
      ch = ch - 'a' + 'A';
    }
    putchar(ch);
  }
}
