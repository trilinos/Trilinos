#include <stdio.h>
/*
Little function that takes a vector and
spits out the associated diagonal matrix.
*/
main() {

  int i;
  double val;

  while ( scanf("%d%lf",&i,&val) != EOF ) {
    printf("%d %d %e\n",i,i,val);
  }
}
