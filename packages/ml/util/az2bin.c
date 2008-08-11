#include <stdio.h>
main(argc, argv)

char *argv[];
int   argc;
{
  FILE  *fp,*fp2;
  int i,total,itemp;
  double dtemp;

  if (argc != 3) {
     fprintf(stderr,"Usage: az2bin input_file output_file\n");
     exit(1);
  }
  if ( (fp = fopen(argv[1],"r")) == NULL) {
     fprintf(stderr,"Could not open file: %s\n",argv[1]);
     exit(1);
  }
  if ( (fp2 = fopen(argv[2],"wb")) == NULL) {
     fprintf(stderr,"Could not open file: %s\n",argv[2]);
     exit(1);
  }

  fscanf(fp,"%d",&total);
  fwrite(&total, sizeof(int), 1, fp2);

  for (i = 0 ; i < total ; i++ ) {
     fscanf(fp,"%d", &itemp);
     fwrite(&itemp, sizeof(int), 1, fp2);
     while (itemp != -1) {
        fscanf(fp,"%lf",&dtemp);
        fwrite(&dtemp, sizeof(double), 1, fp2);
        fscanf(fp,"%d", &itemp);
        fwrite(&itemp, sizeof(int), 1, fp2);
     }
  }
  fclose(fp); fclose(fp2);
}

