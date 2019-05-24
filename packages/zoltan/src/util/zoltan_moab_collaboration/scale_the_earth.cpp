// Scales the coordinates of a mesh to a sphere with radius one.
// Allows atmosphere and ocean to be partitioned together with rcb two-weight

#include <stdio.h>
#include <math.h>

typedef double point[3];

#define ABS(x) ((x) >= 0. ? (x) : -(x))

int getNumCoords(FILE *fp)
{
  int nc;
  fscanf(fp, "%d", &nc);
  return nc;
}

point *scaled_coordinates(int nc, FILE *fp_coords)
{
  double max[3] = {0., 0., 0.};
  point *xyz = new point[nc];

  // read coords and find max magnitude
  for (int i = 0; i < nc; i++) {
    fscanf(fp_coords, "%lf %lf %lf\n", &xyz[i][0], &xyz[i][1], &xyz[i][2]);
    if (ABS(xyz[i][0]) > max[0]) max[0] = ABS(xyz[i][0]);
    if (ABS(xyz[i][1]) > max[1]) max[1] = ABS(xyz[i][1]);
    if (ABS(xyz[i][2]) > max[2]) max[2] = ABS(xyz[i][2]);
  }

  printf("KDD %lf %lf %lf\n", max[0], max[1], max[2]);

  const double epsilon = 0.002;
  double upper = 1. + epsilon;
  double lower = 1. - epsilon;
  for (int i = 0; i < nc; i++) {
    xyz[i][0] /= max[0];
    xyz[i][1] /= max[1];
    xyz[i][2] /= max[2];
    double radius = 
           sqrt(xyz[i][0]*xyz[i][0]+xyz[i][1]*xyz[i][1]+xyz[i][2]*xyz[i][2]);
    if (radius < lower || radius > upper) 
      printf("BAD %lf %lf %lf: %lf\n", xyz[i][0], xyz[i][1], xyz[i][2], radius);
  }
  return xyz;
}

void write_new(int nc, point *xyz, const char *base)
{
  char filename[256];
  FILE *fp;

  sprintf(filename, "%s.graph", base);
  fp = fopen(filename, "w");
  fprintf(fp, "%d 0 000\n", nc);
  fclose(fp);

  sprintf(filename, "%s.coords", base);
  fp = fopen(filename, "w");
  for (int i = 0; i < nc; i++)
    fprintf(fp, "%lf %lf %lf\n", xyz[i][0], xyz[i][1], xyz[i][2]);
  fclose(fp);
}

void write_both(int nOne, point *xyzOne, int nTwo, point *xyzTwo) 
{
  FILE *fpg = fopen("both.graph", "w");
  FILE *fpc = fopen("both.coords", "w");
  int nc = nOne + nTwo;

  fprintf(fpg, "%d 0 010 2\n", nc);

  for (int i = 0; i < nOne; i++) {
    fprintf(fpc, "%lf %lf %lf\n", xyzOne[i][0], xyzOne[i][1], xyzOne[i][2]);
    fprintf(fpg, "1 0\n");
  }

  for (int i = 0; i < nTwo; i++) {
    fprintf(fpc, "%lf %lf %lf\n", xyzTwo[i][0], xyzTwo[i][1], xyzTwo[i][2]);
    fprintf(fpg, "0 1\n");
  }

  fclose(fpg);
  fclose(fpc);
}

int main()
{
  FILE *atmCoords = fopen("atm.coords", "r");
  FILE *atmGraph = fopen("atm.graph", "r");
  FILE *ocnCoords = fopen("ocn.coords", "r");
  FILE *ocnGraph = fopen("ocn.graph", "r");

  int natm = getNumCoords(atmGraph);
  int nocn = getNumCoords(ocnGraph);
  printf("KDD natm = %d nocn = %d\n", natm, nocn);

  point *scaledAtm = scaled_coordinates(natm, atmCoords);
  point *scaledOcn = scaled_coordinates(nocn, ocnCoords);

  write_new(natm, scaledAtm, "atm_new");
  write_new(nocn, scaledOcn, "ocn_new");

  write_both(natm, scaledAtm, nocn, scaledOcn);

  fclose(atmCoords); 
  fclose(atmGraph); 
  fclose(ocnCoords); 
  fclose(ocnGraph); 
  
  return 0;
}
