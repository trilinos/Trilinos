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

point *spherical_coordinates(int nc, FILE *fp_coords)
{
  double max[3] = {0., 0., 0.};
  point *xyz = new point[nc];

  // read coords and convert to spherical coordinates
  const double epsilon = 0.002;
  double upper = 1. + epsilon;
  double lower = 1. - epsilon;
  double degPerRad = 180. / M_PI;

  for (int i = 0; i < nc; i++) {
    double x, y, z;
    fscanf(fp_coords, "%lf %lf %lf\n", &x, &y, &z);
    double radius = sqrt(x*x+y*y+z*z);
    xyz[i][0] = atan2(y, x) * degPerRad;
    xyz[i][1] = asin(z/radius) * degPerRad;
    xyz[i][2] = radius;
  }

  return xyz;
}

void write_new2D(int nc, point *xyz, const char *base)
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
    fprintf(fp, "%lf %lf\n", xyz[i][0], xyz[i][1]);
  fclose(fp);
}

void write_both2D(int nOne, point *xyzOne, int nTwo, point *xyzTwo) 
{
  FILE *fpg = fopen("both_spherical.graph", "w");
  FILE *fpc = fopen("both_spherical.coords", "w");
  int nc = nOne + nTwo;

  fprintf(fpg, "%d 0 010 2\n", nc);

  for (int i = 0; i < nOne; i++) {
    fprintf(fpc, "%lf %lf\n", xyzOne[i][0], xyzOne[i][1]);
    fprintf(fpg, "1 0\n");
  }

  for (int i = 0; i < nTwo; i++) {
    fprintf(fpc, "%lf %lf\n", xyzTwo[i][0], xyzTwo[i][1]);
    fprintf(fpg, "0 1\n");
  }

  fclose(fpg);
  fclose(fpc);
}

int main()
{
  FILE *atmCoords = fopen("atm_new.coords", "r");
  FILE *atmGraph = fopen("atm_new.graph", "r");
  FILE *ocnCoords = fopen("ocn_new.coords", "r");
  FILE *ocnGraph = fopen("ocn_new.graph", "r");

  int natm = getNumCoords(atmGraph);
  int nocn = getNumCoords(ocnGraph);
  printf("KDD natm = %d nocn = %d\n", natm, nocn);

  point *sphericalAtm = spherical_coordinates(natm, atmCoords);
  point *sphericalOcn = spherical_coordinates(nocn, ocnCoords);

  write_new2D(natm, sphericalAtm, "atm_spherical");
  write_new2D(nocn, sphericalOcn, "ocn_spherical");

  write_both2D(natm, sphericalAtm, nocn, sphericalOcn);

  fclose(atmCoords); 
  fclose(atmGraph); 
  fclose(ocnCoords); 
  fclose(ocnGraph); 
  
  return 0;
}
