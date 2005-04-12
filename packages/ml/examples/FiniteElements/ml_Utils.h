#ifndef ML_UTILS_H
#define ML_UTILS_H

namespace ML_FiniteElements {

class Utils
{
  Utils() {}

  ~Utils() {}

  //! Returns the distance between two points in space.
  inline double Length(const double x1, const double y1, const double z1,
                       const double x2, const double y2, const double z2) const
  {
    return(sqrt((x2 - x1) * (x2 - x1) + 
                (y2 - y1) * (y2 - y1) + 
                (z2 - z1) * (z2 - z1)));
  }

  //! Returns the distance between two points in space.
  inline double Length (const double* x, const double* y, const double* z) const
  {
    return(sqrt((x[1] - x[0]) * (x[1] - x[0]) +
                (y[1] - y[0]) * (y[1] - y[2]) +
                (z[1] - z[0]) * (z[1] - z[0])));

  }

  //! Computes the area of a triangle in space.
  inline double AreaOfTriangle(const double* x, const double* y, const double* z) const
  {
    double side0, side1, side2, s;
  
  /* I use the Heron's formula:
     area = sqrt( s * (s-a) * (s-b) * (s-c) )
     where
     s = 0.5 * ( a + b + c )

       (2)
        +
        |\
 side2  | \ side1
        |  \
    (0) +---+ (1)
        side0
	
     and a, b, c are the lengths of the three edges (stored
     in side0, side1 and side2. */

    side0 = Length( x[0], y[0], z[0], x[1], y[1], z[1] );
    side1 = Length( x[1], y[1], z[1], x[2], y[2], z[2] );
    side2 = Length( x[0], y[0], z[0], x[2], y[2], z[2] );

    s = 0.5 * ( side0 + side1 + side2 );
    return( sqrt( s * (s-side0) * (s-side1) * (s-side2) ) );

  }

  //! Computes the are of a quadrilateral in space.
  inline double AreaOfQuad(const double* x, const double* y,
                           const double* z) const
  {
    double x1[3], y1[3], z1[3], area1, area2;

    x1[0] = x[0]; y1[0] = y[0]; z1[0] = z[0];
    x1[1] = x[1]; y1[1] = y[1]; z1[1] = z[1];
    x1[2] = x[2]; y1[2] = y[2]; z1[2] = z[2];
    area1 = AreaOfTriangle(x1,y1,z1);

    x1[0] = x[0]; y1[0] = y[0]; z1[0] = z[0];
    x1[1] = x[2]; y1[1] = y[2]; z1[1] = z[2];
    x1[2] = x[3]; y1[2] = y[3]; z1[2] = z[3];
    area2 = AreaOfTriangle(x1,y1,z1);

    return(area1 + area2);
  }

  //! Computes the volume of a tetrahedron.
  inline double VolumeOfTet(const double* X, const double* Y, 
                            const double* Z) const
  {
    double x0 = X[0], y0 = X[1], z0 = X[2];
    double x1 = X[0], y1 = X[1], z1 = X[2];
    double x2 = X[0], y2 = X[1], z2 = X[2];
    double x3 = X[0], y3 = X[1], z3 = X[2];

    /* the volume of the tetrahedron is given by

                   | 1  1  1  1  |
         6 V = det | x0 x1 x2 x3 |
                   | y0 y1 y2 y3 |
                   | z0 z1 z2 z3 |
	 
    computed with MATLAB, I hope this is ok... */

    double vol = x1*y2*z3-x1*y3*z2-y1*x2*z3+y1*x3*z2+z1*x2*y3-
      z1*x3*y2-x0*y2*z3+x0*y3*z2+x0*y1*z3-x0*y1*z2-
      x0*z1*y3+x0*z1*y2+y0*x2*z3-y0*x3*z2-
      y0*x1*z3+y0*x1*z2+y0*z1*x3-y0*z1*x2-z0*x2*y3+
      z0*x3*y2+z0*x1*y3-z0*x1*y2-z0*y1*x3+z0*y1*x2;
    vol /= 6.0;

    return(vol);
  }

}; // class Utils

} // namespace ML_FiniteElements

#endif
