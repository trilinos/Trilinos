#ifndef VECTOR3D
#define VECTOR3D

  class vector3d
  {
  public:
    // construction
    vector3d ();
    vector3d (double X, double Y, double Z);
    explicit vector3d (double location[3]);
    vector3d (const vector3d& from);

    double x, y, z;

    vector3d& operator= (const vector3d& from);
    bool operator== (const vector3d& from) const;
    bool operator!= (const vector3d& from) const;
    void set(double X, double Y, double Z);
    void set(double location[3]);
    vector3d& reverse();

    vector3d operator- () const;

    vector3d& operator+= (const vector3d& from);
    vector3d& operator-= (const vector3d& from);
    vector3d& operator*= (double scalar);
    vector3d& operator/= (double scalar);

    double length () const;
    double squared_length () const;
    double dot (const vector3d& from) const;
    double normalize (double tolerance = 1e-06);
    vector3d cross (const vector3d& from) const;
    static vector3d plane_normal(const vector3d &v1, const vector3d &v2, const vector3d &v3);
  };

  vector3d operator* (double scalar, const vector3d& vec);
  vector3d operator* (const vector3d& vec, double scalar);
  vector3d operator/ (const vector3d& vec, double scalar);

  vector3d operator+ (const vector3d& vec1, const vector3d& vec2);
  vector3d operator- (const vector3d& vec1, const vector3d& vec2);

  //----------------------------------------------------------------------------
  inline vector3d vector3d::cross (const vector3d& from) const
  {
    return vector3d(y * from.z - z * from.y,
		    z * from.x - x * from.z,
		    x * from.y - y * from.x);
  }
  //----------------------------------------------------------------------------
  inline vector3d& vector3d::operator+= (const vector3d& from)
  {
    x += from.x;
    y += from.y;
    z += from.z;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline vector3d& vector3d::operator-= (const vector3d& from)
  {
    x -= from.x;
    y -= from.y;
    z -= from.z;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline vector3d& vector3d::operator*= (double scalar)
  {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
  }

#endif
