#ifndef VECTOR_3D
#define VECTOR_3D

  class vector3d
  {
  public:
    // construction
    vector3d () : x(0.0), y(0.0), z(0.0)
    {}

    vector3d (double X, double Y, double Z) : x(X), y(Y), z(Z)
    {}

    explicit vector3d (double location[3])
      : x(location[0]), y(location[1]), z(location[2])
    {}

    vector3d (const vector3d& from) : x(from.x), y(from.y), z(from.z)
    {}

    double x, y, z;

    vector3d& operator= (const vector3d& from)
    {
      x = from.x;
      y = from.y;
      z = from.z;
      return *this;
    }

    bool operator== (const vector3d& from) const {
      return ( x == from.x && y == from.y && z == from.z );
    }

    bool operator!= (const vector3d& from) const {
      return ( x != from.x || y != from.y || z != from.z );
    }

    void set(double X, double Y, double Z) {
      x = X;
      y = Y;
      z = Z;
    }

    void set(double location[3]) {
      x = location[0];
      y = location[1];
      z = location[2];
    }

    vector3d& reverse() {
      x = -x;
      y = -y;
      z = -z;
      return *this;
    }

    vector3d operator- () const {
      vector3d tmp(x, y, z);
      return tmp *= -1.0;
    }

    inline vector3d& operator+= (const vector3d& from);
    inline vector3d& operator-= (const vector3d& from);
    inline vector3d& operator*= (double scalar);

    vector3d& operator/= (double scalar)
    {
      if ( scalar != 0.0 ) {
        x /= scalar;
        y /= scalar;
        z /= scalar;
      } else {
        x = HUGE_VAL;
        y = HUGE_VAL;
        z = HUGE_VAL;
      }

      return *this;
    }
    double length () const {
      return sqrt(x*x + y*y + z*z);
    }
    double squared_length () const {
      return x*x + y*y + z*z;
    }
    double dot (const vector3d& from) const {
      return x*from.x + y*from.y + z*from.z;
    }
    double normalize (double tolerance = 1e-06) {
      double mylength = length();
      if ( mylength > tolerance ) {
        x /= mylength;
        y /= mylength;
        z /= mylength;
      } else {
        mylength = 0.0;
      }

      return mylength;
    }
    inline vector3d cross (const vector3d& from) const;

    static vector3d plane_normal(const vector3d &v1, const vector3d &v2, const vector3d &v3) {
      vector3d v32 = v3;   v32 -= v2;
      vector3d v12 = v1;   v12 -= v2;
      return v32.cross(v12);
    }
  };

vector3d operator* (const vector3d& lhs, double scalar) {
  vector3d tmp(lhs);
  return tmp *= scalar;
}

vector3d operator* (double scalar, const vector3d& from)
{
  vector3d tmp(from);
  return tmp *= scalar;
}

vector3d operator/ (const vector3d& lhs, double scalar)
{
  if ( scalar != 0.0 ) {
    vector3d tmp(lhs);
    return tmp /= scalar;
  } else {
    return vector3d(HUGE_VAL, HUGE_VAL, HUGE_VAL);
  }
}
  vector3d operator- (const vector3d& lhs, const vector3d& rhs) {
    vector3d tmp(lhs);
    return tmp -= rhs;
  }


  vector3d operator+ (const vector3d& lhs, const vector3d& rhs) {
    vector3d tmp(lhs);
    return tmp += rhs;
  }

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
