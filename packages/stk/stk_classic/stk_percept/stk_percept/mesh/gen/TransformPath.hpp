#ifndef stk_percept_mesh_gen_TransformPath_hpp
#define stk_percept_mesh_gen_TransformPath_hpp

#include "SweepMesher.hpp"
#include <cmath>

namespace stk_classic
{
  namespace percept
  {
    namespace util
    {

      class TransformPath : public Transform
      {
        typedef boost::array<double,3> Coord;
        Coord m_from;
        Coord m_from_dir;
        Coord m_to;
        Coord m_to_dir;
        Coord m_rotation_axis;
        Coord m_rotation_origin;
        double m_theta;
        Coord m_diff;

        void normalize(Coord& vec)
        {
          double norm = std::sqrt(vec[0]*vec[0]+
                             vec[1]*vec[1]+
                             vec[2]*vec[2]);
          vec[0] /= norm;
          vec[1] /= norm;
          vec[2] /= norm;
        }
        void cross(const Coord& a, const Coord& b, Coord& axb)
        {
          axb[0] = (a[1]*b[2]-a[2]*b[1]);
          axb[1] = -(a[0]*b[2]-a[2]*b[0]);
          axb[2] = (a[0]*b[1]-a[1]*b[0]);
        }

        // from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
        void rotatePointAboutAxis(const Coord& rotation_axis, const Coord& rotation_origin, const Coord& point_to_rotate, const double angleInRadians,
                                  Coord& rotated_point)
        {
          double a = rotation_origin[0], b= rotation_origin[1], c= rotation_origin[2];
          double u = rotation_axis[0], v= rotation_axis[1], w= rotation_axis[2];
          double x = point_to_rotate[0], y = point_to_rotate[1], z = point_to_rotate[2];
          double cosAngle = std::cos(angleInRadians);
          double sinAngle = std::sin(angleInRadians);
          double u2 = u*u;
          double v2 = v*v;
          double w2 = w*w;
          double l2 = u2 + v2 + w2;
          double l = sqrt(l2);
          rotated_point[0] = (a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
                              + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosAngle
                              + l*(-c*v + b*w - w*y + v*z)*sinAngle)/l2;
          rotated_point[1] = (b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
                              + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosAngle
                              + l*(c*u - a*w + w*x - u*z)*sinAngle)/l2;
          rotated_point[2] = (c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
                              + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) +  (u2 + v2)*z)*cosAngle
                              + l*(-b*u + a*v - v*x + u*y)*sinAngle)/l2;

        }
        void rotatePointAboutAxis(const Coord& point_to_rotate, Coord& rotated_point)
        {
          rotatePointAboutAxis(m_rotation_axis, m_rotation_origin, point_to_rotate, m_theta, rotated_point);
        }


      public:

        /// Given some points on a plane ("from"), we want to move them to a new plane by rotation and translation
        /// The initial plane is defined by an origin ([in] from), and it's normal ([in] from_dir).
        /// The final plane is defined by origin/normal: [in] to, to_dir
        /// The translation delta is just the vector {to - from}.
        /// The algorithm rotates the points in the initial plane to the new plane's direction, then does the translate.
        /// Note that the points on the initial plane don't have to lie in a plane, just easier to visualize that way.
        TransformPath(const Coord& from, const Coord& from_dir, 
                  const Coord& to,   const Coord& to_dir) : m_from(from), m_from_dir(from_dir),
                                                            m_to(to) , m_to_dir(to_dir)
        {
          m_diff[0] = m_to[0]-m_from[0];
          m_diff[1] = m_to[1]-m_from[1];
          m_diff[2] = m_to[2]-m_from[2];

          // normalize
          normalize(m_from_dir);
          normalize(m_to_dir);
          // angle between the two direction vectors
          m_theta = acos(m_from_dir[0]*m_to_dir[0]+
                         m_from_dir[1]*m_to_dir[1]+
                         m_from_dir[2]*m_to_dir[2]);
          //std::cout<< "m_theta= " << m_theta << std::endl;

          cross(m_from_dir, m_to_dir, m_rotation_axis);
          m_rotation_origin = m_from;
        }

        using Transform::operator();
        virtual Coord  operator()(const Coord& x)
        {
          Coord y;
          operator()(x, y);
          return y;
        }
        virtual void operator()(const Coord& x, Coord& y) 
        {
          rotatePointAboutAxis(x, y);
          y[0] += m_diff[0];
          y[1] += m_diff[1];
          y[2] += m_diff[2];
        }
      };

    }//namespace utils
  }//namespace percept
}//namespace stk_classic

#endif
