// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#if HAVE_OPENNURBS

#include <percept/mesh/geometry/stk_geom/SplineFit.hpp>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <utility>
#include <string>

  namespace geom {

    // for Mathematica, you can choose to output with precision specified 
    //   - see BSplineFitting.nb in this directory
    //static std::string s_mm_prec="`15";
    static std::string s_mm_prec="";
    bool SplineFit::s_debug_print = false;

    // for cout
    static int s_precision=15;

    inline void replace(std::string &str, const std::string &find_what, const std::string &replace_with)
    {
      std::string::size_type pos = 0;
      while((pos = str.find(find_what, pos)) != std::string::npos)
        {
          str.erase(pos, find_what.length());
          str.insert(pos, replace_with);
          pos += replace_with.length();
        }
    }

    inline std::string convert_to_mm(double d)
    {
      std::ostringstream str;
      str << std::setprecision(s_precision);
      str << d;
      std::string ret = str.str();
      replace(ret, "e", "*^");
      replace(ret, "E", "*^");
      ret = ret + s_mm_prec;
      return ret;
    }

    std::ostream& operator<<(std::ostream& out,  const Vectors2D& pts)
    {
      out << std::setprecision(s_precision);
      out << "{";
      for (unsigned i = 0; i < pts.size(); i++)
        {
          out << " " << pts[i] << (i < pts.size()-1?", ":"");
          if ((i+1) % 8 == 0) out << "\n";
        }
      out << "};";
      return out;
    }
    std::ostream& operator<<(std::ostream& out,  const Vectors3D& pts)
    {
      out << std::setprecision(s_precision);
      out << "{";
      for (unsigned i = 0; i < pts.size(); i++)
        {
          out << " " << pts[i] << (i < pts.size()-1?", ":"");
          if ((i+1) % 8 == 0) out << "\n";
        }
      out << "};";
      return out;
    }
    std::ostream& operator<<(std::ostream& out,  const Points3D& pts)
    {
      out << std::setprecision(s_precision);
      out << "{";
      for (unsigned i = 0; i < pts.size(); i++)
        {
          out << " " << pts[i] << (i < pts.size()-1?", ":"");
          if ((i+1) % 8 == 0) out << "\n";
        }
      out << "};";
      return out;
    }
    std::ostream& operator<<(std::ostream& out,  const Vector2D& pt)
    {
      out << std::setprecision(s_precision);
      out << " {" << convert_to_mm(pt[0])  << ", " << convert_to_mm(pt[1])  << "}";
      return out;
    }
    std::ostream& operator<<(std::ostream& out,  const Point3D& pt)
    {
      out << std::setprecision(s_precision);
      out << " {" << convert_to_mm(pt[0])  << ", " << convert_to_mm(pt[1])  << ", " << convert_to_mm(pt[2])  << "}";
      return out;
    }
    std::ostream& operator<<(std::ostream& out,  const Vector3D& pt)
    {
      out << std::setprecision(s_precision);
      out << " {" << convert_to_mm(pt[0])  << ", " << convert_to_mm(pt[1])  << ", " << convert_to_mm(pt[2])  << "}";
      return out;
    }
    std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec)
    {
      out << std::setprecision(s_precision);
      out << "{";
      for (unsigned i = 0; i < vec.size(); i++)
        {
          out << " " << convert_to_mm(vec[i])  << (i<vec.size()-1?",  ":"");
          if ((i+1) % 8 == 0) out << "\n";
        }
      out << "};";
      return out;
    }

  }

#endif
