// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_COORD_TRANSFORM_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_COORD_TRANSFORM_HPP_

#include "stk_search/SearchInterface.hpp"
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_search_util/ChangeOfBasis.hpp"
#include "stk_topology/topology.hpp"
#include "stk_expreval/Eval.hpp"

namespace stk::search {

using CoordTransformInterface = stk::search::CoordTransformBaseInterface<spmd::EntityKeyPair>;


//-----------------------------------------------------------------------------
// Basic transfers
class CoordTransformIdentity : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    void transform(const EntityKey&, std::vector<double>& /*coords*/) override {};
};

class CoordTransformAddZ : public CoordTransformInterface
{
  public:

    CoordTransformAddZ(const std::string& expression = "")
    {
      m_eval.setExpression(expression);

      m_eval.parse();

      STK_ThrowRequireMsg(m_eval.getParseStatus(),"Expression provided for CoordTransformAddZ was not parsed successfully!");
    }

    using CoordTransformInterface::EntityKey;

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      STK_ThrowAssertMsg(coords.size() == 2, "expected 2D coordinates");
      
      bind_variables(coords);

      coords.push_back(m_eval.evaluate());
    };

  private:

  void bind_variables(std::vector<double>& coords)
  {
    std::vector<std::string> dependentVars = m_eval.get_dependent_variable_names();

    if (dependentVars.size() > 1) {
      STK_ThrowErrorMsg("For CoordTransformAddZ, z is the only permissible dependent variable!");
    }
   
    if (dependentVars.size() == 1) {
      STK_ThrowRequireMsg(m_eval.is_dependent_variable("z"), "For CoordTransformAddZ expression, z is the only permissible dependent variable!");
    }
    
    std::vector<std::string> independentVars = m_eval.get_independent_variable_names();

    if (independentVars.size() > 2) {
      STK_ThrowErrorMsg("For CoordTransformAddZ, x and y are the only permissible independent variables!");
    }

    if (independentVars.size() == 2) {
      bool correctIndepentVariables = m_eval.is_independent_variable("x") && m_eval.is_independent_variable("y");
      STK_ThrowRequireMsg(correctIndepentVariables, "For CoordTransformAddZ, x and y are the only permissible independent variables!");

      m_eval.bindVariable("x", coords[0]);
      m_eval.bindVariable("y", coords[1]);

    }

    if (independentVars.size() == 1) {
      if (m_eval.is_independent_variable("x")) {
        m_eval.bindVariable("x", coords[0]);
      }
      else if (m_eval.is_independent_variable("y")) {
        m_eval.bindVariable("y", coords[1]);
      }
      else {
        STK_ThrowErrorMsg("For CoordTransformAddZ, x and y are the only permissible independent variables!");
      }
    }
  }

    stk::expreval::Eval m_eval;
};

class CoordTransformRemoveZ : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      STK_ThrowAssertMsg(coords.size() == 3, "expected 3D coordinates");
      coords.resize(2);
    };
};

class CoordTransformTranslate : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    CoordTransformTranslate(const std::array<double, 3>& translate) :
      m_translate(translate)
    {}

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      STK_ThrowAssertMsg(coords.size() == 3 || (coords.size() == 2 && m_translate[2] == 0),
                         "Applying translation to 2D coordinates but 3D translation was specified");
      for (size_t i=0; i < coords.size(); ++i)
      {
        coords[i] += m_translate[i];
      }
    };

  private:
    std::array<double, 3> m_translate;
};

class CoordTransformScale : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    CoordTransformScale(const std::array<double, 3>& scale) :
      m_scale(scale)
    {}

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      STK_ThrowAssertMsg(coords.size() == 3 || (coords.size() == 2 && m_scale[2] == 0),
                         "Applying scale to 2D coordinates but 3D scale was specified");
      for (size_t i=0; i < coords.size(); ++i)
      {
        coords[i] *= m_scale[i];
      }
    };

  private:
    std::array<double, 3> m_scale;
};

class CoordTransformPermute : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    CoordTransformPermute(const std::array<unsigned, 3>& dests) :
      m_dests(dests)
    {
      STK_ThrowRequireMsg(dests[0] != dests[1] && dests[1] != dests[2] && dests[1] != dests[2],
                          "CoordTransformPermute must have unique destination indices");
      STK_ThrowRequireMsg((dests[0] < 3) &&
                          (dests[1] < 3) &&
                          (dests[2] == std::numeric_limits<unsigned>::max() || (dests[2] < 3)),
                          "CoordTransformPermute destination indices must be in range [0, 2]");
    }

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      STK_ThrowAssertMsg(coords.size() == 3 || (coords.size() == 2 && m_dests[2] == std::numeric_limits<unsigned>::max()),
                         "Applying coord permutation to 2D coordinates but a 3D permutation was specified");

      std::array<double, 3> coords_copy = {coords[0], coords[1], coords.size() == 3 ? coords[2] : 0};
      for (size_t i=0; i < coords.size(); ++i)
      {
        STK_ThrowAssertMsg(m_dests[i] < coords.size(), "destination index out of range in CoordTransformPermute");
        coords[m_dests[i]] = coords_copy[i];
      }
    };

  private:
    std::array<unsigned, 3> m_dests;
};

class CoordTransformRotateZ : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    CoordTransformRotateZ(const std::array<double, 3>& new_z, bool forward=true) :
      m_change_of_basis(impl::compute_rotation_matrix(new_z)),
      m_forward(forward)
    {}

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      STK_ThrowAssertMsg(coords.size() == 3,
                         "RotateZ transform can only be applied to 3D coordinates");
      std::array<double, 3> coords_copy = {coords[0], coords[1], coords[2]};
      if (m_forward)
      {
        coords_copy = m_change_of_basis.change_basis_forward(coords_copy);
      } else
      {
        coords_copy = m_change_of_basis.change_basis_reverse(coords_copy);
      }
      coords[0] = coords_copy[0];
      coords[1] = coords_copy[1];
      coords[2] = coords_copy[2];
    };

  private:
    impl::ChangeOfBasis m_change_of_basis;
    bool m_forward;
};

class CoordTransformCylindricalCoordinates : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      double x = coords[0];
      double y = coords[1];

      coords[0] = std::sqrt(x*x + y*y);
      coords[1] = std::atan2(y, x);
    };
};

class CoordTransformCartesianCoordinates : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      double r = coords[0];
      double theta = coords[1];

      coords[0] = r*std::cos(theta);
      coords[1] = r*std::sin(theta);
    };
};

class CoordTransformPeriodicSymmetry : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    CoordTransformPeriodicSymmetry(int dim, double xmin, double xmax) :
      m_dim(dim),
      m_xmin(xmin),
      m_xmax(xmax)
    {}

    void transform(const EntityKey&, std::vector<double>& coords) override
    {
      double rng = m_xmax - m_xmin;
      double x = coords[m_dim];

      double x_prime = x - m_xmin;
      int period = std::floor(x_prime / rng);
      double delta_x = x_prime - period*rng;

      if (period % 2 == 0)
      {
        coords[m_dim] = m_xmin + delta_x;

      } else
      {
        coords[m_dim] = m_xmax - delta_x;
      }
    };

  private:
    int m_dim;
    double m_xmin;
    double m_xmax;
};

//-----------------------------------------------------------------------------
// All-in-one transfers that combine the basic transfers above

// For transferring to a axisymmetric body from a portion of that body defined by theta \in  [theta_min, theta_max].
// For example, the source geometry can be a wedge and the destination mesh can be the shape obtained by revolving
// the body about the point of the wedge.  In this case, theta_min and theta_max describe the angle of the
// wedge on the source mesh (when converted to cylindrical coordinates), axis is a vector along the tip of the wedge,
// and pt is used to translate the destination mesh, which can be useful if the axis of the destination mesh is offset
// from the source mesh.
// Another example is a source geometry that is a quarter ellipse and the destination geometry a full ellipse.
// In this case, theta_min and theta_max describe the portion of the cylindrical coordinate space where the
// quarter ellipse is, axis is a vector along the axis of the cylinder, and pt is used to translate the destination mesh,
// to align the axis of its ellipse with that of the source geometry.
class CoordTransformAxisymmetric3D : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    CoordTransformAxisymmetric3D(double theta_min, double theta_max,
                                 const std::array<double, 3>& axis = {0, 0, 1},
                                 const std::array<double, 3>& pt = {0, 0, 0}) :
      m_translate(pt),
      m_rotate(axis),
      m_cartToCyl(),
      m_periodic(1, theta_min, theta_max),
      m_cylToCart()
    {}

    void transform(const EntityKey& key, std::vector<double>& coords) override
    {
      m_translate.transform(key, coords);
      m_rotate.transform(key, coords);
      m_cartToCyl.transform(key, coords);
      m_periodic.transform(key, coords);
      m_cylToCart.transform(key, coords);
    };

  private:
    CoordTransformTranslate m_translate;
    CoordTransformRotateZ m_rotate;
    CoordTransformCylindricalCoordinates m_cartToCyl;
    CoordTransformPeriodicSymmetry m_periodic;
    CoordTransformCartesianCoordinates m_cylToCart;
};


// For transferring between a 2D mesh to a 3D mesh that is axisymmetric, and the 2D
// mesh is a slice of the 3D geometry at some theta = theta_star
// It does this by translating and rotating the 3D mesh so its axis is aligned with
// the 2D mesh axis, then converting to cylindrical coordinates and removing theta
// so only (r, z) remain.
class CoordTransformAxisymmetric2D : public CoordTransformInterface
{
  public:
    using CoordTransformInterface::EntityKey;

    CoordTransformAxisymmetric2D(const std::array<double, 3>& axis = {0, 0, 1},
                                 const std::array<double, 3>& pt = {0, 0, 0}) :
      m_translate(pt),
      m_rotate(axis),
      m_cartToCyl(),
      m_permute({0, 2, 1}),
      m_removeZ()
    {}

    void transform(const EntityKey& key, std::vector<double>& coords) override
    {
      m_translate.transform(key, coords);
      m_rotate.transform(key, coords);
      m_cartToCyl.transform(key, coords);
      m_permute.transform(key, coords);
      m_removeZ.transform(key, coords);
    };

  private:
    CoordTransformTranslate m_translate;
    CoordTransformRotateZ m_rotate;
    CoordTransformCylindricalCoordinates m_cartToCyl;
    CoordTransformPermute m_permute;
    CoordTransformRemoveZ m_removeZ;
};

}

#endif
