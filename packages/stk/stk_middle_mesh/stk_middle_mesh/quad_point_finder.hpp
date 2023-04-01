#ifndef STK_MIDDLE_MESH_QUAD_POINT_FINDER
#define STK_MIDDLE_MESH_QUAD_POINT_FINDER

#include "mesh.hpp"
#include "gauss_newton.hpp"


namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class QuadPointFinder
{
  public:
    QuadPointFinder() :
      m_gaussNewton(3, 2, 1e-13, 100),
      m_xi0(2)
    {}

    // computes the xi coordinates of the point on the quad closest to ptXyz.  Uses ptXiGuess as
    // the initial guess for the iterative method.
    utils::Point compute_nearest_point_on_quad(mesh::MeshEntityPtr el, const utils::Point& ptXyz, const utils::Point& ptXiGuess)
    {
      bool allowFailure = false;
      bool didFail = false;
      return compute_nearest_point_on_quad(el, ptXyz, ptXiGuess, allowFailure, didFail);
    }

    utils::Point compute_nearest_point_on_quad(mesh::MeshEntityPtr el, const utils::Point& ptXyz, const utils::Point& ptXiGuess,
                                               bool allowFailure, bool& didFail)
    {
      assert(el->get_type() == mesh::MeshEntityType::Quad);

      std::array<mesh::MeshEntityPtr, 4> verts;
      mesh::get_downward(el, 0, verts.data());

      std::array<utils::Point, 4> pts;
      for (int i=0; i < 4; ++i)
        pts[i] = verts[i]->get_point_orig(0);

      auto f = [&](const std::vector<double>& quadXi, std::vector<double>& residuals)
      {
        //TODO: this requires calling getDownward every time
        utils::Point currPoint = compute_quad_coords_from_xi_3d(el, utils::Point(quadXi[0], quadXi[1]));

        for (int i=0; i < 3; ++i)
          residuals[i] = ptXyz[i] - currPoint[i];
      };

      auto jacFunc = [&](const std::vector<double>& quadXi, utils::impl::Matrix<double>& jac)
      {
        std::array<double, 2> xiVals, xiDerivs, etaVals, etaDerivs;
        mesh::compute_lagrange_vals(quadXi[0], xiVals.data());
        mesh::compute_lagrange_derivs(quadXi[0], xiDerivs.data());

        mesh::compute_lagrange_vals(quadXi[1], etaVals.data());
        mesh::compute_lagrange_derivs(quadXi[1], etaDerivs.data());

        for (int i=0; i < 3; ++i)
        {
          jac(i, 0) = -(xiDerivs[0]*etaVals[0]*pts[0][i] + xiDerivs[1]*etaVals[0]*pts[1][i] +
                        xiDerivs[1]*etaVals[1]*pts[2][i] + xiDerivs[0]*etaVals[1]*pts[3][i]);
          jac(i, 1) = -(xiVals[0]*etaDerivs[0]*pts[0][i] + xiVals[1]*etaDerivs[0]*pts[1][i] +
                        xiVals[1]*etaDerivs[1]*pts[2][i] + xiVals[0]*etaDerivs[1]*pts[3][i]);
        }
      };


      m_xi0[0] = ptXiGuess.x;
      m_xi0[1] = ptXiGuess.y;
      didFail = !(m_gaussNewton.solve(f, jacFunc, m_xi0, allowFailure));

      if (didFail && !allowFailure)
        throw std::runtime_error("Gauss-Newton failed to converged when finding closest point on quad");
      
      return {m_xi0[0], m_xi0[1]};
    }

  private:
    utils::impl::GaussNewton m_gaussNewton;
    std::vector<double> m_xi0;  

};

}
}
}
}

#endif