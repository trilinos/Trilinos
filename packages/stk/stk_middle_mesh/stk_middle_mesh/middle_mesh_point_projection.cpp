#include "middle_mesh_point_projection.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "mesh_entity.hpp"


namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {


mesh::FieldPtr<utils::Point> MiddleMeshPointProjection::projection_onto_mesh1(std::shared_ptr<XiCoordinates> xiCoords)
{
  int npts = std::max(xiCoords->get_xi_coords(mesh::MeshEntityType::Triangle).size(),
                      xiCoords->get_xi_coords(mesh::MeshEntityType::Quad).size());
  auto projectedXiCoordsPtr = mesh::create_field<utils::Point>(m_meshIn, mesh::FieldShape(0, 0, npts), 1);

  auto& projectedXiCoords   = *projectedXiCoordsPtr;
  auto& meshInToMesh1Elements = *(m_relationalData->meshInElementsToMesh1Elements);
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);
  for (auto elIn : m_meshIn->get_elements())
    if (elIn)
    {
      mesh::MeshEntityPtr el1 = meshInToMesh1Elements(elIn, 0, 0);
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
      int nverts = mesh::get_downward(elIn, 0, verts.data());

      std::array<utils::Point, mesh::MAX_DOWN> vertXiCoordsOnMesh1, vertXYZCoordsOnMesh1;

      for (int i=0; i < nverts; ++i)
      {
        predicates::impl::PointRecord vertRecord = vertsInClassOnMesh1(verts[i], 0, 0);
        if (vertRecord.el != el1)
        {
          vertRecord = m_classifier->classify_onto(vertRecord, el1);
        }

        vertXiCoordsOnMesh1[i] = m_classifier->compute_xi_coords(vertRecord);
        if (el1->get_type() == mesh::MeshEntityType::Quad)
          vertXYZCoordsOnMesh1[i] = m_classifier->compute_xyz_coords(vertRecord);
      }

      const auto& xiCoordsEl = xiCoords->get_xi_coords(elIn->get_type());
      for (size_t i=0; i < xiCoordsEl.size(); ++i)
      {
        std::pair<double, double> xiRangeIn = xiCoords->get_xi_coord_range(elIn->get_type());
        std::pair<double, double> xiRange1  = xiCoords->get_xi_coord_range(el1->get_type());
        utils::Point ptXi = mesh::convert_xi_coords_from_range(xiRangeIn.first, xiRangeIn.second, xiCoordsEl[i]);
        if (el1->get_type() == mesh::MeshEntityType::Triangle)
        {
          utils::Point ptXi1 = mesh::interpolate(elIn->get_type(), vertXiCoordsOnMesh1.data(), ptXi);
          projectedXiCoords(elIn, i, 0) = mesh::convert_xi_coords_to_range(xiRange1.first, xiRange1.second, ptXi1);

        } else if (el1->get_type() == mesh::MeshEntityType::Quad)
        {
          utils::Point ptXyz     = mesh::interpolate(elIn->get_type(), vertXYZCoordsOnMesh1.data(), ptXi);
          utils::Point ptXiGuess = mesh::interpolate(elIn->get_type(), vertXiCoordsOnMesh1.data(), ptXi);
          utils::Point ptXi1     = compute_nearest_point_on_quad(el1, ptXyz, ptXiGuess);          
          projectedXiCoords(elIn, i, 0) = mesh::convert_xi_coords_to_range(xiRange1.first, xiRange1.second, ptXi1);
        }
      }
    }

  return projectedXiCoordsPtr;
}


mesh::FieldPtr<utils::Point> MiddleMeshPointProjection::projection_onto_mesh2(std::shared_ptr<XiCoordinates> xiCoords)
{
  int npts = std::max(xiCoords->get_xi_coords(mesh::MeshEntityType::Triangle).size(),
                      xiCoords->get_xi_coords(mesh::MeshEntityType::Quad).size());
  auto projectedXiCoordsPtr = mesh::create_field<utils::Point>(m_meshIn, mesh::FieldShape(0, 0, npts), 1);

  auto& projectedXiCoords   = *projectedXiCoordsPtr;
  auto& meshInToMesh2Elements = *(m_relationalData->meshInElementsToMesh2Elements);
  for (auto elIn : m_meshIn->get_elements())
    if (elIn)
    {
      mesh::MeshEntityPtr el2 = meshInToMesh2Elements(elIn, 0, 0);
      const auto& xiCoordsEl = xiCoords->get_xi_coords(elIn->get_type());
      for (size_t i=0; i < xiCoordsEl.size(); ++i)
      {
        std::pair<double, double> xiRangeIn = xiCoords->get_xi_coord_range(elIn->get_type());
        std::pair<double, double> xiRange2  = xiCoords->get_xi_coord_range(el2->get_type());
        utils::Point ptXiIn = mesh::convert_xi_coords_from_range(xiRangeIn.first, xiRangeIn.second, xiCoordsEl[i]);

        utils::Point ptXyzIn = mesh::compute_coords_from_xi_3d(elIn, ptXiIn);
        predicates::impl::PointRecord record = m_classifier->classify_reverse(el2, ptXyzIn);

        if (el2->get_type() == mesh::MeshEntityType::Triangle)
        {
          utils::Point ptXi2 = m_classifier->compute_xi_coords(record);
          projectedXiCoords(elIn, i, 0) = mesh::convert_xi_coords_to_range(xiRange2.first, xiRange2.second, ptXi2);
        } else
        {
          utils::Point ptXyzOnTriangle  = m_classifier->compute_xyz_coords(record);
          utils::Point ptXiGuess        = m_classifier->compute_xi_coords(record);
          utils::Point ptXi2 = compute_nearest_point_on_quad(el2, ptXyzOnTriangle, ptXiGuess);
          projectedXiCoords(elIn, i, 0) = mesh::convert_xi_coords_to_range(xiRange2.first, xiRange2.second, ptXi2);
        }
      }
    }

  return projectedXiCoordsPtr;
}

utils::Point MiddleMeshPointProjection::compute_nearest_point_on_quad(mesh::MeshEntityPtr el, const utils::Point& ptXyz, const utils::Point& ptXiGuess)
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
  m_gaussNewton.solve(f, jacFunc, m_xi0);
  utils::Point ptXi(m_xi0[0], m_xi0[1]);
  
  return {m_xi0[0], m_xi0[1]};
}




} // namespace 
} // namespace 
} // namespace 
} // namespace 