#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_search/BoundingBox.hpp>
#include "stk_search/CoarseSearch.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/quad_point_finder.hpp"
#include "stk_middle_mesh/quad_metrics.hpp"
#include "stk_middle_mesh_util/create_stk_mesh.hpp"
#include "stk_middle_mesh_util/exodus_writer.hpp"

#include <string>

using namespace stk::middle_mesh;

using EntityIdentifier = stk::search::IdentProc<int, int>;
using BoundingBoxAndId = std::pair<stk::search::Box<double>, EntityIdentifier>;
using BoundingPointAndId = std::pair<stk::search::Point<double>, EntityIdentifier>;

BoundingBoxAndId create_element_bounding_box(mesh::MeshEntityPtr entity)
{
  assert(mesh::get_type_dimension(entity->get_type()) == 2);

  double expansionFactor = 0.15;

  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  int nverts = mesh::get_downward(entity, 0, verts.data());
  double minVal = std::numeric_limits<double>::min();
  double maxVal = std::numeric_limits<double>::max();
  std::array<double, 3> upperBound = {minVal, minVal, minVal}, lowerBound = {maxVal, maxVal, maxVal};
  for (int i=0; i < nverts; ++i)
  {
    utils::Point pt = verts[i]->get_point_orig(0);

    for (int j=0; j < 3; ++j)
    {
      upperBound[j] = std::max(upperBound[j], pt[j]);
      lowerBound[j] = std::min(lowerBound[j], pt[j]);
    }
  }

  for (int i=0; i < 3; ++i)
  {
    double dist = expansionFactor*(upperBound[i] - lowerBound[i]);
    upperBound[i] += dist/2;
    lowerBound[i] -= dist/2;
  }

  stk::search::Point<double> minCorner(lowerBound[0], lowerBound[1], lowerBound[2]),
                             maxCorner(upperBound[0], upperBound[1], upperBound[2]);

  return std::make_pair(stk::search::Box(minCorner, maxCorner), EntityIdentifier(entity->get_id(), 0));
}

std::vector<BoundingBoxAndId> get_element_bounding_boxes(std::shared_ptr<mesh::Mesh> mesh)
{
  std::vector<BoundingBoxAndId> boxes;
  for (auto el : mesh->get_elements())
    if (el)
    {
      boxes.push_back(create_element_bounding_box(el));
    }

  return boxes;
}

std::vector<BoundingPointAndId> get_vert_bounding_spheres(std::shared_ptr<mesh::Mesh> mesh)
{
  std::vector<BoundingPointAndId> spheres;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      spheres.emplace_back(stk::search::Point<double>(pt.x, pt.y, pt.z),
                           EntityIdentifier(vert->get_id(), 0));;                           
    }

  return spheres;
}

std::vector<std::pair<EntityIdentifier, EntityIdentifier>> compute_possible_intersections(std::shared_ptr<mesh::Mesh> mesh1,
                                                                std::shared_ptr<mesh::Mesh> mesh2)
{
  auto senderBoxes = get_element_bounding_boxes(mesh1);
  auto receiverSpheres = get_vert_bounding_spheres(mesh2);
  std::vector<std::pair<EntityIdentifier, EntityIdentifier>> intersections;

  stk::search::coarse_search(senderBoxes, receiverSpheres, stk::search::KDTREE, MPI_COMM_WORLD, intersections);

  auto cmp = [](const std::pair<EntityIdentifier, EntityIdentifier>& lhs, const std::pair<EntityIdentifier, EntityIdentifier>& rhs)
  {
    if (lhs.second != rhs.second)
      return lhs.second.id() < rhs.second.id();
    else
      return lhs.first.id() < rhs.first.id();
  };

  std::sort(intersections.begin(), intersections.end(), cmp);

  return intersections;
}

void check_all_points_are_matched(std::shared_ptr<mesh::Mesh> recvMesh,
                                  const std::vector<std::pair<EntityIdentifier, EntityIdentifier>>& intersections)
{
  std::vector<int> recvPoints(intersections.size());
  for (size_t i=0; i < intersections.size(); ++i)
    recvPoints[i] = intersections[i].second.id();

  for (auto& vert : recvMesh->get_vertices())
    if (vert)
    {
      if (!std::binary_search(recvPoints.begin(), recvPoints.end(), vert->get_id()))
        throw std::runtime_error("unmatched receive point");
    }
}

struct PointRelation
{
  int senderId;
  int destId;
  utils::Point senderPtXi;
};


class UniqueIntersectionFinder
{
  public:
    UniqueIntersectionFinder(const std::vector<std::pair<EntityIdentifier, EntityIdentifier>>& intersections,
                             std::shared_ptr<mesh::Mesh> meshSend,
                             std::shared_ptr<mesh::Mesh> meshRecv) :
      m_intersections(intersections),
      m_meshSend(meshSend),
      m_meshRecv(meshRecv)
    {}

    std::vector<PointRelation> find()
    {
      std::vector<PointRelation> uniqueIntersections;

      int startIdx = 0, endIdx = 0;
      while (endIdx != int(m_intersections.size()))
      {
        endIdx = get_end_idx(startIdx);
        PointRelation relation = get_unique_intersection(startIdx, endIdx);
        uniqueIntersections.push_back(relation); 
        startIdx = endIdx;   
      }

      return uniqueIntersections;
    }

  private:
    int get_end_idx(int startIdx)
    {
      int recvVertIdx = m_intersections[startIdx].second.id();

      int endIdx = m_intersections.size();
      for (int i=startIdx; i < int(m_intersections.size()); ++i)
        if (m_intersections[i].second.id() != recvVertIdx)
        {
          endIdx = i;
          break;
        }

      return endIdx; 
    }

    bool is_point_inside_element(const utils::Point& ptXi, double tol)
    {
      return ptXi.x >= -tol && ptXi.x <= (1 + tol) &&
             ptXi.y >= -tol && ptXi.y <= (1 + tol);           
    }

    PointRelation get_unique_intersection(int startIdx, int endIdx)
    {
      double tol = 1e-6;
      mesh::MeshEntityPtr vert = m_meshRecv->get_vertices()[m_intersections[startIdx].second.id()];
      mesh::impl::QuadPointFinder finder;
      std::vector<PointRelation> possibleIntersections;

      while (possibleIntersections.size() == 0)
      {
        for (int i=startIdx; i < endIdx; ++i)
        {
          mesh::MeshEntityPtr quad = m_meshSend->get_elements()[m_intersections[i].first.id()];
          utils::Point ptRecv      = vert->get_point_orig(0);

          bool allowFailure = true;
          bool didFail = false;
          utils::Point ptXi   = finder.compute_nearest_point_on_quad(quad, ptRecv, {0.5, 0.5}, allowFailure, didFail);

          if (!didFail && is_point_inside_element(ptXi, tol))
          {
            possibleIntersections.push_back(PointRelation{quad->get_id(), vert->get_id(), ptXi});      
          }
        }

        tol *= 10;

        if (tol > 10)
          throw std::runtime_error("could not find any matches");
      }

      if (possibleIntersections.size() == 0)
        throw std::runtime_error("could not find any possible intersections");


      if (possibleIntersections.size() == 1)
        return possibleIntersections[0];
      else
        return choose_closest_point(possibleIntersections);
    }

    PointRelation choose_closest_point(const std::vector<PointRelation>& possibleIntersections)
    {

      double minDist = std::numeric_limits<double>::max();
      int minIdx = -1;
      for (size_t i=0; i < possibleIntersections.size(); ++i)
      {
        PointRelation relation = possibleIntersections[i];
        mesh::MeshEntityPtr quad = m_meshSend->get_elements()[relation.senderId];
        mesh::MeshEntityPtr vert = m_meshRecv->get_vertices()[relation.destId];

        utils::Point ptSend = mesh::compute_coords_from_xi_3d(quad, relation.senderPtXi);
        utils::Point disp   = ptSend - vert->get_point_orig(0);
        double distSquared  = dot(disp, disp);

        if (distSquared < minDist)
        {
          minDist = distSquared;
          minIdx = i;
        }
      }

      assert(minIdx >= 0);
      return possibleIntersections[minIdx]; 
    }

    const std::vector<std::pair<EntityIdentifier, EntityIdentifier>>& m_intersections;
    std::shared_ptr<mesh::Mesh> m_meshSend;
    std::shared_ptr<mesh::Mesh> m_meshRecv;
};


void interpolate_field(mesh::FieldPtr<double> sendFieldPtr, mesh::FieldPtr<double> recvFieldPtr)
{
  std::shared_ptr<mesh::Mesh> sendMesh = sendFieldPtr->get_mesh();
  std::shared_ptr<mesh::Mesh> recvMesh = recvFieldPtr->get_mesh();

  auto possibleIntersections = compute_possible_intersections(sendMesh, recvMesh);
  check_all_points_are_matched(recvMesh, possibleIntersections);
  UniqueIntersectionFinder finder(possibleIntersections, sendMesh, recvMesh);
  auto uniqueIntersections = finder.find();

  auto& sendField = *sendFieldPtr;
  auto& recvField = *recvFieldPtr;
  for (const PointRelation& relation : uniqueIntersections)
  {
    mesh::MeshEntityPtr quad = sendMesh->get_elements()[relation.senderId];
    mesh::MeshEntityPtr vert = recvMesh->get_vertices()[relation.destId];

    std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
    std::array<double, 4> vertVals;

    mesh::get_downward(quad, 0, verts.data());
    for (int i=0; i < 4; ++i)
      vertVals[i] = sendField(verts[i], 0, 0);

    recvField(vert, 0, 0) = mesh::interpolate(quad->get_type(), vertVals.data(), relation.senderPtXi);
  }
}

template <typename Tfunc>
mesh::FieldPtr<double> set_field(mesh::FieldPtr<double> fieldPtr, Tfunc func)
{
  auto& field = *fieldPtr;
  for (auto& vert : fieldPtr->get_mesh()->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      field(vert, 0, 0) = func(pt);
    }

  return fieldPtr;
}

template <typename Tfunc>
void check_field(mesh::FieldPtr<double> fieldPtr, Tfunc func)
{
  auto& field = *fieldPtr;
  for (auto& vert : field.get_mesh()->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      double valExpected = func(pt);
      //std::cout << "field = " << field(vert, 0, 0) << ", val expected = " << valExpected << ", error = " << std::abs(field(vert, 0, 0) - valExpected) << std::endl;
      if (std::abs(field(vert, 0, 0) - valExpected) > 1e-12)
        throw std::runtime_error("field transfer was not exact");
    }
}

const double one_over_rad3 = 1.0/std::sqrt(3.0);

const std::vector<utils::Point> quadPts = { mesh::convert_xi_coords_from_range(-1, 1, {-one_over_rad3, -one_over_rad3, 0}),
                                            mesh::convert_xi_coords_from_range(-1, 1, { one_over_rad3, -one_over_rad3, 0}),
                                            mesh::convert_xi_coords_from_range(-1, 1, { one_over_rad3,  one_over_rad3, 0}),
                                            mesh::convert_xi_coords_from_range(-1, 1, {-one_over_rad3,  one_over_rad3, 0})};
const std::vector<double> quadWeights = {0.25, 0.25, 0.25, 0.25};


double integrateField(mesh::FieldPtr<double> fieldPtr)
{
  std::vector<double> detJ(quadPts.size());
  mesh::impl::QuadMetrics metrics;  
  double val = 0.0;
  auto& field = *fieldPtr;
  for (auto& quad : fieldPtr->get_mesh()->get_elements())
    if (quad)
    {
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
      std::array<double, 4> vertVals;

      mesh::get_downward(quad, 0, verts.data());
      for (int i=0; i < 4; ++i)
        vertVals[i] = field(verts[i], 0, 0);

      metrics.compute_det_jacobian(quad, quadPts, detJ);

      for (size_t i=0; i < quadPts.size(); ++i)
        val += mesh::interpolate(quad->get_type(), vertVals.data(), quadPts[i]) * quadWeights[i] * detJ[i];
    }

  return val;
}

std::function<double(const utils::Point&)> function_factory(const std::string& functionName)
{
  if (functionName == "constant")
  {
    return [](const utils::Point& /*pt*/) { return 1; };  
  } else if (functionName == "linear")
  {
    return [](const utils::Point& pt) { return pt.x + 2*pt.y + 3*pt.z; };
  } else if (functionName == "quadratic")
  {
    return [](const utils::Point& pt) { return pt.x*pt.x + 2*pt.y*pt.y + 3*pt.z; };
  } else if (functionName == "exponential")
  {
    return [](const utils::Point& pt) { return std::exp(pt.x + pt.y + pt.z ); };
  } else
    throw std::runtime_error("unrecognized function name: " + functionName);
}

void write_output(mesh::FieldPtr<double> field1, mesh::FieldPtr<double> field2, const std::string& functionName)
{
  std::shared_ptr<mesh::Mesh> inputMesh1 = field1->get_mesh();
  std::shared_ptr<mesh::Mesh> inputMesh2 = field2->get_mesh();

  auto func = function_factory(functionName);
  auto field1Exact = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
  auto field2Exact = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
  set_field(field1Exact, func);
  set_field(field2Exact, func);
  
  auto field1Adaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field1, "field");
  auto field1ExactAdaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field1Exact, "fieldExact");
  stk_interface::impl::ExodusWriter writer1(inputMesh1, {field1Adaptor, field1ExactAdaptor});
  writer1.write("mesh1_nonconservative.exo");

  auto field2Adaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field2, "field");
  auto field2ExactAdaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field2Exact, "fieldExact");
  stk_interface::impl::ExodusWriter writer2(inputMesh2, {field2Adaptor, field2ExactAdaptor});
  writer2.write("mesh2_nonconservative.exo");   
}


int main(int argc, char* argv[])
{
  stk::initialize(&argc, &argv);

  {
    std::string defaultFileName1 = "generated:3x3x1|sideset:Z|bbox:0,0,0,1,1,1";
    std::string defaultFileName2 = "generated:4x4x1|sideset:z|bbox:0,0,1,1,1,2";

    std::string mesh1FileName = stk::get_command_line_option(argc, argv, "mesh1", defaultFileName1);
    std::string mesh2FileName = stk::get_command_line_option(argc, argv, "mesh2", defaultFileName2);

    std::string defaultPartName1 = "surface_1";
    std::string defaultPartName2 = "surface_1";
    std::string partName1 = stk::get_command_line_option(argc, argv, "part-name1", defaultPartName1);
    std::string partName2 = stk::get_command_line_option(argc, argv, "part-name2", defaultPartName2);

    std::string defaultFunctionName = "linear";
    std::string functionName = stk::get_command_line_option(argc, argv, "function-name", defaultFunctionName);

    int defaultNumIters = 64;
    int numIters = stk::get_command_line_option(argc, argv, "num-iters", defaultNumIters);

    stk_interface::StkMeshCreator creator1(mesh1FileName);
    std::shared_ptr<mesh::Mesh> inputMesh1 = creator1.create_mesh_from_part(partName1).mesh;

    stk_interface::StkMeshCreator creator2(mesh2FileName);
    std::shared_ptr<mesh::Mesh> inputMesh2 = creator2.create_mesh_from_part(partName2).mesh;


    auto func = function_factory(functionName);
    auto field1 = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
    auto field2 = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);

    set_field(field1, func);

    for (int i=0; i < numIters; ++i)
    {
      std::cout << "\ntransfer iteration " << i << std::endl;
      interpolate_field(field1, field2);
      interpolate_field(field2, field1);

      double sendIntegral = integrateField(field1);
      double recvIntegral = integrateField(field2);
      std::cout << "sendIntegral = " << sendIntegral << ", recvIntegral = " << recvIntegral
                << ", diff = " << std::abs(sendIntegral - recvIntegral) << std::endl;      
    }

    if (functionName == "linear")
      check_field(field2, func);

    write_output(field1, field2, functionName);
  }
  
  stk::finalize();

  return 0;
}
