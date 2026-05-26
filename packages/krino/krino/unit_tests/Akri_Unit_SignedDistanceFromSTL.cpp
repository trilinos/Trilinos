#include <Akri_BoundingBox.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_OutputUtils.hpp>
#include <gtest/gtest.h>
#include <Akri_Unit_CreateFacetedSphere.hpp>

#include <iostream>
#include <fstream>

krino::BoundingBox
compute_point_bbox(std::vector<stk::math::Vector3d> points)
{
  krino::BoundingBox bbox;

  for ( auto && x : points )
    bbox.accommodate( x );

  return bbox;
}

stk::math::Vector3d
random_point_in_bbox(const krino::BoundingBox bbox)
{
  // TODO: importance sampling
  const double randX =  bbox.get_min()[0] + (static_cast<double>(rand()) / RAND_MAX) * (bbox.get_max()[0] - bbox.get_min()[0]);
  const double randY =  bbox.get_min()[1] + (static_cast<double>(rand()) / RAND_MAX) * (bbox.get_max()[1] - bbox.get_min()[1]);
  const double randZ =  bbox.get_min()[2] + (static_cast<double>(rand()) / RAND_MAX) * (bbox.get_max()[2] - bbox.get_min()[2]);
  return stk::math::Vector3d(randX, randY, randZ);
}

std::vector<double>
compute_surface_signed_distance_values(krino::FacetedSurfaceBase & distanceSurface, const std::vector<stk::math::Vector3d> & queryPoints)
{
  const krino::BoundingBox pointBbox = compute_point_bbox(queryPoints);
  distanceSurface.prepare_to_compute(0., pointBbox, 0.);

  std::vector<double> signedDistances;
  for (auto & x : queryPoints)
  {
    const double signedDist = distanceSurface.point_distance(x, 0., 0., true);
    signedDistances.push_back(signedDist);
  }
  return signedDistances;
}

std::vector<stk::math::Vector3d>
compute_random_locations_in_box(const krino::BoundingBox & boundingBox, const unsigned numQueryPoints)
{
  std::vector<stk::math::Vector3d> queryPoints;
  queryPoints.reserve(numQueryPoints);
  for (unsigned i = 0; i < numQueryPoints; ++i)
    queryPoints.push_back(random_point_in_bbox(boundingBox));

  return queryPoints;
}


TEST(ComputeDistanceToSTL, randomQueryPointsBothInsideAndOutsideSphereInSerial_getExpectedSignedDistance)
{
  stk::ParallelMachine comm(MPI_COMM_WORLD);

  const double radius = 0.75;
  const double meshSize = 0.15;
  const std::string filename = "sphere.stl";
  krino::write_stl_for_sphere(filename, radius, meshSize, comm);

  krino::STLSurface distanceSurface(filename);

  krino::BoundingBox boundingBox = distanceSurface.get_bounding_box();
  boundingBox.scale(2.0);

  const std::vector<stk::math::Vector3d> queryPoints = compute_random_locations_in_box(boundingBox, 20);
  const std::vector<double> signedDistances = compute_surface_signed_distance_values(distanceSurface, queryPoints);

  for (size_t i = 0; i < queryPoints.size(); ++i)
  {
    const double goldSignedDist = queryPoints[i].length() - radius;
    EXPECT_NEAR(signedDistances[i], goldSignedDist, 0.02);
  }
}

void compute_training_data(const std::string filename, const unsigned numQueryPoints)
{
  krino::STLSurface distanceSurface(filename);

  krino::BoundingBox boundingBox = distanceSurface.get_bounding_box();
  boundingBox.scale(1.2); // TODO: play with this scaling factor

  const std::vector<stk::math::Vector3d> queryPoints = compute_random_locations_in_box(boundingBox, numQueryPoints);
  const std::vector<double> signedDistances = compute_surface_signed_distance_values(distanceSurface, queryPoints);

  // Write the data to file
  // TODO: use a binary file format
  std::ofstream file("data.csv");
  file << "x,y,z,dist\n";
  for (unsigned i = 0; i < numQueryPoints; ++i)
  {
    file << queryPoints[i][0] << ","
         << queryPoints[i][1] << ","
         << queryPoints[i][2] << ","
         << signedDistances[i] << "\n";
  }
  file.close();
}

int get_num_samples()
{
  const char* p = std::getenv("KRINO_SDF_NUM_SAMPLES");
  return p ? std::stoi(p) : 5000;          // default = 5000 if not set
}

std::string get_stl_filename()
{
  const char* p = std::getenv("KRINO_STL_FILENAME");
  return p ? std::string(p) : "vehicle.stl";   // default = "vehicle.stl" if not set
}

TEST(ComputeDistanceToSTL, computeTrainingData)
{
  stk::ParallelMachine comm(MPI_COMM_WORLD);
  const std::string filename = get_stl_filename();

  if (1 == stk::parallel_machine_size(comm) && krino::does_file_exist(filename))
  {
    int num_data_points = get_num_samples();
    std::cout<<"Computing " << num_data_points << " samples from " << filename << std::endl;
    compute_training_data(filename, num_data_points);
  }
}
