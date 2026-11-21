// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_Transformation.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <Akri_Unit_SurfaceDecomposition.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

// A 10x10x10 cube centered at the origin, rotated about x, then y, then z, and decomposed with
// snapping (sharp features cut before the surface, as by default) on a BCC tet4 mesh of size 2.5,
// 4 cells across the cube, as ngs meshes its Brick10x10x10AtOrigin (UnitTestPGeomMakeMovies.cpp).
// The domain is the same for every rotation, and holds the cube at any of them; unrotated, the
// cube's faces are on mesh planes.  The sweeps also run the cuboid of RotatedCuboidDecomposition.
// RotatedCubeStudy.* (sweeps, movie) are long runs, left out of the default runs by their filters, that
// write a table of the measures of each rotation.

namespace krino
{

namespace
{

using DegreesAboutXYZ = std::array<double, 3>;

// The rotation by degrees[0] about x, then degrees[1] about y, then degrees[2] about z (fixed axes),
// as ngs's RotationAngles, as a rotation vector (see Quaternion).
stk::math::Vector3d rotation_vector_about_x_then_y_then_z(const DegreesAboutXYZ & degrees)
{
  Quaternion rotation;
  for (int axis = 0; axis < 3; ++axis)
  {
    stk::math::Vector3d aboutAxisVector = stk::math::Vector3d::ZERO;
    aboutAxisVector[axis] = degrees[axis] * M_PI / 180.;
    Quaternion aboutAxis;
    aboutAxis.set_from_rotation_vector(aboutAxisVector);
    rotation = aboutAxis.after(rotation);
  }
  return rotation.rotation_vector();
}

// A box (cube or cuboid) of the study, rotated about its center, and the mesh it is decomposed on: a
// BCC tet4 mesh of the given size over a cube, centered at the origin, that holds the box at any rotation.
struct RotatedBoxStudyGeometry
{
  stk::math::Vector3d dimensions;
  stk::math::Vector3d center;
  double meshSize;
  double domainHalfWidth;

  double volume() const { return dimensions[0] * dimensions[1] * dimensions[2]; }
  double area() const
  {
    return 2. * (dimensions[0] * dimensions[1] + dimensions[0] * dimensions[2] + dimensions[1] * dimensions[2]);
  }
  BoundingBoxMesh::BoundingBoxType domain() const
  {
    return BoundingBoxMesh::BoundingBoxType(stk::math::Vector3d(-domainHalfWidth, -domainHalfWidth, -domainHalfWidth),
        stk::math::Vector3d(domainHalfWidth, domainHalfWidth, domainHalfWidth));
  }
  Cuboid box(const DegreesAboutXYZ & degrees) const
  {
    return Cuboid(center, dimensions, rotation_vector_about_x_then_y_then_z(degrees));
  }
};

// The cube: 5 cells from the center to each side of the domain, past the cube's circumscribed sphere
// (radius 8.66) by more than a cell, with the unrotated cube's faces, at +-5, on mesh planes.
const RotatedBoxStudyGeometry rotatedCubeGeometry{stk::math::Vector3d(10., 10., 10.), stk::math::Vector3d::ZERO, 2.5, 12.5};

// The cuboid of RotatedCuboidDecomposition (Akri_Unit_RotatedCuboidDecomposition.cpp): unequal sides and an
// offset center, so no face is on a mesh plane by symmetry; its circumscribed sphere (radius 0.73 about a
// center 0.06 from the origin) is more than a cell inside the domain.
const RotatedBoxStudyGeometry rotatedCuboidStudyGeometry{
    stk::math::Vector3d(1.0, 0.7, 0.8), stk::math::Vector3d(0.03, -0.02, 0.05), 0.2, 1.0};

DecompositionMeasures decompose_rotated_cube(const DegreesAboutXYZ & degrees)
{
  SurfaceDecomposition decomposition(rotatedCubeGeometry.domain(), rotatedCubeGeometry.meshSize);
  return decomposition.decompose_and_measure(rotatedCubeGeometry.box(degrees));
}

// The largest over the processors of the seconds since start.
double max_seconds_since(const std::chrono::steady_clock::time_point start)
{
  const double local = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  double global = local;
  stk::all_reduce_max(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

// The measures of each rotation, one row each, written by processor 0 as comma-separated values.
class RotatedCubeTable
{
public:
  RotatedCubeTable(const std::string & fileName, const RotatedBoxStudyGeometry & geometry)
      : myExactVolume(geometry.volume()), myExactArea(geometry.area())
  {
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) != 0) return;
    myFile.open(fileName);
    myFile << "x_degrees,y_degrees,z_degrees,elements,relative_volume_error,relative_area_error,"
              "worst_mean_ratio,worst_scaled_jacobian";
    for (const double lowMeanRatio : lowMeanRatios)
      myFile << ",elements_below_mean_ratio_" << lowMeanRatio;
    myFile << ",decomposition_seconds,smoothing_seconds\n";
  }

  void add(const DegreesAboutXYZ & degrees,
      const DecompositionMeasures & measures,
      const double decompositionSeconds,
      const double smoothingSeconds)
  {
    const double volumeError = std::abs(measures.negativeVolume - myExactVolume) / myExactVolume;
    const double areaError = std::abs(measures.interfaceArea - myExactArea) / myExactArea;
    myMaxVolumeError = std::max(myMaxVolumeError, volumeError);
    myMaxAreaError = std::max(myMaxAreaError, areaError);
    if (measures.worstQuality < myWorstQuality)
    {
      myWorstQuality = measures.worstQuality;
      myWorstQualityDegrees = degrees;
    }
    if (!myFile.is_open()) return;
    myFile.precision(10);
    myFile << degrees[0] << "," << degrees[1] << "," << degrees[2] << "," << measures.numElements << ","
           << volumeError << "," << areaError << "," << measures.worstQuality << ","
           << measures.worstScaledJacobian;
    for (const size_t count : measures.numElementsBelowLowMeanRatios)
      myFile << "," << count;
    myFile << "," << decompositionSeconds << "," << smoothingSeconds << std::endl;
  }

  void print_summary(const std::string & name) const
  {
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) != 0) return;
    std::cout << name << ": max relative volume error " << myMaxVolumeError
              << ", max relative area error " << myMaxAreaError << ", worst mean ratio "
              << myWorstQuality << " at (" << myWorstQualityDegrees[0] << ", "
              << myWorstQualityDegrees[1] << ", " << myWorstQualityDegrees[2] << ") degrees"
              << std::endl;
  }

private:
  double myExactVolume;
  double myExactArea;
  std::ofstream myFile;
  double myMaxVolumeError{0.};
  double myMaxAreaError{0.};
  double myWorstQuality{1.};
  DegreesAboutXYZ myWorstQualityDegrees{{0., 0., 0.}};
};

} // namespace

TEST(RotatedCubeDecomposition, rotationAboutXThenYThenZ_isAboutTheFixedAxesInThatOrder)
{
  Quaternion rotation;
  rotation.set_from_rotation_vector(rotation_vector_about_x_then_y_then_z({{90., 90., 0.}}));
  // 90 about x leaves x and turns y to z; then 90 about y turns x to -z and z to x.
  expect_eq(stk::math::Vector3d(0., 0., -1.), rotation.rotate_3d_vector(stk::math::Vector3d(1., 0., 0.)));
  expect_eq(stk::math::Vector3d(1., 0., 0.), rotation.rotate_3d_vector(stk::math::Vector3d(0., 1., 0.)));
}

// Unrotated (faces on mesh planes), a general rotation, and the rotation with the worst quality of
// RotatedCubeStudy.sweep: with its corners and edges captured, the negative phase is the cube.
TEST(RotatedCubeDecomposition, snapping_volumeAndAreaOfTheCubeAtSeveralRotations)
{
  for (const DegreesAboutXYZ & degrees :
      {DegreesAboutXYZ{{0., 0., 0.}}, DegreesAboutXYZ{{13., 26., 39.}}, DegreesAboutXYZ{{39., 13., 0.}}})
  {
    const DecompositionMeasures measures = decompose_rotated_cube(degrees);
    EXPECT_NEAR(rotatedCubeGeometry.volume(), measures.negativeVolume, 1.e-10 * rotatedCubeGeometry.volume());
    EXPECT_NEAR(rotatedCubeGeometry.area(), measures.interfaceArea, 1.e-10 * rotatedCubeGeometry.area());
    EXPECT_LT(0., measures.worstQuality);
  }
}

// All rotations about x, y and z of 0, 13, 26 and 39 degrees (64), as ngs's
// tet_testRotatedByIncrementsOf7DegreesAllAxes (UnitTestPGeomMPIThreaded.cpp), of the box of the geometry,
// with the given edge tolerance (CDFEM_Support's default if none), then, if requested, smoothing (the
// nodes off every boundary, or also those on the topological features within them).  The table is <baseName>.csv, with the seconds of decomposing and smoothing
// (the largest over the processors), and each rotation's mesh is written as one file of the series
// <baseName>.e, in the order of the table (time = row).  Writing on more than one processor needs
// parallel IO.
enum class StudySmoothing
{
  NONE,
  KEEPING_BOUNDARY_NODES_FIXED,
  RESPECTING_MESH_BASED_GEOMETRY,
  RESPECTING_MESH_BASED_GEOMETRY_AND_FEATURE_ANGLE
};

// The feature angle of the study: the boxes' 90 degree edges and corners are sharp, their flat faces not.
const double studyCosFeatureAngle = std::cos(135. * M_PI / 180.);

static void sweep_rotated_box(const std::string & baseName,
    const RotatedBoxStudyGeometry & geometry,
    const std::optional<double> edgeTol,
    const StudySmoothing smoothing)
{
  RotatedCubeTable table(baseName + ".csv", geometry);
  int fileIndex = 0;
  for (int x = 0; x < 45; x += 13)
    for (int y = 0; y < 45; y += 13)
      for (int z = 0; z < 45; z += 13)
      {
        const DegreesAboutXYZ degrees{{double(x), double(y), double(z)}};
        SurfaceDecomposition decomposition(geometry.domain(),
            geometry.meshSize,
            SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE,
            std::nullopt,
            stk::topology::TETRAHEDRON_4,
            edgeTol);
        const auto decompositionStart = std::chrono::steady_clock::now();
        DecompositionMeasures measures = decomposition.decompose_and_measure(geometry.box(degrees));
        const double decompositionSeconds = max_seconds_since(decompositionStart);
        double smoothingSeconds = 0.;
        if (smoothing != StudySmoothing::NONE)
        {
          const auto smoothingStart = std::chrono::steady_clock::now();
          if (smoothing == StudySmoothing::KEEPING_BOUNDARY_NODES_FIXED)
            measures = decomposition.smooth_keeping_boundary_nodes_fixed_and_measure();
          else if (smoothing == StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY)
            measures = decomposition.smooth_respecting_mesh_based_geometry_and_measure();
          else
            measures = decomposition.smooth_respecting_mesh_based_geometry_and_measure(studyCosFeatureAngle);
          smoothingSeconds = max_seconds_since(smoothingStart);
        }
        table.add(degrees, measures, decompositionSeconds, smoothingSeconds);
        decomposition.write_active_mesh(baseName + ".e", fileIndex, fileIndex);
        ++fileIndex;
      }
  table.print_summary(baseName);
}

TEST(RotatedCubeStudy, sweep)
{
  sweep_rotated_box("rotated_cube_sweep", rotatedCubeGeometry, std::nullopt, StudySmoothing::NONE);
}

TEST(RotatedCubeStudy, sweepWithEdgeTolerance0_1)
{
  sweep_rotated_box("rotated_cube_sweep_edge_tol_0.1", rotatedCubeGeometry, 0.1, StudySmoothing::NONE);
}

TEST(RotatedCubeStudy, sweepWithSmoothing)
{
  sweep_rotated_box("rotated_cube_sweep_smoothed", rotatedCubeGeometry, std::nullopt, StudySmoothing::KEEPING_BOUNDARY_NODES_FIXED);
}

TEST(RotatedCubeStudy, sweepWithEdgeTolerance0_1AndSmoothing)
{
  sweep_rotated_box("rotated_cube_sweep_edge_tol_0.1_smoothed", rotatedCubeGeometry, 0.1, StudySmoothing::KEEPING_BOUNDARY_NODES_FIXED);
}

TEST(RotatedCubeStudy, sweepWithSmoothingRespectingMeshBasedGeometry)
{
  sweep_rotated_box("rotated_cube_sweep_smoothed_respecting", rotatedCubeGeometry, std::nullopt, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY);
}

TEST(RotatedCubeStudy, sweepWithEdgeTolerance0_1AndSmoothingRespectingMeshBasedGeometry)
{
  sweep_rotated_box("rotated_cube_sweep_edge_tol_0.1_smoothed_respecting", rotatedCubeGeometry, 0.1, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY);
}

TEST(RotatedCubeStudy, cuboidSweepWithSmoothingRespectingMeshBasedGeometry)
{
  sweep_rotated_box("rotated_cuboid_sweep_smoothed_respecting", rotatedCuboidStudyGeometry, std::nullopt, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY);
}

TEST(RotatedCubeStudy, cuboidSweepWithEdgeTolerance0_1AndSmoothingRespectingMeshBasedGeometry)
{
  sweep_rotated_box("rotated_cuboid_sweep_edge_tol_0.1_smoothed_respecting", rotatedCuboidStudyGeometry, 0.1, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY);
}

TEST(RotatedCubeStudy, sweepWithSmoothingRespectingMeshBasedGeometryAndFeatureAngle)
{
  sweep_rotated_box("rotated_cube_sweep_smoothed_respecting_angle", rotatedCubeGeometry, std::nullopt, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY_AND_FEATURE_ANGLE);
}

TEST(RotatedCubeStudy, sweepWithEdgeTolerance0_1AndSmoothingRespectingMeshBasedGeometryAndFeatureAngle)
{
  sweep_rotated_box("rotated_cube_sweep_edge_tol_0.1_smoothed_respecting_angle", rotatedCubeGeometry, 0.1, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY_AND_FEATURE_ANGLE);
}

TEST(RotatedCubeStudy, cuboidSweepWithSmoothingRespectingMeshBasedGeometryAndFeatureAngle)
{
  sweep_rotated_box("rotated_cuboid_sweep_smoothed_respecting_angle", rotatedCuboidStudyGeometry, std::nullopt, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY_AND_FEATURE_ANGLE);
}

TEST(RotatedCubeStudy, cuboidSweepWithEdgeTolerance0_1AndSmoothingRespectingMeshBasedGeometryAndFeatureAngle)
{
  sweep_rotated_box("rotated_cuboid_sweep_edge_tol_0.1_smoothed_respecting_angle", rotatedCuboidStudyGeometry, 0.1, StudySmoothing::RESPECTING_MESH_BASED_GEOMETRY_AND_FEATURE_ANGLE);
}

TEST(RotatedCubeStudy, cuboidSweep)
{
  sweep_rotated_box("rotated_cuboid_sweep", rotatedCuboidStudyGeometry, std::nullopt, StudySmoothing::NONE);
}

TEST(RotatedCubeStudy, cuboidSweepWithEdgeTolerance0_1)
{
  sweep_rotated_box("rotated_cuboid_sweep_edge_tol_0.1", rotatedCuboidStudyGeometry, 0.1, StudySmoothing::NONE);
}

TEST(RotatedCubeStudy, cuboidSweepWithSmoothing)
{
  sweep_rotated_box("rotated_cuboid_sweep_smoothed", rotatedCuboidStudyGeometry, std::nullopt, StudySmoothing::KEEPING_BOUNDARY_NODES_FIXED);
}

TEST(RotatedCubeStudy, cuboidSweepWithEdgeTolerance0_1AndSmoothing)
{
  sweep_rotated_box("rotated_cuboid_sweep_edge_tol_0.1_smoothed", rotatedCuboidStudyGeometry, 0.1, StudySmoothing::KEEPING_BOUNDARY_NODES_FIXED);
}

// The cube turning 2, 1 and 3 degrees about x, y and z each frame until each has turned a full turn
// (360 frames), as ngs's testRotatedAllAxesMakeMovie, writing each frame's mesh as one file of the
// series rotated_cube_movie.e (time = frame) as well as the table.  Writing on more than one
// processor needs parallel IO.
TEST(RotatedCubeStudy, movie)
{
  RotatedCubeTable table("rotated_cube_movie.csv", rotatedCubeGeometry);
  const DegreesAboutXYZ degreesPerFrame{{2., 1., 3.}};
  DegreesAboutXYZ degrees{{0., 0., 0.}};
  for (int frame = 0; degrees[0] < 360. || degrees[1] < 360. || degrees[2] < 360.; ++frame)
  {
    SurfaceDecomposition decomposition(rotatedCubeGeometry.domain(), rotatedCubeGeometry.meshSize);
    const auto decompositionStart = std::chrono::steady_clock::now();
    const DecompositionMeasures measures = decomposition.decompose_and_measure(rotatedCubeGeometry.box(degrees));
    table.add(degrees, measures, max_seconds_since(decompositionStart), 0.);
    decomposition.write_active_mesh("rotated_cube_movie.e", frame, frame);
    for (int axis = 0; axis < 3; ++axis)
      degrees[axis] += degreesPerFrame[axis];
  }
  table.print_summary("RotatedCubeStudy.movie");
}

} // namespace krino
