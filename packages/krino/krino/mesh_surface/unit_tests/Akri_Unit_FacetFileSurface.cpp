// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_BoundingBox.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>
#include <Akri_FacetCoordinates.hpp>
#include <Akri_FacetMesh.hpp>
#include <Akri_FacetFileReader.hpp>
#include <Akri_FacetFileSurface.hpp>
#include <Akri_FacetMeshChecks.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_ScopedTestFile.hpp>
#include <Akri_Unit_CreateFacetedSphere.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace krino
{

namespace
{

// A unit square in z=0 split into two triangles sharing the edge (1,0,0)-(0,1,0).
std::string fac_two_triangle_square()
{
  return "4\n"
         "0 0 0 0\n"
         "1 1 0 0\n"
         "2 0 1 0\n"
         "3 1 1 0\n"
         "2\n"
         "0 0 1 2\n"
         "1 1 3 2\n";
}

std::string ply_two_triangle_square()
{
  return "ply\n"
         "format ascii 1.0\n"
         "element vertex 4\n"
         "property float x\n"
         "property float y\n"
         "property float z\n"
         "element face 2\n"
         "property list uchar int vertex_index\n"
         "end_header\n"
         "0 0 0\n"
         "1 0 0\n"
         "0 1 0\n"
         "1 1 0\n"
         "3 0 1 2\n"
         "3 1 3 2\n";
}

// Writes the contents on processor 0; removed (by the ScopedTestFile) when the test ends.
void write_text_file_on_proc0(const ScopedTestFile & file, const std::string & contents)
{
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
  {
    std::ofstream out(file.name());
    out << contents;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

const BoundingBox queryBox(stk::math::Vector3d(-1, -1, -1), stk::math::Vector3d(2, 2, 2));

std::vector<stk::math::Vector3d> query_points()
{
  std::vector<stk::math::Vector3d> points;
  for (double x : {-0.3, 0.2, 0.55, 0.9, 1.4})
    for (double y : {-0.2, 0.35, 0.7, 1.3})
      for (double z : {-0.6, -0.1, 0.25, 0.8})
        points.emplace_back(x, y, z);
  return points;
}

// Points inside a box, off any symmetry plane of the sphere.
std::vector<stk::math::Vector3d> query_points_in(const BoundingBox & box)
{
  std::vector<stk::math::Vector3d> points;
  for (double fx : {0.13, 0.41, 0.77})
    for (double fy : {0.07, 0.38, 0.62, 0.91})
      for (double fz : {0.19, 0.53, 0.86})
        points.emplace_back(box.get_min()[0] + fx * (box.get_max()[0] - box.get_min()[0]),
            box.get_min()[1] + fy * (box.get_max()[1] - box.get_min()[1]),
            box.get_min()[2] + fz * (box.get_max()[2] - box.get_min()[2]));
  return points;
}

std::unique_ptr<MeshFacetedSurface> build_connected_surface(const std::string & format,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale)
{
  auto surface = std::make_unique<FacetFileSurface>(
      stk::EnvData::parallel_comm(), format, filename, sign, scale);
  surface->prepare_to_compute(0., queryBox, 0.);
  return surface;
}

size_t global_sum_of(const size_t local)
{
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

size_t global_num_vertices_of(const MeshFacetedSurface & surface)
{
  return surface.global_num_vertices();
}

size_t global_num_owned_facets(const MeshFacetedSurface & surface)
{
  const size_t local = surface.num_owned_facets();
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

void expect_same_signed_distances(const Surface & expected, const Surface & actual)
{
  for (auto & x : query_points())
    EXPECT_DOUBLE_EQ(expected.point_signed_distance(x), actual.point_signed_distance(x)) << x;
}

// The facets as a surface of facets each with its own vertices (a soup), all of them on every
// processor: the reference that a surface read from a file is compared to.
std::unique_ptr<MeshFacetedSurface> build_soup_surface(const FacetCoordinates & facets)
{
  return std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, facets.coordinates());
}

// Writes the facets as an EXO file of trishell3 elements, the format the EXO reader reads, on
// processor 0 alone so the file is not split across processors.
void write_exo_on_proc0(const FacetCoordinates & facets, const std::string & baseName)
{
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
  {
    FacetMesh facetMesh(MPI_COMM_SELF, 3, {.sharedVertices = false});
    facetMesh.set_owned_facets_from_coordinates(facets.coordinates());
    write_facets(facetMesh, baseName, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace

TEST(FacetFileSurface, supportedFormats)
{
  EXPECT_TRUE(is_facet_file_format("FAC"));
  EXPECT_TRUE(is_facet_file_format("PLY"));
  EXPECT_TRUE(is_facet_file_format("EXO"));
  EXPECT_TRUE(is_facet_file_format("STL"));
  EXPECT_FALSE(is_facet_file_format("OBJ"));
  EXPECT_ANY_THROW(FacetFileSurface(
      stk::EnvData::parallel_comm(), "OBJ", "unused.obj", 1, stk::math::Vector3d(1., 1., 1.)));
}

// The two-triangle square as disconnected facets, scaled, with the orientation flipped for sign -1.
std::unique_ptr<MeshFacetedSurface>
build_square_reference(const int sign, const stk::math::Vector3d & scale)
{
  const std::array<stk::math::Vector3d, 4> pts{
      {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}}};
  const std::array<std::array<unsigned, 3>, 2> facets{{{{0, 1, 2}}, {{1, 3, 2}}}};
  auto scaled = [&](const unsigned i)
  { return stk::math::Vector3d(pts[i][0] * scale[0], pts[i][1] * scale[1], pts[i][2] * scale[2]); };
  FacetCoordinates referenceFacets(3);
  for (auto & f : facets)
  {
    if (sign == -1)
      referenceFacets.emplace_back_3d(scaled(f[1]), scaled(f[0]), scaled(f[2]));
    else
      referenceFacets.emplace_back_3d(scaled(f[0]), scaled(f[1]), scaled(f[2]));
  }
  std::unique_ptr<MeshFacetedSurface> reference = build_soup_surface(referenceFacets);
  reference->prepare_to_compute(0., queryBox, 0.);
  return reference;
}

TEST(FacetFileSurface, facFile_sharesVerticesAndMatchesDisconnectedFacets)
{
  const ScopedTestFile file("unit_test_connected.fac");
  write_text_file_on_proc0(file, fac_two_triangle_square());
  const stk::math::Vector3d scale(2., 1., 3.);

  for (const int sign : {1, -1})
  {
    std::unique_ptr<MeshFacetedSurface> connected =
        build_connected_surface("FAC", file.name(), sign, scale);
    EXPECT_EQ(2u, global_num_owned_facets(*connected));
    EXPECT_EQ(4u, global_num_vertices_of(*connected)); // the shared edge's vertices are stored once

    expect_same_signed_distances(*build_square_reference(sign, scale), *connected);
  }
}

TEST(FacetFileSurface, plyFile_sharesVerticesAndMatchesDisconnectedFacets)
{
  const ScopedTestFile file("unit_test_connected.ply");
  write_text_file_on_proc0(file, ply_two_triangle_square());
  const stk::math::Vector3d scale(2., 1., 3.);

  for (const int sign : {1, -1})
  {
    std::unique_ptr<MeshFacetedSurface> connected =
        build_connected_surface("PLY", file.name(), sign, scale);
    EXPECT_EQ(2u, global_num_owned_facets(*connected));
    EXPECT_EQ(4u, global_num_vertices_of(*connected));

    expect_same_signed_distances(*build_square_reference(sign, scale), *connected);
  }
}

TEST(FacetFileSurface, exoFile_matchesDisconnectedFacets)
{
  const std::string baseName = "unit_test_connected_exo";
  const ScopedTestFile file(create_file_name(baseName, 0));
  FacetCoordinates facets(3);
  facets.emplace_back_3d(stk::math::Vector3d(0., 0., 0.),
      stk::math::Vector3d(1., 0., 0.),
      stk::math::Vector3d(0., 1., 0.));
  facets.emplace_back_3d(stk::math::Vector3d(1., 0., 0.),
      stk::math::Vector3d(1., 1., 0.),
      stk::math::Vector3d(0., 1., 0.));
  write_exo_on_proc0(facets, baseName);

  const stk::math::Vector3d scale(2., 1., 3.);
  for (const int sign : {1, -1})
  {
    std::unique_ptr<MeshFacetedSurface> connected =
        build_connected_surface("EXO", file.name(), sign, scale);
    EXPECT_EQ(2u, global_num_owned_facets(*connected));

    expect_same_signed_distances(*build_square_reference(sign, scale), *connected);
  }
}

// An EXO block larger than one read batch (on more than one processor, see
// compute_facet_file_read_batch_size) has every facet read once.
TEST(FacetFileSurface, exoBlockLargerThanOneBatch_readsEveryFacetOnce)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> connectivity;
  fill_sphere_vertices_and_connectivity(0.75, 0.05, vertices, connectivity);
  ASSERT_GT(connectivity.size(), 2048u);

  const std::string baseName = "unit_test_exo_many_facets";
  const ScopedTestFile file(create_file_name(baseName, 0));
  FacetCoordinates referenceFacets(3);
  for (auto & facet : connectivity)
    referenceFacets.emplace_back_3d(vertices[facet[0]], vertices[facet[1]], vertices[facet[2]]);
  write_exo_on_proc0(referenceFacets, baseName);
  std::unique_ptr<MeshFacetedSurface> reference = build_soup_surface(referenceFacets);

  const BoundingBox sphereBox(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  reference->prepare_to_compute(0., sphereBox, 0.);

  FacetFileSurface surface(
      stk::EnvData::parallel_comm(), "EXO", file.name(), 1, stk::math::Vector3d(1., 1., 1.));
  surface.prepare_to_compute(0., sphereBox, 0.);
  EXPECT_EQ(connectivity.size(), global_num_owned_facets(surface));
  for (auto & x : query_points())
    EXPECT_DOUBLE_EQ(reference->point_signed_distance(x), surface.point_signed_distance(x)) << x;
}

namespace
{

std::string ply_from_vertices_and_connectivity(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & connectivity)
{
  std::ostringstream out;
  out.precision(17);
  out << "ply\nformat ascii 1.0\nelement vertex " << vertices.size()
      << "\nproperty double x\nproperty double y\nproperty double z\n"
      << "element face " << connectivity.size()
      << "\nproperty list uchar int vertex_index\nend_header\n";
  for (auto & v : vertices)
    out << v[0] << " " << v[1] << " " << v[2] << "\n";
  for (auto & f : connectivity)
    out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
  return out.str();
}

} // namespace

// Each processor queries a disjoint slab in x; the facets read from the file are sent to the first
// processor (>= 1, else 0) whose slab they intersect.
TEST(FacetFileSurface, readsAndDistributesFacetsToTheProcessorsThatQueryThem)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> connectivity;
  fill_sphere_vertices_and_connectivity(radius, 0.05, vertices, connectivity);
  ASSERT_GT(connectivity.size(), 2048u); // several read batches on more than one processor

  const ScopedTestFile file("unit_test_distributed_sphere.ply");
  write_text_file_on_proc0(file, ply_from_vertices_and_connectivity(vertices, connectivity));

  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<BoundingBox> slabs;
  for (int p = 0; p < numProcs; ++p)
  {
    const double width = 2. * radius / numProcs;
    slabs.emplace_back(stk::math::Vector3d(-radius + p * width, -1., -1.),
        stk::math::Vector3d(-radius + (p + 1) * width, 1., 1.));
  }
  const BoundingBox & mySlab = slabs[me];

  size_t expectedNumOwned = 0;
  FacetCoordinates referenceFacets(3);
  for (auto & f : connectivity)
  {
    const Facet3d facet(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
    if (find_destination_proc_for_facet(slabs, facet) == me) ++expectedNumOwned;
    referenceFacets.emplace_back_3d(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
  }
  std::unique_ptr<MeshFacetedSurface> reference = build_soup_surface(referenceFacets);
  reference->prepare_to_compute(0., mySlab, 0.);

  FacetFileSurface surface(
      stk::EnvData::parallel_comm(), "PLY", file.name(), 1, stk::math::Vector3d(1., 1., 1.));
  surface.prepare_to_compute(0., mySlab, 0.);

  EXPECT_EQ(expectedNumOwned, surface.num_owned_facets());
  EXPECT_EQ(connectivity.size(), global_num_owned_facets(surface));
  for (auto & x : query_points_in(mySlab))
    EXPECT_NEAR(reference->point_signed_distance(x), surface.point_signed_distance(x), 1.e-12) << x;

  surface.prepare_to_compute(0., mySlab, 0.); // later preparation does not move facets
  EXPECT_EQ(expectedNumOwned, surface.num_owned_facets());
}

namespace
{

// The sphere with every vertex in z > 0 listed twice, once for the facets with centroid x > 0 and
// once for the others, so the file has a seam: the facets on its two sides do not share those
// vertices.
void split_sphere_vertices_along_seam(std::vector<stk::math::Vector3d> & vertices,
    std::vector<std::array<unsigned, 3>> & connectivity)
{
  const unsigned numOriginal = vertices.size();
  std::vector<unsigned> copyOf(numOriginal, 0);
  for (unsigned v = 0; v < numOriginal; ++v)
  {
    if (vertices[v][2] > 0.)
    {
      copyOf[v] = vertices.size();
      vertices.push_back(vertices[v]);
    }
  }
  for (auto & f : connectivity)
    if (vertices[f[0]][0] + vertices[f[1]][0] + vertices[f[2]][0] > 0.)
      for (auto & v : f)
        if (v < numOriginal && copyOf[v] != 0) v = copyOf[v];
}

} // namespace

TEST(FacetFileSurface, fileWithSeam_mergedVerticesConnectFacetsAcrossIt)
{
  const double radius = 0.75;
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> connectivity;
  fill_sphere_vertices_and_connectivity(radius, 0.15, vertices, connectivity);
  const size_t numSphereVertices = vertices.size();
  FacetCoordinates referenceFacets(3);
  for (auto & f : connectivity)
    referenceFacets.emplace_back_3d(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
  std::unique_ptr<MeshFacetedSurface> reference = build_soup_surface(referenceFacets);

  split_sphere_vertices_along_seam(vertices, connectivity);
  std::vector<unsigned> usedVertices;
  for (auto & f : connectivity)
    usedVertices.insert(usedVertices.end(), f.begin(), f.end());
  stk::util::sort_and_unique(usedVertices);
  const size_t numUsedVerticesWithSeam = usedVertices.size();
  ASSERT_GT(numUsedVerticesWithSeam, numSphereVertices);
  const ScopedTestFile file("unit_test_sphere_with_seam.ply");
  write_text_file_on_proc0(file, ply_from_vertices_and_connectivity(vertices, connectivity));

  const BoundingBox box(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  reference->prepare_to_compute(0., box, 0.);

  // Merged, the vertices along the seam are shared, as on the sphere without it.  Unmerged, they
  // are not, and the sign comes from the facets near the closest point instead; both get the
  // sphere's distances.
  for (const bool merge : {true, false})
  {
    FacetFileSurface surface(stk::EnvData::parallel_comm(),
        "PLY",
        file.name(),
        1,
        stk::math::Vector3d(1., 1., 1.),
        merge);
    surface.prepare_to_compute(0., box, 0.);
    EXPECT_EQ(merge ? numSphereVertices : numUsedVerticesWithSeam, global_num_vertices_of(surface));
    EXPECT_EQ(connectivity.size(), global_num_owned_facets(surface));
    for (auto & x : query_points_in(box))
      EXPECT_DOUBLE_EQ(reference->point_signed_distance(x), surface.point_signed_distance(x))
          << x << " merge " << merge;
  }
}

namespace
{

std::string
ascii_stl_from_vertices_and_connectivity(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & connectivity)
{
  std::ostringstream out;
  out.precision(17);
  out << "solid sphere\n";
  for (auto & f : connectivity)
  {
    const Facet3d facet(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
    const stk::math::Vector3d normal = facet.facet_normal();
    out << "  facet normal " << normal[0] << " " << normal[1] << " " << normal[2]
        << "\n    outer loop\n";
    for (unsigned v : f)
      out << "      vertex " << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2]
          << "\n";
    out << "    endloop\n  endfacet\n";
  }
  out << "endsolid sphere\n";
  return out.str();
}

// The sphere's facets as disconnected facets, with the orientation flipped for sign -1.
std::unique_ptr<MeshFacetedSurface>
build_sphere_reference(const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & connectivity,
    const int sign)
{
  FacetCoordinates referenceFacets(3);
  for (auto & f : connectivity)
  {
    if (sign == -1)
      referenceFacets.emplace_back_3d(vertices[f[0]], vertices[f[2]], vertices[f[1]]);
    else
      referenceFacets.emplace_back_3d(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
  }
  return build_soup_surface(referenceFacets);
}

// Moves each coordinate to a multiple of 2^-20, which for coordinates of magnitude below 1 has at
// most 21 significant bits, so it is exactly representable as a float.  The STL files then hold the
// vertices exactly, and the test needs no float rounding of its own (gcc 12 at -O2 drops the
// rounding of double(float(x)) when the components are stored together).
void make_coordinates_exact_in_float(std::vector<stk::math::Vector3d> & vertices)
{
  for (auto & v : vertices)
    for (unsigned d = 0; d < 3; ++d)
      v[d] = std::ldexp(std::round(std::ldexp(v[d], 20)), -20);
}

void expect_connected_stl_surface_matches_sphere(const std::string & filename,
    const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & connectivity)
{
  const BoundingBox box(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  for (const int sign : {1, -1})
  {
    FacetFileSurface connected(
        stk::EnvData::parallel_comm(), "STL", filename, sign, stk::math::Vector3d(1., 1., 1.));
    std::unique_ptr<MeshFacetedSurface> reference =
        build_sphere_reference(vertices, connectivity, sign);

    BoundingBox connectedBbox, referenceBbox;
    connected.insert_into(connectedBbox);
    for (auto & v : vertices)
      referenceBbox.accommodate(v);
    EXPECT_EQ(referenceBbox.get_min(), connectedBbox.get_min());
    EXPECT_EQ(referenceBbox.get_max(), connectedBbox.get_max());

    connected.prepare_to_compute(0., box, 0.);
    reference->prepare_to_compute(0., box, 0.);
    EXPECT_EQ(connectivity.size(), global_num_owned_facets(connected));
    EXPECT_EQ(vertices.size(), global_num_vertices_of(connected)); // each facet's copies merged
    for (auto & x : query_points_in(box))
      EXPECT_NEAR(reference->point_signed_distance(x), connected.point_signed_distance(x), 1.e-12)
          << x << " sign " << sign;
  }
}

} // namespace

// An STL file lists each facet's vertices; merging the coincident ones connects the facets.
TEST(FacetFileSurface, binaryAndAsciiStl_mergedIntoConnectedSphere)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> connectivity;
  fill_sphere_vertices_and_connectivity(0.75, 0.15, vertices, connectivity);
  make_coordinates_exact_in_float(vertices);

  const ScopedTestFile binaryFile("unit_test_connected_sphere_binary.stl");
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
    write_stl(binaryFile.name(), vertices, connectivity);
  MPI_Barrier(MPI_COMM_WORLD);
  expect_connected_stl_surface_matches_sphere(binaryFile.name(), vertices, connectivity);

  const ScopedTestFile asciiFile("unit_test_connected_sphere_ascii.stl");
  write_text_file_on_proc0(
      asciiFile, ascii_stl_from_vertices_and_connectivity(vertices, connectivity));
  expect_connected_stl_surface_matches_sphere(asciiFile.name(), vertices, connectivity);
}

// Facets with a non-finite normal or no area are skipped.
TEST(FacetFileSurface, stlWithInvalidFacets_skipsThem)
{
  const ScopedTestFile file("unit_test_connected_invalid_facets.stl");
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
  {
    std::ofstream out(file.name(), std::ios::binary);
    const char header[80] = "binary stl with invalid facets";
    out.write(header, 80);
    const unsigned numFacets = 3;
    out.write(reinterpret_cast<const char *>(&numFacets), sizeof(numFacets));
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const std::array<std::array<float, 12>, 3> facets{
        {{{0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0}},      // valid
            {{nan, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0}}, // non-finite normal
            {{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0}}}}; // degenerate
    for (auto & facet : facets)
    {
      out.write(reinterpret_cast<const char *>(facet.data()), sizeof(float) * facet.size());
      const char attribute[2] = {0, 0};
      out.write(attribute, 2);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  FacetFileSurface surface(
      stk::EnvData::parallel_comm(), "STL", file.name(), 1, stk::math::Vector3d(1., 1., 1.));
  surface.prepare_to_compute(0., queryBox, 0.);
  EXPECT_EQ(1u, global_num_owned_facets(surface));
  EXPECT_EQ(3u, global_num_vertices_of(surface));
}

// With remove_coincident_facets, a facet listed twice (as some files do) counts once.
TEST(FacetFileSurface, duplicatedFacetWithRemoveCoincidentFacets_sameAsWithout)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned, 3>> connectivity;
  fill_sphere_vertices_and_connectivity(0.75, 0.15, vertices, connectivity);
  // The duplicate has its own copies of its vertices, so it shares vertex ids with the original
  // only when coincident vertices are merged; otherwise it is found by location.
  std::vector<stk::math::Vector3d> verticesWithCopies = vertices;
  std::vector<std::array<unsigned, 3>> withDuplicate = connectivity;
  std::array<unsigned, 3> duplicate;
  for (unsigned n = 0; n < 3; ++n)
  {
    duplicate[n] = verticesWithCopies.size();
    verticesWithCopies.push_back(vertices[connectivity[7][n]]);
  }
  withDuplicate.push_back(duplicate);

  const ScopedTestFile file("unit_test_connected_duplicated_facet.ply");
  write_text_file_on_proc0(
      file, ply_from_vertices_and_connectivity(verticesWithCopies, withDuplicate));

  const BoundingBox box(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
  for (const bool merge : {true, false})
  {
    FacetFileSurface surface(stk::EnvData::parallel_comm(),
        "PLY",
        file.name(),
        1,
        stk::math::Vector3d(1., 1., 1.),
        merge);
    surface.set_remove_coincident_facets(true);
    surface.set_check_watertight(true);
    surface.prepare_to_compute(0., box, 0.);
    EXPECT_EQ(connectivity.size(), global_num_owned_facets(surface)) << "merge " << merge;
  }
  FacetFileSurface surface(
      stk::EnvData::parallel_comm(), "PLY", file.name(), 1, stk::math::Vector3d(1., 1., 1.));
  surface.set_remove_coincident_facets(true);
  surface.prepare_to_compute(0., box, 0.);
  EXPECT_EQ(
      0u, global_sum_of(find_non_watertight_edges_by_vertex_ids(surface.facet_mesh(), {}).size()));
}

// The STL reader gives each distinct point once when merging, as it reads, else 3 per facet.
TEST(STLFileReader, twoTrianglesSharingAnEdge_mergedAsRead)
{
  const ScopedTestFile file("unit_test_stl_reader_square.stl");
  // -0 and 0 are the same coordinate.
  write_text_file_on_proc0(file,
      "solid s\n"
      "facet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 "
      "0\nendloop\nendfacet\n"
      "facet normal 0 0 1\nouter loop\nvertex 1 -0 0\nvertex 1 1 0\nvertex 0 1 "
      "0\nendloop\nendfacet\n"
      "endsolid s\n");
  if (0 != stk::parallel_machine_rank(MPI_COMM_WORLD)) return;

  for (const bool merge : {true, false})
  {
    auto reader =
        create_facet_file_reader("STL", file.name(), 1, stk::math::Vector3d(1., 1., 1.), merge);
    std::vector<stk::math::Vector3d> coords;
    std::vector<stk::mesh::EntityId> ids;
    reader->read_vertices(coords, ids);
    ASSERT_EQ(2u, reader->num_facets());
    EXPECT_EQ(merge ? 4u : 6u, coords.size());
    EXPECT_EQ(coords.size(), ids.size());
    const std::array<size_t, 3> facet0 = reader->read_facet();
    const std::array<size_t, 3> facet1 = reader->read_facet();
    if (merge)
    {
      EXPECT_EQ(facet0[1], facet1[0]); // (1,0,0) and (1,-0,0)
      EXPECT_EQ(facet0[2], facet1[2]); // (0,1,0)
    }
    else
    {
      EXPECT_EQ((std::array<size_t, 3>{{3, 4, 5}}), facet1);
    }
    EXPECT_EQ(stk::math::Vector3d(1., 1., 0.), coords[facet1[1]]);
  }
}

TEST(FacetFileSurface, missingFile_throwsOnAllProcessors)
{
  FacetFileSurface surface(stk::EnvData::parallel_comm(),
      "PLY",
      "does_not_exist.ply",
      1,
      stk::math::Vector3d(1., 1., 1.));
  EXPECT_ANY_THROW(surface.prepare_to_compute(0., queryBox, 0.));
}

} // namespace krino
