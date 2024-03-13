// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshSurface.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Facet.hpp>
#include <Akri_Transformation.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Ioss_SubSystem.h>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

namespace krino{

void Parallel_Facet_File_Reader::read(const std::string & read_description, const std::function<void(void)> & read_function)
{
  bool ok_locally = true;

  if (0 == stk::EnvData::parallel_rank() )
  {
    my_input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try
    {
      read_function();
    }
    catch (std::ifstream::failure & e)
    {
      ok_locally = false;
    }
  }
  const bool ok_globally = stk::is_true_on_all_procs(stk::EnvData::parallel_comm(), ok_locally);
  if (!ok_globally)
  {
    ThrowRuntimeError("Error " << read_description << " for file " << my_filename);
  }
}

void Parallel_Facet_File_Reader::open_file()
{
  read("opening file", [this](){my_input.open(my_filename.c_str());});
}

void Parallel_Facet_File_Reader::get_batch_size(const int local_num_facets, int & batch_size, int & num_batches)
{
  int num_facets = 0;
  stk::all_reduce_sum( stk::EnvData::parallel_comm(), &local_num_facets, &num_facets, 1 );

  const int min_batch_size = std::min(2048, num_facets);
  batch_size = std::max(min_batch_size, 1+num_facets/stk::EnvData::parallel_size());
  num_batches = 1+num_facets/batch_size;

  if(krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "Reading " << num_facets << " facets, using " << num_batches << " batches of " << batch_size << "." << stk::diag::dendl;
}

Faceted_Surface_From_File::Faceted_Surface_From_File(const std::string & surface_name, const stk::diag::Timer &parent_timer)
: Faceted_Surface<Facet3d>(),
  my_timer("Facet File Reader", parent_timer),
  my_built_local_facets(false)
{
}

void Faceted_Surface_From_File::build_local_facets(const BoundingBox & proc_bbox)
{
  stk::diag::TimeBlock timer_(my_timer);

  if(my_built_local_facets) return;

  std::vector<BoundingBox> proc_bboxes;
  BoundingBox::gather_bboxes( proc_bbox, proc_bboxes );

  read_file(proc_bboxes);
  my_built_local_facets = true;
}

void Faceted_Surface_From_File::insert_into(BoundingBox & bbox) const
{
  // One could reasonably expect that inserting the surface into the bbox is a const method for the surface,
  // but the reading methods used to get the bounding box are non-const.  Ugh.
  const BoundingBox facetBbox = const_cast<Faceted_Surface_From_File *>(this)->get_bounding_box();
  bbox.accommodate(facetBbox);
}

void Parallel_Facet_File_Reader::close_file()
{
  read("closing file", [this](){my_input.close();});
}

FACSurface::FACSurface(const std::string & surface_name,
    const stk::diag::Timer &parent_timer,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale)
 : Faceted_Surface_From_File(surface_name, parent_timer),
   my_reader(filename),
   my_dist_sign(sign),
   my_scale(scale)
{
  // Make sure the file is openable, but close it for now to avoid ulimit restrictions.
  my_reader.open_file();
  my_reader.close_file();
}

BoundingBox
FACSurface::get_bounding_box()
{
  BoundingBox bbox;

  my_reader.open_file();

  std::vector<stk::math::Vector3d> points;
  my_reader.read("reading points", [this, &points](){read_points(points);});

  for (auto && point : points)
  {
    bbox.accommodate(point);
  }

  my_reader.close_file();

  bbox.global_reduce();
  return bbox;
}

void FACSurface::read_file(const std::vector<BoundingBox> & proc_bboxes)
{
  my_reader.open_file();

  std::ifstream & input = my_reader.input();
  std::vector<stk::math::Vector3d> points;
  my_reader.read("reading points", [this, &points](){read_points(points);});

  int num_facets = 0;
  my_reader.read("reading number of facets", [&input, &num_facets](){input >> num_facets; STK_ThrowRequire(num_facets > 0);});

  int batch_size = 0;
  int num_batches = 0;
  my_reader.get_batch_size(num_facets, batch_size, num_batches);

  for (int batch = 0; batch < num_batches; ++batch)
  {
    const int current_batch_size = std::min(batch_size, num_facets-batch*batch_size);
    my_reader.read("reading facets", [this, &points, current_batch_size, num_facets](){read_facets(current_batch_size, num_facets, points);});
    parallel_distribute_facets(current_batch_size, proc_bboxes);
  }

  my_reader.close_file();
}

void FACSurface::read_points(std::vector<stk::math::Vector3d> & points)
{
  std::ifstream & input = my_reader.input();
  int num_points;
  input >> num_points;
  STK_ThrowRequire(num_points > 0);
  points.reserve(num_points);
  for ( int i = 0; i < num_points; i++ )
  {
    int id;
    double X, Y, Z;
    input >> id;
    STK_ThrowRequire(id >= 0 && id < num_points);
    input >> X >> Y >> Z;
    points.emplace_back(X*my_scale[0], Y*my_scale[1], Z*my_scale[2]);
  }
}

void FACSurface::read_facets(const int batch_size, const int num_facets, const std::vector<stk::math::Vector3d> & points)
{
  std::ifstream & input = my_reader.input();
  const int num_points = points.size();
  for ( int i = 0; i < batch_size; ++i )
  {
    int id;
    int p1, p2, p3;
    input >> id;
    STK_ThrowRequire(id >= 0 && id < num_facets);
    input >> p1 >> p2 >> p3;

    STK_ThrowRequire(p1 >= 0 && p1 < num_points);
    STK_ThrowRequire(p2 >= 0 && p2 < num_points);
    STK_ThrowRequire(p3 >= 0 && p3 < num_points);

    if (my_dist_sign == -1)
    {
      // permute to flip normal direction
      int tmp = p1;
      p1 = p2;
      p2 = tmp;
    }
    emplace_back_3d( points[p1], points[p2], points[p3] );
  }
}

PLYSurface::PLYSurface(const std::string & surface_name,
    const stk::diag::Timer &parent_timer,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale)
 : Faceted_Surface_From_File(surface_name, parent_timer),
   my_reader(filename),
   my_dist_sign(sign),
   my_scale(scale)
{
  // Make sure the file is openable, but close it for now to avoid ulimit restrictions.
  my_reader.open_file();
  my_reader.close_file();
}

BoundingBox PLYSurface::get_bounding_box()
{
  BoundingBox bbox;

  my_reader.open_file();

  int num_points = 0;
  int num_facets = 0;

  my_reader.read("reading file header", [this, &num_points, &num_facets](){read_header(num_points, num_facets);});

  std::vector<stk::math::Vector3d> points;
  my_reader.read("reading points", [this, num_points, &points](){read_points(num_points, points);});

  for (auto && point : points)
  {
    bbox.accommodate(point);
  }

  my_reader.close_file();

  bbox.global_reduce();

  return bbox;
}

void PLYSurface::read_file(const std::vector<BoundingBox> & proc_bboxes)
{
  my_reader.open_file();

  int num_points = 0;
  int num_facets = 0;

  my_reader.read("reading file header", [this, &num_points, &num_facets](){read_header(num_points, num_facets);});

  std::vector<stk::math::Vector3d> points;
  my_reader.read("reading points", [this, num_points, &points](){read_points(num_points, points);});

  int batch_size = 0;
  int num_batches = 0;
  my_reader.get_batch_size(num_facets, batch_size, num_batches);

  for (int batch = 0; batch < num_batches; ++batch)
  {
    const int current_batch_size = std::min(batch_size, num_facets-batch*batch_size);
    my_reader.read("reading facets",
        [this, &points, current_batch_size]() { read_facets(current_batch_size, points); });
    parallel_distribute_facets(current_batch_size, proc_bboxes);
  }

  my_reader.close_file();
}

void PLYSurface::read_header(int & num_points, int & num_facets)
{
  std::ifstream & input = my_reader.input();
  num_points = 0;
  num_facets = 0;

  // Read in the file identifier
  std::string symbol;
  input >> symbol;
  STK_ThrowRequire(symbol.compare("ply") == 0);

  while (symbol.compare("end_header") != 0)
  {
    STK_ThrowErrorMsgIf(input.eof(), "Problem reading PLY file, reached end of file.");
    input >> symbol;
    if (symbol.compare("element") == 0)
    {
      input >> symbol;
      if (symbol.compare("vertex") == 0)
      {
        input >> num_points;
      }
      else if (symbol.compare("face") == 0)
      {
        input >> num_facets;
      }
    }
  }
}

void PLYSurface::read_points(const int num_points, std::vector<stk::math::Vector3d> & points)
{
  std::ifstream & input = my_reader.input();
  points.clear();
  points.reserve(num_points);
  for ( int i = 0; i < num_points; i++ )
  {
    double X, Y, Z;
    input >> X >> Y >> Z;
    points.emplace_back(X*my_scale[0], Y*my_scale[1], Z*my_scale[2]);
  }
  // Move to start of next line to prepare for reading facets
  std::string line;
  std::getline(input,line);
}

void PLYSurface::read_facets(const int batch_size, const std::vector<stk::math::Vector3d> & points)
{
  std::ifstream & input = my_reader.input();
  const unsigned num_points = points.size();
  std::string line;
  for ( int i = 0; i < batch_size; ++i )
  {
    unsigned num_facet_nodes;
    unsigned p1, p2, p3;
    std::getline(input,line);
    std::stringstream linestream(line);
    linestream >> num_facet_nodes;
    STK_ThrowRequireMsg(num_facet_nodes == 3, "Failed to read face connectivity correctly.");
    linestream >> p1 >> p2 >> p3;
    STK_ThrowRequire(p1 < num_points);
    STK_ThrowRequire(p2 < num_points);
    STK_ThrowRequire(p3 < num_points);

    if (my_dist_sign == -1)
    {
      // permute to flip normal direction
      int tmp = p1;
      p1 = p2;
      p2 = tmp;
    }
    emplace_back_3d( points[p1], points[p2], points[p3] );
  }
}

STLSurface::STLSurface(const std::string & surface_name,
    const stk::diag::Timer &parent_timer,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale)
 : Faceted_Surface_From_File(surface_name, parent_timer),
   my_is_ascii(true),
   my_reader(filename),
   my_dist_sign(sign),
   my_scale(scale)
{
  // Make sure the file is openable, but close it for now to avoid ulimit restrictions.
  my_reader.open_file();
  my_reader.close_file();
}

static void flip_facet_with_negative_sign(std::array<stk::math::Vector3d,3> & pts, const int sign)
{
  if (sign == -1)
  {
    std::swap(pts[1][0], pts[2][0]);
    std::swap(pts[1][1], pts[2][1]);
    std::swap(pts[1][2], pts[2][2]);
  }
}

static bool is_finite_normal(const stk::math::Vector3d normal)
{
  bool finiteFlag = true;
  if (!std::isfinite(normal[0]) ||
      !std::isfinite(normal[1]) ||
      !std::isfinite(normal[2]) )
  {
    finiteFlag = false;
    krinolog << "Krino::STLSurface found an invalid facet and is skipping." << std::endl;
    krinolog << "  facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
  }
  return finiteFlag;
}

void read_ascii_facet(std::ifstream & input, const stk::math::Vector3d & scale, const int sign, std::array<stk::math::Vector3d,3> & facetPts, bool & isFacetValid)
{
  std::string symbol, nx, ny, nz;
  double X, Y, Z;

  // Read in the facet normal
  input >> symbol;
  STK_ThrowRequire(symbol.compare("normal")  == 0);
  input >> nx >> ny >> nz;
  const stk::math::Vector3d normal(std::atof(nx.c_str()), std::atof(ny.c_str()), std::atof(nz.c_str()));

  // Read in the "outer loop" line
  input >> symbol;
  STK_ThrowRequire(symbol.compare("outer") == 0);
  input >> symbol;
  STK_ThrowRequire(symbol.compare("loop") == 0);

  // Read in the vertices
  for (int i=0; i<3; i++) {
    input >> symbol;
    STK_ThrowRequire(symbol.compare("vertex") == 0);
    input >> X >> Y >> Z;

    facetPts[i] = stk::math::Vector3d(X*scale[0], Y*scale[1], Z*scale[2]);
  }

  // Read in the "endloop" and "endfacet" lines
  input >> symbol;
  STK_ThrowRequire(symbol.compare("endloop") == 0);
  input >> symbol;
  STK_ThrowRequire(symbol.compare("endfacet") == 0);

  isFacetValid = is_finite_normal(normal) && !Facet3d::is_degenerate(facetPts);
  flip_facet_with_negative_sign(facetPts, sign);
}

static bool read_start_of_next_ascii_facet(std::ifstream & input)
{
  // Read the next line
  std::string symbol;
  input >> symbol;
  return symbol.compare("facet") == 0;
}

static float read_binary_float(std::ifstream& input)
{
  float val = 0;
  input.read((char *) &val, sizeof(val));
  return val;
}

static stk::math::Vector3d read_binary_vector(std::ifstream& input)
{
  float x = read_binary_float(input);
  float y = read_binary_float(input);
  float z = read_binary_float(input);
  return stk::math::Vector3d(x, y, z);
}

static void read_binary_facet(std::ifstream & input, const stk::math::Vector3d & scale, const int sign, std::array<stk::math::Vector3d,3> & facetPts, bool & isFacetValid)
{
  const stk::math::Vector3d normal = read_binary_vector(input);

  for (auto && pt : facetPts)
  {
    const stk::math::Vector3d vec = read_binary_vector(input);
    pt = stk::math::Vector3d(vec[0]*scale[0], vec[1]*scale[1], vec[2]*scale[2]);
  }

  char dummy[2];
  input.read(dummy, 2);

  isFacetValid = is_finite_normal(normal) && !Facet3d::is_degenerate(facetPts);
  flip_facet_with_negative_sign(facetPts, sign);
}

BoundingBox
STLSurface::get_bounding_box()
{
  BoundingBox bbox;

  my_reader.open_file();

  my_reader.read("reading file header", [this](){read_header();});
  my_is_ascii = stk::is_true_on_all_procs(stk::EnvData::parallel_comm(), my_is_ascii);

  if (my_is_ascii)
  {
    my_reader.read("reading facets", [this, &bbox](){bbox = get_ascii_facet_bounding_box();});
  }
  else
  {
    unsigned numFacets = 0;
    my_reader.read("reading num facets", [this, &numFacets](){numFacets = read_num_binary_facets();});
    my_reader.read("reading facets", [this, &bbox, numFacets](){bbox = get_binary_facet_bounding_box(numFacets);});
  }

  my_reader.close_file();

  bbox.global_reduce();
  return bbox;
}

void
STLSurface::read_file(const std::vector<BoundingBox> & proc_bboxes)
{
  my_reader.open_file();

  my_reader.read("reading file header", [this](){read_header();});
  my_is_ascii = stk::is_true_on_all_procs(stk::EnvData::parallel_comm(), my_is_ascii);

  unsigned numBinaryFacets = 0;
  if (!my_is_ascii)
    my_reader.read("reading num facets", [this, &numBinaryFacets](){numBinaryFacets = read_num_binary_facets();});

  const unsigned maxBatchSize = 8192;

  unsigned numValidFacetsRead = 0;
  bool done = false;
  while (!done)
  {
    unsigned batchSize = 0;

    if (my_is_ascii)
    {
      my_reader.read("reading facets", [this, &batchSize](){batchSize = read_ascii_facets(maxBatchSize);});
    }
    else
    {
      my_reader.read("reading facets", [this, &batchSize, &numBinaryFacets](){batchSize = read_binary_facets(maxBatchSize, numBinaryFacets);});
    }

    const unsigned localBatchSize = batchSize;
    stk::all_reduce_max( stk::EnvData::parallel_comm(), &localBatchSize, &batchSize, 1 );
    numValidFacetsRead += batchSize;

    parallel_distribute_facets(batchSize, proc_bboxes);

    done = batchSize < maxBatchSize;
  }

  my_reader.close_file();

  krinolog << "Read " << (my_is_ascii ? "ASCII" : "binary") << " STL file " << my_reader.filename() << " with " << numValidFacetsRead << " valid facets." << stk::diag::dendl;
}

void
STLSurface::read_header()
{
  std::ifstream & input = my_reader.input();

  std::string symbol;
  input >> symbol;
  my_is_ascii = symbol.compare("solid") == 0;

  if (my_is_ascii)
  {
    // Read in strings until you get to facet
    while (symbol.compare("facet") != 0)
    {
      STK_ThrowErrorMsgIf(input.eof(), "Problem reading STL file, no facets found.");
      input >> symbol;
    }
    krinolog << "Reading ASCII STL file." << stk::diag::dendl;
  }
  else
  {
    char header_info[80] = "";
    input.clear();
    input.seekg(0);
    input.read(header_info, 80);
    STK_ThrowErrorMsgIf(!input.good(), "Problem reading STL file, cannot read binary file header.");
  }
}

unsigned STLSurface::read_num_binary_facets()
{
  STK_ThrowAssert(!my_is_ascii);
  std::ifstream & input = my_reader.input();
  unsigned numFacets = 0;
  input.read((char *) &numFacets, sizeof(numFacets));
  STK_ThrowErrorMsgIf(!input.good(), "Problem reading STL file, cannot read number of facets.");
  krinolog << "Reading binary STL file with " << numFacets << " facets." << stk::diag::dendl;
  return numFacets;
}

unsigned STLSurface::read_ascii_facets(const unsigned max_batch_size)
{
  std::ifstream & input = my_reader.input();

  bool isFacetValid = false;
  std::array<stk::math::Vector3d,3> facetPts;

  unsigned count = 0;
  bool done = false;
  while (!done)
  {
    read_ascii_facet(input, my_scale, my_dist_sign, facetPts, isFacetValid);
    if (isFacetValid)
    {
      emplace_back_3d(facetPts[0], facetPts[1], facetPts[2]);
      ++count;
    }
    done = (!read_start_of_next_ascii_facet(input) || count >= max_batch_size);
  }
  return count;
}

unsigned STLSurface::read_binary_facets(const unsigned maxBatchSize, unsigned & numRemainingInFile)
{
  std::ifstream & input = my_reader.input();

  bool isFacetValid = false;
  std::array<stk::math::Vector3d,3> facetPts;

  unsigned count = 0;
  bool done = numRemainingInFile == 0;
  while (!done)
  {
    read_binary_facet(input, my_scale, my_dist_sign, facetPts, isFacetValid);
    if (isFacetValid)
    {
      emplace_back_3d(facetPts[0], facetPts[1], facetPts[2]);
      ++count;
    }
    done = (--numRemainingInFile == 0 || count >= maxBatchSize);
  }

  return count;
}

BoundingBox STLSurface::get_ascii_facet_bounding_box()
{
  std::ifstream & input = my_reader.input();

  bool isFacetValid = false;
  std::array<stk::math::Vector3d,3> facetPts;

  BoundingBox bbox;

  bool done = false;
  while (!done)
  {
    read_ascii_facet(input, my_scale, my_dist_sign, facetPts, isFacetValid);
    if (isFacetValid)
    {
      const Facet3d facet(facetPts[0], facetPts[1], facetPts[2]);
      facet.insert_into(bbox);
    }
    done = !read_start_of_next_ascii_facet(input);
  }
  return bbox;
}

BoundingBox STLSurface::get_binary_facet_bounding_box(const unsigned numFacets)
{
  std::ifstream & input = my_reader.input();

  bool isFacetValid = false;
  std::array<stk::math::Vector3d,3> facetPts;

  BoundingBox bbox;

  for (unsigned count=0; count<numFacets; ++count)
  {
    read_binary_facet(input, my_scale, my_dist_sign, facetPts, isFacetValid);
    if (isFacetValid)
    {
      const Facet3d facet(facetPts[0], facetPts[1], facetPts[2]);
      facet.insert_into(bbox);
    }
  }
  return bbox;
}

EXOSurface::EXOSurface(const std::string & surface_name,
    const stk::diag::Timer &parent_timer,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale)
 : Faceted_Surface_From_File(surface_name, parent_timer),
   my_filename(filename),
   my_dist_sign(sign),
   my_scale(scale)
{
}

BoundingBox
EXOSurface::get_bounding_box()
{
  BoundingBox bbox;

  std::vector<double> xyz;
  std::unique_ptr<Ioss::Region> io;
  if ( 0 == stk::EnvData::parallel_rank() )
  {
    /* open file */
    Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodusII", my_filename.c_str(), Ioss::READ_MODEL, MPI_COMM_SELF);
    if ( !db ) {
      ThrowRuntimeError("error reading file " << my_filename);
    }
    io = std::make_unique<Ioss::Region>(db, "EXOSurface IC Region");

    /* exodus description */
    int num_dim   = io->get_property("spatial_dimension").get_int();
    int num_nodes = io->get_property("node_count").get_int();

    STK_ThrowAssert( 3 == num_dim );

    /* generate coordinates */
    xyz.resize(num_nodes*num_dim);

    Ioss::NodeBlock *nb = io->get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", xyz);

    for (int n = 0; n<num_nodes; ++n)
    {
      const stk::math::Vector3d point(&xyz[3*n]);
      bbox.accommodate(point);
    }
  }

  bbox.global_reduce();
  return bbox;
}

void
EXOSurface::read_file(const std::vector<BoundingBox> & proc_bboxes)
{
  const int nodes_per_elem = 3;
  int num_elem_blk = 0;
  std::vector<double> xyz;
  std::vector<int> nmap;
  std::unique_ptr<Ioss::Region> io;
  if ( 0 == stk::EnvData::parallel_rank() )
  {
    /* open file */
    Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodusII", my_filename.c_str(), Ioss::READ_MODEL, MPI_COMM_SELF);
    if ( !db ) {
      ThrowRuntimeError("error reading file " << my_filename);
    }
    io = std::make_unique<Ioss::Region>(db, "EXOSurface IC Region");

    krinolog << "opened file " << my_filename << " for reading..." << std::endl;

    /* exodus description */
    int num_dim   = io->get_property("spatial_dimension").get_int();
    int num_nodes = io->get_property("node_count").get_int();
    int num_elem  = io->get_property("element_count").get_int();
    num_elem_blk = io->get_property("element_block_count").get_int();

    krinolog
        << "num dim = " << num_dim
        << ", num nodes = " << num_nodes
        << ", num elems = " << num_elem
        << ", num elem blks = " << num_elem_blk << std::endl;

    STK_ThrowAssert( 3 == num_dim );

    /* generate coordinates */
    xyz.resize(num_nodes*num_dim);
    nmap.resize(num_nodes);

    Ioss::NodeBlock *nb = io->get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", xyz);
    nb->get_field_data("ids", nmap);
  }

  const int local_num_elem_blk = num_elem_blk;
  stk::all_reduce_max( stk::EnvData::parallel_comm(), &local_num_elem_blk, &num_elem_blk, 1 );


  for (int blk = 0; blk < num_elem_blk; blk++ )
  {
    int num_elem_in_blk = 0;
    std::vector<int> conn;

    if ( 0 == stk::EnvData::parallel_rank() )
    {
      Ioss::ElementBlockContainer ebs = io->get_element_blocks();
      Ioss::ElementBlock *eb = ebs[blk];
      num_elem_in_blk = eb->get_property("entity_count").get_int();

      std::string eb_name = eb->name();
      krinolog
        << "Reading elem blk #" << blk+1
        << " with name " << eb_name << "..." << std::endl;

      int nodes_per_elem_in_blk = eb->get_property("topology_node_count").get_int();
      STK_ThrowAssert( nodes_per_elem_in_blk == nodes_per_elem );

      conn.resize(num_elem_in_blk*nodes_per_elem_in_blk);

      krinolog
        << "  num elem in blk = " << num_elem_in_blk
        << ", nodes per elem in blk = " << nodes_per_elem_in_blk << std::endl;

      eb->get_field_data("connectivity", conn);
    }

    int batch_size = 0;
    int num_batches = 0;
    Parallel_Facet_File_Reader::get_batch_size(num_elem_in_blk, batch_size, num_batches);

    for (int batch = 0; batch < num_batches; ++batch)
    {
      const int current_batch_size = std::min(batch_size, num_elem_in_blk-batch*batch_size);
      if ( 0 == stk::EnvData::parallel_rank() )
      {
        stk::math::Vector3d pt[3];
        for ( int e = 0; e < current_batch_size; e++ )
        {
          for ( int j = 0; j < nodes_per_elem; j++ )
          {
            const int nodeid = conn[e*nodes_per_elem + j];
            const int n = std::distance(nmap.begin(), std::find(nmap.begin(), nmap.end(), nodeid));
            double nodeX = xyz[3*n+0];
            double nodeY = xyz[3*n+1];
            double nodeZ = xyz[3*n+2];
            pt[j] = stk::math::Vector3d( nodeX*my_scale[0], nodeY*my_scale[1], nodeZ*my_scale[2] );
          }

          unsigned p0 = 0;
          unsigned p1 = 1;
          unsigned p2 = 2;
          if (my_dist_sign == -1)
          {
            // permute to flip normal direction
            p1 = 2;
            p2 = 1;
          }

          emplace_back_3d( pt[p0], pt[p1], pt[p2] );
        }
      }
      parallel_distribute_facets(current_batch_size, proc_bboxes);
    }
  }
}

std::unique_ptr<FacetedSurfaceBase> build_mesh_surface(const stk::mesh::MetaData & meta,
    const stk::mesh::Field<double>& coordsField,
    const stk::mesh::Selector & surfaceSelector,
    const int sign)
{
  if (meta.spatial_dimension() == 2)
    return std::make_unique<MeshSurface<Facet2d>>(meta, coordsField, surfaceSelector, sign);
  return std::make_unique<MeshSurface<Facet3d>>(meta, coordsField, surfaceSelector, sign);
}

template <class FACET>
MeshSurface<FACET>::MeshSurface(const stk::mesh::MetaData & meta,
              const stk::mesh::Field<double>& coord_ref,
              const stk::mesh::Selector & surface_selector,
              const int sign)
  : Faceted_Surface<FACET>(),
    my_sign(sign),
    my_mesh_meta(meta),
    my_coord_ref(coord_ref),
    my_surface_selector(surface_selector)
{
}

template <class FACET>
void MeshSurface<FACET>::build_local_facets(const BoundingBox & proc_bbox)
{
  /* %TRACE[ON]% */ Trace trace__("krino::MeshSurface::build_local_facets()"); /* %TRACE% */
  const stk::mesh::BulkData & mesh = my_mesh_meta.mesh_bulk_data();
  stk::mesh::Selector active_locally_owned_selector =
      (NULL != my_mesh_meta.get_part("ACTIVE_CONTEXT_BIT")) ?
      (*my_mesh_meta.get_part("ACTIVE_CONTEXT_BIT") & my_mesh_meta.locally_owned_part()) :
      stk::mesh::Selector(my_mesh_meta.locally_owned_part());

  stk::mesh::Selector active_locally_owned_part_selector = active_locally_owned_selector & my_surface_selector;

  stk::mesh::BucketVector const& buckets = mesh.get_buckets( my_mesh_meta.side_rank(), active_locally_owned_part_selector);

  Faceted_Surface<FACET>::clear();

  stk::mesh::BucketVector::const_iterator ib = buckets.begin();
  stk::mesh::BucketVector::const_iterator ib_end = buckets.end();

  for ( ; ib != ib_end ; ++ib )
  {
    const stk::mesh::Bucket & b = **ib;
    const unsigned length = b.size();

    for ( unsigned iSide = 0; iSide < length; ++iSide )
    {
      stk::mesh::Entity side = b[iSide];

      stk::topology side_topology = mesh.bucket(side).topology();
      if (stk::topology::TRI_3 == side_topology)
      {
        add_facet3d(mesh,side,0,1,2);
      }
      else if (stk::topology::QUAD_4 == side_topology)
      {
        add_quad3d(mesh,side,0,1,2,3);
      }
      else if (stk::topology::TRI_6 == side_topology)
      {
        add_facet3d(mesh,side,0,3,5);
        add_facet3d(mesh,side,1,4,3);
        add_facet3d(mesh,side,2,5,4);
        add_facet3d(mesh,side,3,4,5);
      }
      else if (stk::topology::LINE_2 == side_topology)
      {
        add_facet2d(mesh,side,0,1);
      }
      else if (stk::topology::LINE_3 == side_topology)
      {
        add_facet2d(mesh,side,0,2);
        add_facet2d(mesh,side,2,1);
      }
      else
      {
        ThrowRuntimeError("Elements with side topology " << side_topology.name() << " not supported for mesh surface initialization.");
      }
    }
  }
}

template <class FACET>
void MeshSurface<FACET>::add_facet2d(const stk::mesh::BulkData& mesh, stk::mesh::Entity side, unsigned p0, unsigned p1)
{
  STK_ThrowAssert(2 == my_mesh_meta.spatial_dimension());
  stk::math::Vector3d pt[2];
  
  if (my_sign == -1)
  {
    // permute to flip normal direction
    const unsigned tmp = p0;
    p0 = p1;
    p1 = tmp;
  }
  
  const stk::mesh::Entity* side_nodes = mesh.begin_nodes(side);

  const double * pt0 = stk::mesh::field_data(my_coord_ref, side_nodes[p0]);
  const double * pt1 = stk::mesh::field_data(my_coord_ref, side_nodes[p1]);
  
  pt[0][0] = pt0[0];
  pt[0][1] = pt0[1];
  pt[0][2] = 0.0;
  pt[1][0] = pt1[0];
  pt[1][1] = pt1[1];
  pt[1][2] = 0.0;
  
  add_facet2d(pt[0], pt[1]);
}

template <class FACET>
void MeshSurface<FACET>::add_facet3d(const stk::mesh::BulkData& mesh, stk::mesh::Entity side, unsigned p0, unsigned p1, unsigned p2)
{
  STK_ThrowAssert(3 == my_mesh_meta.spatial_dimension());
  stk::math::Vector3d pt[3];
  
  if (my_sign == -1)
  {
    // permute to flip normal direction
    const unsigned tmp = p1;
    p1 = p2;
    p2 = tmp;
  }
  
  const stk::mesh::Entity* side_nodes = mesh.begin_nodes(side);

  pt[0] = stk::math::Vector3d(stk::mesh::field_data(my_coord_ref, side_nodes[p0]));
  pt[1] = stk::math::Vector3d(stk::mesh::field_data(my_coord_ref, side_nodes[p1]));
  pt[2] = stk::math::Vector3d(stk::mesh::field_data(my_coord_ref, side_nodes[p2]));
  
  add_facet3d(pt[0], pt[1], pt[2]);
}

template <class FACET>
void MeshSurface<FACET>::add_quad3d(const stk::mesh::BulkData& mesh, stk::mesh::Entity side, unsigned p0, unsigned p1, unsigned p2, unsigned p3)
{
  STK_ThrowAssert(3 == my_mesh_meta.spatial_dimension());
  stk::math::Vector3d pt[5];

  if (my_sign == -1)
  {
    // permute to flip normal direction
    const unsigned tmp = p1;
    p1 = p3;
    p3 = tmp;
  }

  const stk::mesh::Entity* side_nodes = mesh.begin_nodes(side);

  pt[0] = stk::math::Vector3d(stk::mesh::field_data(my_coord_ref, side_nodes[p0]));
  pt[1] = stk::math::Vector3d(stk::mesh::field_data(my_coord_ref, side_nodes[p1]));
  pt[2] = stk::math::Vector3d(stk::mesh::field_data(my_coord_ref, side_nodes[p2]));
  pt[3] = stk::math::Vector3d(stk::mesh::field_data(my_coord_ref, side_nodes[p3]));

  // scale is the RMS of the diagonal lengths, sqr_scale is scale*scale, or the average of the lengths squared
  const double sqr_scale = 0.5*((pt[0]-pt[2]).length_squared()+(pt[1]-pt[3]).length_squared());
  // tol is 1e-4, sqr_tol is 1e-8
  const double sqr_tol = 1.e-8;
  // Test is that the distance between the midpoints of the 2 diagonals is less than tol*scale
  const bool is_planar = (0.5*(pt[0]+pt[2])-0.5*(pt[1]+pt[3])).length_squared() < sqr_tol*sqr_scale;
  //const bool is_planar = false; // conservative approach, forces 4 facets per quad

  if (is_planar)
  {
    add_facet3d(pt[0], pt[1], pt[2]);
    add_facet3d(pt[0], pt[2], pt[3]);
  }
  else
  {
    pt[4] = 0.25 * (pt[0]+pt[1]+pt[2]+pt[3]);
    add_facet3d(pt[0], pt[1], pt[4]);
    add_facet3d(pt[1], pt[2], pt[4]);
    add_facet3d(pt[2], pt[3], pt[4]);
    add_facet3d(pt[3], pt[0], pt[4]);
  }
}

template <class FACET>
void MeshSurface<FACET>::add_facet2d(stk::math::Vector3d & p0, stk::math::Vector3d & p1)
{
  STK_ThrowAssert(2 == my_mesh_meta.spatial_dimension());
  Faceted_Surface<FACET>::emplace_back_2d( p0, p1 );
}

template <class FACET>
void MeshSurface<FACET>::add_facet3d(stk::math::Vector3d & p0, stk::math::Vector3d & p1, stk::math::Vector3d & p2)
{
  STK_ThrowAssert(3 == my_mesh_meta.spatial_dimension());
  Faceted_Surface<FACET>::emplace_back_3d( p0, p1, p2 );
}

// Explicit template instantiation
template class MeshSurface<Facet3d>;
template class MeshSurface<Facet2d>;

} // namespace krino
