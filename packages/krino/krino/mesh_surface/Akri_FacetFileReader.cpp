// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_FacetFileReader.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_STLFileParsing.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Ioss_SubSystem.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <unordered_map>

namespace krino
{

namespace
{

// FAC file: the number of points, then "id x y z" for each, with ids 0, 1, ... in order, then the
// number of facets, then "id p1 p2 p3" for each, with p1..p3 point ids.
class FACFileReader : public FacetFileReader
{
public:
  FACFileReader(const std::string & filename, const int sign, const stk::math::Vector3d & scale)
      : myFilename(filename), mySign(sign), myScale(scale)
  {
  }

  void read_vertices(
      std::vector<stk::math::Vector3d> & coords, std::vector<stk::mesh::EntityId> & ids) override
  {
    myInput.open(myFilename);
    STK_ThrowRequireMsg(myInput, "Cannot open file " << myFilename);
    STK_ThrowRequireMsg((myInput >> myNumPoints) && myNumPoints > 0,
        "FAC file " << myFilename << ": cannot read the number of points.");
    coords.clear();
    coords.reserve(myNumPoints);
    for (size_t i = 0; i < myNumPoints; i++)
    {
      size_t id;
      double X, Y, Z;
      STK_ThrowRequireMsg((myInput >> id) && id == i,
          "FAC file " << myFilename << ": point " << i << " is not listed in order.");
      STK_ThrowRequire(myInput >> X >> Y >> Z);
      coords.emplace_back(X * myScale[0], Y * myScale[1], Z * myScale[2]);
    }
    ids.resize(coords.size());
    for (size_t i = 0; i < ids.size(); ++i)
      ids[i] = i + 1;
    STK_ThrowRequireMsg((myInput >> myNumFacets) && myNumFacets > 0,
        "FAC file " << myFilename << ": cannot read the number of facets.");
  }

  size_t num_facets() const override { return myNumFacets; }

  std::array<size_t, 3> read_facet() override
  {
    size_t id;
    size_t p1, p2, p3;
    STK_ThrowRequire((myInput >> id) && id < myNumFacets);
    STK_ThrowRequire(myInput >> p1 >> p2 >> p3);

    STK_ThrowRequire(p1 < myNumPoints);
    STK_ThrowRequire(p2 < myNumPoints);
    STK_ThrowRequire(p3 < myNumPoints);

    if (mySign == -1) return {{p2, p1, p3}}; // reverse the orientation
    return {{p1, p2, p3}};
  }

private:
  std::string myFilename;
  int mySign;
  stk::math::Vector3d myScale;
  std::ifstream myInput;
  size_t myNumPoints{0};
  size_t myNumFacets{0};
};

// PLY file: a header giving the vertex and face counts, then "x y z" for each vertex, then
// "3 p1 p2 p3" for each face, with p1..p3 indices into the vertices.
class PLYFileReader : public FacetFileReader
{
public:
  PLYFileReader(const std::string & filename, const int sign, const stk::math::Vector3d & scale)
      : myFilename(filename), mySign(sign), myScale(scale)
  {
  }

  void read_vertices(
      std::vector<stk::math::Vector3d> & coords, std::vector<stk::mesh::EntityId> & ids) override
  {
    myInput.open(myFilename);
    STK_ThrowRequireMsg(myInput, "Cannot open file " << myFilename);
    read_header();

    coords.clear();
    coords.reserve(myNumPoints);
    for (size_t i = 0; i < myNumPoints; i++)
    {
      double X, Y, Z;
      STK_ThrowRequire(myInput >> X >> Y >> Z);
      coords.emplace_back(X * myScale[0], Y * myScale[1], Z * myScale[2]);
    }
    // Move to start of next line to prepare for reading facets
    std::string line;
    std::getline(myInput, line);

    ids.resize(coords.size());
    for (size_t i = 0; i < ids.size(); ++i)
      ids[i] = i + 1;
  }

  size_t num_facets() const override { return myNumFacets; }

  std::array<size_t, 3> read_facet() override
  {
    std::string line;
    unsigned numFacetNodes = 0;
    size_t p1, p2, p3;
    std::getline(myInput, line);
    std::stringstream linestream(line);
    linestream >> numFacetNodes;
    STK_ThrowRequireMsg(
        numFacetNodes == 3, "PLY file " << myFilename << ": a face is not a triangle.");
    STK_ThrowRequire(linestream >> p1 >> p2 >> p3);
    STK_ThrowRequire(p1 < myNumPoints);
    STK_ThrowRequire(p2 < myNumPoints);
    STK_ThrowRequire(p3 < myNumPoints);

    if (mySign == -1) return {{p2, p1, p3}}; // reverse the orientation
    return {{p1, p2, p3}};
  }

private:
  void read_header()
  {
    myNumPoints = 0;
    myNumFacets = 0;

    // Read in the file identifier
    std::string symbol;
    STK_ThrowRequire((myInput >> symbol) && symbol.compare("ply") == 0);

    while (symbol.compare("end_header") != 0)
    {
      STK_ThrowErrorMsgIf(myInput.eof(), "Problem reading PLY file, reached end of file.");
      STK_ThrowRequire(myInput >> symbol);
      if (symbol.compare("element") == 0)
      {
        STK_ThrowRequire(myInput >> symbol);
        if (symbol.compare("vertex") == 0)
        {
          STK_ThrowRequire(myInput >> myNumPoints);
        }
        else if (symbol.compare("face") == 0)
        {
          STK_ThrowRequire(myInput >> myNumFacets);
        }
      }
    }
  }

  std::string myFilename;
  int mySign;
  stk::math::Vector3d myScale;
  std::ifstream myInput;
  size_t myNumPoints{0};
  size_t myNumFacets{0};
};

// EXO file: a 3D node block, then element blocks of 3-node (trishell3) elements whose connectivity
// gives the file's node ids.
class EXOFileReader : public FacetFileReader
{
public:
  EXOFileReader(const std::string & filename, const int sign, const stk::math::Vector3d & scale)
      : myFilename(filename), mySign(sign), myScale(scale)
  {
  }

  void read_vertices(
      std::vector<stk::math::Vector3d> & coords, std::vector<stk::mesh::EntityId> & ids) override
  {
    // 64-bit ids and connectivity, so files of either integer size are read, including those
    // write_facets(const FacetMesh &) writes.
    Ioss::PropertyManager properties;
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    Ioss::DatabaseIO * db = Ioss::IOFactory::create(
        "exodusII", myFilename.c_str(), Ioss::READ_MODEL, MPI_COMM_SELF, properties);
    STK_ThrowRequireMsg(db, "error reading file " << myFilename);
    myRegion = std::make_unique<Ioss::Region>(db, "facet file region");

    krinolog << "opened file " << myFilename << " for reading..." << stk::diag::dendl;

    const int num_dim = myRegion->get_property("spatial_dimension").get_int();
    const int num_nodes = myRegion->get_property("node_count").get_int();
    const int num_elem = myRegion->get_property("element_count").get_int();
    const int num_elem_blk = myRegion->get_property("element_block_count").get_int();

    krinolog << "num dim = " << num_dim << ", num nodes = " << num_nodes
             << ", num elems = " << num_elem << ", num elem blks = " << num_elem_blk
             << stk::diag::dendl;

    STK_ThrowRequireMsg(3 == num_dim,
        "Exodus surface " << myFilename << " must be 3D, but has dimension " << num_dim << ".");
    STK_ThrowRequireMsg(num_nodes >= 0,
        "Exodus surface " << myFilename << " has an invalid node count " << num_nodes << ".");

    std::vector<double> xyz(size_t(num_nodes) * num_dim);
    std::vector<int64_t> nodeIds(num_nodes);
    Ioss::NodeBlock * nb = myRegion->get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", xyz);
    nb->get_field_data("ids", nodeIds);

    coords.resize(num_nodes);
    ids.resize(num_nodes);
    myIndexOfNodeId.reserve(num_nodes);
    for (int n = 0; n < num_nodes; ++n)
    {
      coords[n] = stk::math::Vector3d(
          xyz[3 * n] * myScale[0], xyz[3 * n + 1] * myScale[1], xyz[3 * n + 2] * myScale[2]);
      ids[n] = nodeIds[n];
      myIndexOfNodeId[nodeIds[n]] = n;
    }

    myNumFacets = 0;
    for (auto * eb : myRegion->get_element_blocks())
      myNumFacets += eb->get_property("entity_count").get_int();
  }

  size_t num_facets() const override { return myNumFacets; }

  std::array<size_t, 3> read_facet() override
  {
    while (myNextElemInBlock >= myBlockConn.size() / nodesPerElem)
      read_next_block();

    const size_t e = myNextElemInBlock++;
    std::array<size_t, 3> facet;
    for (unsigned j = 0; j < nodesPerElem; ++j)
    {
      const auto iter = myIndexOfNodeId.find(myBlockConn[e * nodesPerElem + j]);
      STK_ThrowRequireMsg(iter != myIndexOfNodeId.end(),
          "Exodus surface " << myFilename << " element uses unknown node "
                            << myBlockConn[e * nodesPerElem + j] << ".");
      facet[j] = iter->second;
    }
    if (mySign == -1) std::swap(facet[1], facet[2]); // reverse the orientation
    return facet;
  }

private:
  void read_next_block()
  {
    const Ioss::ElementBlockContainer & blocks = myRegion->get_element_blocks();
    STK_ThrowRequireMsg(myNextBlock < blocks.size(),
        "Exodus surface " << myFilename << ": read past the last facet.");
    Ioss::ElementBlock * eb = blocks[myNextBlock++];
    const int num_elem_in_blk = eb->get_property("entity_count").get_int();
    const int nodes_per_elem_in_blk = eb->get_property("topology_node_count").get_int();

    krinolog << "Reading elem blk #" << myNextBlock << " with name " << eb->name() << "..."
             << stk::diag::dendl;

    STK_ThrowRequireMsg(nodes_per_elem_in_blk == int(nodesPerElem),
        "Exodus surface " << myFilename << " block " << eb->name() << " has "
                          << nodes_per_elem_in_blk << " nodes per element; expected "
                          << nodesPerElem << ".");
    STK_ThrowRequireMsg(num_elem_in_blk >= 0,
        "Exodus surface " << myFilename << " block " << eb->name()
                          << " has an invalid element count " << num_elem_in_blk << ".");

    krinolog << "  num elem in blk = " << num_elem_in_blk
             << ", nodes per elem in blk = " << nodes_per_elem_in_blk << stk::diag::dendl;

    myBlockConn.assign(size_t(num_elem_in_blk) * nodesPerElem, 0);
    eb->get_field_data("connectivity", myBlockConn);
    myNextElemInBlock = 0;
  }

  static constexpr unsigned nodesPerElem = 3;
  std::string myFilename;
  int mySign;
  stk::math::Vector3d myScale;
  std::unique_ptr<Ioss::Region> myRegion;
  std::unordered_map<int64_t, size_t> myIndexOfNodeId;
  size_t myNumFacets{0};
  size_t myNextBlock{0};
  std::vector<int64_t> myBlockConn;
  size_t myNextElemInBlock{0};
};

// STL file: each facet lists its own three vertices.  The vertices are given in facet order, 3 per
// valid facet, so connections between facets come only from merging coincident vertices.
class STLFileReader : public FacetFileReader
{
public:
  STLFileReader(const std::string & filename,
      const int sign,
      const stk::math::Vector3d & scale,
      const bool mergeCoincidentVertices)
      : myFilename(filename),
        mySign(sign),
        myScale(scale),
        myMergeCoincidentVertices(mergeCoincidentVertices)
  {
  }

  void read_vertices(
      std::vector<stk::math::Vector3d> & coords, std::vector<stk::mesh::EntityId> & ids) override
  {
    std::ifstream input(myFilename, std::ios::in | std::ios::binary);
    STK_ThrowRequireMsg(input, "Cannot open file " << myFilename);
    STLFacetReader stlReader(input, mySign, myScale);
    stlReader.read_header();

    // With merging, each distinct point is stored once (as the file is read, so the copies of a
    // vertex in its facets are never all held), else each facet gets its own three vertices.
    std::unordered_map<stk::math::Vector3d, size_t, ExactCoordinatesHash> indexOfPoint;
    auto vertex_index = [&](const stk::math::Vector3d & point)
    {
      if (myMergeCoincidentVertices)
      {
        const auto [iter, inserted] = indexOfPoint.try_emplace(point, coords.size());
        if (inserted) coords.push_back(point);
        return iter->second;
      }
      coords.push_back(point);
      return coords.size() - 1;
    };

    coords.clear();
    myFacets.clear();
    std::array<stk::math::Vector3d, 3> facetPts;
    while (stlReader.read_next_valid_facet(facetPts))
      myFacets.push_back(
          {{vertex_index(facetPts[0]), vertex_index(facetPts[1]), vertex_index(facetPts[2])}});
    ids.resize(coords.size());
    for (size_t i = 0; i < ids.size(); ++i)
      ids[i] = i + 1;
    krinolog << "Read " << (stlReader.is_ascii() ? "ASCII" : "binary") << " STL file " << myFilename
             << " with " << myFacets.size() << " valid facets." << stk::diag::dendl;
  }

  size_t num_facets() const override { return myFacets.size(); }

  std::array<size_t, 3> read_facet() override
  {
    STK_ThrowRequire(myNextFacet < myFacets.size());
    return myFacets[myNextFacet++]; // the orientation is reversed for sign -1 as the facets are
                                    // read
  }

private:
  // Exact equality of the coordinates, as stk::math::Vec ==, with -0 and 0 hashed alike.
  struct ExactCoordinatesHash
  {
    size_t operator()(const stk::math::Vector3d & x) const
    {
      size_t h = 1469598103934665603ull; // FNV-1a
      for (unsigned d = 0; d < 3; ++d)
      {
        const double value = (x[d] == 0.) ? 0. : x[d];
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        h ^= bits;
        h *= 1099511628211ull;
      }
      return h;
    }
  };

  std::string myFilename;
  int mySign;
  stk::math::Vector3d myScale;
  bool myMergeCoincidentVertices;
  std::vector<std::array<size_t, 3>> myFacets;
  size_t myNextFacet{0};
};

} // namespace

void read_on_proc0_and_check_on_all(const stk::ParallelMachine comm,
    const std::string & description,
    const std::string & filename,
    const std::function<void(void)> & read)
{
  std::string error;
  if (0 == stk::parallel_machine_rank(comm))
  {
    try
    {
      read();
    }
    catch (const std::exception & e)
    {
      error = e.what();
    }
  }
  const bool ok = stk::is_true_on_all_procs(comm, error.empty());
  STK_ThrowRequireMsg(ok,
      "Error " << description << " for file " << filename << (error.empty() ? "" : ": " + error));
}

void compute_facet_file_read_batch_size(const stk::ParallelMachine comm,
    const size_t numFacets,
    size_t & batchSize,
    size_t & numBatches)
{
  size_t globalNumFacets = 0;
  stk::all_reduce_sum(comm, &numFacets, &globalNumFacets, 1);

  const size_t minBatchSize = std::min<size_t>(2048, globalNumFacets);
  batchSize =
      std::max<size_t>(minBatchSize, 1 + globalNumFacets / stk::parallel_machine_size(comm));
  numBatches = 1 + globalNumFacets / batchSize;

  if (krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "Reading " << globalNumFacets << " facets, using " << numBatches << " batches of "
             << batchSize << "." << stk::diag::dendl;
}

bool is_facet_file_format(const std::string & format)
{
  return format == "FAC" || format == "PLY" || format == "EXO" || format == "STL";
}

std::unique_ptr<FacetFileReader> create_facet_file_reader(const std::string & format,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale,
    const bool mergeCoincidentVertices)
{
  if (format == "FAC") return std::make_unique<FACFileReader>(filename, sign, scale);
  if (format == "PLY") return std::make_unique<PLYFileReader>(filename, sign, scale);
  if (format == "EXO") return std::make_unique<EXOFileReader>(filename, sign, scale);
  if (format == "STL")
    return std::make_unique<STLFileReader>(filename, sign, scale, mergeCoincidentVertices);
  STK_ThrowErrorMsg("Facet format "
      << format << " is not supported; supported formats are FAC, PLY, EXO and STL.");
  return nullptr;
}

BoundingBox compute_facet_file_bounding_box(const stk::ParallelMachine comm,
    const std::string & format,
    const std::string & filename,
    const int sign,
    const stk::math::Vector3d & scale)
{
  BoundingBox bbox;
  read_on_proc0_and_check_on_all(comm,
      "reading points",
      filename,
      [&]()
      {
        std::vector<stk::math::Vector3d> coords;
        std::vector<stk::mesh::EntityId> ids;
        create_facet_file_reader(format, filename, sign, scale, true)->read_vertices(coords, ids);
        for (auto && point : coords)
          bbox.accommodate(point);
      });
  // Only processor 0 read the file, so it has the box; share it with comm's processors.
  BoundingBox::VecType globalMin, globalMax;
  stk::all_reduce_min(comm, bbox.get_min().data(), globalMin.data(), 3);
  stk::all_reduce_max(comm, bbox.get_max().data(), globalMax.data(), 3);
  return BoundingBox(globalMin, globalMax);
}

} // namespace krino
