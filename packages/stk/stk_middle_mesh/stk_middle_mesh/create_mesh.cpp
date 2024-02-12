#include "create_mesh.hpp"
#include "utils.hpp"
#include <iostream>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

int MeshGenerator::get_idx(const MeshSpec& spec, const int i, const int j)
{
  return i * (spec.numelY + 1) + j;
}

std::array<int, 2> MeshGenerator::get_process_grid_shape(int commsize)
{
  std::array<int, 2> gridShape{0, 0};
  MPI_Dims_create(commsize, 2, gridShape.data());

  return gridShape;
}

std::array<int, 2> MeshGenerator::compute_local_indices(const std::array<int, 2>& gridShape, int myrank)
{
  assert(myrank < gridShape[0] * gridShape[1]);

  int xIdx = myrank % gridShape[0];
  int yIdx = myrank / gridShape[0];

  return std::array<int, 2>{xIdx, yIdx};
}

int MeshGenerator::wrap_processor_index(int idx, int size)
{
  if (idx == size)
    idx = 0;
  else if (idx == -1)
    idx = size - 1;

  return idx;
}

std::array<int, 2> MeshGenerator::wrap_processor_indices(const std::array<int, 2>& gridIndices)
{
  int yidx = gridIndices[1];
  if (m_meshspecGlobal.yPeriodic)
    yidx = wrap_processor_index(yidx, m_processorGridShape[1]);

  int xidx = gridIndices[0];
  if (m_meshspecGlobal.xPeriodic)
    xidx = wrap_processor_index(xidx, m_processorGridShape[0]);

  return {xidx, yidx};
}

bool MeshGenerator::does_block_exist(const std::array<int, 2>& gridIndices)
{
  auto gridIndicesWrapped = wrap_processor_indices(gridIndices);

  return gridIndicesWrapped[0] >= 0 && gridIndicesWrapped[0] < m_processorGridShape[0] && gridIndicesWrapped[1] >= 0 &&
         gridIndicesWrapped[1] < m_processorGridShape[1];
}

int MeshGenerator::get_process_rank(const std::array<int, 2>& gridIndices)
{
  auto gridIndicesWrapped = wrap_processor_indices(gridIndices);
  assert(gridIndicesWrapped[0] >= 0 && gridIndicesWrapped[0] < m_processorGridShape[0]);
  assert(gridIndicesWrapped[1] >= 0 && gridIndicesWrapped[1] < m_processorGridShape[1]);
  return m_processorGridShape[0] * gridIndicesWrapped[1] + gridIndicesWrapped[0];
}

bool MeshGenerator::is_x_periodic_local()
{
  std::array<int, 2> rightNeighborBlock{m_processIndices[0] + 1, m_processIndices[1]};
  rightNeighborBlock    = wrap_processor_indices(rightNeighborBlock);
  bool isXperiodicLocal = rightNeighborBlock[0] == m_processIndices[0];  

  return isXperiodicLocal;
}

bool MeshGenerator::is_y_periodic_local()
{
  std::array<int, 2> topNeighborBlock{m_processIndices[0], m_processIndices[1] + 1};
  topNeighborBlock      = wrap_processor_indices(topNeighborBlock);
  bool isYperiodicLocal = topNeighborBlock[1] == m_processIndices[1];  

  return isYperiodicLocal;
}

MeshSpec MeshGenerator::compute_local_mesh_spec()
{
  auto valsX = get_num_el_in_direction(m_meshspecGlobal.numelX, m_meshspecGlobal.xmin, m_meshspecGlobal.xmax,
                                       m_processorGridShape[0], m_processIndices[0]);
  auto valsY = get_num_el_in_direction(m_meshspecGlobal.numelY, m_meshspecGlobal.ymin, m_meshspecGlobal.ymax,
                                       m_processorGridShape[1], m_processIndices[1]);

  double deltaX = (m_meshspecGlobal.xmax - m_meshspecGlobal.xmin) / m_meshspecGlobal.numelX;
  double deltaY = (m_meshspecGlobal.ymax - m_meshspecGlobal.ymin) / m_meshspecGlobal.numelY;
  int numelX = valsX.first, numelY = valsY.first;
  double xmin = valsX.second, ymin = valsY.second;
  double xmax = xmin + deltaX * numelX;
  double ymax = ymin + deltaY * numelY;

  bool isXperiodic =
      m_meshspecGlobal.xPeriodic && (m_processIndices[0] == 0 || m_processIndices[0] == m_processorGridShape[0] - 1);

  bool isYperiodic =
      m_meshspecGlobal.yPeriodic && (m_processIndices[1] == 0 || m_processIndices[1] == m_processorGridShape[1] - 1);

  return MeshSpec{numelX, numelY, xmin, xmax, ymin, ymax, isXperiodic, isYperiodic};
}

std::pair<int, double> MeshGenerator::get_num_el_in_direction(int numelGlobal, double coordMin, double coordMax,
                                                              int numProcs, int procIdx)
{
  int numelPerProc                  = numelGlobal / numProcs;
  int numelRemaining                = numelGlobal - numelPerProc * numProcs;
  int numLowerProcsWithExtraElement = numelRemaining;
  if (procIdx < numelRemaining)
  {
    numelPerProc += 1;
    numLowerProcsWithExtraElement = 0;
  }

  double deltaX          = (coordMax - coordMin) / numelGlobal;
  double coordStartLocal = coordMin + deltaX * (numelPerProc * procIdx + numLowerProcsWithExtraElement);

  return {numelPerProc, coordStartLocal};
}

std::vector<MeshEntityPtr> MeshGenerator::create_vertices(const MeshSpec& spec,
                                                          const std::vector<utils::Point>& vertCoords)
{
  assert(spec.numelX > 0);
  assert(spec.numelY > 0);

  std::vector<MeshEntityPtr> verts;
  verts.resize((spec.numelX + 1) * (spec.numelY + 1));

  bool isXperiodicLocal = is_x_periodic_local();
  bool isYperiodicLocal = is_y_periodic_local();

  for (int i = 0; i < spec.numelX + 1; ++i)
    for (int j = 0; j < spec.numelY + 1; ++j)
    {
      int idx = get_idx(spec, i, j);
      if (isXperiodicLocal && spec.xPeriodic && i == spec.numelX)
      {
        verts[idx] = verts[get_idx(spec, 0, j)];
      } else if (isYperiodicLocal && spec.yPeriodic && j == spec.numelY)
      {
        verts[idx] = verts[get_idx(spec, i, 0)];
      } else
      {
        verts[idx] = m_mesh->create_vertex(vertCoords[idx]);
      }
    }

  return verts;
}

MeshGenerator::Edges MeshGenerator::create_edges(const MeshSpec& spec, const std::vector<MeshEntityPtr>& verts,
                                                 bool triangles, bool reverseXEdges)
{
  Edges edges(spec);

  bool isXperiodicLocal = is_x_periodic_local();
  bool isYperiodicLocal = is_y_periodic_local();

  // x direction edges
  for (int i = 0; i < spec.numelX; ++i)
    for (int j = 0; j < spec.numelY + 1; ++j)
    {
      if (isYperiodicLocal && spec.yPeriodic && j == spec.numelY)
      {
        edges.get_x_edge(i, j) = edges.get_x_edge(i, 0);
      } else
      {
        MeshEntityPtr v1         = verts[get_idx(spec, i, j)];
        MeshEntityPtr v2         = verts[get_idx(spec, i + 1, j)];
        EntityOrientation orient = EntityOrientation::Standard;

        if (reverseXEdges && j > 0 && j <= spec.numelY)
        {
          auto vTmp = v2;
          v2        = v1;
          v1        = vTmp;
          orient    = EntityOrientation::Reversed;
        }

        edges.get_x_edge(i, j)             = m_mesh->create_edge(v1, v2);
        edges.get_x_edge_orientation(i, j) = orient;
      }
    }

  // y direction edges
  for (int i = 0; i < spec.numelX + 1; ++i)
    for (int j = 0; j < spec.numelY; ++j)
    {
      if ((isXperiodicLocal && spec.xPeriodic && i == spec.numelX))
      {
        edges.get_y_edge(i, j) = edges.get_y_edge(0, j);
      } else
      {
        MeshEntityPtr v1       = verts[get_idx(spec, i, j)];
        MeshEntityPtr v2       = verts[get_idx(spec, i, j + 1)];
        edges.get_y_edge(i, j) = m_mesh->create_edge(v1, v2);
      }
    }

  if (triangles)
  {
    for (int i = 0; i < spec.numelX; ++i)
      for (int j = 0; j < spec.numelY; ++j)
      {
        MeshEntityPtr v1          = verts[get_idx(spec, i, j)];
        MeshEntityPtr v2          = verts[get_idx(spec, i + 1, j + 1)];
        edges.get_diag_edge(i, j) = m_mesh->create_edge(v1, v2);
      }
  }

  return edges;
}

// create mesh of quads from a vector of utils::Points from createVertices
void MeshGenerator::create_mesh_elements(const MeshSpec& spec, Edges& edges, const bool triangles)
{
  assert(spec.numelX > 0);
  assert(spec.numelY > 0);

  if (triangles)
  {
    for (int i = 0; i < spec.numelX; ++i)
      for (int j = 0; j < spec.numelY; ++j)
      {
        auto e1 = edges.get_x_edge(i, j);
        auto e2 = edges.get_y_edge(i + 1, j);
        auto e3 = edges.get_x_edge(i, j + 1);
        auto e4 = edges.get_y_edge(i, j);
        auto e5 = edges.get_diag_edge(i, j);

        EntityOrientation e1Orient = edges.get_x_edge_orientation(i, j);
        EntityOrientation e3Orient = edges.get_x_edge_orientation(i, j + 1);
        m_mesh->create_triangle(e1, e2, e5, e1Orient);
        m_mesh->create_triangle(e5, e3, e4, e3Orient);
      }
  } else
  {
    for (int i = 0; i < spec.numelX; ++i)
      for (int j = 0; j < spec.numelY; ++j)
      {
        auto e1 = edges.get_x_edge(i, j);
        auto e2 = edges.get_y_edge(i + 1, j);
        auto e3 = edges.get_x_edge(i, j + 1);
        auto e4 = edges.get_y_edge(i, j);

        EntityOrientation e1Orient = edges.get_x_edge_orientation(i, j);
        m_mesh->create_quad(e1, e2, e3, e4, e1Orient);
      }
  }
}

void MeshGenerator::set_remote_shared_entities(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts)
{
  if (utils::impl::comm_size(m_comm) == 1)
    return;

  set_edge_adjacent_block_shared_entities(mesh, verts);
  set_corner_adjacent_block_shared_entities(mesh, verts);
}

class EdgeNeighborData
{
  public:
    EdgeNeighborData(MPI_Comm comm, int destRank, int sendTagVerts, int sendTagEdges, int recvTagVerts,
                     int recvTagEdges)
      : m_comm(comm)
      , m_destRank(destRank)
      , m_sendTagVerts(sendTagVerts)
      , m_sendTagEdges(sendTagEdges)
      , m_recvTagVerts(recvTagVerts)
      , m_recvTagEdges(recvTagEdges)
    {}

    std::vector<MeshEntityPtr> verts;
    std::vector<MeshEntityPtr> edges;

    void start_communication()
    {
      if (m_destRank < 0)
        return;

      for (auto& vert : verts)
        m_vertIdsSendBuf.push_back(vert->get_id());
      m_vertIdsRecvBuf.resize(verts.size());

      for (auto& edge : edges)
        m_edgeIdsSendBuf.push_back(edge->get_id());
      m_edgeIdsRecvBuf.resize(edges.size());

      MPI_Isend(m_vertIdsSendBuf.data(), m_vertIdsSendBuf.size(), MPI_INT, m_destRank, m_sendTagVerts, m_comm,
                &m_vertSendReq);
      MPI_Isend(m_edgeIdsSendBuf.data(), m_edgeIdsSendBuf.size(), MPI_INT, m_destRank, m_sendTagEdges, m_comm,
                &m_edgeSendReq);

      MPI_Irecv(m_vertIdsRecvBuf.data(), m_vertIdsRecvBuf.size(), MPI_INT, m_destRank, m_recvTagVerts, m_comm,
                &m_vertRecvReq);
      MPI_Irecv(m_edgeIdsRecvBuf.data(), m_edgeIdsRecvBuf.size(), MPI_INT, m_destRank, m_recvTagEdges, m_comm,
                &m_edgeRecvReq);
    }

    void finish_communication()
    {
      if (m_destRank < 0)
        return;

      MPI_Wait(&m_vertRecvReq, MPI_STATUS_IGNORE);
      MPI_Wait(&m_edgeRecvReq, MPI_STATUS_IGNORE);

      for (size_t i = 0; i < m_vertIdsRecvBuf.size(); ++i)
      {
        auto vert = verts[i];
        vert->add_remote_shared_entity({m_destRank, m_vertIdsRecvBuf[i]});
      }

      for (size_t i = 0; i < m_edgeIdsRecvBuf.size(); ++i)
      {
        auto edge = edges[i];
        edge->add_remote_shared_entity({m_destRank, m_edgeIdsRecvBuf[i]});
      }

      MPI_Wait(&m_vertSendReq, MPI_STATUS_IGNORE);
      MPI_Wait(&m_edgeSendReq, MPI_STATUS_IGNORE);
    }

  private:
    MPI_Comm m_comm;
    int m_destRank;
    int m_sendTagVerts;
    int m_sendTagEdges;
    int m_recvTagVerts;
    int m_recvTagEdges;
    std::vector<int> m_vertIdsSendBuf;
    std::vector<int> m_edgeIdsSendBuf;
    std::vector<int> m_vertIdsRecvBuf;
    std::vector<int> m_edgeIdsRecvBuf;
    MPI_Request m_vertSendReq;
    MPI_Request m_edgeSendReq;
    MPI_Request m_vertRecvReq;
    MPI_Request m_edgeRecvReq;
};

void MeshGenerator::set_edge_adjacent_block_shared_entities(std::shared_ptr<Mesh> mesh,
                                                            std::vector<MeshEntityPtr>& verts)
{
  int procX = m_processIndices[0], procY = m_processIndices[1];
  std::vector<std::array<int, 2>> neighborBlocks = {
      {procX - 1, procY}, {procX + 1, procY}, {procX, procY - 1}, {procX, procY + 1}};
  std::array<int, 4> edgeMap = {1, 0, 3, 2}; // map edge from receivers perspective to edge from senders perspective
  int vertTagStart = 100, edgeTagStart = 104;
  std::vector<EdgeNeighborData> data;
  for (int i = 0; i < 4; ++i)
  {
    int destRank = does_block_exist(neighborBlocks[i]) ? get_process_rank(neighborBlocks[i]) : -1;
    if (destRank == utils::impl::comm_rank(m_comm))
      destRank = -1;

    int destEdge = edgeMap[i];
    data.emplace_back(mesh->get_comm(), destRank, vertTagStart + i, edgeTagStart + i, vertTagStart + destEdge,
                      edgeTagStart + destEdge);
  }

  get_y_shared_entities(data[0], data[1], verts);
  get_x_shared_entities(data[2], data[3], verts);

  for (int i = 0; i < 4; ++i)
    if (does_block_exist(neighborBlocks[i]))
      data[i].start_communication();

  for (int i = 0; i < 4; ++i)
    if (does_block_exist(neighborBlocks[i]))
      data[i].finish_communication();
}

void MeshGenerator::get_y_shared_entities(EdgeNeighborData& leftEdge, EdgeNeighborData& rightEdge,
                                          std::vector<MeshEntityPtr>& verts)
{
  int numVertsY = m_meshspecLocal.numelY;
  if (is_y_periodic_local())
    numVertsY--;

  for (int j=0; j <= numVertsY; ++j)
  {
    leftEdge.verts.push_back(verts[get_idx(m_meshspecLocal, 0, j)]);
    rightEdge.verts.push_back(verts[get_idx(m_meshspecLocal, m_meshspecLocal.numelX, j)]);
  }

  for (int j=1; j <= m_meshspecLocal.numelY; ++j)
  {
    leftEdge.edges.push_back(get_common_edge(verts[get_idx(m_meshspecLocal, 0, j-1)],
                                            verts[get_idx(m_meshspecLocal, 0, j)]));
    rightEdge.edges.push_back(get_common_edge(verts[get_idx(m_meshspecLocal, m_meshspecLocal.numelX, j-1)],
                                            verts[get_idx(m_meshspecLocal, m_meshspecLocal.numelX, j)]));                                            
  }  
}

void MeshGenerator::get_x_shared_entities(EdgeNeighborData& bottomEdge, EdgeNeighborData& topEdge,
                                          std::vector<MeshEntityPtr>& verts)
{
  int numVertsX = m_meshspecLocal.numelX;
  if (is_x_periodic_local())
    numVertsX--;

  for (int i=0; i <= numVertsX; ++i)
  {
    bottomEdge.verts.push_back(verts[get_idx(m_meshspecLocal, i, 0)]);
    topEdge.verts.push_back(verts[get_idx(m_meshspecLocal, i, m_meshspecLocal.numelY)]);
  }

  for (int i=1; i <= m_meshspecLocal.numelX; ++i)
  {
    bottomEdge.edges.push_back(get_common_edge(verts[get_idx(m_meshspecLocal, i-1, 0)],
                                            verts[get_idx(m_meshspecLocal, i,   0)]));
    topEdge.edges.push_back(get_common_edge(verts[get_idx(m_meshspecLocal, i-1, m_meshspecLocal.numelY)],
                                            verts[get_idx(m_meshspecLocal, i,   m_meshspecLocal.numelY)]));                                            
  }    
}                                    

void MeshGenerator::set_corner_adjacent_block_shared_entities(std::shared_ptr<Mesh> mesh,
                                                              std::vector<MeshEntityPtr>& verts)
{
  if (is_x_periodic_local() || is_y_periodic_local())  
    return;

  int numelX = m_meshspecLocal.numelX, numelY = m_meshspecLocal.numelY;
  MeshEntityPtr topRightId    = verts[get_idx(m_meshspecLocal, numelX, numelY)];
  MeshEntityPtr bottomRightId = verts[get_idx(m_meshspecLocal, numelX, 0)];
  MeshEntityPtr bottomLeftId  = verts[get_idx(m_meshspecLocal, 0, 0)];
  MeshEntityPtr topLeftId     = verts[get_idx(m_meshspecLocal, 0, numelY)];

  int procX = m_processIndices[0], procY = m_processIndices[1];
  std::vector<std::array<int, 2>> neighborBlocks = {
      {procX + 1, procY + 1}, {procX + 1, procY - 1}, {procX - 1, procY - 1}, {procX - 1, procY + 1}};
  std::array<MeshEntityPtr, 4> vertsCorner{topRightId, bottomRightId, bottomLeftId, topLeftId};
  std::array<int, 4> idsSend, idsRecv;
  std::array<MPI_Request, 4> sendReqs, recvReqs;
  int startTag                 = 102;
  std::array<int, 4> cornerMap = {2, 3, 0, 1}; // map corner id from my perpective to sender's perspective

  for (int i = 0; i < 4; ++i)
    if (does_block_exist(neighborBlocks[i]))
    {
      idsSend[i] = vertsCorner[i]->get_id();
      MPI_Isend(&(idsSend[i]), 1, MPI_INT, get_process_rank(neighborBlocks[i]), startTag + i, mesh->get_comm(),
                &(sendReqs[i]));
      MPI_Irecv(&(idsRecv[i]), 1, MPI_INT, get_process_rank(neighborBlocks[i]), startTag + cornerMap[i],
                mesh->get_comm(), &(recvReqs[i]));
    }

  for (int i = 0; i < 4; ++i)
    if (does_block_exist(neighborBlocks[i]))
    {
      MPI_Wait(&(recvReqs[i]), MPI_STATUS_IGNORE);
      vertsCorner[i]->add_remote_shared_entity({get_process_rank(neighborBlocks[i]), idsRecv[i]});
    }

  for (int i = 0; i < 4; ++i)
    if (does_block_exist(neighborBlocks[i]))
      MPI_Wait(&(sendReqs[i]), MPI_STATUS_IGNORE);
}

std::shared_ptr<Mesh> create_eigth_sphere(int numelR, int numelTheta, double rInner, double rOuter, MPI_Comm comm, bool createTriangles)
{
  double pi = std::atan(1) * 4;

  MeshSpec spec;
  spec.numelX = numelR;
  spec.numelY = numelTheta;
  spec.xmin   = rInner;
  spec.xmax   = rOuter;
  spec.ymin   = 0;
  spec.ymax   = pi / 2;

  double radius = spec.xmax;
  auto func     = [&](const utils::Point& pt) {
    // interpret x and y as r and theta
    double r     = pt.x;
    double theta = pt.y;

    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    // project upward onto surface of sphere
    double z = std::sqrt(std::max(radius * radius - x * x - y * y, 0.0));
    // double z = 0.0;
    utils::Point pt2(x, y, z);
    return pt2;
  };

  return create_mesh(spec, func, comm, createTriangles);
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
