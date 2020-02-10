#include "stk_ngp_test/ngp_test.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "Kokkos_Core.hpp"
#include <Kokkos_StaticCrsGraph.hpp>
#include "stk_mesh/base/NgpParallelComm.hpp"
#include "stk_mesh/base/NgpSpaces.hpp"

namespace  {

class NgpParallelComm : public ::ngp_testing::Test
{
protected:
  NgpParallelComm()
  {
  }
};

using NeighborCrsViewType = Kokkos::StaticCrsGraph<int, stk::mesh::MemSpace>;
using ValueViewType = Kokkos::View<double*, stk::mesh::MemSpace>;

class ParallelDataExchangeSymPackUnpackHandler
{
public:
  //STK_FUNCTION ParallelDataExchangeSymPackUnpackHandler(stk::mesh::NgpFieldManager & fieldManager)
  ParallelDataExchangeSymPackUnpackHandler(std::vector<std::vector<int>> & neighbors,
                                           ValueViewType & deviceValues)
    : m_neighbors(neighbors),
      m_deviceNeighbors(Kokkos::create_staticcrsgraph<NeighborCrsViewType>("DeviceNeighborsCrs", neighbors)),
      m_deviceValues(deviceValues)
  {
  }

  ParallelDataExchangeSymPackUnpackHandler(const ParallelDataExchangeSymPackUnpackHandler & rhs) = default;

  void hostSizeMessages(int proc, size_t & numValues) const
  {
    numValues = 0;
    for (size_t i = 0; i < m_neighbors.size(); ++i) {
      for (size_t n = 0; n < m_neighbors[i].size(); ++n) {
        if (m_neighbors[i][n] == proc) {
          ++numValues;
        }
      }
    }
  }

  STK_FUNCTION
  void devicePackMessage(int myProc, int proc, ValueViewType & sendData) const
  {
    size_t bufIdx = 0;
    for (size_t i = 0; i < m_deviceNeighbors.numRows(); ++i) {
      const size_t neighborBegin = m_deviceNeighbors.row_map[i];
      const size_t neighborEnd   = m_deviceNeighbors.row_map[i+1];
      for (size_t n = neighborBegin; n < neighborEnd; ++n) {
        if (m_deviceNeighbors.entries[n] == proc) {
          sendData[bufIdx++] = m_deviceValues[i];
        }
      }
    }
  }

  STK_FUNCTION
  void deviceUnpackMessage(int myProc, int proc, ValueViewType & recvData) const
  {
    size_t valueLocation[3];
    size_t valueLocationIdx = 0;
    size_t numSharedValues = 0;
    for (size_t i = 0; i < m_deviceNeighbors.numRows(); ++i) {
      const size_t neighborBegin = m_deviceNeighbors.row_map[i];
      const size_t neighborEnd   = m_deviceNeighbors.row_map[i+1];
      for (size_t n = neighborBegin; n < neighborEnd; ++n) {
        if (m_deviceNeighbors.entries[n] == proc) {
          ++numSharedValues;
          valueLocation[valueLocationIdx++] = i;
        }
      }
    }

    if (numSharedValues == 1u) {
      // Diagonal processor
      m_deviceValues(valueLocation[0]) += recvData[0];
    }
    else {
      // Face-adjacent processor
      m_deviceValues(valueLocation[0]) += recvData[0];

      bool verticalNeighbor = (myProc == proc+1) || (myProc == proc-1);
      if (verticalNeighbor) {
        m_deviceValues(valueLocation[1]) += recvData[2];
        m_deviceValues(valueLocation[2]) += recvData[1];
      }
      else {
        m_deviceValues(valueLocation[1]) += recvData[1];
        m_deviceValues(valueLocation[2]) += recvData[2];
      }
    }
  }

private:
  std::vector<std::vector<int>> m_neighbors;
  NeighborCrsViewType m_deviceNeighbors;
  ValueViewType m_deviceValues;
};

void check_device_values(int pRank, ValueViewType & deviceValues)
{
  Kokkos::parallel_for(deviceValues.size(), KOKKOS_LAMBDA(size_t i)
  {
    //printf("[p%i] deviceValue(%lu) = %f\n", pRank, i, deviceValues(i));
    NGP_EXPECT_NEAR(1.0, deviceValues(i), 1.e-12);
  });
}

void get_domain_shape(const int size, int & nx, int & ny)
{
  nx = 1;
  ny = size;
  int nxTry = std::sqrt(size);
  while (nxTry >= 1) {
    if (size % nxTry == 0) {
      ny = size / nxTry;
      nx = nxTry;
      break;
    }
    --nxTry;
  }
}

std::vector<std::vector<int>> get_neighbors(int pSize, int pRank)
{
  int nx = 0;
  int ny = 0;
  get_domain_shape(pSize, nx, ny);
  int dims[] = {nx, ny};
  int periods[] = {0, 0};
  int reorder = 0;

  // Put together a Cartesian communicator to organize the processors
  // into a grid with well-defined neighbors.
  MPI_Comm comm2D;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2D);

  // Get my coordinates in this processor grid
  int coords[2];
  MPI_Cart_coords(comm2D, pRank, 2, coords);
  const int x = coords[0];
  const int y = coords[1];

  // Query my neighbors in the positive and negative x and y directions.
  // Negative numbers indicate no neighbor in that direction (at the edge
  // of the grid).
  int neighbor[8];
  std::fill_n(neighbor, 8, -1);
  MPI_Cart_shift(comm2D, 0, 1, neighbor+2, neighbor  );
  MPI_Cart_shift(comm2D, 1, 1, neighbor+3, neighbor+1);

  // Query my neighbors in the diagonal directions
  if (neighbor[0] >= 0 && neighbor[1] >= 0) {
    MPI_Cart_rank(comm2D, std::array<int, 2>{{x+1, y+1}}.data(), neighbor+4);
  }
  if (neighbor[1] >= 0 && neighbor[2] >= 0) {
    MPI_Cart_rank(comm2D, std::array<int, 2>{{x-1, y+1}}.data(), neighbor+5);
  }
  if (neighbor[2] >= 0 && neighbor[3] >= 0) {
    MPI_Cart_rank(comm2D, std::array<int, 2>{{x-1, y-1}}.data(), neighbor+6);
  }
  if (neighbor[3] >= 0 && neighbor[0] >= 0) {
    MPI_Cart_rank(comm2D, std::array<int, 2>{{x+1, y-1}}.data(), neighbor+7);
  }

  // We will store 8 values: 1 for each "face" of my subdomain and 1 for each
  // "corner" of my subdomain.  Put together an array of the processors who I
  // will share values with.  For axis-aligned directions (neighbors on a face),
  // there can be only one sharing processor.
  std::vector<std::vector<int>> neighbors(8);
  for (int n = 0; n < 4; ++n) {
    if (neighbor[n] >= 0) {
        neighbors[n].push_back(neighbor[n]);
    }
  }

  // For the diagonal directions, our "corner" value will be shared with either
  // one other processor (we are at the edge of the domain) or three other
  // processors (we are in the middle of the domain).
  for (int n = 0; n < 4; ++n) {
    const int axisDir1 = n;
    const int axisDir2 = (n+1) % 4;
    const int diagDir  = (n+4);

    if (neighbor[axisDir1] >= 0) neighbors[diagDir].push_back(neighbor[axisDir1]);
    if (neighbor[axisDir2] >= 0) neighbors[diagDir].push_back(neighbor[axisDir2]);
    if (neighbor[axisDir1] >= 0 && neighbor[axisDir2] >= 0) {
      neighbors[diagDir].push_back(neighbor[diagDir]);
    }
  }

  return neighbors;
}

NGP_TEST_F(NgpParallelComm, symmetricPackUnpack)
{
  const int pSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int pRank = stk::parallel_machine_rank(MPI_COMM_WORLD);

  std::vector<std::vector<int>> neighbors = get_neighbors(pSize, pRank);

  ValueViewType deviceValues("value", neighbors.size());
  ValueViewType::HostMirror hostValues = Kokkos::create_mirror_view(deviceValues);
  for (size_t n = 0; n < neighbors.size(); ++n) {
    hostValues(n) = 1.0 / (neighbors[n].size() + 1);
  }
  Kokkos::deep_copy(deviceValues, hostValues);

  std::vector<int> commProcs;
  for (size_t i = 0; i < neighbors.size(); ++i) {
    for (const int & n : neighbors[i]) {
      commProcs.push_back(n);
    }
  }
  stk::util::sort_and_unique(commProcs);

  ParallelDataExchangeSymPackUnpackHandler exchangeHandler(neighbors, deviceValues);

  const bool deterministic = false;
  stk::mesh::ngp_parallel_data_exchange_sym_pack_unpack<double>(MPI_COMM_WORLD,
                                                                commProcs,
                                                                exchangeHandler,
                                                                deterministic);

  Kokkos::fence();
  check_device_values(pRank, deviceValues);
}

}
