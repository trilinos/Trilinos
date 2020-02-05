
#ifndef _ZOLTAN2_COORD_PARTITIONMAPPING_HPP_
#define _ZOLTAN2_COORD_PARTITIONMAPPING_HPP_

#include <fstream>
#include <ctime>
#include <vector>
#include <set>
#include <tuple>

#include "Zoltan2_Standards.hpp"
#include "Zoltan2_AlgMultiJagged.hpp"
#include "Teuchos_ArrayViewDecl.hpp"
#include "Zoltan2_PartitionMapping.hpp"
#include "Zoltan2_MachineRepresentation.hpp"
#include "Teuchos_ReductionOp.hpp"

#include "Zoltan2_MappingSolution.hpp"

#include "Zoltan2_GraphModel.hpp"

#include <Zoltan2_TPLTraits.hpp>

#include "Teuchos_Comm.hpp"
#ifdef HAVE_ZOLTAN2_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_ZOLTAN2_MPI
#include <Teuchos_DefaultSerialComm.hpp>

//#define gnuPlot
#include "Zoltan2_XpetraMultiVectorAdapter.hpp"
#include <Zoltan2_Directory_Impl.hpp>

namespace Teuchos{

/*! \brief Zoltan2_ReduceBestMapping Class, reduces the minimum cost
 *  mapping, ties breaks with minimum proc id.
 */
template <typename Ordinal, typename T>
class Zoltan2_ReduceBestMapping  : public ValueTypeReductionOp<Ordinal,T>
{

private:
  T epsilon;

public:
  /*! \brief Default Constructor
   */
  Zoltan2_ReduceBestMapping() : epsilon(std::numeric_limits<T>::epsilon()) {}

  /*! \brief Implement Teuchos::ValueTypeReductionOp interface
   */
  void reduce(const Ordinal count,
              const T inBuffer[],
              T inoutBuffer[]) const {

    for (Ordinal i = 0; i < count; i++) {
      if (inBuffer[0] - inoutBuffer[0] < -epsilon) {
        inoutBuffer[0] = inBuffer[0];
        inoutBuffer[1] = inBuffer[1];
      } else if(inBuffer[0] - inoutBuffer[0] < epsilon &&
          inBuffer[1] - inoutBuffer[1] < epsilon) {
        inoutBuffer[0] = inBuffer[0];
        inoutBuffer[1] = inBuffer[1];
      }
    }
  }
};

} // namespace Teuchos


namespace Zoltan2{

template <typename it>
inline it z2Fact(it x) {
  return (x == 1 ? x : x * z2Fact<it>(x - 1));
}

template <typename gno_t, typename part_t>
class GNO_LNO_PAIR{
public:
  gno_t gno;
  part_t part;
};

//returns the ith permutation indices.
template <typename IT>
void ithPermutation(const IT n, IT i, IT *perm)
{
  IT j, k = 0;
  IT *fact = new IT[n];


  // compute factorial numbers
  fact[k] = 1;
  while (++k < n)
    fact[k] = fact[k - 1] * k;

  // compute factorial code
  for (k = 0; k < n; ++k)
  {
    perm[k] = i / fact[n - 1 - k];
    i = i % fact[n - 1 - k];
  }

  // readjust values to obtain the permutation
  // start from the end and check if preceding values are lower
  for (k = n - 1; k > 0; --k)
    for (j = k - 1; j >= 0; --j)
      if (perm[j] <= perm[k])
        perm[k]++;

  delete [] fact;
}

template <typename part_t>
void getGridCommunicationGraph(part_t taskCount,
                               part_t *&task_comm_xadj,
                               part_t *&task_comm_adj,
                               std::vector<int> grid_dims) {
  int dim = grid_dims.size();
  int neighborCount = 2 * dim;
  task_comm_xadj = new part_t[taskCount + 1];
  task_comm_adj = new part_t[taskCount * neighborCount];

  part_t neighBorIndex = 0;
  task_comm_xadj[0] = 0;
  for (part_t i = 0; i < taskCount; ++i) {
    part_t prevDimMul = 1;
    for (int j = 0; j < dim; ++j) {
      part_t lNeighbor = i - prevDimMul;
      part_t rNeighbor = i + prevDimMul;
      prevDimMul *= grid_dims[j];
      if (lNeighbor >= 0 &&
          lNeighbor/ prevDimMul == i / prevDimMul &&
          lNeighbor < taskCount) {
        task_comm_adj[neighBorIndex++] = lNeighbor;
      }
      if (rNeighbor >= 0 &&
          rNeighbor/ prevDimMul == i / prevDimMul &&
          rNeighbor < taskCount) {
        task_comm_adj[neighBorIndex++] = rNeighbor;
      }
    }
    task_comm_xadj[i + 1] = neighBorIndex;
  }

}
//returns the center of the parts.
template <typename Adapter, typename scalar_t, typename part_t>
void getSolutionCenterCoordinates(
    const Environment *envConst,
    const Teuchos::Comm<int> *comm,
    const Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *coords,
    //const Zoltan2::PartitioningSolution<Adapter> *soln_,
    const part_t *parts,
    int coordDim,
    part_t ntasks,
    scalar_t **partCenters) {

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;

  typedef StridedData<lno_t, scalar_t> input_t;
  ArrayView<const gno_t> gnos;
  ArrayView<input_t>     xyz;
  ArrayView<input_t>     wgts;
  coords->getCoordinates(gnos, xyz, wgts);

  //local and global num coordinates.
  lno_t numLocalCoords = coords->getLocalNumCoordinates();
  //gno_t numGlobalCoords = coords->getGlobalNumCoordinates();

  //local number of points in each part.
  gno_t *point_counts = new gno_t[ntasks];
  memset(point_counts, 0, sizeof(gno_t) * ntasks);

  //global number of points in each part.
  gno_t *global_point_counts = new gno_t[ntasks];

  scalar_t **multiJagged_coordinates = new scalar_t*[coordDim];

  ArrayRCP<ArrayRCP<const scalar_t>> ar = arcp(new ArrayRCP<const scalar_t>[coordDim], 0, coordDim);
  for (int dim=0; dim < coordDim; dim++){
    xyz[dim].getInputArray(ar[dim]);
    //multiJagged coordinate values assignment
    multiJagged_coordinates[dim] = (scalar_t *)ar[dim].getRawPtr();
    memset(partCenters[dim], 0, sizeof(scalar_t) * ntasks);
  }

  //get parts with parallel gnos.
  //const part_t *parts = soln_->getPartListView();

  envConst->timerStart(MACRO_TIMERS, "Mapping - Center Calculation");

  for (lno_t i=0; i < numLocalCoords; i++) {
    part_t p = parts[i];
    //add up all coordinates in each part.
    for(int j = 0; j < coordDim; ++j) {
      scalar_t c = multiJagged_coordinates[j][i];
      partCenters[j][p] += c;
    }
    ++point_counts[p];
  }

  //get global number of points in each part.
  reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_SUM,
                        ntasks, point_counts, global_point_counts);

  for(int j = 0; j < coordDim; ++j) {
    for (part_t i=0; i < ntasks; ++i) {
      if (global_point_counts[i] > 0)
        partCenters[j][i] /= global_point_counts[i];
    }
  }

  scalar_t *tmpCoords = new scalar_t[ntasks];
  for(int j = 0; j < coordDim; ++j) {
    reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM,
                             ntasks, partCenters[j], tmpCoords);

    scalar_t *tmp = partCenters[j];
    partCenters[j] = tmpCoords;
    tmpCoords = tmp;
  }

  envConst->timerStop(MACRO_TIMERS, "Mapping - Center Calculation");

  delete [] point_counts;
  delete [] global_point_counts;

  delete [] tmpCoords;
  delete [] multiJagged_coordinates;
}

//returns the coarsend part graph.
template <typename Adapter, typename scalar_t, typename part_t>
void getCoarsenedPartGraph(
    const Environment *envConst,
    const Teuchos::Comm<int> *comm,
    const Zoltan2::GraphModel<typename Adapter::base_adapter_t> *graph,
    //const Zoltan2::PartitioningSolution<Adapter> *soln_,
    part_t np,
    const part_t *parts,
    ArrayRCP<part_t> &g_part_xadj,
    ArrayRCP<part_t> &g_part_adj,
    ArrayRCP<scalar_t> &g_part_ew) {

  typedef typename Adapter::lno_t t_lno_t;
  typedef typename Adapter::gno_t t_gno_t;
  typedef typename Adapter::scalar_t t_scalar_t;
  typedef typename Adapter::offset_t t_offset_t;
  typedef typename Zoltan2::GraphModel<
    typename Adapter::base_adapter_t>::input_t t_input_t;

  //int numRanks = comm->getSize();
  //int myRank = comm->getRank();

  //get parts with parallel gnos.
  /*
  const part_t *parts = soln_->getPartListView();

  part_t np = soln_->getActualGlobalNumberOfParts();
  if (part_t(soln_->getTargetGlobalNumberOfParts()) > np) {
    np = soln_->getTargetGlobalNumberOfParts();
  }
  */


  t_lno_t localNumVertices = graph->getLocalNumVertices();
  t_lno_t localNumEdges = graph->getLocalNumEdges();

  //get the vertex global ids, and weights
  ArrayView<const t_gno_t> Ids;
  ArrayView<t_input_t> v_wghts;
  graph->getVertexList(Ids, v_wghts);

  //get the edge ids, and weights
  ArrayView<const t_gno_t> edgeIds;
  ArrayView<const t_offset_t> offsets;
  ArrayView<t_input_t> e_wgts;
  graph->getEdgeList(edgeIds, offsets, e_wgts);

  std::vector<t_scalar_t> edge_weights;
  int numWeightPerEdge = graph->getNumWeightsPerEdge();

  if (numWeightPerEdge > 0) {
    edge_weights =  std::vector<t_scalar_t>(localNumEdges);
    for (t_lno_t i = 0; i < localNumEdges; ++i) {
      edge_weights[i] = e_wgts[0][i];
    }
  }

  //create a zoltan dictionary to get the parts of the vertices
  //at the other end of edges
  std::vector<part_t> e_parts(localNumEdges);
#ifdef HAVE_ZOLTAN2_MPI
  if (comm->getSize() > 1) {
    const bool bUseLocalIDs = false;  // Local IDs not needed
    typedef Zoltan2_Directory_Simple<t_gno_t,t_lno_t,part_t> directory_t;
    int debug_level = 0;
    const RCP<const Comm<int> > rcp_comm(comm,false);
    directory_t directory(rcp_comm, bUseLocalIDs, debug_level);
    directory.update(localNumVertices, &Ids[0], NULL, &parts[0],
      NULL, directory_t::Update_Mode::Replace);
    directory.find(localNumEdges, &edgeIds[0], NULL, &e_parts[0],
      NULL, NULL, false);
  } else
#endif
  {

    /*
    std::cout << "localNumVertices:" << localNumVertices
              << " np:" << np
              << " globalNumVertices:" << graph->getGlobalNumVertices()
              << " localNumEdges:" << localNumEdges << std::endl;
              */

    for (t_lno_t i = 0; i < localNumEdges; ++i) {
      t_gno_t ei = edgeIds[i];
      part_t p = parts[ei];
      e_parts[i] = p;
    }

    //get the vertices in each part in my part.
    std::vector<t_lno_t> part_begins(np, -1);
    std::vector<t_lno_t> part_nexts(localNumVertices, -1);

    //cluster vertices according to their parts.
    //create local part graph.
    for (t_lno_t i = 0; i < localNumVertices; ++i) {
      part_t ap = parts[i];
      part_nexts[i] = part_begins[ap];
      part_begins[ap] = i;
    }


    g_part_xadj = ArrayRCP<part_t>(np + 1);
    g_part_adj = ArrayRCP<part_t>(localNumEdges);
    g_part_ew = ArrayRCP<t_scalar_t>(localNumEdges);
    part_t nindex = 0;
    g_part_xadj[0] = 0;
    std::vector<part_t> part_neighbors(np);
    std::vector<t_scalar_t> part_neighbor_weights(np, 0);
    std::vector<t_scalar_t> part_neighbor_weights_ordered(np);

    //coarsen for all vertices in my part in order with parts.
    for (t_lno_t i = 0; i < np; ++i) {
      part_t num_neighbor_parts = 0;
      t_lno_t v = part_begins[i];
      //get part i, and first vertex in this part v.
      while (v != -1) {
        //now get the neightbors of v.
        for (t_offset_t j = offsets[v]; j < offsets[v+1]; ++j) {
          //get the part of the second vertex.
          part_t ep = e_parts[j];

          t_scalar_t ew = 1;
          if (numWeightPerEdge > 0) {
            ew = edge_weights[j];
          }

//          std::cout << "part:" << i << " v:" << v
//            << " part2:" << ep  << " v2:" << edgeIds[j]
//            << " w:" << ew << std::endl;

          if (part_neighbor_weights[ep] < 0.00001) {
            part_neighbors[num_neighbor_parts++] = ep;
          }

          part_neighbor_weights[ep] += ew;
        }

        v = part_nexts[v];
      }


      //now get the part list.
      for (t_lno_t j = 0; j < num_neighbor_parts; ++j) {
        part_t neighbor_part = part_neighbors[j];
        g_part_adj[nindex] = neighbor_part;
        g_part_ew[nindex++] = part_neighbor_weights[neighbor_part];
        part_neighbor_weights[neighbor_part] = 0;
      }
      g_part_xadj[i + 1] = nindex;
    }
  }

#ifdef HAVE_ZOLTAN2_MPI
  if(comm->getSize() > 1) { // Otherwise it's already handled above and done

    // struct for directory data - note more extensive comments in
    // Zoltan2_GraphMetricsUtility.hpp which I didn't want to duplicate
    // because this is in progress. Similar concept there.
    struct part_info {
      part_info() : weight(0) {
      }

      const part_info & operator+=(const part_info & src) {
        this->weight += src.weight;
        return *this;
      }

      bool operator>(const part_info & src) {
        return (destination_part > src.destination_part);
      }

      bool operator==(const part_info & src) {
        return (destination_part == src.destination_part);
      }

      part_t destination_part;
      scalar_t weight;
    };

    bool bUseLocalIDs = false;
    const int debug_level = 0;
    typedef Zoltan2_Directory_Vector<part_t,int,std::vector<part_info>>
        directory_t;
    const RCP<const Comm<int> > rcp_comm(comm,false);
    directory_t directory(rcp_comm, bUseLocalIDs, debug_level);
    std::vector<part_t> part_data;
    std::vector<std::vector<part_info>> user_data;

    envConst->timerStart(MACRO_TIMERS, "GRAPHCREATE Coarsen");
    {
      //get the vertices in each part in my part.
      std::vector<t_lno_t> part_begins(np, -1);
      std::vector<t_lno_t> part_nexts(localNumVertices, -1);

      //cluster vertices according to their parts.
      //create local part graph.
      for (t_lno_t i = 0; i < localNumVertices; ++i) {
        part_t ap = parts[i];
        part_nexts[i] = part_begins[ap];
        part_begins[ap] = i;
      }

      std::vector<part_t> part_neighbors(np);
      std::vector<t_scalar_t> part_neighbor_weights(np, 0);
      std::vector<t_scalar_t> part_neighbor_weights_ordered(np);

      //coarsen for all vertices in my part in order with parts.
      for (t_lno_t i = 0; i < np; ++i) {
        part_t num_neighbor_parts = 0;
        t_lno_t v = part_begins[i];

        //get part i, and first vertex in this part v.
        while (v != -1) {
          //now get the neightbors of v.
          for (t_offset_t j = offsets[v]; j < offsets[v+1]; ++j) {
            //get the part of the second vertex.
            part_t ep = e_parts[j];

            t_scalar_t ew = 1;
            if (numWeightPerEdge > 0) {
              ew = edge_weights[j];
            }

            //add it to my local part neighbors for part i.
            if (part_neighbor_weights[ep] < 0.00001) {
              part_neighbors[num_neighbor_parts++] = ep;
            }

            part_neighbor_weights[ep] += ew;
          }
          
          v = part_nexts[v];
        }

        //now get the part list.
        for (t_lno_t j = 0; j < num_neighbor_parts; ++j) {
          part_t neighbor_part = part_neighbors[j];
          part_neighbor_weights_ordered[j] =
            part_neighbor_weights[neighbor_part];
          part_neighbor_weights[neighbor_part] = 0;
        }

        //insert it to tpetra crsmatrix.
        if (num_neighbor_parts > 0) {
          part_data.push_back(i); // TODO: optimize to avoid push_back
          std::vector<part_info> new_user_data(num_neighbor_parts);
          for(int q = 0; q < num_neighbor_parts; ++q) {
            part_info & info = new_user_data[q];
            info.weight = part_neighbor_weights_ordered[q];
            info.destination_part = part_neighbors[q];
          }
          // TODO: optimize to avoid push_back
          user_data.push_back(new_user_data);
        }
      }
    }
    envConst->timerStop(MACRO_TIMERS, "GRAPHCREATE Coarsen");

    std::vector<part_t> part_indices(np);

    for (part_t i = 0; i < np; ++i) part_indices[i] = i;

    envConst->timerStart(MACRO_TIMERS, "GRAPHCREATE directory update");
    directory.update(part_data.size(), &part_data[0], NULL, &user_data[0],
      NULL, directory_t::Update_Mode::AggregateAdd);
    envConst->timerStop(MACRO_TIMERS, "GRAPHCREATE directory update");

    envConst->timerStart(MACRO_TIMERS, "GRAPHCREATE directory find");
    std::vector<std::vector<part_info>> find_user_data(part_indices.size());
    directory.find(find_user_data.size(), &part_indices[0], NULL,
      &find_user_data[0], NULL, NULL, false);
    envConst->timerStop(MACRO_TIMERS, "GRAPHCREATE directory find");

    // Now reconstruct the output data from the directory find data
    // This code was designed to reproduce the exact format of the original
    // setup for g_part_xadj, g_part_adj, and g_part_ew but before making any
    // further changes I wanted to verify if this formatting should be
    // preserved or potentially changed further.

    // first thing is get the total number of elements
    int get_total_length = 0;
    for (size_t n = 0; n < find_user_data.size(); ++n) {
      get_total_length += find_user_data[n].size();
    }

    // setup data space
    g_part_xadj = ArrayRCP<part_t>(np + 1);
    g_part_adj = ArrayRCP<part_t>(get_total_length);
    g_part_ew = ArrayRCP<t_scalar_t>(get_total_length);

    // loop through again and fill to match the original formatting
    int track_insert_index = 0;
    for(size_t n = 0; n < find_user_data.size(); ++n) {
      g_part_xadj[n] = track_insert_index;
      const std::vector<part_info> & user_data_vector = find_user_data[n];
      for(size_t q = 0; q < user_data_vector.size(); ++q) {
        const part_info & info = user_data_vector[q];
        g_part_adj[track_insert_index] = info.destination_part;
        g_part_ew[track_insert_index] = info.weight;
        ++track_insert_index;
      }
    }
    g_part_xadj[np] = get_total_length; // complete the series
  }
#endif // HAVE_ZOLTAN2_MPI
}


/*! \brief KmeansHeap Class, max heap, but holds the minimum values.
 */
template <class IT, class WT>
class KmeansHeap{
 
  IT heapSize;
  IT *indices;
  WT *values;
  WT epsilon;


public:
  ~KmeansHeap(){
    delete [] this->indices;
    delete [] this->values;
  }

  void setHeapsize(IT heapsize_){
    this->heapSize = heapsize_;
    this->indices = new IT[heapsize_];
    this->values = new WT[heapsize_];
    this->epsilon = std::numeric_limits<WT>::epsilon();
  }

  void addPoint(IT index, WT distance) {
    WT maxVal = this->values[0];
    //add only the distance is smaller than the maximum distance.
//    std::cout << "indeX:" << index
//      << " distance:" << distance
//      << " maxVal:" << maxVal << endl;
    if (distance >= maxVal)
      return;
    else {
      this->values[0] = distance;
      this->indices[0] = index;
      this->push_down(0);
    }
  }

  //heap push down operation
  void push_down(IT index_on_heap) {
    IT child_index1 = 2 * index_on_heap + 1;
    IT child_index2 = 2 * index_on_heap + 2;

    IT biggerIndex = -1;
    if(child_index1 < this->heapSize && child_index2 < this->heapSize) {

      if (this->values[child_index1] < this->values[child_index2]) {
        biggerIndex = child_index2;
      }
      else {
        biggerIndex = child_index1;
      }
    }
    else if(child_index1 < this->heapSize) {
      biggerIndex = child_index1;

    }
    else if(child_index2 < this->heapSize) {
      biggerIndex = child_index2;
    }
    if (biggerIndex >= 0 &&
        this->values[biggerIndex] > this->values[index_on_heap]) {
      WT tmpVal = this->values[biggerIndex];
      this->values[biggerIndex] = this->values[index_on_heap];
      this->values[index_on_heap] = tmpVal;

      IT tmpIndex = this->indices[biggerIndex];
      this->indices[biggerIndex] = this->indices[index_on_heap];
      this->indices[index_on_heap] = tmpIndex;
      this->push_down(biggerIndex);
    }
  }

  void initValues() {
    WT MAXVAL = std::numeric_limits<WT>::max();
    
    for(IT i = 0; i < this->heapSize; ++i) {
      this->values[i] = MAXVAL;
      this->indices[i] = -1;
    }
  }

  //returns the total distance to center in the cluster.
  WT getTotalDistance() {

    WT nc = 0;
    for(IT j = 0; j < this->heapSize; ++j) {
      nc += this->values[j];

//      std::cout << "index:" << this->indices[j]
//        << " distance:" << this->values[j] << endl;
    }
    return nc;
  }

  //returns the new center of the cluster.
  bool getNewCenters(WT *center, WT **coords, int dimension) {
    bool moved = false;
    
    for(int i = 0; i < dimension; ++i) {
      WT nc = 0;
      
      for(IT j = 0; j < this->heapSize; ++j) {
        IT k = this->indices[j];
//        std::cout << "i:" << i
//          << " dim:" << dimension
//          << " k:" << k
//          << " heapSize:" << heapSize << endl;        nc += coords[i][k];
        nc += coords[i][k];
      }

      nc /= this->heapSize;
      moved = (std::abs(center[i] - nc) > this->epsilon || moved );
      center[i] = nc;
    }
    return moved;
  }

  void copyCoordinates(IT *permutation) {
    for (IT i = 0; i < this->heapSize; ++i) {
      permutation[i] = this->indices[i];
    }
  }
};

/*! \brief
 */
template <class IT, class WT>
class KMeansCluster{

  int dimension;
  KmeansHeap<IT,WT> closestPoints;

public:
  
  WT *center;

  ~KMeansCluster() {
    delete [] center;
  }

  void setParams(int dimension_, int heapsize) {
    this->dimension = dimension_;
    this->center = new WT[dimension_];
    this->closestPoints.setHeapsize(heapsize);
  }

  void clearHeap(){
    this->closestPoints.initValues();
  }

  bool getNewCenters( WT **coords) {
    return this->closestPoints.getNewCenters(center, coords, dimension);
  }

  //returns the distance of the coordinate to the center.
  //also adds it to the heap.
  WT getDistance(IT index, WT **elementCoords) {
    WT distance = 0;

    for (int i = 0; i < this->dimension; ++i) {
      WT d = (center[i] - elementCoords[i][index]);
      distance += d * d;
    }
    distance = pow(distance, WT(1.0 / this->dimension));
    closestPoints.addPoint(index, distance);
    
    return distance;
  }

  WT getDistanceToCenter() {
    return closestPoints.getTotalDistance();
  }

  void copyCoordinates(IT *permutation) {
    closestPoints.copyCoordinates(permutation);
  }
};

/*! \brief KMeansAlgorithm Class that performs clustering of the
 *  coordinates, and returns the closest set of coordinates.
 *
 *  Useful to filter the processors, when there are more processors
 *  than needed.
 */
template <class IT, class WT>
class KMeansAlgorithm{

  int dim;
  IT numElements;
  WT **elementCoords;
  IT numClusters;
  IT required_elements;
  KMeansCluster<IT,WT> *clusters;
  WT *maxCoordinates;
  WT *minCoordinates;
public:
  
  ~KMeansAlgorithm(){
    delete [] clusters;
    delete [] maxCoordinates;
    delete [] minCoordinates;
  }

  /*! \brief KMeansAlgorithm Constructor
   */
  KMeansAlgorithm(
      int dim_ ,
      IT numElements_,
      WT **elementCoords_,
      IT required_elements_):
        dim(dim_),
        numElements(numElements_),
        elementCoords(elementCoords_),
        numClusters((1 << dim_) + 1),
        required_elements(required_elements_) {
    this->clusters = new KMeansCluster<IT,WT>[this->numClusters];
    
    //set dimension and the number of required elements for all clusters.
    for (int i = 0; i < numClusters; ++i) {
      this->clusters[i].setParams(this->dim, this->required_elements);
    }

    this->maxCoordinates = new WT[this->dim];
    this->minCoordinates = new WT[this->dim];

    //obtain the min and max coordinates for each dimension.
    for (int j = 0; j < dim; ++j) {
      this->minCoordinates[j] = this->elementCoords[j][0];
      this->maxCoordinates[j] = this->elementCoords[j][0];
      
      for(IT i = 1; i < numElements; ++i) {
        WT t = this->elementCoords[j][i];
        if(t > this->maxCoordinates[j]){
          this->maxCoordinates[j] = t;
        }

        if (t < minCoordinates[j]) {
          this->minCoordinates[j] = t;
        }
      }
    }


    //assign initial cluster centers.
    for (int j = 0; j < dim; ++j) {
      int mod = (1 << (j+1));
      for (int i = 0; i < numClusters - 1; ++i) {
        WT c = 0;

        if ( (i % mod) < mod / 2) {
          c = this->maxCoordinates[j];
//          std::cout << "i:" << i << " j:" << j
//            << " setting max:" << c << endl;
        }
        else {
          c = this->minCoordinates[j];
        }

        this->clusters[i].center[j] = c;
      }
    }

    //last cluster center is placed to middle.
    for (int j = 0; j < dim; ++j) {
      this->clusters[numClusters - 1].center[j] =
        (this->maxCoordinates[j] + this->minCoordinates[j]) / 2;
    }


/*
  for (int i = 0; i < numClusters; ++i) {
//    std::cout << endl << "cluster:" << i << endl << "\t";

    for (int j = 0; j < dim; ++j) {
      std::cout << this->clusters[i].center[j] << " ";
    }
  }
*/
  }

  // Performs kmeans clustering of coordinates.
  void kmeans() {
    for (int it = 0; it < 10; ++it) {
//      std::cout << "it:" << it << endl;

      for (IT j = 0; j < this->numClusters; ++j) {
        this->clusters[j].clearHeap();
      }
      
      for (IT i = 0; i < this->numElements; ++i) {
//        std::cout << "i:" << i << " numEl:" << this->numElements << endl;
        
        for (IT j = 0; j < this->numClusters; ++j) {
//          std::cout << "j:" << j
//            << " numClusters:" << this->numClusters << endl;
          
          this->clusters[j].getDistance(i,this->elementCoords);
        }
      }
      
      bool moved = false;
      
      for (IT j = 0; j < this->numClusters; ++j) {
        moved =
          (this->clusters[j].getNewCenters(this->elementCoords) || moved );
      }
      if (!moved) {
        break;
      }
    }
  }

  // Finds the cluster in which the coordinates are the closest to each
  // other.
  void getMinDistanceCluster(IT *procPermutation) {

    WT minDistance = this->clusters[0].getDistanceToCenter();
    IT minCluster = 0;

//    std::cout << "j:" << 0 << " minDistance:" << minDistance
//      << " minTmpDistance:" << minDistance
//      << " minCluster:" << minCluster << endl;
    
    for (IT j = 1; j < this->numClusters; ++j) {
      WT minTmpDistance = this->clusters[j].getDistanceToCenter();

//      std::cout << "j:" << j << " minDistance:" << minDistance
//        << " minTmpDistance:" << minTmpDistance
//        << " minCluster:" << minCluster << endl;

      if (minTmpDistance < minDistance) {
        minDistance = minTmpDistance;
        minCluster = j;
      }
    }

//    std::cout << "minCluster:" << minCluster << endl;
    this->clusters[minCluster].copyCoordinates(procPermutation);
  }
};


#define MINOF(a,b) (((a)<(b))?(a):(b))

/*! \brief fillContinousArray function
 *
 *  \param arr   array to be filled in with values.
 *  \param arrSize the size of the array.
 *  \param val    the pointer to the value to be filled. if given NULL,
 *         the filling performs arr[i] = i.
 */
template <typename T>
void fillContinousArray(T *arr, size_t arrSize, T *val) {
  if (val == NULL) {

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < arrSize; ++i) {
      arr[i] = i;
    }

  }
  else {
    T v = *val;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < arrSize; ++i) {
//      std::cout << "writing to i:" << i << " arr:" << arrSize << endl;
      arr[i] = v;
    }
  }
}

/*! \brief CommunicationModel Base Class that performs mapping between
 *  the coordinate partitioning result.
 */
template <typename part_t, typename pcoord_t, typename node_t>
class CommunicationModel{

protected:
  double commCost;

public:

  // Number of processors and number of tasks
  part_t no_procs;
  part_t no_tasks;


  CommunicationModel(): commCost(),no_procs(0), no_tasks(0) {}
  CommunicationModel(part_t no_procs_, part_t no_tasks_):
    commCost(),
    no_procs(no_procs_),
    no_tasks(no_tasks_) {}

  virtual ~CommunicationModel() {}

  part_t getNProcs() const{
    return this->no_procs;
  }

  part_t getNTasks()const{
    return this->no_tasks;
  }

  void calculateCommunicationCost(
      part_t *task_to_proc,
      part_t *task_communication_xadj,
      part_t *task_communication_adj,
      pcoord_t *task_communication_edge_weight) {

    double totalCost = 0;

    part_t commCount = 0;
    for (part_t task = 0; task < this->no_tasks; ++task) {
      int assigned_proc = task_to_proc[task];

      part_t task_adj_begin = task_communication_xadj[task];
      part_t task_adj_end = task_communication_xadj[task + 1];

      commCount += task_adj_end - task_adj_begin;

      for (part_t task2 = task_adj_begin; task2 < task_adj_end; ++task2) {

        part_t neighborTask = task_communication_adj[task2];
        int neighborProc = task_to_proc[neighborTask];
        double distance = getProcDistance(assigned_proc, neighborProc);

        if (task_communication_edge_weight == NULL) {
          totalCost += distance ;
        }
        else {
          totalCost += distance * task_communication_edge_weight[task2];
        }
      }
    }

    // commCount
    this->commCost = totalCost;
  }

  double getCommunicationCostMetric() {
    return this->commCost;
  }

  virtual double getProcDistance(int procId1, int procId2) const = 0;

  /*! \brief Function is called whenever nprocs > no_task.
   *  Function returns only the subset of processors that are closest to
   *  each other.
   *
   *  \param proc_to_task_xadj holds the pointer to the task array
   *  \param proc_to_task_adj holds the indices of tasks wrt to
   *         proc_to_task_xadj array.
   *  \param task_to_proc holds the processors mapped to tasks.
   */
  virtual void getMapping(
      int myRank,
      const RCP<const Environment> &env,
      ArrayRCP <part_t> &proc_to_task_xadj,
      ArrayRCP <part_t> &proc_to_task_adj,
      ArrayRCP <part_t> &task_to_proc,
      const Teuchos::RCP <const Teuchos::Comm<int> > comm_
  ) const = 0;
};


/*! \brief CoordinateModelInput Class that performs mapping between the coordinate partitioning result and mpi ranks
 * base on the coordinate results and mpi physical coordinates.
 */
template <typename pcoord_t,  typename tcoord_t, typename part_t, typename node_t>
class CoordinateCommunicationModel :
  public CommunicationModel<part_t, pcoord_t, node_t> {
public:
  //private:
  
  // Dimension of the processors
  int proc_coord_dim;
  // Processor coordinates (allocated outside of the class)
  pcoord_t **proc_coords;
  // Dimension of the tasks coordinates.
  int task_coord_dim;
  // Task coordinates (allocated outside of the class)
  tcoord_t **task_coords;

  // TODO: Perhaps delete this and just reference the view size?
  // Need to check the handling of size -1 versus size 0
  int partArraySize;

  Kokkos::View<part_t *, Kokkos::HostSpace> kokkos_partNoArray;

  int *machine_extent;
  bool *machine_extent_wrap_around;
  const MachineRepresentation<pcoord_t,part_t> *machine;

  int num_ranks_per_node;
  bool divide_to_prime_first;

  //public:
  CoordinateCommunicationModel():
    CommunicationModel<part_t, pcoord_t, node_t>(),
    proc_coord_dim(0),
    proc_coords(0),
    task_coord_dim(0),
    task_coords(0),
    partArraySize(-1),
    machine_extent(NULL),
    machine_extent_wrap_around(NULL),
    machine(NULL),
    num_ranks_per_node(1),
    divide_to_prime_first(false){}

  virtual ~CoordinateCommunicationModel() {}

  /*! \brief Class Constructor:
   *  \param pcoord_dim_ the dimension of the processors
   *  \param pcoords_   the processor coordinates. allocated outside of
             the class.
   *  \param tcoord_dim_   dimension of the tasks coordinates.
   *  \param tcoords_   the task coordinates allocated outside of the
             class.
   *  \param no_procs_   the number of processors
   *  \param no_tasks_   the number of taks.
   */
  CoordinateCommunicationModel(
      int pcoord_dim_,
      pcoord_t **pcoords_,
      int tcoord_dim_,
      tcoord_t **tcoords_,
      part_t no_procs_,
      part_t no_tasks_,
      int *machine_extent_,
      bool *machine_extent_wrap_around_,
      const MachineRepresentation<pcoord_t,part_t> *machine_ = NULL
  ):
    CommunicationModel<part_t, pcoord_t, node_t>(no_procs_, no_tasks_),
    proc_coord_dim(pcoord_dim_), proc_coords(pcoords_),
    task_coord_dim(tcoord_dim_), task_coords(tcoords_),
    partArraySize(-1),
    machine_extent(machine_extent_),
    machine_extent_wrap_around(machine_extent_wrap_around_),
    machine(machine_),
    num_ranks_per_node(1),
    divide_to_prime_first(false) {
  }


  void setPartArraySize(int psize) {
    this->partArraySize = psize;
  }

  void setPartArray(Kokkos::View<part_t *, Kokkos::HostSpace> pNo) {
    this->kokkos_partNoArray = pNo;
  }

  /*! \brief Function is called whenever nprocs > no_task.
   * Function returns only the subset of processors that are closest to each other.
   *
   *  \param proc_permutation holds the indices of the processors that are chosen.
   *  \param nprocs the number of processors.
   *  \param ntasks the number of taks.
   */
  void getClosestSubset(part_t *proc_permutation, part_t nprocs, part_t ntasks) const{
    //currently returns a random subset.

    part_t minCoordDim = MINOF(this->task_coord_dim, this->proc_coord_dim);
    KMeansAlgorithm<part_t, pcoord_t > kma(
        minCoordDim, nprocs,
        this->proc_coords, ntasks);

    kma.kmeans();
    kma.getMinDistanceCluster(proc_permutation);

    for(int i = ntasks; i < nprocs; ++i) {
      proc_permutation[i] = -1;
    }
    /*
  //fill array.
  fillContinousArray<part_t>(proc_permutation, nprocs, NULL);
  '
  int _u_umpa_seed = 847449649;
  srand (time(NULL));
  
  int a = rand() % 1000 + 1;
  _u_umpa_seed -= a;
  
  //permute array randomly.
  update_visit_order(proc_permutation, nprocs,_u_umpa_seed, 1);
     */
  }

  // Temporary, necessary for random permutation.
  static part_t umpa_uRandom(part_t l, int &_u_umpa_seed)
  {
    int   a = 16807;
    int   m = 2147483647;
    int   q = 127773;
    int   r = 2836;
    int   lo, hi, test;
    double d;

    lo = _u_umpa_seed % q;
    hi = _u_umpa_seed / q;
    test = (a * lo) - (r * hi);
    if (test>0)
      _u_umpa_seed = test;
    else
      _u_umpa_seed = test + m;
    
    d = (double) ((double) _u_umpa_seed / (double) m);
    
    return (part_t) (d*(double)l);
  }

  virtual double getProcDistance(int procId1, int procId2) const {
    pcoord_t distance = 0;
    if (machine == NULL) {
      for (int i = 0 ; i < this->proc_coord_dim; ++i) {
        double d =
          std::abs(proc_coords[i][procId1] - proc_coords[i][procId2]);
        if (machine_extent_wrap_around && machine_extent_wrap_around[i]) {
          if (machine_extent[i] - d < d) {
            d = machine_extent[i] - d;
          }
        }
        distance += d;
      }
    }
    else {
      this->machine->getHopCount(procId1, procId2, distance);
    }

    return distance;
  }


  // Temporary, does random permutation.
  void update_visit_order(part_t* visitOrder, part_t n,
                          int &_u_umpa_seed, part_t rndm) {
    part_t *a = visitOrder;

    if (rndm) {
      part_t i, u, v, tmp;

      if (n <= 4)
        return;

//    srand ( time(NULL) );
//    _u_umpa_seed = _u_umpa_seed1 - (rand()%100);
      
      for (i = 0; i < n; i += 16) {
        u = umpa_uRandom(n-4, _u_umpa_seed);
        v = umpa_uRandom(n-4, _u_umpa_seed);

        // FIXME (mfh 30 Sep 2015) This requires including
        // Zoltan2_AlgMultiJagged.hpp.

        ZOLTAN2_ALGMULTIJAGGED_SWAP(a[v], a[u], tmp);
        ZOLTAN2_ALGMULTIJAGGED_SWAP(a[v + 1], a[u + 1], tmp);
        ZOLTAN2_ALGMULTIJAGGED_SWAP(a[v + 2], a[u + 2], tmp);
        ZOLTAN2_ALGMULTIJAGGED_SWAP(a[v + 3], a[u + 3], tmp);
      }
    }
    else {
      part_t i, end = n / 4;

      for (i = 1; i < end; i++) {
        part_t j = umpa_uRandom(n - i, _u_umpa_seed);
        part_t t = a[j];
        a[j] = a[n - i];
        a[n - i] = t;
      }
    }
//    PermuteInPlace(visitOrder, n);
  }



  /*! \brief Function is called whenever nprocs > no_task.
   * Function returns only the subset of processors that are closest to
     each other.
   *  \param proc_to_task_xadj holds the pointer to the task array
   *  \param proc_to_task_xadj holds the indices of tasks wrt to
             proc_to_task_xadj array.
   *  \param task_to_proc holds the processors mapped to tasks.
   */
  virtual void getMapping(
      int myRank,
      const RCP<const Environment> &env,
      ArrayRCP <part_t> &rcp_proc_to_task_xadj,
      ArrayRCP <part_t> &rcp_proc_to_task_adj,
      ArrayRCP <part_t> &rcp_task_to_proc,
      const Teuchos::RCP <const Teuchos::Comm<int> > comm_
  ) const {

    rcp_proc_to_task_xadj = ArrayRCP <part_t>(this->no_procs + 1);
    rcp_proc_to_task_adj = ArrayRCP <part_t>(this->no_tasks);
    rcp_task_to_proc = ArrayRCP <part_t>(this->no_tasks);

    // Holds the pointer to the task array
    part_t *proc_to_task_xadj = rcp_proc_to_task_xadj.getRawPtr();

    // Holds the indices of task wrt to proc_to_task_xadj array.
    part_t *proc_to_task_adj = rcp_proc_to_task_adj.getRawPtr();

    // Holds the processors mapped to tasks.
    part_t *task_to_proc = rcp_task_to_proc.getRawPtr();

    part_t invalid = 0;
    fillContinousArray<part_t> (proc_to_task_xadj, this->no_procs + 1, &invalid);

    // Obtain the number of parts that should be divided.
    part_t num_parts = MINOF(this->no_procs, this->no_tasks);
    
    // Obtain the min coordinate dim.
    //  No more want to do min coord dim. If machine dimension > task_dim,
    //  we end up with a long line.
//    part_t minCoordDim = MINOF(this->task_coord_dim, this->proc_coord_dim);

    int recursion_depth = partArraySize;

//    if (partArraySize < minCoordDim)
//      recursion_depth = minCoordDim;
    if (partArraySize == -1) {

      if (divide_to_prime_first) {
        // It is difficult to estimate the number of steps in this case
        // as each branch will have different depth.
        // The worst case happens when all prime factors are 3s.
        // P = 3^n, n recursion depth will divide parts to 2x and x
        // and n recursion depth with divide 2x into x and x.
        // Set it to upperbound here.
        // We could calculate the exact value here as well, but the
        // partitioning algorithm skips further ones anyways.        recursion_depth = log(float(this->no_procs)) / log(2.0) * 2 + 1;
      }
      else {
        recursion_depth = log(float(this->no_procs)) / log(2.0) + 1;
      }
    }

    // Number of permutations for tasks and processors
    int taskPerm = 1;
    int procPerm = 1;

    // Get number of different permutations for task dimension ordering
    if (this->task_coord_dim <= 8)
      taskPerm = z2Fact<int>(this->task_coord_dim);
    // Prevent overflow
    else
      taskPerm = z2Fact<int>(8);

    // Get number of different permutations for proc dimension ordering
    if (this->proc_coord_dim <= 8)
      procPerm = z2Fact<int>(this->proc_coord_dim);
    // Prevent overflow
    else
      procPerm = z2Fact<int>(8);


    // Total number of permutations (both task and proc permuted)
    int permutations =  taskPerm * procPerm;

    // Add permutations where we divide the processors with longest
    // dimension but permute tasks.
    permutations += taskPerm;

    // Add permutations where we divide the tasks with longest
    // dimension but permute procs.
    permutations += procPerm;

    // Add permutation with both tasks and procs divided by longest
    // dimension
    permutations += 1;
    //add one also that partitions based the longest dimension.

    // Holds the pointers to proc_adjList
    part_t *proc_xadj = new part_t[num_parts+1];

    // Holds the processors in parts according to the result of
    // partitioning algorithm.
    // The processors assigned to part x is at
    // proc_adjList[ proc_xadj[x] : proc_xadj[x + 1] ]
    part_t *proc_adjList = new part_t[this->no_procs];


    part_t used_num_procs = this->no_procs;
    if (this->no_procs > this->no_tasks) {
      // Obtain the subset of the processors that are closest to each
      // other.
      this->getClosestSubset(proc_adjList, this->no_procs,
                             this->no_tasks);
      used_num_procs = this->no_tasks;
    }
    else {
      fillContinousArray<part_t>(proc_adjList,this->no_procs, NULL);
    }

    // Index of the permutation
    int myPermutation = myRank % permutations;
    bool task_partition_along_longest_dim = false;
    bool proc_partition_along_longest_dim = false;


    int myProcPerm = 0;
    int myTaskPerm = 0;

    if (myPermutation == 0) {
      task_partition_along_longest_dim = true;
      proc_partition_along_longest_dim = true;
    }
    else {
      --myPermutation;
      if (myPermutation < taskPerm) {
        proc_partition_along_longest_dim = true;
        // Index of the task permutation
        myTaskPerm  = myPermutation;
      }
      else {
        myPermutation -= taskPerm;
        if (myPermutation < procPerm) {
          task_partition_along_longest_dim = true;
          // Index of the task permutation
          myProcPerm  = myPermutation;
        }
        else {
          myPermutation -= procPerm;
          // Index of the proc permutation
         myProcPerm = myPermutation % procPerm;
          // Index of the task permutation
          myTaskPerm = myPermutation / procPerm;
        }
      }
    }

/*
    if (task_partition_along_longest_dim &&
        proc_partition_along_longest_dim) {
      std::cout <<"me:" << myRank << " task:longest proc:longest"
        << " numPerms:" << permutations << std::endl;
    }
    else if (proc_partition_along_longest_dim) {
      std::cout <<"me:" << myRank << " task:" <<  myTaskPerm
        << " proc:longest" << " numPerms:" << permutations << std::endl;
    }
    else if (task_partition_along_longest_dim) {
      std::cout <<"me:" << myRank << " task: longest" << " proc:"
        <<  myProcPerm  << " numPerms:" << permutations << std::endl;
    }
    else {
      std::cout <<"me:" << myRank << " task:" <<  myTaskPerm << " proc:"
        <<  myProcPerm  << " numPerms:" << permutations << std::endl;
    }
*/

    int *permutation = new int[(this->proc_coord_dim > this->task_coord_dim)
        ? this->proc_coord_dim : this->task_coord_dim];

    // Get the permutation order from the proc permutation index.
    if (this->proc_coord_dim <= 8)
      ithPermutation<int>(this->proc_coord_dim, myProcPerm, permutation);
    else
      ithPermutation<int>(8, myProcPerm, permutation);

/*
    // Reorder the coordinate dimensions.
    pcoord_t **pcoords = allocMemory<pcoord_t *>(this->proc_coord_dim);
    for (int i = 0; i < this->proc_coord_dim; ++i) {
      pcoords[i] = this->proc_coords[permutation[i]];
//      std::cout << permutation[i] << " ";
    }
*/

    int procdim  = this->proc_coord_dim;
    pcoord_t **pcoords = this->proc_coords;

/*
    int procdim  = this->proc_coord_dim;
    procdim  = 6;
    //reorder the coordinate dimensions.
    pcoord_t **pcoords = allocMemory<pcoord_t *>(procdim);
    for (int i = 0; i < procdim; ++i) {
      pcoords[i] = new pcoord_t[used_num_procs] ;
//      this->proc_coords[permutation[i]];
    }

    for (int k = 0; k < used_num_procs ; k++) {
      pcoords[0][k] = (int (this->proc_coords[0][k]) / 2) * 64;
      pcoords[3][k] = (int (this->proc_coords[0][k]) % 2) * 8 ;

      pcoords[1][k] = (int (this->proc_coords[1][k])  / 2) * 8 * 2400;
      pcoords[4][k] = (int (this->proc_coords[1][k])  % 2) * 8;
      pcoords[2][k] = ((int (this->proc_coords[2][k])) / 8) * 160;
      pcoords[5][k] = ((int (this->proc_coords[2][k])) % 8) * 5;

      //if (this->proc_coords[0][k] == 40 &&
      //    this->proc_coords[1][k] == 8 &&
      //    this->proc_coords[2][k] == 48) {
      if (this->proc_coords[0][k] == 5 &&
          this->proc_coords[1][k] == 0 &&
          this->proc_coords[2][k] == 10) {
        std::cout << "pcoords[0][k]:" << pcoords[0][k]
          << "pcoords[1][k]:" << pcoords[1][k]
          << "pcoords[2][k]:" << pcoords[2][k]
          << "pcoords[3][k]:" << pcoords[3][k]
          << "pcoords[4][k]:" << pcoords[4][k]
          << "pcoords[5][k]:" << pcoords[5][k] << std::endl;
      }
      else if (pcoords[0][k] == 64 &&
               pcoords[1][k] == 0 &&
               pcoords[2][k] == 160 &&
               pcoords[3][k] == 16 &&
               pcoords[4][k] == 0 &&
               pcoords[5][k] == 10) {
        std::cout << "this->proc_coords[0][k]:" << this->proc_coords[0][k]
          << "this->proc_coords[1][k]:" << this->proc_coords[1][k]
          << "this->proc_coords[2][k]:" << this->proc_coords[2][k]
          << std::endl;
      }
    }
*/

//    if (partNoArray == NULL)
//      std::cout << "partNoArray is null" << std::endl;
//    std::cout << "recursion_depth:" << recursion_depth
//      << " partArraySize:" << partArraySize << std::endl;

    // Optimization for Dragonfly Networks, First Level of partitioning
    // is imbalanced to ensure procs are divided by first RCA
    // coord (a.k.a. group).
    part_t num_group_count = 1;
    part_t *group_count = NULL;

    if (machine != NULL)
      num_group_count = machine->getNumUniqueGroups();
 
    if (num_group_count > 1) {
      group_count = new part_t[num_group_count];
      memset(group_count, 0, sizeof(part_t) * num_group_count);
 
      machine->getGroupCount(group_count);
    }

    // Do the partitioning and renumber the parts.
    env->timerStart(MACRO_TIMERS, "Mapping - Proc Partitioning");
    // Partitioning of Processors
    AlgMJ<pcoord_t, part_t, part_t, part_t, node_t> mj_partitioner;

    typedef typename node_t::device_type device_t;
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<pcoord_t**, Kokkos::LayoutLeft, device_t>
      kokkos_pcoords("pcoords", this->no_procs, procdim);
    auto host_kokkos_pcoords = Kokkos::create_mirror_view(kokkos_pcoords);
    for(int i = 0; i < procdim; ++i) {
      for(int j = 0; j < this->no_procs; ++j) {
        host_kokkos_pcoords(j,i) = pcoords[i][j];
      }
    }
    Kokkos::deep_copy(kokkos_pcoords, host_kokkos_pcoords);

    Kokkos::View<part_t*, device_t> initial_selected_coords_output_permutation_pcoords(
      "initial_selected_coords_output_permutation_pcoords", this->no_procs);
    typename Kokkos::View<part_t*, device_t>::HostMirror
      host_initial_selected_coords_output_permutation_pcoords =
        Kokkos::create_mirror_view(initial_selected_coords_output_permutation_pcoords);
    for(int n = 0; n < this->no_procs; ++n) {
      host_initial_selected_coords_output_permutation_pcoords(n) =
        proc_adjList[n];
    }
    Kokkos::deep_copy(initial_selected_coords_output_permutation_pcoords,
      host_initial_selected_coords_output_permutation_pcoords);

    // Note num_group_count = 1 when group_count = NULL - perhaps could change
    Kokkos::View<part_t *, Kokkos::HostSpace> kokkos_group_count(
      "kokkos_group_count", group_count ? num_group_count : 0);
    if(group_count) {
      for(int n = 0; n < num_group_count; ++n) {
        kokkos_group_count(n) = group_count[n];
      }
    }

    mj_partitioner.sequential_task_partitioning(
        env,
        this->no_procs,
        used_num_procs,
        num_parts,
        procdim,
        //minCoordDim,
        kokkos_pcoords,
        initial_selected_coords_output_permutation_pcoords,
        proc_xadj,
        recursion_depth,
        kokkos_partNoArray,
        proc_partition_along_longest_dim, // false
        num_ranks_per_node,
        divide_to_prime_first,
        num_group_count,
        kokkos_group_count);
    env->timerStop(MACRO_TIMERS, "Mapping - Proc Partitioning");
//    comm_->barrier();
//    std::cout << "mj_partitioner.for procs over" << std::endl;
//    freeArray<pcoord_t *>(pcoords);

    Kokkos::deep_copy(host_initial_selected_coords_output_permutation_pcoords,
      initial_selected_coords_output_permutation_pcoords);
    for(int n = 0; n < this->no_procs; ++n) {
      proc_adjList[n] =
        host_initial_selected_coords_output_permutation_pcoords(n);
    }

    part_t *task_xadj = new part_t[num_parts + 1];
    part_t *task_adjList = new part_t[this->no_tasks];

    // Fill task_adjList st: task_adjList[i] <- i.
    fillContinousArray<part_t>(task_adjList,this->no_tasks, NULL);

    // Get the permutation order from the task permutation index.
    if (this->task_coord_dim <= 8)
      ithPermutation<int>(this->task_coord_dim, myTaskPerm, permutation);
    else
      ithPermutation<int>(8, myTaskPerm, permutation);

    // Reorder task coordinate dimensions.
    tcoord_t **tcoords = new tcoord_t*[this->task_coord_dim];
    for(int i = 0; i < this->task_coord_dim; ++i) {
      tcoords[i] = this->task_coords[permutation[i]];
    }

    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<tcoord_t**, Kokkos::LayoutLeft, device_t>
      kokkos_tcoords("tcoords", this->no_tasks, this->task_coord_dim);
    auto host_kokkos_tcoords = Kokkos::create_mirror_view(kokkos_tcoords);
    for(int i = 0; i < this->task_coord_dim; ++i) {
      for(int j = 0; j < this->no_tasks; ++j) {
        host_kokkos_tcoords(j,i) = tcoords[i][j];
      }
    }
    Kokkos::deep_copy(kokkos_tcoords, host_kokkos_tcoords);

    env->timerStart(MACRO_TIMERS, "Mapping - Task Partitioning");

    Kokkos::View<part_t*, device_t> initial_selected_coords_output_permutation_tcoords(
      "initial_selected_coords_output_permutation_tcoords", this->no_tasks);
    typename Kokkos::View<part_t*, device_t>::HostMirror
      host_initial_selected_coords_output_permutation_tcoords =
        Kokkos::create_mirror_view(initial_selected_coords_output_permutation_tcoords);
    for(int n = 0; n < this->no_tasks; ++n) {
      host_initial_selected_coords_output_permutation_tcoords(n) =
        task_adjList[n];
    }
    Kokkos::deep_copy(initial_selected_coords_output_permutation_tcoords,
      host_initial_selected_coords_output_permutation_tcoords);

    //partitioning of tasks
    mj_partitioner.sequential_task_partitioning(
        env,
        this->no_tasks,
        this->no_tasks,
        num_parts,
        this->task_coord_dim,
        //minCoordDim,
        kokkos_tcoords,
        initial_selected_coords_output_permutation_tcoords,
        task_xadj,
        recursion_depth,
        kokkos_partNoArray,
        task_partition_along_longest_dim,
        num_ranks_per_node,
        divide_to_prime_first,
        num_group_count,
        kokkos_group_count);
    env->timerStop(MACRO_TIMERS, "Mapping - Task Partitioning");

    Kokkos::deep_copy(host_initial_selected_coords_output_permutation_tcoords,
      initial_selected_coords_output_permutation_tcoords);
    for(int n = 0; n < this->no_tasks; ++n) {
      task_adjList[n] =
        host_initial_selected_coords_output_permutation_tcoords(n);
    }

//    std::cout << "myrank:" << myRank << std::endl;
//    comm_->barrier();
//    std::cout << "mj_partitioner.sequential_task_partitioning over"
//      << std::endl;

    delete [] tcoords;
    delete [] permutation;


    //filling proc_to_task_xadj, proc_to_task_adj, task_to_proc arrays.
    for(part_t i = 0; i < num_parts; ++i) {

      part_t proc_index_begin = proc_xadj[i];
      part_t task_begin_index = task_xadj[i];
      part_t proc_index_end = proc_xadj[i + 1];
      part_t task_end_index = task_xadj[i + 1];


      if(proc_index_end - proc_index_begin != 1) {
        std::cerr << "Error at partitioning of processors" << std::endl;
        std::cerr << "PART:" << i << " is assigned to "
          << proc_index_end - proc_index_begin << " processors."
          << std::endl;
        std::terminate();
      }
      part_t assigned_proc = proc_adjList[proc_index_begin];
      proc_to_task_xadj[assigned_proc] = task_end_index - task_begin_index;
    }

    //holds the pointer to the task array
    //convert proc_to_task_xadj to CSR index array
    part_t *proc_to_task_xadj_work = new part_t[this->no_procs];
    part_t sum = 0;
    for(part_t i = 0; i < this->no_procs; ++i) {
      part_t tmp = proc_to_task_xadj[i];
      proc_to_task_xadj[i] = sum;
      sum += tmp;
      proc_to_task_xadj_work[i] = sum;
    }
    proc_to_task_xadj[this->no_procs] = sum;

    for(part_t i = 0; i < num_parts; ++i){

      part_t proc_index_begin = proc_xadj[i];
      part_t task_begin_index = task_xadj[i];
      part_t task_end_index = task_xadj[i + 1];

      part_t assigned_proc = proc_adjList[proc_index_begin];

      for (part_t j = task_begin_index; j < task_end_index; ++j) {
        part_t taskId = task_adjList[j];

        task_to_proc[taskId] = assigned_proc;

        proc_to_task_adj [--proc_to_task_xadj_work[assigned_proc]] = taskId;
      }
    }

/*
    if (myPermutation == 0) {
      std::ofstream gnuPlotCode ("mymapping.out", std::ofstream::out);

      for (part_t i = 0; i < num_parts; ++i) {

        part_t proc_index_begin = proc_xadj[i];
        part_t proc_index_end = proc_xadj[i + 1];

        if (proc_index_end - proc_index_begin != 1) {
          std::cerr << "Error at partitioning of processors" << std::endl;
          std::cerr << "PART:" << i << " is assigned to "
            << proc_index_end - proc_index_begin << " processors."
            << std::endl;
          exit(1);
        }

        part_t assigned_proc = proc_adjList[proc_index_begin];
        gnuPlotCode << "Rank:" << i << " "
          << this->proc_coords[0][assigned_proc] << " "
          << this->proc_coords[1][assigned_proc] << " "
          << this->proc_coords[2][assigned_proc] << " "
          << pcoords[0][assigned_proc] << " "
          << pcoords[1][assigned_proc] << " "
          << pcoords[2][assigned_proc] << " "
          << pcoords[3][assigned_proc] << std::endl;
      }

      gnuPlotCode << "Machine Extent:" << std::endl;
      //filling proc_to_task_xadj, proc_to_task_adj, task_to_proc arrays.
      for (part_t i = 0; i < num_parts; ++i) {

        part_t proc_index_begin = proc_xadj[i];
        part_t proc_index_end = proc_xadj[i + 1];

        if (proc_index_end - proc_index_begin != 1) {
          std::cerr << "Error at partitioning of processors" << std::endl;
          std::cerr << "PART:" << i << " is assigned to "
            << proc_index_end - proc_index_begin << " processors."
            << std::endl;
          exit(1);
        }

        part_t assigned_proc = proc_adjList[proc_index_begin];
        gnuPlotCode << "Rank:" << i << " "
        << this->proc_coords[0][assigned_proc] << " "
        << this->proc_coords[1][assigned_proc] << " "
        << this->proc_coords[2][assigned_proc] << std::endl;
      }
      gnuPlotCode.close();
    }
*/

    delete [] proc_to_task_xadj_work;
    delete [] task_xadj;
    delete [] task_adjList;
    delete [] proc_xadj;
    delete [] proc_adjList;
  }

};

template <typename Adapter, typename part_t>
class CoordinateTaskMapper:public PartitionMapping<Adapter>{
protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  typedef typename Adapter::scalar_t pcoord_t;
  typedef typename Adapter::scalar_t tcoord_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::lno_t lno_t;

#ifdef KOKKOS_ENABLE_CUDA
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<
      Kokkos::Cuda, Kokkos::CudaSpace>  node_t;
#else
  typedef typename Adapter::node_t node_t;
#endif

#endif

//  RCP<const Environment> env;
 
  // Holds the pointer to the task array
  ArrayRCP<part_t> proc_to_task_xadj;
//    = allocMemory<part_t> (this->no_procs + 1);
  
  // Holds the indices of tasks wrt to proc_to_task_xadj array.
  ArrayRCP<part_t> proc_to_task_adj;
//    = allocMemory<part_t>(this->no_tasks);
  
  // Holds the processors mapped to tasks.
  ArrayRCP<part_t> task_to_proc;
//    = allocMemory<part_t>(this->no_procs);

  // Holds the processors mapped to tasks.
  ArrayRCP<part_t> local_task_to_rank;
//    = allocMemory<part_t>(this->no_procs);

  bool isOwnerofModel;
  CoordinateCommunicationModel<pcoord_t,tcoord_t,part_t,node_t> *proc_task_comm;
  part_t nprocs;
  part_t ntasks;
  ArrayRCP<part_t> task_communication_xadj;
  ArrayRCP<part_t> task_communication_adj;
  ArrayRCP<scalar_t> task_communication_edge_weight;


  /*! \brief doMapping function, calls getMapping function of
   *  communicationModel object.
   */
  void doMapping(int myRank,
                 const Teuchos::RCP<const Teuchos::Comm<int> > comm_) {

    if (this->proc_task_comm) {
      this->proc_task_comm->getMapping(
          myRank,
          this->env,
          // Holds the pointer to the task array
          this->proc_to_task_xadj,
          // Holds the indices of task wrt to proc_to_task_xadj array
          this->proc_to_task_adj,
          // Holds the processors mapped to tasks
          this->task_to_proc,
          comm_
      );
    }
    else {
      std::cerr << "communicationModel is not specified in the Mapper"
        << std::endl;
      exit(1);
    }
  }


  /*! \brief creates and returns the subcommunicator for the
   *  processor group.
   */
  RCP<Comm<int> > create_subCommunicator() {
    int procDim = this->proc_task_comm->proc_coord_dim;
    int taskDim = this->proc_task_comm->task_coord_dim;

    // Number of permutations for tasks and processors
    int taskPerm = 1;
    int procPerm = 1;

    // Get number of different permutations for task dimension ordering
    if (taskDim <= 8)
      taskPerm = z2Fact<int>(taskDim);
    // Prevent overflow
    else
      taskPerm = z2Fact<int>(8);

    // Get number of different permutations for proc dimension ordering
    if (procDim <= 8)
      procPerm = z2Fact<int>(procDim);
    // Prevent overflow
    else
      procPerm = z2Fact<int>(8);

    // Total number of permutations
    int idealGroupSize =  taskPerm * procPerm;

    // For the one that does longest dimension partitioning.
    idealGroupSize += taskPerm + procPerm + 1;

    int myRank = this->comm->getRank();
    int commSize = this->comm->getSize();

    int myGroupIndex = myRank / idealGroupSize;

    int prevGroupBegin = (myGroupIndex - 1)* idealGroupSize;
    if (prevGroupBegin < 0) prevGroupBegin = 0;
    int myGroupBegin = myGroupIndex * idealGroupSize;
    int myGroupEnd = (myGroupIndex + 1) * idealGroupSize;
    int nextGroupEnd = (myGroupIndex + 2)* idealGroupSize;

    if (myGroupEnd > commSize) {
      myGroupBegin = prevGroupBegin;
      myGroupEnd = commSize;
    }
    if (nextGroupEnd > commSize) {
      myGroupEnd = commSize;
    }
    int myGroupSize = myGroupEnd - myGroupBegin;

    part_t *myGroup = new part_t[myGroupSize];
    for (int i = 0; i < myGroupSize; ++i) {
      myGroup[i] = myGroupBegin + i;
    }
//    std::cout << "me:" << myRank << " myGroupBegin:" << myGroupBegin
//      << " myGroupEnd:" << myGroupEnd << endl;

    ArrayView<const part_t> myGroupView(myGroup, myGroupSize);

    RCP<Comm<int> > subComm =
      this->comm->createSubcommunicator(myGroupView);
    delete [] myGroup;
    return subComm;
  }


  /*! \brief finds the lowest cost mapping and broadcasts solution to
   *  everyone.
   */
  void getBestMapping() {
    //create the sub group.
    RCP<Comm<int> > subComm = this->create_subCommunicator();
    //calculate cost.
    double myCost = this->proc_task_comm->getCommunicationCostMetric();
//    std::cout << "me:" << this->comm->getRank() << " myCost:"
//      << myCost << std::endl;
    double localCost[2], globalCost[2];

    localCost[0] = myCost;
    localCost[1] = double(subComm->getRank());

    globalCost[1] = globalCost[0] = std::numeric_limits<double>::max();
    Teuchos::Zoltan2_ReduceBestMapping<int,double> reduceBest;
    reduceAll<int, double>(*subComm, reduceBest,
        2, localCost, globalCost);

    int sender = int(globalCost[1]);

/*
    if ( this->comm->getRank() == 0) {
      std::cout << "me:" << localCost[1] <<
            " localcost:" << localCost[0]<<
            " bestcost:" << globalCost[0] <<
            " Sender:" << sender <<
            " procDim" << proc_task_comm->proc_coord_dim <<
            " taskDim:" << proc_task_comm->task_coord_dim << std::endl;
    }
*/

//    std::cout << "me:" << localCost[1] << " localcost:" << localCost[0]
//      << " bestcost:" << globalCost[0] << endl;
//    std::cout << "me:" << localCost[1] << " proc:" << globalCost[1]
//      << endl;
    broadcast(*subComm, sender, this->ntasks,
              this->task_to_proc.getRawPtr());
    broadcast(*subComm, sender, this->nprocs,
              this->proc_to_task_xadj.getRawPtr());
    broadcast(*subComm, sender, this->ntasks,
              this->proc_to_task_adj.getRawPtr());
  }

  //write mapping to gnuPlot code to visualize.
  void writeMapping() {
    std::ofstream gnuPlotCode("gnuPlot.plot", std::ofstream::out);

    int mindim = MINOF(proc_task_comm->proc_coord_dim,
                       proc_task_comm->task_coord_dim);
    std::string ss = "";
    for (part_t i = 0; i < this->nprocs; ++i) {

      std::string procFile = Teuchos::toString<int>(i) + "_mapping.txt";
      if (i == 0) {
        gnuPlotCode << "plot \"" << procFile << "\"\n";
      }
      else {
        gnuPlotCode << "replot \"" << procFile << "\"\n";
      }

      std::ofstream inpFile(procFile.c_str(), std::ofstream::out);

      std::string gnuPlotArrow = "set arrow from ";
      for (int j = 0; j <  mindim; ++j) {
        if (j == mindim - 1) {
          inpFile << proc_task_comm->proc_coords[j][i];
          gnuPlotArrow +=
            Teuchos::toString<float>(proc_task_comm->proc_coords[j][i]);

        }
        else {
          inpFile << proc_task_comm->proc_coords[j][i] << " ";
          gnuPlotArrow +=
            Teuchos::toString<float>(proc_task_comm->
                proc_coords[j][i]) + ",";
        }
      }
      gnuPlotArrow += " to ";

      inpFile << std::endl;
      ArrayView<part_t> a = this->getAssignedTasksForProc(i);
      
      for (int k = 0; k <  a.size(); ++k) {
        int j = a[k];
//      std::cout << "i:" << i << " j:"
        std::string gnuPlotArrow2 = gnuPlotArrow;
        for (int z = 0; z <  mindim; ++z) {
          if (z == mindim - 1) {

//          std::cout << "z:" << z << " j:" << j << " "
//            << proc_task_comm->task_coords[z][j] << endl;
            inpFile << proc_task_comm->task_coords[z][j];
            gnuPlotArrow2 +=
              Teuchos::toString<float>(proc_task_comm->task_coords[z][j]);
          }
          else {
            inpFile << proc_task_comm->task_coords[z][j] << " ";
            gnuPlotArrow2 +=
              Teuchos::toString<float>(proc_task_comm->
                  task_coords[z][j]) + ",";
          }
        }
        ss += gnuPlotArrow2 + "\n";
        inpFile << std::endl;
      }
      inpFile.close();
    }
    gnuPlotCode << ss;
    gnuPlotCode << "\nreplot\n pause -1 \n";
    gnuPlotCode.close();
  }

  //write mapping to gnuPlot code to visualize.
  void writeMapping2(int myRank) {
    std::string rankStr = Teuchos::toString<int>(myRank);
    std::string gnuPlots = "gnuPlot", extentionS = ".plot";
    std::string outF = gnuPlots + rankStr+ extentionS;
    std::ofstream gnuPlotCode(outF.c_str(), std::ofstream::out);

    CoordinateCommunicationModel<pcoord_t, tcoord_t, part_t, node_t>
      *tmpproc_task_comm =
        static_cast <CoordinateCommunicationModel<
          pcoord_t, tcoord_t, part_t, node_t> * > (
          proc_task_comm);

//    int mindim = MINOF(tmpproc_task_comm->proc_coord_dim,
//                       tmpproc_task_comm->task_coord_dim);
    int mindim = tmpproc_task_comm->proc_coord_dim;
    if (mindim != 3) {
      std::cerr << "Mapping Write is only good for 3 dim" << std::endl;
      return;
    }
    std::string ss = "";
    std::string procs = "";

    std::set < std::tuple<int,int,int,int,int,int> > my_arrows;
    for (part_t origin_rank = 0; origin_rank < this->nprocs; ++origin_rank) {
      ArrayView<part_t> a = this->getAssignedTasksForProc(origin_rank);
      if (a.size() == 0) {
        continue;
      }

      std::string gnuPlotArrow = "set arrow from ";
      for (int j = 0; j <  mindim; ++j) {
        if (j == mindim - 1) {
          gnuPlotArrow +=
            Teuchos::toString<float>(tmpproc_task_comm->
                proc_coords[j][origin_rank]);
          procs +=
            Teuchos::toString<float>(tmpproc_task_comm->
                proc_coords[j][origin_rank]);

        }
        else {
          gnuPlotArrow +=
            Teuchos::toString<float>(tmpproc_task_comm->
                proc_coords[j][origin_rank]) + ",";
          procs +=
            Teuchos::toString<float>(tmpproc_task_comm->
              proc_coords[j][origin_rank])+ " ";
        }
      }
      procs += "\n";

      gnuPlotArrow += " to ";


      for (int k = 0; k < a.size(); ++k) {
        int origin_task = a[k];

        for (int nind = task_communication_xadj[origin_task];
             nind < task_communication_xadj[origin_task + 1]; ++nind) {
          int neighbor_task = task_communication_adj[nind];

          bool differentnode = false;

          int neighbor_rank = this->getAssignedProcForTask(neighbor_task);

          for (int j = 0; j <  mindim; ++j) {
            if (int(tmpproc_task_comm->proc_coords[j][origin_rank]) !=
                int(tmpproc_task_comm->proc_coords[j][neighbor_rank])) {
              differentnode = true; break;
            }
          }
          std::tuple<int,int,int, int, int, int> foo(
              int(tmpproc_task_comm->proc_coords[0][origin_rank]),
              int(tmpproc_task_comm->proc_coords[1][origin_rank]),
              int(tmpproc_task_comm->proc_coords[2][origin_rank]),
              int(tmpproc_task_comm->proc_coords[0][neighbor_rank]),
              int(tmpproc_task_comm->proc_coords[1][neighbor_rank]),
              int(tmpproc_task_comm->proc_coords[2][neighbor_rank]));


          if (differentnode && my_arrows.find(foo) == my_arrows.end()) {
            my_arrows.insert(foo);

            std::string gnuPlotArrow2 = "";
            for (int j = 0; j <  mindim; ++j) {
              if (j == mindim - 1) {
                gnuPlotArrow2 +=
                  Teuchos::toString<float>(tmpproc_task_comm->
                      proc_coords[j][neighbor_rank]);
              }
              else {
                gnuPlotArrow2 +=
                  Teuchos::toString<float>(tmpproc_task_comm->
                      proc_coords[j][neighbor_rank]) + ",";
              }
            }
            ss += gnuPlotArrow + gnuPlotArrow2 + " nohead\n";
          }
        }
      }
    }

    std::ofstream procFile("procPlot.plot", std::ofstream::out);
    procFile << procs << "\n";
    procFile.close();

    //gnuPlotCode << ss;
    if (mindim == 2) {
      gnuPlotCode << "plot \"procPlot.plot\" with points pointsize 3\n";
    } else {
      gnuPlotCode << "splot \"procPlot.plot\" with points pointsize 3\n";
    }

    gnuPlotCode << ss << "\nreplot\n pause -1 \n";
    gnuPlotCode.close();
  }


// KDD Need to provide access to algorithm for getPartBoxes
#ifdef gnuPlot
  void writeGnuPlot(
      const Teuchos::Comm<int> *comm_,
      const Zoltan2::PartitioningSolution<Adapter> *soln_,
      int coordDim,
      tcoord_t **partCenters) {

    std::string file = "gggnuPlot";
    std::string exten = ".plot";
    std::ofstream mm("2d.txt");
    file += Teuchos::toString<int>(comm_->getRank()) + exten;
    std::ofstream ff(file.c_str());
//    ff.seekg(0, ff.end);
    std::vector<Zoltan2::coordinateModelPartBox <tcoord_t, part_t> >
      outPartBoxes =
        ((Zoltan2::PartitioningSolution<Adapter> *)soln_)->
          getPartBoxesView();

    for (part_t i = 0; i < this->ntasks;++i) {
      outPartBoxes[i].writeGnuPlot(ff, mm);
    }
    if (coordDim == 2) {
      ff << "plot \"2d.txt\"" << std::endl;
//      ff << "\n pause -1" << endl;
    }
    else {
      ff << "splot \"2d.txt\"" << std::endl;
//      ff << "\n pause -1" << endl;
    }
    mm.close();

    ff << "set style arrow 5 nohead size screen 0.03,15,135 ls 1"
       << std::endl;
    
    for (part_t i = 0; i < this->ntasks; ++i) {
      part_t pb = task_communication_xadj[i];
      part_t pe = task_communication_xadj[i + 1];

      for (part_t p = pb; p < pe; ++p) {
        part_t n = task_communication_adj[p];

//        std::cout << "i:" << i << " n:" << n << endl;
        std::string arrowline = "set arrow from ";
        for (int j = 0; j < coordDim - 1; ++j) {
          arrowline +=
            Teuchos::toString<tcoord_t>(partCenters[j][n]) + ",";
        }
        arrowline +=
          Teuchos::toString<tcoord_t>(partCenters[coordDim - 1][n]) +
            " to ";

        for (int j = 0; j < coordDim - 1; ++j) {
          arrowline +=
            Teuchos::toString<tcoord_t>(partCenters[j][i]) + ",";
        }
        arrowline +=
          Teuchos::toString<tcoord_t>(partCenters[coordDim - 1][i]) +
            " as 5\n";

//        std::cout << "arrow:" << arrowline << endl;
        ff << arrowline;
      }
    }

    ff << "replot\n pause -1" << std::endl;
    ff.close();
  }
#endif // gnuPlot

public:

  void getProcTask(part_t* &proc_to_task_xadj_,
                   part_t* &proc_to_task_adj_) {
    proc_to_task_xadj_ = this->proc_to_task_xadj.getRawPtr();
    proc_to_task_adj_ = this->proc_to_task_adj.getRawPtr();
  }

  virtual void map(const RCP<MappingSolution<Adapter> > &mappingsoln) {

    // Mapping was already computed in the constructor; we need to store it
    // in the solution.
    mappingsoln->setMap_RankForLocalElements(local_task_to_rank);

    // KDDKDD TODO:  Algorithm is also creating task_to_proc, which maybe
    // KDDKDD is not needed once we use MappingSolution to answer queries
    // KDDKDD instead of this algorithm.
    // KDDKDD Ask Mehmet:  what is the most efficient way to get the answer
    // KDDKDD out of CoordinateTaskMapper and into the MappingSolution?
  }


  virtual ~CoordinateTaskMapper() {
    //freeArray<part_t>(proc_to_task_xadj);
    //freeArray<part_t>(proc_to_task_adj);
    //freeArray<part_t>(task_to_proc);
    if(this->isOwnerofModel) {
      delete this->proc_task_comm;
    }
  }

  void create_local_task_to_rank(
      const lno_t num_local_coords,
      const part_t *local_coord_parts,
      const ArrayRCP<part_t> task_to_proc_) {
    local_task_to_rank = ArrayRCP <part_t>(num_local_coords);

    for (lno_t i = 0; i < num_local_coords; ++i) {
      part_t local_coord_part = local_coord_parts[i];
      part_t rank_index = task_to_proc_[local_coord_part];
      local_task_to_rank[i] = rank_index;
    }
  }



  /*! \brief Constructor.
   *  When this constructor is called, in order to calculate the
   *  communication metric, the task adjacency graph is created based on
   *  the coordinate model input and partitioning of it. If the
   *  communication graph is already calculated, use the other constructors.
   *
   *  \param comm_ is the communication object.
   *  \param machine_ is the machineRepresentation object. Stores the
   *         coordinates of machines.
   *  \param model_ is the input adapter.
   *  \param soln_ is the solution object. Holds the assignment of points.
   *  \param envConst_ is the environment object.
   */
  CoordinateTaskMapper(
          const Teuchos::RCP <const Teuchos::Comm<int> > comm_,
          const Teuchos::RCP <const MachineRepresentation<pcoord_t, part_t> >
            machine_,
          const Teuchos::RCP <const Adapter> input_adapter_,
          const Teuchos::RCP <const Zoltan2::PartitioningSolution<Adapter> >
            soln_,
          const Teuchos::RCP <const Environment> envConst,
          bool is_input_adapter_distributed = true,
          int num_ranks_per_node = 1,
          bool divide_to_prime_first = false,
          bool reduce_best_mapping = true):
      PartitionMapping<Adapter>(comm_, machine_, input_adapter_,
                                soln_, envConst),
      proc_to_task_xadj(0),
      proc_to_task_adj(0),
      task_to_proc(0),
      isOwnerofModel(true),
      proc_task_comm(0),
      task_communication_xadj(0),
      task_communication_adj(0),
      task_communication_edge_weight(0) {

    using namespace Teuchos;
    typedef typename Adapter::base_adapter_t ctm_base_adapter_t;

    RCP<Zoltan2::GraphModel<ctm_base_adapter_t> > graph_model_;
    RCP<Zoltan2::CoordinateModel<ctm_base_adapter_t> > coordinateModel_ ;

    RCP<const Teuchos::Comm<int> > rcp_comm = comm_;
    RCP<const Teuchos::Comm<int> > ia_comm = rcp_comm;
    if (!is_input_adapter_distributed) {
      ia_comm =  Teuchos::createSerialComm<int>();
    }

    RCP<const Environment> envConst_ = envConst;

    RCP<const ctm_base_adapter_t> baseInputAdapter_(
        rcp(dynamic_cast<const ctm_base_adapter_t *>(
            input_adapter_.getRawPtr()), false));

    modelFlag_t coordFlags_, graphFlags_;

    //create coordinate model
    //since this is coordinate task mapper,
    //the adapter has to have the coordinates
    coordinateModel_ = rcp(new CoordinateModel<ctm_base_adapter_t>(
          baseInputAdapter_, envConst_, ia_comm, coordFlags_));

    //if the adapter has also graph model, we will use graph model
    //to calculate the cost mapping.
    BaseAdapterType inputType_ = input_adapter_->adapterType();
    if (inputType_ == MatrixAdapterType ||
        inputType_ == GraphAdapterType ||
        inputType_ == MeshAdapterType)
    {
      graph_model_ = rcp(new GraphModel<ctm_base_adapter_t>(
          baseInputAdapter_, envConst_, ia_comm,
            graphFlags_));
    }

    if (!machine_->hasMachineCoordinates()) {
      throw std::runtime_error("Existing machine does not provide "
                               "coordinates for coordinate task mapping");
    }

    //if mapping type is 0 then it is coordinate mapping
    int procDim = machine_->getMachineDim();
    this->nprocs = machine_->getNumRanks();

    //get processor coordinates.
    pcoord_t **procCoordinates = NULL;
    if (!machine_->getAllMachineCoordinatesView(procCoordinates)) {
      throw std::runtime_error("Existing machine does not implement "
                               "getAllMachineCoordinatesView");
    }

    //get the machine extent.
    //if we have machine extent,
    //if the machine has wrap-around links, we would like to shift the
    //coordinates, so that the largest hap would be the wrap-around.
    std::vector<int> machine_extent_vec(procDim);
    //std::vector<bool> machine_extent_wrap_around_vec(procDim, 0);
    int *machine_extent = &(machine_extent_vec[0]);
    bool *machine_extent_wrap_around = new bool[procDim];
    for (int i = 0; i < procDim; ++i)
      machine_extent_wrap_around[i] = false;

    bool haveWrapArounds = machine_->getMachineExtentWrapArounds(machine_extent_wrap_around);

    // KDDKDD ASK MEHMET:  SHOULD WE GET AND USE machine_dimension HERE IF IT
    // KDDKDD ASK MEHMET:  IS PROVIDED BY THE MACHINE REPRESENTATION?
    // KDDKDD ASK MEHMET:  IF NOT HERE, THEN WHERE?
    // MD: Yes, I ADDED BELOW:
    if (machine_->getMachineExtent(machine_extent) &&
        haveWrapArounds) {

      procCoordinates =
          this->shiftMachineCoordinates(
              procDim,
              machine_extent,
              machine_extent_wrap_around,
              this->nprocs,
              procCoordinates);
    }

    //get the tasks information, such as coordinate dimension,
    //number of parts.
    int coordDim = coordinateModel_->getCoordinateDim();

//    int coordDim = machine_->getMachineDim();

    this->ntasks = soln_->getActualGlobalNumberOfParts();
    if (part_t(soln_->getTargetGlobalNumberOfParts()) > this->ntasks) {
      this->ntasks = soln_->getTargetGlobalNumberOfParts();
    }
    this->solution_parts = soln_->getPartListView();

    //we need to calculate the center of parts.
    tcoord_t **partCenters = new tcoord_t*[coordDim];
    for (int i = 0; i < coordDim; ++i) {
      partCenters[i] = new tcoord_t[this->ntasks];
    }

    typedef typename Adapter::scalar_t t_scalar_t;


    envConst->timerStart(MACRO_TIMERS, "Mapping - Solution Center");
 
    //get centers for the parts.
    getSolutionCenterCoordinates<Adapter, t_scalar_t,part_t>(
        envConst.getRawPtr(),
        ia_comm.getRawPtr(),
        coordinateModel_.getRawPtr(),
        this->solution_parts,
//        soln_->getPartListView();
//        this->soln.getRawPtr(),
        coordDim,
        ntasks,
        partCenters);

    envConst->timerStop(MACRO_TIMERS, "Mapping - Solution Center");

    //create the part graph
    if (graph_model_.getRawPtr() != NULL) {
      getCoarsenedPartGraph<Adapter, t_scalar_t, part_t>(
          envConst.getRawPtr(),
          ia_comm.getRawPtr(),
          graph_model_.getRawPtr(),
          this->ntasks,
          this->solution_parts,
//          soln_->getPartListView(),
//          this->soln.getRawPtr(),
          task_communication_xadj,
          task_communication_adj,
          task_communication_edge_weight
      );
    }
  
    //create coordinate communication model.
    this->proc_task_comm = new Zoltan2::CoordinateCommunicationModel<
      pcoord_t, tcoord_t, part_t, node_t>(
            procDim,
            procCoordinates,
            coordDim,
            partCenters,
            this->nprocs,
            this->ntasks,
            machine_extent,
            machine_extent_wrap_around,
            machine_.getRawPtr());

    int myRank = comm_->getRank();
    this->proc_task_comm->num_ranks_per_node = num_ranks_per_node ;
    this->proc_task_comm->divide_to_prime_first = divide_to_prime_first;

    envConst->timerStart(MACRO_TIMERS, "Mapping - Processor Task map");
    this->doMapping(myRank, comm_);
    envConst->timerStop(MACRO_TIMERS, "Mapping - Processor Task map");

    envConst->timerStart(MACRO_TIMERS, "Mapping - Communication Graph");

/*
    soln_->getCommunicationGraph(task_communication_xadj,
                                 task_communication_adj);
*/

    envConst->timerStop(MACRO_TIMERS, "Mapping - Communication Graph");
  #ifdef gnuPlot1
    if (comm_->getRank() == 0) {

      part_t taskCommCount = task_communication_xadj.size();
      std::cout << " TotalComm:"
        << task_communication_xadj[taskCommCount] << std::endl;
      part_t maxN = task_communication_xadj[0];
      for (part_t i = 1; i <= taskCommCount; ++i) {
        part_t nc =
          task_communication_xadj[i] - task_communication_xadj[i-1];
        if (maxN < nc)
          maxN = nc;
      }
      std::cout << " maxNeighbor:" << maxN << std::endl;
    }

    this->writeGnuPlot(comm_, soln_, coordDim, partCenters);
  #endif

    envConst->timerStart(MACRO_TIMERS, "Mapping - Communication Cost");
   
    if (reduce_best_mapping && task_communication_xadj.getRawPtr() &&
        task_communication_adj.getRawPtr()) {
      this->proc_task_comm->calculateCommunicationCost(
          task_to_proc.getRawPtr(),
          task_communication_xadj.getRawPtr(),
          task_communication_adj.getRawPtr(),
          task_communication_edge_weight.getRawPtr()
      );
    }

//    std::cout << "me: " << comm_->getRank() << " cost:"
//      << this->proc_task_comm->getCommunicationCostMetric() << std::endl;

    envConst->timerStop(MACRO_TIMERS, "Mapping - Communication Cost");

    //processors are divided into groups of size procDim! * coordDim!
    //each processor in the group obtains a mapping with a different
    //rotation and best one is broadcasted all processors.
    this->getBestMapping();
    this->create_local_task_to_rank(
        coordinateModel_->getLocalNumCoordinates(),
        this->solution_parts,
        this->task_to_proc);
/*
    {
      if (task_communication_xadj.getRawPtr() &&
          task_communication_adj.getRawPtr())
        this->proc_task_comm->calculateCommunicationCost(
            task_to_proc.getRawPtr(),
            task_communication_xadj.getRawPtr(),
            task_communication_adj.getRawPtr(),
            task_communication_edge_weight.getRawPtr()
        );
      std::cout << "me: " << comm_->getRank() << " cost:"
        << this->proc_task_comm->getCommunicationCostMetric() << std::endl;
    }
*/


  #ifdef gnuPlot
    this->writeMapping2(comm_->getRank());
  #endif

    delete [] machine_extent_wrap_around;
   
    if (machine_->getMachineExtent(machine_extent) &&
        haveWrapArounds) {
      for (int i = 0; i < procDim; ++i) {
        delete [] procCoordinates[i];
      }
      delete [] procCoordinates;
    }

    for (int i = 0; i < coordDim; ++i) {
      delete [] partCenters[i];
    }
    delete [] partCenters;


  }


  /*! \brief Constructor. Instead of Solution we have two parameters,
   *    numparts
   *
   *  When this constructor is called, in order to calculate the
   *  communication metric, the task adjacency graph is created based on
   *  the coordinate model input and partitioning of it. If the
   *  communication graph is already calculated, use the other
   *  constructors.
   *
   *  \param comm_ is the communication object.
   *  \param machine_ is the machineRepresentation object. Stores the
   *         coordinates of machines.
   *  \param model_ is the input adapter.
   *  \param soln_ is the solution object. Holds the assignment of points.
   *  \param envConst_ is the environment object.
   */
  CoordinateTaskMapper(
          const Teuchos::RCP <const Teuchos::Comm<int> > comm_,
          const Teuchos::RCP <const MachineRepresentation<pcoord_t,part_t> >
            machine_,
          const Teuchos::RCP <const Adapter> input_adapter_,
          const part_t num_parts_,
          const part_t *result_parts,
          const Teuchos::RCP <const Environment> envConst,
          bool is_input_adapter_distributed = true,
          int num_ranks_per_node = 1,
          bool divide_to_prime_first = false,
          bool reduce_best_mapping = true):
      PartitionMapping<Adapter>(comm_, machine_, input_adapter_,
                                num_parts_, result_parts, envConst),
      proc_to_task_xadj(0),
      proc_to_task_adj(0),
      task_to_proc(0),
      isOwnerofModel(true),
      proc_task_comm(0),
      task_communication_xadj(0),
      task_communication_adj(0),
      task_communication_edge_weight(0) {

    using namespace Teuchos;
    typedef typename Adapter::base_adapter_t ctm_base_adapter_t;

    RCP<Zoltan2::GraphModel<ctm_base_adapter_t> > graph_model_;
    RCP<Zoltan2::CoordinateModel<ctm_base_adapter_t> > coordinateModel_ ;

    RCP<const Teuchos::Comm<int> > rcp_comm = comm_;
    RCP<const Teuchos::Comm<int> > ia_comm = rcp_comm;
    if (!is_input_adapter_distributed) {
      ia_comm =  Teuchos::createSerialComm<int>();
    }
    RCP<const Environment> envConst_ = envConst;

    RCP<const ctm_base_adapter_t> baseInputAdapter_(
        rcp(dynamic_cast<const ctm_base_adapter_t *>(
            input_adapter_.getRawPtr()), false));

    modelFlag_t coordFlags_, graphFlags_;

    //create coordinate model
    //since this is coordinate task mapper,
    //the adapter has to have the coordinates
    coordinateModel_ = rcp(new CoordinateModel<ctm_base_adapter_t>(
          baseInputAdapter_, envConst_, ia_comm, coordFlags_));

    //if the adapter has also graph model, we will use graph model
    //to calculate the cost mapping.
    BaseAdapterType inputType_ = input_adapter_->adapterType();
    if (inputType_ == MatrixAdapterType ||
        inputType_ == GraphAdapterType ||
        inputType_ == MeshAdapterType)
    {
      graph_model_ = rcp(new GraphModel<ctm_base_adapter_t>(
          baseInputAdapter_, envConst_, ia_comm,
            graphFlags_));
    }

    if (!machine_->hasMachineCoordinates()) {
      throw std::runtime_error("Existing machine does not provide "
                               "coordinates for coordinate task mapping.");
    }

    //if mapping type is 0 then it is coordinate mapping
    int procDim = machine_->getMachineDim();
    this->nprocs = machine_->getNumRanks();

    //get processor coordinates.
    pcoord_t **procCoordinates = NULL;
    if (!machine_->getAllMachineCoordinatesView(procCoordinates)) {
      throw std::runtime_error("Existing machine does not implement "
                               "getAllMachineCoordinatesView");
    }

    //get the machine extent.
    //if we have machine extent,
    //if the machine has wrap-around links, we would like to shift the
    //coordinates,
    //so that the largest hap would be the wrap-around.
    std::vector<int> machine_extent_vec(procDim);
//    std::vector<bool> machine_extent_wrap_around_vec(procDim, 0);
    int *machine_extent = &(machine_extent_vec[0]);
    bool *machine_extent_wrap_around = new bool[procDim];
    bool haveWrapArounds = machine_->getMachineExtentWrapArounds(machine_extent_wrap_around);

    // KDDKDD ASK MEHMET:  SHOULD WE GET AND USE machine_dimension HERE IF IT
    // KDDKDD ASK MEHMET:  IS PROVIDED BY THE MACHINE REPRESENTATION?
    // KDDKDD ASK MEHMET:  IF NOT HERE, THEN WHERE?
    // MD: Yes, I ADDED BELOW:
    if (machine_->getMachineExtent(machine_extent) &&
        haveWrapArounds) {
      procCoordinates =
          this->shiftMachineCoordinates(
              procDim,
              machine_extent,
              machine_extent_wrap_around,
              this->nprocs,
              procCoordinates);
    }

    //get the tasks information, such as coordinate dimension,
    //number of parts.
    int coordDim = coordinateModel_->getCoordinateDim();

//    int coordDim = machine_->getMachineDim();


    this->ntasks = num_parts_;
    this->solution_parts = result_parts;

    //we need to calculate the center of parts.
    tcoord_t **partCenters = new tcoord_t*[coordDim];
    for (int i = 0; i < coordDim; ++i) {
      partCenters[i] = new tcoord_t[this->ntasks];
    }

    typedef typename Adapter::scalar_t t_scalar_t;


    envConst->timerStart(MACRO_TIMERS, "Mapping - Solution Center");

    //get centers for the parts.
    getSolutionCenterCoordinates<Adapter, t_scalar_t,part_t>(
        envConst.getRawPtr(),
        ia_comm.getRawPtr(),
        coordinateModel_.getRawPtr(),
        this->solution_parts,
//        soln_->getPartListView();
//        this->soln.getRawPtr(),
        coordDim,
        ntasks,
        partCenters);

    envConst->timerStop(MACRO_TIMERS, "Mapping - Solution Center");

    envConst->timerStart(MACRO_TIMERS, "GRAPHCREATE");
    //create the part graph
    if (graph_model_.getRawPtr() != NULL) {
      getCoarsenedPartGraph<Adapter, t_scalar_t, part_t>(
          envConst.getRawPtr(),
          ia_comm.getRawPtr(),
          graph_model_.getRawPtr(),
          this->ntasks,
          this->solution_parts,
//          soln_->getPartListView(),
//          this->soln.getRawPtr(),
          task_communication_xadj,
          task_communication_adj,
          task_communication_edge_weight
      );
    }
    envConst->timerStop(MACRO_TIMERS, "GRAPHCREATE");

    envConst->timerStart(MACRO_TIMERS,
                         "CoordinateCommunicationModel Create");
    //create coordinate communication model.
    this->proc_task_comm = new Zoltan2::CoordinateCommunicationModel<
      pcoord_t, tcoord_t, part_t, node_t>(
            procDim,
            procCoordinates,
            coordDim,
            partCenters,
            this->nprocs,
            this->ntasks,
            machine_extent,
            machine_extent_wrap_around,
            machine_.getRawPtr());

    envConst->timerStop(MACRO_TIMERS,
                        "CoordinateCommunicationModel Create");


    this->proc_task_comm->num_ranks_per_node = num_ranks_per_node;
    this->proc_task_comm->divide_to_prime_first = divide_to_prime_first;

    int myRank = comm_->getRank();


    envConst->timerStart(MACRO_TIMERS, "Mapping - Processor Task map");
    this->doMapping(myRank, comm_);
    envConst->timerStop(MACRO_TIMERS, "Mapping - Processor Task map");


    envConst->timerStart(MACRO_TIMERS, "Mapping - Communication Graph");

/*
    soln_->getCommunicationGraph(task_communication_xadj,
                                 task_communication_adj);
*/

    envConst->timerStop(MACRO_TIMERS, "Mapping - Communication Graph");
  #ifdef gnuPlot1
    if (comm_->getRank() == 0) {

      part_t taskCommCount = task_communication_xadj.size();
      std::cout << " TotalComm:"
        << task_communication_xadj[taskCommCount] << std::endl;
      part_t maxN = task_communication_xadj[0];
      for (part_t i = 1; i <= taskCommCount; ++i) {
        part_t nc =
          task_communication_xadj[i] - task_communication_xadj[i - 1];
        if (maxN < nc)
          maxN = nc;
      }
      std::cout << " maxNeighbor:" << maxN << std::endl;
    }

    this->writeGnuPlot(comm_, soln_, coordDim, partCenters);
  #endif

    envConst->timerStart(MACRO_TIMERS, "Mapping - Communication Cost");

    if (reduce_best_mapping && task_communication_xadj.getRawPtr() &&
        task_communication_adj.getRawPtr()) {
      this->proc_task_comm->calculateCommunicationCost(
          task_to_proc.getRawPtr(),
          task_communication_xadj.getRawPtr(),
          task_communication_adj.getRawPtr(),
          task_communication_edge_weight.getRawPtr()
      );
    }

//    std::cout << "me: " << comm_->getRank() << " cost:"
//    << this->proc_task_comm->getCommunicationCostMetric() << std::endl;

    envConst->timerStop(MACRO_TIMERS, "Mapping - Communication Cost");

    //processors are divided into groups of size procDim! * coordDim!
    //each processor in the group obtains a mapping with a different rotation
    //and best one is broadcasted all processors.
    this->getBestMapping();

    this->create_local_task_to_rank(
         coordinateModel_->getLocalNumCoordinates(),
         this->solution_parts,
         this->task_to_proc);
/*
    {
      if (task_communication_xadj.getRawPtr() &&
          task_communication_adj.getRawPtr())
        this->proc_task_comm->calculateCommunicationCost(
            task_to_proc.getRawPtr(),
            task_communication_xadj.getRawPtr(),
            task_communication_adj.getRawPtr(),
            task_communication_edge_weight.getRawPtr()
        );
      std::cout << "me: " << comm_->getRank() << " cost:"
        << this->proc_task_comm->getCommunicationCostMetric() << std::endl;
    }
*/



  #ifdef gnuPlot
    this->writeMapping2(comm_->getRank());
  #endif

    delete [] machine_extent_wrap_around;
   
    if (machine_->getMachineExtent(machine_extent) &&
        haveWrapArounds) {
      for (int i = 0; i < procDim; ++i) {
        delete [] procCoordinates[i];
      }
      delete [] procCoordinates;
    }

    for (int i = 0; i < coordDim; ++i) {
      delete [] partCenters[i];
    }
    delete [] partCenters;
  }

  /*! \brief Constructor
   * The mapping constructor which will also perform the mapping operation.
   * The result mapping can be obtained by
   *        --getAssignedProcForTask function: which returns the assigned
   *          processor id for the given task
   *        --getPartsForProc: which returns the assigned tasks with the
   *          number of tasks.
   *
   *      -task_comm_xadj, task_comm_adj, task_communication_edge_weight_
   *        can be provided NULL. In this case all processors will
   *        calculate the same mapping.
   *      -If task_comm_xadj, task_comm_adj and provided, algorithm will
   *        perform rotations and processors will calculate different
   *        mappings, and best one will be reduced.
   *      -If task_communication_edge_weight_ is provided with
   *        task_comm_xadj, task_comm_adj this will be used when cost is
   *        calculated.
   *      -recursion_depth is a mandatory argument. In the case
   *        part_no_array is not null, this parameter should represent the
   *        length of part_no_array.
   *        If part_no_array is given as NULL, then this will give the
   *        recursion depth for the algorith,
   *        Maximum number is ceil(log_2(min(num_processors, num_tasks))),
   *        and providing a higher number will be equivalant to this.
   *        Partitioning algorithm will work as RCB when maximum number is
   *        given, which performs the best mapping results.
   *      -part_no_array: The best results are obtained when this parameter
   *        is given as NULL. But if this is provided, partitioning will
   *        use this array for partitioning each dimension to the given
   *        numbers.
   *        The multiplication of these numbers should be equal to
   *        min(num_processors, num_tasks).
   *      -machine_dimensions: This can be NULL, but if provided the
   *        algorithm will perform shift of the machine coords so that
   *        the largest gap is treated as wrap-around link.
   *
   *  \param env_const_ the environment object.
   *  \param problemComm is the communication object.
   *  \param proc_dim dimensions of the processor coordinates.
   *  \param num_processors is the number of processors
   *  \param machine_coords is the coordinates of the processors.
   *
   *  \param task_dim is the dimension of the tasks.
   *  \param num_tasks is the number of tasks.
   *  \param task_coords is the coordinates of the tasks.
   *  \param task_comm_xadj is the task communication graphs xadj array.
   *        (task i adjacency is between task_comm_xadj[i] and
   *         task_comm_xadj[i + 1])
   *  \param task_comm_adj is task communication graphs adj array.
   *  \param task_communication_edge_weight_ is the weight of the
   *         communication in task graph.
   *  \param recursion_depth is the recursion depth that will be applied to
   *         partitioning.
   *         If part_no_array is provided, then it is the length of this
   *         array.
   *  \param part_no_array if part_no_array is provided, partitioning
   *         algorithm will be forced to use this array for partitioning.
   *         However, the multiplication of each entries in this array
   *         should be equal to min(num_processors, num_tasks).
   *  \param *machine_dimensions: the dimensions of the machine network.
   *         For example for hopper 17x8x24
   *         This can be NULL, but if provided the algorithm will perform
   *         shift of the machine coords so that the largest gap is treated
   *         as wrap-around link.
   */
  CoordinateTaskMapper(
          const Environment *env_const_,
          const Teuchos::Comm<int> *problemComm,
          int proc_dim,
          int num_processors,
          pcoord_t **machine_coords,
          int task_dim,
          part_t num_tasks,
          tcoord_t **task_coords,
          ArrayRCP<part_t>task_comm_xadj,
          ArrayRCP<part_t>task_comm_adj,
          pcoord_t *task_communication_edge_weight_,
          int recursion_depth,
          Kokkos::View<part_t *, Kokkos::HostSpace> part_no_array,
          const part_t *machine_dimensions,
          int num_ranks_per_node = 1,
          bool divide_to_prime_first = false,
          bool reduce_best_mapping = true):
      PartitionMapping<Adapter>(
        Teuchos::rcpFromRef<const Teuchos::Comm<int> >(*problemComm),
      Teuchos::rcpFromRef<const Environment>(*env_const_)),
      proc_to_task_xadj(0),
      proc_to_task_adj(0),
      task_to_proc(0),
      isOwnerofModel(true),
      proc_task_comm(0),
      task_communication_xadj(task_comm_xadj),
      task_communication_adj(task_comm_adj) {

    //if mapping type is 0 then it is coordinate mapping
    pcoord_t ** virtual_machine_coordinates  = machine_coords;
    bool *wrap_arounds = new bool [proc_dim];
    for (int i = 0; i < proc_dim; ++i) wrap_arounds[i] = true;

    if (machine_dimensions) {
      virtual_machine_coordinates =
          this->shiftMachineCoordinates(
              proc_dim,
              machine_dimensions,
              wrap_arounds,
              num_processors,
              machine_coords);
    }

    this->nprocs = num_processors;

    int coordDim = task_dim;
    this->ntasks = num_tasks;

    //alloc memory for part centers.
    tcoord_t **partCenters = task_coords;

    //create coordinate communication model.
    this->proc_task_comm =
        new Zoltan2::CoordinateCommunicationModel<pcoord_t,tcoord_t,part_t,node_t>(
            proc_dim,
            virtual_machine_coordinates,
            coordDim,
            partCenters,
            this->nprocs,
            this->ntasks, NULL, NULL
        );

    this->proc_task_comm->num_ranks_per_node = num_ranks_per_node;
    this->proc_task_comm->divide_to_prime_first = divide_to_prime_first;

    this->proc_task_comm->setPartArraySize(recursion_depth);
    this->proc_task_comm->setPartArray(part_no_array);

    int myRank = problemComm->getRank();

    this->doMapping(myRank, this->comm);
#ifdef gnuPlot
    this->writeMapping2(myRank);
#endif

  // MDM added this edge case - for example if NX = 1 NY = 1 NZ = 1
  // That would pass on original develop so updated this so now it will also pass.
    if (reduce_best_mapping && task_communication_xadj.size() &&
        task_communication_adj.size()) {
      this->proc_task_comm->calculateCommunicationCost(
          task_to_proc.getRawPtr(),
          task_communication_xadj.getRawPtr(),
          task_communication_adj.getRawPtr(),
          task_communication_edge_weight_
      );


      this->getBestMapping();

/*
      if (myRank == 0) {
        this->proc_task_comm->calculateCommunicationCost(
            task_to_proc.getRawPtr(),
            task_communication_xadj.getRawPtr(),
            task_communication_adj.getRawPtr(),
            task_communication_edge_weight_
        );
        cout << "me: " << problemComm->getRank() << " cost:"
             << this->proc_task_comm->getCommunicationCostMetric() << endl;
      }
*/

    }

    delete [] wrap_arounds;

    if (machine_dimensions) {
      for (int i = 0; i < proc_dim; ++i) {
        delete [] virtual_machine_coordinates[i];
      }
      delete [] virtual_machine_coordinates;
    }
#ifdef gnuPlot
    if (problemComm->getRank() == 0)
      this->writeMapping2(-1);
#endif
  }


  /*
  double getCommunicationCostMetric() {
    return this->proc_task_comm->getCommCost();
  }
  */

  /*! \brief Returns the number of parts to be assigned to this process.
   */
  virtual size_t getLocalNumberOfParts() const{
    return 0;
  }

  /*! \brief Using the machine dimensions provided, create virtual machine
   * coordinates by assigning the largest gap to be as the wrap around
   * link.
   * \param machine_dim: the number of dimensions in the machine network.
   * \param machine_dimensions: the dimension of the machien network. For
   * example for hopper, 17,8,24
   *
   * \param numProcs: the number of allocated processors.
   * \param mCoords: allocated machine coordinates.
   */
  pcoord_t** shiftMachineCoordinates(
      int machine_dim,
      const part_t *machine_dimensions,
      bool *machine_extent_wrap_around,
      part_t numProcs,
      pcoord_t **mCoords) {
    
    pcoord_t **result_machine_coords = NULL;
    result_machine_coords = new pcoord_t*[machine_dim];
    
    for (int i = 0; i < machine_dim; ++i) {
      result_machine_coords[i] = new pcoord_t [numProcs];
    }

    for (int i = 0; i < machine_dim; ++i) {
      part_t numMachinesAlongDim = machine_dimensions[i];

      part_t *machineCounts = new part_t[numMachinesAlongDim];
      memset(machineCounts, 0, sizeof(part_t) * numMachinesAlongDim);

      int *filledCoordinates = new int[numMachinesAlongDim];

      pcoord_t *coords = mCoords[i];

      for (part_t j = 0; j < numProcs; ++j) {
        part_t mc = (part_t) coords[j];
        ++machineCounts[mc];
      }

      part_t filledCoordinateCount = 0;
      for (part_t j = 0; j < numMachinesAlongDim; ++j) {
        if (machineCounts[j] > 0) {
          filledCoordinates[filledCoordinateCount++] = j;
        }
      }

      part_t firstProcCoord = filledCoordinates[0];
      part_t firstProcCount = machineCounts[firstProcCoord];

      part_t lastProcCoord = filledCoordinates[filledCoordinateCount - 1];
      part_t lastProcCount = machineCounts[lastProcCoord];

      part_t firstLastGap =
        numMachinesAlongDim - lastProcCoord + firstProcCoord;
      part_t firstLastGapProc = lastProcCount + firstProcCount;

      part_t leftSideProcCoord = firstProcCoord;
      part_t leftSideProcCount = firstProcCount;
      part_t biggestGap = 0;
      part_t biggestGapProc = numProcs;

      part_t shiftBorderCoordinate = -1;
      for (part_t j = 1; j < filledCoordinateCount; ++j) {
        part_t rightSideProcCoord= filledCoordinates[j];
        part_t rightSideProcCount = machineCounts[rightSideProcCoord];

        part_t gap = rightSideProcCoord - leftSideProcCoord;
        part_t gapProc = rightSideProcCount + leftSideProcCount;

        // Pick the largest gap in this dimension. Use fewer process on
        // either side of the largest gap to break the tie. An easy
        // addition to this would be to weight the gap by the number of
        // processes.
        if (gap > biggestGap ||
           (gap == biggestGap && biggestGapProc > gapProc)) {
          shiftBorderCoordinate = rightSideProcCoord;
          biggestGapProc = gapProc;
          biggestGap = gap;
        }
        leftSideProcCoord = rightSideProcCoord;
        leftSideProcCount = rightSideProcCount;
      }


      if (!(biggestGap > firstLastGap ||
         (biggestGap == firstLastGap &&
          biggestGapProc < firstLastGapProc))) {
        shiftBorderCoordinate = -1;
      }

      for (part_t j = 0; j < numProcs; ++j) {

        if (machine_extent_wrap_around[i] &&
            coords[j] < shiftBorderCoordinate) {
          result_machine_coords[i][j] = coords[j] + numMachinesAlongDim;

        }
        else {
          result_machine_coords[i][j] = coords[j];
        }
      }
      delete [] machineCounts;
      delete [] filledCoordinates;
    }

    return result_machine_coords;

  }

  /*! \brief getAssignedProcForTask function,
   * returns the assigned tasks with the number of tasks.
   *  \param procId procId being queried.
   *  \param numProcs (output), the number of processor the part is
   *    assigned to.
   *  \param procs (output), the list of processors assigned to given part..
   */
  virtual void getProcsForPart(part_t taskId, part_t &numProcs,
                               part_t *&procs) const {
    numProcs = 1;
    procs = this->task_to_proc.getRawPtr() + taskId;
  }

  /*! \brief getAssignedProcForTask function, returns the assigned
   *    processor id for the given task
   *  \param taskId taskId being queried.
   */
  inline part_t getAssignedProcForTask(part_t taskId) {
    return this->task_to_proc[taskId];
  }

  /*! \brief getAssignedProcForTask function,
   * returns the assigned tasks with the number of tasks.
   *  \param procId procId being queried.
   *  \param numParts (output), the number of parts the processor is
   *    assigned to.
   *  \param parts (output), the list of parts assigned to given processor..
   */
  virtual void getPartsForProc(int procId, part_t &numParts,
                               part_t *&parts) const {

    part_t task_begin = this->proc_to_task_xadj[procId];
    part_t taskend = this->proc_to_task_xadj[procId + 1];

    parts = this->proc_to_task_adj.getRawPtr() + task_begin;
    numParts = taskend - task_begin;
  }

  ArrayView<part_t> getAssignedTasksForProc(part_t procId) {
    part_t task_begin = this->proc_to_task_xadj[procId];
    part_t taskend = this->proc_to_task_xadj[procId + 1];

/*
  std::cout << "part_t:" << procId << " taskCount:"
    << taskend - task_begin << std::endl;

  for (part_t i = task_begin; i < taskend; ++i) {
    std::cout << "part_t:" << procId << " task:"
      << proc_to_task_adj[i] << endl;
  }
*/
    if (taskend - task_begin > 0) {
      ArrayView <part_t> assignedParts(
          this->proc_to_task_adj.getRawPtr() + task_begin,
          taskend - task_begin);
      
      return assignedParts;
    }
    else {
      ArrayView <part_t> assignedParts;
      
      return assignedParts;
    }
  }

};

/*! \brief Constructor
 * The interface function that calls CoordinateTaskMapper which will also
 * perform the mapping operation.
 * The result mapping can be obtained by
 *    -proc_to_task_xadj: Holds the beginning and end indices of
 *     tasks on proc_to_task_adj that is assigned to a processor.
 *     The tasks assigned to processor i are between proc_to_task_xadj[i]
 *     and proc_to_task_xadj[i + 1] on proc_to_task_adj.
 *
 *    -proc_to_task_adj: holds the task adj array.
 *
 *    -task_comm_xadj, task_comm_adj, task_communication_edge_weight_
 *     can be provided NULL.
 *     In this case all processors will calculate the same mapping.
 *
 *    -If task_comm_xadj, task_comm_adj and provided, algorithm will
 *     perform rotations, and processors will calculate different
 *     mappings, and best one will be reduced.
 *
 *    -If task_communication_edge_weight_ is provided with
 *     task_comm_xadj, task_comm_adj, this will be used when cost is
 *     calculated.
 *
 *    -recursion_depth is a mandatory argument. In the case part_no_array
 *     is not null, this parameter
 *     should represent the length of part_no_array.
 *     If part_no_array is given as NULL, then this will give the
 *     recursion depth for the algorithm,
 *     Maximum number is ceil(log_2(min(num_processors, num_tasks))),
 *     and providing a higher number will
 *     be equivalant to this. Partitioning algorithm will work as RCB
 *     when maximum number is given, which performs the best mapping
 *     results.
 *
 *    -part_no_array: The best results are obtained when this parameter
 *     is given as NULL. But if this is provided, partitioning will use
 *     this array for partitioning each dimension to the given numbers.
 *     The multiplication of these numbers should be equal to
 *     min(num_processors, num_tasks).
 *
 *    -machine_dimensions: This can be NULL, but if provided the algorithm
 *     will perform shift of the machine coords so that
 *     the largest gap is treated as wrap-around link.
 *
 *  \param problemComm is the communication object.
 *  \param proc_dim dimensions of the processor coordinates.
 *  \param num_processors is the number of processors
 *  \param machine_coords is the coordinates of the processors.
 *
 *  \param task_dim is the dimension of the tasks.
 *  \param num_tasks is the number of tasks.
 *  \param task_coords is the coordinates of the tasks.
 *  \param task_comm_xadj is the task communication graphs xadj array.
 *        (task i's adjacency is between task_comm_xadj[i] and
 *         task_comm_xadj[i + 1])
 *  \param task_comm_adj is task communication graphs adj array.
 *  \param task_communication_edge_weight_ is the weight of the
 *         communication in task graph.
 *  \param proc_to_task_xadj is is the output for tasks showing which proc
 *         has the which parts.
 *         (proc-i will own the tasks from proc_to_task_xadj[i] to
 *         proc_to_task_xadj[i + 1])
 *  \param proc_to_task_adj is the ouput list of tasks pointed by
 *         proc_to_task_xadj
 *  \param recursion_depth is the recursion depth that will be applied to
 *         partitioning.
 *         If part_no_array is provided, then it is the length of this
 *         array.
 *  \param part_no_array if part_no_array is provided, partitioning
 *         algorithm will be forced to use this array for partitioning.
 *         However, the multiplication of each entries in this array
 *         should be equal to min(num_processors, num_tasks).
 *  \param *machine_dimensions: the dimensions of the machine network. For
 *         example for hopper 17x8x24
 *         This can be NULL, but if provided the algorithm will perform
 *         shift of the machine coords so that
 *         the largest gap is treated as wrap-around link.
 */
template <typename part_t, typename pcoord_t, typename tcoord_t>
void coordinateTaskMapperInterface(
    RCP<const Teuchos::Comm<int> > problemComm,
    int proc_dim,
    int num_processors,
    pcoord_t **machine_coords,
    int task_dim,
    part_t num_tasks,
    tcoord_t **task_coords,
    part_t *task_comm_xadj,
    part_t *task_comm_adj,
    // float-like, same size with task_communication_adj_ weight of the
    // corresponding edge.
    pcoord_t *task_communication_edge_weight_,
    part_t *proc_to_task_xadj, /*output*/
    part_t *proc_to_task_adj, /*output*/
    int recursion_depth,
    Kokkos::View<part_t *, Kokkos::HostSpace> part_no_array,
    const part_t *machine_dimensions,
    int num_ranks_per_node = 1,
    bool divide_to_prime_first = false) {

  const Environment *envConst_ = new Environment(problemComm);

  // mfh 03 Mar 2015: It's OK to omit the Node template
  // parameter in Tpetra, if you're just going to use the
  // default Node.
  typedef Tpetra::MultiVector<tcoord_t, part_t, part_t> tMVector_t;

  Teuchos::ArrayRCP<part_t> task_communication_xadj(
      task_comm_xadj, 0, num_tasks + 1, false);

  Teuchos::ArrayRCP<part_t> task_communication_adj;
  if (task_comm_xadj) {
    Teuchos::ArrayRCP<part_t> tmp_task_communication_adj(
        task_comm_adj, 0, task_comm_xadj[num_tasks], false);
    task_communication_adj = tmp_task_communication_adj;
  }


  CoordinateTaskMapper<XpetraMultiVectorAdapter<tMVector_t>, part_t> *ctm =
      new CoordinateTaskMapper<XpetraMultiVectorAdapter<tMVector_t>, part_t>(
      envConst_,
      problemComm.getRawPtr(),
      proc_dim,
      num_processors,
      machine_coords,
//      machine_coords_,

      task_dim,
      num_tasks,
      task_coords,

      task_communication_xadj,
      task_communication_adj,
      task_communication_edge_weight_,
      recursion_depth,
      part_no_array,
      machine_dimensions,
      num_ranks_per_node,
      divide_to_prime_first);


  part_t* proc_to_task_xadj_;
  part_t* proc_to_task_adj_;

  ctm->getProcTask(proc_to_task_xadj_, proc_to_task_adj_);

  for (part_t i = 0; i <= num_processors; ++i) {
    proc_to_task_xadj[i] = proc_to_task_xadj_[i];
  }

  for (part_t i = 0; i < num_tasks; ++i) {
    proc_to_task_adj[i] = proc_to_task_adj_[i];
  }
  
  delete ctm;
  delete envConst_;
}

template <typename proc_coord_t, typename v_lno_t>
inline void visualize_mapping(int myRank,
    const int machine_coord_dim,
    const int num_ranks,
    proc_coord_t **machine_coords,
    const v_lno_t num_tasks,
    const v_lno_t *task_communication_xadj,
    const v_lno_t *task_communication_adj,
    const int *task_to_rank) {

  std::string rankStr = Teuchos::toString<int>(myRank);
  std::string gnuPlots = "gnuPlot", extentionS = ".plot";
  std::string outF = gnuPlots + rankStr+ extentionS;
  std::ofstream gnuPlotCode( outF.c_str(), std::ofstream::out);

  if (machine_coord_dim != 3) {
    std::cerr << "Mapping Write is only good for 3 dim" << std::endl;
    return;
  }
  std::string ss = "";
  std::string procs = "";

  std::set<std::tuple<int, int, int, int, int, int> > my_arrows;
  
  for (v_lno_t origin_task = 0; origin_task < num_tasks; ++origin_task) {
    int origin_rank = task_to_rank[origin_task];
    std::string gnuPlotArrow = "set arrow from ";

    for (int j = 0; j <  machine_coord_dim; ++j) {
      if (j == machine_coord_dim - 1) {
        gnuPlotArrow +=
          Teuchos::toString<proc_coord_t>(machine_coords[j][origin_rank]);
        procs +=
          Teuchos::toString<proc_coord_t>(machine_coords[j][origin_rank]);

      }
      else {
        gnuPlotArrow +=
          Teuchos::toString<proc_coord_t>(machine_coords[j][origin_rank])
            + ",";
        procs +=
          Teuchos::toString<proc_coord_t>(machine_coords[j][origin_rank])
            + " ";
      }
    }
    procs += "\n";

    gnuPlotArrow += " to ";


    for (int nind = task_communication_xadj[origin_task];
         nind < task_communication_xadj[origin_task + 1]; ++nind) {

      int neighbor_task = task_communication_adj[nind];

      bool differentnode = false;
      int neighbor_rank = task_to_rank[neighbor_task];

      for (int j = 0; j <  machine_coord_dim; ++j) {
        if (int(machine_coords[j][origin_rank]) !=
            int(machine_coords[j][neighbor_rank])) {
          differentnode = true; break;
        }
      }
      
      std::tuple<int,int,int, int, int, int> foo(
          (int)(machine_coords[0][origin_rank]),
          (int)(machine_coords[1][origin_rank]),
          (int)(machine_coords[2][origin_rank]),
          (int)(machine_coords[0][neighbor_rank]),
          (int)(machine_coords[1][neighbor_rank]),
          (int)(machine_coords[2][neighbor_rank]));

      if (differentnode && my_arrows.find(foo) == my_arrows.end()) {
        my_arrows.insert(foo);

        std::string gnuPlotArrow2 = "";
        for (int j = 0; j <  machine_coord_dim; ++j) {
          if (j == machine_coord_dim - 1) {
            gnuPlotArrow2 +=
              Teuchos::toString<float>(machine_coords[j][neighbor_rank]);
          }
          else {
            gnuPlotArrow2 +=
              Teuchos::toString<float>(machine_coords[j][neighbor_rank])
                + ",";
          }
        }
        ss += gnuPlotArrow + gnuPlotArrow2 + " nohead\n";
      }
    }
  }

  std::ofstream procFile("procPlot.plot", std::ofstream::out);
  procFile << procs << "\n";
  procFile.close();

  //gnuPlotCode << ss;
  if (machine_coord_dim == 2) {
    gnuPlotCode << "plot \"procPlot.plot\" with points pointsize 3\n";
  }
  else {
    gnuPlotCode << "splot \"procPlot.plot\" with points pointsize 3\n";
  }

  gnuPlotCode << ss << "\nreplot\n pause -1\npause -1";
  gnuPlotCode.close();
}

} // namespace Zoltan2

#endif
