#ifndef _ZOLTAN2_COORD_PARTITIONMAPPING_HPP_
#define _ZOLTAN2_COORD_PARTITIONMAPPING_HPP_

#include <fstream>
#include <ctime>
#include <vector>
#include "Zoltan2_AlgPQJagged.hpp"
#include "Teuchos_ArrayViewDecl.hpp"
#include "Zoltan2_PartitionMapping.hpp"
#include "Zoltan2_MachineRepresentation.hpp"
#include "Teuchos_ReductionOp.hpp"
#include "Zoltan2_XpetraMultiVectorAdapter.hpp"

#include "Teuchos_ConfigDefs.hpp" // define HAVE_MPI
#include "Teuchos_Comm.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI

//#define gnuPlot

namespace Teuchos{

/*! \brief Zoltan2_ReduceBestMapping Class, reduces the minimum cost mapping, ties breaks with minimum proc id.
 */
template <typename Ordinal, typename T>
class Zoltan2_ReduceBestMapping  : public ValueTypeReductionOp<Ordinal,T>
{
private:
    T _EPSILON;

public:
    /*! \brief Default Constructor
     */
    Zoltan2_ReduceBestMapping ():_EPSILON (std::numeric_limits<T>::epsilon()){}

    /*! \brief Implement Teuchos::ValueTypeReductionOp interface
     */
    void reduce( const Ordinal count, const T inBuffer[], T inoutBuffer[]) const
    {

        for (Ordinal i=0; i < count; i++){
            if (inBuffer[0] - inoutBuffer[0] < -_EPSILON){
                inoutBuffer[0] = inBuffer[0];
                inoutBuffer[1] = inBuffer[1];
            } else if(
                    inBuffer[0] - inoutBuffer[0] < _EPSILON &&
                    inBuffer[1] - inoutBuffer[1] < _EPSILON){
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

template <typename gno_t, typename partId_t>
class GNO_LNO_PAIR{
public:
    gno_t gno;
    partId_t part;
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

template <typename partId_t>
void getGridCommunicationGraph(partId_t taskCount, partId_t *&task_comm_xadj, partId_t *&task_comm_adj, vector <int> grid_dims){
    int dim = grid_dims.size();
    int neighborCount = 2 * dim;
    task_comm_xadj = allocMemory< partId_t>(taskCount);
    task_comm_adj = allocMemory <partId_t>(taskCount * neighborCount);

    partId_t neighBorIndex = 0;
    for (partId_t i = 0; i < taskCount; ++i){
        partId_t prevDimMul = 1;
        for (int j = 0; j < dim; ++j){
            partId_t lNeighbor = i - prevDimMul;
            partId_t rNeighbor = i + prevDimMul;
            prevDimMul *= grid_dims[j];
            if (lNeighbor >= 0 &&  lNeighbor/ prevDimMul == i / prevDimMul && lNeighbor < taskCount){
                task_comm_adj[neighBorIndex++] = lNeighbor;
            }
            if (rNeighbor >= 0 && rNeighbor/ prevDimMul == i / prevDimMul && rNeighbor < taskCount){
                task_comm_adj[neighBorIndex++] = rNeighbor;
            }
        }
        task_comm_xadj[i] = neighBorIndex;
    }

}
//returns the center of the parts.
template <typename Adapter, typename scalar_t, typename partId_t>
void getSolutionCenterCoordinates(
        const Environment *envConst,
        const Teuchos::Comm<int> *comm,
        const Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *coords,
        const Zoltan2::PartitioningSolution<Adapter> *soln_,
        int coordDim,
        partId_t ntasks,
        scalar_t **partCenters){

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
    gno_t *point_counts = allocMemory<gno_t>(ntasks);
    memset(point_counts, 0, sizeof(gno_t) * ntasks);

    //global number of points in each part.
    gno_t *global_point_counts = allocMemory<gno_t>(ntasks);


    scalar_t **pqJagged_coordinates = allocMemory<scalar_t *>(coordDim);

    for (int dim=0; dim < coordDim; dim++){
        ArrayRCP<const scalar_t> ar;
        xyz[dim].getInputArray(ar);
        //pqJagged coordinate values assignment
        pqJagged_coordinates[dim] =  (scalar_t *)ar.getRawPtr();
        memset(partCenters[dim], 0, sizeof(scalar_t) * ntasks);
    }

    //get parts with parallel gnos.
    const partId_t *parts = soln_->getPartList();
    const gno_t *soln_gnos = soln_->getIdList();
/*
    for (lno_t i=0; i < numLocalCoords; i++){
        cout << "me:" << comm->getRank() << " gno:" << soln_gnos[i] << " tmp.part :" << parts[i]<< endl;
    }
    */


    envConst->timerStart(MACRO_TIMERS, "Mapping - Hashing Creation");
    //hash vector
    vector< vector <GNO_LNO_PAIR<gno_t, partId_t> > > hash(numLocalCoords);

    //insert each point in solution to hash.
    for (lno_t i=0; i < numLocalCoords; i++){
        GNO_LNO_PAIR<gno_t, partId_t> tmp;
        tmp.gno = soln_gnos[i];
        tmp.part = parts[i];
        //cout << "gno:" << tmp.gno << " tmp.part :" << tmp.part << endl;
        //count the local number of points in each part.
        ++point_counts[tmp.part];
        lno_t hash_index = tmp.gno % numLocalCoords;
        hash[hash_index].push_back(tmp);
    }

    envConst->timerStop(MACRO_TIMERS, "Mapping - Hashing Creation");
    //get global number of points in each part.
    reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_SUM,
            ntasks, point_counts, global_point_counts
    );



    envConst->timerStart(MACRO_TIMERS, "Mapping - Hashing Search");
    //add up all coordinates in each part.
    for (lno_t i=0; i < numLocalCoords; i++){
        gno_t g = gnos[i];
        lno_t hash_index = g % numLocalCoords;
        lno_t hash_size = hash[hash_index].size();
        partId_t p = -1;
        for (int j =0; j < hash_size; ++j){
            if (hash[hash_index][j].gno == g){
                p = hash[hash_index][j].part;
                break;
            }
        }
        if(p == -1) {
            cerr << "ERROR AT HASHING FOR GNO:"<< g << " LNO:" << i << endl;
        }
        //add uo all coordinates in each part.
        for(int j = 0; j < coordDim; ++j){
            scalar_t c = pqJagged_coordinates[j][i];
            partCenters[j][p] += c;
        }
    }
    envConst->timerStop(MACRO_TIMERS, "Mapping - Hashing Search");

    for(int j = 0; j < coordDim; ++j){
        for (partId_t i=0; i < ntasks; ++i){
            partCenters[j][i] /= global_point_counts[i];
        }
    }

    scalar_t *tmpCoords = allocMemory<scalar_t>(ntasks);
    for(int j = 0; j < coordDim; ++j){
        reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM,
                ntasks, partCenters[j], tmpCoords
        );

        scalar_t *tmp = partCenters[j];
        partCenters[j] = tmpCoords;
        tmpCoords = tmp;
    }

    freeArray<gno_t> (point_counts);
    freeArray<gno_t> (global_point_counts);

    freeArray<scalar_t> (tmpCoords);
    freeArray<scalar_t *>(pqJagged_coordinates);
}


/*! \brief KmeansHeap Class, max heap, but holds the minimum values.
 */
template <class IT, class WT>
class KmeansHeap{
    IT heapSize;
    IT *indices;
    WT *values;
    WT _EPSILON;


public:
    void setHeapsize(IT heapsize_){
        this->heapSize = heapsize_;
        this->indices = allocMemory<IT>(heapsize_ );
        this->values = allocMemory<WT>(heapsize_ );
        this->_EPSILON = numeric_limits<WT>::epsilon();
    }

    ~KmeansHeap(){
        freeArray<IT>(this->indices);
        freeArray<WT>(this->values);
    }


    void addPoint(IT index, WT distance){
        WT maxVal = this->values[0];
        //add only the distance is smaller than the maximum distance.
        //cout << "indeX:" << index << "distance:" <<distance << " maxVal:" << maxVal << endl;
        if (distance >= maxVal) return;
        else {
            this->values[0] = distance;
            this->indices[0] = index;
            this->push_down(0);
        }
    }

    //heap push down operation
    void push_down(IT index_on_heap){
        IT child_index1 = 2 * index_on_heap + 1;
        IT child_index2 = 2 * index_on_heap + 2;

        IT biggerIndex = -1;
        if(child_index1 < this->heapSize && child_index2 < this->heapSize){

            if (this->values[child_index1] < this->values[child_index2]){
                biggerIndex = child_index2;
            }
            else {
                biggerIndex = child_index1;
            }
        }
        else if(child_index1 < this->heapSize){
            biggerIndex = child_index1;

        }
        else if(child_index2 < this->heapSize){
            biggerIndex = child_index2;
        }
        if (biggerIndex >= 0 && this->values[biggerIndex] > this->values[index_on_heap]){
            WT tmpVal = this->values[biggerIndex];
            this->values[biggerIndex] = this->values[index_on_heap];
            this->values[index_on_heap] = tmpVal;

            IT tmpIndex = this->indices[biggerIndex];
            this->indices[biggerIndex] = this->indices[index_on_heap];
            this->indices[index_on_heap] = tmpIndex;
            this->push_down(biggerIndex);
        }
    }

    void initValues(){
        WT MAXVAL = numeric_limits<WT>::max();
        for(IT i = 0; i < this->heapSize; ++i){
            this->values[i] = MAXVAL;
            this->indices[i] = -1;
        }
    }

    //returns the total distance to center in the cluster.
    WT getTotalDistance(){

        WT nc = 0;
        for(IT j = 0; j < this->heapSize; ++j){
            nc += this->values[j];

            //cout << "index:" << this->indices[j] << " distance:" << this->values[j] << endl;
        }
        return nc;
    }

    //returns the new center of the cluster.
    bool getNewCenters(WT *center, WT **coords, int dimension){
        bool moved = false;
        for(int i = 0; i < dimension; ++i){
            WT nc = 0;
            for(IT j = 0; j < this->heapSize; ++j){
                IT k = this->indices[j];
                //cout << "i:" << i << " dim:" << dimension << " k:" << k << " heapSize:" << heapSize << endl;
                nc += coords[i][k];
            }
            nc /= this->heapSize;
            moved = (ABS(center[i] - nc) > this->_EPSILON || moved );
            center[i] = nc;

        }
        return moved;
    }

    void copyCoordinates(IT *permutation){
        for(IT i = 0; i < this->heapSize; ++i){
            permutation[i] = this->indices[i];
        }
    }
};

/*! \brief KMeansCluster Class
 */
template <class IT, class WT>
class KMeansCluster{

    int dimension;
    KmeansHeap<IT,WT> closestPoints;

public:
    WT *center;
    ~KMeansCluster(){
        freeArray<WT>(center);
    }

    void setParams(int dimension_, int heapsize){
        this->dimension = dimension_;
        this->center = allocMemory<WT>(dimension_);
        this->closestPoints.setHeapsize(heapsize);
    }

    void clearHeap(){
        this->closestPoints.initValues();
    }

    bool getNewCenters( WT **coords){
        return this->closestPoints.getNewCenters(center, coords, dimension);
    }

    //returns the distance of the coordinate to the center.
    //also adds it to the heap.
    WT getDistance(IT index, WT **elementCoords){
        WT distance = 0;
        for (int i = 0; i < this->dimension; ++i){
            WT d = (center[i] - elementCoords[i][index]);
            distance += d * d;
        }
        distance = pow(distance, WT(1.0 / this->dimension));
        closestPoints.addPoint(index, distance);
        return distance;
    }

    WT getDistanceToCenter(){
        return closestPoints.getTotalDistance();
    }

    void copyCoordinates(IT *permutation){
        closestPoints.copyCoordinates(permutation);
    }
};

/*! \brief KMeansAlgorithm Class that performs clustering of the coordinates, and returns the closest set of coordinates.
 * Useful to filter the processors, when there are more processors than needed.
 */
template <class IT, class WT>
class KMeansAlgorithm{

    int dim;
    IT numElements;
    WT **elementCoords;
    IT numClusters;
    IT required_elements;
    KMeansCluster <IT,WT> *clusters;
    WT *maxCoordinates;
    WT *minCoordinates;
public:
    ~KMeansAlgorithm(){
        freeArray<KMeansCluster <IT,WT> >(clusters);
        freeArray<WT>(maxCoordinates);
        freeArray<WT>(minCoordinates);
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
                numClusters ((1 << dim_) + 1),
                required_elements(required_elements_)
    {
        this->clusters  = allocMemory<KMeansCluster <IT,WT> >(this->numClusters);
        //set dimension and the number of required elements for all clusters.
        for (int i = 0; i < numClusters; ++i){
            this->clusters[i].setParams(this->dim, this->required_elements);
        }

        this->maxCoordinates = allocMemory <WT> (this->dim);
        this->minCoordinates = allocMemory <WT> (this->dim);

        //obtain the min and max coordiantes for each dimension.
        for (int j = 0; j < dim; ++j){
            this->minCoordinates[j] = this->maxCoordinates[j] = this->elementCoords[j][0];
            for(IT i = 1; i < numElements; ++i){
                WT t = this->elementCoords[j][i];
                if(t > this->maxCoordinates[j]){
                    this->maxCoordinates[j] = t;
                }
                if (t < minCoordinates[j]){
                    this->minCoordinates[j] = t;
                }
            }
        }


        //assign initial cluster centers.
        for (int j = 0; j < dim; ++j){
            int mod = (1 << (j+1));
            for (int i = 0; i < numClusters - 1; ++i){
                WT c = 0;
                if ( (i % mod) < mod / 2){
                    c = this->maxCoordinates[j];
                    //cout << "i:" << i << " j:" << j << " setting max:" << c << endl;
                }
                else {
                    c = this->minCoordinates[j];
                }
                this->clusters[i].center[j] = c;
            }
        }

        //last cluster center is placed to middle.
        for (int j = 0; j < dim; ++j){
            this->clusters[numClusters - 1].center[j] = (this->maxCoordinates[j] + this->minCoordinates[j]) / 2;
        }


        /*
        for (int i = 0; i < numClusters; ++i){
            //cout << endl << "cluster:" << i << endl << "\t";
            for (int j = 0; j < dim; ++j){
                cout << this->clusters[i].center[j] << " ";
            }
        }
         */
    }

    //performs kmeans clustering of coordinates.
    void kmeans(){
        for(int it = 0; it < 10; ++it){
            //cout << "it:" << it << endl;
            for (IT j = 0; j < this->numClusters; ++j){
                this->clusters[j].clearHeap();
            }
            for (IT i = 0; i < this->numElements; ++i){
                //cout << "i:" << i << " numEl:" << this->numElements << endl;
                for (IT j = 0; j < this->numClusters; ++j){
                    //cout << "j:" << j << " numClusters:" << this->numClusters << endl;
                    this->clusters[j].getDistance(i,this->elementCoords);
                }
            }
            bool moved = false;
            for (IT j = 0; j < this->numClusters; ++j){
                moved =(this->clusters[j].getNewCenters(this->elementCoords) || moved );
            }
            if (!moved){
                break;
            }
        }


    }

    //finds the cluster in which the coordinates are the closest to each other.
    void getMinDistanceCluster(IT *procPermutation){

        WT minDistance = this->clusters[0].getDistanceToCenter();
        IT minCluster = 0;
        //cout << "j:" << 0 << " minDistance:" << minDistance << " minTmpDistance:" << minDistance<< " minCluster:" << minCluster << endl;
        for (IT j = 1; j < this->numClusters; ++j){
            WT minTmpDistance = this->clusters[j].getDistanceToCenter();
            //cout << "j:" << j << " minDistance:" << minDistance << " minTmpDistance:" << minTmpDistance<< " minCluster:" << minCluster << endl;
            if(minTmpDistance < minDistance){
                minDistance = minTmpDistance;
                minCluster = j;
            }
        }

        //cout << "minCluster:" << minCluster << endl;
        this->clusters[minCluster].copyCoordinates(procPermutation);
    }
};



#define MINOF(a,b) (((a)<(b))?(a):(b))

/*! \brief fillContinousArray function
 *
 *  \param arr   array to be filled in with values.
 *  \param arrSize the size of the array.
 *  \param val    the pointer to the value to be filled. if given NULL, the filling performs arr[i] = i.
 */
template <typename T>
void fillContinousArray(T *arr, size_t arrSize, T *val){
    if(val == NULL){

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(size_t i = 0; i < arrSize; ++i){
            arr[i] = i;
        }

    }
    else {
        T v = *val;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for(size_t i = 0; i < arrSize; ++i){
            //cout << "writing to i:" << i << " arr:" << arrSize << endl;
            arr[i] = v;
        }
    }
}

/*! \brief CommunicationModel Base Class that performs mapping between the coordinate partitioning result.
 */
template <typename procId_t, typename pcoord_t>
class CommunicationModel{
protected:
    double commCost;
public:

    procId_t no_procs; //the number of processors
    procId_t no_tasks;  //the number of taks.
    CommunicationModel(): commCost(),no_procs(0), no_tasks(0){}
    CommunicationModel(procId_t no_procs_, procId_t no_tasks_):
        commCost(),
        no_procs(no_procs_),
        no_tasks(no_tasks_){}
    virtual ~CommunicationModel(){}
    procId_t getNProcs() const{
        return this->no_procs;
    }
    procId_t getNTasks()const{
        return this->no_tasks;
    }
    /*
    void calculateCommunicationCost2(
            procId_t *task_to_proc,
            procId_t *task_communication_xadj,
            procId_t *task_communication_adj){

        double totalCost = 0;

        procId_t commCount = 0;
        for (procId_t task = 0; task < this->no_tasks; ++task){
            int assigned_proc = task_to_proc[task];
            procId_t task_adj_begin = 0;
            //cout << "task:" << task << endl;
            procId_t task_adj_end = task_communication_xadj[task];
            if (task > 0) task_adj_begin = task_communication_xadj[task - 1];

            commCount += task_adj_end - task_adj_begin;
            //cout << "task:" << task << " proc:" << assigned_proc << endl;
            for (procId_t task2 = task_adj_begin; task2 < task_adj_end; ++task2){

                //cout << "task2:" << task2 << endl;
                procId_t neighborTask = task_communication_adj[task2];
                //cout << "neighborTask :" << neighborTask  << endl;
                int neighborProc = task_to_proc[neighborTask];
                double distance = getProcDistance(assigned_proc, neighborProc);
                cout << assigned_proc << " " << neighborProc << " " << distance << endl;

                totalCost += distance ;
            }
        }

        this->commCost = totalCost;// commCount;
    }
    */

    void calculateCommunicationCost(
            procId_t *task_to_proc,
            procId_t *task_communication_xadj,
            procId_t *task_communication_adj,
            pcoord_t *task_communication_edge_weight){

        double totalCost = 0;

        procId_t commCount = 0;
        for (procId_t task = 0; task < this->no_tasks; ++task){
            int assigned_proc = task_to_proc[task];
            procId_t task_adj_begin = 0;
            //cout << "task:" << task << endl;
            procId_t task_adj_end = task_communication_xadj[task];
            if (task > 0) task_adj_begin = task_communication_xadj[task - 1];

            commCount += task_adj_end - task_adj_begin;
            //cout << "task:" << task << " proc:" << assigned_proc << endl;
            for (procId_t task2 = task_adj_begin; task2 < task_adj_end; ++task2){

                //cout << "task2:" << task2 << endl;
                procId_t neighborTask = task_communication_adj[task2];
                //cout << "neighborTask :" << neighborTask  << endl;
                int neighborProc = task_to_proc[neighborTask];
                double distance = getProcDistance(assigned_proc, neighborProc);
                //cout << "assigned_proc:" << assigned_proc << " neighborProc:" << neighborProc << " d:" << distance << endl;
  

                if (task_communication_edge_weight == NULL){
                  totalCost += distance ;
                } 
                else {
                  totalCost += distance * task_communication_edge_weight[task2];
                }
            }
        }

        this->commCost = totalCost;// commCount;
    }

    double getCommunicationCostMetric(){
        return this->commCost;
    }

    virtual double getProcDistance(int procId1, int procId2) const = 0;

    /*! \brief Function is called whenever nprocs > no_task.
     * Function returns only the subset of processors that are closest to each other.
     *  \param proc_to_task_xadj holds the pointer to the task array
     *  \param proc_to_task_xadj holds the indices of tasks wrt to proc_to_task_xadj array.
     *  \param task_to_proc holds the processors mapped to tasks.
     */
    virtual void getMapping(
            int myRank,
            RCP<const Environment> env,
            ArrayRCP <procId_t> &proc_to_task_xadj, //  = allocMemory<procId_t> (this->no_procs); //holds the pointer to the task array
            ArrayRCP <procId_t> &proc_to_task_adj, // = allocMemory<procId_t>(this->no_tasks); //holds the indices of tasks wrt to proc_to_task_xadj array.
            ArrayRCP <procId_t> &task_to_proc //allocMemory<procId_t>(this->no_tasks); //holds the processors mapped to tasks.
    ) const = 0;
};
/*! \brief CoordinateModelInput Class that performs mapping between the coordinate partitioning result and mpi ranks
 * base on the coordinate results and mpi physical coordinates.
 */
template <typename pcoord_t,  typename tcoord_t, typename procId_t>
class CoordinateCommunicationModel:public CommunicationModel<procId_t, pcoord_t> {
public:
    //private:
    int proc_coord_dim; //dimension of the processors
    pcoord_t **proc_coords; //the processor coordinates. allocated outside of the class.
    int task_coord_dim; //dimension of the tasks coordinates.
    tcoord_t **task_coords; //the task coordinates allocated outside of the class.
    int partArraySize;
    procId_t *partNoArray;

    //public:
    CoordinateCommunicationModel():
        CommunicationModel<procId_t, pcoord_t>(),
        proc_coord_dim(0),
        proc_coords(0),
        task_coord_dim(0),
        task_coords(0),
        partArraySize(-1),
        partNoArray(NULL){}

    virtual ~CoordinateCommunicationModel(){}

    /*! \brief Class Constructor:
     *  \param pcoord_dim_ the dimension of the processors
     *  \param pcoords_   the processor coordinates. allocated outside of the class.
     *  \param tcoord_dim_   dimension of the tasks coordinates.
     *  \param tcoords_   the task coordinates allocated outside of the class.
     *  \param no_procs_   the number of processors
     *  \param no_tasks_   the number of taks.
     */
    CoordinateCommunicationModel(
            int pcoord_dim_,
            pcoord_t **pcoords_,
            int tcoord_dim_,
            tcoord_t **tcoords_,
            procId_t no_procs_,
            procId_t no_tasks_
            ):
                CommunicationModel<procId_t, pcoord_t>(no_procs_, no_tasks_),
                proc_coord_dim(pcoord_dim_), proc_coords(pcoords_),
                task_coord_dim(tcoord_dim_), task_coords(tcoords_),
                partArraySize(min(tcoord_dim_, pcoord_dim_)),
                partNoArray(NULL){
    }


    void setPartArraySize(int psize){
        this->partArraySize = psize;
    }
    void setPartArray(procId_t *pNo){
        this->partNoArray = pNo;
    }

    /*! \brief Function is called whenever nprocs > no_task.
     * Function returns only the subset of processors that are closest to each other.
     *  \param proc_permutation holds the indices of the processors that are chosen.
     *  \param nprocs the number of processors.
     *  \param ntasks the number of taks.
     */
    void getClosestSubset(procId_t *proc_permutation, procId_t nprocs, procId_t ntasks) const{
        //currently returns a random subset.

        procId_t minCoordDim = MINOF(this->task_coord_dim, this->proc_coord_dim);
        KMeansAlgorithm<procId_t, pcoord_t > kma(
                minCoordDim, nprocs,
                this->proc_coords, ntasks);

        kma.kmeans();
        kma.getMinDistanceCluster(proc_permutation);

        for(int i = ntasks; i < nprocs; ++i){
            proc_permutation[i] = -1;
        }
        /*
        //fill array.
        fillContinousArray<procId_t>(proc_permutation, nprocs, NULL);
        int _u_umpa_seed = 847449649;
        srand (time(NULL));
        int a = rand() % 1000 + 1;
        _u_umpa_seed -= a;
        //permute array randomly.
        update_visit_order(proc_permutation, nprocs,_u_umpa_seed, 1);
         */
    }

    //temporary, necessary for random permutation.
    static procId_t umpa_uRandom(procId_t l, int &_u_umpa_seed)
    {
        int   a = 16807;
        int   m = 2147483647;
        int   q = 127773;
        int   r = 2836;
        int   lo, hi, test;
        double d;

        lo = _u_umpa_seed % q;
        hi = _u_umpa_seed / q;
        test = (a*lo)-(r*hi);
        if (test>0)
            _u_umpa_seed = test;
        else
            _u_umpa_seed = test + m;
        d = (double) ((double) _u_umpa_seed / (double) m);
        return (procId_t) (d*(double)l);
    }

    virtual double getProcDistance(int procId1, int procId2) const{
        double distance = 0;
        for (int i = 0 ; i < this->proc_coord_dim; ++i){
            distance += ABS(proc_coords[i][procId1] - proc_coords[i][procId2]);
        }
        return distance;
    }


    //temporary, does random permutation.
    void update_visit_order(procId_t* visitOrder, procId_t n, int &_u_umpa_seed, procId_t rndm) {
        procId_t *a = visitOrder;


        if (rndm){
            procId_t i, u, v, tmp;

            if (n <= 4)
                return;

            //srand ( time(NULL) );

            //_u_umpa_seed = _u_umpa_seed1 - (rand()%100);
            for (i=0; i<n; i+=16)
            {
                u = umpa_uRandom(n-4, _u_umpa_seed);
                v = umpa_uRandom(n-4, _u_umpa_seed);
                SWAP(a[v], a[u], tmp);
                SWAP(a[v+1], a[u+1], tmp);
                SWAP(a[v+2], a[u+2], tmp);
                SWAP(a[v+3], a[u+3], tmp);
            }
        }
        else {
            procId_t i, end = n / 4;

            for (i=1; i<end; i++)
            {
                procId_t j=umpa_uRandom(n-i, _u_umpa_seed);
                procId_t t=a[j];
                a[j] = a[n-i];
                a[n-i] = t;
            }
        }
        //PermuteInPlace(visitOrder, n);
    }


    /*! \brief Function is called whenever nprocs > no_task.
     * Function returns only the subset of processors that are closest to each other.
     *  \param proc_to_task_xadj holds the pointer to the task array
     *  \param proc_to_task_xadj holds the indices of tasks wrt to proc_to_task_xadj array.
     *  \param task_to_proc holds the processors mapped to tasks.
     */
    virtual void getMapping(
            int myRank,
            RCP<const Environment> env,
            ArrayRCP <procId_t> &rcp_proc_to_task_xadj, //  = allocMemory<procId_t> (this->no_procs); //holds the pointer to the task array
            ArrayRCP <procId_t> &rcp_proc_to_task_adj, // = allocMemory<procId_t>(this->no_tasks); //holds the indices of tasks wrt to proc_to_task_xadj array.
            ArrayRCP <procId_t> &rcp_task_to_proc //allocMemory<procId_t>(this->no_tasks); //holds the processors mapped to tasks.
    ) const{

        rcp_proc_to_task_xadj = ArrayRCP <procId_t> (this->no_procs);
        rcp_proc_to_task_adj = ArrayRCP <procId_t> (this->no_tasks);
        rcp_task_to_proc = ArrayRCP <procId_t> (this->no_tasks);

        procId_t *proc_to_task_xadj = rcp_proc_to_task_xadj.getRawPtr(); //holds the pointer to the task array
        procId_t *proc_to_task_adj = rcp_proc_to_task_adj.getRawPtr(); //holds the indices of tasks wrt to proc_to_task_xadj array.
        procId_t *task_to_proc = rcp_task_to_proc.getRawPtr(); //holds the processors mapped to tasks.);


        procId_t invalid = 0;
        fillContinousArray<procId_t> (proc_to_task_xadj, this->no_procs, &invalid);

        //obtain the number of parts that should be divided.
        procId_t num_parts = MINOF(this->no_procs, this->no_tasks);
        //obtain the min coordinate dim.
        procId_t minCoordDim = MINOF(this->task_coord_dim, this->proc_coord_dim);

        int recursion_depth = partArraySize;
        if(partArraySize < minCoordDim) recursion_depth = minCoordDim;

        int taskPerm = z2Fact<int>(this->task_coord_dim); //get the number of different permutations for task dimension ordering
        int procPerm = z2Fact<int>(this->proc_coord_dim); //get the number of different permutations for proc dimension ordering
        int permutations =  taskPerm * procPerm; //total number of permutations

        //holds the pointers to proc_adjList
        procId_t *proc_xadj = allocMemory<procId_t> (num_parts);
        //holds the processors in parts according to the result of partitioning algorithm.
        //the processors assigned to part x is at proc_adjList[ proc_xadj[x - 1] : proc_xadj[x] ]
        procId_t *proc_adjList = allocMemory<procId_t>(this->no_procs);


        procId_t used_num_procs = this->no_procs;
        if(this->no_procs > this->no_tasks){
            //obtain the subset of the processors that are closest to each other.
            this->getClosestSubset(proc_adjList, this->no_procs, this->no_tasks);
            used_num_procs = this->no_tasks;
        }
        else {
            fillContinousArray<procId_t>(proc_adjList,this->no_procs, NULL);
        }

        int myPermutation = myRank % permutations; //the index of the permutation

        int myProcPerm=  myPermutation % procPerm; // the index of the proc permutation
        int myTaskPerm  = myPermutation / procPerm; // the index of the task permutation

        int *permutation = allocMemory<int> ((this->proc_coord_dim > this->task_coord_dim) 
                                                   ? this->proc_coord_dim : this->task_coord_dim);

        //get the permutation order from the proc permutation index.
        ithPermutation<int>(this->proc_coord_dim, myProcPerm, permutation);
        //reorder the coordinate dimensions.
        pcoord_t **pcoords = allocMemory<pcoord_t *> (this->proc_coord_dim);
        for(int i = 0; i < this->proc_coord_dim; ++i){
            pcoords[i] = this->proc_coords[permutation[i]];
            //cout << permutation[i] << " ";
        }


        //do the partitioning and renumber the parts.
        env->timerStart(MACRO_TIMERS, "Mapping - Proc Partitioning");
        sequentialTaskPartitioning<pcoord_t, procId_t, procId_t>(
                env,
                this->no_procs,
                used_num_procs,
                num_parts,
                minCoordDim,
                pcoords,//this->proc_coords,
                proc_adjList,
                proc_xadj,
                recursion_depth,
                partNoArray
                //,"proc_partitioning"
        );
        env->timerStop(MACRO_TIMERS, "Mapping - Proc Partitioning");
        freeArray<pcoord_t *> (pcoords);


        procId_t *task_xadj = allocMemory<procId_t> (num_parts);
        procId_t *task_adjList = allocMemory<procId_t>(this->no_tasks);
        //fill task_adjList st: task_adjList[i] <- i.
        fillContinousArray<procId_t>(task_adjList,this->no_tasks, NULL);

        //get the permutation order from the task permutation index.
        ithPermutation<int>(this->task_coord_dim, myTaskPerm, permutation);

        //reorder task coordinate dimensions.
        tcoord_t **tcoords = allocMemory<tcoord_t *> (this->task_coord_dim);
        for(int i = 0; i < this->task_coord_dim; ++i){
            tcoords[i] = this->task_coords[permutation[i]];
        }

        env->timerStart(MACRO_TIMERS, "Mapping - Task Partitioning");
        //partitioning of tasks
        sequentialTaskPartitioning<tcoord_t, procId_t, procId_t>(
                env,
                this->no_tasks,
                this->no_tasks,
                num_parts,
                minCoordDim,
                tcoords, //this->task_coords,
                task_adjList,
                task_xadj,
                recursion_depth,
                partNoArray
                //,"task_partitioning"
        );
        env->timerStop(MACRO_TIMERS, "Mapping - Task Partitioning");
        freeArray<pcoord_t *> (tcoords);
        freeArray<int> (permutation);


        //filling proc_to_task_xadj, proc_to_task_adj, task_to_proc arrays.
        for(procId_t i = 0; i < num_parts; ++i){

            procId_t proc_index_begin = 0;
            procId_t task_begin_index = 0;

            if (i > 0) {
                proc_index_begin = proc_xadj[i - 1];
                task_begin_index = task_xadj[i - 1];
            }
            procId_t proc_index_end = proc_xadj[i];
            procId_t task_end_index = task_xadj[i];


            if(proc_index_end - proc_index_begin != 1){
                cerr << "Error at partitioning of processors" << endl;
                cerr << "PART:" << i << " is assigned to " << proc_index_end - proc_index_begin << " processors." << std::endl;
                exit(1);
            }
            procId_t assigned_proc = proc_adjList[proc_index_begin];
            proc_to_task_xadj[assigned_proc] = task_end_index - task_begin_index;
        }


        //holds the pointer to the task array
        procId_t *proc_to_task_xadj_work = allocMemory<procId_t> (this->no_procs);
        proc_to_task_xadj_work[0] = proc_to_task_xadj[0];
        for(procId_t i = 1; i < this->no_procs; ++i){
            proc_to_task_xadj[i] += proc_to_task_xadj[i - 1];
            proc_to_task_xadj_work[i] = proc_to_task_xadj[i];
        }

        for(procId_t i = 0; i < num_parts; ++i){

            procId_t proc_index_begin = 0;
            procId_t task_begin_index = 0;

            if (i > 0) {
                proc_index_begin = proc_xadj[i - 1];
                task_begin_index = task_xadj[i - 1];
            }
            procId_t task_end_index = task_xadj[i];

            procId_t assigned_proc = proc_adjList[proc_index_begin];

            for (procId_t j = task_begin_index; j < task_end_index; ++j){
                procId_t taskId = task_adjList[j];

                task_to_proc[taskId] = assigned_proc;

                proc_to_task_adj [ --proc_to_task_xadj_work[assigned_proc] ] = taskId;
            }
        }

        freeArray<procId_t>(proc_to_task_xadj_work);
        freeArray<procId_t>(task_xadj);
        freeArray<procId_t>(task_adjList);
        freeArray<procId_t>(proc_xadj);
        freeArray<procId_t>(proc_adjList);
    }

};

template <typename Adapter, typename procId_t>
class CoordinateTaskMapper:public PartitionMapping<Adapter>{
protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

    typedef typename Adapter::scalar_t pcoord_t;
    typedef typename Adapter::scalar_t tcoord_t;

#endif

    //RCP<const Environment> env;
    ArrayRCP<procId_t> proc_to_task_xadj; //  = allocMemory<procId_t> (this->no_procs); //holds the pointer to the task array
    ArrayRCP<procId_t> proc_to_task_adj; // = allocMemory<procId_t>(this->no_tasks); //holds the indices of tasks wrt to proc_to_task_xadj array.
    ArrayRCP<procId_t> task_to_proc; //allocMemory<procId_t>(this->no_procs); //holds the processors mapped to tasks.
    bool isOwnerofModel;
    CoordinateCommunicationModel<pcoord_t,tcoord_t,procId_t> *proc_task_comm;
    procId_t nprocs;
    procId_t ntasks;
    ArrayRCP<procId_t>task_communication_xadj;
    ArrayRCP<procId_t>task_communication_adj;


public:

    void getProcTask(procId_t* &proc_to_task_xadj_, procId_t* &proc_to_task_adj_){
        proc_to_task_xadj_ = this->proc_to_task_xadj.getRawPtr();
        proc_to_task_adj_ = this->proc_to_task_adj.getRawPtr();
    }


    virtual ~CoordinateTaskMapper(){
        //freeArray<procId_t> (proc_to_task_xadj);
        //freeArray<procId_t> (proc_to_task_adj);
        //freeArray<procId_t> (task_to_proc);
        if(this->isOwnerofModel){
            delete this->proc_task_comm;
        }
    }
    /*! \brief Constructor.
     * When this constructor is called, in order to calculate the communication metric,
     * the task adjacency graph is created based on the coordinate model input and partitioning of it.
     * if the communication graph is already calculated, use the other constructors.
     *  \param comm_ is the communication object.
     *  \param machine_ is the machineRepresentation object. Stores the coordinates of machines.
     *  \param model_ is the input adapter.
     *  \param soln_ is the solution object. Holds the assignment of points.
     *  \param envConst_ is the environment object.
     */
    CoordinateTaskMapper(
            const Teuchos::Comm<int> *comm_,
            const MachineRepresentation<pcoord_t> *machine_,
            const Zoltan2::Model<typename Adapter::base_adapter_t> *model_,
            const Zoltan2::PartitioningSolution<Adapter> *soln_,
            const Environment *envConst
    ):  PartitionMapping<Adapter> (comm_, machine_, model_, soln_, envConst),
            proc_to_task_xadj(0),
            proc_to_task_adj(0),
            task_to_proc(0),
            isOwnerofModel(true),
            proc_task_comm(0),
            task_communication_xadj(0),
            task_communication_adj(0){

        pcoord_t *task_communication_edge_weight_ = NULL;
        //if mapping type is 0 then it is coordinate mapping
        int procDim = machine_->getProcDim();
        this->nprocs = machine_->getNumProcs();
        //get processor coordinates.
        pcoord_t **procCoordinates = machine_->getProcCoords();

        int coordDim = ((Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *)model_)->getCoordinateDim();
        this->ntasks = soln_->getActualGlobalNumberOfParts();
        if (procId_t (soln_->getTargetGlobalNumberOfParts()) > this->ntasks){
            this->ntasks = soln_->getTargetGlobalNumberOfParts();
        }
        //cout << "actual: " << this->ntasks << endl;

        //alloc memory for part centers.
        tcoord_t **partCenters = NULL;
        partCenters = allocMemory<tcoord_t *>(coordDim);
        for (int i = 0; i < coordDim; ++i){
            partCenters[i] = allocMemory<tcoord_t>(this->ntasks);
        }


        envConst->timerStart(MACRO_TIMERS, "Mapping - Solution Center");
        //get centers for the parts.
        getSolutionCenterCoordinates<Adapter, typename Adapter::scalar_t,procId_t>(
                envConst,
                comm_,
                ((Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *)model_),
                this->soln,
                coordDim,
                ntasks,
                partCenters);

        envConst->timerStop(MACRO_TIMERS, "Mapping - Solution Center");


        //create coordinate communication model.
        this->proc_task_comm =
                new Zoltan2::CoordinateCommunicationModel<pcoord_t,tcoord_t,procId_t>(
                        procDim,
                        procCoordinates,
                        coordDim,
                        partCenters,
                        this->nprocs,
                        this->ntasks
                );

        int myRank = comm_->getRank();


        envConst->timerStart(MACRO_TIMERS, "Mapping - Processor Task map");
        this->doMapping(myRank);
        envConst->timerStop(MACRO_TIMERS, "Mapping - Processor Task map");


        envConst->timerStart(MACRO_TIMERS, "Mapping - Communication Graph");
        ((Zoltan2::PartitioningSolution<Adapter> *)soln_)->getCommunicationGraph(
                comm_,
                task_communication_xadj,
                task_communication_adj
         );

        /*
        if (myRank == 0){
            procId_t maxComm = 0;
            procId_t totalComm = 0;
            for (procId_t i = 0; i < this->ntasks;++i){
                procId_t tBegin = 0;
                if (i > 0){
                    tBegin = task_communication_xadj[i - 1];
                }
                procId_t tEnd = task_communication_xadj[i];
                if (tEnd - tBegin > maxComm) maxComm = tEnd - tBegin;
            }
            totalComm = task_communication_xadj[this->ntasks - 1];
            cout << "commMax:" << maxComm << " totalComm:"<< totalComm << endl;
        }
        */




        envConst->timerStop(MACRO_TIMERS, "Mapping - Communication Graph");
#ifdef gnuPlot
        if (comm_->getRank() == 0){

            procId_t taskCommCount = task_communication_xadj.size();
            std::cout << " TotalComm:" << task_communication_xadj[taskCommCount - 1] << std::endl;
            procId_t maxN = task_communication_xadj[0];
            for (procId_t i = 1; i < taskCommCount; ++i){
                procId_t nc = task_communication_xadj[i] - task_communication_xadj[i - 1];
                if (maxN < nc) maxN = nc;
            }
            std::cout << " maxNeighbor:" << maxN << std::endl;
        }

        this->writeGnuPlot(comm_, soln_, coordDim, partCenters);
        /*
        std::string file = "gggnuPlot";
        std::string exten = ".plot";
        ofstream mm("2d.txt");
        file += toString<int>(comm_->getRank()) + exten;
        std::ofstream ff(file.c_str());
        //ff.seekg (0, ff.end);
        RCP < vector <Zoltan2::coordinateModelPartBox <scalar_t, partId_t> > > outPartBoxes = ((Zoltan2::PartitioningSolution<Adapter> *)soln_)->getPartBoxes();

        for (partId_t i = 0; i < this->ntasks;++i){
            (*outPartBoxes)[i].writeGnuPlot(ff, mm);
        }
        if (coordDim == 2){
        ff << "plot \"2d.txt\"" << endl;
        //ff << "\n pause -1" << endl;
        }
        else {
            ff << "splot \"2d.txt\"" << endl;
            //ff << "\n pause -1" << endl;
        }
        mm.close();

        ff << "set style arrow 5 nohead size screen 0.03,15,135 ls 1" << endl;
        for (partId_t i = 0; i < this->ntasks;++i){
            procId_t pb = 0;
            if (i > 0) pb = task_communication_xadj[i -1];
            procId_t pe = task_communication_xadj[i];
            for (procId_t p = pb; p < pe; ++p){
                procId_t n = task_communication_adj[p];

                //cout << "i:" << i << " n:" << n << endl;
                std::string arrowline = "set arrow from ";
                for (int j = 0; j < coordDim - 1; ++j){
                    arrowline += toString<scalar_t>(partCenters[j][n]) + ",";
                }
                arrowline += toString<scalar_t>(partCenters[coordDim -1][n]) + " to ";


                for (int j = 0; j < coordDim - 1; ++j){
                    arrowline += toString<scalar_t>(partCenters[j][i]) + ",";
                }
                arrowline += toString<scalar_t>(partCenters[coordDim -1][i]) + " as 5\n";

                //cout << "arrow:" << arrowline << endl;
                ff << arrowline;
            }
        }

        ff << "replot\n pause -1" << endl;
        ff.close();
        */
#endif

        envConst->timerStart(MACRO_TIMERS, "Mapping - Communication Cost");
        this->proc_task_comm->calculateCommunicationCost(
                task_to_proc.getRawPtr(),
                task_communication_xadj.getRawPtr(),
                task_communication_adj.getRawPtr(),
                task_communication_edge_weight_
        );

        //cout << "me: " << comm_->getRank() << " cost:" << this->proc_task_comm->getCommunicationCostMetric() << endl;

        envConst->timerStop(MACRO_TIMERS, "Mapping - Communication Cost");

        //cout << "me: " << comm_->getRank() << " cost:" << this->proc_task_comm->getCommunicationCostMetric() << endl;
        //processors are divided into groups of size numProc! * numTasks!
        //each processor in the group obtains a mapping with a different rotation
        //and best one is broadcasted all processors.
        this->getBestMapping();
#ifdef gnuPlot
        this->writeMapping2(comm_->getRank());
#endif

        for (int i = 0; i < coordDim; ++i){
            freeArray<tcoord_t>(partCenters[i]);
        }
        freeArray<tcoord_t *>(partCenters);

    }

    /*! \brief Constructor
     * When this constructor is called, it is assumed that the communication graph is already computed.
     * It is assumed that the communication graph is given in task_communication_xadj_ and task_communication_adj_ parameters.
     *  \param machine_ is the machineRepresentation object. Stores the coordinates of machines.
     *  \param model_ is the input adapter.
     *  \param soln_ is the solution object. Holds the assignment of points.
     *  \param envConst_ is the environment object.
     *  \param task_communication_xadj_ 
     *  \param task_communication_adj_ 
     */
    CoordinateTaskMapper(
            const Teuchos::Comm<int> *comm_,
            const MachineRepresentation<pcoord_t> *machine_,
            const Zoltan2::Model<typename Adapter::base_adapter_t> *model_,
            const Zoltan2::PartitioningSolution<Adapter> *soln_,
            const Environment *envConst_,
            ArrayRCP<procId_t>task_communication_xadj_,
            ArrayRCP<procId_t>task_communication_adj_,
            pcoord_t *task_communication_edge_weight_
    ):  PartitionMapping<Adapter> (comm_, machine_, model_, soln_, envConst_),
            proc_to_task_xadj(0),
            proc_to_task_adj(0),
            task_to_proc(0),
            isOwnerofModel(true),
            proc_task_comm(0),
            task_communication_xadj(task_communication_xadj_),
            task_communication_adj(task_communication_adj_){

        //if mapping type is 0 then it is coordinate mapping
        int procDim = machine_->getProcDim();
        this->nprocs = machine_->getNumProcs();
        //get processor coordinates.
        pcoord_t **procCoordinates = machine_->getProcCoords();

        int coordDim = ((Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *)model_)->getCoordinateDim();
        this->ntasks = soln_->getActualGlobalNumberOfParts();
        //cout << "actual: " << this->ntasks << endl;

        //alloc memory for part centers.
        tcoord_t **partCenters = NULL;
        partCenters = allocMemory<tcoord_t *>(coordDim);
        for (int i = 0; i < coordDim; ++i){
            partCenters[i] = allocMemory<tcoord_t>(this->ntasks);
        }
        //get centers for the parts.
        getSolutionCenterCoordinates<Adapter, typename Adapter::scalar_t,procId_t>(
                envConst_,
                comm_,
                ((Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *)model_),
                this->soln,
                coordDim,
                ntasks,
                partCenters);

        //create coordinate communication model.
        this->proc_task_comm =
                new Zoltan2::CoordinateCommunicationModel<pcoord_t,tcoord_t,procId_t>(
                        procDim,
                        procCoordinates,
                        coordDim,
                        partCenters,
                        this->nprocs,
                        this->ntasks
                );

        int myRank = comm_->getRank();
        this->doMapping(myRank);
        this->proc_task_comm->calculateCommunicationCost(
                task_to_proc.getRawPtr(),
                task_communication_xadj.getRawPtr(),
                task_communication_adj.getRawPtr(),
                task_communication_edge_weight_
                );


        //cout << "me: " << comm_->getRank() << " cost:" << this->proc_task_comm->getCommunicationCostMetric() << endl;

        //processors are divided into groups of size numProc! * numTasks!
        //each processor in the group obtains a mapping with a different rotation
        //and best one is broadcasted all processors.
        this->getBestMapping();
#ifdef gnuPlot
        this->writeMapping2(comm_->getRank());
#endif

        for (int i = 0; i < coordDim; ++i){
            freeArray<tcoord_t>(partCenters[i]);
        }
        freeArray<tcoord_t *>(partCenters);
    }



    /*! \brief Constructor
     * Constructor called for fortran interface.
     *  \param envConst_ is the environment object.
     *  \param task_communication_xadj_
     *  \param task_communication_adj_
     */
    CoordinateTaskMapper(
            const Teuchos::Comm<int> *comm_,
            int procDim,
            int numProcessors,
            pcoord_t **machine_coords_,

            int taskDim,
            procId_t numTasks,
            tcoord_t **task_coords,
            const Environment *envConst_,
            ArrayRCP<procId_t>task_communication_xadj_,
            ArrayRCP<procId_t>task_communication_adj_,
            pcoord_t *task_communication_edge_weight_,
            int partArraySize,
            procId_t *partNoArray
    ):  PartitionMapping<Adapter>(comm_, NULL, NULL, NULL, envConst_),
            proc_to_task_xadj(0),
            proc_to_task_adj(0),
            task_to_proc(0),
            isOwnerofModel(true),
            proc_task_comm(0),
            task_communication_xadj(task_communication_xadj_),
            task_communication_adj(task_communication_adj_){

        //if mapping type is 0 then it is coordinate mapping

        this->nprocs = numProcessors;
        //get processor coordinates.
        pcoord_t **procCoordinates = machine_coords_;

        int coordDim = taskDim;
        this->ntasks = numTasks;

        //alloc memory for part centers.
        tcoord_t **partCenters = task_coords;

        //create coordinate communication model.
        this->proc_task_comm =
                new Zoltan2::CoordinateCommunicationModel<pcoord_t,tcoord_t,procId_t>(
                        procDim,
                        procCoordinates,
                        coordDim,
                        partCenters,
                        this->nprocs,
                        this->ntasks
                );
        this->proc_task_comm->setPartArraySize(partArraySize);
        this->proc_task_comm->setPartArray(partNoArray);

        int myRank = comm_->getRank();

        this->doMapping(myRank);
#ifdef gnuPlot
        this->writeMapping2(myRank);
#endif

        this->proc_task_comm->calculateCommunicationCost(
                task_to_proc.getRawPtr(),
                task_communication_xadj.getRawPtr(),
                task_communication_adj.getRawPtr(),
                task_communication_edge_weight_
                );
        //cout << "me: " << comm_->getRank() << " cost:" << this->proc_task_comm->getCommunicationCostMetric() << endl;

        this->getBestMapping();
	/*
	if (myRank == 0){
        	this->proc_task_comm->calculateCommunicationCost2(
                	task_to_proc.getRawPtr(),
               		task_communication_xadj.getRawPtr(),
	                task_communication_adj.getRawPtr()
                );
	}
	*/


#ifdef gnuPlot
        if(comm_->getRank() == 0)
        this->writeMapping2(-1);
#endif
    }


    /*! \brief Constructor. Sequential Constructor for test.
     *  \param env_ Environment object.
     *  \param proc_task_comm_ is the template parameter for which the mapping will be obtained with getMapping() function.
     */
    CoordinateTaskMapper(
            const Environment *env_,
            CoordinateCommunicationModel<pcoord_t,tcoord_t,procId_t> *proc_task_comm_
    ):PartitionMapping<Adapter>(env_),
        proc_to_task_xadj(0),
        proc_to_task_adj(0),
        task_to_proc(0),
        isOwnerofModel(false),
        proc_task_comm(proc_task_comm_),
        nprocs(proc_task_comm_->getNProcs()),
        ntasks(proc_task_comm_->getNTasks()),
        task_communication_xadj(0),
        task_communication_adj(0){
        this->doMapping(0);
    }

    /*! \brief Constructor. Given the machine object and task centers, performs mapping.
     *  \param env_ Environment object.
     *  \param proc_task_comm_ is the template parameter for which the mapping will be obtained with getMapping() function.
     */
    CoordinateTaskMapper(
            const Environment *env_,
            const Teuchos::Comm<int> *comm_,
            const MachineRepresentation<pcoord_t> *machine_,
            int taskDim,
            procId_t taskCount,
            tcoord_t **taskCoords
    ):PartitionMapping<Adapter>(env_, comm_, machine_),
        proc_to_task_xadj(0),
        proc_to_task_adj(0),
        task_to_proc(0),
        isOwnerofModel(true),
        proc_task_comm(0),
        nprocs(comm_->getSize()),
        ntasks(taskCount),
        task_communication_xadj(0),
        task_communication_adj(0){

        this->proc_task_comm = new CoordinateCommunicationModel<pcoord_t,tcoord_t,procId_t> (
                machine_->getProcDim(),
                machine_->getProcCoords(),
                taskDim, taskCoords,
                this->nprocs, taskCount
                );
        this->doMapping(comm_->getRank());
    }


    /*! \brief doMapping function, calls getMapping function of communicationModel object.
     */
    void doMapping(int myRank){

        if(this->proc_task_comm){
            this->proc_task_comm->getMapping(
                    myRank,
                    Teuchos::RCP<const Environment>(this->env, false),
                    this->proc_to_task_xadj, //  = allocMemory<procId_t> (this->no_procs); //holds the pointer to the task array
                    this->proc_to_task_adj, // = allocMemory<procId_t>(this->no_tasks); //holds the indices of tasks wrt to proc_to_task_xadj array.
                    this->task_to_proc //allocMemory<procId_t>(this->no_procs); //holds the processors mapped to tasks.);
            );
        }
        else {
            std::cerr << "communicationModel is not specified in the Mapper" << endl;
            exit(1);
        }
    }

    double getCommunicationCostMetric(){
        return this->proc_task_comm->getCommCost();
    }

    /*! \brief Returns the number of parts to be assigned to this process.
     */
    virtual size_t getLocalNumberOfParts() const{
        return 0;
    }


    /*! \brief creates and returns the subcommunicator for the processor group.
     */
    RCP<Comm<int> > create_subCommunicatior(){
        int procDim = this->proc_task_comm->proc_coord_dim;
        int taskDim = this->proc_task_comm->task_coord_dim;

        int taskPerm = z2Fact<int>(procDim); //get the number of different permutations for task dimension ordering
        int procPerm = z2Fact<int>(taskDim); //get the number of different permutations for proc dimension ordering
        int idealGroupSize =  taskPerm * procPerm; //total number of permutations

        int myRank = this->comm->getRank();
        int commSize = this->comm->getSize();

        int myGroupIndex = myRank / idealGroupSize;

        int prevGroupBegin = (myGroupIndex - 1)* idealGroupSize;
        if (prevGroupBegin < 0) prevGroupBegin = 0;
        int myGroupBegin = myGroupIndex * idealGroupSize;
        int myGroupEnd = (myGroupIndex + 1) * idealGroupSize;
        int nextGroupEnd = (myGroupIndex + 2)* idealGroupSize;

        if (myGroupEnd > commSize){
            myGroupBegin = prevGroupBegin;
            myGroupEnd = commSize;
        }
        if (nextGroupEnd > commSize){
            myGroupEnd = commSize;
        }
        int myGroupSize = myGroupEnd - myGroupBegin;

        int *myGroup = allocMemory<int>(myGroupSize);
        for (int i = 0; i < myGroupSize; ++i){
            myGroup[i] = myGroupBegin + i;
        }
        //cout << "me:" << myRank << " myGroupBegin:" << myGroupBegin << " myGroupEnd:" << myGroupEnd << endl;

        ArrayView<const partId_t> myGroupView(myGroup, myGroupSize);

        RCP<Comm<int> > subComm = this->comm->createSubcommunicator(myGroupView);
        freeArray<int>(myGroup);
        return subComm;
    }


    /*! \brief finds the lowest cost mapping and broadcasts solution to everyone.
     */
    void getBestMapping(){
        //create the sub group.
        RCP<Comm<int> > subComm = this->create_subCommunicatior();
        //calculate cost.
        double myCost = this->proc_task_comm->getCommunicationCostMetric();
        //cout << "me:" << this->comm->getRank() << " myCost:" << myCost << endl;
        double localCost[2], globalCost[2];

        localCost[0] = myCost;
        localCost[1] = double(subComm->getRank());

        globalCost[1] = globalCost[0] = std::numeric_limits<double>::max();
        Teuchos::Zoltan2_ReduceBestMapping<int,double> reduceBest;
        reduceAll<int, double>(*subComm, reduceBest,
                2, localCost, globalCost);

        int sender = int(globalCost[1]);

        //cout << "me:" << localCost[1] << " localcost:" << localCost[0]<< " bestcost:" << globalCost[0] << endl;
        //cout << "me:" << localCost[1] << " proc:" << globalCost[1] << endl;
        broadcast (*subComm, sender, this->ntasks, this->task_to_proc.getRawPtr());
        broadcast (*subComm, sender, this->nprocs, this->proc_to_task_xadj.getRawPtr());
        broadcast (*subComm, sender, this->ntasks, this->proc_to_task_adj.getRawPtr());
    }

    /*! \brief getAssignedProcForTask function,
     * returns the assigned tasks with the number of tasks.
     *  \param procId procId being queried.
     *  \param numProcs (output), the number of processor the part is assigned to.
     *  \param procs (output), the list of processors assigned to given part..
     */
    virtual void getProcsForPart(partId_t taskId, int &numProcs, int *procs) const{
        numProcs = 1;
        procs = this->task_to_proc.getRawPtr() + taskId;
    }
    /*! \brief getAssignedProcForTask function, returns the assigned processor id for the given task
     *  \param taskId taskId being queried.
     */
    inline procId_t getAssignedProcForTask(procId_t taskId){
        return this->task_to_proc[taskId];
    }
    /*! \brief getAssignedProcForTask function,
     * returns the assigned tasks with the number of tasks.
     *  \param procId procId being queried.
     *  \param numParts (output), the number of parts the processor is assigned to.
     *  \param parts (output), the list of parts assigned to given processor..
     */
    virtual void getPartsForProc(int procId, partId_t &numParts, partId_t *parts) const{

        procId_t task_begin = 0;
        if (procId > 0) task_begin = this->proc_to_task_xadj[procId - 1];
        procId_t taskend = this->proc_to_task_xadj[procId];
        parts = this->proc_to_task_adj.getRawPtr() + task_begin;
        numParts = taskend - task_begin;
    }

    ArrayView<procId_t> getAssignedTaksForProc(procId_t procId){
        procId_t task_begin = 0;
        if (procId > 0) task_begin = this->proc_to_task_xadj[procId - 1];
        procId_t taskend = this->proc_to_task_xadj[procId];

        /*
        cout << "procId_t:" << procId << " taskCount:" << taskend - task_begin << endl;
        for(procId_t i = task_begin; i < taskend; ++i){
            cout << "procId_t:" << procId << " task:" << proc_to_task_adj[i] << endl;
        }
         */
        if (taskend - task_begin > 0){
            ArrayView <procId_t> assignedParts(this->proc_to_task_adj.getRawPtr() + task_begin, taskend - task_begin);
            return assignedParts;
        }
        else {
            ArrayView <procId_t> assignedParts;
            return assignedParts;
        }
    }

    //write mapping to gnuPlot code to visualize.
    void writeMapping(){
        std::ofstream gnuPlotCode ("gnuPlot.plot", std::ofstream::out);

        int mindim = MINOF(proc_task_comm->proc_coord_dim, proc_task_comm->task_coord_dim);
        string ss = "";
        for(procId_t i = 0; i < this->nprocs; ++i){

            std::string procFile = toString<int>(i) + "_mapping.txt";
            if (i == 0){
                gnuPlotCode << "plot \"" << procFile << "\"\n";
            }
            else {
                gnuPlotCode << "replot \"" << procFile << "\"\n";
            }

            std::ofstream inpFile (procFile.c_str(), std::ofstream::out);

            string gnuPlotArrow = "set arrow from ";
            for(int j = 0; j <  mindim; ++j){
                if (j == mindim - 1){
                    inpFile << proc_task_comm->proc_coords[j][i];
                    gnuPlotArrow += toString<float>(proc_task_comm->proc_coords[j][i]);

                }
                else {
                    inpFile << proc_task_comm->proc_coords[j][i] << " ";
                    gnuPlotArrow += toString<float>(proc_task_comm->proc_coords[j][i]) +",";
                }
            }
            gnuPlotArrow += " to ";


            inpFile << std::endl;
            ArrayView<procId_t> a = this->getAssignedTaksForProc(i);
            for(int k = 0; k <  a.size(); ++k){
                int j = a[k];
                //cout << "i:" << i << " j:"
                string gnuPlotArrow2 = gnuPlotArrow;
                for(int z = 0; z <  mindim; ++z){
                    if(z == mindim - 1){

                        //cout << "z:" << z << " j:" <<  j << " " << proc_task_comm->task_coords[z][j] << endl;
                        inpFile << proc_task_comm->task_coords[z][j];
                        gnuPlotArrow2 += toString<float>(proc_task_comm->task_coords[z][j]);
                    }
                    else{
                        inpFile << proc_task_comm->task_coords[z][j] << " ";
                        gnuPlotArrow2 += toString<float>(proc_task_comm->task_coords[z][j]) +",";
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
    void writeMapping2(int myRank){

        std::string rankStr = toString<int>(myRank);
        std::string gnuPlots = "gnuPlot", extentionS = ".plot";
        std::string outF = gnuPlots + rankStr+ extentionS;
        std::ofstream gnuPlotCode ( outF.c_str(), std::ofstream::out);

        CoordinateCommunicationModel<pcoord_t, tcoord_t, procId_t> *tmpproc_task_comm =
                static_cast <CoordinateCommunicationModel<pcoord_t, tcoord_t, procId_t> * > (proc_task_comm);
        int mindim = MINOF(tmpproc_task_comm->proc_coord_dim, tmpproc_task_comm->task_coord_dim);
        string ss = "";
        string procs = "", parts = "";
        for(procId_t i = 0; i < this->nprocs; ++i){

            //inpFile << std::endl;
            ArrayView<procId_t> a = this->getAssignedTaksForProc(i);
            if (a.size() == 0){
                continue;
            }

            //std::ofstream inpFile (procFile.c_str(), std::ofstream::out);

            string gnuPlotArrow = "set arrow from ";
            for(int j = 0; j <  mindim; ++j){
                if (j == mindim - 1){
                    //inpFile << proc_task_comm->proc_coords[j][i];
                    gnuPlotArrow += toString<float>(tmpproc_task_comm->proc_coords[j][i]);
                    procs += toString<float>(tmpproc_task_comm->proc_coords[j][i]);

                }
                else {
                    //inpFile << proc_task_comm->proc_coords[j][i] << " ";
                    gnuPlotArrow += toString<float>(tmpproc_task_comm->proc_coords[j][i]) +",";
                    procs += toString<float>(tmpproc_task_comm->proc_coords[j][i])+ " ";
                }
            }
            procs += "\n";

            gnuPlotArrow += " to ";


            for(int k = 0; k <  a.size(); ++k){
                int j = a[k];
                //cout << "i:" << i << " j:"
                string gnuPlotArrow2 = gnuPlotArrow;
                for(int z = 0; z <  mindim; ++z){
                    if(z == mindim - 1){

                        //cout << "z:" << z << " j:" <<  j << " " << proc_task_comm->task_coords[z][j] << endl;
                        //inpFile << proc_task_comm->task_coords[z][j];
                        gnuPlotArrow2 += toString<float>(tmpproc_task_comm->task_coords[z][j]);
                        parts += toString<float>(tmpproc_task_comm->task_coords[z][j]);
                    }
                    else{
                        //inpFile << proc_task_comm->task_coords[z][j] << " ";
                        gnuPlotArrow2 += toString<float>(tmpproc_task_comm->task_coords[z][j]) +",";
                        parts += toString<float>(tmpproc_task_comm->task_coords[z][j]) + " ";
                    }
                }
                parts += "\n";
                ss += gnuPlotArrow2 + " nohead\n";
                //inpFile << std::endl;
            }
            //inpFile.close();

        }


        std::ofstream procFile ("procPlot.plot", std::ofstream::out);
        procFile << procs << "\n";
        procFile.close();

        std::ofstream partFile ("partPlot.plot", std::ofstream::out);
        partFile << parts<< "\n";
        partFile.close();

        std::ofstream extraProcFile ("allProc.plot", std::ofstream::out);

        for(procId_t j = 0; j < this->nprocs; ++j){
            for(int i = 0; i <  mindim; ++i){
                extraProcFile << tmpproc_task_comm->proc_coords[i][j] <<  " ";
            }
            extraProcFile << endl;
        }

        extraProcFile.close();

        gnuPlotCode << ss;
        if(mindim == 2){
            gnuPlotCode << "plot \"procPlot.plot\" with points pointsize 3\n";
        } else {
            gnuPlotCode << "splot \"procPlot.plot\" with points pointsize 3\n";
        }
        gnuPlotCode << "replot \"partPlot.plot\" with points pointsize 3\n";
        gnuPlotCode << "replot \"allProc.plot\" with points pointsize 0.65\n";
        gnuPlotCode << "\nreplot\n pause -1 \n";
        gnuPlotCode.close();

    }


    void writeGnuPlot(
            const Teuchos::Comm<int> *comm_,
            const Zoltan2::PartitioningSolution<Adapter> *soln_,
            int coordDim,
            tcoord_t **partCenters
            ){
        std::string file = "gggnuPlot";
        std::string exten = ".plot";
        ofstream mm("2d.txt");
        file += toString<int>(comm_->getRank()) + exten;
        std::ofstream ff(file.c_str());
        //ff.seekg (0, ff.end);
        RCP < vector <Zoltan2::coordinateModelPartBox <tcoord_t, partId_t> > > outPartBoxes = ((Zoltan2::PartitioningSolution<Adapter> *)soln_)->getPartBoxes();

        for (partId_t i = 0; i < this->ntasks;++i){
            (*outPartBoxes)[i].writeGnuPlot(ff, mm);
        }
        if (coordDim == 2){
        ff << "plot \"2d.txt\"" << endl;
        //ff << "\n pause -1" << endl;
        }
        else {
            ff << "splot \"2d.txt\"" << endl;
            //ff << "\n pause -1" << endl;
        }
        mm.close();

        ff << "set style arrow 5 nohead size screen 0.03,15,135 ls 1" << endl;
        for (partId_t i = 0; i < this->ntasks;++i){
            procId_t pb = 0;
            if (i > 0) pb = task_communication_xadj[i -1];
            procId_t pe = task_communication_xadj[i];
            for (procId_t p = pb; p < pe; ++p){
                procId_t n = task_communication_adj[p];

                //cout << "i:" << i << " n:" << n << endl;
                std::string arrowline = "set arrow from ";
                for (int j = 0; j < coordDim - 1; ++j){
                    arrowline += toString<tcoord_t>(partCenters[j][n]) + ",";
                }
                arrowline += toString<tcoord_t>(partCenters[coordDim -1][n]) + " to ";


                for (int j = 0; j < coordDim - 1; ++j){
                    arrowline += toString<tcoord_t>(partCenters[j][i]) + ",";
                }
                arrowline += toString<tcoord_t>(partCenters[coordDim -1][i]) + " as 5\n";

                //cout << "arrow:" << arrowline << endl;
                ff << arrowline;
            }
        }

        ff << "replot\n pause -1" << endl;
        ff.close();
    }
};


template <typename procId_t,  typename pcoord_t>
pcoord_t **shiftMachineCoordinates(int machine_dim, procId_t *machine_dimensions, procId_t numProcs, pcoord_t **mCoords){
    pcoord_t **result_machine_coords = NULL;
    result_machine_coords = new pcoord_t*[machine_dim];
    for (int i = 0; i < machine_dim; ++i){
        result_machine_coords[i] = new pcoord_t [numProcs];
    }

    for (int i = 0; i < machine_dim; ++i){
        procId_t numMachinesAlongDim = machine_dimensions[i];
        procId_t *machineCounts= new procId_t[numMachinesAlongDim];
        memset(machineCounts, 0, sizeof(procId_t) *numMachinesAlongDim);

        int *filledCoordinates= new int[numMachinesAlongDim];

        pcoord_t *coords = mCoords[i];
        for(procId_t j = 0; j < numProcs; ++j){
            procId_t mc = (procId_t) coords[j];
            ++machineCounts[mc];
        }

        procId_t filledCoordinateCount = 0;
        for(procId_t j = 0; j < numMachinesAlongDim; ++j){
            if (machineCounts[j] > 0){
                filledCoordinates[filledCoordinateCount++] = j;
            }
        }

        procId_t firstProcCoord = filledCoordinates[0];
        procId_t firstProcCount = machineCounts[firstProcCoord];

        procId_t lastProcCoord = filledCoordinates[filledCoordinateCount - 1];
        procId_t lastProcCount = machineCounts[lastProcCoord];

        procId_t firstLastGap = numMachinesAlongDim - lastProcCoord + firstProcCoord;
        procId_t firstLastGapProc = lastProcCount + firstProcCount;

        procId_t leftSideProcCoord = firstProcCoord;
        procId_t leftSideProcCount = firstProcCount;
        procId_t biggestGap = 0;
        procId_t biggestGapProc = numProcs;

        procId_t shiftBorderCoordinate = -1;
        for(procId_t j = 1; j < filledCoordinateCount; ++j){
            procId_t rightSideProcCoord= filledCoordinates[j];
            procId_t rightSideProcCount = machineCounts[rightSideProcCoord];

            procId_t gap = rightSideProcCoord - leftSideProcCoord;
            procId_t gapProc = rightSideProcCount + leftSideProcCount;

            /* Pick the largest gap in this dimension. Use fewer process on either side
               of the largest gap to break the tie. An easy addition to this would
               be to weight the gap by the number of processes. */
            if (gap > biggestGap || (gap == biggestGap && biggestGapProc > gapProc)){
                shiftBorderCoordinate = rightSideProcCoord;
                biggestGapProc = gapProc;
                biggestGap = gap;
            }
            leftSideProcCoord = rightSideProcCoord;
            leftSideProcCount = rightSideProcCount;
        }


        if (!(biggestGap > firstLastGap || (biggestGap == firstLastGap && biggestGapProc < firstLastGapProc))){
            shiftBorderCoordinate = -1;
        }

/*
        for(procId_t j = 0; j < filledCoordinateCount; ++j){
            cout << "dim:" << i << " coord:" << filledCoordinates[j] ;

            if (filledCoordinates[j] < shiftBorderCoordinate){
                cout << " will be shifted to " << filledCoordinates[j] + numMachinesAlongDim << endl;
            }
            else {
                cout << endl;
            }
        }
*/

        for(procId_t j = 0; j < numProcs; ++j){

            if (coords[j] < shiftBorderCoordinate){
                result_machine_coords[i][j] = coords[j] + numMachinesAlongDim;

            }
            else {
                result_machine_coords[i][j] = coords[j];
            }
            //cout << "I:" << i << "j:" << j << " coord:" << coords[j] << " now:" << result_machine_coords[i][j] << endl;
        }
        delete [] machineCounts;
        delete [] filledCoordinates;
    }

    return result_machine_coords;

}

// KDDKDD TODO:  This interface should go away or move to MiniGhost;
// KDDKDD TODO:  it is only a convenience for the MiniGhost experiments.
template <typename procId_t, typename pcoord_t, typename tcoord_t>
void coordinateTaskMapperInterface(
  RCP<const Teuchos::Comm<int> > comm_,
  int procDim,
  int numProcessors,
  pcoord_t **machine_coords_,
  int taskDim,
  procId_t numTasks,
  tcoord_t **task_coords,
  procId_t *task_communication_xadj_,
  procId_t *task_communication_adj_,
  pcoord_t *task_communication_edge_weight_, /*float-like, same size with task_communication_adj_ weight of the corresponding edge.*/
  procId_t *proc_to_task_xadj, /*output*/
  procId_t *proc_to_task_adj, /*output*/
  int partArraySize,
  procId_t *partNoArray,
  procId_t *machine_dimensions
)
{

    const Environment *envConst_ = new Environment();
    //RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();

    typedef Tpetra::MultiVector<tcoord_t, procId_t,procId_t, KokkosClassic::DefaultNode::DefaultNodeType> tMVector_t;



    //cout << "numProcessors:" << numProcessors << endl;
    //cout << "task_communication_xadj_[numProcessors]:" << task_communication_xadj_[numProcessors - 1] << endl;
    Teuchos::ArrayRCP<procId_t> task_communication_xadj (task_communication_xadj_, 0, numProcessors, false);
    Teuchos::ArrayRCP<procId_t> task_communication_adj (task_communication_adj_, 0, task_communication_xadj_[numProcessors -1 /* KDDKDD OK for MEHMET's ODD LAYOUT; WRONG FOR TRADITIONAL */], false);
    /*
    int machine_dimensions[3];
    machine_dimensions[0] = 17;
    machine_dimensions[1] = 8;
    machine_dimensions[2] = 24;

     */
    pcoord_t ** updatedMachine  = machine_coords_;
    if (machine_dimensions){
        updatedMachine =
            shiftMachineCoordinates <procId_t, pcoord_t>(
                            procDim,
                            machine_dimensions,
                            numProcessors,
                            machine_coords_);
    }
    CoordinateTaskMapper<XpetraMultiVectorAdapter <tMVector_t>, procId_t> *ctm = new CoordinateTaskMapper<XpetraMultiVectorAdapter <tMVector_t>, procId_t>(
                comm_.getRawPtr(),
                procDim,
                numProcessors,
                updatedMachine,//machine_coords_,

                taskDim,
                numTasks,
                task_coords,

                envConst_,
                task_communication_xadj,
                task_communication_adj,
                task_communication_edge_weight_,
                partArraySize,
                partNoArray
        );

    if (machine_dimensions){
        for (int i = 0; i < procDim; ++i){
            delete [] updatedMachine[i];
        }
        delete [] updatedMachine;
    }

    procId_t* proc_to_task_xadj_;
    procId_t* proc_to_task_adj_;

    ctm->getProcTask(proc_to_task_xadj_, proc_to_task_adj_);

    for (procId_t i = 0; i < numProcessors; ++i){
        //cout << "i:" << i << " proc_to_task_xadj_[i]:" << proc_to_task_xadj_[i] << endl;
        proc_to_task_xadj[i] = proc_to_task_xadj_[i];
    }
    for (procId_t i = 0; i < numTasks; ++i){
        //cout << "i:" << i << " proc_to_task_adj_[i]:" << proc_to_task_adj_[i] << endl;
        proc_to_task_adj[i] = proc_to_task_adj_[i];
    }
    //cout << "done 3" << endl;
    delete ctm;
    delete envConst_;
    //cout << "done 4" << endl;

}


// KDDKDD TODO:  This interface should go away or move to MiniGhost;
// KDDKDD TODO:  it is only a convenience for the MiniGhost experiments.
template <typename procId_t, typename pcoord_t, typename tcoord_t>
void coordinateTaskMapperInterface_Fortran(
  int *comm_World,
  int procDim,
  int numProcessors,
  pcoord_t **machine_coords_,
  int taskDim,
  procId_t numTasks,
  tcoord_t **task_coords,
  procId_t *task_communication_xadj_,
  procId_t *task_communication_adj_,
  procId_t *proc_to_task_xadj, /*output*/
  procId_t *proc_to_task_adj, /*output*/
  int partArraySize,
  procId_t *partNoArray,
  int *machineDimensions
)
{

#ifdef HAVE_MPI
  MPI_Comm cComm = MPI_Comm_f2c((MPI_Fint) *comm_World);
  RCP<const Teuchos::Comm<int> > tcomm = RCP<const Teuchos::Comm<int> > (new Teuchos::MpiComm<int> (cComm));
#else
  RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
#endif
  coordinateTaskMapperInterface<procId_t, pcoord_t, tcoord_t>(
            tcomm,
	    procDim,
            numProcessors,
            machine_coords_,
            taskDim,
            numTasks,
            task_coords,
            task_communication_xadj_,
            task_communication_adj_,
            (pcoord_t *) NULL,
            proc_to_task_xadj, /*output*/
            proc_to_task_adj, /*output*/
            partArraySize,
            partNoArray,
            machineDimensions);
}

}// namespace Zoltan2

#endif
