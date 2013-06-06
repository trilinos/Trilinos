#include <fstream>
#include <ctime>
#include <vector>
#include "Zoltan2_AlgPQJagged.hpp"
#include "Teuchos_ArrayViewDecl.hpp"

namespace Zoltan2{

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
        distance = pow(distance, 1.0f / this->dimension);
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
                numClusters (pow(2, dim_) + 1),
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
            int mod = pow(2,j + 1);
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


/*! \brief CoordinateModelInput Class that performs mapping between the coordinate partitioning result and mpi ranks
 * base on the coordinate results and mpi physical coordinates.
 */
template <typename pcoord_t,  typename tcoord_t, typename procId_t>
class CoordinateModelInput{
public:
//private:
    int proc_coord_dim; //dimension of the processors
    pcoord_t **proc_coords; //the processor coordinates. allocated outside of the class.
    int task_coord_dim; //dimension of the tasks coordinates.
    tcoord_t **task_coords; //the task coordinates allocated outside of the class.
    procId_t no_procs; //the number of processors
    procId_t no_tasks;  //the number of taks.
//public:
    CoordinateModelInput():proc_coord_dim(0), proc_coords(0),
                           task_coord_dim(0), task_coords(0),
                           no_procs(0), no_tasks(0){}
    ~CoordinateModelInput(){}
    procId_t getNProcs(){
        return this->no_procs;
    }
    procId_t getNTasks(){
        return this->no_tasks;
    }
    /*! \brief Class Constructor:
     *  \param pcoord_dim_ the dimension of the processors
     *  \param pcoords_   the processor coordinates. allocated outside of the class.
     *  \param tcoord_dim_   dimension of the tasks coordinates.
     *  \param tcoords_   the task coordinates allocated outside of the class.
     *  \param no_procs_   the number of processors
     *  \param no_tasks_   the number of taks.
     */
    CoordinateModelInput(int pcoord_dim_, pcoord_t **pcoords_,
                         int tcoord_dim_, tcoord_t **tcoords_,
                         procId_t no_procs_, procId_t no_tasks_):
                             proc_coord_dim(pcoord_dim_), proc_coords(pcoords_),
                             task_coord_dim(tcoord_dim_), task_coords(tcoords_),
                             no_procs(no_procs_), no_tasks(no_tasks_){
    }


    /*! \brief Function is called whenever nprocs > no_task.
     * Function returns only the subset of processors that are closest to each other.
     *  \param proc_permutation holds the indices of the processors that are chosen.
     *  \param nprocs the number of processors.
     *  \param ntasks the number of taks.
     */
    void getClosestSubset(procId_t *proc_permutation, procId_t nprocs, procId_t ntasks){
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
    void getMapping(
                    RCP<const Environment> env,
                    procId_t *&proc_to_task_xadj, //  = allocMemory<procId_t> (this->no_procs); //holds the pointer to the task array
                    procId_t *&proc_to_task_adj, // = allocMemory<procId_t>(this->no_tasks); //holds the indices of tasks wrt to proc_to_task_xadj array.
                    procId_t *&task_to_proc //allocMemory<procId_t>(this->no_tasks); //holds the processors mapped to tasks.
                    ){


        proc_to_task_xadj = allocMemory<procId_t> (this->no_procs); //holds the pointer to the task array
        proc_to_task_adj = allocMemory<procId_t>(this->no_tasks); //holds the indices of tasks wrt to proc_to_task_xadj array.
        task_to_proc = allocMemory<procId_t>(this->no_tasks); //holds the processors mapped to tasks.);

        procId_t invalid = 0;
        fillContinousArray<procId_t> (proc_to_task_xadj, this->no_procs, &invalid);

        //obtain the number of parts that should be divided.
        procId_t num_parts = MINOF(this->no_procs, this->no_tasks);
        //obtain the min coordinate dim.
        procId_t minCoordDim = MINOF(this->task_coord_dim, this->proc_coord_dim);

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

/*
        for(int i = 0; i < no_procs; ++i){
            cout << "permute:" << i << " proc:"  << proc_adjList[i] << endl;
        }
*/
        //partitioning of processors
        sequentialTaskPartitioning<pcoord_t, procId_t, procId_t>(
                env,
                this->no_procs,
                used_num_procs,
                num_parts,
                minCoordDim,
                this->proc_coords,
                proc_adjList,
                proc_xadj,
                "proc_partitioning"
                );

        procId_t *task_xadj = allocMemory<procId_t> (num_parts);
        procId_t *task_adjList = allocMemory<procId_t>(this->no_tasks);
        //fill task_adjList st: task_adjList[i] <- i.
        fillContinousArray<procId_t>(task_adjList,this->no_tasks, NULL);

        //partitioning of tasks
        sequentialTaskPartitioning<tcoord_t, procId_t, procId_t>(
                        env,
                        this->no_tasks,
                        this->no_tasks,
                        num_parts,
                        minCoordDim,
                        this->task_coords,
                        task_adjList,
                        task_xadj,
                        "task_partitioning"
                        );


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
                cerr << "PART:" << i << " is assigned to " << proc_index_end - proc_index_begin << " processors." << endl;
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

/*
        for(procId_t i = 0; i < this->no_procs; ++i){
            cout << " i: " << i << " "<< proc_to_task_xadj[i] << endl;
        }
*/
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

template <typename communicationModel, typename procId_t>
class TaskMapper{
protected:

    RCP<const Environment> env;
    procId_t *proc_to_task_xadj; //  = allocMemory<procId_t> (this->no_procs); //holds the pointer to the task array
    procId_t *proc_to_task_adj; // = allocMemory<procId_t>(this->no_tasks); //holds the indices of tasks wrt to proc_to_task_xadj array.
    procId_t *task_to_proc; //allocMemory<procId_t>(this->no_procs); //holds the processors mapped to tasks.
    communicationModel *proc_task_comm;
    procId_t nprocs;
    procId_t ntasks;

public:
/*
    TaskMapper():
                env(),
                proc_to_task_xadj(0),
                proc_to_task_adj(0),
                task_to_proc(0),
                proc_task_comm(0),
                nprocs(0),
                ntasks(0){}
*/
    ~TaskMapper(){
        freeArray<procId_t> (proc_to_task_xadj);
        freeArray<procId_t> (proc_to_task_adj);
        freeArray<procId_t> (task_to_proc);
    }
    /*! \brief Constructor
     *  \param env_ Environment object.
     *  \param proc_task_comm_ is the templated parameter for which the mapping will be obtained with getMapping() function.
     */
    TaskMapper(RCP<const Environment> env_, communicationModel *proc_task_comm_):
        env(env_),
        proc_to_task_xadj(0),
        proc_to_task_adj(0),
        task_to_proc(0),
        proc_task_comm(proc_task_comm_),
        nprocs(proc_task_comm->getNProcs()),
        ntasks(proc_task_comm->getNTasks()){
        //calls doMapping function
        this->doMapping();

    }

    /*! \brief doMapping function, calls getMapping function of communicationModel object.
     */
    void doMapping(){

        if(this->proc_task_comm){
            this->proc_task_comm->getMapping(
                    this->env,
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

    /*! \brief getAssignedProcForTask function, returns the assigned processor id for the given task
     *  \param taskId taskId being queried.
     */
    inline procId_t getAssignedProcForTask(procId_t taskId){
        return this->task_to_proc[taskId];
    }
    /*! \brief getAssignedProcForTask function, returns the assigned tasks in ArrayView format for the given processor.
     *  \param procId procId being queried.
     */
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
            ArrayView <procId_t> assignedParts(this->proc_to_task_adj + task_begin, taskend - task_begin);
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
    void writeMapping2(){
/*
        for(procId_t i = 0; i < this->nprocs; ++i){
            cout << "Proc:" << i << " assignedParts:" << this->getAssignedTaksForProc(i) << endl;
        }
        //cout << "parts" << endl;
        for(procId_t i = 0; i < this->ntasks; ++i){
            cout << "Part:" << i << " assignedProcs:" << this->getAssignedProcForTask(i) << endl;
        }
*/
        std::ofstream gnuPlotCode ("gnuPlot2.plot", std::ofstream::out);

        int mindim = MINOF(proc_task_comm->proc_coord_dim, proc_task_comm->task_coord_dim);
        string ss = "";
        string procs = "", parts = "";
        for(procId_t i = 0; i < this->nprocs; ++i){

            //inpFile << std::endl;
            ArrayView<procId_t> a = this->getAssignedTaksForProc(i);
            if (a.size() == 0){
                continue;
            }
            /*
            std::string procFile = toString<int>(i) + "_mapping.txt";
            if (i == 0){
                gnuPlotCode << "plot \"" << procFile << "\"\n";
            }
            else {
                gnuPlotCode << "replot \"" << procFile << "\"\n";
            }
            */


            //std::ofstream inpFile (procFile.c_str(), std::ofstream::out);

            string gnuPlotArrow = "set arrow from ";
            for(int j = 0; j <  mindim; ++j){
                if (j == mindim - 1){
                    //inpFile << proc_task_comm->proc_coords[j][i];
                    gnuPlotArrow += toString<float>(proc_task_comm->proc_coords[j][i]);
                    procs += toString<float>(proc_task_comm->proc_coords[j][i]);

                }
                else {
                    //inpFile << proc_task_comm->proc_coords[j][i] << " ";
                    gnuPlotArrow += toString<float>(proc_task_comm->proc_coords[j][i]) +",";
                    procs += toString<float>(proc_task_comm->proc_coords[j][i])+ " ";
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
                        gnuPlotArrow2 += toString<float>(proc_task_comm->task_coords[z][j]);
                        parts += toString<float>(proc_task_comm->task_coords[z][j]);
                    }
                    else{
                        //inpFile << proc_task_comm->task_coords[z][j] << " ";
                        gnuPlotArrow2 += toString<float>(proc_task_comm->task_coords[z][j]) +",";
                        parts += toString<float>(proc_task_comm->task_coords[z][j]) + " ";
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
                extraProcFile << this->proc_task_comm->proc_coords[i][j] <<  " ";
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
        gnuPlotCode << "replot \"allProc.plot\" with points pointsize 0.5\n";
        gnuPlotCode << "\nreplot\n pause -1 \n";
        gnuPlotCode.close();

    }

};
}// namespace Zoltan2

