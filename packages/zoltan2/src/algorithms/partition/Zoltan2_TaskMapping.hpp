#include <fstream>
#include <ctime>
#include <vector>
#include "Zoltan2_AlgPQJagged.hpp"
#include "Teuchos_ArrayViewDecl.hpp"

namespace Zoltan2{

#define MINOF(a,b) (((a)<(b))?(a):(b))
/*
template <typename T>
T *allocMemory(size_t size){
    if (size > 0){
        T * a = new T[size];
        if (a == NULL) {
            throw  "cannot allocate memory";
        }
        return a;
    }
    else {
        return NULL;
    }
}

template <typename T>
void freeArray(T *&array){
    if(array != NULL){
        delete [] array;
        array = NULL;
    }
}
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
            arr[i] = v;
        }
    }
}



template <typename pcoord_t,  typename tcoord_t, typename procId_t>
class CoordinateModelInput{
public:
//private:
    int proc_coord_dim;
    pcoord_t **proc_coords;
    int task_coord_dim;
    tcoord_t **task_coords;
    procId_t no_procs;
    procId_t no_tasks;
//public:
    CoordinateModelInput():proc_coord_dim(0), proc_coords(0),
                           task_coord_dim(0), task_coords(0),
                           no_procs(0), no_tasks(0){}
    ~CoordinateModelInput(){}
    CoordinateModelInput(int pcoord_dim_, pcoord_t **pcoords_,
                         int tcoord_dim_, tcoord_t **tcoords_,
                         procId_t no_procs_, procId_t no_tasks_):
                             proc_coord_dim(pcoord_dim_), proc_coords(pcoords_),
                             task_coord_dim(tcoord_dim_), task_coords(tcoords_),
                             no_procs(no_procs_), no_tasks(no_tasks_){
/*
        for(int j = 0; j <  no_tasks; ++j){
            for(int z = 0; z <  task_coord_dim; ++z){
                cout << "z:" << z << " j:" <<  j << " " << this->task_coords[z][j] << endl;
            }

        }
*/
    }


    void getClosestSubset(procId_t *permutation, procId_t nprocs, procId_t ntasks){

        fillContinousArray<procId_t>(permutation, nprocs, NULL);
        int _u_umpa_seed = 847449649;
        srand (time(NULL));
        int a = rand() % 1000 + 1;
        _u_umpa_seed -= a;
        update_visit_order(permutation, nprocs,_u_umpa_seed, 1);

    }

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


    void getMapping(
                    RCP<const Environment> env,
                    procId_t *&processor_to_part,
                    procId_t *&task_to_part,
                    procId_t *&proc_permutation,
                    procId_t *&proc_begins,
                    procId_t *&assigned_task_permutation,
                    procId_t *&assigned_task_begins,
                    procId_t &nprocs,
                    procId_t &ntasks
                    ){
/*
        for(int j = 0; j <  no_tasks; ++j){
            for(int z = 0; z <  task_coord_dim; ++z){
                cout << "z:" << z << " j:" <<  j << " " << this->task_coords[z][j] << endl;
            }

        }
*/
        nprocs = this->no_procs;
        ntasks = this->no_tasks;
        procId_t num_parts = MINOF(this->no_procs, this->no_tasks);

        procId_t minCoordDim = MINOF(this->task_coord_dim, this->proc_coord_dim);
        proc_permutation = NULL;
        proc_begins = allocMemory<procId_t> (num_parts);
        assigned_task_permutation = NULL;
        assigned_task_begins = allocMemory<procId_t> (num_parts);
        /*
        procId_t *proc_permutation = NULL;
        procId_t *proc_begins = allocMemory<procId_t> (num_parts);

        procId_t *assigned_task_permutation = NULL;
        procId_t *assigned_task_begins = allocMemory<procId_t> (num_parts);
        */
        procId_t used_num_procs = this->no_procs;
        if(this->no_procs > this->no_tasks){
            proc_permutation = allocMemory<procId_t>(this->no_procs);
            assigned_task_permutation = allocMemory<procId_t>(this->no_tasks);
            this->getClosestSubset(proc_permutation, this->no_procs, this->no_tasks);
            fillContinousArray<procId_t>(assigned_task_permutation,this->no_tasks, NULL);
            used_num_procs = this->no_tasks;
        }
        else {
            proc_permutation = allocMemory<procId_t>(this->no_procs);
            assigned_task_permutation = allocMemory<procId_t>(this->no_tasks);

            fillContinousArray<procId_t>(proc_permutation,this->no_procs, NULL);
            fillContinousArray<procId_t>(assigned_task_permutation,this->no_tasks, NULL);
        }

        for (int i = 0; i < used_num_procs; ++i){
            cout << "i:" << i << " proc_permutation:" << proc_permutation[i] << endl;
        }
        cout << "seq1" << endl;
        sequentialTaskPartitioning<pcoord_t, procId_t, procId_t>(
                env,
                this->no_procs,
                used_num_procs,
                num_parts,
                minCoordDim,
                this->proc_coords,
                proc_permutation,
                proc_begins,
                "proc_partitioning"
                );

        cout << "seq2" << endl;
        sequentialTaskPartitioning<tcoord_t, procId_t, procId_t>(
                        env,
                        this->no_tasks,
                        this->no_tasks,
                        num_parts,
                        minCoordDim,
                        this->task_coords,
                        assigned_task_permutation,
                        assigned_task_begins,
                        "task_partitioning"
                        );

        procId_t invalid = -1;
        processor_to_part = allocMemory<procId_t> (this->no_procs);
        fillContinousArray<procId_t> (processor_to_part, this->no_procs, &invalid);

        for(procId_t i = 0; i < num_parts; ++i){

            procId_t proc_index_begin = 0;
            if (i > 0) proc_index_begin = proc_begins[i - 1];
            procId_t proc_index_end = proc_begins[i];

            for (procId_t j = proc_index_begin; j < proc_index_end; ++j){

                procId_t assigned_proc = proc_permutation[j];
                processor_to_part[assigned_proc] = i;

                //cout << " part:" << i << " proc:" << assigned_proc << endl;
            }
        }

        task_to_part = allocMemory<procId_t> (this->no_tasks);
        fillContinousArray<procId_t> (task_to_part, this->no_tasks, &invalid);
        for(procId_t i = 0; i < num_parts; ++i){

            procId_t task_begin_index = 0;
            if (i > 0) task_begin_index = assigned_task_begins[i - 1];
            procId_t task_index_end = assigned_task_begins[i];

            for (procId_t j = task_begin_index; j < task_index_end; ++j){
                procId_t assigned_task = assigned_task_permutation[j];
                task_to_part[assigned_task] = i;
                //cout << " part:" << i << " assigned_task:" << assigned_task << endl;
            }
        }
        /*
        for(int j = 0; j <  no_tasks; ++j){
            for(int z = 0; z <  task_coord_dim; ++z){
                cout << "z:" << z << " j:" <<  j << " " << this->task_coords[z][j] << endl;
            }

        }
        */
    }

};

template <typename communicationModel, typename procId_t>
class TaskMapper{
protected:

    RCP<const Environment> env;
    procId_t *processor_to_part;
    procId_t *task_to_part;
    procId_t *proc_permutation;
    procId_t *proc_begins;
    procId_t *assigned_task_permutation;
    procId_t *assigned_task_begins;
    communicationModel *proc_task_comm;
    procId_t nprocs;
    procId_t ntasks;

public:
    TaskMapper():
                env(new Environment()),
                processor_to_part(0),
                task_to_part(0),
                proc_permutation(0),
                proc_begins(0),
                assigned_task_permutation(0),
                assigned_task_begins(0),
                proc_task_comm(0),
                nprocs(0),
                ntasks(0){}

    ~TaskMapper(){
        freeArray<procId_t> (processor_to_part);
        freeArray<procId_t> (task_to_part);
        freeArray<procId_t> (proc_permutation);
        freeArray<procId_t> (proc_begins);
        freeArray<procId_t> (assigned_task_permutation);
        freeArray<procId_t> (assigned_task_begins);
    }

    TaskMapper(communicationModel *proc_task_comm_):
        env(new Environment()),
        processor_to_part(0),
        task_to_part(0),
        proc_permutation(0),
        proc_begins(0),
        assigned_task_permutation(0),
        assigned_task_begins(0),
        proc_task_comm(proc_task_comm_),
        nprocs(0),
        ntasks(0){
        this->doMapping();

    }

    void doMapping(){
        if(this->proc_task_comm){
        this->proc_task_comm->getMapping(
                this->env,
                this->processor_to_part,
                this->task_to_part,
                this->proc_permutation,
                this->proc_begins,
                this->assigned_task_permutation,
                this->assigned_task_begins,
                this->nprocs,
                this->ntasks);
        }
        else {
            std::cerr << "communicationModel is not specified in the Mapper" << endl;
            exit(1);
        }
    }

    ArrayView<procId_t> getAssignedProcsForTask(procId_t taskId){
        if(taskId < this->ntasks){
            //cout << "t:" << taskId << endl;
            procId_t mapped_part =  this->task_to_part[taskId];
            //cout << "mapped_part:" << mapped_part << endl;

            procId_t partbegin = 0;
            if (mapped_part > 0) partbegin = this->proc_begins[mapped_part - 1];
            procId_t partend = this->proc_begins[mapped_part];
            //cout << "partbegin:" << partbegin << endl;
            //cout << "partend:" << partend << endl;
            //create arrayview here.
            ArrayView <procId_t> assignedProcs(this->proc_permutation + partbegin, partend - partbegin);
            return assignedProcs;
        }
        else {
            ArrayView <procId_t> a(NULL, 0);
            return a;
        }
    }

    ArrayView<procId_t> getAssignedTaksForProc(procId_t procId){
        //cout << "proc id:" << procId << endl;
        procId_t mapped_part =  this->processor_to_part[procId];
        if(mapped_part != -1){
            procId_t partbegin = 0;
            if (mapped_part > 0) partbegin = this->assigned_task_begins[mapped_part - 1];
            procId_t partend = this->assigned_task_begins[mapped_part];
            //create and return arrayview here.
            //cout << "task for proc:" << endl;
            ArrayView <procId_t> assignedParts(this->assigned_task_permutation + partbegin, partend - partbegin);
            return assignedParts;
        } else {
            //cout << "return empty" << endl;
            ArrayView <procId_t> a(NULL,0);
            //cout << "return empty done" << endl;
            return a;
        }
    }

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


    void writeMapping2(){

        for(procId_t i = 0; i < this->nprocs; ++i){
            cout << "Proc:" << i << " assignedParts:" << this->getAssignedTaksForProc(i) << endl;
        }
        //cout << "parts" << endl;
        for(procId_t i = 0; i < this->ntasks; ++i){
            cout << "Part:" << i << " assignedProcs:" << this->getAssignedProcsForTask(i) << endl;
        }

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

        gnuPlotCode << ss;
        if(mindim == 2){
            gnuPlotCode << "plot \"procPlot.plot\" with points pointsize 3\n";
        } else {
            gnuPlotCode << "splot \"procPlot.plot\" with points pointsize 3\n";
        }
        gnuPlotCode << "replot \"partPlot.plot\" with points pointsize 3\n";
        gnuPlotCode << "\nreplot\n pause -1 \n";
        gnuPlotCode.close();

    }

};
}// namespace Zoltan2

