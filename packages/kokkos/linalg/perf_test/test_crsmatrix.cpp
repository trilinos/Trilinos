#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>

#ifdef _OPENMP
#include <KokkosArray_OpenMP.hpp>
#else
#include <KokkosArray_Host.hpp>
#endif
#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_MultiVector.hpp>
#include <KokkosArray_CRSMatrix.hpp>
#ifndef DEVICE
#define DEVICE 1
#endif
#if DEVICE==1
#ifdef _OPENMP
typedef KokkosArray::OpenMP device_type;
#else
typedef KokkosArray::Host device_type;
#endif
#define KokkosArrayHost(a) a
#define KokkosArrayCUDA(a)
#else
typedef KokkosArray::Cuda device_type;
#define KokkosArrayHost(a)
#define KokkosArrayCUDA(a) a
#endif

template< typename ScalarType , typename OrdinalType>
int SparseMatrix_generate(OrdinalType nrows, OrdinalType ncols, OrdinalType &nnz, OrdinalType varianz_nel_row, OrdinalType width_row, ScalarType* &values, OrdinalType* &rowPtr, OrdinalType* &colInd)
{
  rowPtr = new OrdinalType[nrows+1];

  OrdinalType elements_per_row = nnz/nrows;
  srand(13721);
  rowPtr[0] = 0;
  for(int row=0;row<nrows;row++)
  {
    int varianz = (1.0*rand()/INT_MAX-0.5)*varianz_nel_row;
    rowPtr[row+1] = rowPtr[row] + elements_per_row+varianz;
  }
  nnz = rowPtr[nrows];
  values = new ScalarType[nnz];
  colInd = new OrdinalType[nnz];
  for(int row=0;row<nrows;row++)
  {
	 for(int k=rowPtr[row];k<rowPtr[row+1];k++)
	 {
		int pos = (1.0*rand()/INT_MAX-0.5)*width_row+row;
		if(pos<0) pos+=ncols;
		if(pos>=ncols) pos-=ncols;
		colInd[k]= pos;
		values[k] = 100.0*rand()/INT_MAX-50.0;
	 }
  }
  return nnz;
}

template<typename Scalar>
int test_crs_matrix_test(int numRows, int numCols, int nnz, int numVecs, int test) {
	typedef KokkosArray::CrsMatrix<Scalar,int,device_type> matrix_type ;
	typedef typename KokkosArray::MultiVectorDynamic<Scalar,device_type>::type mv_type;
	typedef typename KokkosArray::MultiVectorDynamic<Scalar,device_type>::random_read_type mv_random_read_type;
	typedef typename mv_type::HostMirror h_mv_type;

	Scalar* val = NULL;
	int* row = NULL;
	int* col = NULL;

	srand(17312837);
	nnz = SparseMatrix_generate<Scalar,int>(numRows,numCols,nnz,nnz/numRows*0.2,numRows*0.01,val,row,col);

	matrix_type A("CRS::A",numRows,numCols,nnz,val,row,col,false);
	mv_type x("X",numCols,numVecs);
	mv_random_read_type t_x(x);
	mv_type y("Y",numRows,numVecs);
	h_mv_type h_x = KokkosArray::create_mirror_view(x);
	h_mv_type h_y = KokkosArray::create_mirror_view(y);
	h_mv_type h_y_compare = KokkosArray::create_mirror(y);
    typename matrix_type::CrsArrayType::HostMirror h_graph = KokkosArray::create_mirror(A.graph);
    typename matrix_type::values_type::HostMirror h_values = KokkosArray::create_mirror_view(A.values);

    //KokkosArray::deep_copy(h_graph.row_map,A.graph.row_map);
    for(int k=0;k<numVecs;k++){
	  //h_a(k) = (Scalar) (1.0*(rand()%40)-20.);
	  for(int i=0; i<numCols;i++) {
		  h_x(i,k) = (Scalar) (1.0*(rand()%40)-20.);
		  h_y(i,k) = (Scalar) (1.0*(rand()%40)-20.);
	  }
    }
	for(int i=0;i<numRows;i++) {
		int start = h_graph.row_map(i);
		int end = h_graph.row_map(i+1);
		for(int j=start;j<end;j++) {
		   h_graph.entries(j) = (i+j-start)%numCols;
		   h_values(j) = h_graph.entries(j) + i;
		}
		for(int k = 0; k<numVecs; k++)
		  h_y_compare(i,k) = 0;
		for(int j=start;j<end;j++) {
		   Scalar val = h_graph.entries(j) + i;
		   int idx = h_graph.entries(j);
		   for(int k = 0; k<numVecs; k++)
			   h_y_compare(i,k)+=val*h_x(idx,k);
		}
	}

	KokkosArray::deep_copy(x,h_x);
	KokkosArray::deep_copy(y,h_y);
	KokkosArray::deep_copy(A.graph.entries,h_graph.entries);
	KokkosArray::deep_copy(A.values,h_values);
	/*for(int i=0;i<numRows;i++)
		for(int k = 0; k<numVecs; k++) {
          //error[k]+=(h_y_compare(i,k)-h_y(i,k))*(h_y_compare(i,k)-h_y(i,k));
          printf("%i %i %lf %lf %lf\n",i,k,h_y_compare(i,k),h_y(i,k),h_x(i,k));
		}*/

	KokkosArray::MV_Multiply(0.0,y,1.0,A,x);
	device_type::fence();
	KokkosArray::deep_copy(h_y,y);
	Scalar error[numVecs];
	Scalar sum[numVecs];
	for(int k = 0; k<numVecs; k++) {
		error[k] = 0;
		sum[k] = 0;
	}
	for(int i=0;i<numRows;i++)
		for(int k = 0; k<numVecs; k++) {
          error[k]+=(h_y_compare(i,k)-h_y(i,k))*(h_y_compare(i,k)-h_y(i,k));
          sum[k] += h_y_compare(i,k)*h_y_compare(i,k);
         // printf("%i %i %lf %lf %lf\n",i,k,h_y_compare(i,k),h_y(i,k),h_x(i,k));
		}

	//for(int i=0;i<A.nnz;i++) printf("%i %lf\n",h_graph.entries(i),h_values(i));
    int num_errors = 0;
    double total_error = 0;
    double total_sum = 0;
	for(int k = 0; k<numVecs; k++) {
		num_errors += (error[k]/(sum[k]==0?1:sum[k]))>1e-5?1:0;
		total_error += error[k];
		total_sum += sum[k];
	}

    int loop = 100;
	timespec starttime,endtime;
    clock_gettime(CLOCK_REALTIME,&starttime);
	for(int i=0;i<loop;i++)
		KokkosArray::MV_Multiply(0.0,y,1.0,A,t_x);
	device_type::fence();
	clock_gettime(CLOCK_REALTIME,&endtime);
	double time = endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

	double matrix_size = 1.0*((nnz*(sizeof(Scalar)+sizeof(int)) + numRows*sizeof(int)))/1024/1024;
	double vector_size = 2.0*numRows*numVecs*sizeof(Scalar)/1024/1024;
	double vector_readwrite = 2.0*nnz*numVecs*sizeof(Scalar)/1024/1024;

	double problem_size = matrix_size+vector_size;
    printf("%6.2lf MB %6.2lf GB/s %6.2lf s %i\n",problem_size,(matrix_size+vector_readwrite)/time*loop/1024, time/loop*1000, num_errors);
	return (int)total_error;
}


int test_crs_matrix_type(int numrows, int numcols, int nnz, int numVecs, int type, int test) {
  return test_crs_matrix_test<double>(numrows,numcols,nnz,numVecs,test);
}

int main(int argc, char **argv)
{
 long long int size = 110503; // a prime number
 int numVecs = 4;
 int threads=1;
 int device = 0;
 int numa=1;
 int test=-1;
 int type=-1;

 for(int i=0;i<argc;i++)
 {
  if((strcmp(argv[i],"-d")==0)) {device=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-s")==0)) {size=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-v")==0)) {numVecs=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--threads")==0)) {threads=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--numa")==0)) {numa=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--test")==0)) {test=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--type")==0)) {type=atoi(argv[++i]); continue;}
 }


 KokkosArrayCUDA(
   KokkosArray::Cuda::SelectDevice select_device(device);
   KokkosArray::Cuda::initialize( select_device );
 )

#ifdef _OPENMP
   omp_set_num_threads(numa*threads);
   KokkosArray::OpenMP::initialize( numa);
#pragma message "Compile OpenMP"
#else
   KokkosArray::Host::initialize( numa , threads );
#pragma message "Compile PThreads"
#endif

 int numVecsList[10] = {1, 2, 3, 4, 5, 8, 11, 15, 16, 17};
 int maxNumVecs = numVecs==-1?17:numVecs;
 int numVecIdx = 0;
 if(numVecs == -1) numVecs = numVecsList[numVecIdx++];

 int total_errors = 0;
 while(numVecs<=maxNumVecs) {
   total_errors += test_crs_matrix_type(size,size,size*10,numVecs,type,test);
   if(numVecs<maxNumVecs) numVecs = numVecsList[numVecIdx++];
   else numVecs++;
 }

 if(total_errors == 0)
   printf("Kokkos::MultiVector Test: Passed\n");
 else
   printf("Kokkos::MultiVector Test: Failed\n");


 KokkosArrayCUDA(
#ifdef _OPENMP
 KokkosArray::OpenMP::finalize();
#else
 KokkosArray::Host::finalize();
#endif
 )
 device_type::finalize();
}
