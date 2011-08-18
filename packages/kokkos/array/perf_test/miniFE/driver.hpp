#include <cstddef>
#include <SparseOutput.hpp>


namespace Test {

template<class DeviceType >
void run_kernel(int, int, int, double*);

template<>
void run_kernel<KOKKOS_MACRO_DEVICE>(int x, int y, int z, double* times) 
{
	typedef double Scalar;
	typedef KOKKOS_MACRO_DEVICE     									device_type;

	typedef Kokkos::MDArrayView<Scalar,device_type>::HostView			HostView_scalar;
	typedef Kokkos::MDArrayView<int,device_type>::HostView				HostView_int;

  	typedef Kokkos::MDArrayView<Scalar,device_type> 					scalar_array_d;
  	typedef Kokkos::MDArrayView<int,device_type> 						int_array_d;  	

	typedef Kokkos::MultiVectorView<Scalar , device_type>				scalar_vector_d;
	typedef Kokkos::MultiVectorView<int , device_type>					int_vector_d;

	typedef Kokkos::MultiVectorView<Scalar , Kokkos::DeviceHost>		scalar_vector_h;
	typedef Kokkos::MultiVectorView<int , Kokkos::DeviceHost>			int_vector_h;

	int nelems = x * y * z;
	int nnodes = (x + 1) * (y + 1) * (z + 1);

//	Host Data Structures
	HostView_int 	elem_nodeIDs_h, node_elemIDs_h, elems_per_node_h;
	HostView_scalar 	elem_coords_h;

	int_vector_h 	A_row_h , A_col_h;


//	Device Data Structures
	int_array_d 	elem_nodeIDs_d, node_elemIDs_d, elems_per_node_d;
	scalar_array_d 	elem_coords_d, elem_stiffness, elem_load;

	int_vector_d A_row_d , A_col_d;
	scalar_vector_d A, b, X;
	
	timeval start,stop,result;
	double time = 0.0;
	double total_time = 0.0;


	gettimeofday(&start, NULL);

	init_mesh_and_crsgraph(
			elem_coords_h, 
			elem_nodeIDs_h, 
			node_elemIDs_h, 
			elems_per_node_h, 
			A_row_h, 
			A_col_h, 
			x, y, z);

	gettimeofday(&stop, NULL);
	timersub(&stop, &start, &result);
	time = (result.tv_sec + result.tv_usec/1000000.0);

	total_time += time;	
	times[0] = time;

//	copy host data to device
	elem_coords_d = 	Kokkos::create_mdarray< scalar_array_d > (nelems, 3, 8);
	elem_nodeIDs_d = 	Kokkos::create_mdarray< int_array_d >(nelems, 8);
	node_elemIDs_d = 	Kokkos::create_mdarray< int_array_d >(node_elemIDs_h.dimension(0), 2);
	elems_per_node_d = 	Kokkos::create_mdarray< int_array_d >(elems_per_node_h.dimension(0));

	A_row_d = 			Kokkos::create_labeled_multivector< int_array_d >("A_row_d",A_row_h.length());
	A_col_d = 			Kokkos::create_labeled_multivector< int_array_d >("A_col_d",A_col_h.length());

	elem_stiffness =	Kokkos::create_mdarray< scalar_array_d > (nelems, 8, 8);
	elem_load =			Kokkos::create_mdarray< scalar_array_d > (nelems, 8);

	A = 				Kokkos::create_labeled_multivector< scalar_vector_d > ("A",A_col_h.length());	
	b =					Kokkos::create_labeled_multivector< scalar_vector_d > ("b",nelems, 8);	
	X =					Kokkos::create_labeled_multivector< scalar_vector_d > ("X",nnodes);

	gettimeofday(&start, NULL);

	Kokkos::deep_copy(elem_coords_d, 	elem_coords_h);
	Kokkos::deep_copy(elem_nodeIDs_d, 	elem_nodeIDs_h);
	Kokkos::deep_copy(node_elemIDs_d, 	node_elemIDs_h);
	Kokkos::deep_copy(elems_per_node_d, elems_per_node_h);
	Kokkos::deep_copy(A_row_d, 			A_row_h);
	Kokkos::deep_copy(A_col_d, 			A_col_h);

	gettimeofday(&stop, NULL);
	timersub(&stop, &start, &result);
	time = (result.tv_sec + result.tv_usec/1000000.0);

	total_time += time;	


	Kokkos::parallel_for(nelems, assembleFE<Scalar, device_type>(	
		elem_stiffness, 
		elem_load,
		elem_coords_d,
		x, 
		y, 
		z), time);

	total_time += time;	
	times[1] = time;


	Kokkos::parallel_for(nnodes, CRSMatrixGatherFill<Scalar, device_type>(	
		A,
		b,
		A_row_d,
		A_col_d,
		node_elemIDs_d,
		elem_nodeIDs_d,
		elems_per_node_d,
		elem_stiffness,
		elem_load), time);

	total_time += time;	
	times[2] = time;

	Kokkos::parallel_for(nnodes , Dirichlet<Scalar , device_type>(A, A_row_d ,A_col_d, b,x+1,y+1,z+1,1.0) , time);
	total_time += time;
	times[3] = time;

//	printSparse< scalar_vector_d , int_vector_d>("A.txt",A,A_row_d,A_col_d);

	time = CG_Solve<Scalar, device_type>::run(A , A_row_d, A_col_d , b , X ,times );
	total_time += time;
	
	times[6] = total_time;

//	printGLUT<Scalar , scalar_vector_d , HostView_scalar , HostView_int>
//			("X.txt", X , elem_coords_h , elem_nodeIDs_h,x,y,z);
}

} //Test
