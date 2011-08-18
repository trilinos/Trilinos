template<class Scalar , class DeviceType >
struct Dirichlet;

template<class Scalar >
struct Dirichlet<Scalar , KOKKOS_MACRO_DEVICE >
{
	
	typedef KOKKOS_MACRO_DEVICE 									device_type;
	typedef device_type::size_type									size_type;
	typedef typename Kokkos::MultiVectorView<Scalar, device_type> 	scalar_vector;	
	typedef typename Kokkos::MultiVectorView<int, device_type> 		int_vector;

  	scalar_vector A;
	int_vector A_row;
	int_vector A_col;
  	scalar_vector b;
  	double value;

	int nx, ny, nz;
  
  	Dirichlet(	scalar_vector & arg_A, 
				int_vector & arg_A_row,
				int_vector & arg_A_col, 
				scalar_vector & arg_b, 
				int nodes_x,
				int nodes_y,
				int nodes_z,
				double arg_value) : 

				A(arg_A),
				A_row(arg_A_row),
				A_col(arg_A_col), 
				b(arg_b)
	{
  		value = arg_value;
	
		nx = nodes_x;
		ny = nodes_y;
		nz = nodes_z;

	}
  	
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(int irow) const {
	//	For each node with x coordinate == 0,
	//	Apply a dirichlet boundary condition:
	//		temperature = value

	//	to maintain the symmetry of the original 
	//	global stiffness matrix, zero out the columns
	//	that correspond to boundary conditions, and
	//	adjust the load vector accordingly

	//	if irow is _NOT_ a boundary condition
		if( (irow % (nx * ny)) % ny != 0){
 
 			int stop = A_row(irow+1);

		//	iterate over the columns in irow
			for(int i = A_row(irow); i < stop; i++){

			//	find any columns that are
			//	boundary conditions. Clear them,
			//	and adjust the load vector
				if((A_col(i) % (nx * ny)) % ny == 0){
					b(irow) -= value * A( i);
					A( i) = 0;
				}
			}
		}
	

	//	if irow _IS_ a boundary condition
		else{

		//	set the load vector equal to a specified value
			b(irow) = value;
	
			int stop = A_row(irow+1);

		//	zero each value on the row, and leave a one
		//	on the diagonal
			for(int i = A_row(irow); i < stop; i++){

				if(A_col(i) == irow)
					A(i) = 1;
				else 
					A(i) = 0;
				

			}


		}

  	}

};
