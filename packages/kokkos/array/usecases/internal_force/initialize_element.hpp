template< typename Scalar , class DeviceType >
struct initialize_element;

template<typename Scalar>
struct initialize_element<Scalar, KOKKOS_MACRO_DEVICE>
{
	typedef KOKKOS_MACRO_DEVICE     device_type ;
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;
  typedef typename Kokkos::MDArrayView<int,device_type>    int_array_type ;

  enum { NumNodePerElement = 8 };

  initialize_element(
      const int_array_type & arg_elem_node_connectivity,
      const array_type     & arg_model_coords,
      const array_type     & arg_elem_mass,
      const array_type     & arg_stretch,
      const array_type     & arg_rotation,
      const Scalar           arg_density
      )
    : elem_node_connectivity(arg_elem_node_connectivity)
    , model_coords(arg_model_coords)
    , elem_mass(arg_elem_mass)
    , stretch(arg_stretch)
    , rotation(arg_rotation)
    , density(arg_density)
  {}


  KOKKOS_MACRO_DEVICE_FUNCTION
  void get_nodes( int ielem, int * nodes) const
  {
    for( int i=0; i<NumNodePerElement; ++i) {
      nodes[i] = elem_node_connectivity(ielem,i);
    }
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int ielem )const {
    int nodes[NumNodePerElement];
    get_nodes(ielem,nodes);

    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    const int K_XX = 0;
    const int K_YY = 1;
    const int K_ZZ = 2;

    stretch(ielem,K_XX) = 1;
    stretch(ielem,K_YY) = 1;
    stretch(ielem,K_ZZ) = 1;

    rotation(ielem,K_XX,0) = 1;
    rotation(ielem,K_YY,0) = 1;
    rotation(ielem,K_ZZ,0) = 1;

    rotation(ielem,K_XX,1) = 1;
    rotation(ielem,K_YY,1) = 1;
    rotation(ielem,K_ZZ,1) = 1;

    Scalar x[NumNodePerElement];
    Scalar y[NumNodePerElement];
    Scalar z[NumNodePerElement];

    for( int i=0; i<NumNodePerElement; ++i) {
      z[i] = model_coords(nodes[i], Z);
    }

    //   calc Z difference vectors
    Scalar R42=(z[3] - z[1]);
    Scalar R52=(z[4] - z[1]);
    Scalar R54=(z[4] - z[3]);

    Scalar R63=(z[5] - z[2]);
    Scalar R83=(z[7] - z[2]);
    Scalar R86=(z[7] - z[5]);

    Scalar R31=(z[2] - z[0]);
    Scalar R61=(z[5] - z[0]);
    Scalar R74=(z[6] - z[3]);

    Scalar R72=(z[6] - z[1]);
    Scalar R75=(z[6] - z[4]);
    Scalar R81=(z[7] - z[0]);


    Scalar t1=(R63 + R54);
    Scalar t2=(R61 + R74);
    Scalar t3=(R72 + R81);

    Scalar t4 =(R86 + R42);
    Scalar t5 =(R83 + R52);
    Scalar t6 =(R75 + R31);

    for( int i=0; i<NumNodePerElement; ++i) {
      y[i] = model_coords(nodes[i], Y);
    }

    Scalar grad_x[NumNodePerElement];

    grad_x[0] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5) + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54);
    grad_x[1] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1) - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
    grad_x[2] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2) - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
    grad_x[3] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3) + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
    grad_x[4] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2) - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
    grad_x[5] = (y[6] *  t5) - (y[4] *  t3) - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
    grad_x[6] = (y[7] *  t1) - (y[5] *  t5) - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
    grad_x[7] = (y[4] *  t2) - (y[6] *  t1) + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);

    for( int i=0; i<NumNodePerElement; ++i) {
      x[i] = model_coords(nodes[i], X);
    }

    elem_mass(ielem) = (1.0/12.0) * density *
                      ( x[0] * grad_x[0] +
                        x[1] * grad_x[1] +
                        x[2] * grad_x[2] +
                        x[3] * grad_x[3] +
                        x[4] * grad_x[4] +
                        x[5] * grad_x[5] +
                        x[6] * grad_x[6] +
                        x[7] * grad_x[7] );

  }

  const int_array_type elem_node_connectivity;
	const array_type  model_coords;
	const array_type  elem_mass;
	const array_type  stretch;
	const array_type  rotation;
  const Scalar      density;
};
