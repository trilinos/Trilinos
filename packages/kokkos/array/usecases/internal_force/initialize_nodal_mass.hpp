template< typename Scalar , class DeviceType >
struct initialize_nodal_mass;

template<typename Scalar>
struct initialize_nodal_mass<Scalar, KOKKOS_MACRO_DEVICE>
{
	typedef KOKKOS_MACRO_DEVICE     device_type ;
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;
  typedef typename Kokkos::MDArrayView<int,device_type>    int_array_type ;

  enum { NumNodePerElement = 8 };

  initialize_nodal_mass(
      const int_array_type & arg_node_elem_offsets,
      const int_array_type & arg_node_elem_ids,
      const array_type     & arg_elem_mass,
      const array_type     & arg_nodal_mass)
    : node_elem_offset(arg_node_elem_offsets)
    , node_elem_ids(arg_node_elem_ids)
    , elem_mass(arg_elem_mass)
    , nodal_mass(arg_nodal_mass)
  {}


  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int inode )const {

    const int begin = node_elem_offset(inode);
    const int end = node_elem_offset(inode+1);

    Scalar node_mass = 0;
    for(int i = begin; i != end; ++i) {
      node_mass += elem_mass(node_elem_ids(i)) / NumNodePerElement;
    }

    nodal_mass(inode) = node_mass;

  }

  const int_array_type node_elem_offset;
  const int_array_type node_elem_ids;
	const array_type  elem_mass;
	const array_type  nodal_mass;
};
