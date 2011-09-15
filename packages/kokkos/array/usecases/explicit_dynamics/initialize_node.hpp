template< typename Scalar , class DeviceType >
struct initialize_node;

template<typename Scalar>
struct initialize_node<Scalar, KOKKOS_MACRO_DEVICE>
{
	typedef KOKKOS_MACRO_DEVICE     device_type ;
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;
  typedef typename Kokkos::MDArrayView<int,device_type>    int_array_type ;

  enum { NumNodePerElement = 8 };

  typedef Region<Scalar,device_type> MyRegion;

  initialize_node( const MyRegion & arg_region )
    : region(arg_region)
  {}


  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int inode )const {

    const int begin = region.node_elem_offset(inode);
    const int end = region.node_elem_offset(inode+1);

    Scalar node_mass = 0;
    for(int i = begin; i != end; ++i) {
      node_mass += region.elem_mass(region.node_elem_ids(i,0)) / NumNodePerElement;
    }

    region.nodal_mass(inode) = node_mass;

  }

  const MyRegion region;
};
