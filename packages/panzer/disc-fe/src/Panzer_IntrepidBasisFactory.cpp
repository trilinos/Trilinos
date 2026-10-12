#include "Panzer_IntrepidBasisFactory.hpp"

#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"
#include "Intrepid2_HVOL_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"
#include "Intrepid2_HVOL_TET_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"

#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"
#include "Intrepid2_HCURL_TRI_In_FEM.hpp"

#include "Intrepid2_HCURL_TET_I1_FEM.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"

#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"

#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"
#include "Intrepid2_HDIV_TRI_In_FEM.hpp"

#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"

#include "Intrepid2_HDIV_TET_I1_FEM.hpp"
#include "Intrepid2_HDIV_TET_In_FEM.hpp"

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"


namespace panzer {


  /** \brief Creates an Intrepid2::Basis object given the basis, order and cell topology.

      \param[in] basis_type The name of the basis.
      \param[in] basis_order The order of the polynomial used to construct the basis.
      \param[in] cell_topology Cell topology for the basis.  Taken from shards::CellTopology::getName()
                               after trimming the extended basis suffix.

      To be backwards compatible, this method takes deprecated
      descriptions and transform it into a valid type and order.  For
      example "Q1" is transformed to basis_type="HGrad",basis_order=1.

      \returns A newly allocated panzer::Basis object.
  */
  template <typename DeviceType, typename OutputValueType, typename PointValueType>
  Teuchos::RCP<Intrepid2::Basis<DeviceType,OutputValueType,PointValueType> >
  createIntrepid2Basis(const std::string basis_type, int basis_order,
                       const shards::CellTopology & cell_topology)
  {
    // Shards supports extended topologies so the names have a "size"
    // associated with the number of nodes.  We prune the size to
    // avoid combinatorial explosion of checks.
    std::string cell_topology_type = cell_topology.getName();
    std::size_t end_position = 0;
    end_position = cell_topology_type.find("_");
    std::string cell_type = cell_topology_type.substr(0,end_position);

    Teuchos::RCP<Intrepid2::Basis<DeviceType,OutputValueType,PointValueType> > basis;

    // high order point distribution type;
    // for now equispaced only; to get a permutation map with different orientation,
    // orientation coeff matrix should have the same point distribution.
    const Intrepid2::EPointType point_type = Intrepid2::POINTTYPE_EQUISPACED;

    if ( (basis_type == "Const") && (basis_order == 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_C0_FEM<DeviceType,OutputValueType,PointValueType>(cell_topology) );

    else if ( (basis_type == "HVol") && (basis_order == 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_C0_FEM<DeviceType,OutputValueType,PointValueType>(cell_topology) );

    else if ( (basis_type == "HVol") && (cell_type == "Quadrilateral") && (basis_order > 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_QUAD_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HVol") && (cell_type == "Triangle") && (basis_order > 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_TRI_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HVol") && (cell_type == "Hexahedron") )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_HEX_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HVol") && (cell_type == "Tetrahedron") )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_TET_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HGrad") && (cell_type == "Hexahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Hexahedron") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_C2_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Hexahedron") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HCurl") && (cell_type == "Hexahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_HEX_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Hexahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_HEX_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HDiv") && (cell_type == "Hexahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_HEX_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HDiv") && (cell_type == "Hexahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_HEX_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HGrad") && (cell_type == "Tetrahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TET_C1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Tetrahedron") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TET_C2_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Tetrahedron") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TET_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HCurl") && (cell_type == "Tetrahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TET_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Tetrahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TET_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HDiv") && (cell_type == "Tetrahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TET_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HDiv") && (cell_type == "Tetrahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TET_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HGrad") && (cell_type == "Quadrilateral") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Quadrilateral") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Quadrilateral") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HCurl") && (cell_type == "Quadrilateral") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_QUAD_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Quadrilateral") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_QUAD_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HDiv") && (cell_type == "Quadrilateral") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_QUAD_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HDiv") && (cell_type == "Quadrilateral") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_QUAD_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HGrad") && (cell_type == "Triangle") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TRI_C1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Triangle") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TRI_C2_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Triangle") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TRI_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HCurl") && (cell_type == "Triangle") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TRI_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Triangle") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TRI_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HDiv") && (cell_type == "Triangle") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TRI_I1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HDiv") && (cell_type == "Triangle") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TRI_In_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HGrad") && (cell_type == "Line") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_LINE_C1_FEM<DeviceType,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Line") && (basis_order >= 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<DeviceType,OutputValueType,PointValueType>(basis_order, point_type) );

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
                               "Failed to create the requestedbasis with basis_type=\"" << basis_type <<
                               "\", basis_order=\"" << basis_order << "\", and cell_type=\"" << cell_type << "\"!\n");

    // we compare that the base topologies are the same
    // we do so using the NAME. This avoids the ugly task of getting the
    // cell topology data and constructing a new cell topology object since you cant
    // just get the baseCellTopology directly from a shards cell topology
    TEUCHOS_TEST_FOR_EXCEPTION(cell_topology.getBaseName()!=basis->getBaseCellTopology().getName(),
                               std::runtime_error,
                               "Failed to create basis.  Intrepid2 basis base topology does not match mesh cell base topology!");

    return basis;
  }

  /** \brief Creates an Intrepid2::Basis object given the basis, order and cell topology.

      \param[in] basis_type The name of the basis.
      \param[in] basis_order The order of the polynomial used to construct the basis.
      // TODO BWR Is the cell_topology documentation below correct?
      \param[in] cell_topology Cell topology for the basis.  Taken from shards::CellTopology::getName() after
                               trimming the extended basis suffix.

      To be backwards compatible, this method takes deprecated
      descriptions and transform it into a valid type and order.  For
      example "Q1" is transformed to basis_type="HGrad",basis_order=1.

      \returns A newly allocated panzer::Basis object.
  */
  template <typename DeviceType, typename OutputValueType, typename PointValueType>
  Teuchos::RCP<Intrepid2::Basis<DeviceType,OutputValueType,PointValueType> >
  createIntrepid2Basis(const std::string basis_type, int basis_order,
                      const Teuchos::RCP<const shards::CellTopology> & cell_topology)
  {
    return createIntrepid2Basis<DeviceType,OutputValueType,PointValueType>(basis_type,basis_order,*cell_topology);
  }

}

// Instantiate for the device of every execution space Kokkos has enabled.
// What callers ask for has to be what is instantiated here, and that depends
// on what PHX::Device is.
//
// The DEVICE form is always needed.  Intrepid2's first parameter is a device,
// so Basis<Kokkos::Serial> and Basis<Kokkos::Device<Serial,HostSpace>> are
// unrelated types; callers pass PHX::Device, and code that spells a device
// explicitly -- a host device for a host-only path, say -- asks for one in
// either configuration.  SPACE::device_type is
// Kokkos::Device<SPACE, SPACE::memory_space>, spelled that way because it
// carries no comma and so survives macro argument splitting.
//
// The bare EXECUTION SPACE form is needed only while PHX::Device is still an
// execution space, so it is added just for that build.  SPACE and
// SPACE::device_type are always distinct types, so the two can never collide.
//
// Each backend macro names a distinct execution space, hence a distinct
// device, so no two of these can be the same specialization.  Naming
// PHX::Device directly instead would require knowing which backend it is,
// which the preprocessor cannot work out.  PHX::Device is among this list only
// while its memory space is the execution space's own; a shared memory space
// makes it distinct from all of them, which is what the guarded instantiation
// at the end of this file is for.
#define PANZER_INSTANTIATE_INTREPID2_BASIS_FOR(DEV)                           \
  template Teuchos::RCP<Intrepid2::Basis<DEV,double,double> >                 \
  panzer::createIntrepid2Basis<DEV,double,double>(                            \
      const std::string, int, const shards::CellTopology &);                  \
  template Teuchos::RCP<Intrepid2::Basis<DEV,double,double> >                 \
  panzer::createIntrepid2Basis<DEV,double,double>(                            \
      const std::string, int, const Teuchos::RCP<const shards::CellTopology> &);

#if defined(PHX_DEPRECATED_DEVICE_AS_EXECUTION_SPACE)
#define PANZER_INSTANTIATE_INTREPID2_BASIS(SPACE)                             \
  PANZER_INSTANTIATE_INTREPID2_BASIS_FOR(SPACE::device_type)                  \
  PANZER_INSTANTIATE_INTREPID2_BASIS_FOR(SPACE)
#else
#define PANZER_INSTANTIATE_INTREPID2_BASIS(SPACE)                             \
  PANZER_INSTANTIATE_INTREPID2_BASIS_FOR(SPACE::device_type)
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
PANZER_INSTANTIATE_INTREPID2_BASIS(Kokkos::Serial)
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
PANZER_INSTANTIATE_INTREPID2_BASIS(Kokkos::OpenMP)
#endif
#if defined(KOKKOS_ENABLE_THREADS)
PANZER_INSTANTIATE_INTREPID2_BASIS(Kokkos::Threads)
#endif
#if defined(KOKKOS_ENABLE_CUDA)
PANZER_INSTANTIATE_INTREPID2_BASIS(Kokkos::Cuda)
#endif
#if defined(KOKKOS_ENABLE_HIP)
PANZER_INSTANTIATE_INTREPID2_BASIS(Kokkos::HIP)
#endif
#if defined(KOKKOS_ENABLE_SYCL)
PANZER_INSTANTIATE_INTREPID2_BASIS(Kokkos::SYCL)
#endif

// The list above covers PHX::Device only while its memory space is its
// execution space's own, which is what SPACE::device_type means.  A shared
// memory space breaks that: Kokkos::Device<Cuda,CudaUVMSpace> is not
// Kokkos::Device<Cuda,CudaSpace>, so PHX::Device needs its own instantiation.
//
// Phalanx computes PHX_MEMORY_SPACE_IS_NOT_EXECUTION_SPACE_MEMORY_SPACE for
// exactly this question, and because the memory space is always either the
// execution space's own or Kokkos::SharedSpace, it is exact.  Not needed in
// the deprecated mode, where PHX::Device is an execution space and so already
// instantiated above.  The static_assert is a backstop: nothing should be able
// to reach it.
#if defined(PHX_MEMORY_SPACE_IS_NOT_EXECUTION_SPACE_MEMORY_SPACE) &&          \
    !defined(PHX_DEPRECATED_DEVICE_AS_EXECUTION_SPACE)
static_assert(!std::is_same<PHX::MemorySpace,
                            PHX::ExecutionSpace::memory_space>::value,
              "panzer: PHX::Device would repeat one of the instantiations "
              "above.  Narrow the condition on this block -- the configured "
              "memory space is the execution space's own, so PHX::Device is "
              "already covered by its backend's SPACE::device_type.");
PANZER_INSTANTIATE_INTREPID2_BASIS_FOR(PHX::Device)
#endif

#undef PANZER_INSTANTIATE_INTREPID2_BASIS
#undef PANZER_INSTANTIATE_INTREPID2_BASIS_FOR
