#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGP_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGP_H_

#include <Kokkos_Core.hpp>

#include <ngp/NgpSpaces.hpp>
#include <ngp/NgpMesh.hpp>
#include <ngp/NgpField.hpp>
#include <ngp/NgpAtomics.hpp>
#include <ngp/NgpForEachEntity.hpp>

namespace ngp {

#ifdef KOKKOS_HAVE_CUDA
typedef ngp::StaticMesh StkNgpMesh;
typedef ngp::StaticField<double> StkNgpField;
#else
typedef ngp::WrapperMesh StkNgpMesh;
//typedef ngp::StaticField<double> StkNgpField;
typedef ngp::WrapperField<double> StkNgpField;
#endif

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGP_H_ */
