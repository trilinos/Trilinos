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
using StkNgpMesh = StaticMesh;
template <typename T> using StkNgpField = ngp::StaticField<T>;
#else
using StkNgpMesh = ngp::WrapperMesh;
template <typename T> using StkNgpField = ngp::WrapperField<T>;
//template <typename T> using StkNgpField = ngp::StaticField<T>;
#endif

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGP_H_ */
