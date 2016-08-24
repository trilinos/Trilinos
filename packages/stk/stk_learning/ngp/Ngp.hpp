#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGP_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGP_H_

#include <Kokkos_Core.hpp>

#include <ngp/NgpSpaces.hpp>
#include <ngp/NgpMesh.hpp>
#include <ngp/NgpField.hpp>
#include <ngp/NgpAtomics.hpp>
#include <ngp/NgpForEachEntity.hpp>
#include <ngp/NgpReductions.hpp>

namespace ngp {

#ifdef KOKKOS_HAVE_CUDA
using Mesh = StaticMesh;
template <typename T> using Field = ngp::StaticField<T>;
#else
using Mesh = ngp::StkMeshAdapter;
template <typename T> using Field = ngp::StkFieldAdapter<T>;
//template <typename T> using Field = ngp::StaticField<T>;
#endif

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGP_H_ */
