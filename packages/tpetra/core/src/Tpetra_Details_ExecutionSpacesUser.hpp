// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_EXECUTIONSPACESUSER_HPP
#define TPETRA_DETAILS_EXECUTIONSPACESUSER_HPP

#include <iostream>
#include <sstream>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Teuchos_RCP.hpp"

#include "Tpetra_Details_ExecutionSpacesSlot.hpp"

namespace Tpetra {
namespace Details {
namespace Spaces {

/*! @brief A class can inherit from this if it wants to use Tpetra managed
 * spaces

   Adds a Tpetra::Details::SpaceSlot for each enabled Kokkos execution space,
   and provides a ::space_instance() member function that class members can
   call to retrieve the owned spaces, if they want to use them.
 */
class User {
public:
#ifdef KOKKOS_ENABLE_CUDA

  /**
   * @brief Returns a const smart pointer to an execution space instance.
   *
   * @tparam ExecSpace the execution space type
   * @tparam priority the execution space priority, default value: medium
   *
   * @return a const smart pointer to an execution space instance
   */
  template <typename ExecSpace,
            Spaces::Priority priority = Spaces::Priority::medium,
            Spaces::IsCuda<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace> space_instance() const {
    return cudaSlot.space_instance<priority>();
  }

  /**
   * @brief Returns a const smart pointer to an execution space instance with a
   * specific priority.
   *
   * @tparam ExecSpace the execution space type
   * @param priority the execution space priority
   *
   * @return a const smart pointer to an execution space instance
   */
  template <typename ExecSpace, Spaces::IsCuda<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace>
  space_instance(const Spaces::Priority &priority) const {
    return cudaSlot.space_instance(priority);
  }
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SERIAL
  template <typename ExecSpace,
            Spaces::Priority priority = Spaces::Priority::medium,
            Spaces::IsSerial<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace> space_instance() const {
    return serialSlot.space_instance<priority>();
  }
  template <typename ExecSpace, Spaces::IsSerial<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace>
  space_instance(const Spaces::Priority &priority) const {
    return serialSlot.space_instance(priority);
  }
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_OPENMP
  template <typename ExecSpace,
            Spaces::Priority priority = Spaces::Priority::medium,
            Spaces::IsOpenMP<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace> space_instance() const {
    return openMPSlot.space_instance<priority>();
  }
  template <typename ExecSpace, Spaces::IsOpenMP<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace>
  space_instance(const Spaces::Priority &priority) const {
    return openMPSlot.space_instance(priority);
  }
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_HIP
  template <typename ExecSpace,
            Spaces::Priority priority = Spaces::Priority::medium,
            Spaces::IsHIP<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace> space_instance() const {
    return HIPSlot.space_instance<priority>();
  }
  template <typename ExecSpace, Spaces::IsHIP<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace>
  space_instance(const Spaces::Priority &priority) const {
    return HIPSlot.space_instance(priority);
  }
#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL
  template <typename ExecSpace,
            Spaces::Priority priority = Spaces::Priority::medium,
            Spaces::IsSYCL<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace> space_instance() const {
    return SYCLSlot.space_instance<priority>();
  }
  template <typename ExecSpace, Spaces::IsSYCL<ExecSpace> = true>
  Teuchos::RCP<const ExecSpace>
  space_instance(const Spaces::Priority &priority) const {
    return SYCLSlot.space_instance(priority);
  }
#endif // KOKKOS_ENABLE_SYCL

#ifdef KOKKOS_ENABLE_SERIAL
  Slot<Kokkos::Serial> serialSlot;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  Slot<Kokkos::OpenMP> openMPSlot;
#endif
#ifdef KOKKOS_ENABLE_CUDA
  Slot<Kokkos::Cuda> cudaSlot;
#endif
#ifdef KOKKOS_ENABLE_HIP
  Slot<Kokkos::HIP> HIPSlot;
#endif
#ifdef KOKKOS_ENABLE_SYCL
  Slot<Kokkos::Experimental::SYCL> SYCLSlot;
#endif

}; // User
} // namespace Spaces
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_EXECUTIONSPACESUSER_HPP
