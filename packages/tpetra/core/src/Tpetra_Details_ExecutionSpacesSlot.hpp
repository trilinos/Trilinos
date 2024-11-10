// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_EXEUTIONSPACESSLOT_HPP
#define TPETRA_DETAILS_EXEUTIONSPACESSLOT_HPP

#include <iostream>
#include <sstream>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_ExecutionSpaces.hpp"
#include <Kokkos_Core.hpp>

namespace Tpetra {
namespace Details {
namespace Spaces {

/**
 * @brief Lazily acquires and stores Kokkos Execution Spaces
 *
 * @tparam ExecSpace the type of the execution space to be wrapped.
 */
template <typename ExecSpace> class Slot {
public:
  using execution_space = ExecSpace;

  /**
   * @brief Default constructor that creates instances for each possible space
   * priority level.
   */
  Slot() {
    for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
      // retrieve an RCP from the global instance manager
      instances_[i] = Spaces::space_instance<execution_space>(
          static_cast<Spaces::Priority>(i));
    }
  }

  /**
   * @brief Get a specific execution space instance based on the given priority
   * level.
   *
   * @tparam priority the priority level of the desired execution space
   * instance. Default value is medium.
   * @return Teuchos::RCP<const execution_space> a smart pointer to a const
   * execution_space object.
   * @note This template method is used to get the execution space instance
   * based on the template parameter priority.
   */
  template <Spaces::Priority priority = Spaces::Priority::medium>
  Teuchos::RCP<const execution_space> space_instance() const {
    return instances_[static_cast<int>(priority)];
  }

  /**
   * @brief Get a specific execution space instance based on the given priority
   * level.
   *
   * @param priority the priority level of the desired execution space instance.
   * @return Teuchos::RCP<const execution_space> a smart pointer to a const
   * execution_space object.
   * @note This non-template method is used to get the execution space instance
   * based on the runtime parameter priority. It returns the corresponding
   * execution space instance by calling the corresponding template method.
   * Throws a runtime error if the given priority is not valid.
   */
  Teuchos::RCP<const execution_space>
  space_instance(const Spaces::Priority &priority) const {
    switch (priority) {
    case Spaces::Priority::high:
      return space_instance<Spaces::Priority::high>();
    case Spaces::Priority::medium:
      return space_instance<Spaces::Priority::medium>();
    case Spaces::Priority::low:
      return space_instance<Spaces::Priority::low>();
    default:
      throw std::runtime_error("unexpected Tpetra Space priority");
    }
  }

private:
  /**
   * @brief Array that contains instances of different execution space
   * prioritized by priority levels.
   */
  Teuchos::RCP<const execution_space>
      instances_[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
}; // Slot

} // namespace Spaces
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_EXEUTIONSPACESSLOT_HPP
