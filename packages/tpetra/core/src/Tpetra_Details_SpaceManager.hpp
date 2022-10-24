#ifndef TPETRA_SPACEMANAGER_HPP
#define TPETRA_SPACEMANAGER_HPP

#include <vector>
#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_Spaces.hpp"
#include "Teuchos_RCP.hpp"

namespace Tpetra {
namespace Details {



/* Store RCPs for a particular execution space, acquiring them from the instance
   manager lazily
*/
template <typename ExecSpace>
class SpaceSlot {
public:
    using execution_space = ExecSpace;

    SpaceSlot() {
        for (int i = 0; i < static_cast<int>(Spaces::Priority::NUM_LEVELS); ++i) {
            // retrieve an RCP from the global instance manager
            instances_[i] = Spaces::space_instance<execution_space>(
                static_cast<Spaces::Priority>(i)
            );
        }
    }

    template<Spaces::Priority priority = Spaces::Priority::medium>
    Teuchos::RCP<const execution_space> space_instance() const {
        return instances_[static_cast<int>(priority)];
    }

    Teuchos::RCP<const execution_space> space_instance(const Spaces::Priority &priority) const {
        switch(priority) {
            case Spaces::Priority::high: return space_instance<Spaces::Priority::high>();
            case Spaces::Priority::medium: return space_instance<Spaces::Priority::medium>();
            case Spaces::Priority::low: return space_instance<Spaces::Priority::low>();
            default: throw std::runtime_error("unexpected Tpetra Space priority");
        }
    }

private:
    // an instance of each possible priority
    Teuchos::RCP<const execution_space> instances_[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
}; // SpaceSlot

/*! \brief A class can inherit from this if it wants to use Tpetra managed spaces 
*/
class SpaceUser {
public:

    /* Get execution space i for a given priority

        In some algorithms, the desired priority of independent operations
        may be statically known

        Catch-all when we don't implement priority for spaces (non-CUDA)
    */
#ifdef KOKKOS_ENABLE_CUDA
    template <
      typename ExecSpace, 
      Spaces::Priority priority = Spaces::Priority::medium, 
      Spaces::IsCuda<ExecSpace> = true
    >
    Teuchos::RCP<const ExecSpace> space_instance() const {
        return cudaSlot.space_instance<priority>();
    }
    template <
      typename ExecSpace, 
      Spaces::IsCuda<ExecSpace> = true 
    >
    Teuchos::RCP<const ExecSpace> space_instance(const Spaces::Priority &priority) const {
        return cudaSlot.space_instance(priority);
    }
#endif // KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_SERIAL
    template <
      typename ExecSpace, 
      Spaces::Priority priority = Spaces::Priority::medium, 
      Spaces::IsSerial<ExecSpace> = true
    >
    Teuchos::RCP<const ExecSpace> space_instance() const {
        return serialSlot.space_instance<priority>();
    }
    template <
      typename ExecSpace, 
      Spaces::IsSerial<ExecSpace> = true 
    >
    Teuchos::RCP<const ExecSpace> space_instance(const Spaces::Priority &priority) const {
        return serialSlot.space_instance(priority);
    }
#endif // KOKKOS_ENABLE_SERIAL
#ifdef KOKKOS_ENABLE_OPENMP
    template <
      typename ExecSpace, 
      Spaces::Priority priority = Spaces::Priority::medium, 
      Spaces::IsOpenMP<ExecSpace> = true
    >
    Teuchos::RCP<const ExecSpace> space_instance() const {
        return openMPSlot.space_instance<priority>();
    }
    template <
      typename ExecSpace, 
      Spaces::IsOpenMP<ExecSpace> = true 
    >
    Teuchos::RCP<const ExecSpace> space_instance(const Spaces::Priority &priority) const {
        return openMPSlot.space_instance(priority);
    }
#endif // KOKKOS_ENABLE_OPENMP


#ifdef KOKKOS_ENABLE_CUDA
    SpaceSlot<Kokkos::Cuda> cudaSlot;
#endif
#ifdef KOKKOS_ENABLE_SERIAL
    SpaceSlot<Kokkos::Serial> serialSlot;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    SpaceSlot<Kokkos::OpenMP> openMPSlot;
#endif


}; // SpaceUser
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_SPACEMANAGER_HPP
