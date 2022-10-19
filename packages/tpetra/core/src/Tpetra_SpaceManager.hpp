#ifndef TPETRA_SPACEMANAGER_HPP
#define TPETRA_SPACEMANAGER_HPP

#include <vector>
#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Spaces.hpp"

namespace Tpetra {

/* A class can inherit from this if it wants to own execution space instances
   constructed at compile-time, and accesses are const
*/
class SpaceManager2 {
public:
    SpaceManager2();

    /* get Space i with priority prio
    */
    template <typename ExecSpace>
    const ExecSpace &space_instance(const Spaces::Priority &prio) const {
        switch(prio) {
            case Spaces::Priority::high: return space_instance<ExecSpace, Spaces::Priority::high>();
            case Spaces::Priority::medium: return space_instance<ExecSpace, Spaces::Priority::medium>();
            case Spaces::Priority::low: return space_instance<ExecSpace, Spaces::Priority::low>();
            default: throw std::runtime_error("unexpected Tpetra Space priority");
        }
    }

    /* Get execution space i for a given priority

        In some algorithms, the desired priority of independent operations
        may be statically known

        Catch-all when we don't implement priority for spaces (non-CUDA)
    */
#ifdef KOKKOS_ENABLE_CUDA
    template <typename ExecSpace, Spaces::Priority priority = Spaces::Priority::medium, 
    Spaces::detail::IsCuda<ExecSpace> = true >
    const ExecSpace &space_instance() const {
        return cudaSpace[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_SERIAL
    template <typename ExecSpace, Spaces::Priority priority = Spaces::Priority::medium, 
    Spaces::detail::IsSerial<ExecSpace> = true >
    const ExecSpace &space_instance() const {
        return serialSpace[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_SERIAL
#ifdef KOKKOS_ENABLE_OPENMP
    template <typename ExecSpace, Spaces::Priority priority = Spaces::Priority::medium, 
    Spaces::detail::IsOpenMP<ExecSpace> = true >
    const ExecSpace &space_instance() const {
        return openMPSpace[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_OPENMP


#ifdef KOKKOS_ENABLE_CUDA
    Kokkos::Cuda cudaSpace[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_SERIAL
    Kokkos::Serial serialSpace[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    Kokkos::OpenMP openMPSpace[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
#endif


}; // SpaceManager

/* A class can inherit from this if it wants to own execution space instances
   Constructed on-demand, which means that spaces cannot be accessed from const
   objects
*/
class SpaceManager {
public:
    ~SpaceManager();

    /* get Space i with priority prio
    */
    template <typename ExecSpace>
    ExecSpace &space_instance(int i, const Spaces::Priority &prio) {
        switch(prio) {
            case Spaces::Priority::high: return space_instance<ExecSpace, Spaces::Priority::high>(i);
            case Spaces::Priority::medium: return space_instance<ExecSpace, Spaces::Priority::medium>(i);
            case Spaces::Priority::low: return space_instance<ExecSpace, Spaces::Priority::low>(i);
            default: throw std::runtime_error("unexpected Tpetra Space priority");
        }
    }

    /* Get execution space i for a given priority

        In some algorithms, the desired priority of independent operations
        may be statically known

        Catch-all when we don't implement priority for spaces (non-CUDA)
    */
    template <typename ExecSpace, Spaces::Priority priority = Spaces::Priority::medium>
    ExecSpace &space_instance(int i = 0) {
        if (i < 0) {
            throw std::runtime_error("requested exec space < 0 from Spaces::get");
        }
        if (i >= Details::Behavior::spacesIdWarnLimit()) {
            std::cerr << "WARNING: requested space " << i << std::endl;
        }

        while (spaces<ExecSpace, priority>().size() <= size_t(i)) {
            spaces<ExecSpace, priority>().push_back(Spaces::make_instance<ExecSpace, priority>());
        }
        return spaces<ExecSpace, priority>()[i];
    }

protected:

#ifdef KOKKOS_ENABLE_CUDA
    template <typename ExecSpace, Spaces::Priority priority, 
    Spaces::detail::IsCuda<ExecSpace> = true >
    std::vector<ExecSpace> &spaces() {
        return cudaSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_SERIAL
    template <typename ExecSpace, Spaces::Priority priority, 
    Spaces::detail::IsSerial<ExecSpace> = true >
    std::vector<ExecSpace> &spaces() {
        return serialSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_SERIAL
#ifdef KOKKOS_ENABLE_OPENMP
    template <typename ExecSpace, Spaces::Priority priority, 
    Spaces::detail::IsOpenMP<ExecSpace> = true >
    std::vector<ExecSpace> &spaces() {
        return openMPSpaces[static_cast<int>(priority)];
    }
#endif // KOKKOS_ENABLE_OPENMP


#ifdef KOKKOS_ENABLE_CUDA
    std::vector<Kokkos::Cuda> cudaSpaces[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_SERIAL
    std::vector<Kokkos::Serial> serialSpaces[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    std::vector<Kokkos::OpenMP> openMPSpaces[static_cast<int>(Spaces::Priority::NUM_LEVELS)];
#endif


}; // SpaceManager



} // namespace Tpetra

#endif // TPETRA_SPACEMANAGER_HPP
