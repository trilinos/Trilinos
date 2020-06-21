#ifndef _COMPADRE_PARALLELMANAGER_HPP_
#define _COMPADRE_PARALLELMANAGER_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"

namespace Compadre {

//!  Parallel Manager
/*!
*  This class sets and manages thread / teams levels, scratch memory sizes, and kernel executions.
*  ex:
*  Compadre::ConvertLayoutLeftToRight clr;
*  Compadre::ParallelManager pm;
*  // no tag specified
*  pm.CallFunctorWithTeamThreads(100, clr, "MyFunctorName");
*  // some tag specified
*  pm.CallFunctorWithTeamThreads<DefaultTag>(100, clr);
*/
class ParallelManager {
protected:

    //! lowest level memory for Kokkos::parallel_for for team access memory
    int _scratch_team_level_a;
    int _team_scratch_size_a;

    //! higher (slower) level memory for Kokkos::parallel_for for team access memory
    int _scratch_thread_level_a;
    int _thread_scratch_size_a;

    //! lowest level memory for Kokkos::parallel_for for thread access memory
    int _scratch_team_level_b;
    int _team_scratch_size_b;

    //! higher (slower) level memory for Kokkos::parallel_for for thread access memory
    int _scratch_thread_level_b;
    int _thread_scratch_size_b;

    //! calculated number of threads per team
    int _threads_per_team;
    int _vector_lanes_per_thread;


/** @name Private Modifiers
 *  Private function because information lives on the device
 */
///@{
///@}

/** @name Private Accessors
 *  Private function because information lives on the device
 */
///@{
///@}

/** @name Private Utility
 *  
 */
///@{
///@}

public:

/** @name Instantiation / Destruction
 *  
 */
///@{

    ParallelManager() : _team_scratch_size_a(0), _thread_scratch_size_a(0), 
            _team_scratch_size_b(0), _thread_scratch_size_b(0) {

#ifdef COMPADRE_USE_CUDA
        _scratch_team_level_a = 0;
        _scratch_thread_level_a = 1;
        _scratch_team_level_b = 1;
        _scratch_thread_level_b = 1;
        _threads_per_team = 128;
        _vector_lanes_per_thread = 8;
#else
        _scratch_team_level_a = 0;
        _scratch_thread_level_a = 0;
        _scratch_team_level_b = 0;
        _scratch_thread_level_b = 0;
        _threads_per_team = 1;
        _vector_lanes_per_thread = 1;
#endif
    }

///@}

/** @name Public Utility
 *  
 */
///@{
///@}

/** @name Accessors
 *  Retrieve member variables through public member functions
 */
///@{

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each vector lane executing code be default
    template<typename Tag, class C>
    void CallFunctorWithTeamThreadsAndVectors(const global_index_type batch_size, const int threads_per_team, 
            const int vector_lanes_per_thread, C functor) const {

        if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                // all levels of each type need specified separately
                Kokkos::parallel_for(
                    Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor, typeid(Tag).name());
        } else if (_scratch_team_level_a != _scratch_team_level_b) {
            // scratch thread levels are the same
            Kokkos::parallel_for(
                Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                functor, typeid(Tag).name());
        } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
            // scratch team levels are the same
            Kokkos::parallel_for(
                Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                functor, typeid(Tag).name());
        } else {
            // scratch team levels and thread levels are the same
            Kokkos::parallel_for(
                Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                functor, typeid(Tag).name());
        }

    }

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each vector lane executing code be default
    template<class C>
    void CallFunctorWithTeamThreadsAndVectors(const global_index_type batch_size, const int threads_per_team, 
            const int vector_lanes_per_thread, C functor, std::string functor_name = typeid(C).name()) const {

        if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
            // all levels of each type need specified separately
            Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                functor, functor_name);
        } else if (_scratch_team_level_a != _scratch_team_level_b) {
            // scratch thread levels are the same
            Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                functor, functor_name);
        } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
            // scratch team levels are the same
            Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                functor, functor_name);
        } else {
            // scratch team levels and thread levels are the same
            Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                functor, functor_name);
        }

    }

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each thread executing code be default
    template<typename Tag, class C>
    void CallFunctorWithTeamThreads(const global_index_type batch_size, C functor) const {
        // calls breakout over vector lanes with vector lane size of 1
        CallFunctorWithTeamThreadsAndVectors<Tag,C>(batch_size, this->getThreadsPerTeam(), 1, functor);
    }

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each thread executing code be default
    template<class C>
    void CallFunctorWithTeamThreads(const global_index_type batch_size, C functor, std::string functor_name = typeid(C).name()) const {
        // calls breakout over vector lanes with vector lane size of 1
        CallFunctorWithTeamThreadsAndVectors<C>(batch_size, this->getThreadsPerTeam(), 1, functor, functor_name);
    }

    KOKKOS_INLINE_FUNCTION
    int getTeamScratchLevel(const int level) const {
        if (level == 0) {
            return _scratch_team_level_a;
        } else {
            return _scratch_team_level_b;
        }
    }

    KOKKOS_INLINE_FUNCTION
    int getThreadScratchLevel(const int level) const {
        if (level == 0) {
            return _scratch_thread_level_a;
        } else {
            return _scratch_thread_level_b;
        }
    }

    KOKKOS_INLINE_FUNCTION
    int getTeamScratchSize(const int level) const {
        if (level == 0) {
            return _team_scratch_size_a;
        } else {
            return _team_scratch_size_b;
        }
    }

    KOKKOS_INLINE_FUNCTION
    int getThreadScratchSize(const int level) const {
        if (level == 0) {
            return _thread_scratch_size_a;
        } else {
            return _thread_scratch_size_b;
        }
    }
    
    KOKKOS_INLINE_FUNCTION
    int getThreadsPerTeam(const int vector_lanes_per_thread = 1) const {
        return _threads_per_team / vector_lanes_per_thread;
    }

    KOKKOS_INLINE_FUNCTION
    int getVectorLanesPerThread() const {
        return _vector_lanes_per_thread;
    }

///@}


/** @name Modifiers
 *  Changed member variables through public member functions
 */
///@{

    void setTeamScratchLevel(const int level, const int value) {
        if (level == 0) {
            _scratch_team_level_a = value;
        } else {
            _scratch_team_level_b = value;
        }
    }

    void setThreadScratchLevel(const int level, const int value) {
        if (level == 0) {
            _scratch_thread_level_a = value;
        } else {
            _scratch_thread_level_b = value;
        }
    }

    void setTeamScratchSize(const int level, const int value) {
        if (level == 0) {
            _team_scratch_size_a = value;
        } else {
            _team_scratch_size_b = value;
        }
    }

    void setThreadScratchSize(const int level, const int value) {
        if (level == 0) {
            _thread_scratch_size_a = value;
        } else {
            _thread_scratch_size_b = value;
        }
    }

    void clearScratchSizes() {
        _team_scratch_size_a = 0;
        _team_scratch_size_b = 0;
        _thread_scratch_size_a = 0;
        _thread_scratch_size_b = 0;
    }

    void setThreadsPerTeam(const int value) {
        _threads_per_team = value;
    }

    void setVectorLanesPerThread(const int value) {
        _vector_lanes_per_thread = value;
    }

///@}


}; // ParallelManager Class
} // Compadre

#endif


