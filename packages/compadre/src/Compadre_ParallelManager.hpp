// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
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
*  pm.CallFunctorWithTeamThreads(clr, 100, "MyFunctorName");
*  // some tag specified
*  pm.CallFunctorWithTeamThreads<DefaultTag>(clr, 100);
*/
class ParallelManager {
public:

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

    //! number of threads and vector lanes
    //! if not specified through environment variables THREADS and VECTORLANES,
    //! then set to Kokkos::AUTO
    int _environment_threads;
    int _environment_vector_lanes;


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
        _scratch_thread_level_a = 0;
        _scratch_team_level_b = 1;
        _scratch_thread_level_b = 1;

        _environment_threads = -1;
        _environment_vector_lanes = -1;
#else
        _scratch_team_level_a = 0;
        _scratch_thread_level_a = 0;
        _scratch_team_level_b = 0;
        _scratch_thread_level_b = 0;

        _environment_threads = -1;
        _environment_vector_lanes = -1;
#endif
        if (const char* env_threads = std::getenv("THREADS")) {
            _environment_threads = std::atoi(env_threads);
        }
        if (const char* env_vector_lanes = std::getenv("VECTORLANES")) {
            _environment_vector_lanes = std::atoi(env_vector_lanes);
        }
#ifdef COMPADRE_EXTREME_DEBUG
        printf("ENVIRONMENT:: threads per team: %d, vector lanes per team: %d\n", _environment_threads, _environment_vector_lanes);
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

    //! Creates a team policy for a parallel_for
    //! parallel_for will break out over loops over teams with each vector lane executing code be default
    Kokkos::TeamPolicy<device_execution_space> 
        TeamPolicyThreadsAndVectors(const global_index_type batch_size, const int threads_per_team = -1, 
            const int vector_lanes_per_thread = -1) const {

        // first create object
        Kokkos::TeamPolicy<device_execution_space> tp;
        if (threads_per_team>0 && vector_lanes_per_thread>0) {
            tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, vector_lanes_per_thread);
        } else if (threads_per_team>0) {
            if (_environment_vector_lanes>0) {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, _environment_vector_lanes);
            } else {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, Kokkos::AUTO);
            }
        } else if (vector_lanes_per_thread>0) {
            if (_environment_threads>0) {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, _environment_threads, vector_lanes_per_thread);
            } else {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, Kokkos::AUTO, vector_lanes_per_thread);
            }
        } else {
            if (_environment_vector_lanes>0 && _environment_threads>0) {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, _environment_threads, _environment_vector_lanes);
            } else if (_environment_vector_lanes>0) {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, Kokkos::AUTO, _environment_vector_lanes);
            } else if (_environment_threads>0) {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, _environment_threads, Kokkos::AUTO);
            } else {
                tp = Kokkos::TeamPolicy<device_execution_space>(batch_size, Kokkos::AUTO, Kokkos::AUTO);
            }
        }
        
        if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
            tp.set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
              .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
              .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
              .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
        } else if (_scratch_team_level_a != _scratch_team_level_b) {
            tp.set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
              .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
              .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
        } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
            tp.set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
              .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
              .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
        } else {
            tp.set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
              .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
        }
        return tp;
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

///@}


}; // ParallelManager Class
} // Compadre

#endif


