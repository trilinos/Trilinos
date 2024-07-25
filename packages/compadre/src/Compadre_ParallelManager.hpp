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

    //! largest team size 
    int _default_threads;
    int _default_vector_lanes;


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

        _default_threads = 16;
        _default_vector_lanes = 8;
#else
        _scratch_team_level_a = 0;
        _scratch_thread_level_a = 0;
        _scratch_team_level_b = 0;
        _scratch_thread_level_b = 0;

        _default_threads = 1;
        _default_vector_lanes = 1;
#endif
        if (const char* env_threads = std::getenv("THREADS")) {
            _default_threads = std::atoi(env_threads);
        }
        if (const char* env_vector_lanes = std::getenv("VECTORLANES")) {
            _default_vector_lanes = std::atoi(env_vector_lanes);
        }
#ifdef COMPADRE_EXTREME_DEBUG
        printf("threads per team: %d, vector lanes per team: %d\n", _default_threads, _default_vector_lanes);
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

        if (threads_per_team>0 && vector_lanes_per_thread>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                // all levels of each type need specified separately
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else {
                // scratch team levels and thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            }
        } else if (threads_per_team>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                // all levels of each type need specified separately
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else {
                // scratch team levels and thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            }
        } else if (vector_lanes_per_thread>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                // all levels of each type need specified separately
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else {
                // scratch team levels and thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            }
        } else {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                // all levels of each type need specified separately
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b));
            } else {
                // scratch team levels and thread levels are the same
                return Kokkos::TeamPolicy<device_execution_space>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b));
            }
        }
    }

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each vector lane executing code be default
    template<typename Tag, class C>
    void CallFunctorWithTeamThreadsAndVectors(C functor, const global_index_type batch_size, const int threads_per_team = -1, 
            const int vector_lanes_per_thread = -1) const {

        if (threads_per_team>0 && vector_lanes_per_thread>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                    // all levels of each type need specified separately
                    Kokkos::parallel_for(
                        typeid(Tag).name(),
                        Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                        .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                        .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                        .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                        .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                        functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        } else if (threads_per_team>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                    // all levels of each type need specified separately
                    Kokkos::parallel_for(
                        typeid(Tag).name(),
                        Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, _default_vector_lanes)
                        .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                        .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                        .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                        .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                        functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        } else if (vector_lanes_per_thread>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                    // all levels of each type need specified separately
                    Kokkos::parallel_for(
                        typeid(Tag).name(),
                        Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, vector_lanes_per_thread)
                        .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                        .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                        .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                        .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                        functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        } else {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                    // all levels of each type need specified separately
                    Kokkos::parallel_for(
                        typeid(Tag).name(),
                        Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, _default_vector_lanes)
                        .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                        .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                        .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                        .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                        functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    typeid(Tag).name(),
                    Kokkos::TeamPolicy<Tag>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        }
    }

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each vector lane executing code be default
    template<class C>
    void CallFunctorWithTeamThreadsAndVectors(C functor, const global_index_type batch_size, const int threads_per_team = -1, 
            const int vector_lanes_per_thread = -1, std::string functor_name = typeid(C).name()) const {

        if (threads_per_team>0 && vector_lanes_per_thread>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                // all levels of each type need specified separately
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        } else if (threads_per_team>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                // all levels of each type need specified separately
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, threads_per_team, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        } else if (vector_lanes_per_thread>0) {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, vector_lanes_per_thread)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        } else {
            if ( (_scratch_team_level_a != _scratch_team_level_b) && (_scratch_thread_level_a != _scratch_thread_level_b) ) {
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else if (_scratch_team_level_a != _scratch_team_level_b) {
                // scratch thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a))
                    .set_scratch_size(_scratch_team_level_b, Kokkos::PerTeam(_team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            } else if (_scratch_thread_level_a != _scratch_thread_level_b) {
                // scratch team levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a))
                    .set_scratch_size(_scratch_thread_level_b, Kokkos::PerThread(_thread_scratch_size_b)),
                    functor);
            } else {
                // scratch team levels and thread levels are the same
                Kokkos::parallel_for(
                    functor_name,
                    Kokkos::TeamPolicy<>(batch_size, _default_threads, _default_vector_lanes)
                    .set_scratch_size(_scratch_team_level_a, Kokkos::PerTeam(_team_scratch_size_a + _team_scratch_size_b))
                    .set_scratch_size(_scratch_thread_level_a, Kokkos::PerThread(_thread_scratch_size_a + _thread_scratch_size_b)),
                    functor);
            }
        }
    }

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each thread executing code be default
    template<typename Tag, class C>
    void CallFunctorWithTeamThreads(C functor, const global_index_type batch_size) const {
        // calls breakout over vector lanes with vector lane size of 1
        CallFunctorWithTeamThreadsAndVectors<Tag,C>(functor, batch_size, _default_threads, 1);
    }

    //! Calls a parallel_for
    //! parallel_for will break out over loops over teams with each thread executing code be default
    template<class C>
    void CallFunctorWithTeamThreads(C functor, const global_index_type batch_size, std::string functor_name = typeid(C).name()) const {
        // calls breakout over vector lanes with vector lane size of 1
        CallFunctorWithTeamThreadsAndVectors<C>(functor, batch_size, _default_threads, 1, functor_name);
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


