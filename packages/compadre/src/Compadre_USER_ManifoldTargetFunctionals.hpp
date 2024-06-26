// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// this file picks up at the beginning of the computeTargetFunctionalsOnManifold function
#ifndef _USER_MANIFOLD_TARGET_FUNCTIONALS_HPP_
#define _USER_MANIFOLD_TARGET_FUNCTIONALS_HPP_

bool some_conditions_for_a_user_defined_operation = false;
bool some_conditions_for_another_user_defined_operation = false;

// hint: look in Compadre_GMLS_Target.hpp for examples

if (some_conditions_for_a_user_defined_operation) {
    // these operations are being called at the Team level,
    // so we call single to only perform the operation on one thread
    Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
        // user definition for a target functional goes here


    });
} else if (some_conditions_for_another_user_defined_operation) {
    // these operations are being called at the Team level,
    // so we call single to only perform the operation on one thread
    Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
        // user definition for a different target functional goes here


    });
} else {
    // if the operation was not caught by any user defined TargetFunctional,
    // then it is returned to the toolkit to try to handle the operation
    operation_handled = false;
}

#endif
