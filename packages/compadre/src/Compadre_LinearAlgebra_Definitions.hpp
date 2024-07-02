// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_LINEAR_ALGEBRA_DEFINITIONS_HPP_
#define _COMPADRE_LINEAR_ALGEBRA_DEFINITIONS_HPP_

#include "Compadre_LinearAlgebra_Declarations.hpp"

namespace Compadre {
namespace GMLS_LinearAlgebra {

KOKKOS_INLINE_FUNCTION
void largestTwoEigenvectorsThreeByThreeSymmetric(const member_type& teamMember, scratch_matrix_right_type V, scratch_matrix_right_type PtP, const int dimensions, pool_type& random_number_pool) {

    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

        double maxRange = 100;

        generator_type rand_gen = random_number_pool.get_state();
        // put in a power method here and a deflation by first found eigenvalue
        double eigenvalue_relative_tolerance = 1e-6; // TODO: use something smaller, but really anything close is acceptable for this manifold

        
        double v[3] = {rand_gen.drand(maxRange),rand_gen.drand(maxRange),rand_gen.drand(maxRange)};
        double v_old[3] = {v[0], v[1], v[2]};

        double error = 1;
        double norm_v;

        while (error > eigenvalue_relative_tolerance) {

            double tmp1 = v[0];
            v[0] = PtP(0,0)*tmp1 + PtP(0,1)*v[1];
            if (dimensions>2) v[0] += PtP(0,2)*v[2];

            double tmp2 = v[1];
            v[1] = PtP(1,0)*tmp1 + PtP(1,1)*tmp2;
            if (dimensions>2) v[1] += PtP(1,2)*v[2];

            if (dimensions>2)
                v[2] = PtP(2,0)*tmp1 + PtP(2,1)*tmp2 + PtP(2,2)*v[2];

            norm_v = v[0]*v[0] + v[1]*v[1];
            if (dimensions>2) norm_v += v[2]*v[2];
            norm_v = std::sqrt(norm_v);

            v[0] = v[0] / norm_v;
            v[1] = v[1] / norm_v;
            if (dimensions>2) v[2] = v[2] / norm_v;

            error = (v[0]-v_old[0])*(v[0]-v_old[0]) + (v[1]-v_old[1])*(v[1]-v_old[1]);
            if (dimensions>2) error += (v[2]-v_old[2])*(v[2]-v_old[2]);
            error = std::sqrt(error);
            error /= norm_v;


            v_old[0] = v[0];
            v_old[1] = v[1];
            if (dimensions>2) v_old[2] = v[2];
        }

        double dot_product;
        double norm;

        // if 2D, orthonormalize second vector
        if (dimensions==2) {

            for (int i=0; i<2; ++i) {
                V(0,i) = v[i];
            }

            // orthonormalize second eigenvector against first
            V(1,0) = 1.0; V(1,1) = 1.0;
            dot_product = V(0,0)*V(1,0) + V(0,1)*V(1,1);
            V(1,0) -= dot_product*V(0,0);
            V(1,1) -= dot_product*V(0,1);

            norm = std::sqrt(V(1,0)*V(1,0) + V(1,1)*V(1,1));
            V(1,0) /= norm;
            V(1,1) /= norm;

        } else { // otherwise, work on second eigenvalue

            for (int i=0; i<3; ++i) {
                V(0,i) = v[i];
                for (int j=0; j<3; ++j) {
                    PtP(i,j) -= norm_v*v[i]*v[j];
                }
            }

            error = 1;
            v[0] = rand_gen.drand(maxRange); v[1] = rand_gen.drand(maxRange); v[2] = rand_gen.drand(maxRange);
            v_old[0] = v[0]; v_old[1] = v[1]; v_old[2] =v[2];
            while (error > eigenvalue_relative_tolerance) {

                double tmp1 = v[0];
                v[0] = PtP(0,0)*tmp1 + PtP(0,1)*v[1] + PtP(0,2)*v[2];

                double tmp2 = v[1];
                v[1] = PtP(1,0)*tmp1 + PtP(1,1)*tmp2 + PtP(1,2)*v[2];

                v[2] = PtP(2,0)*tmp1 + PtP(2,1)*tmp2 + PtP(2,2)*v[2];

                norm_v = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

                v[0] = v[0] / norm_v;
                v[1] = v[1] / norm_v;
                v[2] = v[2] / norm_v;

                error = std::sqrt((v[0]-v_old[0])*(v[0]-v_old[0]) + (v[1]-v_old[1])*(v[1]-v_old[1]) + (v[2]-v_old[2])*(v[2]-v_old[2])) / norm_v;

                v_old[0] = v[0];
                v_old[1] = v[1];
                v_old[2] = v[2];
            }

            for (int i=0; i<3; ++i) {
                V(1,i) = v[i];
            }

            // orthonormalize second eigenvector against first
            dot_product = V(0,0)*V(1,0) + V(0,1)*V(1,1) + V(0,2)*V(1,2);

            V(1,0) -= dot_product*V(0,0);
            V(1,1) -= dot_product*V(0,1);
            V(1,2) -= dot_product*V(0,2);

            norm = std::sqrt(V(1,0)*V(1,0) + V(1,1)*V(1,1) + V(1,2)*V(1,2));
            V(1,0) /= norm;
            V(1,1) /= norm;
            V(1,2) /= norm;

            // get normal from cross product
            V(2,0) = V(0,1)*V(1,2) - V(1,1)*V(0,2);
            V(2,1) = V(1,0)*V(0,2) - V(0,0)*V(1,2);
            V(2,2) = V(0,0)*V(1,1) - V(1,0)*V(0,1);

            // orthonormalize third eigenvector (just to be sure)
            norm = std::sqrt(V(2,0)*V(2,0) + V(2,1)*V(2,1) + V(2,2)*V(2,2));
            V(2,0) /= norm;
            V(2,1) /= norm;
            V(2,2) /= norm;

        }

        random_number_pool.free_state(rand_gen);
    });

}

} // GMLS_LinearAlgebra
} // Compadre

#endif

