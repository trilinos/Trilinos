// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_HIERARCHIC_PARALLELISM_HPP
#define PANZER_HIERARCHIC_PARALLELISM_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

  /// Singleton class for accessing kokkos hierarchical parallelism parameters.
  class HP {
    bool use_auto_team_size_; /// If true, the team size is set with Kokkos::AUTO()
    int team_size_;           /// User specified team size.
    int vector_size_;         /// Default vector size for non-AD types.
    int fad_vector_size_;     /// FAD vector size.
    bool use_shared_memory_;     /// Use shared memory kokkos kernels for non-fad types.
    bool fad_use_shared_memory_; /// Use shared memory kokkos kernels for fad types.

    HP(); /// Private ctor.

  public:
    /// Return singleton instance of this class.
    static HP& inst();

    /** Allows the user to override the Kokkos default team and vector
        sizes for kernel dispatch. The values will be capped by
        hardware limits and rounded down to the nearest power of two.

        The final variable will force the values input to be set explicity and
        not round down to the nearest power of two or hardware maximum.

        @param team_size Team size requested for hierarchic kernel
        @param vector_size Vector size requested for hierarchic kernel for non-FAD scalar types
        @param fad_vector_size Vector size requested for hierarchic kernel for FAD scalar types
        @param force_override_safety Ignore the power of two and other checks
    */
    void overrideSizes(const int& team_size,
		       const int& vector_size,
		       const int& fad_vector_size, 
                       const bool force_override_safety=false);

    /// \brief Reset the sizes to default
    void resetSizes() { use_auto_team_size_=true;}

    /** \brief Returns the vector size. Specialized for AD scalar types.

	NOTE: For hierarchic parallelism, if we use the same code for
	both Residual and Jacobian (as we do in most evaluators), the
	loop over vector level is missing for Residual. The loop is
	implemented internally in the AD types for Jacobian where on
	CUDA the warp parallelizes over the derivative dimension. To
	prevent incorrect code, we need to force the vector size to 1
	for non-AD scalar types. Eventual workaround is to use SIMD
	data type with similar hidden vector loop for Residual. In the
	mean time, this function will set correct vector_size of one.
     */
    template<typename Scalar>
    int vectorSize() const
    {
      return Sacado::IsADType<Scalar>::value ? fad_vector_size_ : vector_size_;
    }

    /** \brief Tell kokkos kernels if they should use shared memory. This is very problem dependent.

	If a panzer hierarchic kernel can use shared memory to speed
	the calculation, then it carries a second implementation that
	takes advantage of shared memory. Shared memory on the GPU is
	very limited. On some of the example problems, the shared
	memory runs out if the basis is greated than order 2 on a hex
	mesh. This is also very dependent on the size of the
	derivative array. A large derivative array uses up memory much
	quicker. The default is that for non-fad types, we always
	enable shared memory. For fad types the default is to disable
	use of shared memory, but this function can override for
	specific problems. For example, the
	adapters-stk/examples/MixedPoission problem can use shared
	memory for fad types for basis order 2 or less. It will call
	this function based on the basis order to improve performance.
     */
    void setUseSharedMemory(const bool& use_shared_memory,
			    const bool& fad_use_shared_memory);

    template<typename Scalar>
    bool useSharedMemory() const
    {
      return Sacado::IsADType<Scalar>::value ? fad_use_shared_memory_ : use_shared_memory_;
    }

    /// Returns a TeamPolicy for hierarchic parallelism.
    template<typename ScalarT, typename ... TeamPolicyProperties>
    Kokkos::TeamPolicy<TeamPolicyProperties...> teamPolicy(const int& league_size)
    {
      const int tmp_vector_size = this->template vectorSize<ScalarT>();

      if (use_auto_team_size_)
        return Kokkos::TeamPolicy<TeamPolicyProperties...>(league_size,Kokkos::AUTO(),
                                                           tmp_vector_size);

      return Kokkos::TeamPolicy<TeamPolicyProperties...>(league_size,team_size_,tmp_vector_size);
    }

    /// Returns a TeamPolicy for hierarchic parallelism using an exec_space instance (for cuda streams).
    template<typename ScalarT, typename ... TeamPolicyProperties, typename ExecSpace>
    Kokkos::TeamPolicy<ExecSpace, TeamPolicyProperties...> teamPolicy(ExecSpace exec_space, const int& league_size)
    {
      const int tmp_vector_size = this->template vectorSize<ScalarT>();

      if (use_auto_team_size_)
        return Kokkos::TeamPolicy<TeamPolicyProperties...>(exec_space,league_size,Kokkos::AUTO(),
                                                           tmp_vector_size);

      return Kokkos::TeamPolicy<TeamPolicyProperties...>(exec_space,league_size,team_size_,tmp_vector_size);
    }
  };

}

#endif
