// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_POINTINTERPOLATION_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_POINTINTERPOLATION_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/FieldState.hpp"  // for FieldState, StateNP1, StateNM1
#include "stk_mesh/base/Types.hpp"       // for EntityRank
#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_search_util/MasterElementProvider.hpp"
#include "stk_transfer/TransferTypes.hpp"
#include "stk_transfer_util/FieldUtility.hpp"
#include "stk_transfer_util/LeastSquaresInterpolation.hpp"
#include "stk_transfer_util/Patch.hpp"                      // for LinearPatch
#include "stk_transfer_util/LeastSquares.hpp"               // for Geometric...
#include "stk_transfer_util/LeastSquaresInterpolation.hpp"  // for EntityInt...
#include "stk_transfer_util/RecoverField.hpp"               // for RecoverField
#include "stk_util/util/string_case_compare.hpp"

#include <functional>                    // for function
#include <limits>                        // for numeric_limits
#include <string>                        // for string, operator==, basic_st...
#include <utility>                       // for pair
#include <vector>                        // for vector
namespace stk::mesh { class BulkData; }
namespace stk::mesh { class FieldBase; }
namespace stk::mesh { class MetaData; }
namespace stk::mesh { class Selector; }
namespace stk::mesh { struct Entity; }
namespace stk::mesh { class Part; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

namespace impl {
const stk::mesh::Selector& universal_selector();
}

void nodal_interpolation(stk::mesh::Entity entity,
                         const stk::search::SearchField& field, const std::vector<double>& paramCoords,
                         const stk::search::MasterElementProviderInterface& masterElement, const double lowerBound,
                         const double upperBound, std::vector<double>& fieldDataScratch, std::vector<double>& result);

void nodal_interpolation(stk::mesh::Entity entity,
                         const stk::mesh::FieldBase& field, const std::vector<double>& paramCoords,
                         const stk::search::MasterElementProviderInterface& masterElement, const double lowerBound,
                         const double upperBound, std::vector<double>& fieldDataScratch, std::vector<double>& result);

void nodal_interpolation(stk::mesh::Entity entity,
                         const stk::mesh::FieldBase& field, const std::vector<double>& paramCoords,
                         const stk::search::MasterElementProviderInterface& masterElement, const double lowerBound,
                         const double upperBound, std::vector<double>& result);

template<typename PATCHFILTER>
bool element_centroid_linear(const PATCHFILTER& patchFilter,
                             const stk::mesh::Selector& activeSelector,
                             const stk::transfer::EntityInterpolationData& interpData,
                             const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  stk::transfer::LinearPatch<PATCHFILTER> patch(interpData.sendBulk, interpData.sendEntity, patchFilter, activeSelector);

  const int numSampPtsInPatch = patch.size();
  int nComponents = 1;
  const auto basisSize = static_cast<unsigned>(stk::transfer::RecoverField::TRILINEAR);

  stk::transfer::LeastSquares leastSquaresCalculator(nComponents, numSampPtsInPatch, basisSize);

  return stk::transfer::least_squares_linear_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                           interpData, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_linear(const PATCHFILTER& patchFilter,
                             const stk::transfer::EntityInterpolationData& interpData,
                             const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  return element_centroid_linear<PATCHFILTER>(patchFilter, impl::universal_selector(), interpData, evalPoint, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_linear_mls(const PATCHFILTER& patchFilter,
                                 const stk::mesh::Selector& activeSelector,
                                 const stk::transfer::EntityInterpolationData& interpData,
                                 const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  stk::transfer::LinearPatch<PATCHFILTER> patch(interpData.sendBulk, interpData.sendEntity, patchFilter, activeSelector);

  const int numSampPtsInPatch = patch.size();
  int nComponents = 1;
  const auto basisSize = static_cast<unsigned>(stk::transfer::RecoverField::TRILINEAR);

  stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(nComponents, numSampPtsInPatch, basisSize,
                                                                    interpData.sendBulk, patch.get_patch_entities(),
                                                                    interpData.sendNodalCoordField, evalPoint);

  return stk::transfer::least_squares_linear_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                           interpData, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_linear_mls(const PATCHFILTER& patchFilter,
                                 const stk::transfer::EntityInterpolationData& interpData,
                                 const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  return element_centroid_linear_mls<PATCHFILTER>(patchFilter, impl::universal_selector(), interpData, evalPoint, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_quadratic(const PATCHFILTER& patchFilter,
                                const stk::mesh::Selector& activeSelector,
                                const stk::transfer::EntityInterpolationData& interpData,
                                const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  stk::transfer::QuadraticPatch<PATCHFILTER> patch(interpData.sendBulk, interpData.sendEntity, patchFilter, activeSelector);

  const int numSampPtsInPatch = patch.size();
  int nComponents = 1;
  const auto basisSize = static_cast<unsigned>(stk::transfer::RecoverField::TRIQUADRATIC);

  stk::transfer::LeastSquares leastSquaresCalculator(nComponents, numSampPtsInPatch, basisSize);

  return stk::transfer::least_squares_quadratic_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                              interpData, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_quadratic(const PATCHFILTER& patchFilter,
                                const stk::transfer::EntityInterpolationData& interpData,
                                const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  return element_centroid_quadratic<PATCHFILTER>(patchFilter, impl::universal_selector(), interpData, evalPoint, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_quadratic_mls(const PATCHFILTER& patchFilter,
                                    const stk::mesh::Selector& activeSelector,
                                    const stk::transfer::EntityInterpolationData& interpData,
                                    const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  stk::transfer::QuadraticPatch<PATCHFILTER> patch(interpData.sendBulk, interpData.sendEntity, patchFilter, activeSelector);

  const int numSampPtsInPatch = patch.size();
  int nComponents = 1;
  const auto basisSize = static_cast<unsigned>(stk::transfer::RecoverField::TRIQUADRATIC);

  stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(nComponents, numSampPtsInPatch, basisSize,
                                                                    interpData.sendBulk, patch.get_patch_entities(),
                                                                    interpData.sendNodalCoordField, evalPoint);

  return stk::transfer::least_squares_quadratic_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                              interpData, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_quadratic_mls(const PATCHFILTER& patchFilter,
                                    const stk::transfer::EntityInterpolationData& interpData,
                                    const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  return element_centroid_quadratic_mls<PATCHFILTER>(patchFilter, impl::universal_selector(), interpData, evalPoint, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_cubic(const PATCHFILTER& patchFilter,
                            const stk::mesh::Selector& activeSelector,
                            const stk::transfer::EntityInterpolationData& interpData,
                            const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  stk::transfer::CubicPatch<PATCHFILTER> patch(interpData.sendBulk, interpData.sendEntity, patchFilter, activeSelector);

  const int numSampPtsInPatch = patch.size();
  int nComponents = 1;
  const auto basisSize = static_cast<unsigned>(stk::transfer::RecoverField::TRICUBIC);

  stk::transfer::LeastSquares leastSquaresCalculator(nComponents, numSampPtsInPatch, basisSize);

  return stk::transfer::least_squares_cubic_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                          interpData, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_cubic(const PATCHFILTER& patchFilter,
                            const stk::transfer::EntityInterpolationData& interpData,
                            const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  return element_centroid_cubic<PATCHFILTER>(patchFilter, impl::universal_selector(), interpData, evalPoint, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_cubic_mls(const PATCHFILTER& patchFilter,
                                const stk::mesh::Selector& activeSelector,
                                const stk::transfer::EntityInterpolationData& interpData,
                                const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  stk::transfer::CubicPatch<PATCHFILTER> patch(interpData.sendBulk, interpData.sendEntity, patchFilter, activeSelector);

  const int numSampPtsInPatch = patch.size();
  int nComponents = 1;
  const auto basisSize = static_cast<unsigned>(stk::transfer::RecoverField::TRICUBIC);

  stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(nComponents, numSampPtsInPatch, basisSize,
                                                                    interpData.sendBulk, patch.get_patch_entities(),
                                                                    interpData.sendNodalCoordField, evalPoint);

  return stk::transfer::least_squares_cubic_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                              interpData, toFieldSize, toFieldPtr);
}

template<typename PATCHFILTER>
bool element_centroid_cubic_mls(const PATCHFILTER& patchFilter,
                                const stk::transfer::EntityInterpolationData& interpData,
                                const std::vector<double>& evalPoint, unsigned toFieldSize, double* toFieldPtr)
{
  return element_centroid_cubic_mls<PATCHFILTER>(patchFilter, impl::universal_selector(), interpData, evalPoint, toFieldSize, toFieldPtr);
}

} // namespace transfer
} // namespace stk

#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_POINTINTERPOLATION_HPP_ */
