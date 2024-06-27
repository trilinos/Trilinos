// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_CHEBYSHEV_WEIGHTS_HPP
#define IFPACK2_DETAILS_CHEBYSHEV_WEIGHTS_HPP

#include <stdexcept>
#include <vector>
#include "Teuchos_StandardParameterEntryValidators.hpp"

/// \file Ifpack2_Details_Chebyshev_Weights.hpp
/// \brief Definition of Chebyshev implementation
/// \author Malachi Phillips
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It defines the weights needed for fourth kind Chebyshev polynomials
/// with optimal weights. See: https://arxiv.org/pdf/2202.08830.pdf for details.


namespace Ifpack2 {
namespace Details {

  /// \brief Generate optimal weights for using the fourth kind Chebyshev polynomials
  ///   see:
  ///    https://arxiv.org/pdf/2202.08830.pdf
  ///   for details.
  ///
  /// \pre 0 < chebyOrder <= 16 -- this is an arbitrary distinction,
  ///      but the weights are currently hard-coded from the
  ///      MATLAB scripts from https://arxiv.org/pdf/2202.08830.pdf.
  ///
  /// \param chebyOrder [in] Number of Chebyshev iterations.
  template<typename ScalarType>
  std::vector<ScalarType>
  optimalWeightsImpl (const int chebyOrder)
  {
    if (chebyOrder == 0){
      return {};
    }
    if (chebyOrder == 1){
      return {
          1.12500000000000,
      };
    }

    if (chebyOrder == 2){
      return {
          1.02387287570313,
          1.26408905371085,
      };
    }

    if (chebyOrder == 3){
      return {
          1.00842544782028,
          1.08867839208730,
          1.33753125909618,
      };
    }

    if (chebyOrder == 4){
      return {
          1.00391310427285,
          1.04035811188593,
          1.14863498546254,
          1.38268869241000,
      };
    }

    if (chebyOrder == 5){
      return {
          1.00212930146164,
          1.02173711549260,
          1.07872433192603,
          1.19810065292663,
          1.41322542791682,
      };
    }

    if (chebyOrder == 6){
      return {
          1.00128517255940,
          1.01304293035233,
          1.04678215124113,
          1.11616489419675,
          1.23829020218444,
          1.43524297106744,
      };
    }

    if (chebyOrder == 7){
      return {
          1.00083464397912,
          1.00843949430122,
          1.03008707768713,
          1.07408384092003,
          1.15036186707366,
          1.27116474046139,
          1.45186658649364,
      };
    }

    if (chebyOrder == 8){
      return {
          1.00057246631197,
          1.00577427662415,
          1.02050187922941,
          1.05019803444565,
          1.10115572984941,
          1.18086042806856,
          1.29838585382576,
          1.46486073151099,
      };
    }

    if (chebyOrder == 9){
      return {
          1.00040960072832,
          1.00412439506106,
          1.01460212148266,
          1.03561113626671,
          1.07139972529194,
          1.12688273710962,
          1.20785219140729,
          1.32121930716746,
          1.47529642820699,
      };
    }

    if (chebyOrder == 10){
      return {
          1.00030312229652,
          1.00304840660796,
          1.01077022715387,
          1.02619011597640,
          1.05231724933755,
          1.09255743207549,
          1.15083376663972,
          1.23172250870894,
          1.34060802024460,
          1.48386124407011,
      };
    }

    if (chebyOrder == 11){
      return {
          1.00023058595209,
          1.00231675024028,
          1.00817245396304,
          1.01982986566342,
          1.03950210235324,
          1.06965042700541,
          1.11305754295742,
          1.17290876275564,
          1.25288300576792,
          1.35725579919519,
          1.49101672564139,
      };
    }

    if (chebyOrder == 12){
      return {
          1.00017947200828,
          1.00180189139619,
          1.00634861907307,
          1.01537864566306,
          1.03056942830760,
          1.05376019693943,
          1.08699862592072,
          1.13259183097913,
          1.19316273358172,
          1.27171293675110,
          1.37169337969799,
          1.49708418575562,
      };
    }

    if (chebyOrder == 13){
      return {
          1.00014241921559,
          1.00142906932629,
          1.00503028986298,
          1.01216910518495,
          1.02414874342792,
          1.04238158880820,
          1.06842008128700,
          1.10399010936759,
          1.15102748242645,
          1.21171811910125,
          1.28854264865128,
          1.38432619380991,
          1.50229418757368,
      };
    }

    if (chebyOrder == 14){
      return {
          1.00011490538261,
          1.00115246376914,
          1.00405357333264,
          1.00979590573153,
          1.01941300472994,
          1.03401425035436,
          1.05480599606629,
          1.08311420301813,
          1.12040891660892,
          1.16833095655446,
          1.22872122288238,
          1.30365305707817,
          1.39546814053678,
          1.50681646209583,
      };
    }

    if (chebyOrder == 15){
      return {
          1.00009404750752,
          1.00094291696343,
          1.00331449056444,
          1.00800294833816,
          1.01584236259140,
          1.02772083317705,
          1.04459535422831,
          1.06750761206125,
          1.09760092545889,
          1.13613855366157,
          1.18452361426236,
          1.24432087304475,
          1.31728069083392,
          1.40536543893560,
          1.51077872501845,
      };
    }

    if (chebyOrder == 16){
      return {
          1.00007794828179,
          1.00078126847253,
          1.00274487974401,
          1.00662291017015,
          1.01309858836971,
          1.02289448329337,
          1.03678321409983,
          1.05559875719896,
          1.08024848405560,
          1.11172607131497,
          1.15112543431072,
          1.19965584614973,
          1.25865841744946,
          1.32962412656664,
          1.41421360695576,
          1.51427891730346,
      };
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "Ifpack2::Details::optimalWeightsImpl::"
        "Requested Chebyshev order exceeds maximum of 16."
      );
  }

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CHEBYSHEV_WEIGHTS_HPP
