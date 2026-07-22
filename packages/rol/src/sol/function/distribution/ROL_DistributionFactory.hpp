// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DISTRIBUTIONFACTORY_HPP
#define ROL_DISTRIBUTIONFACTORY_HPP

#include "ROL_ParameterList.hpp"

#include "ROL_Dirac.hpp"
#include "ROL_Gaussian.hpp"
#include "ROL_TruncatedGaussian.hpp"
#include "ROL_Uniform.hpp"
#include "ROL_Logistic.hpp"
#include "ROL_Triangle.hpp"
#include "ROL_Parabolic.hpp"
#include "ROL_RaisedCosine.hpp"
#include "ROL_Laplace.hpp"
#include "ROL_Cauchy.hpp"
#include "ROL_Smale.hpp"
#include "ROL_Arcsine.hpp"
#include "ROL_Kumaraswamy.hpp"
#include "ROL_Exponential.hpp"
#include "ROL_TruncatedExponential.hpp"
#include "ROL_Gamma.hpp"
#include "ROL_Beta.hpp"

namespace ROL {

  enum EDistribution {
    DISTRIBUTION_ARCSINE = 0,
    DISTRIBUTION_BETA,
    DISTRIBUTION_CAUCHY,
    DISTRIBUTION_DIRAC, 
    DISTRIBUTION_EXPONENTIAL,
    DISTRIBUTION_GAMMA,
    DISTRIBUTION_GAUSSIAN,
    DISTRIBUTION_KUMARASWAMY,
    DISTRIBUTION_LAPLACE,
    DISTRIBUTION_LOGISTIC,
    DISTRIBUTION_PARABOLIC,
    DISTRIBUTION_RAISEDCOSINE,
    DISTRIBUTION_SMALE,
    DISTRIBUTION_TRIANGLE,
    DISTRIBUTION_TRUNCATEDEXPONENTIAL,
    DISTRIBUTION_TRUNCATEDGAUSSIAN,
    DISTRIBUTION_UNIFORM, 
    DISTRIBUTION_LAST
  };

  inline std::string EDistributionToString(EDistribution ed) {
    std::string retString;
    switch(ed) {
      case DISTRIBUTION_ARCSINE:              retString = "Arcsine";               break;
      case DISTRIBUTION_BETA:                 retString = "Beta";                  break;
      case DISTRIBUTION_CAUCHY:               retString = "Cauchy";                break;
      case DISTRIBUTION_DIRAC:                retString = "Dirac";                 break;
      case DISTRIBUTION_EXPONENTIAL:          retString = "Exponential";           break;
      case DISTRIBUTION_GAMMA:                retString = "Gamma";                 break;
      case DISTRIBUTION_GAUSSIAN:             retString = "Gaussian";              break;
      case DISTRIBUTION_KUMARASWAMY:          retString = "Kumaraswamy";           break;
      case DISTRIBUTION_LAPLACE:              retString = "Laplace";               break;
      case DISTRIBUTION_LOGISTIC:             retString = "Logistic";              break;
      case DISTRIBUTION_PARABOLIC:            retString = "Parabolic";             break;
      case DISTRIBUTION_RAISEDCOSINE:         retString = "Raised Cosine";         break;
      case DISTRIBUTION_SMALE:                retString = "Smale";                 break;
      case DISTRIBUTION_TRIANGLE:             retString = "Triangle";              break;
      case DISTRIBUTION_TRUNCATEDEXPONENTIAL: retString = "Truncated Exponential"; break;
      case DISTRIBUTION_TRUNCATEDGAUSSIAN:    retString = "Truncated Gaussian";    break;
      case DISTRIBUTION_UNIFORM:              retString = "Uniform";               break;
      case DISTRIBUTION_LAST:                 retString = "Last Type (Dummy)";     break;
      default:                                retString = "INVALID EDistribution"; break;
    }
    return retString;
  }

  inline int isValidDistribution(EDistribution ed) {
    return( (ed == DISTRIBUTION_DIRAC) ||
            (ed == DISTRIBUTION_BETA) ||
            (ed == DISTRIBUTION_GAMMA) ||
            (ed == DISTRIBUTION_GAUSSIAN) ||
            (ed == DISTRIBUTION_TRUNCATEDGAUSSIAN) ||
            (ed == DISTRIBUTION_UNIFORM) ||
            (ed == DISTRIBUTION_LOGISTIC) ||
            (ed == DISTRIBUTION_TRIANGLE) ||
            (ed == DISTRIBUTION_PARABOLIC) ||
            (ed == DISTRIBUTION_RAISEDCOSINE) ||
            (ed == DISTRIBUTION_LAPLACE) ||
            (ed == DISTRIBUTION_CAUCHY) ||
            (ed == DISTRIBUTION_SMALE) ||
            (ed == DISTRIBUTION_ARCSINE) ||
            (ed == DISTRIBUTION_KUMARASWAMY) || 
            (ed == DISTRIBUTION_EXPONENTIAL) ||
            (ed == DISTRIBUTION_TRUNCATEDEXPONENTIAL) );
  }

  inline EDistribution & operator++(EDistribution &type) {
    return type = static_cast<EDistribution>(type+1);
  }

  inline EDistribution operator++(EDistribution &type, int) {
    EDistribution oldval = type;
    ++type;
    return oldval;
  }

  inline EDistribution & operator--(EDistribution &type) {
    return type = static_cast<EDistribution>(type-1);
  }

  inline EDistribution operator--(EDistribution &type, int) {
    EDistribution oldval = type;
    --type;
    return oldval;
  }

  inline EDistribution StringToEDistribution(std::string s) {
    s = removeStringFormat(s);
    for ( EDistribution tr = DISTRIBUTION_ARCSINE; tr < DISTRIBUTION_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(EDistributionToString(tr))) ) {
        return tr;
      }
    }
    return DISTRIBUTION_UNIFORM;
  }

  template<class Real>
  inline ROL::Ptr<Distribution<Real> > DistributionFactory(ROL::ParameterList &parlist) {
    std::string dist;
    ROL::ParameterList sollist;
    if ( parlist.isSublist("SOL") ) {
      dist.assign(parlist.sublist("SOL").sublist("Distribution").get("Name","Dirac"));
      sollist = parlist;
    }
    else {
      dist.assign(parlist.sublist("Distribution").get("Name","Dirac"));
      sollist.sublist("SOL") = parlist;
    }
    EDistribution ed = StringToEDistribution(dist);
    switch(ed) {
      case DISTRIBUTION_ARCSINE:              return ROL::makePtr<Arcsine<Real>>(sollist);
      case DISTRIBUTION_BETA:                 return ROL::makePtr<Beta<Real>>(sollist);
      case DISTRIBUTION_CAUCHY:               return ROL::makePtr<Cauchy<Real>>(sollist);
      case DISTRIBUTION_DIRAC:                return ROL::makePtr<Dirac<Real>>(sollist);
      case DISTRIBUTION_EXPONENTIAL:          return ROL::makePtr<Exponential<Real>>(sollist);
      case DISTRIBUTION_GAMMA:                return ROL::makePtr<Gamma<Real>>(sollist);
      case DISTRIBUTION_GAUSSIAN:             return ROL::makePtr<Gaussian<Real>>(sollist);
      case DISTRIBUTION_KUMARASWAMY:          return ROL::makePtr<Kumaraswamy<Real>>(sollist);
      case DISTRIBUTION_LAPLACE:              return ROL::makePtr<Laplace<Real>>(sollist);
      case DISTRIBUTION_LOGISTIC:             return ROL::makePtr<Logistic<Real>>(sollist);
      case DISTRIBUTION_PARABOLIC:            return ROL::makePtr<Parabolic<Real>>(sollist);
      case DISTRIBUTION_RAISEDCOSINE:         return ROL::makePtr<RaisedCosine<Real>>(sollist);
      case DISTRIBUTION_SMALE:                return ROL::makePtr<Smale<Real>>(sollist);
      case DISTRIBUTION_TRIANGLE:             return ROL::makePtr<Triangle<Real>>(sollist);
      case DISTRIBUTION_TRUNCATEDEXPONENTIAL: return ROL::makePtr<TruncatedExponential<Real>>(sollist);
      case DISTRIBUTION_TRUNCATEDGAUSSIAN:    return ROL::makePtr<TruncatedGaussian<Real>>(sollist);
      case DISTRIBUTION_UNIFORM:              return ROL::makePtr<Uniform<Real>>(sollist);
      default:                                return ROL::nullPtr;
    }
  }
}
#endif
