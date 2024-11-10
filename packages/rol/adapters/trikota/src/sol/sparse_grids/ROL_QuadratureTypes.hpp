// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATURETYPES_HPP
#define ROL_QUADRATURETYPES_HPP

namespace ROL {

  /** \enum  ROL::EQuadrature
      \brief Enumerations of integration rules provided in Dakota.
  */
  enum EQuadrature {
    QUAD_CHEBYSHEV1 = 0,
    QUAD_CHEBYSHEV2,
    QUAD_CLENSHAWCURTIS,
    QUAD_FEJER2,
    QUAD_LEGENDRE,
    QUAD_PATTERSON,
    QUAD_TRAPEZOIDAL,
    QUAD_HERMITE,
    QUAD_GENZKEISTER,
    QUAD_LAGUERRE,
    QUAD_LAST
  };

  inline std::string EQuadratureToString(EQuadrature rule) {
    std::string retString;
    switch(rule) {
      case QUAD_CHEBYSHEV1:       retString = "Gauss-Chebyshev Type 1";    break;
      case QUAD_CHEBYSHEV2:       retString = "Gauss-Chebyshev Type 2";    break;
      case QUAD_CLENSHAWCURTIS:   retString = "Clenshaw-Curtis";           break;
      case QUAD_FEJER2:           retString = "Fejer Type 2";              break;
      case QUAD_LEGENDRE:         retString = "Gauss-Legendre";            break;
      case QUAD_PATTERSON:        retString = "Gauss-Patterson";           break;
      case QUAD_TRAPEZOIDAL:      retString = "Trapezoidal Rule";          break;
      case QUAD_HERMITE:          retString = "Gauss-Hermite";             break;
      case QUAD_GENZKEISTER:      retString = "Hermite-Genz-Keister";      break;
      case QUAD_LAGUERRE:         retString = "Gauss-Laguerre";            break;
      case QUAD_LAST:             retString = "Last Type (Dummy)";         break;
      default:                    retString = "INVALID EQuadrature";
    }
    return retString;
  }

  inline int isValidQuadrature(EQuadrature rule) {
    return( (rule == QUAD_CHEBYSHEV1) ||
            (rule == QUAD_CHEBYSHEV2) ||
            (rule == QUAD_CLENSHAWCURTIS) ||
            (rule == QUAD_FEJER2) ||
            (rule == QUAD_LEGENDRE) ||
            (rule == QUAD_PATTERSON) ||
            (rule == QUAD_TRAPEZOIDAL) ||
            (rule == QUAD_HERMITE) ||
            (rule == QUAD_GENZKEISTER) ||
            (rule == QUAD_LAGUERRE) );
  }

  inline EQuadrature & operator++(EQuadrature &type) {
    return type = static_cast<EQuadrature>(type+1);
  }

  inline EQuadrature operator++(EQuadrature &type, int) {
    EQuadrature oldval = type;
    ++type;
    return oldval;
  }

  inline EQuadrature & operator--(EQuadrature &type) {
    return type = static_cast<EQuadrature>(type-1);
  }

  inline EQuadrature operator--(EQuadrature &type, int) {
    EQuadrature oldval = type;
    --type;
    return oldval;
  }

  inline EQuadrature StringToEQuadrature(std::string s) {
    s = removeStringFormat(s);
    for ( EQuadrature q = QUAD_CHEBYSHEV1; q < QUAD_LAST; q++ ) {
      if ( !s.compare(removeStringFormat(EQuadratureToString(q))) ) {
        return q;
      }
    }
    return QUAD_CHEBYSHEV1;
  }

  enum EGrowth {
    GROWTH_DEFAULT = 0,
    GROWTH_SLOWLIN,
    GROWTH_SLOWLINODD,
    GROWTH_MODLIN,
    GROWTH_SLOWEXP,
    GROWTH_MODEXP,
    GROWTH_FULLEXP,
    GROWTH_LAST
  };

  inline int isValidGrowth(EGrowth rule) {
    return( (rule == GROWTH_DEFAULT) ||
            (rule == GROWTH_SLOWLIN) ||
            (rule == GROWTH_SLOWLINODD) ||
            (rule == GROWTH_MODLIN) ||
            (rule == GROWTH_SLOWEXP) ||
            (rule == GROWTH_MODEXP) ||
            (rule == GROWTH_FULLEXP) );
  }

  inline std::string EGrowthToString(EGrowth rule) {
    std::string retString;
    switch(rule) {
      case GROWTH_DEFAULT:    retString = "Default";              break;
      case GROWTH_SLOWLIN:    retString = "Slow Linear";          break;
      case GROWTH_SLOWLINODD: retString = "Slow Linear Odd";      break;
      case GROWTH_MODLIN:     retString = "Moderate Linear";      break;
      case GROWTH_SLOWEXP:    retString = "Slow Exponential";     break;
      case GROWTH_MODEXP:     retString = "Moderate Exponential"; break;
      case GROWTH_FULLEXP:    retString = "Fully Exponential";    break;
      case GROWTH_LAST:       retString = "Last Type (Dummy)";    break;
      default:                retString = "INVALID EGrowth";
    }
    return retString;
  }

  inline EGrowth & operator++(EGrowth &type) {
    return type = static_cast<EGrowth>(type+1);
  }

  inline EGrowth operator++(EGrowth &type, int) {
    EGrowth oldval = type;
    ++type;
    return oldval;
  }

  inline EGrowth & operator--(EGrowth &type) {
    return type = static_cast<EGrowth>(type-1);
  }

  inline EGrowth operator--(EGrowth &type, int) {
    EGrowth oldval = type;
    --type;
    return oldval;
  }

  inline EGrowth StringToEGrowth(std::string s) {
    s = removeStringFormat(s);
    for ( EGrowth q = GROWTH_DEFAULT; q < GROWTH_LAST; q++ ) {
      if ( !s.compare(removeStringFormat(EGrowthToString(q))) ) {
        return q;
      }
    }
    return GROWTH_DEFAULT;
  }

  struct QuadratureInfo {
    int dim;
    int maxLevel;
    std::vector<EQuadrature> rule1D;
    std::vector<EGrowth> growth1D;
    bool normalized;
    bool adaptive;
    bool print;
    std::string name;
    QuadratureInfo(void)
      : maxLevel(1), normalized(true), adaptive(false),
        print(false), name("Default") {}
  };

  inline int growthRule1D(int index, EGrowth growth, EQuadrature rule) {
    //
    //  Compute the growth sequence for 1D quadrature rules according to growth.
    //  For more information on growth rules, see 
    //  
    //  J. Burkardt. 1D Quadrature Rules For Sparse Grids.
    //  http://people.sc.fsu.edu/~jburkardt/presentations/sgmga_1d_rules.pdf.
    //  
    //  Drew P. Kouri
    //  Sandia National Laboratories - CSRI
    //  May 27, 2011
    //
  
    int level = index-1;
    //int level = index;
    if (rule==QUAD_CLENSHAWCURTIS) { // Clenshaw-Curtis
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        if (level==0) {
          return 1;
        }
        else { 
          int o = 2;
          while(o<2*level+1) {
            o = 2*(o-1)+1;
          }
          return o;
        }
      }
      else if (growth==GROWTH_MODEXP) {
        if (level==0) {
          return 1;
        }
        else {
          int o = 2;
          while (o<4*level+1) {
            o = 2*(o-1)+1;
          }
          return o;
        }
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        if (level==0) {
          return 1;
        }
        else {
          return std::pow(2,level)+1;
        }
      }
    }
    else if (rule==QUAD_FEJER2) { // Fejer Type 2
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        int o = 1;
        while (o<2*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_MODEXP) {
        int o = 1;
        while (o<4*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        return std::pow(2,level+1)-1;
      }
    }
 
    else if (rule==QUAD_PATTERSON) { // Gauss-Patterson
      if (growth==GROWTH_SLOWLIN||
          growth==GROWTH_SLOWLINODD||
          growth==GROWTH_MODLIN) {
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>> (ROL::growthRule1D): Specified Growth Rule Not Permitted!");
        return 0;
      }
      else if (growth==GROWTH_SLOWEXP) {
        if (level==0) {
          return 1;
        }
        else {
        int p = 5;
        int o = 3;
        while (p<2*level+1) {
          p = 2*p+1;
          o = 2*o+1;
        }
        return o;
        }
      }
      else if (growth==GROWTH_MODEXP) {
        if (level==0) {
          return 1;
        }
        else {
          int p = 5;
          int o = 3;
          while (p<4*level+1) {
            p = 2*p+1;
            o = 2*o+1;
          }
          return o;
        }
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        return std::pow(2,level+1)-1;
      }
    }
  
    else if (rule==QUAD_LEGENDRE) { // Gauss-Legendre
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        int o = 1;
        while (2*o-1<2*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_MODEXP) {
        int o = 1;
        while (2*o-1<4*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        return std::pow(2,level+1)-1;
      }
    }
  
    else if (rule==QUAD_HERMITE) { // Gauss-Hermite
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        int o = 1;
        while (2*o-1<2*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_MODEXP) {
        int o = 1;
        while (2*o-1<4*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        return std::pow(2,level+1)-1;
      }
    }
    
    else if (rule==QUAD_LAGUERRE) { // Gauss-Laguerre
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        int o = 1;
        while (2*o-1<2*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_MODEXP) {
        int o = 1;
        while (2*o-1<4*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        return std::pow(2,level+1)-1;
      }
    }
  
    else if (rule==QUAD_CHEBYSHEV1) { // Gauss-Chebyshev Type 1
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        int o = 1;
        while (2*o-1<2*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_MODEXP) {
        int o = 1;
        while (2*o-1<4*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        return std::pow(2,level+1)-1;
      }
    }
  
  
    else if (rule==QUAD_CHEBYSHEV2) { // Gauss-Chebyshev Type 2
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        int o = 1;
        while (2*o-1<2*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_MODEXP) {
        int o = 1;
        while (2*o-1<4*level+1) {
          o = 2*o+1;
        }
        return o;
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        return std::pow(2,level+1)-1;
      }
    }
    
    else if (rule==QUAD_GENZKEISTER) { // Hermite-Genz-Keister  
      static int o_hgk[5] = { 1, 3, 9, 19, 35 };
      static int p_hgk[5] = { 1, 5, 15, 29, 51 };
      if (growth==GROWTH_SLOWLIN||
          growth==GROWTH_SLOWLINODD||
          growth==GROWTH_MODLIN) {
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>> (ROL::growthRule1D): Specified Growth Rule Not Permitted!");
        return 0;
      }
      else if (growth==GROWTH_SLOWEXP) { 
        int l = 0, p = p_hgk[l], o = o_hgk[l];
        while (p<2*level+1 && l<4) {
          l++;
          p = p_hgk[l];
          o = o_hgk[l];
        }
        return o;
      }
      else if (growth==GROWTH_MODEXP) {
        int l = 0, p = p_hgk[l], o = o_hgk[l];
        while (p<4*level+1 && l<4) {
          l++;
          p = p_hgk[l];
          o = o_hgk[l];
        }
        return o;
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        int l = level; l = std::max(l,0); l = std::min(l,4);
        return o_hgk[l];
      }
    }  
  
    else if (rule==QUAD_TRAPEZOIDAL) { // Trapezoidal
      if (growth==GROWTH_SLOWLIN) {
        return level+1;
      }
      else if (growth==GROWTH_SLOWLINODD) {
        return 2*((level+1)/2)+1;
      }
      else if (growth==GROWTH_MODLIN) {
        return 2*level+1;
      }
      else if (growth==GROWTH_SLOWEXP) {
        if (level==0) {
          return 1;
        }
        else { 
          int o = 2;
          while(o<2*level+1) {
            o = 2*(o-1)+1;
          }
          return o;
        }
      }
      else if (growth==GROWTH_MODEXP) {
        if (level==0) {
          return 1;
        }
        else {
          int o = 2;
          while (o<4*level+1) {
            o = 2*(o-1)+1;
          }
          return o;
        }
      }
      else if (growth==GROWTH_FULLEXP||growth==GROWTH_DEFAULT) {
        if (level==0) {
          return 1;
        }
        else {
          return std::pow(2,level)+1;
        }
      }
    }
    return 0;
  } // end growthRule1D
} // end ROL namespace

#endif
