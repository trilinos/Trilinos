//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_HEXELEMENT_HPP
#define KOKKOS_HEXELEMENT_HPP

namespace Kokkos {
namespace Example {

template< unsigned NodeCount >
class HexElement_TensorData ;

template< unsigned NodeCount , class Device >
class HexElement_TensorEval ;

//----------------------------------------------------------------------------
/** \brief  Evaluate Hex element on interval [-1,1]^3 */
template<>
class HexElement_TensorData< 8 > {
public:

  static const unsigned element_node_count    = 8 ;
  static const unsigned spatial_dimension     = 3 ;
  static const unsigned integration_count_1d  = 2 ;
  static const unsigned function_count_1d     = 2 ;

  double values_1d [ function_count_1d ][ integration_count_1d ];
  double derivs_1d [ function_count_1d ][ integration_count_1d ];
  double weights_1d[ integration_count_1d ];

  unsigned char eval_map[ element_node_count ][4] ;

  static double eval_value_1d( const unsigned jf , const double x )
  {
    return 0 == jf ? 0.5 * ( 1.0 - x ) : (
           1 == jf ? 0.5 * ( 1.0 + x ) : 0 );
  }

  static double eval_deriv_1d( const unsigned jf , const double )
  {
    return 0 == jf ? -0.5 : (
           1 == jf ?  0.5 : 0 );
  }

  HexElement_TensorData()
  {
    const unsigned char tmp_map[ element_node_count ][ spatial_dimension ] =
      { { 0 , 0 , 0 },
        { 1 , 0 , 0 },
        { 1 , 1 , 0 },
        { 0 , 1 , 0 },
        { 0 , 0 , 1 },
        { 1 , 0 , 1 },
        { 1 , 1 , 1 },
        { 0 , 1 , 1 } };

    weights_1d[0] = 1 ;
    weights_1d[1] = 1 ;

    const double points_1d[ integration_count_1d ] =
      { -0.577350269189623 , 0.577350269189623 };

    for ( unsigned i = 0 ; i < element_node_count ; ++i ) {
      eval_map[i][0] = tmp_map[i][0];
      eval_map[i][1] = tmp_map[i][1];
      eval_map[i][2] = tmp_map[i][2];
    }

    for ( unsigned xp = 0 ; xp < integration_count_1d ; ++xp ) {
    for ( unsigned xf = 0 ; xf < function_count_1d ; ++xf ) {
      values_1d[xp][xf] = eval_value_1d( xf , points_1d[xp] );
      derivs_1d[xp][xf] = eval_deriv_1d( xf , points_1d[xp] );
    }}
  }
};

//----------------------------------------------------------------------------

template<>
class HexElement_TensorData< 27 > {
public:

  static const unsigned element_node_count    = 27 ;
  static const unsigned spatial_dimension     = 3 ;
  static const unsigned integration_count_1d  = 3 ;
  static const unsigned function_count_1d     = 3 ;

  double values_1d [ function_count_1d ][ integration_count_1d ];
  double derivs_1d [ function_count_1d ][ integration_count_1d ];
  double weights_1d[ integration_count_1d ];

  unsigned char eval_map[ element_node_count ][4] ;

  // sizeof(EvaluateElementHex) = 111 bytes =
  //   sizeof(double) * 9 +
  //   sizeof(double) * 9 +
  //   sizeof(double) * 3 +
  //   sizeof(char) * 27 

  static double eval_value_1d( const unsigned jf , const double p )
  {
    return 0 == jf ? 0.5 * p * ( p - 1 ) : (
           1 == jf ? 1.0 - p * p : (
           2 == jf ? 0.5 * p * ( p + 1 ) : 0 ));
  }

  static double eval_deriv_1d( const unsigned jf , const double p )
  {
    return 0 == jf ? p - 0.5 : (
           1 == jf ? -2.0 * p : (
           2 == jf ? p + 0.5 : 0 ));
  }

  HexElement_TensorData()
  {
    const unsigned char tmp_map[ element_node_count ][ spatial_dimension ] =
      { { 0 , 0 , 0 },
        { 2 , 0 , 0 },
        { 2 , 2 , 0 },
        { 0 , 2 , 0 },
        { 0 , 0 , 2 },
        { 2 , 0 , 2 },
        { 2 , 2 , 2 },
        { 0 , 2 , 2 },
        { 1 , 0 , 0 },
        { 2 , 1 , 0 },
        { 1 , 2 , 0 },
        { 0 , 1 , 0 },
        { 0 , 0 , 1 },
        { 2 , 0 , 1 },
        { 2 , 2 , 1 },
        { 0 , 2 , 1 },
        { 1 , 0 , 2 },
        { 2 , 1 , 2 },
        { 1 , 2 , 2 },
        { 0 , 1 , 2 },
        { 1 , 1 , 1 },
        { 1 , 1 , 0 },
        { 1 , 1 , 2 },
        { 0 , 1 , 1 },
        { 2 , 1 , 1 },
        { 1 , 0 , 1 },
        { 1 , 2 , 1 } };

    // Interval [-1,1]

    weights_1d[0] = 0.55555555555556 ;
    weights_1d[1] = 0.88888888888889 ;
    weights_1d[2] = 0.55555555555556 ;

    const double points_1d[3] = { -0.774596669241483 ,
                                   0.000000000000000 ,
                                   0.774596669241483 };

    for ( unsigned i = 0 ; i < element_node_count ; ++i ) {
      eval_map[i][0] = tmp_map[i][0];
      eval_map[i][1] = tmp_map[i][1];
      eval_map[i][2] = tmp_map[i][2];
    }

    for ( unsigned xp = 0 ; xp < integration_count_1d ; ++xp ) {
    for ( unsigned xf = 0 ; xf < function_count_1d ; ++xf ) {
      values_1d[xp][xf] = eval_value_1d( xf , points_1d[xp] );
      derivs_1d[xp][xf] = eval_deriv_1d( xf , points_1d[xp] );
    }}
  }
};

//----------------------------------------------------------------------------

template< unsigned NodeCount >
class HexElement_Data {
public:
  static const unsigned spatial_dimension   = 3 ;
  static const unsigned element_node_count  = NodeCount ;
  static const unsigned integration_count   = NodeCount ;
  static const unsigned function_count      = NodeCount ;

  double weights[   integration_count ] ;
  double values[    integration_count ][ function_count ];
  double gradients[ integration_count ][ spatial_dimension ][ function_count ];

  HexElement_Data()
  {
    HexElement_TensorData< NodeCount > tensor_data ;

    for ( unsigned ip = 0 ; ip < integration_count ; ++ip ) {

      const unsigned ipx = tensor_data.eval_map[ip][0] ;
      const unsigned ipy = tensor_data.eval_map[ip][1] ;
      const unsigned ipz = tensor_data.eval_map[ip][2] ;

      weights[ip] = tensor_data.weights_1d[ ipx ] *
                    tensor_data.weights_1d[ ipy ] *
                    tensor_data.weights_1d[ ipz ] ;

      for ( unsigned jf = 0 ; jf < function_count ; ++jf ) {

        const unsigned jfx = tensor_data.eval_map[jf][0] ;
        const unsigned jfy = tensor_data.eval_map[jf][1] ;
        const unsigned jfz = tensor_data.eval_map[jf][2] ;

        values[ip][jf] = tensor_data.values_1d[ ipx ][ jfx ] *
                         tensor_data.values_1d[ ipy ][ jfy ] *
                         tensor_data.values_1d[ ipz ][ jfz ] ;

        gradients[ip][0][jf] = tensor_data.derivs_1d[ ipx ][ jfx ] *
                               tensor_data.values_1d[ ipy ][ jfy ] *
                               tensor_data.values_1d[ ipz ][ jfz ] ;

        gradients[ip][1][jf] = tensor_data.values_1d[ ipx ][ jfx ] *
                               tensor_data.derivs_1d[ ipy ][ jfy ] *
                               tensor_data.values_1d[ ipz ][ jfz ] ;

        gradients[ip][2][jf] = tensor_data.values_1d[ ipx ][ jfx ] *
                               tensor_data.values_1d[ ipy ][ jfy ] *
                               tensor_data.derivs_1d[ ipz ][ jfz ] ;
      }
    }
  }
};

//----------------------------------------------------------------------------

} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_HEXELEMENT_HPP */


