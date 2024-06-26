// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_QUADRATUREMANAGER_HPP_
#define _COMPADRE_QUADRATUREMANAGER_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"
#include "Compadre_Utilities.hpp"

namespace Compadre {

// ENUM for quadrature types
enum QuadratureType : int {INVALID, LINE, TRI, QUAD, TET, HEX};

//!  Quadrature
/*!
*  This class sets and manages quadrature orders, rules, etc...
*/
class Quadrature {
protected:

    bool _is_initialized;
    int _number_of_quadrature_points;
    int _order_of_quadrature_points;
    int _dimension_of_quadrature_points;
    Kokkos::View<double*, layout_right> _quadrature_weights;
    Kokkos::View<double**, layout_right> _parameterized_quadrature_sites;
    QuadratureType _qt;

/** @name Private Modifiers
 *  Private function because information lives on the device
 */
///@{

    //! Generates quadrature
    void generateQuadrature(const int order, const int dimension) {

        std::vector<int> oneDNumPointsFromOrder = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5};
        std::vector<int> twoDNumPointsFromOrder = {1, 1, 3, 4, 6, 7, 12, 13, 16, 19, 24, 27, 32, 36, 42, 55, 61, 66, 73, 78};

        if (dimension == 0) {
            _number_of_quadrature_points = 0;
        } else if (dimension == 1) {
            _number_of_quadrature_points = oneDNumPointsFromOrder[order];
        } else if (dimension == 2) {
            _number_of_quadrature_points = twoDNumPointsFromOrder[order];
        } else {
            compadre_assert_release(false && "Higher than 2D quadrature not currently supported.");
        }

        _quadrature_weights = decltype(_quadrature_weights)("quadrature weights", _number_of_quadrature_points);
        _parameterized_quadrature_sites = decltype(_parameterized_quadrature_sites)("quadrature sites", _number_of_quadrature_points, dimension);

        decltype(_quadrature_weights)::HostMirror quadrature_weights = create_mirror_view(_quadrature_weights);
        decltype(_parameterized_quadrature_sites)::HostMirror parameterized_quadrature_sites = create_mirror_view(_parameterized_quadrature_sites);

        if (dimension == 1) {

            switch (_number_of_quadrature_points) {
            case 1:
                quadrature_weights(0) = 1.0;
                parameterized_quadrature_sites(0,0) = 0.5;
                break;
            case 2:
                quadrature_weights(0) = 0.5;
                quadrature_weights(1) = 0.5;
                parameterized_quadrature_sites(0,0) = 0.5*(1-std::sqrt(1./3.));
                parameterized_quadrature_sites(1,0) = 0.5*(1+std::sqrt(1./3.));
                break;
            case 3:
                quadrature_weights(0) = 2.5/9.;
                quadrature_weights(1) = 2.5/9;
                quadrature_weights(2) = 4./9.;
                parameterized_quadrature_sites(0,0) = 0.5*(1-std::sqrt(3./5.));
                parameterized_quadrature_sites(1,0) = 0.5*(1+std::sqrt(3./5.));
                parameterized_quadrature_sites(2,0) = 0.5;
                break;
            case 4:
                quadrature_weights(0) = (18+std::sqrt(30))/72.;
                quadrature_weights(1) = (18+std::sqrt(30))/72.;
                quadrature_weights(2) = (18-std::sqrt(30))/72.;
                quadrature_weights(3) = (18-std::sqrt(30))/72.;
                parameterized_quadrature_sites(0,0) = 0.5*(1-std::sqrt(3./7.-2./7.*std::sqrt(6./5.)));
                parameterized_quadrature_sites(1,0) = 0.5*(1+std::sqrt(3./7.-2./7.*std::sqrt(6./5.)));
                parameterized_quadrature_sites(2,0) = 0.5*(1-std::sqrt(3./7.+2./7.*std::sqrt(6./5.)));
                parameterized_quadrature_sites(3,0) = 0.5*(1+std::sqrt(3./7.+2./7.*std::sqrt(6./5.)));
                break;
            case 5:
                quadrature_weights(0) = 128./450.;
                quadrature_weights(1) = (322+13*std::sqrt(70))/1800.;
                quadrature_weights(2) = (322+13*std::sqrt(70))/1800.;
                quadrature_weights(3) = (322-13*std::sqrt(70))/1800.;
                quadrature_weights(4) = (322-13*std::sqrt(70))/1800.;
                parameterized_quadrature_sites(0,0) = 0.5;
                parameterized_quadrature_sites(1,0) = 0.5*(1+1./3.*std::sqrt(5-2*std::sqrt(10./7.)));
                parameterized_quadrature_sites(2,0) = 0.5*(1-1./3.*std::sqrt(5-2*std::sqrt(10./7.)));
                parameterized_quadrature_sites(3,0) = 0.5*(1+1./3.*std::sqrt(5+2*std::sqrt(10./7.)));
                parameterized_quadrature_sites(4,0) = 0.5*(1-1./3.*std::sqrt(5+2*std::sqrt(10./7.)));
                break;
            default:
                compadre_assert_release(false && "Unsupported number of quadrature points.");
            }

        } else if (dimension == 2) {
            switch (_number_of_quadrature_points) {
            case 1:
                quadrature_weights(0) = 0.5;
                parameterized_quadrature_sites(0,0) = 1./3.;
                parameterized_quadrature_sites(0,1) = 1./3.;
                break;
            case 3:
                quadrature_weights(0) = 1./6.;
                quadrature_weights(1) = 1./6.;
                quadrature_weights(2) = 1./6.;
                //{{1./6., 1./6., 0.0},
                parameterized_quadrature_sites(0,0) = 1./6.;
                parameterized_quadrature_sites(0,1) = 1./6.;
                //{1./6., 2./3., 0.0},
                parameterized_quadrature_sites(1,0) = 1./6.;
                parameterized_quadrature_sites(1,1) = 2./3.;
                //{2./3., 1./6., 0.0}},
                parameterized_quadrature_sites(2,0) = 2./3.;
                parameterized_quadrature_sites(2,1) = 1./6.;
                break;
            case 4:
                quadrature_weights(0) = -9./32.;
                quadrature_weights(1) = 25./96.;
                quadrature_weights(2) = 25./96.;
                quadrature_weights(3) = 25./96.;

                //{1./3., 1./3., 0.0},
                parameterized_quadrature_sites(0,0) = 1./3.;
                parameterized_quadrature_sites(0,1) = 1./3.;
                //{1./5., 1./5., 0.0},
                parameterized_quadrature_sites(1,0) = 1./5.;
                parameterized_quadrature_sites(1,1) = 1./5.;
                //{1./5., 3./5., 0.0},
                parameterized_quadrature_sites(2,0) = 1./5.;
                parameterized_quadrature_sites(2,1) = 3./5.;
                //{3./5., 1./5., 0.0}},
                parameterized_quadrature_sites(3,0) = 3./5.;
                parameterized_quadrature_sites(3,1) = 1./5.;
                break;
            case 6:
                quadrature_weights(0) = 1.1169079483900573053953942803038657817028e-1;
                quadrature_weights(1) = 1.1169079483900573053953942803038657817028e-1;
                quadrature_weights(2) = 1.1169079483900573053953942803038657817028e-1;
                quadrature_weights(3) = 5.4975871827660936395761590801109078322738e-2;
                quadrature_weights(4) = 5.4975871827660936395761590801109078322738e-2;
                quadrature_weights(5) = 5.4975871827660936395761590801109078322738e-2;
                //{4.4594849091596487577332043252695176084800e-1,  4.4594849091596487577332043252695176084800e-1, 0.0},
                parameterized_quadrature_sites(0,0) = 4.4594849091596487577332043252695176084800e-1;
                parameterized_quadrature_sites(0,1) = 4.4594849091596487577332043252695176084800e-1;
                //{4.4594849091596487577332043252695176084800e-1,  1.0810301816807024845335913494609647830400e-1, 0.0},
                parameterized_quadrature_sites(1,0) = 4.4594849091596487577332043252695176084800e-1;
                parameterized_quadrature_sites(1,1) = 1.0810301816807024845335913494609647830400e-1;
                //{1.0810301816807024845335913494609647830400e-1,  4.4594849091596487577332043252695176084800e-1, 0.0},
                parameterized_quadrature_sites(2,0) = 1.0810301816807024845335913494609647830400e-1;
                parameterized_quadrature_sites(2,1) = 4.4594849091596487577332043252695176084800e-1;
                //{9.1576213509770745141706146463816673424873e-2,  9.1576213509770745141706146463816673424873e-2, 0.0},
                parameterized_quadrature_sites(3,0) = 9.1576213509770745141706146463816673424873e-2;
                parameterized_quadrature_sites(3,1) = 9.1576213509770745141706146463816673424873e-2;
                //{9.1576213509770745141706146463816673424873e-2,  8.1684757298045850971658770707236665315025e-1, 0.0},
                parameterized_quadrature_sites(4,0) = 9.1576213509770745141706146463816673424873e-2;
                parameterized_quadrature_sites(4,1) = 8.1684757298045850971658770707236665315025e-1;
                //{8.1684757298045850971658770707236665315025e-1,  9.1576213509770745141706146463816673424873e-2, 0.0}},
                parameterized_quadrature_sites(5,0) = 8.1684757298045850971658770707236665315025e-1;
                parameterized_quadrature_sites(5,1) = 9.1576213509770745141706146463816673424873e-2;
                break;
			case 7:
			    parameterized_quadrature_sites(0,0) = 3.3333333333333333333333333333333333333333e-1;
			    parameterized_quadrature_sites(0,1) = 3.3333333333333333333333333333333333333333e-1;
			    parameterized_quadrature_sites(1,0) = 1.0128650732345633880098736191512382805558e-1;
			    parameterized_quadrature_sites(1,1) = 1.0128650732345633880098736191512382805558e-1;
			    parameterized_quadrature_sites(2,0) = 7.9742698535308732239802527616975234388885e-1;
			    parameterized_quadrature_sites(2,1) = 1.0128650732345633880098736191512382805558e-1;
			    parameterized_quadrature_sites(3,0) = 1.0128650732345633880098736191512382805558e-1;
			    parameterized_quadrature_sites(3,1) = 7.9742698535308732239802527616975234388885e-1;
			    parameterized_quadrature_sites(4,0) = 4.7014206410511508977044120951344760051585e-1;
			    parameterized_quadrature_sites(4,1) = 4.7014206410511508977044120951344760051585e-1;
			    parameterized_quadrature_sites(5,0) = 5.9715871789769820459117580973104798968293e-2;
			    parameterized_quadrature_sites(5,1) = 4.7014206410511508977044120951344760051585e-1;
			    parameterized_quadrature_sites(6,0) = 4.7014206410511508977044120951344760051585e-1;
			    parameterized_quadrature_sites(6,1) = 5.9715871789769820459117580973104798968293e-2;
			    quadrature_weights(0) = 1.1250000000000000000000000000000000000000e-1;
			    quadrature_weights(1) = 6.2969590272413576297841972750090666828820e-2;
			    quadrature_weights(2) = 6.2969590272413576297841972750090666828820e-2;
			    quadrature_weights(3) = 6.2969590272413576297841972750090666828820e-2;
			    quadrature_weights(4) = 6.6197076394253090368824693916575999837847e-2;
			    quadrature_weights(5) = 6.6197076394253090368824693916575999837847e-2;
			    quadrature_weights(6) = 6.6197076394253090368824693916575999837847e-2;
			    break;
			case 12:
			    parameterized_quadrature_sites(0,0) = 6.3089014491502228340331602870819157341003e-2;
			    parameterized_quadrature_sites(0,1) = 6.3089014491502228340331602870819157341003e-2;
			    parameterized_quadrature_sites(1,0) = 6.3089014491502228340331602870819157341003e-2;
			    parameterized_quadrature_sites(1,1) = 8.7382197101699554331933679425836168531799e-1;
			    parameterized_quadrature_sites(2,0) = 8.7382197101699554331933679425836168531799e-1;
			    parameterized_quadrature_sites(2,1) = 6.3089014491502228340331602870819157341003e-2;
			    parameterized_quadrature_sites(3,0) = 2.4928674517091042129163855310701907608796e-1;
			    parameterized_quadrature_sites(3,1) = 2.4928674517091042129163855310701907608796e-1;
			    parameterized_quadrature_sites(4,0) = 2.4928674517091042129163855310701907608796e-1;
			    parameterized_quadrature_sites(4,1) = 5.0142650965817915741672289378596184782407e-1;
			    parameterized_quadrature_sites(5,0) = 5.0142650965817915741672289378596184782407e-1;
			    parameterized_quadrature_sites(5,1) = 2.4928674517091042129163855310701907608796e-1;
			    parameterized_quadrature_sites(6,0) = 3.1035245103378441286759566723265117140082e-1;
			    parameterized_quadrature_sites(6,1) = 5.3145049844816939902261738355299128796183e-2;
			    parameterized_quadrature_sites(7,0) = 6.3650249912139865667831578657006182134850e-1;
			    parameterized_quadrature_sites(7,1) = 3.1035245103378439596843454179854003165774e-1;
			    parameterized_quadrature_sites(8,0) = 5.3145049844816930454088546197287007250682e-2;
			    parameterized_quadrature_sites(8,1) = 6.3650249912139866412930371984616083954608e-1;
			    parameterized_quadrature_sites(9,0) = 5.3145049844816939902261738355299128796183e-2;
			    parameterized_quadrature_sites(9,1) = 3.1035245103378441286759566723265117140082e-1;
			    parameterized_quadrature_sites(10,0) = 6.3650249912139866412930371984616083954608e-1;
			    parameterized_quadrature_sites(10,1) = 5.3145049844816930454088546197287007250682e-2;
			    parameterized_quadrature_sites(11,0) = 3.1035245103378439596843454179854003165774e-1;
			    parameterized_quadrature_sites(11,1) = 6.3650249912139865667831578657006182134850e-1;
			    quadrature_weights(0) = 2.5422453185103408460468404553434492023395e-2;
			    quadrature_weights(1) = 2.5422453185103408460468404553434492023395e-2;
			    quadrature_weights(2) = 2.5422453185103408460468404553434492023395e-2;
			    quadrature_weights(3) = 5.8393137863189683012644805692789720663043e-2;
			    quadrature_weights(4) = 5.8393137863189683012644805692789720663043e-2;
			    quadrature_weights(5) = 5.8393137863189683012644805692789720663043e-2;
			    quadrature_weights(6) = 4.1425537809186787596776728210221226990114e-2;
			    quadrature_weights(7) = 4.1425537809186787596776728210221226990114e-2;
			    quadrature_weights(8) = 4.1425537809186787596776728210221226990114e-2;
			    quadrature_weights(9) = 4.1425537809186787596776728210221226990114e-2;
			    quadrature_weights(10) = 4.1425537809186787596776728210221226990114e-2;
			    quadrature_weights(11) = 4.1425537809186787596776728210221226990114e-2;
			    break;
			case 13:
			    parameterized_quadrature_sites(0,0) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(0,1) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(1,0) = 2.60345966079040e-01;
			    parameterized_quadrature_sites(1,1) = 2.60345966079040e-01;
			    parameterized_quadrature_sites(2,0) = 2.60345966079040e-01;
			    parameterized_quadrature_sites(2,1) = 4.79308067841920e-01;
			    parameterized_quadrature_sites(3,0) = 4.79308067841920e-01;
			    parameterized_quadrature_sites(3,1) = 2.60345966079040e-01;
			    parameterized_quadrature_sites(4,0) = 6.51301029022160e-02;
			    parameterized_quadrature_sites(4,1) = 6.51301029022160e-02;
			    parameterized_quadrature_sites(5,0) = 6.51301029022160e-02;
			    parameterized_quadrature_sites(5,1) = 8.69739794195568e-01;
			    parameterized_quadrature_sites(6,0) = 8.69739794195568e-01;
			    parameterized_quadrature_sites(6,1) = 6.51301029022160e-02;
			    parameterized_quadrature_sites(7,0) = 3.12865496004874e-01;
			    parameterized_quadrature_sites(7,1) = 6.38444188569810e-01;
			    parameterized_quadrature_sites(8,0) = 6.38444188569810e-01;
			    parameterized_quadrature_sites(8,1) = 4.86903154253160e-02;
			    parameterized_quadrature_sites(9,0) = 4.86903154253160e-02;
			    parameterized_quadrature_sites(9,1) = 3.12865496004874e-01;
			    parameterized_quadrature_sites(10,0) = 3.12865496004874e-01;
			    parameterized_quadrature_sites(10,1) = 4.86903154253160e-02;
			    parameterized_quadrature_sites(11,0) = 6.38444188569810e-01;
			    parameterized_quadrature_sites(11,1) = 3.12865496004874e-01;
			    parameterized_quadrature_sites(12,0) = 4.86903154253160e-02;
			    parameterized_quadrature_sites(12,1) = 6.38444188569810e-01;
			    quadrature_weights(0) = -7.47850222338410e-02;
			    quadrature_weights(1) = 8.78076287166040e-02;
			    quadrature_weights(2) = 8.78076287166040e-02;
			    quadrature_weights(3) = 8.78076287166040e-02;
			    quadrature_weights(4) = 2.66736178044190e-02;
			    quadrature_weights(5) = 2.66736178044190e-02;
			    quadrature_weights(6) = 2.66736178044190e-02;
			    quadrature_weights(7) = 3.85568804451285e-02;
			    quadrature_weights(8) = 3.85568804451285e-02;
			    quadrature_weights(9) = 3.85568804451285e-02;
			    quadrature_weights(10) = 3.85568804451285e-02;
			    quadrature_weights(11) = 3.85568804451285e-02;
			    quadrature_weights(12) = 3.85568804451285e-02;
			    break;
			case 16:
			    parameterized_quadrature_sites(0,0) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(0,1) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(1,0) = 4.59292588292723e-01;
			    parameterized_quadrature_sites(1,1) = 4.59292588292723e-01;
			    parameterized_quadrature_sites(2,0) = 4.59292588292723e-01;
			    parameterized_quadrature_sites(2,1) = 8.14148234145540e-02;
			    parameterized_quadrature_sites(3,0) = 8.14148234145540e-02;
			    parameterized_quadrature_sites(3,1) = 4.59292588292723e-01;
			    parameterized_quadrature_sites(4,0) = 1.70569307751760e-01;
			    parameterized_quadrature_sites(4,1) = 1.70569307751760e-01;
			    parameterized_quadrature_sites(5,0) = 1.70569307751760e-01;
			    parameterized_quadrature_sites(5,1) = 6.58861384496480e-01;
			    parameterized_quadrature_sites(6,0) = 6.58861384496480e-01;
			    parameterized_quadrature_sites(6,1) = 1.70569307751760e-01;
			    parameterized_quadrature_sites(7,0) = 5.05472283170310e-02;
			    parameterized_quadrature_sites(7,1) = 5.05472283170310e-02;
			    parameterized_quadrature_sites(8,0) = 5.05472283170310e-02;
			    parameterized_quadrature_sites(8,1) = 8.98905543365938e-01;
			    parameterized_quadrature_sites(9,0) = 8.98905543365938e-01;
			    parameterized_quadrature_sites(9,1) = 5.05472283170310e-02;
			    parameterized_quadrature_sites(10,0) = 2.63112829634638e-01;
			    parameterized_quadrature_sites(10,1) = 7.28492392955404e-01;
			    parameterized_quadrature_sites(11,0) = 7.28492392955404e-01;
			    parameterized_quadrature_sites(11,1) = 8.39477740995798e-03;
			    parameterized_quadrature_sites(12,0) = 8.39477740995798e-03;
			    parameterized_quadrature_sites(12,1) = 2.63112829634638e-01;
			    parameterized_quadrature_sites(13,0) = 2.63112829634638e-01;
			    parameterized_quadrature_sites(13,1) = 8.39477740995798e-03;
			    parameterized_quadrature_sites(14,0) = 7.28492392955404e-01;
			    parameterized_quadrature_sites(14,1) = 2.63112829634638e-01;
			    parameterized_quadrature_sites(15,0) = 8.39477740995798e-03;
			    parameterized_quadrature_sites(15,1) = 7.28492392955404e-01;
			    quadrature_weights(0) = 7.21578038388935e-02;
			    quadrature_weights(1) = 4.75458171336425e-02;
			    quadrature_weights(2) = 4.75458171336425e-02;
			    quadrature_weights(3) = 4.75458171336425e-02;
			    quadrature_weights(4) = 5.16086852673590e-02;
			    quadrature_weights(5) = 5.16086852673590e-02;
			    quadrature_weights(6) = 5.16086852673590e-02;
			    quadrature_weights(7) = 1.62292488115990e-02;
			    quadrature_weights(8) = 1.62292488115988e-02;
			    quadrature_weights(9) = 1.62292488115990e-02;
			    quadrature_weights(10) = 1.36151570872175e-02;
			    quadrature_weights(11) = 1.36151570872175e-02;
			    quadrature_weights(12) = 1.36151570872175e-02;
			    quadrature_weights(13) = 1.36151570872175e-02;
			    quadrature_weights(14) = 1.36151570872175e-02;
			    quadrature_weights(15) = 1.36151570872175e-02;
			    break;
			case 19:
			    parameterized_quadrature_sites(0,0) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(0,1) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(1,0) = 4.89682519198738e-01;
			    parameterized_quadrature_sites(1,1) = 4.89682519198738e-01;
			    parameterized_quadrature_sites(2,0) = 4.89682519198738e-01;
			    parameterized_quadrature_sites(2,1) = 2.06349616025250e-02;
			    parameterized_quadrature_sites(3,0) = 2.06349616025250e-02;
			    parameterized_quadrature_sites(3,1) = 4.89682519198738e-01;
			    parameterized_quadrature_sites(4,0) = 4.37089591492937e-01;
			    parameterized_quadrature_sites(4,1) = 4.37089591492937e-01;
			    parameterized_quadrature_sites(5,0) = 4.37089591492937e-01;
			    parameterized_quadrature_sites(5,1) = 1.25820817014127e-01;
			    parameterized_quadrature_sites(6,0) = 1.25820817014127e-01;
			    parameterized_quadrature_sites(6,1) = 4.37089591492937e-01;
			    parameterized_quadrature_sites(7,0) = 1.88203535619033e-01;
			    parameterized_quadrature_sites(7,1) = 1.88203535619033e-01;
			    parameterized_quadrature_sites(8,0) = 1.88203535619033e-01;
			    parameterized_quadrature_sites(8,1) = 6.23592928761935e-01;
			    parameterized_quadrature_sites(9,0) = 6.23592928761935e-01;
			    parameterized_quadrature_sites(9,1) = 1.88203535619033e-01;
			    parameterized_quadrature_sites(10,0) = 4.47295133944530e-02;
			    parameterized_quadrature_sites(10,1) = 4.47295133944530e-02;
			    parameterized_quadrature_sites(11,0) = 4.47295133944530e-02;
			    parameterized_quadrature_sites(11,1) = 9.10540973211095e-01;
			    parameterized_quadrature_sites(12,0) = 9.10540973211095e-01;
			    parameterized_quadrature_sites(12,1) = 4.47295133944530e-02;
			    parameterized_quadrature_sites(13,0) = 2.21962989160766e-01;
			    parameterized_quadrature_sites(13,1) = 7.41198598784498e-01;
			    parameterized_quadrature_sites(14,0) = 7.41198598784498e-01;
			    parameterized_quadrature_sites(14,1) = 3.68384120547360e-02;
			    parameterized_quadrature_sites(15,0) = 3.68384120547360e-02;
			    parameterized_quadrature_sites(15,1) = 2.21962989160766e-01;
			    parameterized_quadrature_sites(16,0) = 2.21962989160766e-01;
			    parameterized_quadrature_sites(16,1) = 3.68384120547360e-02;
			    parameterized_quadrature_sites(17,0) = 7.41198598784498e-01;
			    parameterized_quadrature_sites(17,1) = 2.21962989160766e-01;
			    parameterized_quadrature_sites(18,0) = 3.68384120547360e-02;
			    parameterized_quadrature_sites(18,1) = 7.41198598784498e-01;
			    quadrature_weights(0) = 4.85678981413995e-02;
			    quadrature_weights(1) = 1.56673501135695e-02;
			    quadrature_weights(2) = 1.56673501135695e-02;
			    quadrature_weights(3) = 1.56673501135695e-02;
			    quadrature_weights(4) = 3.89137705023870e-02;
			    quadrature_weights(5) = 3.89137705023870e-02;
			    quadrature_weights(6) = 3.89137705023870e-02;
			    quadrature_weights(7) = 3.98238694636050e-02;
			    quadrature_weights(8) = 3.98238694636050e-02;
			    quadrature_weights(9) = 3.98238694636050e-02;
			    quadrature_weights(10) = 1.27888378293490e-02;
			    quadrature_weights(11) = 1.27888378293490e-02;
			    quadrature_weights(12) = 1.27888378293490e-02;
			    quadrature_weights(13) = 2.16417696886445e-02;
			    quadrature_weights(14) = 2.16417696886445e-02;
			    quadrature_weights(15) = 2.16417696886445e-02;
			    quadrature_weights(16) = 2.16417696886445e-02;
			    quadrature_weights(17) = 2.16417696886445e-02;
			    quadrature_weights(18) = 2.16417696886445e-02;
			    break;
			case 24:
			    parameterized_quadrature_sites(0,0) = 3.8102570854643002e-03;
			    parameterized_quadrature_sites(0,1) = 8.6854386943076545e-01;
			    parameterized_quadrature_sites(1,0) = 8.3865349500109043e-01;
			    parameterized_quadrature_sites(1,1) = 1.6134650499890960e-01;
			    parameterized_quadrature_sites(2,0) = 0.0e-00;
			    parameterized_quadrature_sites(2,1) = 3.9366774470722010e-01;
			    parameterized_quadrature_sites(3,0) = 7.7757518429429107e-01;
			    parameterized_quadrature_sites(3,1) = 0.0e-00;
			    parameterized_quadrature_sites(4,0) = 4.7768381772022403e-02;
			    parameterized_quadrature_sites(4,1) = 9.2899486985787905e-01;
			    parameterized_quadrature_sites(5,0) = 1.0939142057119900e-02;
			    parameterized_quadrature_sites(5,1) = 1.7690730625559031e-01;
			    parameterized_quadrature_sites(6,0) = 4.6374383867430541e-01;
			    parameterized_quadrature_sites(6,1) = 0.0e-00;
			    parameterized_quadrature_sites(7,0) = 9.3049846900263089e-01;
			    parameterized_quadrature_sites(7,1) = 2.9553592846822900e-02;
			    parameterized_quadrature_sites(8,0) = 3.9099745550423302e-02;
			    parameterized_quadrature_sites(8,1) = 3.5319656252586103e-02;
			    parameterized_quadrature_sites(9,0) = 4.8798437805397499e-01;
			    parameterized_quadrature_sites(9,1) = 5.0365825075943971e-01;
			    parameterized_quadrature_sites(10,0) = 1.9305903224251941e-01;
			    parameterized_quadrature_sites(10,1) = 3.0573404093099301e-02;
			    parameterized_quadrature_sites(11,0) = 2.2376358774275851e-01;
			    parameterized_quadrature_sites(11,1) = 7.4726591728868819e-01;
			    parameterized_quadrature_sites(12,0) = 3.6036266787907702e-02;
			    parameterized_quadrature_sites(12,1) = 6.3491832379200652e-01;
			    parameterized_quadrature_sites(13,0) = 7.6777680170023954e-01;
			    parameterized_quadrature_sites(13,1) = 1.0614642990290001e-01;
			    parameterized_quadrature_sites(14,0) = 1.0954959855585469e-01;
			    parameterized_quadrature_sites(14,1) = 7.5329402776254240e-01;
			    parameterized_quadrature_sites(15,0) = 6.4203365318662664e-01;
			    parameterized_quadrature_sites(15,1) = 2.9530445535851102e-01;
			    parameterized_quadrature_sites(16,0) = 1.0999439055630450e-01;
			    parameterized_quadrature_sites(16,1) = 1.6929927488966459e-01;
			    parameterized_quadrature_sites(17,0) = 3.3947290311800549e-01;
			    parameterized_quadrature_sites(17,1) = 9.5379208487721703e-02;
			    parameterized_quadrature_sites(18,0) = 8.4198522115543697e-02;
			    parameterized_quadrature_sites(18,1) = 3.8729657913960353e-01;
			    parameterized_quadrature_sites(19,0) = 5.7966325105486349e-01;
			    parameterized_quadrature_sites(19,1) = 8.0491894656105595e-02;
			    parameterized_quadrature_sites(20,0) = 3.6419744430339263e-01;
			    parameterized_quadrature_sites(20,1) = 5.2433682558924433e-01;
			    parameterized_quadrature_sites(21,0) = 2.7586334089315973e-01;
			    parameterized_quadrature_sites(21,1) = 2.6481531651496770e-01;
			    parameterized_quadrature_sites(22,0) = 2.0776116575484829e-01;
			    parameterized_quadrature_sites(22,1) = 5.0550507373529086e-01;
			    parameterized_quadrature_sites(23,0) = 4.8123289062464247e-01;
			    parameterized_quadrature_sites(23,1) = 2.7542385024412980e-01;
			    quadrature_weights(0) = 5.3333456379563004e-03;
			    quadrature_weights(1) = 5.5334733911235499e-03;
			    quadrature_weights(2) = 5.9501524071177998e-03;
			    quadrature_weights(3) = 6.6995403929284002e-03;
			    quadrature_weights(4) = 7.3041646115608498e-03;
			    quadrature_weights(5) = 7.4233118233402503e-03;
			    quadrature_weights(6) = 7.5781309587828498e-03;
			    quadrature_weights(7) = 7.6015470130063496e-03;
			    quadrature_weights(8) = 8.8279222964811003e-03;
			    quadrature_weights(9) = 9.7797061087609508e-03;
			    quadrature_weights(10) = 1.5707351180334426e-02;
			    quadrature_weights(11) = 1.6151301011728625e-02;
			    quadrature_weights(12) = 2.1481164696039650e-02;
			    quadrature_weights(13) = 2.6268329551205824e-02;
			    quadrature_weights(14) = 2.6365190704419499e-02;
			    quadrature_weights(15) = 2.7457802196927176e-02;
			    quadrature_weights(16) = 2.8987070299349601e-02;
			    quadrature_weights(17) = 3.0871122933097676e-02;
			    quadrature_weights(18) = 3.1343962650456775e-02;
			    quadrature_weights(19) = 3.3093067237629975e-02;
			    quadrature_weights(20) = 3.7039823668389947e-02;
			    quadrature_weights(21) = 4.2207221279855073e-02;
			    quadrature_weights(22) = 4.3362019313832399e-02;
			    quadrature_weights(23) = 4.7633278635674972e-02;
			    break;
			case 27:
			    parameterized_quadrature_sites(0,0) = 3.7802163891336921e-01;
			    parameterized_quadrature_sites(0,1) = 6.1948431533135195e-01;
			    parameterized_quadrature_sites(1,0) = 3.2899822292186298e-02;
			    parameterized_quadrature_sites(1,1) = 9.3614893514675623e-01;
			    parameterized_quadrature_sites(2,0) = 9.3551434285897095e-01;
			    parameterized_quadrature_sites(2,1) = 3.3268560622678398e-02;
			    parameterized_quadrature_sites(3,0) = 3.4222771841359197e-02;
			    parameterized_quadrature_sites(3,1) = 3.2916403878999703e-02;
			    parameterized_quadrature_sites(4,0) = 1.4354532010930900e-02;
			    parameterized_quadrature_sites(4,1) = 3.9659731669586501e-01;
			    parameterized_quadrature_sites(5,0) = 2.2120535196161799e-02;
			    parameterized_quadrature_sites(5,1) = 1.6892970982290231e-01;
			    parameterized_quadrature_sites(6,0) = 8.1562969693268217e-01;
			    parameterized_quadrature_sites(6,1) = 2.6807150626772601e-02;
			    parameterized_quadrature_sites(7,0) = 2.7719522918618601e-02;
			    parameterized_quadrature_sites(7,1) = 8.1626233715968810e-01;
			    parameterized_quadrature_sites(8,0) = 1.7400571673032261e-01;
			    parameterized_quadrature_sites(8,1) = 2.5252704638304500e-02;
			    parameterized_quadrature_sites(9,0) = 3.8913981113319362e-01;
			    parameterized_quadrature_sites(9,1) = 2.2592651051306600e-02;
			    parameterized_quadrature_sites(10,0) = 8.0364834053903877e-01;
			    parameterized_quadrature_sites(10,1) = 1.6655614492060569e-01;
			    parameterized_quadrature_sites(11,0) = 1.6429286715713459e-01;
			    parameterized_quadrature_sites(11,1) = 8.0454974747615537e-01;
			    parameterized_quadrature_sites(12,0) = 6.1758873171277151e-01;
			    parameterized_quadrature_sites(12,1) = 2.5660186833052399e-02;
			    parameterized_quadrature_sites(13,0) = 2.6297199713764201e-02;
			    parameterized_quadrature_sites(13,1) = 6.1924873232110123e-01;
			    parameterized_quadrature_sites(14,0) = 5.9895439629934211e-01;
			    parameterized_quadrature_sites(14,1) = 3.7272769861629101e-01;
			    parameterized_quadrature_sites(15,0) = 8.1721404855381805e-02;
			    parameterized_quadrature_sites(15,1) = 3.2719878157552901e-01;
			    parameterized_quadrature_sites(16,0) = 1.3035453031942690e-01;
			    parameterized_quadrature_sites(16,1) = 1.3667083534390509e-01;
			    parameterized_quadrature_sites(17,0) = 7.1027868107761583e-01;
			    parameterized_quadrature_sites(17,1) = 1.3828000204292321e-01;
			    parameterized_quadrature_sites(18,0) = 1.4118119730952799e-01;
			    parameterized_quadrature_sites(18,1) = 7.0099267949645228e-01;
			    parameterized_quadrature_sites(19,0) = 5.3141960154079959e-01;
			    parameterized_quadrature_sites(19,1) = 1.2417148586801489e-01;
			    parameterized_quadrature_sites(20,0) = 3.4992914334288650e-01;
			    parameterized_quadrature_sites(20,1) = 5.6938486195327997e-01;
			    parameterized_quadrature_sites(21,0) = 3.1909737814681871e-01;
			    parameterized_quadrature_sites(21,1) = 1.1698976413323441e-01;
			    parameterized_quadrature_sites(22,0) = 1.2454405910544100e-01;
			    parameterized_quadrature_sites(22,1) = 5.1353143433447235e-01;
			    parameterized_quadrature_sites(23,0) = 4.1132499178904658e-01;
			    parameterized_quadrature_sites(23,1) = 2.6677168071577739e-01;
			    parameterized_quadrature_sites(24,0) = 5.3634228112084714e-01;
			    parameterized_quadrature_sites(24,1) = 3.2081957909482989e-01;
			    parameterized_quadrature_sites(25,0) = 2.2789955884347499e-01;
			    parameterized_quadrature_sites(25,1) = 2.8790310224819649e-01;
			    parameterized_quadrature_sites(26,0) = 2.9133859436942361e-01;
			    parameterized_quadrature_sites(26,1) = 4.6494564773693992e-01;
			    quadrature_weights(0) = 5.6375356078552001e-03;
			    quadrature_weights(1) = 6.5529674260441251e-03;
			    quadrature_weights(2) = 6.6871546432313500e-03;
			    quadrature_weights(3) = 7.2740543402810501e-03;
			    quadrature_weights(4) = 9.1794421950682249e-03;
			    quadrature_weights(5) = 1.0688904048956826e-02;
			    quadrature_weights(6) = 1.1556190176810576e-02;
			    quadrature_weights(7) = 1.1821029782882325e-02;
			    quadrature_weights(8) = 1.2035150250304100e-02;
			    quadrature_weights(9) = 1.3361393587366349e-02;
			    quadrature_weights(10) = 1.3378318081626600e-02;
			    quadrature_weights(11) = 1.3797950075089925e-02;
			    quadrature_weights(12) = 1.4360888683763475e-02;
			    quadrature_weights(13) = 1.4717279053216750e-02;
			    quadrature_weights(14) = 1.6298436572453948e-02;
			    quadrature_weights(15) = 2.1298852449057524e-02;
			    quadrature_weights(16) = 2.2769129560075176e-02;
			    quadrature_weights(17) = 2.4595497454187175e-02;
			    quadrature_weights(18) = 2.6030267669060574e-02;
			    quadrature_weights(19) = 2.7060738237627526e-02;
			    quadrature_weights(20) = 2.7485999004421752e-02;
			    quadrature_weights(21) = 2.7495506358711074e-02;
			    quadrature_weights(22) = 2.8322942561349675e-02;
			    quadrature_weights(23) = 2.9305884203704700e-02;
			    quadrature_weights(24) = 3.1096575055627549e-02;
			    quadrature_weights(25) = 3.3071222942480799e-02;
			    quadrature_weights(26) = 3.4120689978745525e-02;
			    break;
			case 32:
			    parameterized_quadrature_sites(0,0) = 9.2734897448394982e-01;
			    parameterized_quadrature_sites(0,1) = 0.0e-00;
			    parameterized_quadrature_sites(1,0) = 2.3551733249578700e-02;
			    parameterized_quadrature_sites(1,1) = 9.5526919357006035e-01;
			    parameterized_quadrature_sites(2,0) = 0.0e-00;
			    parameterized_quadrature_sites(2,1) = 8.5815888421533082e-01;
			    parameterized_quadrature_sites(3,0) = 9.4547507322097091e-01;
			    parameterized_quadrature_sites(3,1) = 4.3010560106405499e-02;
			    parameterized_quadrature_sites(4,0) = 1.5406460162685609e-01;
			    parameterized_quadrature_sites(4,1) = 8.4593539837314391e-01;
			    parameterized_quadrature_sites(5,0) = 0.0e-00;
			    parameterized_quadrature_sites(5,1) = 6.2731531923241179e-01;
			    parameterized_quadrature_sites(6,0) = 2.7110971356255800e-02;
			    parameterized_quadrature_sites(6,1) = 2.9754117496841800e-02;
			    parameterized_quadrature_sites(7,0) = 1.4604496167217570e-01;
			    parameterized_quadrature_sites(7,1) = 9.2296909059649008e-03;
			    parameterized_quadrature_sites(8,0) = 2.1152223383121900e-02;
			    parameterized_quadrature_sites(8,1) = 1.5557066896897950e-01;
			    parameterized_quadrature_sites(9,0) = 1.4566514788347000e-02;
			    parameterized_quadrature_sites(9,1) = 3.6384660446077510e-01;
			    parameterized_quadrature_sites(10,0) = 7.8860171922313160e-01;
			    parameterized_quadrature_sites(10,1) = 1.8920633061715941e-01;
			    parameterized_quadrature_sites(11,0) = 7.4918973979067949e-01;
			    parameterized_quadrature_sites(11,1) = 2.3088148766115799e-02;
			    parameterized_quadrature_sites(12,0) = 7.1871496101589105e-02;
			    parameterized_quadrature_sites(12,1) = 8.5431474947580432e-01;
			    parameterized_quadrature_sites(13,0) = 3.3212908394764512e-01;
			    parameterized_quadrature_sites(13,1) = 2.4506286636990001e-02;
			    parameterized_quadrature_sites(14,0) = 3.6118159118967208e-01;
			    parameterized_quadrature_sites(14,1) = 6.1600929617267497e-01;
			    parameterized_quadrature_sites(15,0) = 2.4345813394879970e-01;
			    parameterized_quadrature_sites(15,1) = 9.3448087604440996e-02;
			    parameterized_quadrature_sites(16,0) = 5.8168921474014745e-01;
			    parameterized_quadrature_sites(16,1) = 3.9316510319604808e-01;
			    parameterized_quadrature_sites(17,0) = 5.4444667627192522e-01;
			    parameterized_quadrature_sites(17,1) = 2.5716283623693902e-02;
			    parameterized_quadrature_sites(18,0) = 8.2600331401756000e-01;
			    parameterized_quadrature_sites(18,1) = 7.9955384841381302e-02;
			    parameterized_quadrature_sites(19,0) = 1.1638649906727730e-01;
			    parameterized_quadrature_sites(19,1) = 8.9602705800587407e-02;
			    parameterized_quadrature_sites(20,0) = 2.0376848107772980e-01;
			    parameterized_quadrature_sites(20,1) = 7.1788185898052326e-01;
			    parameterized_quadrature_sites(21,0) = 6.4413220382260494e-02;
			    parameterized_quadrature_sites(21,1) = 7.1008125956836521e-01;
			    parameterized_quadrature_sites(22,0) = 9.5428585810584596e-02;
			    parameterized_quadrature_sites(22,1) = 2.6077068256562902e-01;
			    parameterized_quadrature_sites(23,0) = 2.4498296509349021e-01;
			    parameterized_quadrature_sites(23,1) = 2.1117939909804931e-01;
			    parameterized_quadrature_sites(24,0) = 7.0566724344036796e-02;
			    parameterized_quadrature_sites(24,1) = 4.9732063377796598e-01;
			    parameterized_quadrature_sites(25,0) = 6.1938125736255578e-01;
			    parameterized_quadrature_sites(25,1) = 1.2512299505810390e-01;
			    parameterized_quadrature_sites(26,0) = 6.2768261568031403e-01;
			    parameterized_quadrature_sites(26,1) = 2.5015500335339208e-01;
			    parameterized_quadrature_sites(27,0) = 4.2260565743346001e-01;
			    parameterized_quadrature_sites(27,1) = 1.2953296900433620e-01;
			    parameterized_quadrature_sites(28,0) = 2.1078525939140391e-01;
			    parameterized_quadrature_sites(28,1) = 3.7986021093401962e-01;
			    parameterized_quadrature_sites(29,0) = 4.0896380449124481e-01;
			    parameterized_quadrature_sites(29,1) = 4.6631787462323071e-01;
			    parameterized_quadrature_sites(30,0) = 2.1377743253005960e-01;
			    parameterized_quadrature_sites(30,1) = 5.5802528953120256e-01;
			    parameterized_quadrature_sites(31,0) = 4.0978657777002531e-01;
			    parameterized_quadrature_sites(31,1) = 3.0141709320909299e-01;
			    quadrature_weights(0) = 2.4440249073300249e-03;
			    quadrature_weights(1) = 3.3379500136836750e-03;
			    quadrature_weights(2) = 3.4227673271718501e-03;
			    quadrature_weights(3) = 3.5598757180403499e-03;
			    quadrature_weights(4) = 3.8572461868124248e-03;
			    quadrature_weights(5) = 4.8273543712181750e-03;
			    quadrature_weights(6) = 5.2546633678012501e-03;
			    quadrature_weights(7) = 5.3404218288141498e-03;
			    quadrature_weights(8) = 9.2418429056154005e-03;
			    quadrature_weights(9) = 9.2727402108032757e-03;
			    quadrature_weights(10) = 1.0310002059841049e-02;
			    quadrature_weights(11) = 1.0842542708505725e-02;
			    quadrature_weights(12) = 1.1245373099579100e-02;
			    quadrature_weights(13) = 1.2452036600753851e-02;
			    quadrature_weights(14) = 1.2549586713842551e-02;
			    quadrature_weights(15) = 1.3971867159939951e-02;
			    quadrature_weights(16) = 1.4072779302606675e-02;
			    quadrature_weights(17) = 1.4084827229864975e-02;
			    quadrature_weights(18) = 1.5264586206036251e-02;
			    quadrature_weights(19) = 1.5287638802019450e-02;
			    quadrature_weights(20) = 1.9786802896485999e-02;
			    quadrature_weights(21) = 2.0640943697731274e-02;
			    quadrature_weights(22) = 2.2968921082895850e-02;
			    quadrature_weights(23) = 2.3749787662653600e-02;
			    quadrature_weights(24) = 2.4074402518453650e-02;
			    quadrature_weights(25) = 2.5482462438393926e-02;
			    quadrature_weights(26) = 2.6676041524410474e-02;
			    quadrature_weights(27) = 2.7073436306583775e-02;
			    quadrature_weights(28) = 2.9718916975567701e-02;
			    quadrature_weights(29) = 2.9994853663553051e-02;
			    quadrature_weights(30) = 3.1582273211328352e-02;
			    quadrature_weights(31) = 3.7611031301662198e-02;
			    break;
			case 36:
			    parameterized_quadrature_sites(0,0) = 2.4293535159026700e-02;
			    parameterized_quadrature_sites(0,1) = 9.4930592938464031e-01;
			    parameterized_quadrature_sites(1,0) = 2.6519342772158901e-02;
			    parameterized_quadrature_sites(1,1) = 2.4269513064041098e-02;
			    parameterized_quadrature_sites(2,0) = 9.4921260235510574e-01;
			    parameterized_quadrature_sites(2,1) = 2.6506796643724899e-02;
			    parameterized_quadrature_sites(3,0) = 3.3775763749036999e-03;
			    parameterized_quadrature_sites(3,1) = 4.7673164123630779e-01;
			    parameterized_quadrature_sites(4,0) = 4.7576722981011582e-01;
			    parameterized_quadrature_sites(4,1) = 5.1989218291019390e-01;
			    parameterized_quadrature_sites(5,0) = 5.1907831934706850e-01;
			    parameterized_quadrature_sites(5,1) = 5.5912706202052003e-03;
			    parameterized_quadrature_sites(6,0) = 8.6168397453205303e-01;
			    parameterized_quadrature_sites(6,1) = 1.3399604861818300e-02;
			    parameterized_quadrature_sites(7,0) = 1.2492097599255590e-01;
			    parameterized_quadrature_sites(7,1) = 8.6130543213341393e-01;
			    parameterized_quadrature_sites(8,0) = 1.3856545386105401e-02;
			    parameterized_quadrature_sites(8,1) = 1.2477337173584679e-01;
			    parameterized_quadrature_sites(9,0) = 2.1188706422168500e-02;
			    parameterized_quadrature_sites(9,1) = 8.4384383512226457e-01;
			    parameterized_quadrature_sites(10,0) = 8.4322967872188404e-01;
			    parameterized_quadrature_sites(10,1) = 1.3545636458303650e-01;
			    parameterized_quadrature_sites(11,0) = 1.3542317978649990e-01;
			    parameterized_quadrature_sites(11,1) = 2.1348282065620599e-02;
			    parameterized_quadrature_sites(12,0) = 3.0888535106794068e-01;
			    parameterized_quadrature_sites(12,1) = 2.2191966301360600e-02;
			    parameterized_quadrature_sites(13,0) = 6.6850575951690716e-01;
			    parameterized_quadrature_sites(13,1) = 3.0890128793894273e-01;
			    parameterized_quadrature_sites(14,0) = 2.2654501255714700e-02;
			    parameterized_quadrature_sites(14,1) = 6.6917099433209937e-01;
			    parameterized_quadrature_sites(15,0) = 2.8085154087720221e-01;
			    parameterized_quadrature_sites(15,1) = 6.9247181551062442e-01;
			    parameterized_quadrature_sites(16,0) = 6.9224467490505948e-01;
			    parameterized_quadrature_sites(16,1) = 2.6872334502594599e-02;
			    parameterized_quadrature_sites(17,0) = 2.6861744711943400e-02;
			    parameterized_quadrature_sites(17,1) = 2.8100939732219082e-01;
			    parameterized_quadrature_sites(18,0) = 1.1417784854701600e-01;
			    parameterized_quadrature_sites(18,1) = 7.9735814135857996e-01;
			    parameterized_quadrature_sites(19,0) = 7.9748079220612744e-01;
			    parameterized_quadrature_sites(19,1) = 8.7980650879088101e-02;
			    parameterized_quadrature_sites(20,0) = 8.9280729389424204e-02;
			    parameterized_quadrature_sites(20,1) = 1.1450205611275180e-01;
			    parameterized_quadrature_sites(21,0) = 1.0524878924550450e-01;
			    parameterized_quadrature_sites(21,1) = 6.6869041199220447e-01;
			    parameterized_quadrature_sites(22,0) = 6.6630222807398454e-01;
			    parameterized_quadrature_sites(22,1) = 2.2750516318320271e-01;
			    parameterized_quadrature_sites(23,0) = 2.3078037375469529e-01;
			    parameterized_quadrature_sites(23,1) = 1.0545725612213259e-01;
			    parameterized_quadrature_sites(24,0) = 1.7050591575403051e-01;
			    parameterized_quadrature_sites(24,1) = 5.1740643986577728e-01;
			    parameterized_quadrature_sites(25,0) = 5.0865939730425092e-01;
			    parameterized_quadrature_sites(25,1) = 3.1705238552093218e-01;
			    parameterized_quadrature_sites(26,0) = 3.1418238622808309e-01;
			    parameterized_quadrature_sites(26,1) = 1.8107063616590391e-01;
			    parameterized_quadrature_sites(27,0) = 4.6174608178640142e-01;
			    parameterized_quadrature_sites(27,1) = 4.6785945398040590e-01;
			    parameterized_quadrature_sites(28,0) = 6.9308749608105902e-02;
			    parameterized_quadrature_sites(28,1) = 4.6228560420845410e-01;
			    parameterized_quadrature_sites(29,0) = 4.6519552592682439e-01;
			    parameterized_quadrature_sites(29,1) = 7.2435780566898006e-02;
			    parameterized_quadrature_sites(30,0) = 2.5786258578926041e-01;
			    parameterized_quadrature_sites(30,1) = 6.1313950391771632e-01;
			    parameterized_quadrature_sites(31,0) = 6.1126277667792195e-01;
			    parameterized_quadrature_sites(31,1) = 1.3003608346093859e-01;
			    parameterized_quadrature_sites(32,0) = 1.3051821359335081e-01;
			    parameterized_quadrature_sites(32,1) = 2.5817138288836389e-01;
			    parameterized_quadrature_sites(33,0) = 4.2814379918281070e-01;
			    parameterized_quadrature_sites(33,1) = 2.3620059698167339e-01;
			    parameterized_quadrature_sites(34,0) = 3.3569957837300629e-01;
			    parameterized_quadrature_sites(34,1) = 4.3110263085883421e-01;
			    parameterized_quadrature_sites(35,0) = 2.3054242988361631e-01;
			    parameterized_quadrature_sites(35,1) = 3.4560139493758052e-01;
			    quadrature_weights(0) = 4.1560249689274751e-03;
			    quadrature_weights(1) = 4.1702924944386254e-03;
			    quadrature_weights(2) = 4.1707642266642000e-03;
			    quadrature_weights(3) = 4.3920217520637501e-03;
			    quadrature_weights(4) = 4.6118665461224503e-03;
			    quadrature_weights(5) = 4.9485602546932503e-03;
			    quadrature_weights(6) = 5.0885098963676500e-03;
			    quadrature_weights(7) = 5.1713215985039248e-03;
			    quadrature_weights(8) = 5.2067841521619253e-03;
			    quadrature_weights(9) = 7.9454944569679001e-03;
			    quadrature_weights(10) = 8.0118008810179994e-03;
			    quadrature_weights(11) = 8.0151920286416749e-03;
			    quadrature_weights(12) = 1.0769148979569475e-02;
			    quadrature_weights(13) = 1.0961835383476750e-02;
			    quadrature_weights(14) = 1.0980241818319100e-02;
			    quadrature_weights(15) = 1.1998798092264475e-02;
			    quadrature_weights(16) = 1.2095156518331500e-02;
			    quadrature_weights(17) = 1.2121685584365801e-02;
			    quadrature_weights(18) = 1.3924112200605200e-02;
			    quadrature_weights(19) = 1.4025659108902326e-02;
			    quadrature_weights(20) = 1.4129753092321276e-02;
			    quadrature_weights(21) = 1.7232247266740025e-02;
			    quadrature_weights(22) = 1.7930333402221024e-02;
			    quadrature_weights(23) = 1.8186348024405299e-02;
			    quadrature_weights(24) = 1.9720183418426299e-02;
			    quadrature_weights(25) = 2.0252858637809099e-02;
			    quadrature_weights(26) = 2.0643132476366850e-02;
			    quadrature_weights(27) = 2.1051114183261624e-02;
			    quadrature_weights(28) = 2.1089638332624925e-02;
			    quadrature_weights(29) = 2.1299246712198301e-02;
			    quadrature_weights(30) = 2.2571133201311374e-02;
			    quadrature_weights(31) = 2.2857078587113475e-02;
			    quadrature_weights(32) = 2.2906976635229850e-02;
			    quadrature_weights(33) = 2.5639334372397973e-02;
			    quadrature_weights(34) = 2.5828991535335525e-02;
			    quadrature_weights(35) = 2.5896359179831598e-02;
			    break;
			case 42:
			    parameterized_quadrature_sites(0,0) = 4.88963910362179e-01;
			    parameterized_quadrature_sites(0,1) = 4.88963910362179e-01;
			    parameterized_quadrature_sites(1,0) = 4.88963910362179e-01;
			    parameterized_quadrature_sites(1,1) = 2.20721792756430e-02;
			    parameterized_quadrature_sites(2,0) = 2.20721792756430e-02;
			    parameterized_quadrature_sites(2,1) = 4.88963910362179e-01;
			    parameterized_quadrature_sites(3,0) = 4.17644719340454e-01;
			    parameterized_quadrature_sites(3,1) = 4.17644719340454e-01;
			    parameterized_quadrature_sites(4,0) = 4.17644719340454e-01;
			    parameterized_quadrature_sites(4,1) = 1.64710561319092e-01;
			    parameterized_quadrature_sites(5,0) = 1.64710561319092e-01;
			    parameterized_quadrature_sites(5,1) = 4.17644719340454e-01;
			    parameterized_quadrature_sites(6,0) = 2.73477528308839e-01;
			    parameterized_quadrature_sites(6,1) = 2.73477528308839e-01;
			    parameterized_quadrature_sites(7,0) = 2.73477528308839e-01;
			    parameterized_quadrature_sites(7,1) = 4.53044943382323e-01;
			    parameterized_quadrature_sites(8,0) = 4.53044943382323e-01;
			    parameterized_quadrature_sites(8,1) = 2.73477528308839e-01;
			    parameterized_quadrature_sites(9,0) = 1.77205532412543e-01;
			    parameterized_quadrature_sites(9,1) = 1.77205532412543e-01;
			    parameterized_quadrature_sites(10,0) = 1.77205532412543e-01;
			    parameterized_quadrature_sites(10,1) = 6.45588935174913e-01;
			    parameterized_quadrature_sites(11,0) = 6.45588935174913e-01;
			    parameterized_quadrature_sites(11,1) = 1.77205532412543e-01;
			    parameterized_quadrature_sites(12,0) = 6.17998830908730e-02;
			    parameterized_quadrature_sites(12,1) = 6.17998830908730e-02;
			    parameterized_quadrature_sites(13,0) = 6.17998830908730e-02;
			    parameterized_quadrature_sites(13,1) = 8.76400233818255e-01;
			    parameterized_quadrature_sites(14,0) = 8.76400233818255e-01;
			    parameterized_quadrature_sites(14,1) = 6.17998830908730e-02;
			    parameterized_quadrature_sites(15,0) = 1.93909612487010e-02;
			    parameterized_quadrature_sites(15,1) = 1.93909612487010e-02;
			    parameterized_quadrature_sites(16,0) = 1.93909612487010e-02;
			    parameterized_quadrature_sites(16,1) = 9.61218077502598e-01;
			    parameterized_quadrature_sites(17,0) = 9.61218077502598e-01;
			    parameterized_quadrature_sites(17,1) = 1.93909612487010e-02;
			    parameterized_quadrature_sites(18,0) = 1.72266687821356e-01;
			    parameterized_quadrature_sites(18,1) = 7.70608554774996e-01;
			    parameterized_quadrature_sites(19,0) = 7.70608554774996e-01;
			    parameterized_quadrature_sites(19,1) = 5.71247574036480e-02;
			    parameterized_quadrature_sites(20,0) = 5.71247574036480e-02;
			    parameterized_quadrature_sites(20,1) = 1.72266687821356e-01;
			    parameterized_quadrature_sites(21,0) = 1.72266687821356e-01;
			    parameterized_quadrature_sites(21,1) = 5.71247574036480e-02;
			    parameterized_quadrature_sites(22,0) = 7.70608554774996e-01;
			    parameterized_quadrature_sites(22,1) = 1.72266687821356e-01;
			    parameterized_quadrature_sites(23,0) = 5.71247574036480e-02;
			    parameterized_quadrature_sites(23,1) = 7.70608554774996e-01;
			    parameterized_quadrature_sites(24,0) = 3.36861459796345e-01;
			    parameterized_quadrature_sites(24,1) = 5.70222290846683e-01;
			    parameterized_quadrature_sites(25,0) = 5.70222290846683e-01;
			    parameterized_quadrature_sites(25,1) = 9.29162493569720e-02;
			    parameterized_quadrature_sites(26,0) = 9.29162493569720e-02;
			    parameterized_quadrature_sites(26,1) = 3.36861459796345e-01;
			    parameterized_quadrature_sites(27,0) = 3.36861459796345e-01;
			    parameterized_quadrature_sites(27,1) = 9.29162493569720e-02;
			    parameterized_quadrature_sites(28,0) = 5.70222290846683e-01;
			    parameterized_quadrature_sites(28,1) = 3.36861459796345e-01;
			    parameterized_quadrature_sites(29,0) = 9.29162493569720e-02;
			    parameterized_quadrature_sites(29,1) = 5.70222290846683e-01;
			    parameterized_quadrature_sites(30,0) = 2.98372882136258e-01;
			    parameterized_quadrature_sites(30,1) = 6.86980167808088e-01;
			    parameterized_quadrature_sites(31,0) = 6.86980167808088e-01;
			    parameterized_quadrature_sites(31,1) = 1.46469500556540e-02;
			    parameterized_quadrature_sites(32,0) = 1.46469500556540e-02;
			    parameterized_quadrature_sites(32,1) = 2.98372882136258e-01;
			    parameterized_quadrature_sites(33,0) = 2.98372882136258e-01;
			    parameterized_quadrature_sites(33,1) = 1.46469500556540e-02;
			    parameterized_quadrature_sites(34,0) = 6.86980167808088e-01;
			    parameterized_quadrature_sites(34,1) = 2.98372882136258e-01;
			    parameterized_quadrature_sites(35,0) = 1.46469500556540e-02;
			    parameterized_quadrature_sites(35,1) = 6.86980167808088e-01;
			    parameterized_quadrature_sites(36,0) = 1.18974497696957e-01;
			    parameterized_quadrature_sites(36,1) = 8.79757171370171e-01;
			    parameterized_quadrature_sites(37,0) = 8.79757171370171e-01;
			    parameterized_quadrature_sites(37,1) = 1.26833093287199e-03;
			    parameterized_quadrature_sites(38,0) = 1.26833093287199e-03;
			    parameterized_quadrature_sites(38,1) = 1.18974497696957e-01;
			    parameterized_quadrature_sites(39,0) = 1.18974497696957e-01;
			    parameterized_quadrature_sites(39,1) = 1.26833093287199e-03;
			    parameterized_quadrature_sites(40,0) = 8.79757171370171e-01;
			    parameterized_quadrature_sites(40,1) = 1.18974497696957e-01;
			    parameterized_quadrature_sites(41,0) = 1.26833093287199e-03;
			    parameterized_quadrature_sites(41,1) = 8.79757171370171e-01;
			    quadrature_weights(0) = 1.09417906847145e-02;
			    quadrature_weights(1) = 1.09417906847145e-02;
			    quadrature_weights(2) = 1.09417906847145e-02;
			    quadrature_weights(3) = 1.63941767720625e-02;
			    quadrature_weights(4) = 1.63941767720625e-02;
			    quadrature_weights(5) = 1.63941767720625e-02;
			    quadrature_weights(6) = 2.58870522536460e-02;
			    quadrature_weights(7) = 2.58870522536460e-02;
			    quadrature_weights(8) = 2.58870522536460e-02;
			    quadrature_weights(9) = 2.10812943684965e-02;
			    quadrature_weights(10) = 2.10812943684965e-02;
			    quadrature_weights(11) = 2.10812943684965e-02;
			    quadrature_weights(12) = 7.21684983488850e-03;
			    quadrature_weights(13) = 7.21684983488850e-03;
			    quadrature_weights(14) = 7.21684983488850e-03;
			    quadrature_weights(15) = 2.46170180120000e-03;
			    quadrature_weights(16) = 2.46170180120000e-03;
			    quadrature_weights(17) = 2.46170180120000e-03;
			    quadrature_weights(18) = 1.23328766062820e-02;
			    quadrature_weights(19) = 1.23328766062820e-02;
			    quadrature_weights(20) = 1.23328766062820e-02;
			    quadrature_weights(21) = 1.23328766062820e-02;
			    quadrature_weights(22) = 1.23328766062820e-02;
			    quadrature_weights(23) = 1.23328766062820e-02;
			    quadrature_weights(24) = 1.92857553935305e-02;
			    quadrature_weights(25) = 1.92857553935305e-02;
			    quadrature_weights(26) = 1.92857553935305e-02;
			    quadrature_weights(27) = 1.92857553935305e-02;
			    quadrature_weights(28) = 1.92857553935305e-02;
			    quadrature_weights(29) = 1.92857553935305e-02;
			    quadrature_weights(30) = 7.21815405676700e-03;
			    quadrature_weights(31) = 7.21815405676700e-03;
			    quadrature_weights(32) = 7.21815405676700e-03;
			    quadrature_weights(33) = 7.21815405676700e-03;
			    quadrature_weights(34) = 7.21815405676700e-03;
			    quadrature_weights(35) = 7.21815405676700e-03;
			    quadrature_weights(36) = 2.50511441925050e-03;
			    quadrature_weights(37) = 2.50511441925050e-03;
			    quadrature_weights(38) = 2.50511441925050e-03;
			    quadrature_weights(39) = 2.50511441925050e-03;
			    quadrature_weights(40) = 2.50511441925050e-03;
			    quadrature_weights(41) = 2.50511441925050e-03;
			    break;
			case 55:
			    parameterized_quadrature_sites(0,0) = 1.0e-00;
			    parameterized_quadrature_sites(0,1) = 0.0e-00;
			    parameterized_quadrature_sites(1,0) = 0.0e-00;
			    parameterized_quadrature_sites(1,1) = 1.0e-00;
			    parameterized_quadrature_sites(2,0) = 0.0e-00;
			    parameterized_quadrature_sites(2,1) = 0.0e-00;
			    parameterized_quadrature_sites(3,0) = 9.3988635835771928e-01;
			    parameterized_quadrature_sites(3,1) = 4.9848744634100996e-03;
			    parameterized_quadrature_sites(4,0) = 5.4380668305835503e-02;
			    parameterized_quadrature_sites(4,1) = 9.3864056186166756e-01;
			    parameterized_quadrature_sites(5,0) = 9.3940049163876004e-03;
			    parameterized_quadrature_sites(5,1) = 5.2642446269734702e-02;
			    parameterized_quadrature_sites(6,0) = 1.6434508636240200e-02;
			    parameterized_quadrature_sites(6,1) = 9.4690355173508323e-01;
			    parameterized_quadrature_sites(7,0) = 9.4694872698624577e-01;
			    parameterized_quadrature_sites(7,1) = 3.6337367716693998e-02;
			    parameterized_quadrature_sites(8,0) = 4.2660400576765102e-02;
			    parameterized_quadrature_sites(8,1) = 1.5122454179941101e-02;
			    parameterized_quadrature_sites(9,0) = 1.2226949543872000e-02;
			    parameterized_quadrature_sites(9,1) = 8.6937735106643133e-01;
			    parameterized_quadrature_sites(10,0) = 8.6736965210466677e-01;
			    parameterized_quadrature_sites(10,1) = 1.2049172857742969e-01;
			    parameterized_quadrature_sites(11,0) = 8.4567440213890721e-01;
			    parameterized_quadrature_sites(11,1) = 1.5776396787000199e-02;
			    parameterized_quadrature_sites(12,0) = 1.3957596321026139e-01;
			    parameterized_quadrature_sites(12,1) = 8.4481208703747090e-01;
			    parameterized_quadrature_sites(13,0) = 1.3178217432308281e-01;
			    parameterized_quadrature_sites(13,1) = 1.3500960558402201e-02;
			    parameterized_quadrature_sites(14,0) = 1.5795512630024801e-02;
			    parameterized_quadrature_sites(14,1) = 1.4552749385359881e-01;
			    parameterized_quadrature_sites(15,0) = 7.3654628844363068e-01;
			    parameterized_quadrature_sites(15,1) = 1.5569754090822801e-02;
			    parameterized_quadrature_sites(16,0) = 1.3968843033038900e-02;
			    parameterized_quadrature_sites(16,1) = 7.3798368944501946e-01;
			    parameterized_quadrature_sites(17,0) = 2.5478951860390298e-01;
			    parameterized_quadrature_sites(17,1) = 7.2976156897705524e-01;
			    parameterized_quadrature_sites(18,0) = 7.3163865225549030e-01;
			    parameterized_quadrature_sites(18,1) = 2.5430766833150520e-01;
			    parameterized_quadrature_sites(19,0) = 1.5725372895084501e-02;
			    parameterized_quadrature_sites(19,1) = 2.6962397957906031e-01;
			    parameterized_quadrature_sites(20,0) = 2.6623028436468249e-01;
			    parameterized_quadrature_sites(20,1) = 1.4478395630801300e-02;
			    parameterized_quadrature_sites(21,0) = 8.6735040652140771e-01;
			    parameterized_quadrature_sites(21,1) = 5.9167941040048203e-02;
			    parameterized_quadrature_sites(22,0) = 7.4149366695661204e-02;
			    parameterized_quadrature_sites(22,1) = 8.6347825750608687e-01;
			    parameterized_quadrature_sites(23,0) = 1.5928594836003299e-02;
			    parameterized_quadrature_sites(23,1) = 4.1912389552381862e-01;
			    parameterized_quadrature_sites(24,0) = 1.5606102806777700e-02;
			    parameterized_quadrature_sites(24,1) = 5.8092229211457624e-01;
			    parameterized_quadrature_sites(25,0) = 5.9100948174838852e-01;
			    parameterized_quadrature_sites(25,1) = 1.5925145265094101e-02;
			    parameterized_quadrature_sites(26,0) = 4.0347714968887188e-01;
			    parameterized_quadrature_sites(26,1) = 5.8067003681039198e-01;
			    parameterized_quadrature_sites(27,0) = 5.6947456285259768e-01;
			    parameterized_quadrature_sites(27,1) = 4.1494951463020030e-01;
			    parameterized_quadrature_sites(28,0) = 6.7849370065030001e-02;
			    parameterized_quadrature_sites(28,1) = 7.6121867859137604e-02;
			    parameterized_quadrature_sites(29,0) = 4.2659685902715933e-01;
			    parameterized_quadrature_sites(29,1) = 1.5750969231154401e-02;
			    parameterized_quadrature_sites(30,0) = 6.7098250788970207e-02;
			    parameterized_quadrature_sites(30,1) = 7.7418983124212093e-01;
			    parameterized_quadrature_sites(31,0) = 7.5283102314795158e-01;
			    parameterized_quadrature_sites(31,1) = 8.1911949563924294e-02;
			    parameterized_quadrature_sites(32,0) = 7.7537277835568841e-01;
			    parameterized_quadrature_sites(32,1) = 1.5771284572917341e-01;
			    parameterized_quadrature_sites(33,0) = 1.6890731577873661e-01;
			    parameterized_quadrature_sites(33,1) = 7.5039430997422452e-01;
			    parameterized_quadrature_sites(34,0) = 1.6873358329194171e-01;
			    parameterized_quadrature_sites(34,1) = 7.0831150726781894e-02;
			    parameterized_quadrature_sites(35,0) = 8.2124470843632405e-02;
			    parameterized_quadrature_sites(35,1) = 1.7629966267710759e-01;
			    parameterized_quadrature_sites(36,0) = 6.2887053633447976e-01;
			    parameterized_quadrature_sites(36,1) = 8.0774495331656301e-02;
			    parameterized_quadrature_sites(37,0) = 8.1141301526575199e-02;
			    parameterized_quadrature_sites(37,1) = 3.0543735897757762e-01;
			    parameterized_quadrature_sites(38,0) = 2.9691120650804809e-01;
			    parameterized_quadrature_sites(38,1) = 6.2274859888709300e-01;
			    parameterized_quadrature_sites(39,0) = 7.6754231417057298e-02;
			    parameterized_quadrature_sites(39,1) = 6.2472471495456661e-01;
			    parameterized_quadrature_sites(40,0) = 6.2230223338447721e-01;
			    parameterized_quadrature_sites(40,1) = 3.0114858211656370e-01;
			    parameterized_quadrature_sites(41,0) = 3.1037862880509631e-01;
			    parameterized_quadrature_sites(41,1) = 7.7909836507944599e-02;
			    parameterized_quadrature_sites(42,0) = 8.1921821518658594e-02;
			    parameterized_quadrature_sites(42,1) = 4.6036330383508761e-01;
			    parameterized_quadrature_sites(43,0) = 4.7170226650134689e-01;
			    parameterized_quadrature_sites(43,1) = 8.2155400679671906e-02;
			    parameterized_quadrature_sites(44,0) = 4.5466034152504742e-01;
			    parameterized_quadrature_sites(44,1) = 4.6375650338896440e-01;
			    parameterized_quadrature_sites(45,0) = 1.7010913392369389e-01;
			    parameterized_quadrature_sites(45,1) = 6.4222778081881993e-01;
			    parameterized_quadrature_sites(46,0) = 6.4060043294867430e-01;
			    parameterized_quadrature_sites(46,1) = 1.8982935372556059e-01;
			    parameterized_quadrature_sites(47,0) = 1.9122675837165989e-01;
			    parameterized_quadrature_sites(47,1) = 1.7399556853425760e-01;
			    parameterized_quadrature_sites(48,0) = 1.8853157670702370e-01;
			    parameterized_quadrature_sites(48,1) = 4.7989140704057581e-01;
			    parameterized_quadrature_sites(49,0) = 4.7729299576907452e-01;
			    parameterized_quadrature_sites(49,1) = 3.3483565981193042e-01;
			    parameterized_quadrature_sites(50,0) = 3.1269746217597721e-01;
			    parameterized_quadrature_sites(50,1) = 4.9579721972587398e-01;
			    parameterized_quadrature_sites(51,0) = 4.9612259459456259e-01;
			    parameterized_quadrature_sites(51,1) = 1.9275536689044351e-01;
			    parameterized_quadrature_sites(52,0) = 1.9288053128670610e-01;
			    parameterized_quadrature_sites(52,1) = 3.1610158072607569e-01;
			    parameterized_quadrature_sites(53,0) = 3.3600414538164958e-01;
			    parameterized_quadrature_sites(53,1) = 1.8948928012898231e-01;
			    parameterized_quadrature_sites(54,0) = 3.3372805508479741e-01;
			    parameterized_quadrature_sites(54,1) = 3.3435710218114523e-01;
			    quadrature_weights(0) = 1.5506499627784999e-04;
			    quadrature_weights(1) = 1.5787936779320001e-04;
			    quadrature_weights(2) = 1.7716503897180001e-04;
			    quadrature_weights(3) = 1.3790929042020500e-03;
			    quadrature_weights(4) = 1.5673101913945000e-03;
			    quadrature_weights(5) = 1.9632852206504501e-03;
			    quadrature_weights(6) = 2.3637870966120248e-03;
			    quadrature_weights(7) = 2.4456127817772499e-03;
			    quadrature_weights(8) = 2.4965410872361499e-03;
			    quadrature_weights(9) = 3.4388454704036252e-03;
			    quadrature_weights(10) = 3.5244794510020999e-03;
			    quadrature_weights(11) = 3.7411716084290001e-03;
			    quadrature_weights(12) = 3.9024375902997751e-03;
			    quadrature_weights(13) = 3.9420923337041246e-03;
			    quadrature_weights(14) = 4.3948636595678749e-03;
			    quadrature_weights(15) = 5.1028460067507001e-03;
			    quadrature_weights(16) = 5.2390719653994996e-03;
			    quadrature_weights(17) = 5.2678353249450997e-03;
			    quadrature_weights(18) = 5.4411690050507498e-03;
			    quadrature_weights(19) = 5.5572102174651500e-03;
			    quadrature_weights(20) = 5.6046673420516247e-03;
			    quadrature_weights(21) = 5.7530654248287997e-03;
			    quadrature_weights(22) = 5.9203475624943751e-03;
			    quadrature_weights(23) = 6.4366160841976749e-03;
			    quadrature_weights(24) = 6.4489200402012004e-03;
			    quadrature_weights(25) = 6.4518081902498253e-03;
			    quadrature_weights(26) = 6.5085808014669752e-03;
			    quadrature_weights(27) = 6.6442035402279253e-03;
			    quadrature_weights(28) = 6.6446190457719254e-03;
			    quadrature_weights(29) = 6.6883082309448748e-03;
			    quadrature_weights(30) = 9.3946951660218506e-03;
			    quadrature_weights(31) = 9.5766473548822492e-03;
			    quadrature_weights(32) = 9.6212423756326000e-03;
			    quadrature_weights(33) = 9.7404956463108747e-03;
			    quadrature_weights(34) = 9.8651027886874250e-03;
			    quadrature_weights(35) = 1.0309119452445125e-02;
			    quadrature_weights(36) = 1.2821810962084575e-02;
			    quadrature_weights(37) = 1.2910141048365926e-02;
			    quadrature_weights(38) = 1.2955751056727700e-02;
			    quadrature_weights(39) = 1.3213199704525376e-02;
			    quadrature_weights(40) = 1.3462639325681701e-02;
			    quadrature_weights(41) = 1.3547383232981901e-02;
			    quadrature_weights(42) = 1.4618428661110900e-02;
			    quadrature_weights(43) = 1.4821579209082151e-02;
			    quadrature_weights(44) = 1.4858956918716301e-02;
			    quadrature_weights(45) = 1.5795006396574451e-02;
			    quadrature_weights(46) = 1.5823171128833101e-02;
			    quadrature_weights(47) = 1.6017684044289202e-02;
			    quadrature_weights(48) = 2.0301014897957576e-02;
			    quadrature_weights(49) = 2.0360937838258826e-02;
			    quadrature_weights(50) = 2.0366980031034525e-02;
			    quadrature_weights(51) = 2.0376263702112225e-02;
			    quadrature_weights(52) = 2.0379116623473949e-02;
			    quadrature_weights(53) = 2.0423276490578225e-02;
			    quadrature_weights(54) = 2.3080458363263227e-02;
			    break;
			case 61:
			    parameterized_quadrature_sites(0,0) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(0,1) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(1,0) = 4.97170540556774e-01;
			    parameterized_quadrature_sites(1,1) = 4.97170540556774e-01;
			    parameterized_quadrature_sites(2,0) = 4.97170540556774e-01;
			    parameterized_quadrature_sites(2,1) = 5.65891888645198e-03;
			    parameterized_quadrature_sites(3,0) = 5.65891888645198e-03;
			    parameterized_quadrature_sites(3,1) = 4.97170540556774e-01;
			    parameterized_quadrature_sites(4,0) = 4.82176322624625e-01;
			    parameterized_quadrature_sites(4,1) = 4.82176322624625e-01;
			    parameterized_quadrature_sites(5,0) = 4.82176322624625e-01;
			    parameterized_quadrature_sites(5,1) = 3.56473547507510e-02;
			    parameterized_quadrature_sites(6,0) = 3.56473547507510e-02;
			    parameterized_quadrature_sites(6,1) = 4.82176322624625e-01;
			    parameterized_quadrature_sites(7,0) = 4.50239969020782e-01;
			    parameterized_quadrature_sites(7,1) = 4.50239969020782e-01;
			    parameterized_quadrature_sites(8,0) = 4.50239969020782e-01;
			    parameterized_quadrature_sites(8,1) = 9.95200619584370e-02;
			    parameterized_quadrature_sites(9,0) = 9.95200619584370e-02;
			    parameterized_quadrature_sites(9,1) = 4.50239969020782e-01;
			    parameterized_quadrature_sites(10,0) = 4.00266239377397e-01;
			    parameterized_quadrature_sites(10,1) = 4.00266239377397e-01;
			    parameterized_quadrature_sites(11,0) = 4.00266239377397e-01;
			    parameterized_quadrature_sites(11,1) = 1.99467521245206e-01;
			    parameterized_quadrature_sites(12,0) = 1.99467521245206e-01;
			    parameterized_quadrature_sites(12,1) = 4.00266239377397e-01;
			    parameterized_quadrature_sites(13,0) = 2.52141267970953e-01;
			    parameterized_quadrature_sites(13,1) = 2.52141267970953e-01;
			    parameterized_quadrature_sites(14,0) = 2.52141267970953e-01;
			    parameterized_quadrature_sites(14,1) = 4.95717464058095e-01;
			    parameterized_quadrature_sites(15,0) = 4.95717464058095e-01;
			    parameterized_quadrature_sites(15,1) = 2.52141267970953e-01;
			    parameterized_quadrature_sites(16,0) = 1.62047004658461e-01;
			    parameterized_quadrature_sites(16,1) = 1.62047004658461e-01;
			    parameterized_quadrature_sites(17,0) = 1.62047004658461e-01;
			    parameterized_quadrature_sites(17,1) = 6.75905990683077e-01;
			    parameterized_quadrature_sites(18,0) = 6.75905990683077e-01;
			    parameterized_quadrature_sites(18,1) = 1.62047004658461e-01;
			    parameterized_quadrature_sites(19,0) = 7.58758822607460e-02;
			    parameterized_quadrature_sites(19,1) = 7.58758822607460e-02;
			    parameterized_quadrature_sites(20,0) = 7.58758822607460e-02;
			    parameterized_quadrature_sites(20,1) = 8.48248235478508e-01;
			    parameterized_quadrature_sites(21,0) = 8.48248235478508e-01;
			    parameterized_quadrature_sites(21,1) = 7.58758822607460e-02;
			    parameterized_quadrature_sites(22,0) = 1.56547269678220e-02;
			    parameterized_quadrature_sites(22,1) = 1.56547269678220e-02;
			    parameterized_quadrature_sites(23,0) = 1.56547269678220e-02;
			    parameterized_quadrature_sites(23,1) = 9.68690546064356e-01;
			    parameterized_quadrature_sites(24,0) = 9.68690546064356e-01;
			    parameterized_quadrature_sites(24,1) = 1.56547269678220e-02;
			    parameterized_quadrature_sites(25,0) = 3.34319867363658e-01;
			    parameterized_quadrature_sites(25,1) = 6.55493203809423e-01;
			    parameterized_quadrature_sites(26,0) = 6.55493203809423e-01;
			    parameterized_quadrature_sites(26,1) = 1.01869288269190e-02;
			    parameterized_quadrature_sites(27,0) = 1.01869288269190e-02;
			    parameterized_quadrature_sites(27,1) = 3.34319867363658e-01;
			    parameterized_quadrature_sites(28,0) = 3.34319867363658e-01;
			    parameterized_quadrature_sites(28,1) = 1.01869288269190e-02;
			    parameterized_quadrature_sites(29,0) = 6.55493203809423e-01;
			    parameterized_quadrature_sites(29,1) = 3.34319867363658e-01;
			    parameterized_quadrature_sites(30,0) = 1.01869288269190e-02;
			    parameterized_quadrature_sites(30,1) = 6.55493203809423e-01;
			    parameterized_quadrature_sites(31,0) = 2.92221537796944e-01;
			    parameterized_quadrature_sites(31,1) = 5.72337590532020e-01;
			    parameterized_quadrature_sites(32,0) = 5.72337590532020e-01;
			    parameterized_quadrature_sites(32,1) = 1.35440871671036e-01;
			    parameterized_quadrature_sites(33,0) = 1.35440871671036e-01;
			    parameterized_quadrature_sites(33,1) = 2.92221537796944e-01;
			    parameterized_quadrature_sites(34,0) = 2.92221537796944e-01;
			    parameterized_quadrature_sites(34,1) = 1.35440871671036e-01;
			    parameterized_quadrature_sites(35,0) = 5.72337590532020e-01;
			    parameterized_quadrature_sites(35,1) = 2.92221537796944e-01;
			    parameterized_quadrature_sites(36,0) = 1.35440871671036e-01;
			    parameterized_quadrature_sites(36,1) = 5.72337590532020e-01;
			    parameterized_quadrature_sites(37,0) = 3.19574885423190e-01;
			    parameterized_quadrature_sites(37,1) = 6.26001190286228e-01;
			    parameterized_quadrature_sites(38,0) = 6.26001190286228e-01;
			    parameterized_quadrature_sites(38,1) = 5.44239242905830e-02;
			    parameterized_quadrature_sites(39,0) = 5.44239242905830e-02;
			    parameterized_quadrature_sites(39,1) = 3.19574885423190e-01;
			    parameterized_quadrature_sites(40,0) = 3.19574885423190e-01;
			    parameterized_quadrature_sites(40,1) = 5.44239242905830e-02;
			    parameterized_quadrature_sites(41,0) = 6.26001190286228e-01;
			    parameterized_quadrature_sites(41,1) = 3.19574885423190e-01;
			    parameterized_quadrature_sites(42,0) = 5.44239242905830e-02;
			    parameterized_quadrature_sites(42,1) = 6.26001190286228e-01;
			    parameterized_quadrature_sites(43,0) = 1.90704224192292e-01;
			    parameterized_quadrature_sites(43,1) = 7.96427214974071e-01;
			    parameterized_quadrature_sites(44,0) = 7.96427214974071e-01;
			    parameterized_quadrature_sites(44,1) = 1.28685608336370e-02;
			    parameterized_quadrature_sites(45,0) = 1.28685608336370e-02;
			    parameterized_quadrature_sites(45,1) = 1.90704224192292e-01;
			    parameterized_quadrature_sites(46,0) = 1.90704224192292e-01;
			    parameterized_quadrature_sites(46,1) = 1.28685608336370e-02;
			    parameterized_quadrature_sites(47,0) = 7.96427214974071e-01;
			    parameterized_quadrature_sites(47,1) = 1.90704224192292e-01;
			    parameterized_quadrature_sites(48,0) = 1.28685608336370e-02;
			    parameterized_quadrature_sites(48,1) = 7.96427214974071e-01;
			    parameterized_quadrature_sites(49,0) = 1.80483211648746e-01;
			    parameterized_quadrature_sites(49,1) = 7.52351005937729e-01;
			    parameterized_quadrature_sites(50,0) = 7.52351005937729e-01;
			    parameterized_quadrature_sites(50,1) = 6.71657824135240e-02;
			    parameterized_quadrature_sites(51,0) = 6.71657824135240e-02;
			    parameterized_quadrature_sites(51,1) = 1.80483211648746e-01;
			    parameterized_quadrature_sites(52,0) = 1.80483211648746e-01;
			    parameterized_quadrature_sites(52,1) = 6.71657824135240e-02;
			    parameterized_quadrature_sites(53,0) = 7.52351005937729e-01;
			    parameterized_quadrature_sites(53,1) = 1.80483211648746e-01;
			    parameterized_quadrature_sites(54,0) = 6.71657824135240e-02;
			    parameterized_quadrature_sites(54,1) = 7.52351005937729e-01;
			    parameterized_quadrature_sites(55,0) = 8.07113136795640e-02;
			    parameterized_quadrature_sites(55,1) = 9.04625504095608e-01;
			    parameterized_quadrature_sites(56,0) = 9.04625504095608e-01;
			    parameterized_quadrature_sites(56,1) = 1.46631822248280e-02;
			    parameterized_quadrature_sites(57,0) = 1.46631822248280e-02;
			    parameterized_quadrature_sites(57,1) = 8.07113136795640e-02;
			    parameterized_quadrature_sites(58,0) = 8.07113136795640e-02;
			    parameterized_quadrature_sites(58,1) = 1.46631822248280e-02;
			    parameterized_quadrature_sites(59,0) = 9.04625504095608e-01;
			    parameterized_quadrature_sites(59,1) = 8.07113136795640e-02;
			    parameterized_quadrature_sites(60,0) = 1.46631822248280e-02;
			    parameterized_quadrature_sites(60,1) = 9.04625504095608e-01;
			    quadrature_weights(0) = 1.67185996454015e-02;
			    quadrature_weights(1) = 2.54670772025350e-03;
			    quadrature_weights(2) = 2.54670772025350e-03;
			    quadrature_weights(3) = 2.54670772025350e-03;
			    quadrature_weights(4) = 7.33543226381900e-03;
			    quadrature_weights(5) = 7.33543226381900e-03;
			    quadrature_weights(6) = 7.33543226381900e-03;
			    quadrature_weights(7) = 1.21754391768360e-02;
			    quadrature_weights(8) = 1.21754391768360e-02;
			    quadrature_weights(9) = 1.21754391768360e-02;
			    quadrature_weights(10) = 1.55537754344845e-02;
			    quadrature_weights(11) = 1.55537754344845e-02;
			    quadrature_weights(12) = 1.55537754344845e-02;
			    quadrature_weights(13) = 1.56285556093100e-02;
			    quadrature_weights(14) = 1.56285556093100e-02;
			    quadrature_weights(15) = 1.56285556093100e-02;
			    quadrature_weights(16) = 1.24078271698325e-02;
			    quadrature_weights(17) = 1.24078271698325e-02;
			    quadrature_weights(18) = 1.24078271698325e-02;
			    quadrature_weights(19) = 7.02803653527850e-03;
			    quadrature_weights(20) = 7.02803653527850e-03;
			    quadrature_weights(21) = 7.02803653527850e-03;
			    quadrature_weights(22) = 1.59733808688950e-03;
			    quadrature_weights(23) = 1.59733808688950e-03;
			    quadrature_weights(24) = 1.59733808688950e-03;
			    quadrature_weights(25) = 4.05982765949650e-03;
			    quadrature_weights(26) = 4.05982765949650e-03;
			    quadrature_weights(27) = 4.05982765949650e-03;
			    quadrature_weights(28) = 4.05982765949650e-03;
			    quadrature_weights(29) = 4.05982765949650e-03;
			    quadrature_weights(30) = 4.05982765949650e-03;
			    quadrature_weights(31) = 1.34028711415815e-02;
			    quadrature_weights(32) = 1.34028711415815e-02;
			    quadrature_weights(33) = 1.34028711415815e-02;
			    quadrature_weights(34) = 1.34028711415815e-02;
			    quadrature_weights(35) = 1.34028711415815e-02;
			    quadrature_weights(36) = 1.34028711415815e-02;
			    quadrature_weights(37) = 9.22999660541100e-03;
			    quadrature_weights(38) = 9.22999660541100e-03;
			    quadrature_weights(39) = 9.22999660541100e-03;
			    quadrature_weights(40) = 9.22999660541100e-03;
			    quadrature_weights(41) = 9.22999660541100e-03;
			    quadrature_weights(42) = 9.22999660541100e-03;
			    quadrature_weights(43) = 4.23843426716400e-03;
			    quadrature_weights(44) = 4.23843426716400e-03;
			    quadrature_weights(45) = 4.23843426716400e-03;
			    quadrature_weights(46) = 4.23843426716400e-03;
			    quadrature_weights(47) = 4.23843426716400e-03;
			    quadrature_weights(48) = 4.23843426716400e-03;
			    quadrature_weights(49) = 9.14639838501250e-03;
			    quadrature_weights(50) = 9.14639838501250e-03;
			    quadrature_weights(51) = 9.14639838501250e-03;
			    quadrature_weights(52) = 9.14639838501250e-03;
			    quadrature_weights(53) = 9.14639838501250e-03;
			    quadrature_weights(54) = 9.14639838501250e-03;
			    quadrature_weights(55) = 3.33281600208250e-03;
			    quadrature_weights(56) = 3.33281600208250e-03;
			    quadrature_weights(57) = 3.33281600208250e-03;
			    quadrature_weights(58) = 3.33281600208250e-03;
			    quadrature_weights(59) = 3.33281600208250e-03;
			    quadrature_weights(60) = 3.33281600208250e-03;
			    break;
			case 66:
			    parameterized_quadrature_sites(0,0) = 1.1673105966841200e-02;
			    parameterized_quadrature_sites(0,1) = 9.8125659512890129e-01;
			    parameterized_quadrature_sites(1,0) = 9.8100308583879503e-01;
			    parameterized_quadrature_sites(1,1) = 7.1462504863216000e-03;
			    parameterized_quadrature_sites(2,0) = 1.0696631709169700e-02;
			    parameterized_quadrature_sites(2,1) = 1.1515393337596600e-02;
			    parameterized_quadrature_sites(3,0) = 9.3824769835505051e-01;
			    parameterized_quadrature_sites(3,1) = 4.9557059134064198e-02;
			    parameterized_quadrature_sites(4,0) = 1.2662751841721401e-02;
			    parameterized_quadrature_sites(4,1) = 9.3701236206150318e-01;
			    parameterized_quadrature_sites(5,0) = 5.9810940998380198e-02;
			    parameterized_quadrature_sites(5,1) = 1.2136457892184800e-02;
			    parameterized_quadrature_sites(6,0) = 1.3736329792672100e-02;
			    parameterized_quadrature_sites(6,1) = 6.1278362559696799e-02;
			    parameterized_quadrature_sites(7,0) = 9.2295279594054480e-01;
			    parameterized_quadrature_sites(7,1) = 1.4112827060242099e-02;
			    parameterized_quadrature_sites(8,0) = 6.3310735499269494e-02;
			    parameterized_quadrature_sites(8,1) = 9.2201972917274344e-01;
			    parameterized_quadrature_sites(9,0) = 1.1726510033460201e-02;
			    parameterized_quadrature_sites(9,1) = 1.5005204752290349e-01;
			    parameterized_quadrature_sites(10,0) = 1.5547205873234721e-01;
			    parameterized_quadrature_sites(10,1) = 8.3251471215892492e-01;
			    parameterized_quadrature_sites(11,0) = 8.3432938889821573e-01;
			    parameterized_quadrature_sites(11,1) = 1.2522815875883600e-02;
			    parameterized_quadrature_sites(12,0) = 8.5016380319567597e-01;
			    parameterized_quadrature_sites(12,1) = 1.3719975087357841e-01;
			    parameterized_quadrature_sites(13,0) = 1.2881635052197599e-02;
			    parameterized_quadrature_sites(13,1) = 8.4776270634792006e-01;
			    parameterized_quadrature_sites(14,0) = 1.5108016089587781e-01;
			    parameterized_quadrature_sites(14,1) = 1.3652692403937501e-02;
			    parameterized_quadrature_sites(15,0) = 1.0191787921658400e-02;
			    parameterized_quadrature_sites(15,1) = 5.7704386183448575e-01;
			    parameterized_quadrature_sites(16,0) = 2.8133723993032811e-01;
			    parameterized_quadrature_sites(16,1) = 7.0668537596231984e-01;
			    parameterized_quadrature_sites(17,0) = 7.1243746285009335e-01;
			    parameterized_quadrature_sites(17,1) = 1.2456978098990301e-02;
			    parameterized_quadrature_sites(18,0) = 2.7630252508633668e-01;
			    parameterized_quadrature_sites(18,1) = 1.2174131138564200e-02;
			    parameterized_quadrature_sites(19,0) = 1.0965836856061799e-02;
			    parameterized_quadrature_sites(19,1) = 4.1943067124662847e-01;
			    parameterized_quadrature_sites(20,0) = 4.2891105178839167e-01;
			    parameterized_quadrature_sites(20,1) = 5.5996160674689166e-01;
			    parameterized_quadrature_sites(21,0) = 4.2154205551147350e-01;
			    parameterized_quadrature_sites(21,1) = 1.1647599478465699e-02;
			    parameterized_quadrature_sites(22,0) = 5.7112585904443613e-01;
			    parameterized_quadrature_sites(22,1) = 1.1821831398851200e-02;
			    parameterized_quadrature_sites(23,0) = 5.8268682705109343e-01;
			    parameterized_quadrature_sites(23,1) = 4.0578895811771831e-01;
			    parameterized_quadrature_sites(24,0) = 1.3056780671324699e-02;
			    parameterized_quadrature_sites(24,1) = 2.7250237508679159e-01;
			    parameterized_quadrature_sites(25,0) = 1.3076040096391800e-02;
			    parameterized_quadrature_sites(25,1) = 7.2247125232334730e-01;
			    parameterized_quadrature_sites(26,0) = 7.2634370624067746e-01;
			    parameterized_quadrature_sites(26,1) = 2.6029840192506443e-01;
			    parameterized_quadrature_sites(27,0) = 6.8723006863737404e-02;
			    parameterized_quadrature_sites(27,1) = 6.3141727720962701e-02;
			    parameterized_quadrature_sites(28,0) = 8.6523021015294610e-01;
			    parameterized_quadrature_sites(28,1) = 7.2061183733767895e-02;
			    parameterized_quadrature_sites(29,0) = 6.4859907103736694e-02;
			    parameterized_quadrature_sites(29,1) = 8.5904335439099433e-01;
			    parameterized_quadrature_sites(30,0) = 1.4834949433620581e-01;
			    parameterized_quadrature_sites(30,1) = 7.8887883522396707e-01;
			    parameterized_quadrature_sites(31,0) = 6.2435989839593801e-02;
			    parameterized_quadrature_sites(31,1) = 1.4939354993542750e-01;
			    parameterized_quadrature_sites(32,0) = 7.8713690117350699e-01;
			    parameterized_quadrature_sites(32,1) = 6.5638204275659501e-02;
			    parameterized_quadrature_sites(33,0) = 5.1910492160953101e-02;
			    parameterized_quadrature_sites(33,1) = 5.2556356956052430e-01;
			    parameterized_quadrature_sites(34,0) = 1.5431299274438229e-01;
			    parameterized_quadrature_sites(34,1) = 7.1638392691700595e-02;
			    parameterized_quadrature_sites(35,0) = 2.6178427456029407e-01;
			    parameterized_quadrature_sites(35,1) = 6.2147948528815097e-02;
			    parameterized_quadrature_sites(36,0) = 7.6672578728127994e-01;
			    parameterized_quadrature_sites(36,1) = 1.6582115548313259e-01;
			    parameterized_quadrature_sites(37,0) = 2.5821036766273009e-01;
			    parameterized_quadrature_sites(37,1) = 6.8001197661390189e-01;
			    parameterized_quadrature_sites(38,0) = 6.7906592514742597e-02;
			    parameterized_quadrature_sites(38,1) = 7.5715154377818017e-01;
			    parameterized_quadrature_sites(39,0) = 5.2935782748041971e-01;
			    parameterized_quadrature_sites(39,1) = 4.1215038411072058e-01;
			    parameterized_quadrature_sites(40,0) = 6.6603615048415998e-02;
			    parameterized_quadrature_sites(40,1) = 2.6125130878864999e-01;
			    parameterized_quadrature_sites(41,0) = 5.8567546189943198e-02;
			    parameterized_quadrature_sites(41,1) = 3.9022361145349760e-01;
			    parameterized_quadrature_sites(42,0) = 6.4453536041083406e-02;
			    parameterized_quadrature_sites(42,1) = 6.3736265597609554e-01;
			    parameterized_quadrature_sites(43,0) = 6.7481384291513691e-01;
			    parameterized_quadrature_sites(43,1) = 6.3758334206129100e-02;
			    parameterized_quadrature_sites(44,0) = 3.9146023103687089e-01;
			    parameterized_quadrature_sites(44,1) = 5.5032380905631106e-01;
			    parameterized_quadrature_sites(45,0) = 6.4877014923071408e-01;
			    parameterized_quadrature_sites(45,1) = 2.8367283602629478e-01;
			    parameterized_quadrature_sites(46,0) = 3.9464982204080379e-01;
			    parameterized_quadrature_sites(46,1) = 6.0517552255370900e-02;
			    parameterized_quadrature_sites(47,0) = 5.3901371519333352e-01;
			    parameterized_quadrature_sites(47,1) = 6.1199017693642201e-02;
			    parameterized_quadrature_sites(48,0) = 1.6278950827847499e-01;
			    parameterized_quadrature_sites(48,1) = 6.8613221410348235e-01;
			    parameterized_quadrature_sites(49,0) = 6.8124363226406448e-01;
			    parameterized_quadrature_sites(49,1) = 1.5679683458990931e-01;
			    parameterized_quadrature_sites(50,0) = 1.5428328780201980e-01;
			    parameterized_quadrature_sites(50,1) = 1.6675126240198401e-01;
			    parameterized_quadrature_sites(51,0) = 2.5227277504445078e-01;
			    parameterized_quadrature_sites(51,1) = 2.5048039333948502e-01;
			    parameterized_quadrature_sites(52,0) = 2.5479815324070432e-01;
			    parameterized_quadrature_sites(52,1) = 4.9940906490431908e-01;
			    parameterized_quadrature_sites(53,0) = 1.4855805491943541e-01;
			    parameterized_quadrature_sites(53,1) = 5.7560230960873771e-01;
			    parameterized_quadrature_sites(54,0) = 2.9302396064361819e-01;
			    parameterized_quadrature_sites(54,1) = 5.6568973541618528e-01;
			    parameterized_quadrature_sites(55,0) = 2.8089912723099042e-01;
			    parameterized_quadrature_sites(55,1) = 1.4379215742477949e-01;
			    parameterized_quadrature_sites(56,0) = 4.8209895929708219e-01;
			    parameterized_quadrature_sites(56,1) = 2.5185575358650381e-01;
			    parameterized_quadrature_sites(57,0) = 5.6418782454436134e-01;
			    parameterized_quadrature_sites(57,1) = 1.4629667431525920e-01;
			    parameterized_quadrature_sites(58,0) = 1.3076996443439021e-01;
			    parameterized_quadrature_sites(58,1) = 4.4895775861167753e-01;
			    parameterized_quadrature_sites(59,0) = 1.4796922219475581e-01;
			    parameterized_quadrature_sites(59,1) = 3.0011743868291701e-01;
			    parameterized_quadrature_sites(60,0) = 5.6386842229459166e-01;
			    parameterized_quadrature_sites(60,1) = 2.8137720892975088e-01;
			    parameterized_quadrature_sites(61,0) = 4.3611574287904659e-01;
			    parameterized_quadrature_sites(61,1) = 4.2520534464204729e-01;
			    parameterized_quadrature_sites(62,0) = 3.6032639352854701e-01;
			    parameterized_quadrature_sites(62,1) = 2.5991900048886368e-01;
			    parameterized_quadrature_sites(63,0) = 4.2241883346742481e-01;
			    parameterized_quadrature_sites(63,1) = 1.4532384433026860e-01;
			    parameterized_quadrature_sites(64,0) = 3.7190018330523877e-01;
			    parameterized_quadrature_sites(64,1) = 3.7801227035670099e-01;
			    parameterized_quadrature_sites(65,0) = 2.4136450069284729e-01;
			    parameterized_quadrature_sites(65,1) = 3.8475632849397318e-01;
			    quadrature_weights(0) = 6.2914392466127504e-04;
			    quadrature_weights(1) = 6.3183630018052495e-04;
			    quadrature_weights(2) = 8.3173238332957496e-04;
			    quadrature_weights(3) = 2.0375873031348501e-03;
			    quadrature_weights(4) = 2.1533881435404252e-03;
			    quadrature_weights(5) = 2.1946686544827001e-03;
			    quadrature_weights(6) = 2.4274896390419250e-03;
			    quadrature_weights(7) = 2.5616552978717250e-03;
			    quadrature_weights(8) = 2.7099422085185751e-03;
			    quadrature_weights(9) = 3.2346347543960748e-03;
			    quadrature_weights(10) = 3.4084955895738499e-03;
			    quadrature_weights(11) = 3.4619332036664001e-03;
			    quadrature_weights(12) = 3.4855385026212752e-03;
			    quadrature_weights(13) = 3.6030349991898998e-03;
			    quadrature_weights(14) = 3.8425863883506501e-03;
			    quadrature_weights(15) = 4.0622450563140748e-03;
			    quadrature_weights(16) = 4.2429576070036248e-03;
			    quadrature_weights(17) = 4.2522133105332002e-03;
			    quadrature_weights(18) = 4.2738380168663498e-03;
			    quadrature_weights(19) = 4.3472213639774750e-03;
			    quadrature_weights(20) = 4.3635990609680250e-03;
			    quadrature_weights(21) = 4.4601689321659248e-03;
			    quadrature_weights(22) = 4.4611715969841001e-03;
			    quadrature_weights(23) = 4.4761584388088504e-03;
			    quadrature_weights(24) = 4.5314939050175749e-03;
			    quadrature_weights(25) = 4.6196209720506748e-03;
			    quadrature_weights(26) = 4.6448391092779500e-03;
			    quadrature_weights(27) = 5.0804287944135754e-03;
			    quadrature_weights(28) = 5.3442915452291753e-03;
			    quadrature_weights(29) = 5.7979213524570498e-03;
			    quadrature_weights(30) = 6.8606677714776751e-03;
			    quadrature_weights(31) = 7.2575480585093501e-03;
			    quadrature_weights(32) = 7.3630684626356750e-03;
			    quadrature_weights(33) = 7.4859062907272752e-03;
			    quadrature_weights(34) = 7.6756737029697996e-03;
			    quadrature_weights(35) = 8.1315841465677500e-03;
			    quadrature_weights(36) = 8.1971052126526001e-03;
			    quadrature_weights(37) = 8.2808668797998992e-03;
			    quadrature_weights(38) = 8.6541881718643493e-03;
			    quadrature_weights(39) = 8.6770343494031749e-03;
			    quadrature_weights(40) = 8.6843012350962742e-03;
			    quadrature_weights(41) = 8.7132190613550004e-03;
			    quadrature_weights(42) = 8.7150390296491503e-03;
			    quadrature_weights(43) = 8.8867892493721002e-03;
			    quadrature_weights(44) = 9.0045749095674504e-03;
			    quadrature_weights(45) = 9.0731571460701751e-03;
			    quadrature_weights(46) = 9.5474425520800498e-03;
			    quadrature_weights(47) = 9.8063200029475748e-03;
			    quadrature_weights(48) = 1.2067753147187676e-02;
			    quadrature_weights(49) = 1.2247803039156025e-02;
			    quadrature_weights(50) = 1.2430520846804525e-02;
			    quadrature_weights(51) = 1.2676643424644550e-02;
			    quadrature_weights(52) = 1.2744299851074425e-02;
			    quadrature_weights(53) = 1.3034001591679976e-02;
			    quadrature_weights(54) = 1.3086521873117450e-02;
			    quadrature_weights(55) = 1.3111017088793300e-02;
			    quadrature_weights(56) = 1.3186491120564975e-02;
			    quadrature_weights(57) = 1.3236226593190750e-02;
			    quadrature_weights(58) = 1.3559889862520475e-02;
			    quadrature_weights(59) = 1.3586755085482525e-02;
			    quadrature_weights(60) = 1.3677513715971475e-02;
			    quadrature_weights(61) = 1.3932208647817225e-02;
			    quadrature_weights(62) = 1.4443356605827800e-02;
			    quadrature_weights(63) = 1.4634844540566926e-02;
			    quadrature_weights(64) = 1.5225981266989826e-02;
			    quadrature_weights(65) = 1.5931849111237351e-02;
			    break;
			case 73:
			    parameterized_quadrature_sites(0,0) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(0,1) = 3.33333333333333e-01;
			    parameterized_quadrature_sites(1,0) = 4.89609987073006e-01;
			    parameterized_quadrature_sites(1,1) = 4.89609987073006e-01;
			    parameterized_quadrature_sites(2,0) = 4.89609987073006e-01;
			    parameterized_quadrature_sites(2,1) = 2.07800258539870e-02;
			    parameterized_quadrature_sites(3,0) = 2.07800258539870e-02;
			    parameterized_quadrature_sites(3,1) = 4.89609987073006e-01;
			    parameterized_quadrature_sites(4,0) = 4.54536892697893e-01;
			    parameterized_quadrature_sites(4,1) = 4.54536892697893e-01;
			    parameterized_quadrature_sites(5,0) = 4.54536892697893e-01;
			    parameterized_quadrature_sites(5,1) = 9.09262146042150e-02;
			    parameterized_quadrature_sites(6,0) = 9.09262146042150e-02;
			    parameterized_quadrature_sites(6,1) = 4.54536892697893e-01;
			    parameterized_quadrature_sites(7,0) = 4.01416680649431e-01;
			    parameterized_quadrature_sites(7,1) = 4.01416680649431e-01;
			    parameterized_quadrature_sites(8,0) = 4.01416680649431e-01;
			    parameterized_quadrature_sites(8,1) = 1.97166638701138e-01;
			    parameterized_quadrature_sites(9,0) = 1.97166638701138e-01;
			    parameterized_quadrature_sites(9,1) = 4.01416680649431e-01;
			    parameterized_quadrature_sites(10,0) = 2.55551654403098e-01;
			    parameterized_quadrature_sites(10,1) = 2.55551654403098e-01;
			    parameterized_quadrature_sites(11,0) = 2.55551654403098e-01;
			    parameterized_quadrature_sites(11,1) = 4.88896691193805e-01;
			    parameterized_quadrature_sites(12,0) = 4.88896691193805e-01;
			    parameterized_quadrature_sites(12,1) = 2.55551654403098e-01;
			    parameterized_quadrature_sites(13,0) = 1.77077942152130e-01;
			    parameterized_quadrature_sites(13,1) = 1.77077942152130e-01;
			    parameterized_quadrature_sites(14,0) = 1.77077942152130e-01;
			    parameterized_quadrature_sites(14,1) = 6.45844115695741e-01;
			    parameterized_quadrature_sites(15,0) = 6.45844115695741e-01;
			    parameterized_quadrature_sites(15,1) = 1.77077942152130e-01;
			    parameterized_quadrature_sites(16,0) = 1.10061053227952e-01;
			    parameterized_quadrature_sites(16,1) = 1.10061053227952e-01;
			    parameterized_quadrature_sites(17,0) = 1.10061053227952e-01;
			    parameterized_quadrature_sites(17,1) = 7.79877893544096e-01;
			    parameterized_quadrature_sites(18,0) = 7.79877893544096e-01;
			    parameterized_quadrature_sites(18,1) = 1.10061053227952e-01;
			    parameterized_quadrature_sites(19,0) = 5.55286242518400e-02;
			    parameterized_quadrature_sites(19,1) = 5.55286242518400e-02;
			    parameterized_quadrature_sites(20,0) = 5.55286242518400e-02;
			    parameterized_quadrature_sites(20,1) = 8.88942751496321e-01;
			    parameterized_quadrature_sites(21,0) = 8.88942751496321e-01;
			    parameterized_quadrature_sites(21,1) = 5.55286242518400e-02;
			    parameterized_quadrature_sites(22,0) = 1.26218637772290e-02;
			    parameterized_quadrature_sites(22,1) = 1.26218637772290e-02;
			    parameterized_quadrature_sites(23,0) = 1.26218637772290e-02;
			    parameterized_quadrature_sites(23,1) = 9.74756272445543e-01;
			    parameterized_quadrature_sites(24,0) = 9.74756272445543e-01;
			    parameterized_quadrature_sites(24,1) = 1.26218637772290e-02;
			    parameterized_quadrature_sites(25,0) = 3.95754787356943e-01;
			    parameterized_quadrature_sites(25,1) = 6.00633794794645e-01;
			    parameterized_quadrature_sites(26,0) = 6.00633794794645e-01;
			    parameterized_quadrature_sites(26,1) = 3.61141784841201e-03;
			    parameterized_quadrature_sites(27,0) = 3.61141784841201e-03;
			    parameterized_quadrature_sites(27,1) = 3.95754787356943e-01;
			    parameterized_quadrature_sites(28,0) = 3.95754787356943e-01;
			    parameterized_quadrature_sites(28,1) = 3.61141784841201e-03;
			    parameterized_quadrature_sites(29,0) = 6.00633794794645e-01;
			    parameterized_quadrature_sites(29,1) = 3.95754787356943e-01;
			    parameterized_quadrature_sites(30,0) = 3.61141784841201e-03;
			    parameterized_quadrature_sites(30,1) = 6.00633794794643e-01;
			    parameterized_quadrature_sites(31,0) = 3.07929983880436e-01;
			    parameterized_quadrature_sites(31,1) = 5.57603261588784e-01;
			    parameterized_quadrature_sites(32,0) = 5.57603261588784e-01;
			    parameterized_quadrature_sites(32,1) = 1.34466754530780e-01;
			    parameterized_quadrature_sites(33,0) = 1.34466754530780e-01;
			    parameterized_quadrature_sites(33,1) = 3.07929983880436e-01;
			    parameterized_quadrature_sites(34,0) = 3.07929983880436e-01;
			    parameterized_quadrature_sites(34,1) = 1.34466754530780e-01;
			    parameterized_quadrature_sites(35,0) = 5.57603261588784e-01;
			    parameterized_quadrature_sites(35,1) = 3.07929983880436e-01;
			    parameterized_quadrature_sites(36,0) = 1.34466754530780e-01;
			    parameterized_quadrature_sites(36,1) = 5.57603261588784e-01;
			    parameterized_quadrature_sites(37,0) = 2.64566948406520e-01;
			    parameterized_quadrature_sites(37,1) = 7.20987025817365e-01;
			    parameterized_quadrature_sites(38,0) = 7.20987025817365e-01;
			    parameterized_quadrature_sites(38,1) = 1.44460257761150e-02;
			    parameterized_quadrature_sites(39,0) = 1.44460257761150e-02;
			    parameterized_quadrature_sites(39,1) = 2.64566948406520e-01;
			    parameterized_quadrature_sites(40,0) = 2.64566948406520e-01;
			    parameterized_quadrature_sites(40,1) = 1.44460257761150e-02;
			    parameterized_quadrature_sites(41,0) = 7.20987025817365e-01;
			    parameterized_quadrature_sites(41,1) = 2.64566948406520e-01;
			    parameterized_quadrature_sites(42,0) = 1.44460257761150e-02;
			    parameterized_quadrature_sites(42,1) = 7.20987025817365e-01;
			    parameterized_quadrature_sites(43,0) = 3.58539352205951e-01;
			    parameterized_quadrature_sites(43,1) = 5.94527068955871e-01;
			    parameterized_quadrature_sites(44,0) = 5.94527068955871e-01;
			    parameterized_quadrature_sites(44,1) = 4.69335788381780e-02;
			    parameterized_quadrature_sites(45,0) = 4.69335788381780e-02;
			    parameterized_quadrature_sites(45,1) = 3.58539352205951e-01;
			    parameterized_quadrature_sites(46,0) = 3.58539352205951e-01;
			    parameterized_quadrature_sites(46,1) = 4.69335788381780e-02;
			    parameterized_quadrature_sites(47,0) = 5.94527068955871e-01;
			    parameterized_quadrature_sites(47,1) = 3.58539352205951e-01;
			    parameterized_quadrature_sites(48,0) = 4.69335788381780e-02;
			    parameterized_quadrature_sites(48,1) = 5.94527068955871e-01;
			    parameterized_quadrature_sites(49,0) = 1.57807405968595e-01;
			    parameterized_quadrature_sites(49,1) = 8.39331473680839e-01;
			    parameterized_quadrature_sites(50,0) = 8.39331473680839e-01;
			    parameterized_quadrature_sites(50,1) = 2.86112035056701e-03;
			    parameterized_quadrature_sites(51,0) = 2.86112035056701e-03;
			    parameterized_quadrature_sites(51,1) = 1.57807405968595e-01;
			    parameterized_quadrature_sites(52,0) = 1.57807405968595e-01;
			    parameterized_quadrature_sites(52,1) = 2.86112035056701e-03;
			    parameterized_quadrature_sites(53,0) = 8.39331473680839e-01;
			    parameterized_quadrature_sites(53,1) = 1.57807405968595e-01;
			    parameterized_quadrature_sites(54,0) = 2.86112035056701e-03;
			    parameterized_quadrature_sites(54,1) = 8.39331473680839e-01;
			    parameterized_quadrature_sites(55,0) = 7.50505969759110e-02;
			    parameterized_quadrature_sites(55,1) = 7.01087978926173e-01;
			    parameterized_quadrature_sites(56,0) = 7.01087978926173e-01;
			    parameterized_quadrature_sites(56,1) = 2.23861424097916e-01;
			    parameterized_quadrature_sites(57,0) = 2.23861424097916e-01;
			    parameterized_quadrature_sites(57,1) = 7.50505969759110e-02;
			    parameterized_quadrature_sites(58,0) = 7.50505969759110e-02;
			    parameterized_quadrature_sites(58,1) = 2.23861424097916e-01;
			    parameterized_quadrature_sites(59,0) = 7.01087978926173e-01;
			    parameterized_quadrature_sites(59,1) = 7.50505969759110e-02;
			    parameterized_quadrature_sites(60,0) = 2.23861424097916e-01;
			    parameterized_quadrature_sites(60,1) = 7.01087978926173e-01;
			    parameterized_quadrature_sites(61,0) = 1.42421601113383e-01;
			    parameterized_quadrature_sites(61,1) = 8.22931324069857e-01;
			    parameterized_quadrature_sites(62,0) = 8.22931324069857e-01;
			    parameterized_quadrature_sites(62,1) = 3.46470748167600e-02;
			    parameterized_quadrature_sites(63,0) = 3.46470748167600e-02;
			    parameterized_quadrature_sites(63,1) = 1.42421601113383e-01;
			    parameterized_quadrature_sites(64,0) = 1.42421601113383e-01;
			    parameterized_quadrature_sites(64,1) = 3.46470748167600e-02;
			    parameterized_quadrature_sites(65,0) = 8.22931324069857e-01;
			    parameterized_quadrature_sites(65,1) = 1.42421601113383e-01;
			    parameterized_quadrature_sites(66,0) = 3.46470748167600e-02;
			    parameterized_quadrature_sites(66,1) = 8.22931324069857e-01;
			    parameterized_quadrature_sites(67,0) = 6.54946280829380e-02;
			    parameterized_quadrature_sites(67,1) = 9.24344252620784e-01;
			    parameterized_quadrature_sites(68,0) = 9.24344252620784e-01;
			    parameterized_quadrature_sites(68,1) = 1.01611192962780e-02;
			    parameterized_quadrature_sites(69,0) = 1.01611192962780e-02;
			    parameterized_quadrature_sites(69,1) = 6.54946280829380e-02;
			    parameterized_quadrature_sites(70,0) = 6.54946280829380e-02;
			    parameterized_quadrature_sites(70,1) = 1.01611192962780e-02;
			    parameterized_quadrature_sites(71,0) = 9.24344252620784e-01;
			    parameterized_quadrature_sites(71,1) = 6.54946280829380e-02;
			    parameterized_quadrature_sites(72,0) = 1.01611192962780e-02;
			    parameterized_quadrature_sites(72,1) = 9.24344252620784e-01;
			    quadrature_weights(0) = 1.64531656944595e-02;
			    quadrature_weights(1) = 5.16536594563600e-03;
			    quadrature_weights(2) = 5.16536594563600e-03;
			    quadrature_weights(3) = 5.16536594563600e-03;
			    quadrature_weights(4) = 1.11936236315080e-02;
			    quadrature_weights(5) = 1.11936236315080e-02;
			    quadrature_weights(6) = 1.11936236315080e-02;
			    quadrature_weights(7) = 1.51330629347340e-02;
			    quadrature_weights(8) = 1.51330629347340e-02;
			    quadrature_weights(9) = 1.51330629347340e-02;
			    quadrature_weights(10) = 1.52454839010990e-02;
			    quadrature_weights(11) = 1.52454839010990e-02;
			    quadrature_weights(12) = 1.52454839010990e-02;
			    quadrature_weights(13) = 1.20796063708205e-02;
			    quadrature_weights(14) = 1.20796063708205e-02;
			    quadrature_weights(15) = 1.20796063708205e-02;
			    quadrature_weights(16) = 8.02540179340050e-03;
			    quadrature_weights(17) = 8.02540179340050e-03;
			    quadrature_weights(18) = 8.02540179340050e-03;
			    quadrature_weights(19) = 4.04229013089200e-03;
			    quadrature_weights(20) = 4.04229013089200e-03;
			    quadrature_weights(21) = 4.04229013089200e-03;
			    quadrature_weights(22) = 1.03968101374250e-03;
			    quadrature_weights(23) = 1.03968101374250e-03;
			    quadrature_weights(24) = 1.03968101374250e-03;
			    quadrature_weights(25) = 1.94243845249050e-03;
			    quadrature_weights(26) = 1.94243845249050e-03;
			    quadrature_weights(27) = 1.94243845249050e-03;
			    quadrature_weights(28) = 1.94243845249050e-03;
			    quadrature_weights(29) = 1.94243845249050e-03;
			    quadrature_weights(30) = 1.94243845249050e-03;
			    quadrature_weights(31) = 1.27870803060110e-02;
			    quadrature_weights(32) = 1.27870803060110e-02;
			    quadrature_weights(33) = 1.27870803060110e-02;
			    quadrature_weights(34) = 1.27870803060110e-02;
			    quadrature_weights(35) = 1.27870803060110e-02;
			    quadrature_weights(36) = 1.27870803060110e-02;
			    quadrature_weights(37) = 4.44045178666900e-03;
			    quadrature_weights(38) = 4.44045178666900e-03;
			    quadrature_weights(39) = 4.44045178666900e-03;
			    quadrature_weights(40) = 4.44045178666900e-03;
			    quadrature_weights(41) = 4.44045178666900e-03;
			    quadrature_weights(42) = 4.44045178666900e-03;
			    quadrature_weights(43) = 8.06227338086550e-03;
			    quadrature_weights(44) = 8.06227338086550e-03;
			    quadrature_weights(45) = 8.06227338086550e-03;
			    quadrature_weights(46) = 8.06227338086550e-03;
			    quadrature_weights(47) = 8.06227338086550e-03;
			    quadrature_weights(48) = 8.06227338086550e-03;
			    quadrature_weights(49) = 1.24597090874550e-03;
			    quadrature_weights(50) = 1.24597090874550e-03;
			    quadrature_weights(51) = 1.24597090874550e-03;
			    quadrature_weights(52) = 1.24597090874550e-03;
			    quadrature_weights(53) = 1.24597090874550e-03;
			    quadrature_weights(54) = 1.24597090874550e-03;
			    quadrature_weights(55) = 9.12142005947550e-03;
			    quadrature_weights(56) = 9.12142005947550e-03;
			    quadrature_weights(57) = 9.12142005947550e-03;
			    quadrature_weights(58) = 9.12142005947550e-03;
			    quadrature_weights(59) = 9.12142005947550e-03;
			    quadrature_weights(60) = 9.12142005947550e-03;
			    quadrature_weights(61) = 5.12928186809950e-03;
			    quadrature_weights(62) = 5.12928186809950e-03;
			    quadrature_weights(63) = 5.12928186809950e-03;
			    quadrature_weights(64) = 5.12928186809950e-03;
			    quadrature_weights(65) = 5.12928186809950e-03;
			    quadrature_weights(66) = 5.12928186809950e-03;
			    quadrature_weights(67) = 1.89996442765100e-03;
			    quadrature_weights(68) = 1.89996442765100e-03;
			    quadrature_weights(69) = 1.89996442765100e-03;
			    quadrature_weights(70) = 1.89996442765100e-03;
			    quadrature_weights(71) = 1.89996442765100e-03;
			    quadrature_weights(72) = 1.89996442765100e-03;
			    break;
			case 78:
			    parameterized_quadrature_sites(0,0) = 8.9411337112035999e-03;
			    parameterized_quadrature_sites(0,1) = 8.6983293701984998e-03;
			    parameterized_quadrature_sites(1,0) = 9.7926226298067365e-01;
			    parameterized_quadrature_sites(1,1) = 1.0264413374365100e-02;
			    parameterized_quadrature_sites(2,0) = 1.0547538211187800e-02;
			    parameterized_quadrature_sites(2,1) = 9.7855142025151109e-01;
			    parameterized_quadrature_sites(3,0) = 2.3777061947122002e-03;
			    parameterized_quadrature_sites(3,1) = 6.3655109860361700e-02;
			    parameterized_quadrature_sites(4,0) = 6.3042511579465998e-02;
			    parameterized_quadrature_sites(4,1) = 4.1506347508631003e-03;
			    parameterized_quadrature_sites(5,0) = 9.3084224967299967e-01;
			    parameterized_quadrature_sites(5,1) = 4.8053482262546002e-03;
			    parameterized_quadrature_sites(6,0) = 6.2907655549027400e-02;
			    parameterized_quadrature_sites(6,1) = 9.3167900694812233e-01;
			    parameterized_quadrature_sites(7,0) = 9.3159622463806491e-01;
			    parameterized_quadrature_sites(7,1) = 6.2626488180135900e-02;
			    parameterized_quadrature_sites(8,0) = 6.1951689414552003e-03;
			    parameterized_quadrature_sites(8,1) = 9.2935870585640645e-01;
			    parameterized_quadrature_sites(9,0) = 2.8712581923668101e-02;
			    parameterized_quadrature_sites(9,1) = 3.1020212299716299e-02;
			    parameterized_quadrature_sites(10,0) = 9.2938444783052321e-01;
			    parameterized_quadrature_sites(10,1) = 3.4215296821852897e-02;
			    parameterized_quadrature_sites(11,0) = 3.7545756662128102e-02;
			    parameterized_quadrature_sites(11,1) = 9.2578688846693047e-01;
			    parameterized_quadrature_sites(12,0) = 8.6895739063833997e-03;
			    parameterized_quadrature_sites(12,1) = 1.5849712515099221e-01;
			    parameterized_quadrature_sites(13,0) = 1.5475970539646791e-01;
			    parameterized_quadrature_sites(13,1) = 8.3636066576882862e-01;
			    parameterized_quadrature_sites(14,0) = 8.3310252941849239e-01;
			    parameterized_quadrature_sites(14,1) = 8.9257244824476004e-03;
			    parameterized_quadrature_sites(15,0) = 8.3742310735260950e-01;
			    parameterized_quadrature_sites(15,1) = 1.5291673040783921e-01;
			    parameterized_quadrature_sites(16,0) = 1.5593625052337881e-01;
			    parameterized_quadrature_sites(16,1) = 9.4966240058029002e-03;
			    parameterized_quadrature_sites(17,0) = 9.8599642095236004e-03;
			    parameterized_quadrature_sites(17,1) = 8.3422114935955050e-01;
			    parameterized_quadrature_sites(18,0) = 4.0558737332891631e-01;
			    parameterized_quadrature_sites(18,1) = 7.4389302007913001e-03;
			    parameterized_quadrature_sites(19,0) = 5.9647278986182350e-01;
			    parameterized_quadrature_sites(19,1) = 3.9563308093107152e-01;
			    parameterized_quadrature_sites(20,0) = 8.0747800415767006e-03;
			    parameterized_quadrature_sites(20,1) = 4.0313194259026802e-01;
			    parameterized_quadrature_sites(21,0) = 7.5073977720710996e-03;
			    parameterized_quadrature_sites(21,1) = 5.8516095946805691e-01;
			    parameterized_quadrature_sites(22,0) = 3.9367645192372991e-01;
			    parameterized_quadrature_sites(22,1) = 5.9748965928987985e-01;
			    parameterized_quadrature_sites(23,0) = 5.8465307262122179e-01;
			    parameterized_quadrature_sites(23,1) = 8.7250464968192006e-03;
			    parameterized_quadrature_sites(24,0) = 4.8708041121196383e-01;
			    parameterized_quadrature_sites(24,1) = 2.0212922991194000e-02;
			    parameterized_quadrature_sites(25,0) = 2.6835128117845169e-01;
			    parameterized_quadrature_sites(25,1) = 7.2023400886682198e-01;
			    parameterized_quadrature_sites(26,0) = 7.2239562887479880e-01;
			    parameterized_quadrature_sites(26,1) = 2.6623993664561901e-01;
			    parameterized_quadrature_sites(27,0) = 2.7168267423572212e-01;
			    parameterized_quadrature_sites(27,1) = 1.1288269880823600e-02;
			    parameterized_quadrature_sites(28,0) = 1.1258084204589300e-02;
			    parameterized_quadrature_sites(28,1) = 7.1696959633251023e-01;
			    parameterized_quadrature_sites(29,0) = 1.1503473436974001e-02;
			    parameterized_quadrature_sites(29,1) = 2.7400671101656832e-01;
			    parameterized_quadrature_sites(30,0) = 7.1405259005638033e-01;
			    parameterized_quadrature_sites(30,1) = 1.1351156049706200e-02;
			    parameterized_quadrature_sites(31,0) = 4.9028710531115449e-01;
			    parameterized_quadrature_sites(31,1) = 4.9364918414683351e-01;
			    parameterized_quadrature_sites(32,0) = 2.0142342520930698e-02;
			    parameterized_quadrature_sites(32,1) = 4.8325734596013992e-01;
			    parameterized_quadrature_sites(33,0) = 3.6110746485855001e-02;
			    parameterized_quadrature_sites(33,1) = 9.3567950158201393e-02;
			    parameterized_quadrature_sites(34,0) = 8.6079988198508572e-01;
			    parameterized_quadrature_sites(34,1) = 3.9737906707539197e-02;
			    parameterized_quadrature_sites(35,0) = 1.0058915260013050e-01;
			    parameterized_quadrature_sites(35,1) = 8.5863434193517962e-01;
			    parameterized_quadrature_sites(36,0) = 9.1874071705841595e-02;
			    parameterized_quadrature_sites(36,1) = 3.9551300197337699e-02;
			    parameterized_quadrature_sites(37,0) = 8.6048882961910289e-01;
			    parameterized_quadrature_sites(37,1) = 9.6622405707924700e-02;
			    parameterized_quadrature_sites(38,0) = 4.3984217867325599e-02;
			    parameterized_quadrature_sites(38,1) = 8.5618863491067676e-01;
			    parameterized_quadrature_sites(39,0) = 2.0110176067354310e-01;
			    parameterized_quadrature_sites(39,1) = 7.4491158356262255e-01;
			    parameterized_quadrature_sites(40,0) = 7.4499937262632787e-01;
			    parameterized_quadrature_sites(40,1) = 5.3686563816580400e-02;
			    parameterized_quadrature_sites(41,0) = 5.3218664130983202e-02;
			    parameterized_quadrature_sites(41,1) = 1.9637542759350521e-01;
			    parameterized_quadrature_sites(42,0) = 7.4539846474005178e-01;
			    parameterized_quadrature_sites(42,1) = 1.9820658055500051e-01;
			    parameterized_quadrature_sites(43,0) = 1.9572899328760179e-01;
			    parameterized_quadrature_sites(43,1) = 5.5571383315608597e-02;
			    parameterized_quadrature_sites(44,0) = 1.0925320579875419e-01;
			    parameterized_quadrature_sites(44,1) = 6.1000361824130300e-01;
			    parameterized_quadrature_sites(45,0) = 5.6762570200051501e-02;
			    parameterized_quadrature_sites(45,1) = 7.4091218949591942e-01;
			    parameterized_quadrature_sites(46,0) = 4.8383793347481101e-02;
			    parameterized_quadrature_sites(46,1) = 6.0751356609779783e-01;
			    parameterized_quadrature_sites(47,0) = 1.0806128097601329e-01;
			    parameterized_quadrature_sites(47,1) = 1.1220815104370099e-01;
			    parameterized_quadrature_sites(48,0) = 6.1856059009905007e-01;
			    parameterized_quadrature_sites(48,1) = 2.6987537030349740e-01;
			    parameterized_quadrature_sites(49,0) = 7.7212960134965625e-01;
			    parameterized_quadrature_sites(49,1) = 1.1141173953329921e-01;
			    parameterized_quadrature_sites(50,0) = 6.1157348011327173e-01;
			    parameterized_quadrature_sites(50,1) = 3.3893676779306348e-01;
			    parameterized_quadrature_sites(51,0) = 3.3813261033758418e-01;
			    parameterized_quadrature_sites(51,1) = 4.9469393878745799e-02;
			    parameterized_quadrature_sites(52,0) = 1.1730841282542900e-01;
			    parameterized_quadrature_sites(52,1) = 7.6964513097951825e-01;
			    parameterized_quadrature_sites(53,0) = 2.6745512605961458e-01;
			    parameterized_quadrature_sites(53,1) = 1.1157188081540730e-01;
			    parameterized_quadrature_sites(54,0) = 6.5421001600256889e-01;
			    parameterized_quadrature_sites(54,1) = 1.9065483146999149e-01;
			    parameterized_quadrature_sites(55,0) = 5.3829748115775802e-02;
			    parameterized_quadrature_sites(55,1) = 3.3586168268491179e-01;
			    parameterized_quadrature_sites(56,0) = 1.8488403241167711e-01;
			    parameterized_quadrature_sites(56,1) = 1.5518315238513730e-01;
			    parameterized_quadrature_sites(57,0) = 3.3762671047443338e-01;
			    parameterized_quadrature_sites(57,1) = 6.0814025962944529e-01;
			    parameterized_quadrature_sites(58,0) = 6.0671020344994708e-01;
			    parameterized_quadrature_sites(58,1) = 5.4263279559821201e-02;
			    parameterized_quadrature_sites(59,0) = 4.6126140854956371e-01;
			    parameterized_quadrature_sites(59,1) = 6.8817667072165398e-02;
			    parameterized_quadrature_sites(60,0) = 1.5254653656712561e-01;
			    parameterized_quadrature_sites(60,1) = 6.5102408457488470e-01;
			    parameterized_quadrature_sites(61,0) = 7.0058254354307500e-02;
			    parameterized_quadrature_sites(61,1) = 4.6619043927415987e-01;
			    parameterized_quadrature_sites(62,0) = 4.7042013790318088e-01;
			    parameterized_quadrature_sites(62,1) = 4.6348264553531421e-01;
			    parameterized_quadrature_sites(63,0) = 1.2164616937459330e-01;
			    parameterized_quadrature_sites(63,1) = 2.3814948755156831e-01;
			    parameterized_quadrature_sites(64,0) = 6.3714040527021165e-01;
			    parameterized_quadrature_sites(64,1) = 1.2383993845133670e-01;
			    parameterized_quadrature_sites(65,0) = 2.3799045151187120e-01;
			    parameterized_quadrature_sites(65,1) = 6.3702164523263760e-01;
			    parameterized_quadrature_sites(66,0) = 1.4839298571771459e-01;
			    parameterized_quadrature_sites(66,1) = 4.8941885777801442e-01;
			    parameterized_quadrature_sites(67,0) = 3.5980695715496319e-01;
			    parameterized_quadrature_sites(67,1) = 1.4528808662532389e-01;
			    parameterized_quadrature_sites(68,0) = 4.9414410550951349e-01;
			    parameterized_quadrature_sites(68,1) = 3.6102163838181101e-01;
			    parameterized_quadrature_sites(69,0) = 1.4406306879808209e-01;
			    parameterized_quadrature_sites(69,1) = 3.5135083418870572e-01;
			    parameterized_quadrature_sites(70,0) = 5.0197644400035468e-01;
			    parameterized_quadrature_sites(70,1) = 1.4354916632930600e-01;
			    parameterized_quadrature_sites(71,0) = 3.5554238342982608e-01;
			    parameterized_quadrature_sites(71,1) = 5.0164915995018422e-01;
			    parameterized_quadrature_sites(72,0) = 2.4434395407713269e-01;
			    parameterized_quadrature_sites(72,1) = 2.4060521291041001e-01;
			    parameterized_quadrature_sites(73,0) = 2.4370649893418969e-01;
			    parameterized_quadrature_sites(73,1) = 5.1090172770553444e-01;
			    parameterized_quadrature_sites(74,0) = 5.1222008073208247e-01;
			    parameterized_quadrature_sites(74,1) = 2.4527379735428820e-01;
			    parameterized_quadrature_sites(75,0) = 2.5260383151777532e-01;
			    parameterized_quadrature_sites(75,1) = 3.7003195550936951e-01;
			    parameterized_quadrature_sites(76,0) = 3.7598956528506539e-01;
			    parameterized_quadrature_sites(76,1) = 2.5054066116305501e-01;
			    parameterized_quadrature_sites(77,0) = 3.7290779871441049e-01;
			    parameterized_quadrature_sites(77,1) = 3.7537502775491960e-01;
			    quadrature_weights(0) = 5.4361363496487503e-04;
			    quadrature_weights(1) = 7.2467838162779998e-04;
			    quadrature_weights(2) = 7.7115073341352497e-04;
			    quadrature_weights(3) = 8.6004082758987495e-04;
			    quadrature_weights(4) = 1.0474618002941001e-03;
			    quadrature_weights(5) = 1.1184512874480750e-03;
			    quadrature_weights(6) = 1.1763605203547751e-03;
			    quadrature_weights(7) = 1.2216983937579000e-03;
			    quadrature_weights(8) = 1.2981910842235500e-03;
			    quadrature_weights(9) = 1.8518264745199500e-03;
			    quadrature_weights(10) = 1.9938852575230501e-03;
			    quadrature_weights(11) = 2.0887630727400748e-03;
			    quadrature_weights(12) = 2.4041665215934500e-03;
			    quadrature_weights(13) = 2.4079564462609749e-03;
			    quadrature_weights(14) = 2.4644365189580250e-03;
			    quadrature_weights(15) = 2.5664470075342249e-03;
			    quadrature_weights(16) = 2.5797025777773248e-03;
			    quadrature_weights(17) = 2.6572750407618498e-03;
			    quadrature_weights(18) = 2.6720326723837502e-03;
			    quadrature_weights(19) = 2.6742255252576500e-03;
			    quadrature_weights(20) = 2.7256615428406252e-03;
			    quadrature_weights(21) = 2.7474945893627501e-03;
			    quadrature_weights(22) = 2.8355763807253751e-03;
			    quadrature_weights(23) = 3.0133910732457001e-03;
			    quadrature_weights(24) = 3.4904798455149752e-03;
			    quadrature_weights(25) = 3.5286997884102249e-03;
			    quadrature_weights(26) = 3.5482586761459252e-03;
			    quadrature_weights(27) = 3.6053169067081750e-03;
			    quadrature_weights(28) = 3.6176086713727751e-03;
			    quadrature_weights(29) = 3.6237442468015752e-03;
			    quadrature_weights(30) = 3.6346693923518499e-03;
			    quadrature_weights(31) = 3.6491047731400249e-03;
			    quadrature_weights(32) = 3.6828644616571501e-03;
			    quadrature_weights(33) = 4.1865990826037754e-03;
			    quadrature_weights(34) = 4.2238875114384497e-03;
			    quadrature_weights(35) = 4.2355665721040998e-03;
			    quadrature_weights(36) = 4.3267543023872996e-03;
			    quadrature_weights(37) = 4.3631136623173503e-03;
			    quadrature_weights(38) = 4.4304305539834754e-03;
			    quadrature_weights(39) = 7.0706006005818254e-03;
			    quadrature_weights(40) = 7.1249178122030250e-03;
			    quadrature_weights(41) = 7.1251411634678250e-03;
			    quadrature_weights(42) = 7.5161805869424503e-03;
			    quadrature_weights(43) = 7.5507819270567997e-03;
			    quadrature_weights(44) = 7.5996784019234003e-03;
			    quadrature_weights(45) = 7.6417199018501249e-03;
			    quadrature_weights(46) = 7.6516853250572746e-03;
			    quadrature_weights(47) = 7.7332517050150247e-03;
			    quadrature_weights(48) = 7.7443455208702503e-03;
			    quadrature_weights(49) = 7.8286562636164757e-03;
			    quadrature_weights(50) = 7.8393373347939495e-03;
			    quadrature_weights(51) = 7.8580117321664505e-03;
			    quadrature_weights(52) = 7.8795535973380255e-03;
			    quadrature_weights(53) = 8.1062034496375505e-03;
			    quadrature_weights(54) = 8.6878038096496243e-03;
			    quadrature_weights(55) = 8.7598363731706003e-03;
			    quadrature_weights(56) = 8.7679355077581993e-03;
			    quadrature_weights(57) = 8.8032303833538994e-03;
			    quadrature_weights(58) = 8.8153876245372752e-03;
			    quadrature_weights(59) = 9.1600805085676005e-03;
			    quadrature_weights(60) = 9.1933276917576506e-03;
			    quadrature_weights(61) = 9.2918915734212758e-03;
			    quadrature_weights(62) = 9.3342892901571994e-03;
			    quadrature_weights(63) = 1.0099333664696650e-02;
			    quadrature_weights(64) = 1.0339501015941751e-02;
			    quadrature_weights(65) = 1.0548944796751400e-02;
			    quadrature_weights(66) = 1.2386275100928551e-02;
			    quadrature_weights(67) = 1.2510481528527850e-02;
			    quadrature_weights(68) = 1.2644864677864924e-02;
			    quadrature_weights(69) = 1.3000930254698475e-02;
			    quadrature_weights(70) = 1.3038339197138024e-02;
			    quadrature_weights(71) = 1.3122478808947751e-02;
			    quadrature_weights(72) = 1.4978994062894975e-02;
			    quadrature_weights(73) = 1.4990249935657125e-02;
			    quadrature_weights(74) = 1.4997881803219050e-02;
			    quadrature_weights(75) = 1.5853329586233849e-02;
			    quadrature_weights(76) = 1.5882796527709402e-02;
			    quadrature_weights(77) = 1.5930165141800151e-02;
			    break;
            default:
                compadre_assert_release(false && "Number of quadrature points not supported.");
        }
    } else {
        compadre_assert_release(false && "Dimension not supported for quadrature.");
    }

//  {
//    7,
//    {{3.3333333333333333333333333333333333333333e-1,  3.3333333333333333333333333333333333333333e-1, 0.0},
//     {1.0128650732345633880098736191512382805558e-1,  1.0128650732345633880098736191512382805558e-1, 0.0},
//     {7.9742698535308732239802527616975234388885e-1,  1.0128650732345633880098736191512382805558e-1, 0.0},
//     {1.0128650732345633880098736191512382805558e-1,  7.9742698535308732239802527616975234388885e-1, 0.0},
//     {4.7014206410511508977044120951344760051585e-1,  4.7014206410511508977044120951344760051585e-1, 0.0},
//     {5.9715871789769820459117580973104798968293e-2,  4.7014206410511508977044120951344760051585e-1, 0.0},
//     {4.7014206410511508977044120951344760051585e-1,  5.9715871789769820459117580973104798968293e-2, 0.0}},
//    {1.1250000000000000000000000000000000000000e-1,
//     6.2969590272413576297841972750090666828820e-2,
//     6.2969590272413576297841972750090666828820e-2,
//     6.2969590272413576297841972750090666828820e-2,
//     6.6197076394253090368824693916575999837847e-2,
//     6.6197076394253090368824693916575999837847e-2,
//     6.6197076394253090368824693916575999837847e-2}
//  } // order 5
//  {
//    12,
//    {{6.3089014491502228340331602870819157341003e-2,  6.3089014491502228340331602870819157341003e-2, 0.0},
//     {6.3089014491502228340331602870819157341003e-2,  8.7382197101699554331933679425836168531799e-1, 0.0},
//     {8.7382197101699554331933679425836168531799e-1,  6.3089014491502228340331602870819157341003e-2, 0.0},
//     {2.4928674517091042129163855310701907608796e-1,  2.4928674517091042129163855310701907608796e-1, 0.0},
//     {2.4928674517091042129163855310701907608796e-1,  5.0142650965817915741672289378596184782407e-1, 0.0},
//     {5.0142650965817915741672289378596184782407e-1,  2.4928674517091042129163855310701907608796e-1, 0.0},
//     {3.1035245103378441286759566723265117140082e-1,  5.3145049844816939902261738355299128796183e-2, 0.0},
//     {6.3650249912139865667831578657006182134850e-1,  3.1035245103378439596843454179854003165774e-1, 0.0},
//     {5.3145049844816930454088546197287007250682e-2,  6.3650249912139866412930371984616083954608e-1, 0.0},
//     {5.3145049844816939902261738355299128796183e-2,  3.1035245103378441286759566723265117140082e-1, 0.0},
//     {6.3650249912139866412930371984616083954608e-1,  5.3145049844816930454088546197287007250682e-2, 0.0},
//     {3.1035245103378439596843454179854003165774e-1,  6.3650249912139865667831578657006182134850e-1, 0.0}},
//    {2.5422453185103408460468404553434492023395e-2,
//     2.5422453185103408460468404553434492023395e-2,
//     2.5422453185103408460468404553434492023395e-2,
//     5.8393137863189683012644805692789720663043e-2,
//     5.8393137863189683012644805692789720663043e-2,
//     5.8393137863189683012644805692789720663043e-2,
//     4.1425537809186787596776728210221226990114e-2,
//     4.1425537809186787596776728210221226990114e-2,
//     4.1425537809186787596776728210221226990114e-2,
//     4.1425537809186787596776728210221226990114e-2,
//     4.1425537809186787596776728210221226990114e-2,
//     4.1425537809186787596776728210221226990114e-2}
//  }, // order 6
//  {
//    13,
//    {{3.33333333333333e-01, 3.33333333333333e-01, 0.0},
//     {2.60345966079040e-01, 2.60345966079040e-01, 0.0},
//     {2.60345966079040e-01, 4.79308067841920e-01, 0.0},
//     {4.79308067841920e-01, 2.60345966079040e-01, 0.0},
//     {6.51301029022160e-02, 6.51301029022160e-02, 0.0},
//     {6.51301029022160e-02, 8.69739794195568e-01, 0.0},
//     {8.69739794195568e-01, 6.51301029022160e-02, 0.0},
//     {3.12865496004874e-01, 6.38444188569810e-01, 0.0},
//     {6.38444188569810e-01, 4.86903154253160e-02, 0.0},
//     {4.86903154253160e-02, 3.12865496004874e-01, 0.0},
//     {3.12865496004874e-01, 4.86903154253160e-02, 0.0},
//     {6.38444188569810e-01, 3.12865496004874e-01, 0.0},
//     {4.86903154253160e-02, 6.38444188569810e-01, 0.0}},
//    {-7.47850222338410e-02,
//     8.78076287166040e-02,
//     8.78076287166040e-02,
//     8.78076287166040e-02,
//     2.66736178044190e-02,
//     2.66736178044190e-02,
//     2.66736178044190e-02,
//     3.85568804451285e-02,
//     3.85568804451285e-02,
//     3.85568804451285e-02,
//     3.85568804451285e-02,
//     3.85568804451285e-02,
//     3.85568804451285e-02}
//  }, // order 7
//{
//    16,
//    {{3.33333333333333e-01, 3.33333333333333e-01, 0.0},
//     {4.59292588292723e-01, 4.59292588292723e-01, 0.0},
//     {4.59292588292723e-01, 8.14148234145540e-02, 0.0},
//     {8.14148234145540e-02, 4.59292588292723e-01, 0.0},
//     {1.70569307751760e-01, 1.70569307751760e-01, 0.0},
//     {1.70569307751760e-01, 6.58861384496480e-01, 0.0},
//     {6.58861384496480e-01, 1.70569307751760e-01, 0.0},
//     {5.05472283170310e-02, 5.05472283170310e-02, 0.0},
//     {5.05472283170310e-02, 8.98905543365938e-01, 0.0},
//     {8.98905543365938e-01, 5.05472283170310e-02, 0.0},
//     {2.63112829634638e-01, 7.28492392955404e-01, 0.0},
//     {7.28492392955404e-01, 8.39477740995798e-03, 0.0},
//     {8.39477740995798e-03, 2.63112829634638e-01, 0.0},
//     {2.63112829634638e-01, 8.39477740995798e-03, 0.0},
//     {7.28492392955404e-01, 2.63112829634638e-01, 0.0},
//     {8.39477740995798e-03, 7.28492392955404e-01, 0.0}},
//    {7.21578038388935e-02,
//     4.75458171336425e-02,
//     4.75458171336425e-02,
//     4.75458171336425e-02,
//     5.16086852673590e-02,
//     5.16086852673590e-02,
//     5.16086852673590e-02,
//     1.62292488115990e-02,
//     1.62292488115988e-02,
//     1.62292488115990e-02,
//     1.36151570872175e-02,
//     1.36151570872175e-02,
//     1.36151570872175e-02,
//     1.36151570872175e-02,
//     1.36151570872175e-02,
//     1.36151570872175e-02}
//  }, // order 8
//  {
//    19,
//    {{3.33333333333333e-01, 3.33333333333333e-01, 0.0},
//     {4.89682519198738e-01, 4.89682519198738e-01, 0.0},
//     {4.89682519198738e-01, 2.06349616025250e-02, 0.0},
//     {2.06349616025250e-02, 4.89682519198738e-01, 0.0},
//     {4.37089591492937e-01, 4.37089591492937e-01, 0.0},
//     {4.37089591492937e-01, 1.25820817014127e-01, 0.0},
//     {1.25820817014127e-01, 4.37089591492937e-01, 0.0},
//     {1.88203535619033e-01, 1.88203535619033e-01, 0.0},
//     {1.88203535619033e-01, 6.23592928761935e-01, 0.0},
//     {6.23592928761935e-01, 1.88203535619033e-01, 0.0},
//     {4.47295133944530e-02, 4.47295133944530e-02, 0.0},
//     {4.47295133944530e-02, 9.10540973211095e-01, 0.0},
//     {9.10540973211095e-01, 4.47295133944530e-02, 0.0},
//     {2.21962989160766e-01, 7.41198598784498e-01, 0.0},
//     {7.41198598784498e-01, 3.68384120547360e-02, 0.0},
//     {3.68384120547360e-02, 2.21962989160766e-01, 0.0},
//     {2.21962989160766e-01, 3.68384120547360e-02, 0.0},
//     {7.41198598784498e-01, 2.21962989160766e-01, 0.0},
//     {3.68384120547360e-02, 7.41198598784498e-01, 0.0}},
//    {4.85678981413995e-02,
//     1.56673501135695e-02,
//     1.56673501135695e-02,
//     1.56673501135695e-02,
//     3.89137705023870e-02,
//     3.89137705023870e-02,
//     3.89137705023870e-02,
//     3.98238694636050e-02,
//     3.98238694636050e-02,
//     3.98238694636050e-02,
//     1.27888378293490e-02,
//     1.27888378293490e-02,
//     1.27888378293490e-02,
//     2.16417696886445e-02,
//     2.16417696886445e-02,
//     2.16417696886445e-02,
//     2.16417696886445e-02,
//     2.16417696886445e-02,
//     2.16417696886445e-02}
//  }, // order 9
//{
//    24,
//    {{3.8102570854643002e-03, 8.6854386943076545e-01, 0.0},
//     {8.3865349500109043e-01, 1.6134650499890960e-01, 0.0},
//     {0.0                   , 3.9366774470722010e-01, 0.0},
//     {7.7757518429429107e-01, 0.0                   , 0.0},
//     {4.7768381772022403e-02, 9.2899486985787905e-01, 0.0},
//     {1.0939142057119900e-02, 1.7690730625559031e-01, 0.0},
//     {4.6374383867430541e-01, 0.0                   , 0.0},
//     {9.3049846900263089e-01, 2.9553592846822900e-02, 0.0},
//     {3.9099745550423302e-02, 3.5319656252586103e-02, 0.0},
//     {4.8798437805397499e-01, 5.0365825075943971e-01, 0.0},
//     {1.9305903224251941e-01, 3.0573404093099301e-02, 0.0},
//     {2.2376358774275851e-01, 7.4726591728868819e-01, 0.0},
//     {3.6036266787907702e-02, 6.3491832379200652e-01, 0.0},
//     {7.6777680170023954e-01, 1.0614642990290001e-01, 0.0},
//     {1.0954959855585469e-01, 7.5329402776254240e-01, 0.0},
//     {6.4203365318662664e-01, 2.9530445535851102e-01, 0.0},
//     {1.0999439055630450e-01, 1.6929927488966459e-01, 0.0},
//     {3.3947290311800549e-01, 9.5379208487721703e-02, 0.0},
//     {8.4198522115543697e-02, 3.8729657913960353e-01, 0.0},
//     {5.7966325105486349e-01, 8.0491894656105595e-02, 0.0},
//     {3.6419744430339263e-01, 5.2433682558924433e-01, 0.0},
//     {2.7586334089315973e-01, 2.6481531651496770e-01, 0.0},
//     {2.0776116575484829e-01, 5.0550507373529086e-01, 0.0},
//     {4.8123289062464247e-01, 2.7542385024412980e-01, 0.0}},
//    {5.3333456379563004e-03,
//     5.5334733911235499e-03,
//     5.9501524071177998e-03,
//     6.6995403929284002e-03,
//     7.3041646115608498e-03,
//     7.4233118233402503e-03,
//     7.5781309587828498e-03,
//     7.6015470130063496e-03,
//     8.8279222964811003e-03,
//     9.7797061087609508e-03,
//     1.5707351180334426e-02,
//     1.6151301011728625e-02,
//     2.1481164696039650e-02,
//     2.6268329551205824e-02,
//     2.6365190704419499e-02,
//     2.7457802196927176e-02,
//     2.8987070299349601e-02,
//     3.0871122933097676e-02,
//     3.1343962650456775e-02,
//     3.3093067237629975e-02,
//     3.7039823668389947e-02,
//     4.2207221279855073e-02,
//     4.3362019313832399e-02,
//     4.7633278635674972e-02}
//  }, // order 10
//  {
//    27,
//    {{3.7802163891336921e-01, 6.1948431533135195e-01, 0.0},
//     {3.2899822292186298e-02, 9.3614893514675623e-01, 0.0},
//     {9.3551434285897095e-01, 3.3268560622678398e-02, 0.0},
//     {3.4222771841359197e-02, 3.2916403878999703e-02, 0.0},
//     {1.4354532010930900e-02, 3.9659731669586501e-01, 0.0},
//     {2.2120535196161799e-02, 1.6892970982290231e-01, 0.0},
//     {8.1562969693268217e-01, 2.6807150626772601e-02, 0.0},
//     {2.7719522918618601e-02, 8.1626233715968810e-01, 0.0},
//     {1.7400571673032261e-01, 2.5252704638304500e-02, 0.0},
//     {3.8913981113319362e-01, 2.2592651051306600e-02, 0.0},
//     {8.0364834053903877e-01, 1.6655614492060569e-01, 0.0},
//     {1.6429286715713459e-01, 8.0454974747615537e-01, 0.0},
//     {6.1758873171277151e-01, 2.5660186833052399e-02, 0.0},
//     {2.6297199713764201e-02, 6.1924873232110123e-01, 0.0},
//     {5.9895439629934211e-01, 3.7272769861629101e-01, 0.0},
//     {8.1721404855381805e-02, 3.2719878157552901e-01, 0.0},
//     {1.3035453031942690e-01, 1.3667083534390509e-01, 0.0},
//     {7.1027868107761583e-01, 1.3828000204292321e-01, 0.0},
//     {1.4118119730952799e-01, 7.0099267949645228e-01, 0.0},
//     {5.3141960154079959e-01, 1.2417148586801489e-01, 0.0},
//     {3.4992914334288650e-01, 5.6938486195327997e-01, 0.0},
//     {3.1909737814681871e-01, 1.1698976413323441e-01, 0.0},
//     {1.2454405910544100e-01, 5.1353143433447235e-01, 0.0},
//     {4.1132499178904658e-01, 2.6677168071577739e-01, 0.0},
//     {5.3634228112084714e-01, 3.2081957909482989e-01, 0.0},
//     {2.2789955884347499e-01, 2.8790310224819649e-01, 0.0},
//     {2.9133859436942361e-01, 4.6494564773693992e-01, 0.0}},
//    {5.6375356078552001e-03,
//     6.5529674260441251e-03,
//     6.6871546432313500e-03,
//     7.2740543402810501e-03,
//     9.1794421950682249e-03,
//     1.0688904048956826e-02,
//     1.1556190176810576e-02,
//     1.1821029782882325e-02,
//     1.2035150250304100e-02,
//     1.3361393587366349e-02,
//     1.3378318081626600e-02,
//     1.3797950075089925e-02,
//     1.4360888683763475e-02,
//     1.4717279053216750e-02,
//     1.6298436572453948e-02,
//     2.1298852449057524e-02,
//     2.2769129560075176e-02,
//     2.4595497454187175e-02,
//     2.6030267669060574e-02,
//     2.7060738237627526e-02,
//     2.7485999004421752e-02,
//     2.7495506358711074e-02,
//     2.8322942561349675e-02,
//     2.9305884203704700e-02,
//     3.1096575055627549e-02,
//     3.3071222942480799e-02,
//     3.4120689978745525e-02}
//  },
//  {
//    32,
//    {{9.2734897448394982e-01, 0.0                   , 0.0},
//     {2.3551733249578700e-02, 9.5526919357006035e-01, 0.0},
//     {0.0                   , 8.5815888421533082e-01, 0.0},
//     {9.4547507322097091e-01, 4.3010560106405499e-02, 0.0},
//     {1.5406460162685609e-01, 8.4593539837314391e-01, 0.0},
//     {0.0                   , 6.2731531923241179e-01, 0.0},
//     {2.7110971356255800e-02, 2.9754117496841800e-02, 0.0},
//     {1.4604496167217570e-01, 9.2296909059649008e-03, 0.0},
//     {2.1152223383121900e-02, 1.5557066896897950e-01, 0.0},
//     {1.4566514788347000e-02, 3.6384660446077510e-01, 0.0},
//     {7.8860171922313160e-01, 1.8920633061715941e-01, 0.0},
//     {7.4918973979067949e-01, 2.3088148766115799e-02, 0.0},
//     {7.1871496101589105e-02, 8.5431474947580432e-01, 0.0},
//     {3.3212908394764512e-01, 2.4506286636990001e-02, 0.0},
//     {3.6118159118967208e-01, 6.1600929617267497e-01, 0.0},
//     {2.4345813394879970e-01, 9.3448087604440996e-02, 0.0},
//     {5.8168921474014745e-01, 3.9316510319604808e-01, 0.0},
//     {5.4444667627192522e-01, 2.5716283623693902e-02, 0.0},
//     {8.2600331401756000e-01, 7.9955384841381302e-02, 0.0},
//     {1.1638649906727730e-01, 8.9602705800587407e-02, 0.0},
//     {2.0376848107772980e-01, 7.1788185898052326e-01, 0.0},
//     {6.4413220382260494e-02, 7.1008125956836521e-01, 0.0},
//     {9.5428585810584596e-02, 2.6077068256562902e-01, 0.0},
//     {2.4498296509349021e-01, 2.1117939909804931e-01, 0.0},
//     {7.0566724344036796e-02, 4.9732063377796598e-01, 0.0},
//     {6.1938125736255578e-01, 1.2512299505810390e-01, 0.0},
//     {6.2768261568031403e-01, 2.5015500335339208e-01, 0.0},
//     {4.2260565743346001e-01, 1.2953296900433620e-01, 0.0},
//     {2.1078525939140391e-01, 3.7986021093401962e-01, 0.0},
//     {4.0896380449124481e-01, 4.6631787462323071e-01, 0.0},
//     {2.1377743253005960e-01, 5.5802528953120256e-01, 0.0},
//     {4.0978657777002531e-01, 3.0141709320909299e-01, 0.0}},
//    {2.4440249073300249e-03,
//     3.3379500136836750e-03,
//     3.4227673271718501e-03,
//     3.5598757180403499e-03,
//     3.8572461868124248e-03,
//     4.8273543712181750e-03,
//     5.2546633678012501e-03,
//     5.3404218288141498e-03,
//     9.2418429056154005e-03,
//     9.2727402108032757e-03,
//     1.0310002059841049e-02,
//     1.0842542708505725e-02,
//     1.1245373099579100e-02,
//     1.2452036600753851e-02,
//     1.2549586713842551e-02,
//     1.3971867159939951e-02,
//     1.4072779302606675e-02,
//     1.4084827229864975e-02,
//     1.5264586206036251e-02,
//     1.5287638802019450e-02,
//     1.9786802896485999e-02,
//     2.0640943697731274e-02,
//     2.2968921082895850e-02,
//     2.3749787662653600e-02,
//     2.4074402518453650e-02,
//     2.5482462438393926e-02,
//     2.6676041524410474e-02,
//     2.7073436306583775e-02,
//     2.9718916975567701e-02,
//     2.9994853663553051e-02,
//     3.1582273211328352e-02,
//     3.7611031301662198e-02}
//  },
//  {
//    36,
//    {{2.4293535159026700e-02, 9.4930592938464031e-01, 0.0},
//     {2.6519342772158901e-02, 2.4269513064041098e-02, 0.0},
//     {9.4921260235510574e-01, 2.6506796643724899e-02, 0.0},
//     {3.3775763749036999e-03, 4.7673164123630779e-01, 0.0},
//     {4.7576722981011582e-01, 5.1989218291019390e-01, 0.0},
//     {5.1907831934706850e-01, 5.5912706202052003e-03, 0.0},
//     {8.6168397453205303e-01, 1.3399604861818300e-02, 0.0},
//     {1.2492097599255590e-01, 8.6130543213341393e-01, 0.0},
//     {1.3856545386105401e-02, 1.2477337173584679e-01, 0.0},
//     {2.1188706422168500e-02, 8.4384383512226457e-01, 0.0},
//     {8.4322967872188404e-01, 1.3545636458303650e-01, 0.0},
//     {1.3542317978649990e-01, 2.1348282065620599e-02, 0.0},
//     {3.0888535106794068e-01, 2.2191966301360600e-02, 0.0},
//     {6.6850575951690716e-01, 3.0890128793894273e-01, 0.0},
//     {2.2654501255714700e-02, 6.6917099433209937e-01, 0.0},
//     {2.8085154087720221e-01, 6.9247181551062442e-01, 0.0},
//     {6.9224467490505948e-01, 2.6872334502594599e-02, 0.0},
//     {2.6861744711943400e-02, 2.8100939732219082e-01, 0.0},
//     {1.1417784854701600e-01, 7.9735814135857996e-01, 0.0},
//     {7.9748079220612744e-01, 8.7980650879088101e-02, 0.0},
//     {8.9280729389424204e-02, 1.1450205611275180e-01, 0.0},
//     {1.0524878924550450e-01, 6.6869041199220447e-01, 0.0},
//     {6.6630222807398454e-01, 2.2750516318320271e-01, 0.0},
//     {2.3078037375469529e-01, 1.0545725612213259e-01, 0.0},
//     {1.7050591575403051e-01, 5.1740643986577728e-01, 0.0},
//     {5.0865939730425092e-01, 3.1705238552093218e-01, 0.0},
//     {3.1418238622808309e-01, 1.8107063616590391e-01, 0.0},
//     {4.6174608178640142e-01, 4.6785945398040590e-01, 0.0},
//     {6.9308749608105902e-02, 4.6228560420845410e-01, 0.0},
//     {4.6519552592682439e-01, 7.2435780566898006e-02, 0.0},
//     {2.5786258578926041e-01, 6.1313950391771632e-01, 0.0},
//     {6.1126277667792195e-01, 1.3003608346093859e-01, 0.0},
//     {1.3051821359335081e-01, 2.5817138288836389e-01, 0.0},
//     {4.2814379918281070e-01, 2.3620059698167339e-01, 0.0},
//     {3.3569957837300629e-01, 4.3110263085883421e-01, 0.0},
//     {2.3054242988361631e-01, 3.4560139493758052e-01, 0.0}},
//    {4.1560249689274751e-03,
//     4.1702924944386254e-03,
//     4.1707642266642000e-03,
//     4.3920217520637501e-03,
//     4.6118665461224503e-03,
//     4.9485602546932503e-03,
//     5.0885098963676500e-03,
//     5.1713215985039248e-03,
//     5.2067841521619253e-03,
//     7.9454944569679001e-03,
//     8.0118008810179994e-03,
//     8.0151920286416749e-03,
//     1.0769148979569475e-02,
//     1.0961835383476750e-02,
//     1.0980241818319100e-02,
//     1.1998798092264475e-02,
//     1.2095156518331500e-02,
//     1.2121685584365801e-02,
//     1.3924112200605200e-02,
//     1.4025659108902326e-02,
//     1.4129753092321276e-02,
//     1.7232247266740025e-02,
//     1.7930333402221024e-02,
//     1.8186348024405299e-02,
//     1.9720183418426299e-02,
//     2.0252858637809099e-02,
//     2.0643132476366850e-02,
//     2.1051114183261624e-02,
//     2.1089638332624925e-02,
//     2.1299246712198301e-02,
//     2.2571133201311374e-02,
//     2.2857078587113475e-02,
//     2.2906976635229850e-02,
//     2.5639334372397973e-02,
//     2.5828991535335525e-02,
//     2.5896359179831598e-02}
//  },
//  {
//    42,
//    {{4.88963910362179e-01, 4.88963910362179e-01, 0.0},
//     {4.88963910362179e-01, 2.20721792756430e-02, 0.0},
//     {2.20721792756430e-02, 4.88963910362179e-01, 0.0},
//     {4.17644719340454e-01, 4.17644719340454e-01, 0.0},
//     {4.17644719340454e-01, 1.64710561319092e-01, 0.0},
//     {1.64710561319092e-01, 4.17644719340454e-01, 0.0},
//     {2.73477528308839e-01, 2.73477528308839e-01, 0.0},
//     {2.73477528308839e-01, 4.53044943382323e-01, 0.0},
//     {4.53044943382323e-01, 2.73477528308839e-01, 0.0},
//     {1.77205532412543e-01, 1.77205532412543e-01, 0.0},
//     {1.77205532412543e-01, 6.45588935174913e-01, 0.0},
//     {6.45588935174913e-01, 1.77205532412543e-01, 0.0},
//     {6.17998830908730e-02, 6.17998830908730e-02, 0.0},
//     {6.17998830908730e-02, 8.76400233818255e-01, 0.0},
//     {8.76400233818255e-01, 6.17998830908730e-02, 0.0},
//     {1.93909612487010e-02, 1.93909612487010e-02, 0.0},
//     {1.93909612487010e-02, 9.61218077502598e-01, 0.0},
//     {9.61218077502598e-01, 1.93909612487010e-02, 0.0},
//     {1.72266687821356e-01, 7.70608554774996e-01, 0.0},
//     {7.70608554774996e-01, 5.71247574036480e-02, 0.0},
//     {5.71247574036480e-02, 1.72266687821356e-01, 0.0},
//     {1.72266687821356e-01, 5.71247574036480e-02, 0.0},
//     {7.70608554774996e-01, 1.72266687821356e-01, 0.0},
//     {5.71247574036480e-02, 7.70608554774996e-01, 0.0},
//     {3.36861459796345e-01, 5.70222290846683e-01, 0.0},
//     {5.70222290846683e-01, 9.29162493569720e-02, 0.0},
//     {9.29162493569720e-02, 3.36861459796345e-01, 0.0},
//     {3.36861459796345e-01, 9.29162493569720e-02, 0.0},
//     {5.70222290846683e-01, 3.36861459796345e-01, 0.0},
//     {9.29162493569720e-02, 5.70222290846683e-01, 0.0},
//     {2.98372882136258e-01, 6.86980167808088e-01, 0.0},
//     {6.86980167808088e-01, 1.46469500556540e-02, 0.0},
//     {1.46469500556540e-02, 2.98372882136258e-01, 0.0},
//     {2.98372882136258e-01, 1.46469500556540e-02, 0.0},
//     {6.86980167808088e-01, 2.98372882136258e-01, 0.0},
//     {1.46469500556540e-02, 6.86980167808088e-01, 0.0},
//     {1.18974497696957e-01, 8.79757171370171e-01, 0.0},
//     {8.79757171370171e-01, 1.26833093287199e-03, 0.0},
//     {1.26833093287199e-03, 1.18974497696957e-01, 0.0},
//     {1.18974497696957e-01, 1.26833093287199e-03, 0.0},
//     {8.79757171370171e-01, 1.18974497696957e-01, 0.0},
//     {1.26833093287199e-03, 8.79757171370171e-01, 0.0}},
//    {1.09417906847145e-02,
//     1.09417906847145e-02,
//     1.09417906847145e-02,
//     1.63941767720625e-02,
//     1.63941767720625e-02,
//     1.63941767720625e-02,
//     2.58870522536460e-02,
//     2.58870522536460e-02,
//     2.58870522536460e-02,
//     2.10812943684965e-02,
//     2.10812943684965e-02,
//     2.10812943684965e-02,
//     7.21684983488850e-03,
//     7.21684983488850e-03,
//     7.21684983488850e-03,
//     2.46170180120000e-03,
//     2.46170180120000e-03,
//     2.46170180120000e-03,
//     1.23328766062820e-02,
//     1.23328766062820e-02,
//     1.23328766062820e-02,
//     1.23328766062820e-02,
//     1.23328766062820e-02,
//     1.23328766062820e-02,
//     1.92857553935305e-02,
//     1.92857553935305e-02,
//     1.92857553935305e-02,
//     1.92857553935305e-02,
//     1.92857553935305e-02,
//     1.92857553935305e-02,
//     7.21815405676700e-03,
//     7.21815405676700e-03,
//     7.21815405676700e-03,
//     7.21815405676700e-03,
//     7.21815405676700e-03,
//     7.21815405676700e-03,
//     2.50511441925050e-03,
//     2.50511441925050e-03,
//     2.50511441925050e-03,
//     2.50511441925050e-03,
//     2.50511441925050e-03,
//     2.50511441925050e-03}
//  },
//  {
//    55,
//    {{1.0                   , 0.0                   , 0.0},
//     {0.0                   , 1.0                   , 0.0},
//     {0.0                   , 0.0                   , 0.0},
//     {9.3988635835771928e-01, 4.9848744634100996e-03, 0.0},
//     {5.4380668305835503e-02, 9.3864056186166756e-01, 0.0},
//     {9.3940049163876004e-03, 5.2642446269734702e-02, 0.0},
//     {1.6434508636240200e-02, 9.4690355173508323e-01, 0.0},
//     {9.4694872698624577e-01, 3.6337367716693998e-02, 0.0},
//     {4.2660400576765102e-02, 1.5122454179941101e-02, 0.0},
//     {1.2226949543872000e-02, 8.6937735106643133e-01, 0.0},
//     {8.6736965210466677e-01, 1.2049172857742969e-01, 0.0},
//     {8.4567440213890721e-01, 1.5776396787000199e-02, 0.0},
//     {1.3957596321026139e-01, 8.4481208703747090e-01, 0.0},
//     {1.3178217432308281e-01, 1.3500960558402201e-02, 0.0},
//     {1.5795512630024801e-02, 1.4552749385359881e-01, 0.0},
//     {7.3654628844363068e-01, 1.5569754090822801e-02, 0.0},
//     {1.3968843033038900e-02, 7.3798368944501946e-01, 0.0},
//     {2.5478951860390298e-01, 7.2976156897705524e-01, 0.0},
//     {7.3163865225549030e-01, 2.5430766833150520e-01, 0.0},
//     {1.5725372895084501e-02, 2.6962397957906031e-01, 0.0},
//     {2.6623028436468249e-01, 1.4478395630801300e-02, 0.0},
//     {8.6735040652140771e-01, 5.9167941040048203e-02, 0.0},
//     {7.4149366695661204e-02, 8.6347825750608687e-01, 0.0},
//     {1.5928594836003299e-02, 4.1912389552381862e-01, 0.0},
//     {1.5606102806777700e-02, 5.8092229211457624e-01, 0.0},
//     {5.9100948174838852e-01, 1.5925145265094101e-02, 0.0},
//     {4.0347714968887188e-01, 5.8067003681039198e-01, 0.0},
//     {5.6947456285259768e-01, 4.1494951463020030e-01, 0.0},
//     {6.7849370065030001e-02, 7.6121867859137604e-02, 0.0},
//     {4.2659685902715933e-01, 1.5750969231154401e-02, 0.0},
//     {6.7098250788970207e-02, 7.7418983124212093e-01, 0.0},
//     {7.5283102314795158e-01, 8.1911949563924294e-02, 0.0},
//     {7.7537277835568841e-01, 1.5771284572917341e-01, 0.0},
//     {1.6890731577873661e-01, 7.5039430997422452e-01, 0.0},
//     {1.6873358329194171e-01, 7.0831150726781894e-02, 0.0},
//     {8.2124470843632405e-02, 1.7629966267710759e-01, 0.0},
//     {6.2887053633447976e-01, 8.0774495331656301e-02, 0.0},
//     {8.1141301526575199e-02, 3.0543735897757762e-01, 0.0},
//     {2.9691120650804809e-01, 6.2274859888709300e-01, 0.0},
//     {7.6754231417057298e-02, 6.2472471495456661e-01, 0.0},
//     {6.2230223338447721e-01, 3.0114858211656370e-01, 0.0},
//     {3.1037862880509631e-01, 7.7909836507944599e-02, 0.0},
//     {8.1921821518658594e-02, 4.6036330383508761e-01, 0.0},
//     {4.7170226650134689e-01, 8.2155400679671906e-02, 0.0},
//     {4.5466034152504742e-01, 4.6375650338896440e-01, 0.0},
//     {1.7010913392369389e-01, 6.4222778081881993e-01, 0.0},
//     {6.4060043294867430e-01, 1.8982935372556059e-01, 0.0},
//     {1.9122675837165989e-01, 1.7399556853425760e-01, 0.0},
//     {1.8853157670702370e-01, 4.7989140704057581e-01, 0.0},
//     {4.7729299576907452e-01, 3.3483565981193042e-01, 0.0},
//     {3.1269746217597721e-01, 4.9579721972587398e-01, 0.0},
//     {4.9612259459456259e-01, 1.9275536689044351e-01, 0.0},
//     {1.9288053128670610e-01, 3.1610158072607569e-01, 0.0},
//     {3.3600414538164958e-01, 1.8948928012898231e-01, 0.0},
//     {3.3372805508479741e-01, 3.3435710218114523e-01, 0.0}},
//    {1.5506499627784999e-04,
//     1.5787936779320001e-04,
//     1.7716503897180001e-04,
//     1.3790929042020500e-03,
//     1.5673101913945000e-03,
//     1.9632852206504501e-03,
//     2.3637870966120248e-03,
//     2.4456127817772499e-03,
//     2.4965410872361499e-03,
//     3.4388454704036252e-03,
//     3.5244794510020999e-03,
//     3.7411716084290001e-03,
//     3.9024375902997751e-03,
//     3.9420923337041246e-03,
//     4.3948636595678749e-03,
//     5.1028460067507001e-03,
//     5.2390719653994996e-03,
//     5.2678353249450997e-03,
//     5.4411690050507498e-03,
//     5.5572102174651500e-03,
//     5.6046673420516247e-03,
//     5.7530654248287997e-03,
//     5.9203475624943751e-03,
//     6.4366160841976749e-03,
//     6.4489200402012004e-03,
//     6.4518081902498253e-03,
//     6.5085808014669752e-03,
//     6.6442035402279253e-03,
//     6.6446190457719254e-03,
//     6.6883082309448748e-03,
//     9.3946951660218506e-03,
//     9.5766473548822492e-03,
//     9.6212423756326000e-03,
//     9.7404956463108747e-03,
//     9.8651027886874250e-03,
//     1.0309119452445125e-02,
//     1.2821810962084575e-02,
//     1.2910141048365926e-02,
//     1.2955751056727700e-02,
//     1.3213199704525376e-02,
//     1.3462639325681701e-02,
//     1.3547383232981901e-02,
//     1.4618428661110900e-02,
//     1.4821579209082151e-02,
//     1.4858956918716301e-02,
//     1.5795006396574451e-02,
//     1.5823171128833101e-02,
//     1.6017684044289202e-02,
//     2.0301014897957576e-02,
//     2.0360937838258826e-02,
//     2.0366980031034525e-02,
//     2.0376263702112225e-02,
//     2.0379116623473949e-02,
//     2.0423276490578225e-02,
//     2.3080458363263227e-02}
//  },
//  {
//    55,
//    {{1.0                   , 0.0                   , 0.0},
//     {0.0                   , 1.0                   , 0.0},
//     {0.0                   , 0.0                   , 0.0},
//     {9.3988635835771928e-01, 4.9848744634100996e-03, 0.0},
//     {5.4380668305835503e-02, 9.3864056186166756e-01, 0.0},
//     {9.3940049163876004e-03, 5.2642446269734702e-02, 0.0},
//     {1.6434508636240200e-02, 9.4690355173508323e-01, 0.0},
//     {9.4694872698624577e-01, 3.6337367716693998e-02, 0.0},
//     {4.2660400576765102e-02, 1.5122454179941101e-02, 0.0},
//     {1.2226949543872000e-02, 8.6937735106643133e-01, 0.0},
//     {8.6736965210466677e-01, 1.2049172857742969e-01, 0.0},
//     {8.4567440213890721e-01, 1.5776396787000199e-02, 0.0},
//     {1.3957596321026139e-01, 8.4481208703747090e-01, 0.0},
//     {1.3178217432308281e-01, 1.3500960558402201e-02, 0.0},
//     {1.5795512630024801e-02, 1.4552749385359881e-01, 0.0},
//     {7.3654628844363068e-01, 1.5569754090822801e-02, 0.0},
//     {1.3968843033038900e-02, 7.3798368944501946e-01, 0.0},
//     {2.5478951860390298e-01, 7.2976156897705524e-01, 0.0},
//     {7.3163865225549030e-01, 2.5430766833150520e-01, 0.0},
//     {1.5725372895084501e-02, 2.6962397957906031e-01, 0.0},
//     {2.6623028436468249e-01, 1.4478395630801300e-02, 0.0},
//     {8.6735040652140771e-01, 5.9167941040048203e-02, 0.0},
//     {7.4149366695661204e-02, 8.6347825750608687e-01, 0.0},
//     {1.5928594836003299e-02, 4.1912389552381862e-01, 0.0},
//     {1.5606102806777700e-02, 5.8092229211457624e-01, 0.0},
//     {5.9100948174838852e-01, 1.5925145265094101e-02, 0.0},
//     {4.0347714968887188e-01, 5.8067003681039198e-01, 0.0},
//     {5.6947456285259768e-01, 4.1494951463020030e-01, 0.0},
//     {6.7849370065030001e-02, 7.6121867859137604e-02, 0.0},
//     {4.2659685902715933e-01, 1.5750969231154401e-02, 0.0},
//     {6.7098250788970207e-02, 7.7418983124212093e-01, 0.0},
//     {7.5283102314795158e-01, 8.1911949563924294e-02, 0.0},
//     {7.7537277835568841e-01, 1.5771284572917341e-01, 0.0},
//     {1.6890731577873661e-01, 7.5039430997422452e-01, 0.0},
//     {1.6873358329194171e-01, 7.0831150726781894e-02, 0.0},
//     {8.2124470843632405e-02, 1.7629966267710759e-01, 0.0},
//     {6.2887053633447976e-01, 8.0774495331656301e-02, 0.0},
//     {8.1141301526575199e-02, 3.0543735897757762e-01, 0.0},
//     {2.9691120650804809e-01, 6.2274859888709300e-01, 0.0},
//     {7.6754231417057298e-02, 6.2472471495456661e-01, 0.0},
//     {6.2230223338447721e-01, 3.0114858211656370e-01, 0.0},
//     {3.1037862880509631e-01, 7.7909836507944599e-02, 0.0},
//     {8.1921821518658594e-02, 4.6036330383508761e-01, 0.0},
//     {4.7170226650134689e-01, 8.2155400679671906e-02, 0.0},
//     {4.5466034152504742e-01, 4.6375650338896440e-01, 0.0},
//     {1.7010913392369389e-01, 6.4222778081881993e-01, 0.0},
//     {6.4060043294867430e-01, 1.8982935372556059e-01, 0.0},
//     {1.9122675837165989e-01, 1.7399556853425760e-01, 0.0},
//     {1.8853157670702370e-01, 4.7989140704057581e-01, 0.0},
//     {4.7729299576907452e-01, 3.3483565981193042e-01, 0.0},
//     {3.1269746217597721e-01, 4.9579721972587398e-01, 0.0},
//     {4.9612259459456259e-01, 1.9275536689044351e-01, 0.0},
//     {1.9288053128670610e-01, 3.1610158072607569e-01, 0.0},
//     {3.3600414538164958e-01, 1.8948928012898231e-01, 0.0},
//     {3.3372805508479741e-01, 3.3435710218114523e-01, 0.0}},
//    {1.5506499627784999e-04,
//     1.5787936779320001e-04,
//     1.7716503897180001e-04,
//     1.3790929042020500e-03,
//     1.5673101913945000e-03,
//     1.9632852206504501e-03,
//     2.3637870966120248e-03,
//     2.4456127817772499e-03,
//     2.4965410872361499e-03,
//     3.4388454704036252e-03,
//     3.5244794510020999e-03,
//     3.7411716084290001e-03,
//     3.9024375902997751e-03,
//     3.9420923337041246e-03,
//     4.3948636595678749e-03,
//     5.1028460067507001e-03,
//     5.2390719653994996e-03,
//     5.2678353249450997e-03,
//     5.4411690050507498e-03,
//     5.5572102174651500e-03,
//     5.6046673420516247e-03,
//     5.7530654248287997e-03,
//     5.9203475624943751e-03,
//     6.4366160841976749e-03,
//     6.4489200402012004e-03,
//     6.4518081902498253e-03,
//     6.5085808014669752e-03,
//     6.6442035402279253e-03,
//     6.6446190457719254e-03,
//     6.6883082309448748e-03,
//     9.3946951660218506e-03,
//     9.5766473548822492e-03,
//     9.6212423756326000e-03,
//     9.7404956463108747e-03,
//     9.8651027886874250e-03,
//     1.0309119452445125e-02,
//     1.2821810962084575e-02,
//     1.2910141048365926e-02,
//     1.2955751056727700e-02,
//     1.3213199704525376e-02,
//     1.3462639325681701e-02,
//     1.3547383232981901e-02,
//     1.4618428661110900e-02,
//     1.4821579209082151e-02,
//     1.4858956918716301e-02,
//     1.5795006396574451e-02,
//     1.5823171128833101e-02,
//     1.6017684044289202e-02,
//     2.0301014897957576e-02,
//     2.0360937838258826e-02,
//     2.0366980031034525e-02,
//     2.0376263702112225e-02,
//     2.0379116623473949e-02,
//     2.0423276490578225e-02,
//     2.3080458363263227e-02}
//  },
//  {
//    61,
//    {{3.33333333333333e-01, 3.33333333333333e-01, 0.0},
//     {4.97170540556774e-01, 4.97170540556774e-01, 0.0},
//     {4.97170540556774e-01, 5.65891888645198e-03, 0.0},
//     {5.65891888645198e-03, 4.97170540556774e-01, 0.0},
//     {4.82176322624625e-01, 4.82176322624625e-01, 0.0},
//     {4.82176322624625e-01, 3.56473547507510e-02, 0.0},
//     {3.56473547507510e-02, 4.82176322624625e-01, 0.0},
//     {4.50239969020782e-01, 4.50239969020782e-01, 0.0},
//     {4.50239969020782e-01, 9.95200619584370e-02, 0.0},
//     {9.95200619584370e-02, 4.50239969020782e-01, 0.0},
//     {4.00266239377397e-01, 4.00266239377397e-01, 0.0},
//     {4.00266239377397e-01, 1.99467521245206e-01, 0.0},
//     {1.99467521245206e-01, 4.00266239377397e-01, 0.0},
//     {2.52141267970953e-01, 2.52141267970953e-01, 0.0},
//     {2.52141267970953e-01, 4.95717464058095e-01, 0.0},
//     {4.95717464058095e-01, 2.52141267970953e-01, 0.0},
//     {1.62047004658461e-01, 1.62047004658461e-01, 0.0},
//     {1.62047004658461e-01, 6.75905990683077e-01, 0.0},
//     {6.75905990683077e-01, 1.62047004658461e-01, 0.0},
//     {7.58758822607460e-02, 7.58758822607460e-02, 0.0},
//     {7.58758822607460e-02, 8.48248235478508e-01, 0.0},
//     {8.48248235478508e-01, 7.58758822607460e-02, 0.0},
//     {1.56547269678220e-02, 1.56547269678220e-02, 0.0},
//     {1.56547269678220e-02, 9.68690546064356e-01, 0.0},
//     {9.68690546064356e-01, 1.56547269678220e-02, 0.0},
//     {3.34319867363658e-01, 6.55493203809423e-01, 0.0},
//     {6.55493203809423e-01, 1.01869288269190e-02, 0.0},
//     {1.01869288269190e-02, 3.34319867363658e-01, 0.0},
//     {3.34319867363658e-01, 1.01869288269190e-02, 0.0},
//     {6.55493203809423e-01, 3.34319867363658e-01, 0.0},
//     {1.01869288269190e-02, 6.55493203809423e-01, 0.0},
//     {2.92221537796944e-01, 5.72337590532020e-01, 0.0},
//     {5.72337590532020e-01, 1.35440871671036e-01, 0.0},
//     {1.35440871671036e-01, 2.92221537796944e-01, 0.0},
//     {2.92221537796944e-01, 1.35440871671036e-01, 0.0},
//     {5.72337590532020e-01, 2.92221537796944e-01, 0.0},
//     {1.35440871671036e-01, 5.72337590532020e-01, 0.0},
//     {3.19574885423190e-01, 6.26001190286228e-01, 0.0},
//     {6.26001190286228e-01, 5.44239242905830e-02, 0.0},
//     {5.44239242905830e-02, 3.19574885423190e-01, 0.0},
//     {3.19574885423190e-01, 5.44239242905830e-02, 0.0},
//     {6.26001190286228e-01, 3.19574885423190e-01, 0.0},
//     {5.44239242905830e-02, 6.26001190286228e-01, 0.0},
//     {1.90704224192292e-01, 7.96427214974071e-01, 0.0},
//     {7.96427214974071e-01, 1.28685608336370e-02, 0.0},
//     {1.28685608336370e-02, 1.90704224192292e-01, 0.0},
//     {1.90704224192292e-01, 1.28685608336370e-02, 0.0},
//     {7.96427214974071e-01, 1.90704224192292e-01, 0.0},
//     {1.28685608336370e-02, 7.96427214974071e-01, 0.0},
//     {1.80483211648746e-01, 7.52351005937729e-01, 0.0},
//     {7.52351005937729e-01, 6.71657824135240e-02, 0.0},
//     {6.71657824135240e-02, 1.80483211648746e-01, 0.0},
//     {1.80483211648746e-01, 6.71657824135240e-02, 0.0},
//     {7.52351005937729e-01, 1.80483211648746e-01, 0.0},
//     {6.71657824135240e-02, 7.52351005937729e-01, 0.0},
//     {8.07113136795640e-02, 9.04625504095608e-01, 0.0},
//     {9.04625504095608e-01, 1.46631822248280e-02, 0.0},
//     {1.46631822248280e-02, 8.07113136795640e-02, 0.0},
//     {8.07113136795640e-02, 1.46631822248280e-02, 0.0},
//     {9.04625504095608e-01, 8.07113136795640e-02, 0.0},
//     {1.46631822248280e-02, 9.04625504095608e-01, 0.0}},
//    {1.67185996454015e-02,
//     2.54670772025350e-03,
//     2.54670772025350e-03,
//     2.54670772025350e-03,
//     7.33543226381900e-03,
//     7.33543226381900e-03,
//     7.33543226381900e-03,
//     1.21754391768360e-02,
//     1.21754391768360e-02,
//     1.21754391768360e-02,
//     1.55537754344845e-02,
//     1.55537754344845e-02,
//     1.55537754344845e-02,
//     1.56285556093100e-02,
//     1.56285556093100e-02,
//     1.56285556093100e-02,
//     1.24078271698325e-02,
//     1.24078271698325e-02,
//     1.24078271698325e-02,
//     7.02803653527850e-03,
//     7.02803653527850e-03,
//     7.02803653527850e-03,
//     1.59733808688950e-03,
//     1.59733808688950e-03,
//     1.59733808688950e-03,
//     4.05982765949650e-03,
//     4.05982765949650e-03,
//     4.05982765949650e-03,
//     4.05982765949650e-03,
//     4.05982765949650e-03,
//     4.05982765949650e-03,
//     1.34028711415815e-02,
//     1.34028711415815e-02,
//     1.34028711415815e-02,
//     1.34028711415815e-02,
//     1.34028711415815e-02,
//     1.34028711415815e-02,
//     9.22999660541100e-03,
//     9.22999660541100e-03,
//     9.22999660541100e-03,
//     9.22999660541100e-03,
//     9.22999660541100e-03,
//     9.22999660541100e-03,
//     4.23843426716400e-03,
//     4.23843426716400e-03,
//     4.23843426716400e-03,
//     4.23843426716400e-03,
//     4.23843426716400e-03,
//     4.23843426716400e-03,
//     9.14639838501250e-03,
//     9.14639838501250e-03,
//     9.14639838501250e-03,
//     9.14639838501250e-03,
//     9.14639838501250e-03,
//     9.14639838501250e-03,
//     3.33281600208250e-03,
//     3.33281600208250e-03,
//     3.33281600208250e-03,
//     3.33281600208250e-03,
//     3.33281600208250e-03,
//     3.33281600208250e-03}
//  },
//  {
//    66,
//    {{1.1673105966841200e-02, 9.8125659512890129e-01, 0.0},
//     {9.8100308583879503e-01, 7.1462504863216000e-03, 0.0},
//     {1.0696631709169700e-02, 1.1515393337596600e-02, 0.0},
//     {9.3824769835505051e-01, 4.9557059134064198e-02, 0.0},
//     {1.2662751841721401e-02, 9.3701236206150318e-01, 0.0},
//     {5.9810940998380198e-02, 1.2136457892184800e-02, 0.0},
//     {1.3736329792672100e-02, 6.1278362559696799e-02, 0.0},
//     {9.2295279594054480e-01, 1.4112827060242099e-02, 0.0},
//     {6.3310735499269494e-02, 9.2201972917274344e-01, 0.0},
//     {1.1726510033460201e-02, 1.5005204752290349e-01, 0.0},
//     {1.5547205873234721e-01, 8.3251471215892492e-01, 0.0},
//     {8.3432938889821573e-01, 1.2522815875883600e-02, 0.0},
//     {8.5016380319567597e-01, 1.3719975087357841e-01, 0.0},
//     {1.2881635052197599e-02, 8.4776270634792006e-01, 0.0},
//     {1.5108016089587781e-01, 1.3652692403937501e-02, 0.0},
//     {1.0191787921658400e-02, 5.7704386183448575e-01, 0.0},
//     {2.8133723993032811e-01, 7.0668537596231984e-01, 0.0},
//     {7.1243746285009335e-01, 1.2456978098990301e-02, 0.0},
//     {2.7630252508633668e-01, 1.2174131138564200e-02, 0.0},
//     {1.0965836856061799e-02, 4.1943067124662847e-01, 0.0},
//     {4.2891105178839167e-01, 5.5996160674689166e-01, 0.0},
//     {4.2154205551147350e-01, 1.1647599478465699e-02, 0.0},
//     {5.7112585904443613e-01, 1.1821831398851200e-02, 0.0},
//     {5.8268682705109343e-01, 4.0578895811771831e-01, 0.0},
//     {1.3056780671324699e-02, 2.7250237508679159e-01, 0.0},
//     {1.3076040096391800e-02, 7.2247125232334730e-01, 0.0},
//     {7.2634370624067746e-01, 2.6029840192506443e-01, 0.0},
//     {6.8723006863737404e-02, 6.3141727720962701e-02, 0.0},
//     {8.6523021015294610e-01, 7.2061183733767895e-02, 0.0},
//     {6.4859907103736694e-02, 8.5904335439099433e-01, 0.0},
//     {1.4834949433620581e-01, 7.8887883522396707e-01, 0.0},
//     {6.2435989839593801e-02, 1.4939354993542750e-01, 0.0},
//     {7.8713690117350699e-01, 6.5638204275659501e-02, 0.0},
//     {5.1910492160953101e-02, 5.2556356956052430e-01, 0.0},
//     {1.5431299274438229e-01, 7.1638392691700595e-02, 0.0},
//     {2.6178427456029407e-01, 6.2147948528815097e-02, 0.0},
//     {7.6672578728127994e-01, 1.6582115548313259e-01, 0.0},
//     {2.5821036766273009e-01, 6.8001197661390189e-01, 0.0},
//     {6.7906592514742597e-02, 7.5715154377818017e-01, 0.0},
//     {5.2935782748041971e-01, 4.1215038411072058e-01, 0.0},
//     {6.6603615048415998e-02, 2.6125130878864999e-01, 0.0},
//     {5.8567546189943198e-02, 3.9022361145349760e-01, 0.0},
//     {6.4453536041083406e-02, 6.3736265597609554e-01, 0.0},
//     {6.7481384291513691e-01, 6.3758334206129100e-02, 0.0},
//     {3.9146023103687089e-01, 5.5032380905631106e-01, 0.0},
//     {6.4877014923071408e-01, 2.8367283602629478e-01, 0.0},
//     {3.9464982204080379e-01, 6.0517552255370900e-02, 0.0},
//     {5.3901371519333352e-01, 6.1199017693642201e-02, 0.0},
//     {1.6278950827847499e-01, 6.8613221410348235e-01, 0.0},
//     {6.8124363226406448e-01, 1.5679683458990931e-01, 0.0},
//     {1.5428328780201980e-01, 1.6675126240198401e-01, 0.0},
//     {2.5227277504445078e-01, 2.5048039333948502e-01, 0.0},
//     {2.5479815324070432e-01, 4.9940906490431908e-01, 0.0},
//     {1.4855805491943541e-01, 5.7560230960873771e-01, 0.0},
//     {2.9302396064361819e-01, 5.6568973541618528e-01, 0.0},
//     {2.8089912723099042e-01, 1.4379215742477949e-01, 0.0},
//     {4.8209895929708219e-01, 2.5185575358650381e-01, 0.0},
//     {5.6418782454436134e-01, 1.4629667431525920e-01, 0.0},
//     {1.3076996443439021e-01, 4.4895775861167753e-01, 0.0},
//     {1.4796922219475581e-01, 3.0011743868291701e-01, 0.0},
//     {5.6386842229459166e-01, 2.8137720892975088e-01, 0.0},
//     {4.3611574287904659e-01, 4.2520534464204729e-01, 0.0},
//     {3.6032639352854701e-01, 2.5991900048886368e-01, 0.0},
//     {4.2241883346742481e-01, 1.4532384433026860e-01, 0.0},
//     {3.7190018330523877e-01, 3.7801227035670099e-01, 0.0},
//     {2.4136450069284729e-01, 3.8475632849397318e-01, 0.0}},
//    {6.2914392466127504e-04,
//     6.3183630018052495e-04,
//     8.3173238332957496e-04,
//     2.0375873031348501e-03,
//     2.1533881435404252e-03,
//     2.1946686544827001e-03,
//     2.4274896390419250e-03,
//     2.5616552978717250e-03,
//     2.7099422085185751e-03,
//     3.2346347543960748e-03,
//     3.4084955895738499e-03,
//     3.4619332036664001e-03,
//     3.4855385026212752e-03,
//     3.6030349991898998e-03,
//     3.8425863883506501e-03,
//     4.0622450563140748e-03,
//     4.2429576070036248e-03,
//     4.2522133105332002e-03,
//     4.2738380168663498e-03,
//     4.3472213639774750e-03,
//     4.3635990609680250e-03,
//     4.4601689321659248e-03,
//     4.4611715969841001e-03,
//     4.4761584388088504e-03,
//     4.5314939050175749e-03,
//     4.6196209720506748e-03,
//     4.6448391092779500e-03,
//     5.0804287944135754e-03,
//     5.3442915452291753e-03,
//     5.7979213524570498e-03,
//     6.8606677714776751e-03,
//     7.2575480585093501e-03,
//     7.3630684626356750e-03,
//     7.4859062907272752e-03,
//     7.6756737029697996e-03,
//     8.1315841465677500e-03,
//     8.1971052126526001e-03,
//     8.2808668797998992e-03,
//     8.6541881718643493e-03,
//     8.6770343494031749e-03,
//     8.6843012350962742e-03,
//     8.7132190613550004e-03,
//     8.7150390296491503e-03,
//     8.8867892493721002e-03,
//     9.0045749095674504e-03,
//     9.0731571460701751e-03,
//     9.5474425520800498e-03,
//     9.8063200029475748e-03,
//     1.2067753147187676e-02,
//     1.2247803039156025e-02,
//     1.2430520846804525e-02,
//     1.2676643424644550e-02,
//     1.2744299851074425e-02,
//     1.3034001591679976e-02,
//     1.3086521873117450e-02,
//     1.3111017088793300e-02,
//     1.3186491120564975e-02,
//     1.3236226593190750e-02,
//     1.3559889862520475e-02,
//     1.3586755085482525e-02,
//     1.3677513715971475e-02,
//     1.3932208647817225e-02,
//     1.4443356605827800e-02,
//     1.4634844540566926e-02,
//     1.5225981266989826e-02,
//     1.5931849111237351e-02}
//  },
//  {
//    73,
//    {{3.33333333333333e-01, 3.33333333333333e-01, 0.0},
//     {4.89609987073006e-01, 4.89609987073006e-01, 0.0},
//     {4.89609987073006e-01, 2.07800258539870e-02, 0.0},
//     {2.07800258539870e-02, 4.89609987073006e-01, 0.0},
//     {4.54536892697893e-01, 4.54536892697893e-01, 0.0},
//     {4.54536892697893e-01, 9.09262146042150e-02, 0.0},
//     {9.09262146042150e-02, 4.54536892697893e-01, 0.0},
//     {4.01416680649431e-01, 4.01416680649431e-01, 0.0},
//     {4.01416680649431e-01, 1.97166638701138e-01, 0.0},
//     {1.97166638701138e-01, 4.01416680649431e-01, 0.0},
//     {2.55551654403098e-01, 2.55551654403098e-01, 0.0},
//     {2.55551654403098e-01, 4.88896691193805e-01, 0.0},
//     {4.88896691193805e-01, 2.55551654403098e-01, 0.0},
//     {1.77077942152130e-01, 1.77077942152130e-01, 0.0},
//     {1.77077942152130e-01, 6.45844115695741e-01, 0.0},
//     {6.45844115695741e-01, 1.77077942152130e-01, 0.0},
//     {1.10061053227952e-01, 1.10061053227952e-01, 0.0},
//     {1.10061053227952e-01, 7.79877893544096e-01, 0.0},
//     {7.79877893544096e-01, 1.10061053227952e-01, 0.0},
//     {5.55286242518400e-02, 5.55286242518400e-02, 0.0},
//     {5.55286242518400e-02, 8.88942751496321e-01, 0.0},
//     {8.88942751496321e-01, 5.55286242518400e-02, 0.0},
//     {1.26218637772290e-02, 1.26218637772290e-02, 0.0},
//     {1.26218637772290e-02, 9.74756272445543e-01, 0.0},
//     {9.74756272445543e-01, 1.26218637772290e-02, 0.0},
//     {3.95754787356943e-01, 6.00633794794645e-01, 0.0},
//     {6.00633794794645e-01, 3.61141784841201e-03, 0.0},
//     {3.61141784841201e-03, 3.95754787356943e-01, 0.0},
//     {3.95754787356943e-01, 3.61141784841201e-03, 0.0},
//     {6.00633794794645e-01, 3.95754787356943e-01, 0.0},
//     {3.61141784841201e-03, 6.00633794794643e-01, 0.0},
//     {3.07929983880436e-01, 5.57603261588784e-01, 0.0},
//     {5.57603261588784e-01, 1.34466754530780e-01, 0.0},
//     {1.34466754530780e-01, 3.07929983880436e-01, 0.0},
//     {3.07929983880436e-01, 1.34466754530780e-01, 0.0},
//     {5.57603261588784e-01, 3.07929983880436e-01, 0.0},
//     {1.34466754530780e-01, 5.57603261588784e-01, 0.0},
//     {2.64566948406520e-01, 7.20987025817365e-01, 0.0},
//     {7.20987025817365e-01, 1.44460257761150e-02, 0.0},
//     {1.44460257761150e-02, 2.64566948406520e-01, 0.0},
//     {2.64566948406520e-01, 1.44460257761150e-02, 0.0},
//     {7.20987025817365e-01, 2.64566948406520e-01, 0.0},
//     {1.44460257761150e-02, 7.20987025817365e-01, 0.0},
//     {3.58539352205951e-01, 5.94527068955871e-01, 0.0},
//     {5.94527068955871e-01, 4.69335788381780e-02, 0.0},
//     {4.69335788381780e-02, 3.58539352205951e-01, 0.0},
//     {3.58539352205951e-01, 4.69335788381780e-02, 0.0},
//     {5.94527068955871e-01, 3.58539352205951e-01, 0.0},
//     {4.69335788381780e-02, 5.94527068955871e-01, 0.0},
//     {1.57807405968595e-01, 8.39331473680839e-01, 0.0},
//     {8.39331473680839e-01, 2.86112035056701e-03, 0.0},
//     {2.86112035056701e-03, 1.57807405968595e-01, 0.0},
//     {1.57807405968595e-01, 2.86112035056701e-03, 0.0},
//     {8.39331473680839e-01, 1.57807405968595e-01, 0.0},
//     {2.86112035056701e-03, 8.39331473680839e-01, 0.0},
//     {7.50505969759110e-02, 7.01087978926173e-01, 0.0},
//     {7.01087978926173e-01, 2.23861424097916e-01, 0.0},
//     {2.23861424097916e-01, 7.50505969759110e-02, 0.0},
//     {7.50505969759110e-02, 2.23861424097916e-01, 0.0},
//     {7.01087978926173e-01, 7.50505969759110e-02, 0.0},
//     {2.23861424097916e-01, 7.01087978926173e-01, 0.0},
//     {1.42421601113383e-01, 8.22931324069857e-01, 0.0},
//     {8.22931324069857e-01, 3.46470748167600e-02, 0.0},
//     {3.46470748167600e-02, 1.42421601113383e-01, 0.0},
//     {1.42421601113383e-01, 3.46470748167600e-02, 0.0},
//     {8.22931324069857e-01, 1.42421601113383e-01, 0.0},
//     {3.46470748167600e-02, 8.22931324069857e-01, 0.0},
//     {6.54946280829380e-02, 9.24344252620784e-01, 0.0},
//     {9.24344252620784e-01, 1.01611192962780e-02, 0.0},
//     {1.01611192962780e-02, 6.54946280829380e-02, 0.0},
//     {6.54946280829380e-02, 1.01611192962780e-02, 0.0},
//     {9.24344252620784e-01, 6.54946280829380e-02, 0.0},
//     {1.01611192962780e-02, 9.24344252620784e-01, 0.0}},
//    {1.64531656944595e-02,
//     5.16536594563600e-03,
//     5.16536594563600e-03,
//     5.16536594563600e-03,
//     1.11936236315080e-02,
//     1.11936236315080e-02,
//     1.11936236315080e-02,
//     1.51330629347340e-02,
//     1.51330629347340e-02,
//     1.51330629347340e-02,
//     1.52454839010990e-02,
//     1.52454839010990e-02,
//     1.52454839010990e-02,
//     1.20796063708205e-02,
//     1.20796063708205e-02,
//     1.20796063708205e-02,
//     8.02540179340050e-03,
//     8.02540179340050e-03,
//     8.02540179340050e-03,
//     4.04229013089200e-03,
//     4.04229013089200e-03,
//     4.04229013089200e-03,
//     1.03968101374250e-03,
//     1.03968101374250e-03,
//     1.03968101374250e-03,
//     1.94243845249050e-03,
//     1.94243845249050e-03,
//     1.94243845249050e-03,
//     1.94243845249050e-03,
//     1.94243845249050e-03,
//     1.94243845249050e-03,
//     1.27870803060110e-02,
//     1.27870803060110e-02,
//     1.27870803060110e-02,
//     1.27870803060110e-02,
//     1.27870803060110e-02,
//     1.27870803060110e-02,
//     4.44045178666900e-03,
//     4.44045178666900e-03,
//     4.44045178666900e-03,
//     4.44045178666900e-03,
//     4.44045178666900e-03,
//     4.44045178666900e-03,
//     8.06227338086550e-03,
//     8.06227338086550e-03,
//     8.06227338086550e-03,
//     8.06227338086550e-03,
//     8.06227338086550e-03,
//     8.06227338086550e-03,
//     1.24597090874550e-03,
//     1.24597090874550e-03,
//     1.24597090874550e-03,
//     1.24597090874550e-03,
//     1.24597090874550e-03,
//     1.24597090874550e-03,
//     9.12142005947550e-03,
//     9.12142005947550e-03,
//     9.12142005947550e-03,
//     9.12142005947550e-03,
//     9.12142005947550e-03,
//     9.12142005947550e-03,
//     5.12928186809950e-03,
//     5.12928186809950e-03,
//     5.12928186809950e-03,
//     5.12928186809950e-03,
//     5.12928186809950e-03,
//     5.12928186809950e-03,
//     1.89996442765100e-03,
//     1.89996442765100e-03,
//     1.89996442765100e-03,
//     1.89996442765100e-03,
//     1.89996442765100e-03,
//     1.89996442765100e-03}
//  },
//  {
//    78,
//    {{8.9411337112035999e-03, 8.6983293701984998e-03, 0.0},
//     {9.7926226298067365e-01, 1.0264413374365100e-02, 0.0},
//     {1.0547538211187800e-02, 9.7855142025151109e-01, 0.0},
//     {2.3777061947122002e-03, 6.3655109860361700e-02, 0.0},
//     {6.3042511579465998e-02, 4.1506347508631003e-03, 0.0},
//     {9.3084224967299967e-01, 4.8053482262546002e-03, 0.0},
//     {6.2907655549027400e-02, 9.3167900694812233e-01, 0.0},
//     {9.3159622463806491e-01, 6.2626488180135900e-02, 0.0},
//     {6.1951689414552003e-03, 9.2935870585640645e-01, 0.0},
//     {2.8712581923668101e-02, 3.1020212299716299e-02, 0.0},
//     {9.2938444783052321e-01, 3.4215296821852897e-02, 0.0},
//     {3.7545756662128102e-02, 9.2578688846693047e-01, 0.0},
//     {8.6895739063833997e-03, 1.5849712515099221e-01, 0.0},
//     {1.5475970539646791e-01, 8.3636066576882862e-01, 0.0},
//     {8.3310252941849239e-01, 8.9257244824476004e-03, 0.0},
//     {8.3742310735260950e-01, 1.5291673040783921e-01, 0.0},
//     {1.5593625052337881e-01, 9.4966240058029002e-03, 0.0},
//     {9.8599642095236004e-03, 8.3422114935955050e-01, 0.0},
//     {4.0558737332891631e-01, 7.4389302007913001e-03, 0.0},
//     {5.9647278986182350e-01, 3.9563308093107152e-01, 0.0},
//     {8.0747800415767006e-03, 4.0313194259026802e-01, 0.0},
//     {7.5073977720710996e-03, 5.8516095946805691e-01, 0.0},
//     {3.9367645192372991e-01, 5.9748965928987985e-01, 0.0},
//     {5.8465307262122179e-01, 8.7250464968192006e-03, 0.0},
//     {4.8708041121196383e-01, 2.0212922991194000e-02, 0.0},
//     {2.6835128117845169e-01, 7.2023400886682198e-01, 0.0},
//     {7.2239562887479880e-01, 2.6623993664561901e-01, 0.0},
//     {2.7168267423572212e-01, 1.1288269880823600e-02, 0.0},
//     {1.1258084204589300e-02, 7.1696959633251023e-01, 0.0},
//     {1.1503473436974001e-02, 2.7400671101656832e-01, 0.0},
//     {7.1405259005638033e-01, 1.1351156049706200e-02, 0.0},
//     {4.9028710531115449e-01, 4.9364918414683351e-01, 0.0},
//     {2.0142342520930698e-02, 4.8325734596013992e-01, 0.0},
//     {3.6110746485855001e-02, 9.3567950158201393e-02, 0.0},
//     {8.6079988198508572e-01, 3.9737906707539197e-02, 0.0},
//     {1.0058915260013050e-01, 8.5863434193517962e-01, 0.0},
//     {9.1874071705841595e-02, 3.9551300197337699e-02, 0.0},
//     {8.6048882961910289e-01, 9.6622405707924700e-02, 0.0},
//     {4.3984217867325599e-02, 8.5618863491067676e-01, 0.0},
//     {2.0110176067354310e-01, 7.4491158356262255e-01, 0.0},
//     {7.4499937262632787e-01, 5.3686563816580400e-02, 0.0},
//     {5.3218664130983202e-02, 1.9637542759350521e-01, 0.0},
//     {7.4539846474005178e-01, 1.9820658055500051e-01, 0.0},
//     {1.9572899328760179e-01, 5.5571383315608597e-02, 0.0},
//     {1.0925320579875419e-01, 6.1000361824130300e-01, 0.0},
//     {5.6762570200051501e-02, 7.4091218949591942e-01, 0.0},
//     {4.8383793347481101e-02, 6.0751356609779783e-01, 0.0},
//     {1.0806128097601329e-01, 1.1220815104370099e-01, 0.0},
//     {6.1856059009905007e-01, 2.6987537030349740e-01, 0.0},
//     {7.7212960134965625e-01, 1.1141173953329921e-01, 0.0},
//     {6.1157348011327173e-01, 3.3893676779306348e-01, 0.0},
//     {3.3813261033758418e-01, 4.9469393878745799e-02, 0.0},
//     {1.1730841282542900e-01, 7.6964513097951825e-01, 0.0},
//     {2.6745512605961458e-01, 1.1157188081540730e-01, 0.0},
//     {6.5421001600256889e-01, 1.9065483146999149e-01, 0.0},
//     {5.3829748115775802e-02, 3.3586168268491179e-01, 0.0},
//     {1.8488403241167711e-01, 1.5518315238513730e-01, 0.0},
//     {3.3762671047443338e-01, 6.0814025962944529e-01, 0.0},
//     {6.0671020344994708e-01, 5.4263279559821201e-02, 0.0},
//     {4.6126140854956371e-01, 6.8817667072165398e-02, 0.0},
//     {1.5254653656712561e-01, 6.5102408457488470e-01, 0.0},
//     {7.0058254354307500e-02, 4.6619043927415987e-01, 0.0},
//     {4.7042013790318088e-01, 4.6348264553531421e-01, 0.0},
//     {1.2164616937459330e-01, 2.3814948755156831e-01, 0.0},
//     {6.3714040527021165e-01, 1.2383993845133670e-01, 0.0},
//     {2.3799045151187120e-01, 6.3702164523263760e-01, 0.0},
//     {1.4839298571771459e-01, 4.8941885777801442e-01, 0.0},
//     {3.5980695715496319e-01, 1.4528808662532389e-01, 0.0},
//     {4.9414410550951349e-01, 3.6102163838181101e-01, 0.0},
//     {1.4406306879808209e-01, 3.5135083418870572e-01, 0.0},
//     {5.0197644400035468e-01, 1.4354916632930600e-01, 0.0},
//     {3.5554238342982608e-01, 5.0164915995018422e-01, 0.0},
//     {2.4434395407713269e-01, 2.4060521291041001e-01, 0.0},
//     {2.4370649893418969e-01, 5.1090172770553444e-01, 0.0},
//     {5.1222008073208247e-01, 2.4527379735428820e-01, 0.0},
//     {2.5260383151777532e-01, 3.7003195550936951e-01, 0.0},
//     {3.7598956528506539e-01, 2.5054066116305501e-01, 0.0},
//     {3.7290779871441049e-01, 3.7537502775491960e-01, 0.0}},
//    {5.4361363496487503e-04,
//     7.2467838162779998e-04,
//     7.7115073341352497e-04,
//     8.6004082758987495e-04,
//     1.0474618002941001e-03,
//     1.1184512874480750e-03,
//     1.1763605203547751e-03,
//     1.2216983937579000e-03,
//     1.2981910842235500e-03,
//     1.8518264745199500e-03,
//     1.9938852575230501e-03,
//     2.0887630727400748e-03,
//     2.4041665215934500e-03,
//     2.4079564462609749e-03,
//     2.4644365189580250e-03,
//     2.5664470075342249e-03,
//     2.5797025777773248e-03,
//     2.6572750407618498e-03,
//     2.6720326723837502e-03,
//     2.6742255252576500e-03,
//     2.7256615428406252e-03,
//     2.7474945893627501e-03,
//     2.8355763807253751e-03,
//     3.0133910732457001e-03,
//     3.4904798455149752e-03,
//     3.5286997884102249e-03,
//     3.5482586761459252e-03,
//     3.6053169067081750e-03,
//     3.6176086713727751e-03,
//     3.6237442468015752e-03,
//     3.6346693923518499e-03,
//     3.6491047731400249e-03,
//     3.6828644616571501e-03,
//     4.1865990826037754e-03,
//     4.2238875114384497e-03,
//     4.2355665721040998e-03,
//     4.3267543023872996e-03,
//     4.3631136623173503e-03,
//     4.4304305539834754e-03,
//     7.0706006005818254e-03,
//     7.1249178122030250e-03,
//     7.1251411634678250e-03,
//     7.5161805869424503e-03,
//     7.5507819270567997e-03,
//     7.5996784019234003e-03,
//     7.6417199018501249e-03,
//     7.6516853250572746e-03,
//     7.7332517050150247e-03,
//     7.7443455208702503e-03,
//     7.8286562636164757e-03,
//     7.8393373347939495e-03,
//     7.8580117321664505e-03,
//     7.8795535973380255e-03,
//     8.1062034496375505e-03,
//     8.6878038096496243e-03,
//     8.7598363731706003e-03,
//     8.7679355077581993e-03,
//     8.8032303833538994e-03,
//     8.8153876245372752e-03,
//     9.1600805085676005e-03,
//     9.1933276917576506e-03,
//     9.2918915734212758e-03,
//     9.3342892901571994e-03,
//     1.0099333664696650e-02,
//     1.0339501015941751e-02,
//     1.0548944796751400e-02,
//     1.2386275100928551e-02,
//     1.2510481528527850e-02,
//     1.2644864677864924e-02,
//     1.3000930254698475e-02,
//     1.3038339197138024e-02,
//     1.3122478808947751e-02,
//     1.4978994062894975e-02,
//     1.4990249935657125e-02,
//     1.4997881803219050e-02,
//     1.5853329586233849e-02,
//     1.5882796527709402e-02,
//     1.5930165141800151e-02}
//  }
    

        Kokkos::deep_copy(_quadrature_weights, quadrature_weights);
        Kokkos::deep_copy(_parameterized_quadrature_sites, parameterized_quadrature_sites);
    }
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

    Quadrature() {
        _is_initialized = false;
        _number_of_quadrature_points = 0;
        _order_of_quadrature_points = 0;
        _dimension_of_quadrature_points = 0;
        _qt = QuadratureType::INVALID;
    }

    Quadrature(const int order, const int dimension = 0, std::string quadrature_type = "LINE") {
        _number_of_quadrature_points = 0;
        _order_of_quadrature_points = order;
        _dimension_of_quadrature_points = dimension;
        _qt = this->parseQuadratureType(quadrature_type);
        if (dimension > 0) {
            // populates _number_of_quadrature_points
            this->generateQuadrature(order, dimension);
            _is_initialized = true;
        } else {
            _is_initialized = false;
        }
    }

///@}

/** @name Public Utility
 *  
 */
///@{

    static QuadratureType parseQuadratureType(std::string quadrature_type) { 
        transform(quadrature_type.begin(), quadrature_type.end(), quadrature_type.begin(), ::tolower);
        if (quadrature_type=="line") {
            return QuadratureType::LINE;
        } else if (quadrature_type=="tri" || quadrature_type=="triangle") {
            return QuadratureType::TRI;
        } else if (quadrature_type=="") {
            return QuadratureType::INVALID;
        } else {
            compadre_assert_release(false && "Quadrature type not available.");
        }
    }

///@}

/** @name Accessors
 *  Retrieve member variables through public member functions
 */
///@{

    KOKKOS_INLINE_FUNCTION 
    bool validQuadrature() const {
        return _is_initialized;
    }

    KOKKOS_INLINE_FUNCTION 
    int getNumberOfQuadraturePoints() const {
        return _number_of_quadrature_points;
    }

    KOKKOS_INLINE_FUNCTION 
    int getOrderOfQuadraturePoints() const {
        return _order_of_quadrature_points;
    }

    KOKKOS_INLINE_FUNCTION 
    int getDimensionOfQuadraturePoints() const {
        return _dimension_of_quadrature_points;
    }

    KOKKOS_INLINE_FUNCTION 
    QuadratureType getQuadratureType() const {
        return _qt;
    }

    decltype(_quadrature_weights) getWeights() const {
        return _quadrature_weights;
    }

    decltype(_parameterized_quadrature_sites) getSites() const {
        return _parameterized_quadrature_sites;
    }

    KOKKOS_INLINE_FUNCTION 
    double getWeight(const int index) const {
        return _quadrature_weights(index);
    }

    KOKKOS_INLINE_FUNCTION 
    double getSite(const int index, const int component) const {
        return _parameterized_quadrature_sites(index, component);
    }

///@}


/** @name Modifiers
 *  Changed member variables through public member functions
 */
///@{
///@}


}; // Quadrature Class
} // Compadre

#endif


