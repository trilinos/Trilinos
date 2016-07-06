// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureDirectLineGaussJacobi20Def.hpp
    \brief  Definition file for the Intrepid2::CubatureDirectLineGaussJacobi20 class.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
            Kokkorized by Kyungjoo Kim
*/

namespace Intrepid2 {

  template <typename SpT, typename PT, typename WT>
  CubatureDirectLineGaussJacobi20<SpT,PT,WT>::
  CubatureDirectLineGaussJacobi20(const ordinal_type degree)
    : CubatureDirect<SpT>(degree, 1) {
    
    INTREPID2_TEST_FOR_EXCEPTION( degree < 0 ||
                                  degree > static_cast<ordinal_type>(Parameters::MaxCubatureDegreePyr), std::out_of_range,
                                  ">>> ERROR (CubatureDirectLineGaussJacobi20): No cubature rule implemented for the desired polynomial degree.");

    typedef Kokkos::DynRankView<PT,Kokkos::LayoutRight,Kokkos::HostSpace> pointViewHostType;
    typedef Kokkos::DynRankView<WT,Kokkos::LayoutRight,Kokkos::HostSpace> weightViewHostType;
    
    this->cubatureData_.numPoints_ = cubatureDataStatic_[this->degree_].numPoints_;
    const Kokkos::pair<ordinal_type,ordinal_type> pointRange(0, this->cubatureData_.numPoints_);
    {
      // src
      const pointViewHostType points_host(const_cast<PT*>( &(cubatureDataStatic_[this->degree_].points_[0][0]) ),
                                          pointRange.second,
                                          Parameters::MaxDimension);

      auto points = Kokkos::create_mirror_view(typename SpT::memory_space(), points_host);

      Kokkos::deep_copy(points,points_host);

      // dst
      this->cubatureData_.points_ = Kokkos::DynRankView<PT,SpT>("CubatureDirectLineGaussJacobi20::cubatureData_::points_",
                                                                pointRange.second,
                                                                Parameters::MaxDimension);
      // copy
      Kokkos::deep_copy(this->cubatureData_.points_ , points );
    }
    {
      // src
      const weightViewHostType weights(const_cast<PT*>( &(cubatureDataStatic_[this->degree_].weights_[0]) ),
                                       pointRange.second);

      // dst
      this->cubatureData_.weights_ = Kokkos::DynRankView<WT,SpT>("CubatureDirectLineGaussJacobi20::cubatureData_::weights_",
                                                                 pointRange.second);
      // copy
      Kokkos::deep_copy(Kokkos::subdynrankview(this->cubatureData_.weights_, Kokkos::ALL()) , Kokkos::subdynrankview(weights, Kokkos::ALL()));
    }
  }


  //-------------------------------------------------------------------------------------//
  //                          Definition of cubature templates                           //
  //-------------------------------------------------------------------------------------//

  /*
    Cubature templates for lines are defined the reference cell:

    Line -> (-1,0,0),(1,0,0)
  */

  /*
    This static const member contains templates for GaussJacobi20(-Legendre) rules.
  */
  template<typename SpT, typename PT, typename WT>
  const typename CubatureDirect<SpT,PT,WT>::CubatureDataStatic
  CubatureDirectLineGaussJacobi20<SpT,PT,WT>::
  cubatureDataStatic_[cubatureDataStaticSize] = {
    // Collection of GaussJacobi20 rules on [-1,1]
    // The rule with array index i is exact for polynomials up to order i
    {
      1,
      {{-0.5, 0.0, 0.0}},
      {2.66666666666666666666666666}
    },
    {
      1,
      {{-0.5, 0.0, 0.0}},
      {2.66666666666666666666666666}
    },
    {
      2,
      {{-7.549703546891172e-1, 0.0, 0.0},
       {8.830368802245062e-2, 0.0, 0.0}},
      {1.860379610028064,
       8.062870566386037e-01}
    },
    {
      2,
      {{-7.549703546891172e-1, 0.0, 0.0},
       {8.830368802245062e-2, 0.0, 0.0}},
      {1.860379610028064,
       8.062870566386037e-01}
    },
    {
      3,
      {{-8.540119518537008e-01, 0.0, 0.0},
       {-3.059924679232963e-01, 0.0, 0.0},
       { 4.100044197769969e-01, 0.0, 0.0}},
      {1.257090888519093e+00,
       1.169970154078928e+00,
       2.396056240686456e-01}
    },
    {
      3,
      {{-8.540119518537008e-01, 0.0, 0.0},
       {-3.059924679232963e-01, 0.0, 0.0},
       { 4.100044197769969e-01, 0.0, 0.0}},
      {1.257090888519093e+00,
       1.169970154078928e+00,
       2.396056240686456e-01}
    },
    {
      4,
      {{-9.029989011060054e-01, 0.0, 0.0},
       {-5.227985248962754e-01, 0.0, 0.0},
       {3.409459020873505e-02, 0.0, 0.0},
       {5.917028357935457e-01, 0.0, 0.0}},
      {8.871073248902235e-01,
       1.147670318393715e+00,
       5.490710973833849e-01,
       8.281792599934450e-02}
    },
    {
      4,
      {{-9.029989011060054e-01, 0.0, 0.0},
       {-5.227985248962754e-01, 0.0, 0.0},
       {3.409459020873505e-02, 0.0, 0.0},
       {5.917028357935457e-01, 0.0, 0.0}},
      {8.871073248902235e-01,
       1.147670318393715e+00,
       5.490710973833849e-01,
       8.281792599934450e-02}
    },
    {
      5,
      {{-9.308421201635699e-01, 0.0, 0.0},
       {-6.530393584566085e-01, 0.0, 0.0},
       {-2.202272258689614e-01, 0.0, 0.0},
       {2.686669452617736e-01, 0.0, 0.0},
       {7.021084258940329e-01, 0.0, 0.0}},
      {6.541182742861678e-01,
       1.009591695199292e+00,
       7.136012897727201e-01,
       2.564448057836956e-01,
       3.291060162479211e-02}
    },
    {
      5,
      {{-9.308421201635699e-01, 0.0, 0.0},
       {-6.530393584566085e-01, 0.0, 0.0},
       {-2.202272258689614e-01, 0.0, 0.0},
       {2.686669452617736e-01, 0.0, 0.0},
       {7.021084258940329e-01, 0.0, 0.0}},
      {6.541182742861678e-01,
       1.009591695199292e+00,
       7.136012897727201e-01,
       2.564448057836956e-01,
       3.291060162479211e-02}
    },
    {
      6,
      {{-9.481908898126656e-01, 0.0, 0.0},
       {-7.368721166840297e-01, 0.0, 0.0},
       {-3.951261639542174e-01, 0.0, 0.0},
       {1.807282632950432e-02, 0.0, 0.0},
       {4.313622546234276e-01, 0.0, 0.0},
       {7.736112323551237e-01, 0.0, 0.0}},
      {5.003096218126469e-01,
       8.590119978942462e-01,
       7.566174939883307e-01,
       4.103165690369299e-01,
       1.257623774795603e-01,
       1.464860645495425e-02}
    },
    {
      6,
      {{-9.481908898126656e-01, 0.0, 0.0},
       {-7.368721166840297e-01, 0.0, 0.0},
       {-3.951261639542174e-01, 0.0, 0.0},
       {1.807282632950432e-02, 0.0, 0.0},
       {4.313622546234276e-01, 0.0, 0.0},
       {7.736112323551237e-01, 0.0, 0.0}},
      {5.003096218126469e-01,
       8.590119978942462e-01,
       7.566174939883307e-01,
       4.103165690369299e-01,
       1.257623774795603e-01,
       1.464860645495425e-02}
    } // end GaussJacobi20

  };
    
} // end namespace Intrepid2
