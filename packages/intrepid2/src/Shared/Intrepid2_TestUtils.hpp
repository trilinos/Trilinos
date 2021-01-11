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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_TestUtils.hpp
    \brief  Utility methods for Intrepid2 unit tests.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_TestUtils_h
#define Intrepid2_TestUtils_h

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_DerivedBasisFamily.hpp"
#include "Intrepid2_FunctorIterator.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_Sacado.hpp" // Sacado includes, guarded by the appropriate preprocessor variable
#include "Intrepid2_Utils.hpp"

#include "Teuchos_UnitTestHarness.hpp"

namespace Intrepid2
{
  //! Maximum number of derivatives to track for Fad types in tests.
  constexpr int MAX_FAD_DERIVATIVES_FOR_TESTS = 3;

  //! Use Teuchos small number determination on host; pass this to Intrepid2::relErr() on device.
  template <class Scalar>
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  smallNumber()
  {
    using ST = Teuchos::ScalarTraits<Scalar>;
    return ST::magnitude(Teuchos::RelErrSmallNumber<ST::hasMachineParameters,Scalar>::smallNumber());
  }

  //! Adapted from Teuchos::relErr(); for use in code that may be executed on device.
  template <class Scalar1, class Scalar2>
  KOKKOS_INLINE_FUNCTION
  bool
  relErrMeetsTol( const Scalar1 &s1, const Scalar2 &s2, const typename Teuchos::ScalarTraits< typename std::common_type<Scalar1,Scalar2>::type >::magnitudeType &smallNumber, const double &tol )
  {
    using Scalar        = typename std::common_type<Scalar1,Scalar2>::type;
    const Scalar s1Abs  = fabs(s1);
    const Scalar s2Abs  = fabs(s2);
    const Scalar maxAbs = (s1Abs > s2Abs) ? s1Abs : s2Abs;
    const Scalar relErr = fabs( s1 - s2 ) / ( smallNumber + maxAbs );
    return relErr < tol;
  }

  template <class Scalar1, class Scalar2>
  KOKKOS_INLINE_FUNCTION
  bool
  errMeetsAbsAndRelTol( const Scalar1 &s1, const Scalar2 &s2, const double &relTol, const double &absTol )
  {
    return fabs( s1 - s2 ) <= absTol + fabs(s1) * relTol;
  }

  static const double TEST_TOLERANCE_TIGHT = 1.e2 * std::numeric_limits<double>::epsilon();

  // we use DynRankView for both input points and values
  template<typename ScalarType>
  using ViewType = Kokkos::DynRankView<ScalarType,Kokkos::DefaultExecutionSpace>;

  template<typename ScalarType>
  using FixedRankViewType = Kokkos::View<ScalarType,Kokkos::DefaultExecutionSpace>;

  template<typename ScalarType>
  KOKKOS_INLINE_FUNCTION bool valuesAreSmall(const ScalarType &a, const ScalarType &b, const double &epsilon)
  {
    using std::abs;
    return (abs(a) < epsilon) && (abs(b) < epsilon);
  }

  inline bool approximatelyEqual(double a, double b, double epsilon)
  {
    const double larger_magnitude = (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a));
    return std::abs(a - b) <= larger_magnitude * epsilon;
  }

  inline bool essentiallyEqual(double a, double b, double epsilon)
  {
    const double smaller_magnitude = (std::abs(a) > std::abs(b) ? std::abs(b) : std::abs(a));
    return std::abs(a - b) <= smaller_magnitude * epsilon;
  }

  // conversion from the ref element [0,1] (used by ESEAS) to the ref element [-1,1] (used by Intrepid2)
  KOKKOS_INLINE_FUNCTION double fromZeroOne(double x_zero_one)
  {
    return x_zero_one * 2.0 - 1.0;
  }

  // conversion from the ref element [-1,1] (used by Intrepid2) to the ref element [0,1] (used by ESEAS)
  KOKKOS_INLINE_FUNCTION double toZeroOne(double x_minus_one_one)
  {
    return (x_minus_one_one + 1.0) / 2.0;
  }

  // conversion from the ref element [0,1] (used by ESEAS) to the ref element [-1,1] (used by Intrepid2)
  KOKKOS_INLINE_FUNCTION double fromZeroOne_dx(double dx_zero_one)
  {
    return dx_zero_one / 2.0;
  }

  // conversion from the ref element [-1,1] (used by Intrepid2) to the ref element [0,1] (used by ESEAS)
  KOKKOS_INLINE_FUNCTION double toZeroOne_dx(double dx_minus_one_one)
  {
    return dx_minus_one_one * 2.0;
  }

  template<class DeviceViewType>
  typename DeviceViewType::HostMirror getHostCopy(const DeviceViewType &deviceView)
  {
    typename DeviceViewType::HostMirror hostView = Kokkos::create_mirror(deviceView);
    Kokkos::deep_copy(hostView, deviceView);
    return hostView;
  }

  template<class BasisFamily>
  inline Teuchos::RCP< Intrepid2::Basis<Kokkos::DefaultExecutionSpace,double,double> > getBasisUsingFamily(shards::CellTopology cellTopo, Intrepid2::EFunctionSpace fs,
                                                                                                           int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1)
  {
    using BasisPtr = typename BasisFamily::BasisPtr;
    
    BasisPtr basis;
    using namespace Intrepid2;
    
    if (cellTopo.getBaseKey() == shards::Line<>::key)
    {
      basis = getLineBasis<BasisFamily>(fs,polyOrder_x);
    }
    else if (cellTopo.getBaseKey() == shards::Quadrilateral<>::key)
    {
      INTREPID2_TEST_FOR_EXCEPTION(polyOrder_y < 0, std::invalid_argument, "polyOrder_y must be specified");
      basis = getQuadrilateralBasis<BasisFamily>(fs,polyOrder_x,polyOrder_y);
    }
    else if (cellTopo.getBaseKey() == shards::Triangle<>::key)
    {
      basis = getTriangleBasis<BasisFamily>(fs,polyOrder_x);
    }
    else if (cellTopo.getBaseKey() == shards::Hexahedron<>::key)
    {
      INTREPID2_TEST_FOR_EXCEPTION(polyOrder_y < 0, std::invalid_argument, "polyOrder_y must be specified");
      INTREPID2_TEST_FOR_EXCEPTION(polyOrder_z < 0, std::invalid_argument, "polyOrder_z must be specified");
      basis = getHexahedronBasis<BasisFamily>(fs,polyOrder_x,polyOrder_y,polyOrder_z);
    }
    else if (cellTopo.getBaseKey() == shards::Tetrahedron<>::key)
    {
      basis = getTetrahedronBasis<BasisFamily>(fs, polyOrder_x);
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported cell topology");
    }
    return basis;
  }

  template<bool defineVertexFunctions>
  inline Teuchos::RCP< Intrepid2::Basis<Kokkos::DefaultExecutionSpace,double,double> > getHierarchicalBasis(shards::CellTopology cellTopo, Intrepid2::EFunctionSpace fs,
                                                                                                            int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1)
  {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using Scalar = double;
    using namespace Intrepid2;
    
    using LineBasisGrad = Intrepid2::IntegratedLegendreBasis_HGRAD_LINE<ExecSpace, Scalar, Scalar, defineVertexFunctions, true>;
    using LineBasisVol  = Intrepid2::LegendreBasis_HVOL_LINE< ExecSpace, Scalar, Scalar>;
    using TriangleBasisFamily = Intrepid2::HierarchicalTriangleBasisFamily<ExecSpace, Scalar, Scalar, defineVertexFunctions>;
    using TetrahedronBasisFamily = Intrepid2::HierarchicalTetrahedronBasisFamily<ExecSpace, Scalar, Scalar, defineVertexFunctions>;
    
    using BasisFamily = DerivedBasisFamily<LineBasisGrad, LineBasisVol, TriangleBasisFamily, TetrahedronBasisFamily>;
    
    return getBasisUsingFamily<BasisFamily>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
  }

  template<typename ValueType, class ... DimArgs>
  inline ViewType<ValueType> getView(const std::string &label, DimArgs... dims)
  {
    const bool allocateFadStorage = !std::is_pod<ValueType>::value;
    if (!allocateFadStorage)
    {
      return ViewType<ValueType>(label,dims...);
    }
    else
    {
      return ViewType<ValueType>(label,dims...,MAX_FAD_DERIVATIVES_FOR_TESTS+1);
    }
  }

  template<typename ValueType, class ... DimArgs>
  inline FixedRankViewType< typename RankExpander<ValueType, sizeof...(DimArgs) >::value_type > getFixedRankView(const std::string &label, DimArgs... dims)
  {
    const bool allocateFadStorage = !std::is_pod<ValueType>::value;
    using value_type = typename RankExpander<ValueType, sizeof...(dims) >::value_type;
    if (!allocateFadStorage)
    {
      return FixedRankViewType<value_type>(label,dims...);
    }
    else
    {
      return FixedRankViewType<value_type>(label,dims...,MAX_FAD_DERIVATIVES_FOR_TESTS+1);
    }
  }

/** \brief Returns a DynRankView containing regularly-spaced points on the specified cell topology.
    \param [in] cellTopo - the cell topology on which the points will be defined.
    \param [in] numPoints_1D - the number of points that will be defined along each edge.

The total number of points defined will be a triangular number; if n=numPointsBase, then the point count is the nth triangular number, given by n*(n+1)/2.
*/
  template <typename PointValueType>
  inline ViewType<PointValueType> getInputPointsView(shards::CellTopology &cellTopo, int numPoints_1D)
  {
    const ordinal_type order = numPoints_1D - 1;
    ordinal_type numPoints = PointTools::getLatticeSize(cellTopo, order);
    ordinal_type spaceDim  = cellTopo.getDimension();
    
    ViewType<PointValueType> inputPoints = getView<PointValueType>("input points",numPoints,spaceDim);
    PointTools::getLattice(inputPoints, cellTopo, order, 0, POINTTYPE_EQUISPACED );
    
    return inputPoints;
  }

  template<typename OutputValueType>
  inline ViewType<OutputValueType> getOutputView(Intrepid2::EFunctionSpace fs, Intrepid2::EOperator op, int basisCardinality, int numPoints, int spaceDim)
  {
    switch (fs) {
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        switch (op) {
          case Intrepid2::OPERATOR_VALUE:
            return getView<OutputValueType>("H^1 value output",basisCardinality,numPoints);
          case Intrepid2::OPERATOR_GRAD:
            return getView<OutputValueType>("H^1 derivative output",basisCardinality,numPoints,spaceDim);
          case Intrepid2::OPERATOR_D1:
          case Intrepid2::OPERATOR_D2:
          case Intrepid2::OPERATOR_D3:
          case Intrepid2::OPERATOR_D4:
          case Intrepid2::OPERATOR_D5:
          case Intrepid2::OPERATOR_D6:
          case Intrepid2::OPERATOR_D7:
          case Intrepid2::OPERATOR_D8:
          case Intrepid2::OPERATOR_D9:
          case Intrepid2::OPERATOR_D10:
          {
            const auto dkcard = getDkCardinality(op, spaceDim);
            return getView<OutputValueType>("H^1 derivative output",basisCardinality,numPoints,dkcard);
          }
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
        }
      case Intrepid2::FUNCTION_SPACE_HCURL:
        switch (op) {
          case Intrepid2::OPERATOR_VALUE:
            return getView<OutputValueType>("H(curl) value output",basisCardinality,numPoints,spaceDim);
          case Intrepid2::OPERATOR_CURL:
            if (spaceDim == 2)
              return getView<OutputValueType>("H(curl) derivative output",basisCardinality,numPoints);
            else if (spaceDim == 3)
              return getView<OutputValueType>("H(curl) derivative output",basisCardinality,numPoints,spaceDim);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
        }
      case Intrepid2::FUNCTION_SPACE_HDIV:
        switch (op) {
          case Intrepid2::OPERATOR_VALUE:
            return getView<OutputValueType>("H(div) value output",basisCardinality,numPoints,spaceDim);
          case Intrepid2::OPERATOR_DIV:
            return getView<OutputValueType>("H(div) derivative output",basisCardinality,numPoints);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
        }
        
      case Intrepid2::FUNCTION_SPACE_HVOL:
        switch (op) {
          case Intrepid2::OPERATOR_VALUE:
            return getView<OutputValueType>("H(vol) value output",basisCardinality,numPoints);
          case Intrepid2::OPERATOR_D1:
          case Intrepid2::OPERATOR_D2:
          case Intrepid2::OPERATOR_D3:
          case Intrepid2::OPERATOR_D4:
          case Intrepid2::OPERATOR_D5:
          case Intrepid2::OPERATOR_D6:
          case Intrepid2::OPERATOR_D7:
          case Intrepid2::OPERATOR_D8:
          case Intrepid2::OPERATOR_D9:
          case Intrepid2::OPERATOR_D10:
          {
            const auto dkcard = getDkCardinality(op, spaceDim);
            return getView<OutputValueType>("H(vol) derivative output",basisCardinality,numPoints,dkcard);
          }
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
        }
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
    }
  }

  // ! This returns a vector whose entries are vector<int>s containing 1-3 polynomial orders from 1 up to and including those specified
  // ! Intended for testing bases that support anisotropic polynomial degree, such as the hierarchical bases
  inline std::vector< std::vector<int> > getBasisTestCasesUpToDegree(int spaceDim, int minDegree, int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1)
  {
    std::vector<int> degrees(spaceDim);
    degrees[0] = polyOrder_x;
    if (spaceDim > 1) degrees[1] = polyOrder_y;
    if (spaceDim > 2) degrees[2] = polyOrder_z;
    
    int numCases = degrees[0];
    for (unsigned d=1; d<degrees.size(); d++)
    {
      numCases = numCases * (degrees[d] + 1 - minDegree);
    }
    std::vector< std::vector<int> > subBasisDegreeTestCases(numCases);
    for (int caseOrdinal=0; caseOrdinal<numCases; caseOrdinal++)
    {
      std::vector<int> subBasisDegrees(degrees.size());
      int caseRemainder = caseOrdinal;
      for (int d=degrees.size()-1; d>=0; d--)
      {
        int subBasisDegree = caseRemainder % (degrees[d] + 1 - minDegree);
        caseRemainder = caseRemainder / (degrees[d] + 1 - minDegree);
        subBasisDegrees[d] = subBasisDegree + minDegree;
      }
      subBasisDegreeTestCases[caseOrdinal] = subBasisDegrees;
    }
    return subBasisDegreeTestCases;
  }

  //! Copy the values for the specified functor
  template<class Functor, class Scalar, int rank>
  typename ViewType<Scalar>::HostMirror copyFunctorToHostView(const Functor &deviceFunctor)
  {
    INTREPID2_TEST_FOR_EXCEPTION(rank != getFunctorRank(deviceFunctor), std::invalid_argument, "functor rank must match the template argument");
    
    ViewType<Scalar> view;
    const std::string label = "functor copy";
    const auto &f = deviceFunctor;
    switch (rank)
    {
      case 0:
        view = getView<Scalar>(label);
        break;
      case 1:
        view = getView<Scalar>(label, f.extent_int(0));
        break;
      case 2:
        view = getView<Scalar>(label, f.extent_int(0), f.extent_int(1));
        break;
      case 3:
        view = getView<Scalar>(label, f.extent_int(0), f.extent_int(1), f.extent_int(2));
        break;
      case 4:
        view = getView<Scalar>(label, f.extent_int(0), f.extent_int(1), f.extent_int(2), f.extent_int(3));
        break;
      case 5:
        view = getView<Scalar>(label, f.extent_int(0), f.extent_int(1), f.extent_int(2), f.extent_int(3), f.extent_int(4));
        break;
      case 6:
        view = getView<Scalar>(label, f.extent_int(0), f.extent_int(1), f.extent_int(2), f.extent_int(3), f.extent_int(4), f.extent_int(5));
        break;
      case 7:
        view = getView<Scalar>(label, f.extent_int(0), f.extent_int(1), f.extent_int(2), f.extent_int(3), f.extent_int(4), f.extent_int(5), f.extent_int(6));
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported functor rank");
    }
    
    int entryCount = view.size();
    
    using ExecutionSpace = typename ViewType<Scalar>::execution_space;
    
    using ViewIteratorScalar    = Intrepid2::ViewIterator<ViewType<Scalar>, Scalar>;
    using FunctorIteratorScalar = FunctorIterator<Functor, Scalar, rank>;
    
    Kokkos::RangePolicy < ExecutionSpace > policy(0,entryCount);
    Kokkos::parallel_for( policy,
    KOKKOS_LAMBDA (const int &enumerationIndex )
    {
      ViewIteratorScalar     vi(view);
      FunctorIteratorScalar  fi(f);
      
      vi.setEnumerationIndex(enumerationIndex);
      fi.setEnumerationIndex(enumerationIndex);
      
      vi.set(fi.get());
    }
    );
    
    auto hostView = Kokkos::create_mirror_view(view);
    Kokkos::deep_copy(hostView, view);
    return hostView;
  }

  template<class FunctorType, typename Scalar, int rank>
  void printFunctor(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    auto functorHostCopy = copyFunctorToHostView<FunctorType, Scalar, rank>(functor);
    
    std::string name = (functorName == "") ? "Functor" : functorName;
    
    out << "\n******** " << name << " (rank " << rank << ") ********\n";
    out << "dimensions: (";
    for (int r=0; r<rank; r++)
    {
      out << functor.extent_int(r);
      if (r<rank-1) out << ",";
    }
    out << ")\n";
    
    ViewIterator<decltype(functorHostCopy),Scalar> vi(functorHostCopy);
    
    bool moreEntries = true;
    while (moreEntries)
    {
      Scalar value = vi.get();
      
      auto location = vi.getLocation();
      out << functorName << "(";
      for (ordinal_type i=0; i<rank; i++)
      {
        out << location[i];
        if (i<rank-1)
        {
          out << ",";
        }
      }
      out << ") " << value << std::endl;

      moreEntries = (vi.increment() != -1);
    }
    out << "\n";
  }
  
  template<class FunctorType>
  void printFunctor1(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    using Scalar = typename std::remove_reference<decltype(functor(0))>::type;
    printFunctor<FunctorType, Scalar, 1>(functor, out, functorName);
  }

  template<class FunctorType>
  void printFunctor2(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    using Scalar = typename std::remove_const<typename std::remove_reference<decltype(functor(0,0))>::type>::type;
    printFunctor<FunctorType, Scalar, 2>(functor, out, functorName);
  }
    
  template<class FunctorType>
  void printFunctor3(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    using Scalar = typename std::remove_const<typename std::remove_reference<decltype(functor(0,0,0))>::type>::type;
    printFunctor<FunctorType, Scalar, 3>(functor, out, functorName);
  }
    
  template<class FunctorType>
  void printFunctor4(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    using Scalar = typename std::remove_const<typename std::remove_reference<decltype(functor(0,0,0,0))>::type>::type;
    printFunctor<FunctorType, Scalar, 4>(functor, out, functorName);
  }
    
  template<class FunctorType>
  void printFunctor5(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    using Scalar = typename std::remove_const<typename std::remove_reference<decltype(functor(0,0,0,0,0))>::type>::type;
    printFunctor<FunctorType, Scalar, 5>(functor, out, functorName);
  }

  template<class FunctorType>
  void printFunctor6(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    using Scalar = typename std::remove_const<typename std::remove_reference<decltype(functor(0,0,0,0,0,0))>::type>::type;
    printFunctor<FunctorType, Scalar, 6>(functor, out, functorName);
  }

  template<class FunctorType>
  void printFunctor7(const FunctorType &functor, std::ostream &out, const std::string &functorName = "")
  {
    using Scalar = typename std::remove_const<typename std::remove_reference<decltype(functor(0,0,0,0,0,0,0))>::type>::type;
    printFunctor<FunctorType, Scalar, 7>(functor, out, functorName);
  }

  template<class View>
  void printView(const View &view, std::ostream &out, const std::string &viewName = "")
  {
    using Scalar   = typename View::value_type;
    using HostView = typename View::HostMirror;
    using HostViewIteratorScalar = Intrepid2::ViewIterator<HostView, Scalar>;
    
    auto hostView = getHostCopy(view);
    
    HostViewIteratorScalar vi(hostView);
    
    bool moreEntries = (vi.nextIncrementRank() != -1);
    while (moreEntries)
    {
      Scalar value = vi.get();
      
      auto location = vi.getLocation();
      out << viewName << "(";
      for (unsigned i=0; i<getFunctorRank(view); i++)
      {
        out << location[i];
        if (i<getFunctorRank(view)-1)
        {
          out << ",";
        }
      }
      out << ") " << value << std::endl;

      moreEntries = (vi.increment() != -1);
    }
  }

  template <class FunctorType1, class FunctorType2, int rank, typename Scalar=typename FunctorType1::value_type, class ExecutionSpace = typename FunctorType1::execution_space>
  typename std::enable_if< !supports_rank<FunctorType1,rank>::value, void >::type
  testFloatingEquality(const FunctorType1 &functor1, const FunctorType2 &functor2, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                       std::string functor1Name = "Functor 1", std::string functor2Name = "Functor 2")
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "testFloatingEquality() called on FunctorType1 that does not support the specified rank.\n");
  }

  //! this method assumes both functors are accesible on device.
  template <class FunctorType1, class FunctorType2, int rank, typename Scalar=typename FunctorType1::value_type, class ExecutionSpace = typename FunctorType1::execution_space>
  typename std::enable_if< supports_rank<FunctorType1,rank>::value, void >::type
  testFloatingEquality(const FunctorType1 &functor1, const FunctorType2 &functor2, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                       std::string functor1Name = "Functor 1", std::string functor2Name = "Functor 2")
  {
    static_assert( supports_rank<FunctorType2,rank>::value, "Both Functor1 and Functor2 must support the specified rank through operator().");
    using Functor1IteratorScalar = FunctorIterator<FunctorType1, Scalar, rank>;
    using Functor2IteratorScalar = FunctorIterator<FunctorType2, Scalar, rank>;
    using ViewIteratorBool       = Intrepid2::ViewIterator<ViewType<bool>, bool>;

    // check that rank/size match
    TEUCHOS_TEST_FOR_EXCEPTION(getFunctorRank(functor1) != rank, std::invalid_argument, "functor1 and functor2 must agree in rank"); // Kokkos::View does not provide rank() method; getFunctorRank() provides common interface
    TEUCHOS_TEST_FOR_EXCEPTION(getFunctorRank(functor2) != rank, std::invalid_argument, "functor1 and functor2 must agree in rank"); // Kokkos::View does not provide rank() method; getFunctorRank() provides common interface
    
    int entryCount = 1;
    for (unsigned i=0; i<getFunctorRank(functor1); i++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(functor1.extent_int(i) != functor2.extent_int(i), std::invalid_argument,
                                 "functor1 and functor2 must agree in size in each dimension; functor1 has extent "
                                 + std::to_string(functor1.extent_int(i)) + " in dimension " + std::to_string(i)
                                 + "; functor2 has extent " + std::to_string(functor2.extent_int(i)) );
      entryCount *= functor1.extent_int(i);
    }
    if (entryCount == 0) return; // nothing to test
    
    ViewType<bool> valuesMatch = getView<bool>("valuesMatch", entryCount);
    
    Kokkos::RangePolicy < ExecutionSpace > policy(0,entryCount);
    Kokkos::parallel_for( policy,
    KOKKOS_LAMBDA (const int &enumerationIndex )
    {
      Functor1IteratorScalar vi1(functor1);
      Functor2IteratorScalar vi2(functor2);

      vi1.setEnumerationIndex(enumerationIndex);
      vi2.setEnumerationIndex(enumerationIndex);
      
      const Scalar & value1 = vi1.get();
      const Scalar & value2 = vi2.get();
      
      bool errMeetsTol = errMeetsAbsAndRelTol(value1, value2, relTol, absTol);
      valuesMatch(enumerationIndex) = errMeetsTol;
    }
    );

    int failureCount = 0; 
    Kokkos::RangePolicy<ExecutionSpace > reducePolicy(0, entryCount);
    Kokkos::parallel_reduce( reducePolicy,
    KOKKOS_LAMBDA( const int &enumerationIndex, int &reducedValue )
    {
      reducedValue += valuesMatch(enumerationIndex) ? 0 : 1;
    }, failureCount);
        
    const bool allValuesMatch = (failureCount == 0);
    success = success && allValuesMatch;
    
    if (!allValuesMatch)
    {
      // copy functors to host views
      auto functor1CopyHost = copyFunctorToHostView<FunctorType1,Scalar,rank>(functor1);
      auto functor2CopyHost = copyFunctorToHostView<FunctorType2,Scalar,rank>(functor2);
      
      auto valuesMatchHost = getHostCopy(valuesMatch);
      
      FunctorIterator<decltype(functor1CopyHost),Scalar,rank> vi1(functor1CopyHost);
      FunctorIterator<decltype(functor2CopyHost),Scalar,rank> vi2(functor2CopyHost);
      ViewIteratorBool    viMatch(valuesMatchHost);
      
      bool moreEntries = true;
      while (moreEntries)
      {
        bool errMeetsTol = viMatch.get();
        
        if (!errMeetsTol)
        {
          const Scalar value1 = vi1.get();
          const Scalar value2 = vi2.get();
          
          success = false;
          auto location = vi1.getLocation();
          out << "At location (";
          for (unsigned i=0; i<getFunctorRank(functor1); i++)
          {
            out << location[i];
            if (i<getFunctorRank(functor1)-1)
            {
              out << ",";
            }
          }
          out << ") " << functor1Name << " value != " << functor2Name << " value (";
          out << value1 << " != " << value2 << "); diff is " << std::abs(value1-value2) << std::endl;
        }

        moreEntries =                (vi1.increment() != -1);
        moreEntries = moreEntries && (vi2.increment() != -1);
        moreEntries = moreEntries && (viMatch.increment() != -1);
      }
    }
  }

  template <class ViewType, class FunctorType>
  void testFloatingEquality1(const ViewType &view, const FunctorType &functor, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                             std::string view1Name = "View", std::string view2Name = "Functor")
  {
    testFloatingEquality<ViewType, FunctorType, 1>(view, functor, relTol, absTol, out, success, view1Name, view2Name);
  }

  template <class ViewType, class FunctorType>
  void testFloatingEquality2(const ViewType &view, const FunctorType &functor, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                             std::string view1Name = "View", std::string view2Name = "Functor")
  {
    testFloatingEquality<ViewType, FunctorType, 2>(view, functor, relTol, absTol, out, success, view1Name, view2Name);
  }
    
  template <class ViewType, class FunctorType>
  void testFloatingEquality3(const ViewType &view, const FunctorType &functor, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                             std::string view1Name = "View", std::string view2Name = "Functor")
  {
    testFloatingEquality<ViewType, FunctorType, 3>(view, functor, relTol, absTol, out, success, view1Name, view2Name);
  }
    
  template <class ViewType, class FunctorType>
  void testFloatingEquality4(const ViewType &view, const FunctorType &functor, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                             std::string view1Name = "View", std::string view2Name = "Functor")
  {
    testFloatingEquality<ViewType, FunctorType, 4>(view, functor, relTol, absTol, out, success, view1Name, view2Name);
  }

  template <class ViewType, class FunctorType>
  void testFloatingEquality5(const ViewType &view, const FunctorType &functor, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                             std::string view1Name = "View", std::string view2Name = "Functor")
  {
    testFloatingEquality<ViewType, FunctorType, 5>(view, functor, relTol, absTol, out, success, view1Name, view2Name);
  }

  template <class ViewType, class FunctorType>
  void testFloatingEquality6(const ViewType &view, const FunctorType &functor, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                             std::string view1Name = "View", std::string view2Name = "Functor")
  {
    testFloatingEquality<ViewType, FunctorType, 6>(view, functor, relTol, absTol, out, success, view1Name, view2Name);
  }

  template <class ViewType, class FunctorType>
  void testFloatingEquality7(const ViewType &view, const FunctorType &functor, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                             std::string view1Name = "View", std::string view2Name = "Functor")
  {
    testFloatingEquality<ViewType, FunctorType, 7>(view, functor, relTol, absTol, out, success, view1Name, view2Name);
  }

  template <class ViewType1, class ViewType2>
  void testViewFloatingEquality(const ViewType1 &view1, const ViewType2 &view2, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                                std::string view1Name = "View 1", std::string view2Name = "View 2")
  {
    // check that rank/size match
    TEUCHOS_TEST_FOR_EXCEPTION(view1.rank() != view2.rank(), std::invalid_argument, "views must agree in rank");
    for (unsigned i=0; i<view1.rank(); i++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(view1.extent_int(i) != view2.extent_int(i), std::invalid_argument, "views must agree in size in each dimension");
    }
    
    if (view1.size() == 0) return; // nothing to test

    const int rank = view1.rank();
    switch (rank)
    {
      case 0: testFloatingEquality<ViewType1, ViewType2, 0>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      case 1: testFloatingEquality<ViewType1, ViewType2, 1>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      case 2: testFloatingEquality<ViewType1, ViewType2, 2>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      case 3: testFloatingEquality<ViewType1, ViewType2, 3>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      case 4: testFloatingEquality<ViewType1, ViewType2, 4>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      case 5: testFloatingEquality<ViewType1, ViewType2, 5>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      case 6: testFloatingEquality<ViewType1, ViewType2, 6>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      case 7: testFloatingEquality<ViewType1, ViewType2, 7>(view1, view2, relTol, absTol, out, success, view1Name, view2Name);  break;
      default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported rank");
    }
  }

} // namespace Intrepid2

#ifdef HAVE_INTREPID2_SACADO
using Sacado_Fad_DFadType = Sacado::Fad::DFad<double>;
#define INTREPID2_POINTSCALAR_TEST_INSTANT(GROUP_NAME, TEST_NAME) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( GROUP_NAME, TEST_NAME, double ) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( GROUP_NAME, TEST_NAME, Sacado_Fad_DFadType ) \

#define INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT(GROUP_NAME, TEST_NAME) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( GROUP_NAME, TEST_NAME, double, double ) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( GROUP_NAME, TEST_NAME, Sacado_Fad_DFadType, Sacado_Fad_DFadType ) \

#else
#define INTREPID2_POINTSCALAR_TEST_INSTANT(GROUP_NAME, TEST_NAME) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( GROUP_NAME, TEST_NAME, double ) \

#define INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT(GROUP_NAME, TEST_NAME) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( GROUP_NAME, TEST_NAME, double, double ) \

#endif

#endif /* Intrepid2_TestUtils_h */
