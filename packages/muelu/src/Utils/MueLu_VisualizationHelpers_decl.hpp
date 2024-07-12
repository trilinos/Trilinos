// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_VISUALIZATIONHELPERS_DECL_HPP_
#define MUELU_VISUALIZATIONHELPERS_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VisualizationHelpers_fwd.hpp"
#include "MueLu_LWGraph.hpp"

#include <list>

namespace MueLu {

class Level;
// Utility classes used in convex hull algorithm

class myTriangle {
 public:
  myTriangle()
    : v1(0)
    , v2(0)
    , v3(0) {}
  myTriangle(int v1in, int v2in, int v3in)
    : v1(v1in)
    , v2(v2in)
    , v3(v3in) {}
  ~myTriangle() {}
  bool operator==(const myTriangle& l) {
    if (l.v1 == v1 && l.v2 == v2 && l.v3 == v3)
      return true;
    return false;
  }
  int v1;
  int v2;
  int v3;
};

class myVec3 {
 public:
  myVec3()
    : x(0)
    , y(0)
    , z(0) {}
  myVec3(double xin, double yin, double zin)
    : x(xin)
    , y(yin)
    , z(zin) {}
  ~myVec3() {}
  double x;
  double y;
  double z;
};

class myVec2 {
 public:
  myVec2()
    : x(0)
    , y(0) {}
  myVec2(double xin, double yin)
    : x(xin)
    , y(yin) {}
  ~myVec2() {}
  double x;
  double y;
};

/*!
  @class VisualizationHelpers class.
  @brief Base class providing routines to visualize aggregates and coarsening information.

  This class is the base class for the CoarseningVisualizationFactory as well as the AggregationExporterFactory to
  visualize aggregates or coarsening information from the transfer operators.

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class VisualizationHelpers {
#undef MUELU_VISUALIZATIONHELPERS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  VisualizationHelpers() {}

  //! Destructor.
  virtual ~VisualizationHelpers() {}
  //@}

  RCP<ParameterList> GetValidParameterList() const;

 protected:
  void writeFileVTKOpening(std::ofstream& fout, std::vector<int>& uniqueFine, std::vector<int>& geomSizesFine) const;
  void writeFileVTKNodes(std::ofstream& fout, std::vector<int>& uniqueFine, Teuchos::RCP<const Map>& nodeMap) const;
  void writeFileVTKData(std::ofstream& fout, std::vector<int>& uniqueFine, LocalOrdinal myAggOffset, ArrayRCP<LocalOrdinal>& vertex2AggId, int myRank) const;
  void writeFileVTKCoordinates(std::ofstream& fout, std::vector<int>& uniqueFine, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fx, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fy, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fz, int dim) const;
  void writeFileVTKCells(std::ofstream& fout, std::vector<int>& uniqueFine, std::vector<LocalOrdinal>& vertices, std::vector<LocalOrdinal>& geomSize) const;
  void writeFileVTKClosing(std::ofstream& fout) const;
  void writePVTU(std::ofstream& pvtu, std::string baseFname, int numProcs, bool bFineEdges = false, bool bCoarseEdges = false) const;
  void buildColormap() const;

  std::string getFileName(int numProcs, int myRank, int level, const Teuchos::ParameterList& pL) const;
  std::string getBaseFileName(int numProcs, int level, const Teuchos::ParameterList& pL) const;
  std::string getPVTUFileName(int numProcs, int myRank, int level, const Teuchos::ParameterList& pL) const;

  // move these routines to a common base class for visualization factories?
  static void doPointCloud(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes);
  static void doJacks(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& isRoot, const ArrayRCP<LO>& vertex2AggId);
  static void doConvexHulls2D(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& isRoot, const ArrayRCP<LO>& vertex2AggId, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords);
  static void doConvexHulls3D(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& isRoot, const ArrayRCP<LO>& vertex2AggId, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& zCoords);

  static void doGraphEdges(std::vector<int>& vertices, std::vector<int>& geomSizes, Teuchos::RCP<LWGraph>& G, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fx, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fy, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fz);

  static int ccw(const myVec2& a, const myVec2& b, const myVec2& c);
  static myVec3 crossProduct(myVec3 v1, myVec3 v2);
  static double dotProduct(myVec2 v1, myVec2 v2);
  static double dotProduct(myVec3 v1, myVec3 v2);
  static bool isInFront(myVec3 point, myVec3 inPlane, myVec3 n);
  static double mymagnitude(myVec2 vec);
  static double mymagnitude(myVec3 vec);
  static double distance(myVec2 p1, myVec2 p2);
  static double distance(myVec3 p1, myVec3 p2);
  static myVec2 vecSubtract(myVec2 v1, myVec2 v2);
  static myVec3 vecSubtract(myVec3 v1, myVec3 v2);
  static myVec2 getNorm(myVec2 v);
  static myVec3 getNorm(myVec3 v1, myVec3 v2, myVec3 v3);
  static double pointDistFromTri(myVec3 point, myVec3 v1, myVec3 v2, myVec3 v3);
  static std::vector<myTriangle> processTriangle(std::list<myTriangle>& tris, myTriangle tri, std::list<int>& pointsInFront, myVec3& barycenter, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& zCoords);
  static std::vector<int> giftWrap(std::vector<myVec2>& points, std::vector<int>& nodes, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords);

  std::string replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const;
  std::vector<int> makeUnique(std::vector<int>& vertices) const;  //!< replaces node indices in vertices with compressed unique indices, and returns list of unique points
};                                                                // class VisualizationHelpers
}  // namespace MueLu

#define MUELU_VISUALIZATIONHELPERS_SHORT

#endif /* MUELU_VISUALIZATIONHELPERS_DECL_HPP_ */
