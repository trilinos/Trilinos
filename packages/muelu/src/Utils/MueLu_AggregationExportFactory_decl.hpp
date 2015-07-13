// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_AggregationExportFactory_decl.hpp
 *
 *  Created on: Feb 10, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_
#define MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_AggregationExportFactory_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

  class Level;
  //Utility classes used in convex hull algorithm

  /*!
    @class AggregationExportFactory class.
    @brief Factory for exporting aggregates data

  */
  class Triangle_
  {
    public:
      Triangle_() : v1(0), v2(0), v3(0) {}
      Triangle_(int v1in, int v2in, int v3in) : v1(v1in), v2(v2in), v3(v3in) {}
      ~Triangle_() {}
      bool operator==(const Triangle_& l)
      {
        if(l.v1 == v1 && l.v2 == v2 && l.v3 == v3)
          return true;
        return false;
      }
      int v1;
      int v2;
      int v3;
  };
  class vec3_
  {
    public:
      vec3_() : x(0), y(0), z(0) {}
      vec3_(double xin, double yin, double zin) : x(xin), y(yin), z(zin) {}
      ~vec3_() {}
      double x;
      double y;
      double z;
  };
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class AggregationExportFactory : public TwoLevelFactoryBase {
#undef MUELU_AGGREGATIONEXPORTFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    AggregationExportFactory() { }

    //! Destructor.
    virtual ~AggregationExportFactory() { }
    //@}

    RCP<const ParameterList> GetValidParameterList() const;

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level &fineLevel, Level &coarseLevel) const;

    //@}

  private:
    std::string replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const;
    //Break different viz styles into separate functions for organization:
    void doPointCloud_(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
    void doJacks_(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
    void doJacksPlus_(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
    void doConvexHulls_(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
    void doConvexHulls2D_(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
    void doConvexHulls3D_(std::vector<int>& vertices, std::vector<int>& geomSizes) const;
    void doAlphaHulls_(std::vector<int>& vertices, std::vector<int>& geomSizes) const; //not implemented yet
    void doGraphEdges_(std::ofstream& fout, Teuchos::RCP<Matrix>& A, Teuchos::RCP<GraphBase>& G, bool fine) const; //add geometry to display node connections from a matrix. Connections in graph but not matrix have different color.
    void writeFile_(std::ofstream& fout, std::string styleName, std::vector<int>& vertices, std::vector<int>& geomSizes, std::vector<int>& verticesCoarse, std::vector<int>& geomSizesCoarse) const; //write the local .vtu file with the computed geometry
    void buildColormap_() const;
    void writePVTU_(std::ofstream& pvtu, std::string baseFname, int numProcs) const;
    std::vector<int> makeUnique_(std::vector<int>& vertices) const; //replaces node indices in vertices with compressed unique indices, and returns list of unique points
    //Utility functions for convex hulls
    static vec3_ crossProduct_(vec3_ v1, vec3_ v2);
    static double dotProduct_(vec3_ v1, vec3_ v2);
    static bool isInFront_(vec3_ point, vec3_ inPlane, vec3_ n);
    static double magnitude_(vec3_ vec);
    static double distance_(vec3_ p1, vec3_ p2);
    static vec3_ vecSubtract_(vec3_ v1, vec3_ v2);
    static vec3_ getNorm_(vec3_ v1, vec3_ v2, vec3_ v3);
    static double pointDistFromTri_(vec3_ point, vec3_ v1, vec3_ v2, vec3_ v3);
    //Returns a list of the triangles that were removed and replaced
    std::vector<Triangle_>  processTriangle_(std::list<Triangle_>& tris, Triangle_ tri, std::list<int>& pointsInFront, vec3_& barycenter) const;
    const int CONTRAST_1_ = -1;
    const int CONTRAST_2_ = -2;
    const int CONTRAST_3_ = -3;
    //Data that the different styles need to have available when building geometry
    mutable Teuchos::ArrayRCP<const double> xCoords_; //fine local coordinates
    mutable Teuchos::ArrayRCP<const double> yCoords_;
    mutable Teuchos::ArrayRCP<const double> zCoords_;
    mutable Teuchos::ArrayRCP<const double> cx_; //coarse local coordinates
    mutable Teuchos::ArrayRCP<const double> cy_;
    mutable Teuchos::ArrayRCP<const double> cz_;
    mutable Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId_;
    mutable Teuchos::ArrayRCP<LocalOrdinal> aggSizes_;
    mutable std::vector<bool> isRoot_;
    mutable bool doFineGraphEdges_;
    mutable bool doCoarseGraphEdges_;
    mutable int numNodes_;
    mutable int numAggs_;
    mutable int dims_;
    mutable int myRank_;
    mutable Teuchos::RCP<const Map> nodeMap_; //map used in A and Coordinates to map local ordinals to global ordinals. Need the whole map especially if it's not contiguous.
    mutable Teuchos::RCP<const Map> nodeMapCoarse_; //Map for Ac
    mutable int aggsOffset_;            //in a global list of aggregates, the offset of local aggregate indices
  }; // class AggregationExportFactory
} // namespace MueLu

#define MUELU_AGGREGATIONEXPORTFACTORY_SHORT

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_ */
