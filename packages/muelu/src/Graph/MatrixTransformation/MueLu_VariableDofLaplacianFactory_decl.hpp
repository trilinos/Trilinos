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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DECL_HPP_
#define PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_VariableDofLaplacianFactory_fwd.hpp"
#include "MueLu_Level_fwd.hpp"

#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class VariableDofLaplacianFactory class.
  @brief Factory for building scalar Laplace operator (that is used as fake operator for variable dof size problems)

  Build distance Laplacian associated with input matrix A (which might have a variable number of DOFs per node).
  Coordinates are needed to calculate the distance laplacian values. The user-provided array "DofPresent" stores whether
  an array is present (=1) or not (=0) in the matrix. The length of the array is number of nodes * maxDofPerNode and
  therefore it might be larger or equal than the number of rows in the input matrix.

  The factory produces the distance laplacian matrix A as output (with one dof per node) as well as the coarse version
  of the DofStatus (needed for the next coarser level), containing information about (artificial) Dirichlet rows in the matrix.

  @ingroup MueLuGraphClasses

  ## Input/output of VariableDofLaplacianFactory ##

  ### User parameters of VariableDofLaplacianFactory ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
   A        | Factory | null |   | * | *  | Generating factory of the input matrix A with potentially variable number of DOFs. Might be padded or non-padded. Padded means, that the matrix has additional artificial rows and columns to have a constant number of DOFs per node.
   Coordinates | Factory | null | | * | * | Generating factory for Coordinates needed for building distance laplacian.
   DofPresent | Teuchos::ArrayRCP<LocalOrdinal> | NoFactory | | | (*) | Optional array containing information whether DOF is actually present in matrix or not.
   Advanced Dirichlet: threshold | double | 1e-5 |   | * |   | Drop tolerance for Dirichlet detection
   Variable DOF amalgamation: threshold | double | 1.8e-9 |  | * |  | Drop tolerance for amalgamation process
   maxDofPerNode | int | 1 |  | * |  | Maximum number of DOFs per node


  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see VariableDofLaplacianFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see VariableDofLaplacianFactory::DeclareInput).

  ### Variables provided by VariableDofLaplacianFactory ###

  After TentativePFactory::Build the following data is available (if requested)

  Parameter | generated by | description
  ----------|--------------|------------
  | A       | VariableDofLaplacianFactory   | Laplacian operator
  | DofStatus | VariableDofLaplacianFactory | Status array for next coarse level
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class VariableDofLaplacianFactory : public SingleLevelFactoryBase {
#undef MUELU_VARIABLEDOFLAPLACIANFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  VariableDofLaplacianFactory();

  //! Destructor
  virtual ~VariableDofLaplacianFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  void Build(Level& currentLevel) const;  // Build

 private:
  void buildPaddedMap(const Teuchos::ArrayRCP<const LocalOrdinal>& dofPresent, std::vector<LocalOrdinal>& map, size_t nDofs) const;
  void assignGhostLocalNodeIds(const Teuchos::RCP<const Map>& rowDofMap, const Teuchos::RCP<const Map>& colDofMap, std::vector<LocalOrdinal>& myLocalNodeIds, const std::vector<LocalOrdinal>& dofMap, size_t maxDofPerNode, size_t& nLocalNodes, size_t& nLocalPlusGhostNodes, Teuchos::RCP<const Teuchos::Comm<int> > comm) const;
  void squeezeOutNnzs(Teuchos::ArrayRCP<size_t>& rowPtr, Teuchos::ArrayRCP<LocalOrdinal>& cols, Teuchos::ArrayRCP<Scalar>& vals, const std::vector<bool>& keep) const;
  void buildLaplacian(const Teuchos::ArrayRCP<size_t>& rowPtr, const Teuchos::ArrayRCP<LocalOrdinal>& cols, Teuchos::ArrayRCP<Scalar>& vals, const size_t& numdim, const RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >& ghostedCoords) const;

  template <class listType>
  void MueLu_az_sort(listType list[], size_t N, size_t list2[], Scalar list3[]) const {
    /* local variables */

    listType RR, K;
    size_t l, r, j, i;
    int flag;
    size_t RR2;
    Scalar RR3;

    /*********************** execution begins ******************************/

    if (N <= 1) return;

    l  = N / 2 + 1;
    r  = N - 1;
    l  = l - 1;
    RR = list[l - 1];
    K  = list[l - 1];

    if ((list2 != NULL) && (list3 != NULL)) {
      RR2 = list2[l - 1];
      RR3 = list3[l - 1];
      while (r != 0) {
        j    = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[i - 1]  = list[j - 1];
              list2[i - 1] = list2[j - 1];
              list3[i - 1] = list3[j - 1];
            } else {
              flag = 0;
            }
          }
        }

        list[i - 1]  = RR;
        list2[i - 1] = RR2;
        list3[i - 1] = RR3;

        if (l == 1) {
          RR  = list[r];
          RR2 = list2[r];
          RR3 = list3[r];

          K        = list[r];
          list[r]  = list[0];
          list2[r] = list2[0];
          list3[r] = list3[0];
          r        = r - 1;
        } else {
          l   = l - 1;
          RR  = list[l - 1];
          RR2 = list2[l - 1];
          RR3 = list3[l - 1];
          K   = list[l - 1];
        }
      }

      list[0]  = RR;
      list2[0] = RR2;
      list3[0] = RR3;
    } else if (list2 != NULL) {
      RR2 = list2[l - 1];
      while (r != 0) {
        j    = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[i - 1]  = list[j - 1];
              list2[i - 1] = list2[j - 1];
            } else {
              flag = 0;
            }
          }
        }

        list[i - 1]  = RR;
        list2[i - 1] = RR2;

        if (l == 1) {
          RR  = list[r];
          RR2 = list2[r];

          K        = list[r];
          list[r]  = list[0];
          list2[r] = list2[0];
          r        = r - 1;
        } else {
          l   = l - 1;
          RR  = list[l - 1];
          RR2 = list2[l - 1];
          K   = list[l - 1];
        }
      }

      list[0]  = RR;
      list2[0] = RR2;
    } else if (list3 != NULL) {
      RR3 = list3[l - 1];
      while (r != 0) {
        j    = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[i - 1]  = list[j - 1];
              list3[i - 1] = list3[j - 1];
            } else {
              flag = 0;
            }
          }
        }

        list[i - 1]  = RR;
        list3[i - 1] = RR3;

        if (l == 1) {
          RR  = list[r];
          RR3 = list3[r];

          K        = list[r];
          list[r]  = list[0];
          list3[r] = list3[0];
          r        = r - 1;
        } else {
          l   = l - 1;
          RR  = list[l - 1];
          RR3 = list3[l - 1];
          K   = list[l - 1];
        }
      }

      list[0]  = RR;
      list3[0] = RR3;

    } else {
      while (r != 0) {
        j    = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[i - 1] = list[j - 1];
            } else {
              flag = 0;
            }
          }
        }

        list[i - 1] = RR;

        if (l == 1) {
          RR = list[r];

          K       = list[r];
          list[r] = list[0];
          r       = r - 1;
        } else {
          l  = l - 1;
          RR = list[l - 1];
          K  = list[l - 1];
        }
      }

      list[0] = RR;
    }
  }

};  // class CoalesceDropFactory

}  // namespace MueLu

#define MUELU_VARIABLEDOFLAPLACIANFACTORY_SHORT

#endif /* PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DECL_HPP_ */
