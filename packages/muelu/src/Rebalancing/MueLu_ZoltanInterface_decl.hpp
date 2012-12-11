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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_ZOLTANINTERFACE_DECL_HPP
#define MUELU_ZOLTANINTERFACE_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#include <zoltan_cpp.h>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ZoltanInterface_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

  /*!
    @class ZoltanInterface
    @brief Interface to Zoltan library.

    This interface provides access to partitioning methods in Zoltan.
    Currently, it supports the RCB algorithm only.
  */

  //FIXME: this class should not be templated
  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ZoltanInterface : public SingleLevelFactoryBase {

    typedef double Scalar; // FIXME
#undef MUELU_ZOLTANINTERFACE_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors
    //@{

    //! Constructor
    ZoltanInterface() {}

    //! Destructor
    virtual ~ZoltanInterface() { }
    //@}

    //! @name Input
    //@{
    void DeclareInput(Level & level) const;
    //@}

    //! @name Build methods.
    //@{
    void Build(Level &level) const;

    //@}

    //! @name Query methods (really functions) required by Zoltan.
    //@{

    /*! Callback function that returns the local number of objects. Required by Zoltan.

    In this case, the number of objects is the number of local rows.

    @param data (in) void pointer to an Xpetra::Matrix.
    @param ierr (out) error code.
    */
    static int GetLocalNumberOfRows(void *data, int *ierr);

    /*! Callback function that returns the local number of nonzeros in the matrix. Required by Zoltan.

    FIXME: Note that this will not work properly for non-point matrices.

    @param data (in) void pointer to an Xpetra::Matrix
    @param weights (out) array whose <tt>i</tt><sup>th</sup> entry is the number of nonzeros in local row \c i.
    @param ierr (out) error code
    */
    static void GetLocalNumberOfNonzeros(void *data, int NumGidEntries, int NumLidEntries, ZOLTAN_ID_PTR gids,
                                         ZOLTAN_ID_PTR lids, int wgtDim, float *weights, int *ierr);

    /*! Callback function that returns the problem dimension. Required by Zoltan.

    @param data (in) void pointer to integer dimension
    @param ierr (out) error code
    */
    static int GetProblemDimension(void *data, int *ierr);


    /*! Callback function that returns the problem dimension. Required by Zoltan.

    @param data (in) void pointer to Xpetra::MultiVector.
    @param coordinates (out) array of double coordinates, arranged like so: [x1 y1 z1 x2 y2 z2 ...].
    @param ierr (out) error code

    TODO -- should I return a view of the coordinates instead of copying them?
    */
    static void GetProblemGeometry(void *data, int numGIDEntries, int numLIDEntries, int numObjectIDs,
                                   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int dim, double *coordinates, int *ierr);

    //@}

  private:

    static ArrayRCP<double> coalesceCoordinates(ArrayRCP<double> coord, LocalOrdinal blksize);

  };  //class ZoltanInterface

} //namespace MueLu

#define MUELU_ZOLTANINTERFACE_SHORT
#endif //if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#endif // MUELU_ZOLTANINTERFACE_DECL_HPP
