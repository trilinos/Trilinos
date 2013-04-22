// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_BlockedTpetraLinearObjContainer_hpp__
#define __Panzer_BlockedTpetraLinearObjContainer_hpp__

#include "Panzer_config.hpp"

// Tpetra includes
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Teuchos_RCP.hpp"

#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include <boost/unordered_map.hpp>

namespace panzer {

/** Linear object container for Block operators, this
  * always assumes the matrix is square.
  */
template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=Kokkos::DefaultNode::DefaultNodeType>
class BlockedTpetraLinearObjContainer : public LinearObjContainer
                                      , public ThyraObjContainer<ScalarT> {
public:
   typedef Thyra::VectorBase<ScalarT> VectorType;
   typedef Thyra::LinearOpBase<ScalarT> CrsMatrixType;

   typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;

   //! Make sure row and column spaces match up 
   bool checkCompatibility() const;

   virtual void clear();

   //! Put a particular scalar in the matrix
   void initializeMatrix(ScalarT value);

   void setMapsForBlocks(const std::vector<Teuchos::RCP<const MapType> > & blockMaps)
   { blockMaps_ = blockMaps; }

   Teuchos::RCP<const MapType> getMapForBlock(std::size_t i) const
   { return blockMaps_[i]; }

   inline void set_x(const Teuchos::RCP<VectorType> & in) { set_x_th(in); }
   inline Teuchos::RCP<VectorType> get_x() const { return get_x_th(); }

   inline void set_dxdt(const Teuchos::RCP<VectorType> & in) { set_dxdt_th(in); }
   inline Teuchos::RCP<VectorType> get_dxdt() const { return get_dxdt_th(); }

   inline void set_f(const Teuchos::RCP<VectorType> & in) { set_f_th(in); }
   inline Teuchos::RCP<VectorType> get_f() const { return get_f_th(); }

   inline void set_A(const Teuchos::RCP<CrsMatrixType> & in) { set_A_th(in); }
   inline Teuchos::RCP<CrsMatrixType> get_A() const { return get_A_th(); }

   // Inherited from LinearObjContainer
   virtual void initialize();

   // Inherited from ThyraObjContainer
 
   void set_x_th(const Teuchos::RCP<VectorType> & in) { x = in; }
   Teuchos::RCP<VectorType> get_x_th() const { return x; }

   void set_dxdt_th(const Teuchos::RCP<VectorType> & in) { dxdt = in; }
   Teuchos::RCP<VectorType> get_dxdt_th() const { return dxdt; }

   void set_f_th(const Teuchos::RCP<VectorType> & in) { f = in; }
   Teuchos::RCP<VectorType> get_f_th() const { return f; }

   void set_A_th(const Teuchos::RCP<CrsMatrixType> & in) { A = in; }
   Teuchos::RCP<CrsMatrixType> get_A_th() const { return A; }

private:
   Teuchos::RCP<VectorType> x, dxdt, f;
   Teuchos::RCP<CrsMatrixType> A;

   std::vector<Teuchos::RCP<const MapType> > blockMaps_;
};

}

#include "Panzer_BlockedTpetraLinearObjContainer_impl.hpp"

#endif
