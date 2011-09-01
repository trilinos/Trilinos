// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


#ifndef AMESOS2_MULTIVECADAPTER_DEF_HPP
#define AMESOS2_MULTIVECADAPTER_DEF_HPP

#include "Amesos2_TpetraMultiVecAdapter_def.hpp"
// EpetraMultiVecAdapter_def.hpp not included because the specialization is not a template

#include "Amesos2_Util.hpp"	// for getDistributionMap

namespace Amesos2{

  namespace Util {

    ////////////////////////////
    // Copy-getting utilities //
    ////////////////////////////
    
    /*
     * If the multivector scalar type and the desired scalar tpye are
     * the same, then we can do a simple straight copy.
     */
    template <typename MV>
    void same_type_get_copy<MV>::apply(const Teuchos::Ptr<const MV>& mv,
				       const Teuchos::ArrayView<typename MV::scalar_t>& v,
				       const size_t ldx,
				       Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map )
    {
      mv->get1dCopy(v, ldx, distribution_map);
    }

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    void diff_type_get_copy<MV,S>::apply(const Teuchos::Ptr<const MV>& mv,
					 const Teuchos::ArrayView<S>& v,
					 const size_t& ldx,
					 Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map )
    {
      typedef typename MV::scalar_t mv_scalar_t;
	  
      int vals_length = v.size();
      Teuchos::Array<mv_scalar_t> vals_tmp(vals_length);
	  
      mv->get1dCopy(vals_tmp(), ldx, distribution_map);

      for ( int i = 0; i < vals_length; ++i ){
	v[i] = Teuchos::as<S>(vals_tmp[i]);
      }
    }

    /** \internal
     * 
     * \brief Helper class for getting 1-D copies of multivectors
     *
     * Handles datatype conversion when appropriate.
     */
    template <class MV, typename S>
    void get_1d_copy_helper<MV,S>::do_get(const Teuchos::Ptr<const MV>& mv,
					  const Teuchos::ArrayView<S>& vals,
					  const size_t ldx,
					  Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map )
    {
      // Dispatch to the copy function appropriate for the type
      if_then_else<is_same<typename MV::scalar_t,S>::value,
	same_type_get_copy<MV>,
	diff_type_get_copy<MV,S> >::type::apply(mv, vals, ldx, distribution_map);
    }
    
    template <class MV, typename S>
    void get_1d_copy_helper<MV,S>::do_get(const Teuchos::Ptr<const MV>& mv,
					  const Teuchos::ArrayView<S>& vals,
					  const size_t ldx,
					  EDistribution distribution)
    {
      typedef typename MV::local_ordinal_t lo_t;
      typedef typename MV::global_ordinal_t go_t;
      typedef typename MV::global_size_t gs_t;
      typedef typename MV::node_t node_t;
  
      const Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
	= Amesos2::Util::getDistributionMap<lo_t,go_t,gs_t,node_t>(distribution,
								   mv->getGlobalLength(),
								   mv->getComm());
      do_get(mv, vals, ldx, Teuchos::ptrInArg(*map));
    }

    template <class MV, typename S>
    void get_1d_copy_helper<MV,S>::do_get(const Teuchos::Ptr<const MV>& mv,
					  const Teuchos::ArrayView<S>& vals,
					  const size_t ldx)
    {
      const Teuchos::RCP<const Tpetra::Map<typename MV::local_ordinal_t,
	typename MV::global_ordinal_t,
	typename MV::node_t> > map
	= mv->getMap();
      do_get(mv, vals, ldx, Teuchos::ptrInArg(*map));
    }
    
    
    ///////////////////////////
    // Copy-puting utilities //
    ///////////////////////////
    
    template <typename MV>
    void same_type_data_put<MV>::apply(const Teuchos::Ptr<MV>& mv,
				       const Teuchos::ArrayView<typename MV::scalar_t>& data,
				       const size_t ldx,
				       Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map )
    {
      mv->put1dData(data, ldx, distribution_map);
    }
    
    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    void diff_type_data_put<MV,S>::apply(const Teuchos::Ptr<MV>& mv,
					 const Teuchos::ArrayView<S>& data,
					 const size_t& ldx,
					 Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map )
    {
      typedef typename MV::scalar_t mv_scalar_t;
      
      int vals_length = data.size();
      Teuchos::Array<mv_scalar_t> data_tmp(vals_length);
      
      for ( int i = 0; i < vals_length; ++i ){
	data_tmp[i] = Teuchos::as<mv_scalar_t>(data[i]);
      }
      
      mv->put1dData(data_tmp(), ldx, distribution_map);
    }
    
    
    /** \internal
     * 
     * \brief Helper class for putting 1-D data arrays into multivectors
     *
     * Handles dataype conversion when necessary before putting the data.
     */
    template <class MV, typename S>
    void put_1d_data_helper<MV,S>::do_put(const Teuchos::Ptr<MV>& mv,
					  const Teuchos::ArrayView<S>& data,
					  const size_t ldx,
					  Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map )
    {
      // Dispatch to the copy function appropriate for the type
      if_then_else<is_same<typename MV::scalar_t,S>::value,
	same_type_data_put<MV>,
	diff_type_data_put<MV,S> >::type::apply(mv, data, ldx, distribution_map);
    }
    
    template <class MV, typename S>
    void put_1d_data_helper<MV,S>::do_put(const Teuchos::Ptr<MV>& mv,
					  const Teuchos::ArrayView<S>& data,
					  const size_t ldx,
					  EDistribution distribution)
    {
      typedef typename MV::local_ordinal_t lo_t;
      typedef typename MV::global_ordinal_t go_t;
      typedef typename MV::global_size_t gs_t;
      typedef typename MV::node_t node_t;
      
      const Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
	= Amesos2::Util::getDistributionMap<lo_t,go_t,gs_t,node_t>(distribution,
								   mv->getGlobalLength(),
								   mv->getComm());
      do_put(mv, data, ldx, Teuchos::ptrInArg(*map));
    }
    
    template <class MV, typename S>
    void put_1d_data_helper<MV,S>::do_put(const Teuchos::Ptr<MV>& mv,
					  const Teuchos::ArrayView<S>& data,
					  const size_t ldx)
    {
      const Teuchos::RCP<const Tpetra::Map<typename MV::local_ordinal_t,
	typename MV::global_ordinal_t,
	typename MV::node_t> > map
	= mv->getMap();
      do_put(mv, data, ldx, Teuchos::ptrInArg(*map));
    }
    
  } // end namespace Util
  
} // end namespace Amesos2

#endif	// AMESOS2_EPETRAMULTIVECADAPTER_DEF
