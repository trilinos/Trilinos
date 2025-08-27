// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/PEncr_RecoverUtilities.hpp>
//#include <stk_util/util/Fortran.hpp>
#include <Teuchos_config.h>

  namespace percept {

    // This number is used to determine the effective rank of the least squares
    // system. The effective rank is the order of the largest leading triangular
    // submatrix in the QR factorization with pivoting, whose estimated condition
    // number is less than 1/RCOND.
    static double least_squares_rcond = 5.0e-12;

    // LAPACK routine,
    // Computes the minimum-norm solution to a linear least
    // squares problem using a complete orthogonal factorization.
    extern "C" void F77_FUNC(dgelsy,DGELSY)(
                                          const P_int*  M,
                                          const P_int*  N,
                                          const P_int*  NRHS,
                                          P_Real* A,
                                          const P_int*  LDA,
                                          P_Real* B,
                                          const P_int*  LDB,
                                          P_int*  JPVT,
                                          const P_Real* RCOND,
                                          P_int*  RANK,
                                          P_Real* WORK,
                                          const P_int*  LWORK,
                                          P_int*  INFO );


    void
    RecoverUtilities::getValuesAtNodes(const std::vector<stk::mesh::Entity> &nodes,
                                       const stk::mesh::FieldBase& field,
                                       std::vector<double> &values)
    {
      const P_uint nNodes = nodes.size();
      //const P_uint length = field.field()->length();
      const P_uint length = field.max_size();

      STK_ThrowAssert(values.size()==nNodes*length);

      for (P_uint n=0; n<nNodes; n++)
        {
          const P_Real * data = m_eMesh.field_data(&field, nodes[n]);
          //double * PerceptMesh::field_data(const stk::mesh::FieldBase *field, const stk::mesh::Entity entity, unsigned *stride)
          STK_ThrowAssert(data);

          for (P_uint i=0; i<length; i++)
            values[n*length+i] = data[i];
        }
    }

    void RecoverUtilities::
    object_patch_if(const stk::mesh::BulkData& bulk_data,
                    stk::mesh::Entity obj,
                    const Filter& predicate,
                    std::vector<stk::mesh::Entity>& element_patch,
                    std::vector<stk::mesh::Entity>& nodes_in_patch,
                    const stk::mesh::EntityRank /*patch_obj_type*/)
    {
      element_patch.resize(0);
      nodes_in_patch.resize(0);
      
      std::set<stk::mesh::Entity> elem_patch_set;
      std::set<stk::mesh::Entity> node_patch_set;
      if ( bulk_data.entity_rank(obj) == stk::topology::NODE_RANK)
        {
          const MyPairIterRelation node_elems(*m_eMesh.get_bulk_data(), obj, stk::topology::ELEMENT_RANK) ;
          for (unsigned i=0; i < node_elems.size(); ++i)
            {
              stk::mesh::Entity elem = node_elems[i].entity();
              if (predicate.pass(elem))
                elem_patch_set.insert(elem);

              const MyPairIterRelation elem_nodes(*m_eMesh.get_bulk_data(), elem, stk::topology::NODE_RANK) ;
              for (unsigned j=0; j < elem_nodes.size(); ++j)
                {
                  if (predicate.pass(elem_nodes[j].entity()))
                    node_patch_set.insert(elem_nodes[j].entity());
                }
            }

          element_patch.assign(elem_patch_set.begin(), elem_patch_set.end());
          nodes_in_patch.assign(node_patch_set.begin(), node_patch_set.end());
          return;
        }
      else if ( bulk_data.entity_rank(obj) == stk::topology::ELEMENT_RANK )
        {
          const MyPairIterRelation elem_nodes_0(*m_eMesh.get_bulk_data(), obj, stk::topology::NODE_RANK) ;
          for (unsigned i0=0; i0 < elem_nodes_0.size(); ++i0)
            {
              stk::mesh::Entity node = elem_nodes_0[i0].entity();
              const MyPairIterRelation node_elems(*m_eMesh.get_bulk_data(), node, stk::topology::ELEMENT_RANK) ;
              for (unsigned i=0; i < node_elems.size(); ++i)
                {
                  stk::mesh::Entity elem = node_elems[i].entity();
                  if (predicate.pass(elem))
                    elem_patch_set.insert(elem);

                  const MyPairIterRelation elem_nodes(*m_eMesh.get_bulk_data(), elem, stk::topology::NODE_RANK) ;
                  for (unsigned j=0; j < elem_nodes.size(); ++j)
                    {
                      if (predicate.pass(elem_nodes[j].entity()))
                        node_patch_set.insert(elem_nodes[j].entity());
                    }
                }
            }

          element_patch.assign(elem_patch_set.begin(), elem_patch_set.end());
          nodes_in_patch.assign(node_patch_set.begin(), node_patch_set.end());
          return;
        }
      else
        {
          //error
        }
    }

    void RecoverUtilities::createMeshObjPatch(stk::mesh::Entity obj,
                                              std::vector<stk::mesh::Entity> & elmts,
                                              std::vector<stk::mesh::Entity> & nodes,
                                              //const bool includeConstrained,
                                              const stk::mesh::EntityRank patch_type,
                                              stk::mesh::Selector *meshpart )
    {
      const stk::mesh::BulkData & stk_mesh = *m_eMesh.get_bulk_data();

      elmts.clear();
      nodes.clear();
      //object_patch_if (stk_mesh, obj, ShellFilter(), elmts, nodes, includeConstrained, patch_type);
      object_patch_if (stk_mesh, obj, Filter(), elmts, nodes, patch_type);

      if (patch_type!=stk::topology::ELEMENT_RANK)
        {
          filterMeshObjsFromMeshPart(nodes, meshpart);
          filterMeshObjsFromMeshPart(elmts, meshpart);
        }

#if 0
      if (encorelog.shouldPrint(LOG_RECOVER))
        {
          encorelog << "    num nodes = " << nodes.size() << Diag::dendl;
          encorelog << "    num elmts = " << elmts.size() << Diag::dendl;

          encorelog << "    created patch nodes: ";
          for (P_uint i=0; i<nodes.size(); ++i)
            encorelog << stk_mesh.global_id(nodes[i]) << " ";
          encorelog << Diag::dendl;

          encorelog << "    created patch objs: ";
          for (P_uint i=0; i<elmts.size(); ++i)
            encorelog << stk_mesh.global_id(elmts[i]) << " ";
          encorelog << Diag::dendl;

          encorelog << "    patch boundary nodes: ";

          stk::mesh::Selector exposed_block_selector(
                                                     region->mesh_meta_data().exposed_boundary_part());
          for (P_uint i=0; i<nodes.size(); ++i)
            if (exposed_block_selector(stk_mesh.bucket(nodes[i])))
              encorelog << stk_mesh.global_id(nodes[i]) << " ";
          encorelog << Diag::dendl;

          const sierra::Mmod::HAdapt * hAdapt =
            Mmod::HAdapt::get_hadapt(*(const_cast<Fmwk::Region*>(region)));
          if (hAdapt)
            {
              stk::mesh::Selector constrained_selector(hAdapt->constrained_mesh_part());

              encorelog << "    patch constrained nodes: ";
              for (P_uint i=0; i<nodes.size(); ++i)
                {
                  if (constrained_selector(stk_mesh.bucket(nodes[i])))
                    encorelog << stk_mesh.global_id(nodes[i]) << " ";
                }
              encorelog << Diag::dendl;
            }
        }
#endif
    }

    void RecoverUtilities::filterMeshObjsFromMeshPart(std::vector<stk::mesh::Entity> & objs,
                                                      stk::mesh::Selector *meshpart)
    {
      if (NULL == meshpart) return;
      const stk::mesh::BulkData& stk_mesh = *m_eMesh.get_bulk_data();

      std::vector<stk::mesh::Entity>::iterator iobj = objs.begin();
      while (iobj!=objs.end())
        {
          if ( !(*meshpart)(stk_mesh.bucket(*iobj)) )
            iobj = objs.erase(iobj);
          else
            iobj++;
        }
    }

#if 0
    bool ShellFilter::pass (stk::mesh::Entity obj) const
    {
      /* %TRACE% */ Traceback trace__("sierra::Encr::RecoverUtilities::ShellFilter::pass(stk::mesh::Entity obj) const"); /* %TRACE% */
      return true; //isShellElement(obj.getCellTopology());
    }
#endif

    void RecoverUtilities::fitPolynomialLS(const P_uint spatial_dim,
                                           const P_uint sample_dim,
                                           const P_uint nSamplePts,
                                           const P_uint polyDegree,
                                           const std::vector<double> & sample_coords,
                                           std::vector<double> & function_sample,
                                           const P_uint nEvalPts,
                                           const std::vector<double> & eval_coords,
                                           std::vector<double> & eval_values)
    {
      const P_uint nBasis = simplexPolynomialBasisSize(spatial_dim, polyDegree);

      // evaluate the poly basis at the sampling points.
      std::vector<double> basis_sample (nBasis*nSamplePts, 0.);

      evalSimplexPolynomialBasis(spatial_dim,
                                 nSamplePts,
                                 polyDegree,
                                 sample_coords,
                                 basis_sample);

      // compute the LS fit for the poly coeffs.
      // overwrites function_sample with poly coefficients and basis_sample with LS factorization
      leastSquaresSolve(sample_dim,
                        nBasis,
                        nSamplePts,
                        function_sample,
                        basis_sample);

      // re-evaluate the poly basis at the evaluation points.
      basis_sample.clear();
      basis_sample.resize(nBasis*nEvalPts, 0.);

      evalSimplexPolynomialBasis(spatial_dim,
                                 nEvalPts,
                                 polyDegree,
                                 eval_coords,
                                 basis_sample);

      computePolynomialValues(sample_dim,
                              nEvalPts,
                              nBasis,
                              nSamplePts,
                              eval_values,
                              basis_sample,
                              function_sample);
    }

    void RecoverUtilities::fitSimpleAverage(const P_uint sample_dim,
                                            const P_uint nSamplePts,
                                            std::vector<double> & function_sample,
                                            const P_uint nEvalPts,
                                            const std::vector<double> & /*eval_coords*/,
                                            std::vector<double> & eval_values)
    {
      const double wt = 1.0 / (double) nSamplePts;

      for (P_uint i=0; i<sample_dim; i++)
        {
          eval_values[i] = 0.;

          for (P_uint j=0; j<nSamplePts; j++)
            eval_values[i] += function_sample[i*nSamplePts + j];

          eval_values[i] *= wt;
        }

      // copy first eval values to other points
      if (nEvalPts>1)
        {
          for (P_uint j=1; j<nEvalPts; j++)
            for (P_uint i=0; i<sample_dim; i++)
              eval_values[j*sample_dim + i] = eval_values[i];
        }
    }


    P_uint RecoverUtilities::simplexPolynomialBasisSize(const P_uint spatial_dim,
                                                        const P_uint polyDegree)
    {
      STK_ThrowAssert(polyDegree>0 && (spatial_dim==2 || spatial_dim==3));

      return (spatial_dim==2) ? ( (polyDegree+1)*(polyDegree+2) )/2
        : ( (polyDegree+1)*(polyDegree+2)*(polyDegree+3) )/6;
    }


    void RecoverUtilities::evalSimplexPolynomialBasis(const P_uint spatial_dim,
                                                      const P_uint nSamplePts,
                                                      const P_uint polyDegree,
                                                      const std::vector<double> & sample_coords,
                                                      std::vector<double> & sample_values)
    {
      for (P_uint q=0; q<nSamplePts; ++q)
        {
          const double x = sample_coords[q*spatial_dim];
          const double y = sample_coords[q*spatial_dim + 1];

          sample_values[0*nSamplePts + q] = 1.0;
          sample_values[1*nSamplePts + q] = x;
          sample_values[2*nSamplePts + q] = y;

          if (spatial_dim==2)
            {
              if (polyDegree>=2)
                {
                  sample_values[3*nSamplePts + q] = x*x;
                  sample_values[4*nSamplePts + q] = x*y;
                  sample_values[5*nSamplePts + q] = y*y;
                }

              if (polyDegree>=3)
                {
                  sample_values[6*nSamplePts + q] = x*x*x;
                  sample_values[7*nSamplePts + q] = x*x*y;
                  sample_values[8*nSamplePts + q] = x*y*y;
                  sample_values[9*nSamplePts + q] = y*y*y;
                }
            }
          else if (spatial_dim==3)
            {
              const double z = sample_coords[q*spatial_dim + 2];

              sample_values[3*nSamplePts + q] = z;

              if (polyDegree>=2)
                {
                  sample_values[4*nSamplePts + q] = x*x;
                  sample_values[5*nSamplePts + q] = x*y;
                  sample_values[6*nSamplePts + q] = x*z;
                  sample_values[7*nSamplePts + q] = y*y;
                  sample_values[8*nSamplePts + q] = y*z;
                  sample_values[9*nSamplePts + q] = z*z;
                }

              if (polyDegree>=3)
                {
                  sample_values[10*nSamplePts + q] = x*x*x;
                  sample_values[11*nSamplePts + q] = x*x*y;
                  sample_values[12*nSamplePts + q] = x*x*z;
                  sample_values[13*nSamplePts + q] = x*y*y;
                  sample_values[14*nSamplePts + q] = x*y*z;
                  sample_values[15*nSamplePts + q] = x*z*z;
                  sample_values[16*nSamplePts + q] = y*y*y;
                  sample_values[17*nSamplePts + q] = y*z*y;
                  sample_values[18*nSamplePts + q] = y*z*z;
                  sample_values[19*nSamplePts + q] = z*z*z;
                }
            }
        }
    }

    void RecoverUtilities::
    leastSquaresSolve(const P_uint sample_dim,
                      const P_uint nBasis,
                      const P_uint nSamplePts,
                      std::vector<double> & func_samp,
                      std::vector<double> & basis_samp)
    {
      static std::vector<P_int>  my_iwork;
      static std::vector<double> my_work(1);

      STK_ThrowRequire(sample_dim>=1);
      STK_ThrowRequire(nSamplePts>=1);

      const P_int NRHS = sample_dim;
      const P_int M = nSamplePts;
      P_int N = nBasis;
      P_int LDB;

      // determine leading dimension of the matrix used to hold
      // both the RHS and the solution of the least squares system
      if (M < N) {

        // There are less samples than there are coefficients to solve; i.e.,
        // the system is under-determined.

        LDB = N;

        // We need to enlarge the leading dimension of the RHS matrix
        // so it has room to hold the solution, we do this by expanding the
        // stride.
        func_samp.resize(NRHS*N); // larger than (NRHS*M)

        // now copy data into expanded matrix, column by column
        for (P_int d=sample_dim-1; d>-1; --d)
          {
            for (P_int q=M-1; q>-1; --q)
              func_samp[d*N+q] = func_samp[d*M+q];
          }
        // zero out new rows
        for (P_int q=M; q<N; ++q)
          {
            for (P_int d=0; d<NRHS; ++d)
              func_samp[d*N+q] = 0.;
          }

        // Now reduce the size of the basis we are using.
        N = M;
      }
      else {

        // There are at least as many samples as the basis,
        // or the system is over-determined.

        LDB = M;
      }

      //
      // Query the optimal size of the workspace array.
      //
      P_int lwork = -1;

      // zero out the pivot array(1,...,N)
      my_iwork.resize( N );
      std::fill( my_iwork.begin(), my_iwork.end(), 0 );

      P_int rank = 0;
      P_int ierr = -1;

      F77_FUNC(dgelsy,DGELSY)(
                            &M, &N, &NRHS,
                            NULL, &M,
                            NULL, &LDB,
                            &my_iwork[0],
                            &least_squares_rcond,
                            &rank,
                            &my_work[0],
                            &lwork,
                            &ierr );

      lwork = (int)my_work[0];

      if( ierr != 0)
        throw std::runtime_error("ERROR in least_squares_solve when sizing the workspace array "+toString(ierr));

      my_work.resize( lwork );

      // Compute the solution to the linear least squares problem.
      F77_FUNC(dgelsy,DGELSY) (
                             &M, &N, &NRHS,
                             &basis_samp[0], &M,
                             &func_samp[0], &LDB,
                             &my_iwork[0],
                             &least_squares_rcond,
                             &rank,
                             &my_work[0],
                             &lwork,
                             &ierr );

      if( ierr != 0 )
        throw std::runtime_error("ERROR in least_squares_solve when solving for coefficients: ierr = "+toString(ierr));
    }

    void
    RecoverUtilities::computePolynomialValues(const P_uint sample_dim,
                            const P_uint nEvalPts,
                            const P_uint nBasis,
                            const P_uint nSamplePts,
                            std::vector<double> & eval_values,
                            const std::vector<double> & basis_sample,
                            const std::vector<double> & coeff)
    {
      // need to handle under and over determined cases
      const P_uint M = (nSamplePts < nBasis) ? nBasis : nSamplePts;

      // TODO: use the BLAS Luke ...
      for (P_uint d=0; d<sample_dim; ++d)
        {
          for (P_uint n=0; n<nEvalPts; ++n)
            {
              eval_values[n*sample_dim + d] = 0.;

              for (P_uint j=0; j<nBasis; ++j)
                // NOTE: we use nSamplePts here because the coeff array was written
                // over the function_sample array.
                eval_values[n*sample_dim + d] += basis_sample[j*nEvalPts + n] * coeff[d*M + j];
            }
        }
    }

#if 0
    //------------------------------------------------------------------------
    // This method will return false if the given node is not the vertex of
    // every element that it is USED_BY.
    bool nodeIsVertex(stk::mesh::Entity node,
                      const Fmwk::FieldRef & vertexNodeMarkRef)
    {
      /* %TRACE% */ Traceback trace__("sierra::Encr::RecoverUtilities::nodeIsVertex(stk::mesh::Entity node)"); /* %TRACE% */
      const stk::mesh::BulkData& stk_mesh = vertexNodeMarkRef.field()->registrar().mesh_meta_data().mesh_bulk_data();
      const P_int * node_mark = field_data(vertexNodeMarkRef, node);

      if (encorelog.shouldPrint(LOG_RECOVER))
        {
          const String isIsNot = (*node_mark>0) ? " is " : " is NOT ";
          encorelog << "  node " << stk_mesh.global_id(node)
                    << isIsNot << "a vertex node"
                    << Diag::dendl;
        }

      return (*node_mark>0);
    }
#endif

    std::vector<double> RecoverUtilities::
    computeElementRecoveryValue(
                                stk::mesh::Entity el,
                                const std::vector<double>& point,
                                const stk::mesh::FieldBase& coeff,
                                const stk::mesh::FieldBase& coords,
                                const P_uint coeff_dim)
    {
      std::vector<stk::mesh::Entity> elmts;
      std::vector<stk::mesh::Entity> nodes;

      // TODO: enable adaptivity
      //const bool includeConstrained = false;

      const stk::mesh::BulkData & stk_mesh = *m_eMesh.get_bulk_data();
      createMeshObjPatch( el, elmts, nodes, stk_mesh.entity_rank(el));

      // TODO: enable recovery of quadratic Fields
      const P_uint polyDegree = 2;

      const P_uint spatial_dim = point.size();
      const P_uint nSamplePts = nodes.size();

      std::vector<double> sample_coords (nSamplePts*spatial_dim);
      std::vector<double> function_sample (nSamplePts*coeff_dim);

      getValuesAtNodes(nodes, coords, sample_coords);
      getValuesAtNodes(nodes, coeff, function_sample);

      // We are taking in only a single point for now
      const P_uint nEvalPts = 1;

      std::vector<double> eval_values(coeff_dim,0.0);

      fitPolynomialLS(
                      spatial_dim,
                      coeff_dim,
                      nSamplePts,
                      polyDegree,
                      sample_coords,
                      function_sample,
                      nEvalPts,
                      point,
                      eval_values);

      return eval_values;
    }


  }

