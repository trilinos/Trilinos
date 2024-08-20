// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_INTEGRATION_VALUES2_HPP
#define PANZER_INTEGRATION_VALUES2_HPP

#include "Teuchos_RCP.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_Dimension.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid2_Cubature.hpp"

namespace panzer {

  class SubcellConnectivity;

  template <typename Scalar>
  class IntegrationValues2 {
  public:
    typedef typename ArrayTraits<Scalar,PHX::MDField<Scalar> >::size_type size_type;

    typedef PHX::MDField<Scalar> ArrayDynamic;
    typedef PHX::MDField<double> DblArrayDynamic;

    typedef PHX::MDField<Scalar,IP>               Array_IP;
    typedef PHX::MDField<Scalar,IP,Dim>           Array_IPDim;
    typedef PHX::MDField<Scalar,Point>            Array_Point;
    typedef PHX::MDField<Scalar,Cell,IP>          Array_CellIP;
    typedef PHX::MDField<Scalar,Cell,IP,Dim>      Array_CellIPDim;
    typedef PHX::MDField<Scalar,Cell,IP,Dim,Dim>  Array_CellIPDimDim;
    typedef PHX::MDField<Scalar,Cell,BASIS,Dim>   Array_CellBASISDim;

    typedef PHX::MDField<const Scalar,IP>               ConstArray_IP;
    typedef PHX::MDField<const Scalar,IP,Dim>           ConstArray_IPDim;
    typedef PHX::MDField<const Scalar,Point>            ConstArray_Point;
    typedef PHX::MDField<const Scalar,Cell,IP>          ConstArray_CellIP;
    typedef PHX::MDField<const Scalar,Cell,IP,Dim>      ConstArray_CellIPDim;
    typedef PHX::MDField<const Scalar,Cell,IP,Dim,Dim>  ConstArray_CellIPDimDim;
    typedef PHX::MDField<const Scalar,Cell,BASIS,Dim>   ConstArray_CellBASISDim;

    /**
     * \brief Base constructor
     *
     * \param[in] pre Prefix to apply to all internal field names
     * \param[in] allocArrays (Classic Interface Only) Allocate array data in 'setupArrays' call
     */
    IntegrationValues2(const std::string & pre="",
                       const bool allocArrays=false);


    // =====================================================================================================
    // Classic Interface (DEPRECATED)

    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir);

    /** \brief Evaluate basis values.

        @param cell_node_coordinates [in] Cell node coordinates, not
        basis coordinates.
        @param num_cells [in] (optional) number of cells in the
        workset. This can be less than the workset size. If set to
        zero, extent(0) of the evaluated array is used which equates
        to the workset size.
        @param face_connectivity [in] (optional) connectivity used to
        enforce quadrature alignment for surface integration.
     */
    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & cell_node_coordinates,
                        const int num_cells = -1,
                        const Teuchos::RCP<const SubcellConnectivity> & face_connectivity = Teuchos::null,
                        const int num_virtual_cells = -1);

    /** \brief Match IP.

       Optionally provide IP coordinates for an element 'other' that
       shares the same side. If provided, a permutation of the
       cubature points is calculated so that the integration values
       are ordered according to the other element's. This permutation
       is then applied so that all fields are ordered accordingly in
       their IP dimension.

        @param num_cells [in] (optional) number of cells in the
        workset. This can be less than the workset size. If set to
        zero, extent(0) of the evaluated array is used which equates
        to the workset size.
    */
    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & cell_node_coordinates,
                        const PHX::MDField<Scalar,Cell,IP,Dim> & other_ip_coordinates,
                        const int num_cells = -1);

    // Reference space quantities
    mutable Array_IPDim cub_points;              // <IP,Dim>
    mutable Array_IPDim side_cub_points;         // <IP,Dim> points on face topology (dim-1)
    mutable Array_IP cub_weights;                // <IP>

    // Physical space quantities
    mutable Array_CellBASISDim node_coordinates; // <Cell,BASIS,Dim>
    mutable Array_CellIPDimDim jac;              // <Cell,IP,Dim,Dim>
    mutable Array_CellIPDimDim jac_inv;          // <Cell,IP,Dim,Dim>
    mutable Array_CellIP jac_det;                // <Cell,IP>
    mutable Array_CellIP weighted_measure;       // <Cell,IP>
    mutable Array_CellIPDim weighted_normals;    // <Cell,IP,Dim>
    mutable Array_CellIPDim surface_normals;    // <Cell,IP,Dim>
    mutable Array_CellIPDimDim surface_rotation_matrices;    // <Cell,IP,Dim,Dim>
      // this (appears) is a matrix where the first row is the "normal" direction
      // and the remaining two rows lie in the hyperplane

    // for Shakib stabilization <Cell,IP,Dim,Dim>
    mutable Array_CellIPDimDim covarient;
    mutable Array_CellIPDimDim contravarient;
    mutable Array_CellIP norm_contravarient;

    // integration points
    mutable Array_CellIPDim ip_coordinates;      // <Cell,IP,Dim>
    mutable Array_CellIPDim ref_ip_coordinates;  // <Cell,IP,Dim> for Control Volumes or Surface integrals


    Teuchos::RCP<const panzer::IntegrationRule> int_rule;

    Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>> intrepid_cubature;

    // =====================================================================================================
    // Lazy evaluation interface

    /**
     * The lazy evaluation construction path is designed such that you only allocate and fill arrays on demand,
     * with an option of caching those generated fields. This is useful for when we are worried about a
     * code's memory footprint.
     *
     * The process for setting up one of these objects is to initialize the IntegrationValues2 class as follows:
     *
     * IntegrationValues2<double> values;   // Constructor
     * values.setup(ir, cell_node_coordinates);  // Required - must come before all other calls
     * values.setupPermutations(...);            // Optional - only needed for surface integration and some side integration rules - must come before get*() calls
     * auto array = values.get*();               // Lazy evaluate whatever field you need
     */

    /**
     * \brief Main setup call for the lazy evaluation interface
     *
     * \todo Instead of IntegrationRule, we just need to load the integration descriptor and the cell topology
     *
     * \param[in] ir Integration rule descripting integration scheme
     * \param[in] cell_node_coordinates Node/Vertex <cell, node, dim> coordinates describing cell geometry
     * \param[in] num_cells In case you need to only generate integration values for the first 'num_cells' of the node_coordinates - defaults to all cells
     */
    void
    setup(const Teuchos::RCP<const panzer::IntegrationRule>& ir,
          const PHX::MDField<Scalar,Cell,NODE,Dim> & cell_node_coordinates,
          const int num_cells = -1);

    /**
     * \brief Initialize the permutation arrays given a face connectivity
     *
     * \note REQUIRED FOR SURFACE INTEGRATION
     * \note Must be called AFTER setup
     * \note Virtual cells have a unique way to generate surface normals and rotation matrices, hence we need to know how many are included
     *
     * \param[in] face_connectivity Connectivity describing how sides are connected to cells
     * \param[in] num_virtual_cells Number of virtual cells included in the node coordinates (found at end of node coordinate array)
     */
    void
    setupPermutations(const Teuchos::RCP<const SubcellConnectivity> & face_connectivity,
                      const int num_virtual_cells);

    /**
     * \brief Initialize the permutation arrays given another IntegrationValues2<Scalar>::getCubaturePoints() array
     *
     * \note Required if you want points to line up between two integration values (e.g. for side integration)
     * \note Must be called AFTER setup
     *
     * \param[in] other_ip_coordinates Cubature points to align with
     */
    void
    setupPermutations(const PHX::MDField<Scalar,Cell,IP,Dim> & other_ip_coordinates);

    /**
     * \brief Get the uniform cubature points
     *
     * \note DEPRECATED Please use getCubaturePointsRef call instead
     * \note Option is only supported for volume integration, and even then, may not align with the getCubaturePointsRef call
     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     * \param[in] apply_permutation ADVANCED Do not change this unless you know what it does (it can break things)
     *
     * \return Array <point, dim>
     */
    ConstArray_IPDim
    getUniformCubaturePointsRef(const bool cache = true,
                                const bool force = false,
                                const bool apply_permutation = true) const;

    /**
     * \brief Get the uniform cubature points for a side
     *
     * \note DEPRECATED Please use getCubaturePointsRef call instead
     * \note Option is only supported for side integration, and even then, may not align with the getCubaturePointsRef call
     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     * \param[in] apply_permutation ADVANCED Do not change this unless you know what it does (it can break things)
     *
     * \return Array <point, dim>
     */
    ConstArray_IPDim
    getUniformSideCubaturePointsRef(const bool cache = true,
                                    const bool force = false,
                                    const bool apply_permutation = true) const;

    /**
     * \brief Get the uniform cubature weights
     *
     * \note DEPRECATED Please do not use
     * \note Option is only supported for some integration types
     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     * \param[in] apply_permutation ADVANCED Do not change this unless you know what it does (it can break things)
     *
     * \return Array <point>
     */
    ConstArray_IP
    getUniformCubatureWeightsRef(const bool cache = true,
                                 const bool force = false,
                                 const bool apply_permutation = true) const;


    /**
     * \brief Get the node coordinates describing the geometry of the mesh
     *
     * \return Array <cell, node, dim>
     */
    ConstArray_CellBASISDim
    getNodeCoordinates() const;

    /**
     * \brief Get the Jacobian matrix evaluated at the cubature points
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim, dim>
     */
    ConstArray_CellIPDimDim
    getJacobian(const bool cache = true,
                const bool force = false) const;


    /**
     * \brief Get the inverse of the Jacobian matrix evaluated at the cubature points
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim, dim>
     */
    ConstArray_CellIPDimDim
    getJacobianInverse(const bool cache = true,
                       const bool force = false) const;

    /**
     * \brief Get the determinant of the Jacobian matrix evaluated at the cubature points
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point>
     */
    ConstArray_CellIP
    getJacobianDeterminant(const bool cache = true,
                           const bool force = false) const;

    /**
     * \brief Get the weighted measure (integration weights)
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point>
     */
    ConstArray_CellIP
    getWeightedMeasure(const bool cache = true,
                       const bool force = false) const;

    /**
     * \brief Get the weighted normals
     *
     * \note Support: CV_SIDE
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim>
     */
    ConstArray_CellIPDim
    getWeightedNormals(const bool cache = true,
                       const bool force = false) const;

    /**
     * \brief Get the surface normals
     *
     * \note Support: SURFACE
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim>
     */
    ConstArray_CellIPDim
    getSurfaceNormals(const bool cache = true,
                      const bool force = false) const;

    /**
     * \brief Get the surface rotation matrices
     *
     * \note Support: SURFACE
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, 3, 3>
     */
    ConstArray_CellIPDimDim
    getSurfaceRotationMatrices(const bool cache = true,
                               const bool force = false) const;

    /**
     * \brief Get the covarient matrix
     *
     * cov(i,j) = jacobian(i,k) * jacobian(j,k)
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim, dim>
     */
    ConstArray_CellIPDimDim
    getCovarientMatrix(const bool cache = true,
                       const bool force = false) const;

    /**
     * \brief Get the contravarient matrix
     *
     * contra = (getCovarientMatrix())^{-1}
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim, dim>
     */
    ConstArray_CellIPDimDim
    getContravarientMatrix(const bool cache = true,
                           const bool force = false) const;

    /**
     * \brief Get the contravarient matrix
     *
     * norm = sqrt(\sum_{ij} cov(i,j) * cov(i,j))
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point>
     */
    ConstArray_CellIP
    getNormContravarientMatrix(const bool cache = true,
                               const bool force = false) const;

    /**
     * \brief Get the cubature points in physical space
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim>
     */
    ConstArray_CellIPDim
    getCubaturePoints(const bool cache = true,
                      const bool force = false) const;

    /**
     * \brief Get the cubature points in the reference space
     *
     * \note Support: VOLUME, SURFACE, SIDE, CV_VOLUME, CV_SIDE, CV_BOUNDARY
     *     *
     * \param[in] cache If true, the result will be stored in the IntegrationValues2 class
     * \param[in] force Force the re-evaluation of the array
     *
     * \return Array <cell, point, dim>
     */
    ConstArray_CellIPDim
    getCubaturePointsRef(const bool cache = true,
                         const bool force = false) const;

    /**
     * \brief Returns the IntegrationRule
     *
     * \return panzer::IntegrationRule
     */
    Teuchos::RCP<const panzer::IntegrationRule>
    getIntegrationRule() const
    {return int_rule;}

    // =====================================================================================================

  protected:

    // Reset all the lazy evaluation arrays
    void
    resetArrays();

    // Number of cells in mesh
    int num_cells_;

    // Number of cells in mesh to evaluate
    int num_evaluate_cells_;

    // Number of virtual cells in the mesh - used for surface evaluations
    int num_virtual_cells_;

    // Permutations (used to re-orient arrays similar to orientations in BasisValues2)
    bool requires_permutation_;

    // Array contains the mapping from uniform reference space to permuted space
    PHX::MDField<const int,Cell,IP> permutations_;

    // TODO: There is a way around this, but it will require some work
    // Subcell connectivity is required for surface evaluations (normals and rotation matrices)
    Teuchos::RCP<const SubcellConnectivity> side_connectivity_;

    // Lazy evaluation checks
    mutable bool cub_points_evaluated_;
    mutable bool side_cub_points_evaluated_;
    mutable bool cub_weights_evaluated_;
    mutable bool node_coordinates_evaluated_;
    mutable bool jac_evaluated_;
    mutable bool jac_inv_evaluated_;
    mutable bool jac_det_evaluated_;
    mutable bool weighted_measure_evaluated_;
    mutable bool weighted_normals_evaluated_;
    mutable bool surface_normals_evaluated_;
    mutable bool surface_rotation_matrices_evaluated_;
    mutable bool covarient_evaluated_;
    mutable bool contravarient_evaluated_;
    mutable bool norm_contravarient_evaluated_;
    mutable bool ip_coordinates_evaluated_;
    mutable bool ref_ip_coordinates_evaluated_;

    // Backward compatibility call that evaluates all internal values for CV, surface, side, or volume integration schemes
    void
    evaluateEverything();

  private:

    bool alloc_arrays_;
    std::string prefix_;
    std::vector<PHX::index_size_type> ddims_;

  };

} // namespace panzer

#endif
