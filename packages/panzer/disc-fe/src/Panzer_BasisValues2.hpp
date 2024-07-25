// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BASIS_VALUES2_HPP
#define PANZER_BASIS_VALUES2_HPP

#include "Teuchos_RCP.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_Orientation.hpp"

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_ArrayTraits.hpp"

namespace panzer {

  // For applying orientations with lazy evaluation
  class OrientationsInterface;

  /**
    * \class BasisValues2
    *
    * Data structure that holds all evaluated fields associated
    * with a basis function and integration rule. This class will
    * allocate the memory and evaluate the basis functions. The
    * orientations must be applied using the
    * <code>applyOrientations</code> method.
    *
    *
    * Lazy evaluation path:
    *
    * BasisValues can be defined in standard or uniform forms.
    * Uniform means that the reference space is the same for all cells,
    * while standard means that there is a different reference space per cell.
    * If the standard method is used, all get*Ref calls will throw errors.
    *
    * The construction path for the lazy evaluation form of the BasisValues2 object
    * is as follows:
    *
    * auto basis_values = Teuchos::rcp(new BasisValues<Scalar>("prefix"));
    *
    * # (Required) Main setup call for lazy evaluation
    * basis_values->setup[Uniform](basis, reference_points, jacobian, jacobian_determinant, jacobian_inverse);
    *
    * # (Optional) Some basis (HCurl/HDiv) require orientations
    * basis_values->setOrientations(orientations);
    *
    * # (Optional) If Scalar is a Sacado::Fad type, we need to know the dimension
    * basis_values->setExtendedDimensions(ddims);
    *
    * # (Optional) If you are going to integrate quantities, you'll need the weighted measure (see IntegrationValues2)
    * basis_values->setWeightedMeasure(weighted_measure);
    *
    * # basis_values is now ready for lazy evaluation
    *
    * Note that the above optional/required calls should happen in this order,
    * but all must happen before any 'get' calls can be called.
    *
    */
  template <typename Scalar>
  class BasisValues2 {
  public:
    typedef typename ArrayTraits<Scalar,PHX::MDField<Scalar> >::size_type size_type;
    using IntrepidBasis = Intrepid2::Basis<PHX::Device::execution_space,Scalar,Scalar>;

    typedef PHX::MDField<Scalar>                    ArrayDynamic;
    typedef PHX::MDField<Scalar,BASIS,IP>           Array_BasisIP;
    typedef PHX::MDField<Scalar,BASIS,IP,Dim>       Array_BasisIPDim;
    typedef PHX::MDField<Scalar,BASIS,Dim>          Array_BasisDim;
    typedef PHX::MDField<Scalar,Cell,BASIS,IP>      Array_CellBasisIP;
    typedef PHX::MDField<Scalar,Cell,BASIS,IP,Dim>  Array_CellBasisIPDim;
    typedef PHX::MDField<Scalar,Cell,BASIS,Dim>     Array_CellBasisDim;

    typedef PHX::MDField<const Scalar>                    ConstArrayDynamic;
    typedef PHX::MDField<const Scalar,BASIS,IP>           ConstArray_BasisIP;
    typedef PHX::MDField<const Scalar,BASIS,IP,Dim>       ConstArray_BasisIPDim;
    typedef PHX::MDField<const Scalar,BASIS,Dim>          ConstArray_BasisDim;
    typedef PHX::MDField<const Scalar,Cell,BASIS,IP>      ConstArray_CellBasisIP;
    typedef PHX::MDField<const Scalar,Cell,BASIS,IP,Dim>  ConstArray_CellBasisIPDim;
    typedef PHX::MDField<const Scalar,Cell,BASIS,Dim>     ConstArray_CellBasisDim;

    //=======================================================================================================
    // DEPRECATED Interface

    /**
     * \brief Main constructor
     *
     * \note This call must be followed with setupArrays
     *
     * \param[in] prefix Name to prefix all arrays with
     * \param[in] allocArrays If true we will allocate all arrays in the setupArrays call
     * \param[in] buildWeighted Builds the weighted components for the basis (i.e. multiplied by weighted measure for integration purposes)
     */
    BasisValues2(const std::string & prefix="",
                 const bool allocArrays=false,
                 const bool buildWeighted=false);

    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<const panzer::BasisIRLayout>& basis,
                     bool computeDerivatives=true);

    void evaluateValues(const PHX::MDField<Scalar,IP,Dim> & cub_points,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
                        const PHX::MDField<Scalar,Cell,IP> & jac_det,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
                        const int in_num_cells = -1);

    void evaluateValues(const PHX::MDField<Scalar,IP,Dim> & cub_points,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
                        const PHX::MDField<Scalar,Cell,IP> & jac_det,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
                        const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                        const PHX::MDField<Scalar,Cell,NODE,Dim> & node_coordinates,
                        bool use_node_coordinates=true,
                        const int in_num_cells = -1);

    void evaluateValuesCV(const PHX::MDField<Scalar,Cell,IP,Dim> & cell_cub_points,
                          const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
                          const PHX::MDField<Scalar,Cell,IP> & jac_det,
                          const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv);


    void evaluateValuesCV(const PHX::MDField<Scalar,Cell,IP,Dim> & cell_cub_points,
                          const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
                          const PHX::MDField<Scalar,Cell,IP> & jac_det,
                          const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
                          const PHX::MDField<Scalar,Cell,NODE,Dim> & node_coordinates,
                          bool use_node_coordinates=true,
                          const int in_num_cells = -1);


    void evaluateValues(const PHX::MDField<Scalar,Cell,IP,Dim> & cub_points,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
                        const PHX::MDField<Scalar,Cell,IP> & jac_det,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
                        const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                        const PHX::MDField<Scalar,Cell,NODE,Dim> & node_coordinates,
                        bool use_node_coordinates=true,
                        const int in_num_cells = -1);


    //! Method to apply orientations to a basis values container.
    // some evaluators use this apply orientation (this will be deprecated)
    void applyOrientations(const PHX::MDField<const Scalar,Cell,BASIS> & orientations);

    // this is used in workset factory
    void applyOrientations(const std::vector<Intrepid2::Orientation> & orientations,
                           const int in_num_cells = -1);

    void setExtendedDimensions(const std::vector<PHX::index_size_type> & ddims)
    { ddims_ = ddims; }

    PureBasis::EElementSpace getElementSpace() const;

    mutable Array_BasisIP     basis_ref_scalar;           // <BASIS,IP>
    mutable Array_CellBasisIP basis_scalar;               // <Cell,BASIS,IP>

    mutable Array_BasisIPDim     basis_ref_vector;           // <BASIS,IP,Dim>
    mutable Array_CellBasisIPDim basis_vector;               // <Cell,BASIS,IP,Dim>

    mutable Array_BasisIPDim     grad_basis_ref;             // <BASIS,IP,Dim>
    mutable Array_CellBasisIPDim grad_basis;                 // <Cell,BASIS,IP,Dim>

    mutable Array_BasisIP     curl_basis_ref_scalar;         // <BASIS,IP> - 2D!
    mutable Array_CellBasisIP curl_basis_scalar;             // <Cell,BASIS,IP> - 2D!

    mutable Array_BasisIPDim     curl_basis_ref_vector;      // <BASIS,IP,Dim> - 3D!
    mutable Array_CellBasisIPDim curl_basis_vector;          // <Cell,BASIS,IP,Dim> - 3D!

    mutable Array_BasisIP     div_basis_ref;           // <BASIS,IP>
    mutable Array_CellBasisIP div_basis;               // <Cell,BASIS,IP>

    mutable Array_CellBasisIP weighted_basis_scalar;                  // <Cell,BASIS,IP>
    mutable Array_CellBasisIPDim weighted_basis_vector;               // <Cell,BASIS,IP,Dim>
    mutable Array_CellBasisIPDim weighted_grad_basis;                 // <Cell,BASIS,IP,Dim>
    mutable Array_CellBasisIP weighted_curl_basis_scalar;             // <Cell,BASIS,IP>
    mutable Array_CellBasisIPDim weighted_curl_basis_vector;          // <Cell,BASIS,IP,Dim>
    mutable Array_CellBasisIP weighted_div_basis;                     // <Cell,BASIS,IP>

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    mutable Array_BasisDim basis_coordinates_ref;     // <BASIS,Dim>

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    mutable Array_CellBasisDim basis_coordinates;     // <Cell,BASIS,Dim>

    Teuchos::RCP<const panzer::BasisIRLayout> basis_layout;

    Teuchos::RCP<IntrepidBasis> intrepid_basis;

    bool compute_derivatives;
    bool build_weighted;
    bool alloc_arrays;
    std::string prefix;
    std::vector<PHX::index_size_type> ddims_;

    bool orientationsApplied() const
    {return orientations_applied_;}

    void evaluateBasisCoordinates(const PHX::MDField<Scalar,Cell,NODE,Dim> & node_coordinates,
                                  const int in_num_cells = -1);

  private:

    bool references_evaluated;

    // Have orientations been applied
    bool orientations_applied_;

    //=======================================================================================================
    // New Interface

  protected:

    // Reset internal data structures for lazy evaluation
    void resetArrays();

    // Number of total cells (owned + ghost + virtual)
    int num_cells_;

    // Number of cells to evaluate (HACK for backward compatability)
    int num_evaluate_cells_;

    // True if a uniform point space is used
    bool is_uniform_;

    // Only valid for uniform reference space
    PHX::MDField<const Scalar,IP,Dim>           cubature_points_uniform_ref_;

    // For non-uniform reference space
    PHX::MDField<const Scalar,Cell,IP,Dim>      cubature_points_ref_;

    // Geometry objects that represent cubature/point values
    PHX::MDField<const Scalar,Cell,IP,Dim,Dim>  cubature_jacobian_;
    PHX::MDField<const Scalar,Cell,IP>          cubature_jacobian_determinant_;
    PHX::MDField<const Scalar,Cell,IP,Dim,Dim>  cubature_jacobian_inverse_;
    PHX::MDField<const Scalar,Cell,IP>          cubature_weights_;

    // Coordinates of the mesh nodes
    PHX::MDField<const Scalar,Cell,NODE,Dim> cell_node_coordinates_;

    // Cell topology from the mesh
    Teuchos::RCP<const shards::CellTopology> cell_topology_;

    // Number of cells to apply orientations to (required in situations where virtual cells exist)
    int num_orientations_cells_;

    // Orientations object
    std::vector<Intrepid2::Orientation>   orientations_;

    /// Used to check if arrays have been cached
    mutable bool basis_ref_scalar_evaluated_;
    mutable bool basis_scalar_evaluated_;
    mutable bool basis_ref_vector_evaluated_;
    mutable bool basis_vector_evaluated_;
    mutable bool grad_basis_ref_evaluated_;
    mutable bool grad_basis_evaluated_;
    mutable bool curl_basis_ref_scalar_evaluated_;
    mutable bool curl_basis_scalar_evaluated_;
    mutable bool curl_basis_ref_vector_evaluated_;
    mutable bool curl_basis_vector_evaluated_;
    mutable bool div_basis_ref_evaluated_;
    mutable bool div_basis_evaluated_;
    mutable bool weighted_basis_scalar_evaluated_;
    mutable bool weighted_basis_vector_evaluated_;
    mutable bool weighted_grad_basis_evaluated_;
    mutable bool weighted_curl_basis_scalar_evaluated_;
    mutable bool weighted_curl_basis_vector_evaluated_;
    mutable bool weighted_div_basis_evaluated_;
    mutable bool basis_coordinates_ref_evaluated_;
    mutable bool basis_coordinates_evaluated_;

  public:

    /**
     * \brief Setup for lazy evaluation for non-uniform point layout
     *
     * \note This call is used when there are different reference points per cell
     *
     * \param[in] basis Basis layout - contains information for intrepid
     * \param[in] reference_points Points (e.g. cubature points) in the reference space of the cell
     * \param[in] point_jacobian Cell jacobian evaluated at the reference points
     * \param[in] point_jacobian_determinant Determinant of point_jacobian array
     * \param[in] point_jacobian_inverse Inverse of point_jacobian array
     * \param[in] num_evaluated_cells Used to force evaluation of arrays over subset of cells (default: all cells)
     *
     */
    void
    setup(const Teuchos::RCP<const panzer::BasisIRLayout> & basis,
          PHX::MDField<const Scalar, Cell, IP, Dim>         reference_points,
          PHX::MDField<const Scalar, Cell, IP, Dim, Dim>    point_jacobian,
          PHX::MDField<const Scalar, Cell, IP>              point_jacobian_determinant,
          PHX::MDField<const Scalar, Cell, IP, Dim, Dim>    point_jacobian_inverse,
          const int                                         num_evaluated_cells = -1);

    /**
     * \brief Setup for lazy evaluation for uniform point layout
     *
     * \note This call is used when the same reference points are used for each cell
     *
     * \param[in] basis Basis layout - contains information for intrepid
     * \param[in] reference_points Points (e.g. cubature points) in the reference space of the cell
     * \param[in] point_jacobian Cell jacobian evaluated at the reference points
     * \param[in] point_jacobian_determinant Determinant of point_jacobian array
     * \param[in] point_jacobian_inverse Inverse of point_jacobian array
     * \param[in] num_evaluated_cells Used to force evaluation of arrays over subset of cells (default: all cells)
     *
     */
    void
    setupUniform(const Teuchos::RCP<const panzer::BasisIRLayout> & basis,
                 PHX::MDField<const Scalar, IP, Dim>               reference_points,
                 PHX::MDField<const Scalar, Cell, IP, Dim, Dim>    point_jacobian,
                 PHX::MDField<const Scalar, Cell, IP>              point_jacobian_determinant,
                 PHX::MDField<const Scalar, Cell, IP, Dim, Dim>    point_jacobian_inverse,
                 const int                                         num_evaluated_cells = -1);

    /// Set the orientations object for applying orientations using the lazy evaluation path - required for certain bases
    void
    setOrientations(const std::vector<Intrepid2::Orientation> & orientations,
                    const int num_orientations_cells = -1);

    /// Set the cubature weights (weighted measure) for the basis values object - required to get weighted basis objects
    void
    setWeightedMeasure(PHX::MDField<const Scalar, Cell, IP> weighted_measure);

    /////////// TO BE DEPRECATED.....
    /// Set the cell vertex coordinates (required for getBasisCoordinates())
    void
    setCellVertexCoordinates(PHX::MDField<Scalar,Cell,NODE,Dim> vertex_coordinates);
    ////////// END TO BE DEPRECATED

    /// Set the cell node coordinates (required for getBasisCoordinates())
    void
    setCellNodeCoordinates(PHX::MDField<Scalar,Cell,NODE,Dim> node_coordinates);

    /// Check if reference point space is uniform across all cells (faster evaluation)
    bool
    hasUniformReferenceSpace() const
    {return is_uniform_;}

    /// Return the basis descriptor
    panzer::BasisDescriptor getBasisDescriptor() const;

    /// Get the extended dimensions used by sacado AD allocations
    const std::vector<PHX::index_size_type> &
    getExtendedDimensions() const
    {return ddims_;}

    //=======================================================================================================
    // Reference space evaluations

    /**
     * \brief Get the reference coordinates for basis
     *
     * \throws If reference points are not uniform
     * \throws If basis does not support a coordinate space
     *
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * return Array <basis, dim>
    */
    ConstArray_BasisDim
    getBasisCoordinatesRef(const bool cache = true,
                           const bool force = false) const;

    /**
     * \brief Get the basis values evaluated at reference points
     *
     * \throws If not a scalar basis
     * \throws If reference points are not uniform
     *
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <basis, point>
     */
    ConstArray_BasisIP
    getBasisValuesRef(const bool cache = true,
                      const bool force = false) const;

    /**
     * \brief Get the vector basis values evaluated at reference points
     *
     * \throws If not a vector basis
     * \throws If reference points are not uniform
     *
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <basis, point, dim>
     */
    ConstArray_BasisIPDim
    getVectorBasisValuesRef(const bool cache = true,
                            const bool force = false) const;

    /**
     * \brief Get the gradient of the basis evaluated at reference points
     *
     * \throws If not a scalar basis
     * \throws If reference points are not uniform
     *
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <basis, point, dim>
     */
    ConstArray_BasisIPDim
    getGradBasisValuesRef(const bool cache = true,
                          const bool force = false) const;

    /**
     * \brief Get the curl of a vector basis evaluated at reference points
     *
     * \throws If reference points are not uniform
     * \throws If not a 2D space
     * \throws If not a HCurl basis
     *
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <basis, point>
     */
    ConstArray_BasisIP
    getCurl2DVectorBasisRef(const bool cache = true,
                            const bool force = false) const;

    /**
     * \brief Get the curl of a vector basis evaluated at reference points
     *
     * \throws If reference points are not uniform
     * \throws If not a 3D space
     * \throws If not a HCurl basis
     *
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <basis, point, dim>
     */
    ConstArray_BasisIPDim
    getCurlVectorBasisRef(const bool cache = true,
                          const bool force = false) const;

    /**
     * \brief Get the divergence of a vector basis evaluated at reference points
     *
     * \throws If reference points are not uniform
     * \throws If not a 3D space
     * \throws If not a HDiv basis
     *
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <basis, point>
     */
    ConstArray_BasisIP
    getDivVectorBasisRef(const bool cache = true,
                         const bool force = false) const;

    //=======================================================================================================
    // Mesh evaluations

    /**
     * \brief Carterisan coordinates for basis coefficients in mesh space
     *
     * \throws If basis does not support a coordinate space
     *
     * return Array <cell, basis, dim>
    */
    ConstArray_CellBasisDim
    getBasisCoordinates(const bool cache = true,
                        const bool force = false) const;

    /**
     * \brief Get the basis values evaluated at mesh points
     *
     * \throws If not a scalar basis
     *
     * \param[in] weighted Add cubature weighting for integration purposes
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <cell, basis, point>
     */
    ConstArray_CellBasisIP
    getBasisValues(const bool weighted,
                   const bool cache = true,
                   const bool force = false) const;

    /**
     * \brief Get the vector basis values evaluated at mesh points
     *
     * \throws If not a vector basis
     *
     * \param[in] weighted Add cubature weighting for integration purposes
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <cell, basis, point, dim>
     */
    ConstArray_CellBasisIPDim
    getVectorBasisValues(const bool weighted,
                         const bool cache = true,
                         const bool force = false) const;

    /**
     * \brief Get the gradient of the basis evaluated at mesh points
     *
     * \throws If not a scalar basis
     *
     * \param[in] weighted Add cubature weighting for integration purposes
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <cell, basis, point, dim>
     */
    ConstArray_CellBasisIPDim
    getGradBasisValues(const bool weighted,
                       const bool cache = true,
                       const bool force = false) const;

    /**
     * \brief Get the curl of a 2D vector basis evaluated at mesh points
     *
     * \throws If not a HCurl basis
     * \throws If not a 2D space
     *
     * \param[in] weighted Add cubature weighting for integration purposes
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <cell, basis, point>
     */
    ConstArray_CellBasisIP
    getCurl2DVectorBasis(const bool weighted,
                         const bool cache = true,
                         const bool force = false) const;

    /**
     * \brief Get the curl of a 3D vector basis evaluated at mesh points
     *
     * \throws If not a HCurl basis
     * \throws If not a 3D space
     *
     * \param[in] weighted Add cubature weighting for integration purposes
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <cell, basis, point, dim>
     */
    ConstArray_CellBasisIPDim
    getCurlVectorBasis(const bool weighted,
                       const bool cache = true,
                       const bool force = false) const;

    /**
     * \brief Get the divergence of a vector basis evaluated at mesh points
     *
     * \throws If not a HDiv basis
     * \throws If not a 3D space
     *
     * \param[in] weighted Add cubature weighting for integration purposes
     * \param[in] cache If true, the returned object will be cached for later use
     * \param[in] force Force re-evaluation of cached array
     *
     * \return Array <cell, basis, point>
     */
    ConstArray_CellBasisIP
    getDivVectorBasis(const bool weighted,
                      const bool cache = true,
                      const bool force = false) const;

    //=======================================================================================================

  };

} // namespace panzer

#endif
