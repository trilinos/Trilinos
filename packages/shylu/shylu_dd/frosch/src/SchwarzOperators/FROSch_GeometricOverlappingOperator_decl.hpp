// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_GEOMETRICOVERLAPPINGOPERATOR_DECL_HPP
#define _FROSCH_GEOMETRICOVERLAPPINGOPERATOR_DECL_HPP

#include <FROSch_OverlappingOperator_def.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;
    /**
     * @class GeometricOverlappingOperator
     *
     * The GeometricOverlappingOperator is derived from the OverlappingOperator. 
     * It is used to build a geometric (optimized) Schwarz preconditioner.  
     *
     * In the construction of the preconditioner, the algebraic decomposition of the matrix 
     * itself always results in Dirichlet boundary conditions. The local system matrix must 
     * be assembled separately to apply different boundary conditions in the preconditioner, 
     * e.g., Robin, perfectly matched layer, or other absorbing boundary conditions.  This 
     * is precisely what we do here: we employ the locally assembled system matrix to allow 
     * for different interface conditions.
     *
     * This backbone class is designed for use in the OneLevelOptimizedPreconditioner 
     * or TwoLevelOptimizedPreconditioner Preconditioner class. It handles the creation 
     * of overlapping subproblems (based on some global problem description) and provides 
     * tools to create local subproblems in other finite element software (e.g., deal.II).
     *
     * @tparam SC The scalar type. (default is double)
     * @tparam LO The local ordinal type. (default is int)
     * @tparam GO The global ordinal type. (default is DefaultGlobalOrdinal)
     * @tparam NO The node type. (default is Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
     */
    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class GeometricOverlappingOperator : public OverlappingOperator<SC,LO,GO,NO> {

    protected:
        using CommPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMultiVector            = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr         = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;

        template <typename S>
        using XMultiVectorTemplate    = Xpetra::MultiVector<S,LO,GO,NO>;
        template <typename S>
        using XMultiVectorTemplatePtr = RCP<XMultiVectorTemplate<S>>;

        using GraphPtr                = typename SchwarzOperator<SC,LO,GO,NO>::GraphPtr;

        using XMapPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr            = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;

        using XMatrix                 = typename SchwarzOperator<SC,LO,GO,NO>::XMatrix;
        using XMatrixPtr              = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr         = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using ConstXCrsGraphPtr       = typename SchwarzOperator<SC,LO,GO,NO>::ConstXCrsGraphPtr;

        using ParameterListPtr        = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

    public:

        /**
         * @brief Constructor for the GeometricOverlappingOperator class.
         *
         * This is the default constructor that takes the (global) system_matrix 
         * to construct the GeometricOverlappingOperator.
         *
         * To fine-tune the resulting preconditioner, several options can be controlled 
         * via the provided parameter list:
         * - **Combine Values in Overlap**: The method to combine values in overlap. Options are: Restricted, Averaging, Full. Default Restricted.
         * - **Solver type**: The type of solver to use. The option is: Amesos2. 
         * - **Solver**: The specific solver to use. The options are all direct solvers from Amesos2. Default KLU2
         *
         * @param k The (global) system_matrix.
         * @param dualGraph The dual graph of the grid.
         * @param overlap The overlap for the Schwarz preconditioner.
         * @param parameterList The list of parameters for the GeometricOverlappingOperator.
         */
        GeometricOverlappingOperator(ConstXMatrixPtr  k,
                                     int              overlap,
                                     GraphPtr         dualGraph,
                                     ParameterListPtr parameterList);

        /**
         * Not implemented. 
         */
        virtual int 
        initialize() override;

        /**
         * @brief Initializes the GeometricOverlappingOperator based on the given overlap and the overlappingMap.
         *
         * This function creates the overlapping map for the cells (via calling 
         * OverlappingOperator::buildOverlappingMap()). The map is computed based on the given 
         * overlap and the dual graph, that was provided in the constructor of this class. 
         * Moreover, this function assigns the overlapping DoF map to the internal stores it in 
         * OverlappingMap_. 
         *
         * After assigning the overlapping dof map to this->OverlappingMap_, the other functions, 
         * that are normally called in initialize are called (i.e. 
         * OverlappingOperator::initializeOverlappingOperator() and 
         * OverlappingOperator::updateLocalOverlappingMatrices_Symbolic).
         *
         * @param overlappingMap The overlapping dof map.
         * @return Error code. 0 if successful.
         */
        int 
        initialize(XMapPtr overlappingMap);

        /**
         * Not implemented. 
         */
        int 
        compute() override;

        /**
         * @brief Compute the preconditioner.
         *
         * This function actually computes the preconditioner. Therefore, it computes the 
         * inverse of the subproblem matrix, which is constructed based on the NeumannMatrix 
         * (the local system_matrix) and the RobinMatrix (which describes the local optimized 
         * interface conditions).
         *
         * @warning This function is computationally expensive, as it computes the inverse of 
         * all matrices on the subproblems.
         *
         * @param neumannMatrix The local system matrix.
         * @param robinMatrix The matrix containing the optimized interface conditions.
         * @return Error code. 0 if successful.
         */
        int 
        compute(ConstXMatrixPtr neumannMatrix, ConstXMatrixPtr robinMatrix);

        /**
         * @brief Redistributes the nodes_vector, cell_vector and auxillary_vector onto the overlapping domain.
         *
         * This function takes an Xpetra::MultiVector, that stores a list of nodes
         * (where each node is given by its coordinates), a Xpetra::MultiVector 
         * that contains a list of cells (where each cell is described by its vertices)
         * and an Xpetra::MultiVector auxillary list that contains any further important 
         * infromation.
         *
         * Those vectors redistributes onto the overlapping domain. 
         * Thhose vectors are redistributed based on the overlapping cell map computed 
         * in intialize(). But we have to take special care of the node list, as we
         * can not directly use the overlapping cell map but have to compute a new map
         * based on the overlapping cell map and the content of the cell_vector.
         *
         * 
         * Moreover, we have to take special care of the cell_vector. As explained above, 
         * each cell is described by its vertices. The vertices are supplied as the number 
         * of the vertex. Therefore, this function also takes care to update those vertex 
         * numbers inside the cell_vector to match the the newly created nodes_vector.
         *
         * @param nodeList The nodes_vector (dim, <number vertices>).
         * @param elementList The cell_vector (<vertices_per_cell>, <number cells>).
         * @param auxillaryList The auxillary_vector (<arbitrary>, <number cells>).
         * @param nodeListOverlapping The redistributed nodes_vector.
         * @param elementListOverlapping The redistributed cell_vector.
         * @param auxillaryListOverlapping The redistributed auxillary_vector.
         * @return Error code. 0 if successful.
         */
        int 
        communicateOverlappingTriangulation(
            XMultiVectorPtr                     nodeList,
            XMultiVectorTemplatePtr<long long>  elementList,
            XMultiVectorTemplatePtr<long long>  auxillaryList,
            XMultiVectorPtr                    &nodeListOverlapping,
            XMultiVectorTemplatePtr<long long> &elementListOverlapping,
            XMultiVectorTemplatePtr<long long> &auxillaryListOverlapping);

        /**
         * @brief Prints the description of the GeometricOverlappingOperator to a stream.
         *
         * This function prints a detailed description of the GeometricOverlappingOperator 
         * to the provided output stream. The level of detail of the description is 
         * controlled by the verbLevel parameter.
         *
         * @param out The output stream to which the description is printed.
         * @param verbLevel The verbosity level controlling the detail of the description.
         */
        void 
        describe(FancyOStream          &out,
                 const EVerbosityLevel  verbLevel) const override;

        /**
         * @brief This function returns a string with the name of the operator.
         * 
         * @return The name of the operator.
         */
        string 
        description() const override;

    protected:
      /**
       * @brief Computes the subdomains matrices based on the system_matrix and the overlapping dof map.
       *
       * @param overlap The overlap for the Schwarz preconditioner. Must be >= 1.
       * @return Error code. 0 if successful.
       */
      int 
      buildOverlappingMap(int overlap);
    
      /**
       * @brief Extracts the local subdomain matrices from the system_matrix when it was not done earlier.
       *
       * This function is called by buildOverlappingMap. It also adds the local 
       * robin matrix, that stores the optimized interface conditions onto the 
       * subdomain matrix.
       *
       * @return Error code. 0 if successful.
       */
      int 
      updateLocalOverlappingMatrices() override;
    
      /**
       * @brief An internal function called by updateLocalOverlappingMatrices.
       *
       * @return Error code. 0 if successful.
       */
      int 
      updateLocalOverlappingMatrices_Symbolic();
    
      /**
       * @brief An internal function called by updateLocalOverlappingMatrices_Symbolic.
       */
      void 
      extractLocalSubdomainMatrix_Symbolic();
    
    private:
      /**
       * @brief A RCP<Xpetra::CrsGraph<SC,LO,GO,NO>> which contains the dual graph.
       *
       * If there is an entry in (row i, column j), element i and element j are neighbors.
       */
      GraphPtr 
      DualGraph_;
    
      /**
       * @brief The DualGraph but with overlap.
       */
      ConstXCrsGraphPtr 
      OverlappingGraph_;
    
      /**
       * @brief The column map of the DualGraph with overlap.
       */
      ConstXMapPtr 
      OverlappingElementMap_;
    
      /**
       * @brief Local matrix that contains the local system matrix.
       */
      ConstXMatrixPtr 
      NeumannMatrix_;
    
      /**
       * @brief Local matrix that contains the optimized interface conditions.
       */
      ConstXMatrixPtr 
      RobinMatrix_;

    };

} //namespace FROSch

#endif // _FROSCH_GEOMETRICOVERLAPPINGOPERATOR_DECL_HPP

