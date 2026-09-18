// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_GEOMETRICONELEVELPRECONDITIONER_DECL_HPP
#define _FROSCH_GEOMETRICONELEVELPRECONDITIONER_DECL_HPP

#include <FROSch_SchwarzPreconditioner_decl.hpp>
#include <FROSch_SchwarzPreconditioner_def.hpp>
#include <FROSch_GeometricOverlappingOperator_decl.hpp>
#include <FROSch_GeometricOverlappingOperator_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;    

    /**
     * @class GeometricOneLevelPreconditioner
     *
     * The GeometricOneLevelPreconditioner is a variant of the OneLevelPrecondtioner that 
     * builds on the GeometricOverlappingOperator instead of the AlgebraicOverlappingOperator.
     *
     * In the construction of the preconditioner, the algebraic decomposition of the matrix 
     * itself always results in Dirichlet boundary conditions. The local system matrix must 
     * be assembled separately to apply different boundary conditions in the preconditioner, 
     * e.g., Robin, perfectly matched layer, or other absorbing boundary conditions.  This 
     * is precisely what we do here: we employ the locally assembled system matrix to allow 
     * for different interface conditions.
     *
     * To construct a GeometricOneLevelPreconditioner preconditioner, a list of vertices, a 
     * list of cells (where 2*dim vertices describe each cell), and an auxiliary list of the 
     * same length as the cell list containing additional information, such as material 
     * constants is required. This information must be provided by a finite element program. 
     * Where each rank only stores all vertices, cells, and auxiliary data relevant to itself. 
     * Also, the node graph of how the vertices are connected has to be provided from the 
     * finite element program. Please note that the parallel distributed triangulation is 
     * assumed to correspond to a non-overlapping domain decomposition. Based on the dual 
     * graph, an overlapping list of vertices, cells, and auxiliary information is computed. 
     * Based on the returned overlapping lists, the local matrices must be constructed by the 
     * finite element program in use and provided to this class.
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
    class GeometricOneLevelPreconditioner : public SchwarzPreconditioner<SC,LO,GO,NO> {

    protected:

        using XMapPtr                           = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                      = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtr;

        using GraphPtr                          = RCP<Xpetra::CrsGraph<LO,GO,NO>>;

        using XMatrixPtr                        = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector                      = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVectorPtr;

        template <typename S>
        using XMultiVectorTemplate              = Xpetra::MultiVector<S,LO,GO,NO>;
        template <typename S>
        using XMultiVectorTemplatePtr           = RCP<XMultiVectorTemplate<S>>;

        using ParameterListPtr                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr;

        using SumOperatorPtr                    = typename SchwarzPreconditioner<SC,LO,GO,NO>::SumOperatorPtr;
        using MultiplicativeOperatorPtr         = typename SchwarzPreconditioner<SC,LO,GO,NO>::MultiplicativeOperatorPtr;

        using OverlappingOperatorPtr            = typename SchwarzPreconditioner<SC,LO,GO,NO>::OverlappingOperatorPtr;
        using GeometricOverlappingOperatorPtr   = RCP<GeometricOverlappingOperator<SC,LO,GO,NO>>;

    public:

        /**
         * @brief Constructor for the GeometricOneLevelPreconditioner class.
         *
         * This is the default constructor that takes the (global) system_matrix  and the
         * dual graph to construct the GeometricOverlappingOperator.
         *
         * To fine-tune the resulting preconditioner, several options can be controlled 
         * via the provided parameter list:
         * - **OverlappingOperator Type**: The type of overlapping operator to use. Options are: GeometricOverlappingOperator. 
         * - **Level Combination**: The level combination method. Options are: Additive, Multiplicative. Default Additive.
         * - **Overlap**: The overlap for the Schwarz preconditioner. This should be a positive integer. Default 0.
         *
         * And the Optimized Operator can be fine-tuned by the following options:
         * - **Combine Values in Overlap**: The method to combine values in overlap. Options are: Restricted, Averaging, Full. Default Restricted.
         * - **Solver type**: The type of solver to use. The option is: Amesos2. 
         * - **Solver**: The specific solver to use. The options are all direct solvers from Amesos2. Default KLU2
         *
         * @param k The (global) system_matrix.
         * @param dualGraph The dual graph of the grid.
         * @param overlap The overlap for the Schwarz preconditioner.
         * @param parameterList The list of parameters for the GeometricOverlappingOperator.
         */
        GeometricOneLevelPreconditioner(ConstXMatrixPtr  k,
                                        GraphPtr         dualGraph,
                                        ParameterListPtr parameterList);

        /**
         * Not implemented. 
         */
        virtual int 
        initialize(bool useDefault) override;

        /**
         * @brief Initializes the GeometricOneLevelPreconditioner based on the given overlap and the overlappingMap.
         *
         * This function creates the overlapping map for the cells. The map is computed based 
         * on the given overlap and the dual graph provided in this class's constructor. 
         * Which is done by calling the initialize function from the underlying GeometricOverlappingOperator.
         *
         * @param overlappingMap The overlapping dof map.
         * @return Error code. 0 if successful.
         */
        int 
        initialize(XMapPtr overlappingMap);

        /**
         * Not implemented. 
         */
        virtual int 
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
        compute(ConstXMatrixPtr neumannMatrix, 
                ConstXMatrixPtr robinMatrix);

        
        virtual void 
        apply(const XMultiVector &x,
              XMultiVector &y,
              ETransp       mode  = NO_TRANS,
              SC            alpha = ScalarTraits<SC>::one(),
              SC            beta  = ScalarTraits<SC>::zero()) const override;

        /**
         * @brief Return the domain map.
         * 
         * @return The domain map
         */
        virtual const ConstXMapPtr 
        getDomainMap() const override;

        /**
         * @brief Return the range map.
         * 
         * @return The range map
         */
        virtual const ConstXMapPtr 
        getRangeMap() const override;

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
         * @brief Prints the description of the GeometricOneLevelPreconditioner to a stream.
         *
         * This function prints a detailed description of the GeometricOneLevelPreconditioner 
         * to the provided output stream. The level of detail of the description is 
         * controlled by the verbLevel parameter.
         *
         * @param out The output stream to which the description is printed.
         * @param verbLevel The verbosity level controlling the detail of the description.
         */
        virtual void 
        describe(
            FancyOStream          &out,
            const EVerbosityLevel  verbLevel = Describable::verbLevel_default) const override;

        /**
         * @brief This function returns a string with the name of the operator.
         * 
         * @return The name of the operator.
         */
        virtual string 
        description() const override;

        virtual int 
        resetMatrix(ConstXMatrixPtr &k);

    protected:

        ConstXMatrixPtr           K_;

        SumOperatorPtr            SumOperator_;
        MultiplicativeOperatorPtr MultiplicativeOperator_;
        OverlappingOperatorPtr    OverlappingOperator_;
        bool                      UseMultiplicative_ = false;
    };

}

#endif // _FROSCH_GEOMETRICONElEVELPRECONDITIONER_DECL_HPP


