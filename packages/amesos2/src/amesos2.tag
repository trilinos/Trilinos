<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>Amesos2 - Direct Sparse Solver Interfaces</title>
    <filename>index</filename>
    <docanchor file="index">Amesos2_intro</docanchor>
    <docanchor file="index">Amesos2_outline</docanchor>
    <docanchor file="index">Amesos2_sclasses</docanchor>
    <docanchor file="index">Amesos2_aclasses</docanchor>
    <docanchor file="index">Amesos2_startup</docanchor>
    <docanchor file="index">Amesos2_contributors</docanchor>
  </compound>
  <compound kind="file">
    <name>Amesos2_Control.cpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Control_8cpp</filename>
    <includes id="Amesos2__Control_8hpp" name="Amesos2_Control.hpp" local="yes" imported="no">Amesos2_Control.hpp</includes>
  </compound>
  <compound kind="file">
    <name>Amesos2_Control.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Control_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_EpetraCrsMatrix_MatrixAdapter_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__EpetraCrsMatrix__MatrixAdapter__decl_8hpp</filename>
    <class kind="class">Amesos::ConcreteMatrixAdapter&lt; Epetra_CrsMatrix &gt;</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_EpetraMultiVecAdapter_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__EpetraMultiVecAdapter__decl_8hpp</filename>
    <class kind="class">Amesos::MultiVecAdapter&lt; Epetra_MultiVector &gt;</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_EpetraMultiVecAdapter_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__EpetraMultiVecAdapter__def_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__EpetraRowMatrix__AbstractMatrixAdapter__decl_8hpp</filename>
    <class kind="class">Amesos::AbstractConcreteMatrixAdapter&lt; Epetra_RowMatrix, DerivedMat &gt;</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__EpetraRowMatrix__AbstractMatrixAdapter__def_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Factory_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Factory__decl_8hpp</filename>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>a1bfd24259d3d02f48d770f39663a1416</anchor>
      <arglist>(const char *solverName)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>afcec03d38bbe79c3170ef377cd0db841</anchor>
      <arglist>(const std::string solverName)</arglist>
    </member>
    <docanchor file="Amesos2__Factory__decl_8hpp">usage</docanchor>
  </compound>
  <compound kind="file">
    <name>Amesos2_Factory_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Factory__def_8hpp</filename>
    <member kind="function">
      <type>std::string</type>
      <name>tolower</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>a77f0e3a2a428dc8416850e8f46f241c8</anchor>
      <arglist>(const std::string &amp;s)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>a0d5ff6f140807099a30e6a7f050039e0</anchor>
      <arglist>(Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>a9093713df97505f28528b80feecec3ba</anchor>
      <arglist>(RCP&lt; Matrix &gt; A, RCP&lt; Vector &gt; X, RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>a0415713b95269067132218d80225973a</anchor>
      <arglist>(const char *solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>aec1fc3fa8691356a9e3f3d22a159d3e3</anchor>
      <arglist>(const char *solverName, const RCP&lt; Matrix &gt; A, const RCP&lt; Vector &gt; X, const RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>a18abf3bcfa605087216d8c98f74b5d47</anchor>
      <arglist>(const std::string solverName, Matrix *A)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>ad75304da50137ecefca16143e4d0b4ce</anchor>
      <arglist>(const std::string solverName, const RCP&lt; Matrix &gt; A)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>afac0e3332344c53172545fc7baa57499</anchor>
      <arglist>(const std::string solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>afc97da118080236a3665fc0621be20a6</anchor>
      <arglist>(const std::string solver_name, const RCP&lt; Matrix &gt; A, const RCP&lt; Vector &gt; X, const RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>a1bfd24259d3d02f48d770f39663a1416</anchor>
      <arglist>(const char *solverName)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>query</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>afcec03d38bbe79c3170ef377cd0db841</anchor>
      <arglist>(const std::string solverName)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Amesos2_FunctionMap.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__FunctionMap_8hpp</filename>
    <class kind="struct">Amesos::FunctionMap</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_MultiVecAdapter.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__MultiVecAdapter_8hpp</filename>
    <includes id="Amesos2__TpetraMultiVecAdapter__decl_8hpp" name="Amesos2_TpetraMultiVecAdapter_decl.hpp" local="yes" imported="no">Amesos2_TpetraMultiVecAdapter_decl.hpp</includes>
    <includes id="Amesos2__EpetraMultiVecAdapter__decl_8hpp" name="Amesos2_EpetraMultiVecAdapter_decl.hpp" local="yes" imported="no">Amesos2_EpetraMultiVecAdapter_decl.hpp</includes>
    <includes id="Amesos2__TpetraMultiVecAdapter__def_8hpp" name="Amesos2_TpetraMultiVecAdapter_def.hpp" local="yes" imported="no">Amesos2_TpetraMultiVecAdapter_def.hpp</includes>
    <includes id="Amesos2__EpetraMultiVecAdapter__def_8hpp" name="Amesos2_EpetraMultiVecAdapter_def.hpp" local="yes" imported="no">Amesos2_EpetraMultiVecAdapter_def.hpp</includes>
  </compound>
  <compound kind="file">
    <name>Amesos2_Solver_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Solver__decl_8hpp</filename>
    <includes id="Amesos2__Control_8hpp" name="Amesos2_Control.hpp" local="yes" imported="no">Amesos2_Control.hpp</includes>
    <includes id="Amesos2__Status_8hpp" name="Amesos2_Status.hpp" local="yes" imported="no">Amesos2_Status.hpp</includes>
    <includes id="Amesos2__Timers_8hpp" name="Amesos2_Timers.hpp" local="yes" imported="no">Amesos2_Timers.hpp</includes>
    <class kind="class">Amesos::Solver</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_SolverCore_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__SolverCore__decl_8hpp</filename>
    <includes id="Amesos2__MultiVecAdapter_8hpp" name="Amesos2_MultiVecAdapter.hpp" local="yes" imported="no">Amesos2_MultiVecAdapter.hpp</includes>
    <class kind="class">Amesos::SolverCore</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_SolverCore_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__SolverCore__def_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Status.cpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Status_8cpp</filename>
    <includes id="Amesos2__Status_8hpp" name="Amesos2_Status.hpp" local="yes" imported="no">Amesos2_Status.hpp</includes>
  </compound>
  <compound kind="file">
    <name>Amesos2_Status.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Status_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlu_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlu__decl_8hpp</filename>
    <includes id="Amesos2__Superlu__FunctionMap_8hpp" name="Amesos2_Superlu_FunctionMap.hpp" local="yes" imported="no">Amesos2_Superlu_FunctionMap.hpp</includes>
    <includes id="Amesos2__Superlu__TypeMap_8hpp" name="Amesos2_Superlu_TypeMap.hpp" local="yes" imported="no">Amesos2_Superlu_TypeMap.hpp</includes>
    <class kind="class">Amesos::Superlu</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlu_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlu__def_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlu_FunctionMap.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlu__FunctionMap_8hpp</filename>
    <includes id="Amesos2__FunctionMap_8hpp" name="Amesos2_FunctionMap.hpp" local="yes" imported="no">Amesos2_FunctionMap.hpp</includes>
    <includes id="Amesos2__Superlu__TypeMap_8hpp" name="Amesos2_Superlu_TypeMap.hpp" local="yes" imported="no">Amesos2_Superlu_TypeMap.hpp</includes>
    <class kind="struct">Amesos::FunctionMap&lt; Superlu, Scalar &gt;</class>
    <member kind="typedef">
      <type>int</type>
      <name>int_t</name>
      <anchorfile>namespaceSLU.html</anchorfile>
      <anchor>a63ea4864076bf60a5fe72f332e75b5ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; SLU::SuperMatrix &gt;</type>
      <name>getPermMatrix</name>
      <anchorfile>namespaceSLU.html</anchorfile>
      <anchor>a2cb2fad27c774aeb6456fecbb15f60ad</anchor>
      <arglist>(SLU::superlu_options_t *options, SLU::SuperMatrix *A, int *perm_c, int *etree)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; SLU::SuperMatrix &gt;</type>
      <name>getPermMatrix</name>
      <anchorfile>namespaceSLU.html</anchorfile>
      <anchor>aabd538f07da42f3c66069dfe8472b02e</anchor>
      <arglist>(superlu_options_t *options, SuperMatrix *A, int *perm_c, int *etree)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlu_TypeMap.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlu__TypeMap_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superludist_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superludist__decl_8hpp</filename>
    <includes id="Amesos2__Superludist__FunctionMap_8hpp" name="Amesos2_Superludist_FunctionMap.hpp" local="yes" imported="no">Amesos2_Superludist_FunctionMap.hpp</includes>
    <class kind="class">Amesos::Superludist</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superludist_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superludist__def_8hpp</filename>
    <member kind="function">
      <type>int</type>
      <name>sp_colorder</name>
      <anchorfile>namespaceSLUD.html</anchorfile>
      <anchor>a0e0667a5a4a7c7d5edf5b0a5a1a0f46f</anchor>
      <arglist>(SuperMatrix *, int *, superludist_options_t *, SuperMatrix *)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superludist_FunctionMap.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superludist__FunctionMap_8hpp</filename>
    <includes id="Amesos2__FunctionMap_8hpp" name="Amesos2_FunctionMap.hpp" local="yes" imported="no">Amesos2_FunctionMap.hpp</includes>
    <includes id="Amesos2__Superludist__TypeMap_8hpp" name="Amesos2_Superludist_TypeMap.hpp" local="yes" imported="no">Amesos2_Superludist_TypeMap.hpp</includes>
    <class kind="struct">Amesos::FunctionMap&lt; Superludist, Scalar &gt;</class>
    <member kind="function">
      <type>SLUD::DiagScale_t</type>
      <name>get_diag_scale</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>ad37e7830fcfa9ab9354ff0ae3b64f6bf</anchor>
      <arglist>(char eq)</arglist>
    </member>
    <member kind="function">
      <type>char</type>
      <name>get_equed</name>
      <anchorfile>namespaceAmesos.html</anchorfile>
      <anchor>aa354aab097b75818dc67471a0850f599</anchor>
      <arglist>(SLUD::DiagScale_t ds)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superludist_TypeMap.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superludist__TypeMap_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlumt_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlumt__decl_8hpp</filename>
    <includes id="Amesos2__Superlumt__FunctionMap_8hpp" name="Amesos2_Superlumt_FunctionMap.hpp" local="yes" imported="no">Amesos2_Superlumt_FunctionMap.hpp</includes>
    <class kind="class">Amesos::Superlumt</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlumt_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlumt__def_8hpp</filename>
    <member kind="function">
      <type>int</type>
      <name>sp_colorder</name>
      <anchorfile>namespaceSLUMT.html</anchorfile>
      <anchor>ae5643a218ce061424cb3ae4534697045</anchor>
      <arglist>(SuperMatrix *, int *, superlumt_options_t *, SuperMatrix *)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlumt_FunctionMap.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlumt__FunctionMap_8hpp</filename>
    <includes id="Amesos2__FunctionMap_8hpp" name="Amesos2_FunctionMap.hpp" local="yes" imported="no">Amesos2_FunctionMap.hpp</includes>
    <includes id="Amesos2__Superlumt__TypeMap_8hpp" name="Amesos2_Superlumt_TypeMap.hpp" local="yes" imported="no">Amesos2_Superlumt_TypeMap.hpp</includes>
    <class kind="struct">Amesos::FunctionMap&lt; Superlumt, Scalar &gt;</class>
    <member kind="typedef">
      <type>int</type>
      <name>int_t</name>
      <anchorfile>namespaceSLUMT.html</anchorfile>
      <anchor>ade29f670d027f9d24d1dc29ec09d3eca</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Amesos2_Superlumt_TypeMap.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Superlumt__TypeMap_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Timers.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Timers_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__TpetraCrsMatrix__MatrixAdapter__decl_8hpp</filename>
    <class kind="class">Amesos::ConcreteMatrixAdapter&lt; Tpetra::CrsMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps &gt; &gt;</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_TpetraMultiVecAdapter_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__TpetraMultiVecAdapter__decl_8hpp</filename>
    <class kind="class">Amesos::MultiVecAdapter&lt; Tpetra::MultiVector&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt; &gt;</class>
  </compound>
  <compound kind="file">
    <name>Amesos2_TpetraMultiVecAdapter_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__TpetraMultiVecAdapter__def_8hpp</filename>
  </compound>
  <compound kind="file">
    <name>Amesos2_Util_decl.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Util__decl_8hpp</filename>
    <class kind="struct">Amesos::Util::has_special_impl</class>
    <class kind="struct">Amesos::Util::row_access</class>
    <class kind="struct">Amesos::Util::col_access</class>
    <class kind="struct">Amesos::Util::get_cxs_helper</class>
    <class kind="struct">Amesos::Util::get_ccs_helper</class>
    <class kind="struct">Amesos::Util::get_crs_helper</class>
    <member kind="enumvalue">
      <name>Distributed</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>aad35f806a27308be6fa12ea9acbaaa88a4343856a900ea33811ed7e3c49fe8c54</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Distributed_No_Overlap</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>aad35f806a27308be6fa12ea9acbaaa88a68fc183f38d9d4791403a1b7e99d8599</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Globally_Replicated</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>aad35f806a27308be6fa12ea9acbaaa88af10754c1b1f6de9ba1248441ec06ffb3</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Rooted</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>aad35f806a27308be6fa12ea9acbaaa88acf2b74303756800e4f337d8dccd54643</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Sorted_Indices</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>ab4b7c1000135b83dd295aae5989072bda758a1e31d43985250b31e58dcf9093d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Arbitrary</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>ab4b7c1000135b83dd295aae5989072bda2c668209535846ba0fc51bab4452d2b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Tpetra::Map&lt; LO, GO, Node &gt; &gt;</type>
      <name>epetra_map_to_tpetra_map</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>aed570a2f70ea3395c64a75f1c27bcdca</anchor>
      <arglist>(const Epetra_BlockMap &amp;map)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Epetra_Map &gt;</type>
      <name>tpetra_map_to_epetra_map</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a93e07b55ec4d21a1f9e72e2c2ef30d7b</anchor>
      <arglist>(const Tpetra::Map&lt; LO, GO, Node &gt; &amp;map)</arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>to_teuchos_comm</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a12090631546a38cf01e3f2a1b7f19aa0</anchor>
      <arglist>(RCP&lt; const Epetra_Comm &gt; c)</arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Epetra_Comm &gt;</type>
      <name>to_epetra_comm</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a8c3a8f5b260a9ee1c8869d857e838f63</anchor>
      <arglist>(RCP&lt; const Teuchos::Comm&lt; int &gt; &gt; c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>transpose</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>ad1b59ab1d2ed7ddb1a4e2cc23167e3f7</anchor>
      <arglist>(ArrayView&lt; Scalar &gt; vals, ArrayView&lt; GlobalOrdinal &gt; indices, ArrayView&lt; GlobalSizeT &gt; ptr, ArrayView&lt; Scalar &gt; trans_vals, ArrayView&lt; GlobalOrdinal &gt; trans_indices, ArrayView&lt; GlobalSizeT &gt; trans_ptr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scale</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a8e72c92cd7c3935d1ed8f7224486396b</anchor>
      <arglist>(ArrayView&lt; Scalar1 &gt; vals, size_t l, size_t ld, ArrayView&lt; Scalar2 &gt; s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scale</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a43449c000c2761790f50353837945479</anchor>
      <arglist>(ArrayView&lt; Scalar1 &gt; vals, size_t l, size_t ld, ArrayView&lt; Scalar2 &gt; s, BinaryOp binary_op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeTrueResidual</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a264cdfa2318fd49c83041158fd4af1cf</anchor>
      <arglist>(const RCP&lt; Matrix &gt; &amp;A, const RCP&lt; Vector &gt; &amp;X, const RCP&lt; Vector &gt; &amp;B, const Teuchos::ETransp trans=Teuchos::NO_TRANS, const std::string prefix=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeVectorNorms</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>acb2b90901fcf3b898e045efa50b3d5bf</anchor>
      <arglist>(const RCP&lt; Matrix &gt; X, const RCP&lt; Vector &gt; B, std::string prefix=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printLine</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a68278fd71c063772d61aaa863871ddcf</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Amesos2_Util_def.hpp</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>Amesos2__Util__def_8hpp</filename>
    <member kind="function">
      <type>RCP&lt; Tpetra::Map&lt; LO, GO, Node &gt; &gt;</type>
      <name>epetra_map_to_tpetra_map</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>aed570a2f70ea3395c64a75f1c27bcdca</anchor>
      <arglist>(const Epetra_BlockMap &amp;map)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; Epetra_Map &gt;</type>
      <name>tpetra_map_to_epetra_map</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a93e07b55ec4d21a1f9e72e2c2ef30d7b</anchor>
      <arglist>(const Tpetra::Map&lt; LO, GO, Node &gt; &amp;map)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>to_teuchos_comm</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a10fd40490abfa5ee7f80d0b05f7805b3</anchor>
      <arglist>(Teuchos::RCP&lt; const Epetra_Comm &gt; c)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Epetra_Comm &gt;</type>
      <name>to_epetra_comm</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a87a8a9d612ba535b85857cb9969f8344</anchor>
      <arglist>(Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt; c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>transpose</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>aa99673c67895b46c076bf36ea9ef746c</anchor>
      <arglist>(Teuchos::ArrayView&lt; Scalar &gt; vals, Teuchos::ArrayView&lt; GlobalOrdinal &gt; indices, Teuchos::ArrayView&lt; GlobalSizeT &gt; ptr, Teuchos::ArrayView&lt; Scalar &gt; trans_vals, Teuchos::ArrayView&lt; GlobalOrdinal &gt; trans_indices, Teuchos::ArrayView&lt; GlobalSizeT &gt; trans_ptr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scale</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>ad6d975ce99679e335eee7553ced30259</anchor>
      <arglist>(Teuchos::ArrayView&lt; Scalar1 &gt; vals, size_t l, size_t ld, Teuchos::ArrayView&lt; Scalar2 &gt; s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scale</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a8a5d61e1c5036f91680465532e3e14cb</anchor>
      <arglist>(Teuchos::ArrayView&lt; Scalar1 &gt; vals, size_t l, size_t ld, Teuchos::ArrayView&lt; Scalar2 &gt; s, BinaryOp binary_op)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeTrueResidual</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a72a207610994242ae87eddca1660a089</anchor>
      <arglist>(const Teuchos::RCP&lt; Matrix &gt; &amp;A, const Teuchos::RCP&lt; Vector &gt; &amp;X, const Teuchos::RCP&lt; Vector &gt; &amp;B, const Teuchos::ETransp trans=Teuchos::NO_TRANS, const std::string prefix=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeVectorNorms</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>ab53374443490f66b316ce14502615c5c</anchor>
      <arglist>(const Teuchos::RCP&lt; Matrix &gt; X, const Teuchos::RCP&lt; Vector &gt; B, std::string prefix=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printLine</name>
      <anchorfile>namespaceAmesos_1_1Util.html</anchorfile>
      <anchor>a68278fd71c063772d61aaa863871ddcf</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>amesos2_adapters</name>
    <title>Amesos2 Linear Algebra Object Adapters</title>
    <filename>group__amesos2__adapters.html</filename>
    <subgroup>amesos2_matrix_adapters</subgroup>
    <subgroup>amesos2_multivec_adapters</subgroup>
  </compound>
  <compound kind="group">
    <name>amesos2_matrix_adapters</name>
    <title>Amesos2 Matrix Adapters</title>
    <filename>group__amesos2__matrix__adapters.html</filename>
    <class kind="class">Amesos::ConcreteMatrixAdapter&lt; Epetra_CrsMatrix &gt;</class>
    <class kind="class">Amesos::AbstractConcreteMatrixAdapter&lt; Epetra_RowMatrix, DerivedMat &gt;</class>
    <class kind="class">Amesos::MatrixAdapter</class>
    <class kind="class">Amesos::ConcreteMatrixAdapter&lt; Tpetra::CrsMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps &gt; &gt;</class>
    <class kind="class">Amesos::AbstractConcreteMatrixAdapter&lt; Tpetra::RowMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt;, DerivedMat &gt;</class>
  </compound>
  <compound kind="group">
    <name>amesos2_multivec_adapters</name>
    <title>Amesos2 MultiVector Adapters</title>
    <filename>group__amesos2__multivec__adapters.html</filename>
    <class kind="class">Amesos::MultiVecAdapter&lt; Epetra_MultiVector &gt;</class>
    <class kind="struct">Amesos::MultiVecAdapter</class>
    <class kind="class">Amesos::MultiVecAdapter&lt; Tpetra::MultiVector&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt; &gt;</class>
  </compound>
  <compound kind="group">
    <name>amesos2_solvers</name>
    <title>Amesos2 Solvers</title>
    <filename>group__amesos2__solvers.html</filename>
    <subgroup>amesos2_solver_framework</subgroup>
    <subgroup>amesos2_solver_interfaces</subgroup>
    <subgroup>amesos2_solver_parameters</subgroup>
  </compound>
  <compound kind="group">
    <name>amesos2_solver_framework</name>
    <title>Amesos2 Solver Framework</title>
    <filename>group__amesos2__solver__framework.html</filename>
    <class kind="class">Amesos::Solver</class>
    <class kind="class">Amesos::SolverCore</class>
    <file>Amesos2_Factory_decl.hpp</file>
  </compound>
  <compound kind="group">
    <name>amesos2_solver_interfaces</name>
    <title>Amesos2 Solver Interfaces</title>
    <filename>group__amesos2__solver__interfaces.html</filename>
    <class kind="class">Amesos::Superlu</class>
    <class kind="class">Amesos::Superlumt</class>
  </compound>
  <compound kind="group">
    <name>amesos2_solver_parameters</name>
    <title>Supported Solver Parameters</title>
    <filename>group__amesos2__solver__parameters.html</filename>
    <docanchor file="group__amesos2__solver__parameters">superlu_parameters</docanchor>
    <docanchor file="group__amesos2__solver__parameters">amesos2_solver_parameters</docanchor>
    <docanchor file="group__amesos2__solver__parameters">amesos2_parameters</docanchor>
    <docanchor file="group__amesos2__solver__parameters">superlu_mt_parameters</docanchor>
  </compound>
  <compound kind="group">
    <name>amesos2_util</name>
    <title>Amesos2 Utilities</title>
    <filename>group__amesos2__util.html</filename>
    <class kind="struct">Amesos::Util::get_cxs_helper</class>
    <class kind="struct">Amesos::Util::get_ccs_helper</class>
    <class kind="struct">Amesos::Util::get_crs_helper</class>
  </compound>
  <compound kind="class">
    <name>Amesos::AbstractConcreteMatrixAdapter</name>
    <filename>classAmesos_1_1AbstractConcreteMatrixAdapter.html</filename>
    <templarg></templarg>
    <templarg></templarg>
  </compound>
  <compound kind="class">
    <name>Amesos::ConcreteMatrixAdapter&lt; Epetra_CrsMatrix &gt;</name>
    <filename>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</filename>
    <base>AbstractConcreteMatrixAdapter&lt; Epetra_RowMatrix, Epetra_CrsMatrix &gt;</base>
    <member kind="typedef">
      <type>Epetra_CrsMatrix</type>
      <name>matrix_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a685bed98cd31188329b4b5afc546d3e0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::scalar_t</type>
      <name>scalar_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a88e2fa80d21933e7d9af74c6525eb6c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::local_ordinal_t</type>
      <name>local_ordinal_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a59e6daf7ebf45ebaf22fe42e2d5a93b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::global_ordinal_t</type>
      <name>global_ordinal_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a4d7743178e60baab81cf99416275e79b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::node_t</type>
      <name>node_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a3071c59f6ddb968522a50b88f486f721</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::global_size_t</type>
      <name>global_size_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a093d45e64f4f901b63548aa5dbe89b6c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ConcreteMatrixAdapter&lt; matrix_t &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>ac10bf2cfa03713732f31daff8cfcd94c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ConcreteMatrixAdapter</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a44008b452ad8127a62e458eba897074f</anchor>
      <arglist>(RCP&lt; matrix_t &gt; m)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const MatrixAdapter&lt; matrix_t &gt; &gt;</type>
      <name>get_impl</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a7a808e4fda6f69586a62ae3d1a828faa</anchor>
      <arglist>(const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; map) const </arglist>
    </member>
    <member kind="typedef" protection="private">
      <type>AbstractConcreteMatrixAdapter&lt; Epetra_RowMatrix, Epetra_CrsMatrix &gt;</type>
      <name>super_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>ad146411a9e1d57dd5b8c6fab2e4ca32b</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>MatrixAdapter&lt; Epetra_RowMatrix &gt;</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Epetra__CrsMatrix_01_4.html</anchorfile>
      <anchor>a2774c14861cd340b9aee7c56e3a57a15</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::MultiVecAdapter&lt; Epetra_MultiVector &gt;</name>
    <filename>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</filename>
    <member kind="typedef">
      <type>double</type>
      <name>scalar_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a99abff6d59d8b6171c507f57dea80b71</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>local_ordinal_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>aee1c62eeae6a5f9b934e92a23e304c53</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>global_ordinal_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a1e590c8b072d2638e72461083042d51c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>size_t</type>
      <name>global_size_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a6f445c00fedf9c310784106b8ea97479</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Tpetra::DefaultPlatform::DefaultPlatformType::NodeType</type>
      <name>node_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>ac837afbdb3cebb9ccf2b4b4cadea14c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Epetra_MultiVector</type>
      <name>multivec_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a48b1f3e648f59fdfe3b4fe062b1aae78</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MultiVecAdapter</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a316c884859a8ecf2643cf82beeb6f392</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MultiVecAdapter</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>acb3e78f5560bf3c51a7e9c140f22656b</anchor>
      <arglist>(const MultiVecAdapter&lt; multivec_type &gt; &amp;adapter)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MultiVecAdapter</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a5ba697693b39ea867fae048e01f2b666</anchor>
      <arglist>(const Teuchos::RCP&lt; multivec_type &gt; &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; multivec_type &gt;</type>
      <name>getAdaptee</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>aca38eee6991e56abc35a16754a12699e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MultiVecAdapter&lt; multivec_type &gt; &amp;</type>
      <name>scale</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>aa55606898e80fbba2d0b170f49e28f7e</anchor>
      <arglist>(const scalar_type alpha)</arglist>
    </member>
    <member kind="function">
      <type>MultiVecAdapter&lt; multivec_type &gt; &amp;</type>
      <name>update</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>abbc37d7ce559cd86d2641cb0ec7d8b7a</anchor>
      <arglist>(const scalar_type beta, const MultiVecAdapter&lt; multivec_type &gt; &amp;B, const scalar_type alpha)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isLocal</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>acc1312788158ef21e72975c9dc393e85</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>ab39e02be6eb072cb3e01e6961b663e74</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalLength</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a9d36a7e62ab093171f5f48f024d87478</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumVectors</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>aabd89d3d4805d384c365110cd64e16e2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>global_size_type</type>
      <name>getGlobalLength</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>afe85534858698846e34fa7c3c75f4ea4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getGlobalNumVectors</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a614a00c332793799f14c3b01064569bc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getStride</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a3e2d1b2b944598ec45bbbc9303699a2b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isConstantStride</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a4e855f61b54b51264142aee2577c873a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Tpetra::Vector&lt; scalar_type, local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>getVector</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a0605dc5de54885e2f0bdf38c8039c2fb</anchor>
      <arglist>(size_t j) const </arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Tpetra::Vector&lt; scalar_type, local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>getVectorNonConst</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a74ed7666f236b92b5428b1da46d59a68</anchor>
      <arglist>(size_t j)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get1dCopy</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>ae4401c9411b2ad7d89fd30e6319f4523</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_type &gt; &amp;A, size_t lda, bool global_copy) const </arglist>
    </member>
    <member kind="function">
      <type>Teuchos::ArrayRCP&lt; scalar_type &gt;</type>
      <name>get1dViewNonConst</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>abbb87b03e54e45c587c00a43bd4d7d7c</anchor>
      <arglist>(bool local=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get2dCopy</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>ad2c59a3308c90ce07a932e823c987bc1</anchor>
      <arglist>(Teuchos::ArrayView&lt; const Teuchos::ArrayView&lt; scalar_type &gt; &gt; A)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::ArrayRCP&lt; Teuchos::ArrayRCP&lt; scalar_type &gt; &gt;</type>
      <name>get2dViewNonConst</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>aedabfc314c389bc354a20e9e3245a95c</anchor>
      <arglist>(bool local=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>globalize</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a6094863621a9be6b16912235498651f7</anchor>
      <arglist>(int root=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>globalize</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a3854190e2454115fb12b5d1a51e1b132</anchor>
      <arglist>(const Teuchos::ArrayView&lt; Value_t &gt; &amp;newVals, int root=0)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>ab407d2d5cbaf71ef8197227279834d78</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a28bba365e0d9d6ed1cdac7158070564b</anchor>
      <arglist>(Teuchos::FancyOStream &amp;os, const Teuchos::EVerbosityLevel verbLevel) const </arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const char *</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>ac3636cf08be452a47b8708a94916359f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>localize</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a6c4570daffa3c98da032de0e18401267</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Teuchos::RCP&lt; const Tpetra::Map&lt; local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>getMap</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a711960b7d88e2bca4fa7e1a06f3fd335</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; multivec_type &gt;</type>
      <name>mv_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a584f8ea2b9d8be22849836fa12c6c8e5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; multivec_type &gt;</type>
      <name>l_mv_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a6e54815ead43181d43c22cf51a33fe75</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; multivec_type &gt;</type>
      <name>l_l_mv_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a8228cd2f44451a855dd10ba5e4a827c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Epetra_BlockMap &gt;</type>
      <name>o_map_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a89e11abb2f07a140c1a557c4cf723b7f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Epetra_BlockMap &gt;</type>
      <name>l_map_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>ab03c72fc0f3f52b6c5abcea4400e231e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Epetra_Import &gt;</type>
      <name>importer_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>afb5b0fdc6138156ae2068d263ac94b00</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Epetra_Export &gt;</type>
      <name>exporter_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Epetra__MultiVector_01_4.html</anchorfile>
      <anchor>a0ee52974ca6d994263186d3a7aba3512</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::AbstractConcreteMatrixAdapter&lt; Epetra_RowMatrix, DerivedMat &gt;</name>
    <filename>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</filename>
    <templarg></templarg>
    <base>MatrixAdapter&lt; DerivedMat &gt;</base>
    <member kind="typedef">
      <type>MatrixTraits&lt; Epetra_RowMatrix &gt;::scalar_t</type>
      <name>scalar_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>affddce1920a772c3b8229a0a4ced0d8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; Epetra_RowMatrix &gt;::local_ordinal_t</type>
      <name>local_ordinal_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>ae009b0b3f9913d96ab47a173970d7090</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; Epetra_RowMatrix &gt;::global_ordinal_t</type>
      <name>global_ordinal_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a223c1fd684c762a5f2c2f8cba5d34b65</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; Epetra_RowMatrix &gt;::node_t</type>
      <name>node_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>adbe7a699eeebb39e3cdd20f38d03cf65</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DerivedMat</type>
      <name>matrix_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a525d2cedc9f348b1f4640e50882a35b3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::global_size_t</type>
      <name>global_size_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a584a40f6c7a3b3ddd6844531b6490086</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>AbstractConcreteMatrixAdapter&lt; matrix_t, DerivedMat &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>ad309903f8ef1e93182011d90de5fd020</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Util::no_special_impl</type>
      <name>get_crs_spec</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>ad477ed334bd099f71d17d9b8eb3183dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Util::no_special_impl</type>
      <name>get_ccs_spec</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a298f743d18de12fb71bd84f2363fc80c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; Epetra_RowMatrix &gt;::major_access</type>
      <name>major_access</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>ab4830bbfa6b189298808e35e69daf544</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ConcreteMatrixAdapter&lt; DerivedMat &gt;</type>
      <name>adapter_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a3bfa635f15a7ee77ba78dc64f5822e9c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>AbstractConcreteMatrixAdapter</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>ab178b95ed4f9e2385cdd39bc17cc1d66</anchor>
      <arglist>(RCP&lt; matrix_t &gt; m)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getGlobalRowCopy_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a99e7efa2197a9433662d9355248f1ba0</anchor>
      <arglist>(global_ordinal_t row, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getGlobalColCopy_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a0c28244bac0481602c1b3433e2285f2f</anchor>
      <arglist>(global_ordinal_t col, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const </arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a5265fa69aec138d473f79147e86c6193</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>ab47c9560145f47f2026394c11d3f4545</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMaxRowNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a58b31425bc1b98b230a4b56a9ee79e17</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMaxColNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a7168892e0791cf587ae443f40c7adb69</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getGlobalRowNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a6e97df991c83a3e7dabbdfb20aeaae62</anchor>
      <arglist>(global_ordinal_t row) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalRowNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a920f297d72689dc7bdd6049f9fa62f01</anchor>
      <arglist>(local_ordinal_t row) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getGlobalColNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>adbf1390ae0f9fcdc0b31afcab57fa80f</anchor>
      <arglist>(global_ordinal_t col) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalColNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a41b6838a4f7d8ca4349f9e6a00be88f6</anchor>
      <arglist>(local_ordinal_t col) const </arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>getRowMap_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>aa7a85aad361aaf0dbecd679cee8b7fe2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>getColMap_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a075ff46c1cca66abda570323c7d10d23</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a719631507055da6bbc6da9e40493b5a6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isLocallyIndexed_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>abdde9fd45e5d9ef5479902020fd12b36</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isGloballyIndexed_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a6fd983c3e9429dfccdb210c8bffb808e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const super_t &gt;</type>
      <name>get_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>ac1d3a20e1706792602544e498d455fbf</anchor>
      <arglist>(const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; map) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a7da98dde7decdcd4006c4203193a1edb</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; rowmap, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a9b9abaec0efdf184ac502bcf0ef08b3c</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a54a1fb9c3f9203c6a7ec4ba9c87dae25</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; colmap, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad85223e6476499e4ab135661bf8339aa</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aae644193c73b59c573720132be955431</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNumRows</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ac729469ca361ac261590f143624fa6b4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNumCols</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a2c1f4eed533cf858db94067149c70c0b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aa144b89709fd20ed4627d5966997132e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumRows</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aa72c5a28011f01c5f09f1be50df5de2e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumCols</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ada217a4504703e97a0ed2e165ec9b766</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ac1fcfa87de3617bc6fee93f884ca6bf7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t &gt; &gt;</type>
      <name>getRowMap</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a079a384379ea6ba42f82f254ec322ab2</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t &gt; &gt;</type>
      <name>getColMap</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad095ea1b4a87c53258bc95344540d369</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aac6aec881feab1fa3d30c4c3a7d3a2f1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a0e6d8d17887b32364e946403a2617365</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getGlobalRowCopy</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a82a04efccdb9fc4c677f6971b0f2a9ff</anchor>
      <arglist>(global_ordinal_t row, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getGlobalColCopy</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a1930fc8e940fa82ae8933c8fbd61c58d</anchor>
      <arglist>(global_ordinal_t col, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getMaxRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aec3505d353a2b1bfac631222aa6e582d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getMaxColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>af4fca7c867756723a39d7251eb5e214e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getGlobalRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a8d5fb3fc6c2db1015f75399f880d8e1e</anchor>
      <arglist>(global_ordinal_t row) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getLocalRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a6286f7e449913ec698831cbe9f855ed0</anchor>
      <arglist>(local_ordinal_t row) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getGlobalColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a139320c33a3e6954e1263ac083c56fe9</anchor>
      <arglist>(global_ordinal_t col) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getLocalColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a59b1d89868a0044cdb5d6442c505ba21</anchor>
      <arglist>(local_ordinal_t col) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>isLocallyIndexed</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a37bec34548dee2f1bfb9ab42f5c872a6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>isGloballyIndexed</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a760657acc1d3e7abfcb35b689ed7c3a5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>RCP&lt; const type &gt;</type>
      <name>get</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ab2e42b619d1ecc2796bce408615bd15d</anchor>
      <arglist>(const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; map) const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const RCP&lt; const DerivedMat &gt;</type>
      <name>mat_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aaa7c67b31929c0252ceb73e3a5fc8abf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>row_map_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a91680cc808fdc40b64a01ab7ca6ca554</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>col_map_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a06617c51632049ccc3dee8fe7cdcf4f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>comm_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a898d99aef23d15a31a92fc54ea9a5116</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef" protection="private">
      <type>MatrixAdapter&lt; DerivedMat &gt;</type>
      <name>super_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>aebcc1559a498de5a058aeb592d26936c</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>MatrixAdapter&lt; DerivedMat &gt;</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Epetra__RowMatrix_00_01DerivedMat_01_4.html</anchorfile>
      <anchor>a0a0a3629b22324189309040873178be5</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>Util::get_cxs_helper</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad17a6004ba9322d05cb84067e15a3016</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Amesos::FunctionMap</name>
    <filename>structAmesos_1_1FunctionMap.html</filename>
    <templarg>ConcreteSolver</templarg>
    <templarg></templarg>
  </compound>
  <compound kind="class">
    <name>Amesos::MatrixAdapter</name>
    <filename>classAmesos_1_1MatrixAdapter.html</filename>
    <templarg>Matrix</templarg>
    <member kind="typedef">
      <type>MatrixTraits&lt; Matrix &gt;::scalar_t</type>
      <name>scalar_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a37e717815af215fdd76eb3c8b9af2888</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; Matrix &gt;::local_ordinal_t</type>
      <name>local_ordinal_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aaec907ef5173c4fa3fd096c2f1563ec9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; Matrix &gt;::global_ordinal_t</type>
      <name>global_ordinal_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a2e387c76fc0fc760deb7ac4030bd820a</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; Matrix &gt;::node_t</type>
      <name>node_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aac35e7fe67eeddbc7fd9ce5237f19746</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Tpetra::global_size_t</type>
      <name>global_size_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a998e3af52656ac6c583043a3be7f9b65</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix</type>
      <name>matrix_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a08c8937cfff1d1c1420d99063c86c7e9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; Matrix &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a4ae7efbcc781289cb68af408e5f19160</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ConcreteMatrixAdapter&lt; Matrix &gt;</type>
      <name>adapter_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a3bfa635f15a7ee77ba78dc64f5822e9c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MatrixAdapter</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ab3d7d119a1beeccfe22490a4791e4d0a</anchor>
      <arglist>(RCP&lt; Matrix &gt; m)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a7da98dde7decdcd4006c4203193a1edb</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; rowmap, EStorage_Ordering ordering=Util::Arbitrary) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a9b9abaec0efdf184ac502bcf0ef08b3c</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Util::Arbitrary) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a54a1fb9c3f9203c6a7ec4ba9c87dae25</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; colmap, EStorage_Ordering ordering=Util::Arbitrary) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad85223e6476499e4ab135661bf8339aa</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Util::Arbitrary) const </arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aae644193c73b59c573720132be955431</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNumRows</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ac729469ca361ac261590f143624fa6b4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNumCols</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a2c1f4eed533cf858db94067149c70c0b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aa144b89709fd20ed4627d5966997132e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumRows</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aa72c5a28011f01c5f09f1be50df5de2e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumCols</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ada217a4504703e97a0ed2e165ec9b766</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ac1fcfa87de3617bc6fee93f884ca6bf7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t &gt; &gt;</type>
      <name>getRowMap</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a079a384379ea6ba42f82f254ec322ab2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t &gt; &gt;</type>
      <name>getColMap</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad095ea1b4a87c53258bc95344540d369</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aac6aec881feab1fa3d30c4c3a7d3a2f1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a0e6d8d17887b32364e946403a2617365</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getGlobalRowCopy</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a82a04efccdb9fc4c677f6971b0f2a9ff</anchor>
      <arglist>(global_ordinal_t row, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getGlobalColCopy</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a1930fc8e940fa82ae8933c8fbd61c58d</anchor>
      <arglist>(global_ordinal_t col, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getMaxRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aec3505d353a2b1bfac631222aa6e582d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getMaxColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>af4fca7c867756723a39d7251eb5e214e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getGlobalRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a8d5fb3fc6c2db1015f75399f880d8e1e</anchor>
      <arglist>(global_ordinal_t row) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getLocalRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a6286f7e449913ec698831cbe9f855ed0</anchor>
      <arglist>(local_ordinal_t row) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getGlobalColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a139320c33a3e6954e1263ac083c56fe9</anchor>
      <arglist>(global_ordinal_t col) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getLocalColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a59b1d89868a0044cdb5d6442c505ba21</anchor>
      <arglist>(local_ordinal_t col) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>isLocallyIndexed</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a37bec34548dee2f1bfb9ab42f5c872a6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>isGloballyIndexed</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a760657acc1d3e7abfcb35b689ed7c3a5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>RCP&lt; const type &gt;</type>
      <name>get</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ab2e42b619d1ecc2796bce408615bd15d</anchor>
      <arglist>(const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; map) const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const RCP&lt; const Matrix &gt;</type>
      <name>mat_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aaa7c67b31929c0252ceb73e3a5fc8abf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>row_map_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a91680cc808fdc40b64a01ab7ca6ca554</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>col_map_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a06617c51632049ccc3dee8fe7cdcf4f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>comm_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a898d99aef23d15a31a92fc54ea9a5116</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>help_getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a516319c4eaf1fe77a80e8aee76a86ae1</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; rowmap, EStorage_Ordering ordering, has_special_impl hsi) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>help_getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>af78df62595bc621c7b0b3f01061f1846</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; rowmap, EStorage_Ordering ordering, no_special_impl nsi) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>do_getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a9bc0fa6fda56e9b8b9243b768bc294ce</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; rowmap, EStorage_Ordering ordering, row_access ra) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>do_getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a71c4e15fae62d88416aa0a29220b0fb5</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; rowmap, EStorage_Ordering ordering, col_access ca) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>help_getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>af3b897379ed6ab2db8dcaa67ab8804c5</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; colmap, EStorage_Ordering ordering, has_special_impl hsi) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>help_getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a335456d757f8a2d6cebcfbb792cded4a</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; colmap, EStorage_Ordering ordering, no_special_impl nsi) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>do_getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a46b37e798269d77885cded85de493c60</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; colmap, EStorage_Ordering ordering, row_access ra) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>do_getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a517096171045f2c97409f4902aee04d5</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; colmap, EStorage_Ordering ordering, col_access ca) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>const Teuchos::RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>getDistributionMap</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a339df63c80cb78a004417574164d1f0c</anchor>
      <arglist>(Util::EDistribution distribution) const </arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>Util::get_cxs_helper</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad17a6004ba9322d05cb84067e15a3016</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Amesos::MatrixHelper</name>
    <filename>structAmesos_1_1MatrixHelper.html</filename>
    <templarg>ConcreteSolver</templarg>
  </compound>
  <compound kind="struct">
    <name>Amesos::MultiVecAdapter</name>
    <filename>structAmesos_1_1MultiVecAdapter.html</filename>
    <templarg>MultiVecType</templarg>
    <member kind="function">
      <type>Teuchos::RCP&lt; MultiVecAdapter&lt; MV &gt; &gt;</type>
      <name>createMultiVecAdapter</name>
      <anchorfile>structAmesos_1_1MultiVecAdapter.html</anchorfile>
      <anchor>aeb4d8aaf14a2b2b2eb8c0cf8ce67d28f</anchor>
      <arglist>(Teuchos::RCP&lt; MV &gt; mv)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::Solver</name>
    <filename>classAmesos_1_1Solver.html</filename>
    <templarg>Matrix</templarg>
    <templarg>Vector</templarg>
    <member kind="typedef">
      <type>Solver&lt; Matrix, Vector &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a8cdf7660b04fbef2ab7e67254e70bfa5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual type &amp;</type>
      <name>preOrdering</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>acd253d0f641b84ec92920492ec83dcaf</anchor>
      <arglist>(void)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual type &amp;</type>
      <name>symbolicFactorization</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a20c17edc5bb2cc59afd29d4fcd23c119</anchor>
      <arglist>(void)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual type &amp;</type>
      <name>numericFactorization</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a758cfc703a91de4d5d0d7599de0c429e</anchor>
      <arglist>(void)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ad5055ac72f596931984160a8f6c83d3d</anchor>
      <arglist>(void)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a1ef4c6ede05860d16811094bc64844b3</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; const Vector &gt; B) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual type &amp;</type>
      <name>setParameters</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a82af5cf1d5624c19d94ca8ed8ebc5a6e</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a8a1923bc53a3e78585d9878075eaa490</anchor>
      <arglist>(void) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setA</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>adfc219139987f8a17133186ebb748617</anchor>
      <arglist>(const Teuchos::RCP&lt; Matrix &gt; a)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>matrixShapeOK</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a8d031d5530cb0387cb5b76a621cb226b</anchor>
      <arglist>(void)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setX</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a3fbf70e39fabe035c4564f5cc5caf448</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; x)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const Teuchos::RCP&lt; Vector &gt;</type>
      <name>getX</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a066fc8d9963f874c384420fa09c2f3e5</anchor>
      <arglist>(void)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setB</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a0b805f4759a96a0d8e39a463cce04da0</anchor>
      <arglist>(const Teuchos::RCP&lt; const Vector &gt; b)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const Teuchos::RCP&lt; const Vector &gt;</type>
      <name>getB</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a6e42153b106323f98c31e5684101c4c4</anchor>
      <arglist>(void)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ae47f7c19ded1827114b3791ccbd650d0</anchor>
      <arglist>(void) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getNumSymbolicFact</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a21b281e1147ae78b6f171c5ec1633401</anchor>
      <arglist>(void) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getNumNumericFact</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ac8c9615752ca229033625a429bcfbded</anchor>
      <arglist>(void) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getNumSolve</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a78a3c4c88d792f702179a27407bd87a3</anchor>
      <arglist>(void) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::string</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>afdfaf68e6c49c92c6451936503c6e294</anchor>
      <arglist>(void) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a95460337927b17ea14852803ec9d5a80</anchor>
      <arglist>(void) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>acd201473e99e7210946067f6f9e62ed2</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>printTiming</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a49d43c33ef101053bc204a40e398ce81</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>getTiming</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a661eec947d0e688ce14afbdd95ef55df</anchor>
      <arglist>(Teuchos::ParameterList &amp;timingParameterList) const =0</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a415501e885f7e753298fa2cf61793932</anchor>
      <arglist>(Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a479b23341a0b171d241cffc87a92bcf0</anchor>
      <arglist>(Teuchos::RCP&lt; Matrix &gt; A, Teuchos::RCP&lt; Vector &gt; X, Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ad8aa0892d3a641a41d18884185d31ab6</anchor>
      <arglist>(const char *solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a53132a9d360914a443142b7101294a07</anchor>
      <arglist>(const char *solverName, const Teuchos::RCP&lt; Matrix &gt; A, const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>abb8acb7cceac1ce640a41dd05cc6024b</anchor>
      <arglist>(const std::string solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>aea67ca5c343fc0f248abb8bd59423098</anchor>
      <arglist>(const std::string solverName, const Teuchos::RCP&lt; Matrix &gt; A, const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ab239b521c543b7cbed4f8823bb84c23f</anchor>
      <arglist>(const std::string solverName, Matrix *A)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>afed8390a758da6367f4f28c1602df807</anchor>
      <arglist>(const std::string solverName, const Teuchos::RCP&lt; Matrix &gt; A)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::SolverCore</name>
    <filename>classAmesos_1_1SolverCore.html</filename>
    <templarg>ConcreteSolver</templarg>
    <templarg>Matrix</templarg>
    <templarg>Vector</templarg>
    <base>Solver&lt; Matrix, Vector &gt;</base>
    <member kind="typedef">
      <type>SolverCore&lt; ConcreteSolver, Matrix, Vector &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e3660d3adce807015ae6b1490fedb16</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Solver&lt; Matrix, Vector &gt;</type>
      <name>super_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae0a58e56723c5dd02253c8e7ea2b05fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ConcreteSolver&lt; Matrix, Vector &gt;</type>
      <name>solver_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a6dacca606d92777840246d455c94f335</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix</type>
      <name>matrix_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>af1cdac3715d3b1d826483444c5a42743</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Vector</type>
      <name>vector_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a68e5678b5631e357bf2102c91e2b9bad</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; matrix_type &gt;::scalar_t</type>
      <name>scalar_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a2813db55d915edf1cea4377f8408c40b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; matrix_type &gt;::local_ordinal_t</type>
      <name>local_ordinal_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a06f1f172d073fc8530219aa5ff3ccca3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; matrix_type &gt;::global_ordinal_t</type>
      <name>global_ordinal_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a6595e103e535d53a15125a4340e46bb0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; matrix_type &gt;::global_size_t</type>
      <name>global_size_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a8e96cb999e58dbbd0527bced160f69cb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>matrixShapeOK</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a99a5d733412d7b9f8e7ff692e6c27c2a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setA</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a72d655870f183deb536d126fe161a3a6</anchor>
      <arglist>(const Teuchos::RCP&lt; Matrix &gt; a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac36fe794be11b6d03547e71dd8f89ac4</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; x)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; Vector &gt;</type>
      <name>getX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeee7773536622735fe0f2b3242521e05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>afc887943fe97b8bbc16809606692ac0f</anchor>
      <arglist>(const Teuchos::RCP&lt; const Vector &gt; b)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Vector &gt;</type>
      <name>getB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeb10617c7d0feab6bfab5c8d4240ac9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a28df22e8d7ea6c53d4f38fe08c253b5b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9ffe894f7d3973a9081b43eb2cd3b974</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a35f2e143409f7fcbf6e1b3c610f4a676</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a684fc9ecf86f1dc99462faa7aece791b</anchor>
      <arglist>(Teuchos::ParameterList &amp;timingParameterList) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aade60205a3b6ece6bdc76f5d1555184e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SolverCore</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a67e4c02187489d7c57c8f0621a003863</anchor>
      <arglist>(Teuchos::RCP&lt; Matrix &gt; A, Teuchos::RCP&lt; Vector &gt; X, Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SolverCore</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a56fdd2b797e376c8572b8c51322812da</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SolverCore</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a463a303468c55cf33364f5b759e2d60b</anchor>
      <arglist>(const solver_type &amp;rhs)</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>operator=</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1535500c05cddba86a0c915be4bc7c5a</anchor>
      <arglist>(const solver_type *rhs)</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>preOrdering</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ace48a1882df1a866f7fba00e321d1abe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>symbolicFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ab0602bd80f5ba358adf710423c013186</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>numericFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a00d8910072b4ff26abbb71e81ed30ef9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a89ed4e401b4836bee75de516b6de12cd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a7ab43f4fd19cb4f3e7b9c1068d78f381</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; const Vector &gt; B) const </arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>setParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a2f812b0093f7491cf86574eb5f5d1117</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac42d1a2b4c037065b2637c06b23b06bd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a53d48aa613aca50ed957272aeeeb9485</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>getNonconstParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bf21b207f1e9eb860fe326e3fb6dcbc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>unsetParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4a62e2cbe4d173534302af0f93088a05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a6b805f6cab10b90c0843f908e00e0344</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSymbolicFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a559b355352ff76b91d7ce2cc94b6242e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumNumericFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a32637ffd666f7af4726811eb7fdaa47a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSolve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a942a1be951a2bb725aef7039219e6381</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>preOrderingDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9f679ee04f9db29ac263dead626a6c76</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>symbolicFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bb9fc2c4313ba09f8b74b3c190f1704</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>numericFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a549a2117cc7b78842f55627be58cffe2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; MatrixAdapter&lt; Matrix &gt; &gt;</type>
      <name>matrixA_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a73862ad4870aaad942cf1d85910e035e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; Vector &gt;</type>
      <name>multiVecX_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a57e4fe1d6563e004f933a390d89f5d6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; const Vector &gt;</type>
      <name>multiVecB_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a8948ea996ec300e9423b3c6977197b0b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumRows_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4f98ce6a877b2f6ece913c05adef6fc4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumCols_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9c9458f37c6bfe6aa23be99213265ebf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumNonZeros_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e824badc831fc6cd4627ea834b3782d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Status</type>
      <name>status_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>af688193f21876bb4a209738d2bfeae3e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Control</type>
      <name>control_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae9cc6813b5540966424eaab3cd9b0e15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Timers</type>
      <name>timers_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aa6a15406034089ec343445003db0a276</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>refreshA</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a73d94a5327526b05c8b6032d7dddfa44</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>EPhase</type>
      <name>last_phase_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a5083cf9c3b368b22f6ff42a996a9d93a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::Superlu</name>
    <filename>classAmesos_1_1Superlu.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <base>Amesos::SolverCore</base>
    <member kind="typedef">
      <type>Superlu&lt; Matrix, Vector &gt;</type>
      <name>solver_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a7fc9203b4723fe06ca2648ff5e3f5ea5</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix</type>
      <name>matrix_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a03e5d53c814c4aec43a2a6b49a50d632</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Vector</type>
      <name>vector_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a0ac170c6769977d954ddebc02dcc3591</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::scalar_t</type>
      <name>scalar_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a3c783b4c6b2fa475195be19cb7d7f089</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TypeMap&lt; Amesos::Superlu, scalar_type &gt;::magnitude_type</type>
      <name>magnitude_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a6ee3035b363c06389c316e0f1142c205</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::local_ordinal_t</type>
      <name>local_ordinal_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>ac29c9c5d05b92e89b707913488c653dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::global_ordinal_t</type>
      <name>global_ordinal_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>ae29b80b4387cbc7b635efc6497b9d234</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; matrix_type &gt;::global_size_t</type>
      <name>global_size_type</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>ad64bf151d92cc5867d2ad550016b5a26</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SolverCore&lt; ConcreteSolver, Matrix, Vector &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e3660d3adce807015ae6b1490fedb16</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Solver&lt; Matrix, Vector &gt;</type>
      <name>super_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae0a58e56723c5dd02253c8e7ea2b05fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>matrixShapeOK</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a99a5d733412d7b9f8e7ff692e6c27c2a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setA</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a72d655870f183deb536d126fe161a3a6</anchor>
      <arglist>(const Teuchos::RCP&lt; Matrix &gt; a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac36fe794be11b6d03547e71dd8f89ac4</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; x)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; Vector &gt;</type>
      <name>getX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeee7773536622735fe0f2b3242521e05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>afc887943fe97b8bbc16809606692ac0f</anchor>
      <arglist>(const Teuchos::RCP&lt; const Vector &gt; b)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Vector &gt;</type>
      <name>getB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeb10617c7d0feab6bfab5c8d4240ac9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a28df22e8d7ea6c53d4f38fe08c253b5b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9ffe894f7d3973a9081b43eb2cd3b974</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a35f2e143409f7fcbf6e1b3c610f4a676</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a684fc9ecf86f1dc99462faa7aece791b</anchor>
      <arglist>(Teuchos::ParameterList &amp;timingParameterList) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aade60205a3b6ece6bdc76f5d1555184e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Superlu</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>af77f52720abefe152903af640f3b2d01</anchor>
      <arglist>(Teuchos::RCP&lt; Matrix &gt; A, Teuchos::RCP&lt; Vector &gt; X, Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Superlu</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>ad58105d7f8f4c2a25aff0a97dd4a74b5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>preOrdering</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ace48a1882df1a866f7fba00e321d1abe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>symbolicFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ab0602bd80f5ba358adf710423c013186</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>numericFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a00d8910072b4ff26abbb71e81ed30ef9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a89ed4e401b4836bee75de516b6de12cd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a7ab43f4fd19cb4f3e7b9c1068d78f381</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; const Vector &gt; B) const </arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>setParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a2f812b0093f7491cf86574eb5f5d1117</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac42d1a2b4c037065b2637c06b23b06bd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a53d48aa613aca50ed957272aeeeb9485</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>getNonconstParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bf21b207f1e9eb860fe326e3fb6dcbc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>unsetParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4a62e2cbe4d173534302af0f93088a05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a6b805f6cab10b90c0843f908e00e0344</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSymbolicFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a559b355352ff76b91d7ce2cc94b6242e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumNumericFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a32637ffd666f7af4726811eb7fdaa47a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSolve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a942a1be951a2bb725aef7039219e6381</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const char *</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>add86bcb019c5e5957ac99a62ac9dfa0f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>preOrderingDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9f679ee04f9db29ac263dead626a6c76</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>symbolicFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bb9fc2c4313ba09f8b74b3c190f1704</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>numericFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a549a2117cc7b78842f55627be58cffe2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; MatrixAdapter&lt; Matrix &gt; &gt;</type>
      <name>matrixA_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a73862ad4870aaad942cf1d85910e035e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; Vector &gt;</type>
      <name>multiVecX_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a57e4fe1d6563e004f933a390d89f5d6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; const Vector &gt;</type>
      <name>multiVecB_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a8948ea996ec300e9423b3c6977197b0b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumRows_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4f98ce6a877b2f6ece913c05adef6fc4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumCols_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9c9458f37c6bfe6aa23be99213265ebf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumNonZeros_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e824badc831fc6cd4627ea834b3782d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Status</type>
      <name>status_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>af688193f21876bb4a209738d2bfeae3e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Control</type>
      <name>control_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae9cc6813b5540966424eaab3cd9b0e15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Timers</type>
      <name>timers_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aa6a15406034089ec343445003db0a276</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>preOrdering_impl</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a216162041a294344d6f63e6b9143c93b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>symbolicFactorization_impl</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>ae9de591b55f71ee19ae74af52020384d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>numericFactorization_impl</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a145e984feb3a8dcaed65abf2bb021f6d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>solve_impl</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a1b299a6e81b41b54969392f773e0a9e8</anchor>
      <arglist>(const Teuchos::Ptr&lt; MultiVecAdapter&lt; Vector &gt; &gt; X, const Teuchos::Ptr&lt; MultiVecAdapter&lt; Vector &gt; &gt; B) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>matrixShapeOK_impl</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a611ac3d50741005c3c737e24b0107944</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>setParameters_impl</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a3958e48e255671b7f2a717538d1ddfe5</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters_impl</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a18651c0b1cbc0de69074de8545132547</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>struct Amesos::Superlu::SLUData</type>
      <name>data_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a7f69e41669b5b6124e721ec3b33acc32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; typename TypeMap&lt; Amesos::Superlu, scalar_type &gt;::type &gt;</type>
      <name>nzvals_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>abf2abc87874f9d312a9b7fb184e10857</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; int &gt;</type>
      <name>rowind_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a344e8a3fb495f6b01b096fa1e9e3189c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; int &gt;</type>
      <name>colptr_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>af6954a800f57d694d9eed245ce7b3cfc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; typename TypeMap&lt; Amesos::Superlu, scalar_type &gt;::type &gt;</type>
      <name>xvals_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>ad5c2edfa805189416a8abe33bed72f70</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>ldx_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a6329763e564cc4f037a141fe13b82884</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; typename TypeMap&lt; Amesos::Superlu, scalar_type &gt;::type &gt;</type>
      <name>bvals_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>af9ed7ae4329d7e28f23d76e806fefc65</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>ldb_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>ae70d0f73b6cb97f2bb7a73676d19dfb7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>factorizationDone_</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a6c046eab00755a4fc5853aedb153a338</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>SolverCore&lt; Amesos::Superlu, Matrix, Vector &gt;</name>
      <anchorfile>classAmesos_1_1Superlu.html</anchorfile>
      <anchor>a15859bd317c4ce5564416b53b3d05b7d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Amesos::FunctionMap&lt; Superlu, Scalar &gt;</name>
    <filename>structAmesos_1_1FunctionMap_3_01Superlu_00_01Scalar_01_4.html</filename>
    <templarg></templarg>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gssvx</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlu_00_01Scalar_01_4.html</anchorfile>
      <anchor>ada2db36117b2bcf6eda88675acdac374</anchor>
      <arglist>(SLU::superlu_options_t *options, SLU::SuperMatrix *A, int *perm_c, int *perm_r, int *etree, char *equed, typename TypeMap&lt; Superlu, Scalar &gt;::magnitude_type *R, typename TypeMap&lt; Superlu, Scalar &gt;::magnitude_type *C, SLU::SuperMatrix *L, SLU::SuperMatrix *U, void *work, int lwork, SLU::SuperMatrix *B, SLU::SuperMatrix *X, typename TypeMap&lt; Superlu, Scalar &gt;::magnitude_type *recip_pivot_growth, typename TypeMap&lt; Superlu, Scalar &gt;::magnitude_type *rcond, typename TypeMap&lt; Superlu, Scalar &gt;::magnitude_type *ferr, typename TypeMap&lt; Superlu, Scalar &gt;::magnitude_type *berr, SLU::mem_usage_t *mem_usage, SLU::SuperLUStat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gstrf</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlu_00_01Scalar_01_4.html</anchorfile>
      <anchor>a39aa5ee8777be9fc466e0209f9c2bfce</anchor>
      <arglist>(SLU::superlu_options_t *options, SLU::SuperMatrix *A, int relax, int panel_size, int *etree, void *work, int lwork, int *perm_c, int *perm_r, SLU::SuperMatrix *L, SLU::SuperMatrix *U, SLU::SuperLUStat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_CompCol_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlu_00_01Scalar_01_4.html</anchorfile>
      <anchor>a149af6127d3b9ab27473c735298ccc41</anchor>
      <arglist>(SLU::SuperMatrix *A, int numrows, int numcols, int nnz, typename TypeMap&lt; Superlu, Scalar &gt;::type *nzval, int *rowind, int *colptr, SLU::Stype_t storage_t, SLU::Dtype_t data_t, SLU::Mtype_t mat_t)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_CompRow_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlu_00_01Scalar_01_4.html</anchorfile>
      <anchor>a02c5822c013d312d91ec34403282ace5</anchor>
      <arglist>(SLU::SuperMatrix *A, int numrows, int numcols, int nnz, typename TypeMap&lt; Superlu, Scalar &gt;::type *nzval, int *rowind, int *colptr, SLU::Stype_t storage_t, SLU::Dtype_t data_t, SLU::Mtype_t mat_t)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_Dense_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlu_00_01Scalar_01_4.html</anchorfile>
      <anchor>abe4ca35c7de2d36a166a677867a8947d</anchor>
      <arglist>(SLU::SuperMatrix *X, int numrows, int numcols, typename TypeMap&lt; Superlu, Scalar &gt;::type *x, int ldx, SLU::Stype_t storage_t, SLU::Dtype_t data_t, SLU::Mtype_t mat_t)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::Superludist</name>
    <filename>classAmesos_1_1Superludist.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <base>SolverCore&lt; Amesos::Superludist, Matrix, Vector &gt;</base>
    <member kind="typedef">
      <type>Superludist&lt; Matrix, Vector &gt;</type>
      <name>solver_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a7f5aed1a08636edf84dc16e15ab0c2c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix</type>
      <name>matrix_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>aa58959fecaa9d94ef3b92b461f9fb8da</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Vector</type>
      <name>vector_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a7a7d9ad803ffea25096d7086569fb23c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::scalar_t</type>
      <name>scalar_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a23619367c18783cea452050a2803fce5</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::local_ordinal_t</type>
      <name>local_ordinal_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a2b413cd40b4972b1c65658b45ba43bce</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::global_ordinal_t</type>
      <name>global_ordinal_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a00f0b854066259fb29419340206667bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::node_t</type>
      <name>node_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a76874516d072c886b799d178f31f1339</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; matrix_type &gt;::global_size_t</type>
      <name>global_size_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a7d6f804b00253f990930569eb3df4330</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TypeMap&lt; Amesos::Superludist, scalar_type &gt;::type</type>
      <name>slu_dist_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ab1b9cae603a405e96183b9cc18d2037b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>FunctionMap&lt; Amesos::Superludist, scalar_type &gt;</type>
      <name>function_map</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a2bcc44a82bb705e9dd126be6261e8e13</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixHelper&lt; Amesos::Superludist &gt;</type>
      <name>matrix_helper</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a813b421329df18c05edf15f0b06d3f47</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SolverCore&lt; Amesos::Superludist, Matrix, Vector &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e3660d3adce807015ae6b1490fedb16</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Solver&lt; Matrix, Vector &gt;</type>
      <name>super_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae0a58e56723c5dd02253c8e7ea2b05fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>matrixShapeOK</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a99a5d733412d7b9f8e7ff692e6c27c2a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setA</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a72d655870f183deb536d126fe161a3a6</anchor>
      <arglist>(const Teuchos::RCP&lt; Matrix &gt; a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac36fe794be11b6d03547e71dd8f89ac4</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; x)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; Vector &gt;</type>
      <name>getX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeee7773536622735fe0f2b3242521e05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>afc887943fe97b8bbc16809606692ac0f</anchor>
      <arglist>(const Teuchos::RCP&lt; const Vector &gt; b)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Vector &gt;</type>
      <name>getB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeb10617c7d0feab6bfab5c8d4240ac9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a28df22e8d7ea6c53d4f38fe08c253b5b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9ffe894f7d3973a9081b43eb2cd3b974</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a35f2e143409f7fcbf6e1b3c610f4a676</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a684fc9ecf86f1dc99462faa7aece791b</anchor>
      <arglist>(Teuchos::ParameterList &amp;timingParameterList) const</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aade60205a3b6ece6bdc76f5d1555184e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Superludist</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a73ffe10366949086d5c43e9cdebd1915</anchor>
      <arglist>(Teuchos::RCP&lt; Matrix &gt; A, Teuchos::RCP&lt; Vector &gt; X, Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Superludist</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a8a5b2d3aeb58029d97f7ca548de36cf7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>preOrdering</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ace48a1882df1a866f7fba00e321d1abe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>symbolicFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ab0602bd80f5ba358adf710423c013186</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>numericFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a00d8910072b4ff26abbb71e81ed30ef9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a89ed4e401b4836bee75de516b6de12cd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a7ab43f4fd19cb4f3e7b9c1068d78f381</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; const Vector &gt; B) const</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>setParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a2f812b0093f7491cf86574eb5f5d1117</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac42d1a2b4c037065b2637c06b23b06bd</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a53d48aa613aca50ed957272aeeeb9485</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>getNonconstParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bf21b207f1e9eb860fe326e3fb6dcbc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>unsetParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4a62e2cbe4d173534302af0f93088a05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a6b805f6cab10b90c0843f908e00e0344</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSymbolicFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a559b355352ff76b91d7ce2cc94b6242e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumNumericFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a32637ffd666f7af4726811eb7fdaa47a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSolve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a942a1be951a2bb725aef7039219e6381</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const char *</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ae66e66d962f94727d9b7545d8fe7daa8</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>preOrderingDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9f679ee04f9db29ac263dead626a6c76</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>symbolicFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bb9fc2c4313ba09f8b74b3c190f1704</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>numericFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a549a2117cc7b78842f55627be58cffe2</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; MatrixAdapter&lt; Matrix &gt; &gt;</type>
      <name>matrixA_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a73862ad4870aaad942cf1d85910e035e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; Vector &gt;</type>
      <name>multiVecX_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a57e4fe1d6563e004f933a390d89f5d6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; const Vector &gt;</type>
      <name>multiVecB_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a8948ea996ec300e9423b3c6977197b0b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumRows_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4f98ce6a877b2f6ece913c05adef6fc4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumCols_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9c9458f37c6bfe6aa23be99213265ebf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumNonZeros_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e824badc831fc6cd4627ea834b3782d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Status</type>
      <name>status_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>af688193f21876bb4a209738d2bfeae3e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Control</type>
      <name>control_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae9cc6813b5540966424eaab3cd9b0e15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Timers</type>
      <name>timers_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aa6a15406034089ec343445003db0a276</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef" protection="private">
      <type>TypeMap&lt; Amesos::Superludist, scalar_type &gt;::magnitude_type</type>
      <name>magnitude_type</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a567ca0120342e36a673d9523a4bbb693</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>preOrdering_impl</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>acb7c22de662b13e6e242c7a5c37fb5c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>symbolicFactorization_impl</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a5978f915f8569b0814d6659d7b89d6d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>numericFactorization_impl</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>aaaf3f59d1f6ffdee5ec098c5215083fc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>solve_impl</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>abfdfba9fc1c9caa1df8c383c608bda21</anchor>
      <arglist>(const Teuchos::Ptr&lt; MultiVecAdapter&lt; Vector &gt; &gt; X, const Teuchos::Ptr&lt; const MultiVecAdapter&lt; Vector &gt; &gt; B)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>matrixShapeOK_impl</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ab823edbb4ad2c26ec6e754ed46c63123</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>setParameters_impl</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a7ba4b666a671e960a84ac5ba20aede3e</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters_impl</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>aba72ca535646011f0b96b9a36da6b40b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>get_default_grid_size</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ab566386bc63131028e72ce0433633653</anchor>
      <arglist>(int nprocs, SLUD::int_t &amp;nprow, SLUD::int_t &amp;npcol) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>struct Amesos::Superludist::SLUData</type>
      <name>data_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a125b8b0589514cbe241c687ae75f5e65</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; slu_dist_type &gt;</type>
      <name>nzvals_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>af65d767b320e2ad7468c588b2ead7dc6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; int &gt;</type>
      <name>rowind_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ad9aece372dd1be6810e89cda107c43cb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; int &gt;</type>
      <name>colptr_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a6a04ff115f1171ee0e52d534f09f679a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; slu_dist_type &gt;</type>
      <name>bxvals_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a77d1730ec10baccaa1490f18748978d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>size_t</type>
      <name>ldbx_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a3cfb3f69390a04b6d0eab8bf136a3f32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>in_grid</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>a6f15074ad903d84e03fede9009c150bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; const Tpetra::Map&lt; local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>superlu_rowmap_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ad85e9237d848f421be6b7e82eafe4df5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; const MatrixAdapter&lt; matrix_type &gt; &gt;</type>
      <name>redist_A_</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ae78230a30d9ad8d2f2d94a42af50bfbd</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>SolverCore&lt; Amesos::Superludist, Matrix, Vector &gt;</name>
      <anchorfile>classAmesos_1_1Superludist.html</anchorfile>
      <anchor>ad28e9f9b31e39feeb0c58d15d7fdf934</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a415501e885f7e753298fa2cf61793932</anchor>
      <arglist>(Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a479b23341a0b171d241cffc87a92bcf0</anchor>
      <arglist>(Teuchos::RCP&lt; Matrix &gt; A, Teuchos::RCP&lt; Vector &gt; X, Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ad8aa0892d3a641a41d18884185d31ab6</anchor>
      <arglist>(const char *solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a53132a9d360914a443142b7101294a07</anchor>
      <arglist>(const char *solverName, const Teuchos::RCP&lt; Matrix &gt; A, const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>abb8acb7cceac1ce640a41dd05cc6024b</anchor>
      <arglist>(const std::string solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>aea67ca5c343fc0f248abb8bd59423098</anchor>
      <arglist>(const std::string solverName, const Teuchos::RCP&lt; Matrix &gt; A, const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ab239b521c543b7cbed4f8823bb84c23f</anchor>
      <arglist>(const std::string solverName, Matrix *A)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>afed8390a758da6367f4f28c1602df807</anchor>
      <arglist>(const std::string solverName, const Teuchos::RCP&lt; Matrix &gt; A)</arglist>
    </member>
    <docanchor file="classAmesos_1_1Superludist">slu_dist_options</docanchor>
  </compound>
  <compound kind="struct">
    <name>Amesos::FunctionMap&lt; Superludist, Scalar &gt;</name>
    <filename>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</filename>
    <templarg></templarg>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>distribute</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>a284a45002dfa75d1f661f12d9b4c9ca6</anchor>
      <arglist>(SLUD::fact_t fact, SLUD::int_t n, SLUD::SuperMatrix *A, SLUD::Glu_freeable_t *glu_freeable, SLUD::LUstruct_t *lu, SLUD::gridinfo_t *grid)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>distribute_loc</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>a7a40ff3df17d75ca1ed340c0fecb3a8d</anchor>
      <arglist>(SLUD::fact_t fact, SLUD::int_t n, SLUD::SuperMatrix *A, SLUD::ScalePermstruct *scale_perm, SLUD::Glu_freeable_t *glu_freeable, SLUD::LUstruct_t *lu, SLUD::gridinfo_t *grid)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gstrf</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>ad56c74fa9d61e2daf669d37a0c97f296</anchor>
      <arglist>(SLUD::superlu_options_t *options, int m, int n, double anorm, SLUD::LUstruct_t *LU, SLUD::gridinfo_t *grid, SLUD::SuperLUStat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gstrs</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>ac5753c3badcadb36650415b9bbaffcba</anchor>
      <arglist>(SLUD::int_t n, SLUD::LUStruct_t *lu_struct, SLUD::ScalePermstruct_t *scale_perm_struct, SLUD::gridinfo_t *grid, typename TypeMap&lt; Superludist, Scalar &gt;::type *B, SLUD::int_t l_numrows, SLUD::int_t fst_global_row, SLUD::int_t ldb, int nrhs, SLUD::SOLVEstruct_t *solve_struct, SLUD::SuperLUStat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gstrs_Bglobal</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>ae80c40d0739b064e0f2247d64ff8d39c</anchor>
      <arglist>(SLUD::int_t n, SLUD::LUStruct_t *lu_struct, SLUD::gridinfo_t *grid, typename TypeMap&lt; Superludist, Scalar &gt;::type *B, SLUD::int_t ldb, int nrhs, SLUD::SuperLUStat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gsrfs</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>a9e606dbd27896d8cb852998d97450428</anchor>
      <arglist>(SLUD::int_t n, SLUD::SuperMatrix *A, double anorm, SLUD::LUStruct_t *lu_struct, SLUD::ScalePermstruct_t *scale_perm, SLUD::gridinfo_t *grid, typename TypeMap&lt; Superludist, Scalar &gt;::type *B, SLUD::int_t ldb, typename TypeMap&lt; Superludist, Scalar &gt;::type *X, SLUD::int_t ldx, int nrhs, SLUD::SOLVEstruct_t *solve_struct, double *berr, SLUD::SuperLUStat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gsrfs_AXBglobal</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>a59b20565434519ba1db6deb29a8918b2</anchor>
      <arglist>(SLUD::int_t n, SLUD::SuperMatrix *A, double anorm, SLUD::LUStruct_t *lu_struct, SLUD::gridinfo_t *grid, typename TypeMap&lt; Superludist, Scalar &gt;::type *B, SLUD::int_t ldb, typename TypeMap&lt; Superludist, Scalar &gt;::type *X, SLUD::int_t ldx, int nrhs, double *berr, SLUD::SuperLUStat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_CompRowLoc_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>aae20f0ea68db2609325f8f428decb886</anchor>
      <arglist>(SLUD::SuperMatrix *A, SLUD::int_t g_numrows, SLUD::int_t g_numcols, SLUD::int_t l_nnz, SLUD::int_t l_numrows, SLUD::int_t fst_global_row, typename TypeMap&lt; Superludist, Scalar &gt;::type *nzval, SLUD::int_t *colind, SLUD::int_t *rowptr, SLUD::Stype_t storage_t, SLUD::Dtype_t data_t, SLUD::Mtype_t mat_t)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_CompCol_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>ae9bb33c0c534f0f1b4a8f02255f8e5e3</anchor>
      <arglist>(SLUD::SuperMatrix *A, SLUD::int_t numrows, SLUD::int_t numcols, SLUD::int_t nnz, typename TypeMap&lt; Superludist, Scalar &gt;::type *nzval, SLUD::int_t *rowind, SLUD::int_t *colptr, SLUD::Stype_t storage_t, SLUD::Dtype_t data_t, SLUD::Mtype_t mat_t)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_Dense_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>a7e9d06ffd32c4adbae9edd4d3e56c0ab</anchor>
      <arglist>(SLUD::SuperMatrix *X, int numrows, int numcols, typename TypeMap&lt; Superludist, Scalar &gt;::type *x, int ldx, SLUD::Stype_t storage_t, SLUD::Dtype_t data_t, SLUD::Mtype_t mat_t)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gsequ_loc</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>a6d78858eb9678637cd4f974cedcbb062</anchor>
      <arglist>(SLUD::SuperMatrix *A, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *r, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *c, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *rowcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *colcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *amax, int *info, SLUD::gridinfo_t *grid)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gsequ</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>a8cd73c22f4aeadb5a780e92957ee50f9</anchor>
      <arglist>(SLUD::SuperMatrix *A, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *r, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *c, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *rowcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *colcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *amax, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>laqgs_loc</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>abb5a24d6025950e15cf5b2239a097a59</anchor>
      <arglist>(SLUD::SuperMatrix *A, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *r, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *c, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type rowcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type colcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type amax, SLUD::DiagScale_t *equed)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>laqgs_loc</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superludist_00_01Scalar_01_4.html</anchorfile>
      <anchor>abb5a24d6025950e15cf5b2239a097a59</anchor>
      <arglist>(SLUD::SuperMatrix *A, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *r, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type *c, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type rowcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type colcnd, typename TypeMap&lt; Superludist, Scalar &gt;::magnitude_type amax, SLUD::DiagScale_t *equed)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::Superlumt</name>
    <filename>classAmesos_1_1Superlumt.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <base>SolverCore&lt; Amesos::Superlumt, Matrix, Vector &gt;</base>
    <member kind="typedef">
      <type>Superlumt&lt; Matrix, Vector &gt;</type>
      <name>solver_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>aa72d6585621e28544c189ce66fcd08d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Matrix</type>
      <name>matrix_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a765b2daf180a22faae9bf9ea75da1281</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Vector</type>
      <name>vector_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a447270e7a37e14e72e32529f8dbd6d05</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::scalar_t</type>
      <name>scalar_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a66f90d82aa46d1d563dc7e83fa710265</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::local_ordinal_t</type>
      <name>local_ordinal_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>ad6e7c39c3b37e2b484d83c62d992720f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::global_ordinal_t</type>
      <name>global_ordinal_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>aa89f626d1e8afa7d7800c9bbbecba712</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixTraits&lt; matrix_type &gt;::node_t</type>
      <name>node_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a9cb38937dd1d8a394ff7d97076279b99</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>MatrixAdapter&lt; matrix_type &gt;::global_size_t</type>
      <name>global_size_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>ad7d60fbebe2ad63bf59ed36f9edfe6c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SolverCore&lt; Amesos::Superlumt, Matrix, Vector &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e3660d3adce807015ae6b1490fedb16</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Solver&lt; Matrix, Vector &gt;</type>
      <name>super_type</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae0a58e56723c5dd02253c8e7ea2b05fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>matrixShapeOK</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a99a5d733412d7b9f8e7ff692e6c27c2a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setA</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a72d655870f183deb536d126fe161a3a6</anchor>
      <arglist>(const Teuchos::RCP&lt; Matrix &gt; a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac36fe794be11b6d03547e71dd8f89ac4</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; x)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; Vector &gt;</type>
      <name>getX</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeee7773536622735fe0f2b3242521e05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>afc887943fe97b8bbc16809606692ac0f</anchor>
      <arglist>(const Teuchos::RCP&lt; const Vector &gt; b)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Vector &gt;</type>
      <name>getB</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aeb10617c7d0feab6bfab5c8d4240ac9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a28df22e8d7ea6c53d4f38fe08c253b5b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9ffe894f7d3973a9081b43eb2cd3b974</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a35f2e143409f7fcbf6e1b3c610f4a676</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getTiming</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a684fc9ecf86f1dc99462faa7aece791b</anchor>
      <arglist>(Teuchos::ParameterList &amp;timingParameterList) const</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aade60205a3b6ece6bdc76f5d1555184e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Superlumt</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a8287f6cf7afd882cdc56ef879c9c19ad</anchor>
      <arglist>(Teuchos::RCP&lt; Matrix &gt; A, Teuchos::RCP&lt; Vector &gt; X, Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Superlumt</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>aef7a5948be3ad71804402c1c5cb0e46e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>preOrdering</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ace48a1882df1a866f7fba00e321d1abe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>symbolicFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ab0602bd80f5ba358adf710423c013186</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>numericFactorization</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a00d8910072b4ff26abbb71e81ed30ef9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a89ed4e401b4836bee75de516b6de12cd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>solve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a7ab43f4fd19cb4f3e7b9c1068d78f381</anchor>
      <arglist>(const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; const Vector &gt; B) const</arglist>
    </member>
    <member kind="function">
      <type>super_type &amp;</type>
      <name>setParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a2f812b0093f7491cf86574eb5f5d1117</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ac42d1a2b4c037065b2637c06b23b06bd</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a53d48aa613aca50ed957272aeeeb9485</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>getNonconstParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bf21b207f1e9eb860fe326e3fb6dcbc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Teuchos::ParameterList &gt;</type>
      <name>unsetParameterList</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4a62e2cbe4d173534302af0f93088a05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a6b805f6cab10b90c0843f908e00e0344</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSymbolicFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a559b355352ff76b91d7ce2cc94b6242e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumNumericFact</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a32637ffd666f7af4726811eb7fdaa47a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumSolve</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a942a1be951a2bb725aef7039219e6381</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const char *</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a04d5cf1ab860dff836a2f3f059ac7563</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>preOrderingDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9f679ee04f9db29ac263dead626a6c76</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>symbolicFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a3bb9fc2c4313ba09f8b74b3c190f1704</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>numericFactorizationDone</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a549a2117cc7b78842f55627be58cffe2</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; MatrixAdapter&lt; Matrix &gt; &gt;</type>
      <name>matrixA_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a73862ad4870aaad942cf1d85910e035e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; Vector &gt;</type>
      <name>multiVecX_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a57e4fe1d6563e004f933a390d89f5d6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Teuchos::RCP&lt; const Vector &gt;</type>
      <name>multiVecB_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a8948ea996ec300e9423b3c6977197b0b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumRows_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a4f98ce6a877b2f6ece913c05adef6fc4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumCols_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a9c9458f37c6bfe6aa23be99213265ebf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>global_size_type</type>
      <name>globalNumNonZeros_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>a1e824badc831fc6cd4627ea834b3782d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Status</type>
      <name>status_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>af688193f21876bb4a209738d2bfeae3e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Control</type>
      <name>control_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>ae9cc6813b5540966424eaab3cd9b0e15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Timers</type>
      <name>timers_</name>
      <anchorfile>classAmesos_1_1SolverCore.html</anchorfile>
      <anchor>aa6a15406034089ec343445003db0a276</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef" protection="private">
      <type>TypeMap&lt; Amesos::Superlumt, scalar_type &gt;::magnitude_type</type>
      <name>magnitude_type</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a35e2a94313ba2582db3a85b6deae631a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>preOrdering_impl</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>ad62d2a3bb7d8984c99aafe48626bf63e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>symbolicFactorization_impl</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a70d29c43e206a87766b2849c8e77ff2c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>numericFactorization_impl</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a1610f79bf3b34ae9d8ad04c746bbb4f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>solve_impl</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a6a104a2dbb4e280a4a9df82cec51d821</anchor>
      <arglist>(const Teuchos::Ptr&lt; MultiVecAdapter&lt; Vector &gt; &gt; X, const Teuchos::Ptr&lt; MultiVecAdapter&lt; Vector &gt; &gt; B) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>matrixShapeOK_impl</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a46d87d71c92e1409bc4dbcafb0da322d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>setParameters_impl</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a3f626a55a7a166f0ed7881846c202531</anchor>
      <arglist>(const Teuchos::RCP&lt; Teuchos::ParameterList &gt; &amp;parameterList)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Teuchos::RCP&lt; const Teuchos::ParameterList &gt;</type>
      <name>getValidParameters_impl</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>af0fabaf93d9d15c078660be7f67bb117</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>struct Amesos::Superlumt::SLUData</type>
      <name>data_</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a7335b8921db86343baa06233f99a9dd0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; typename TypeMap&lt; Amesos::Superlumt, scalar_type &gt;::type &gt;</type>
      <name>nzvals_</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a95c4e4cf9ba3fcda93eea04f3c8a816c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; int &gt;</type>
      <name>rowind_</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>ad2e63c87f259163e5dad00b8bebce2b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::Array&lt; int &gt;</type>
      <name>colptr_</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>af39a59034242ff1c24cef8fbb9efa3ed</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>SolverCore&lt; Amesos::Superlumt, Matrix, Vector &gt;</name>
      <anchorfile>classAmesos_1_1Superlumt.html</anchorfile>
      <anchor>a248935705fdcc5bf341c3cdd2d926ee0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a415501e885f7e753298fa2cf61793932</anchor>
      <arglist>(Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a479b23341a0b171d241cffc87a92bcf0</anchor>
      <arglist>(Teuchos::RCP&lt; Matrix &gt; A, Teuchos::RCP&lt; Vector &gt; X, Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ad8aa0892d3a641a41d18884185d31ab6</anchor>
      <arglist>(const char *solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>a53132a9d360914a443142b7101294a07</anchor>
      <arglist>(const char *solverName, const Teuchos::RCP&lt; Matrix &gt; A, const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>abb8acb7cceac1ce640a41dd05cc6024b</anchor>
      <arglist>(const std::string solverName, Matrix *A, Vector *X, Vector *B)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>aea67ca5c343fc0f248abb8bd59423098</anchor>
      <arglist>(const std::string solverName, const Teuchos::RCP&lt; Matrix &gt; A, const Teuchos::RCP&lt; Vector &gt; X, const Teuchos::RCP&lt; Vector &gt; B)</arglist>
    </member>
    <member kind="function">
      <type>Solver&lt; Matrix, Vector &gt; *</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>ab239b521c543b7cbed4f8823bb84c23f</anchor>
      <arglist>(const std::string solverName, Matrix *A)</arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Solver&lt; Matrix, Vector &gt; &gt;</type>
      <name>create</name>
      <anchorfile>classAmesos_1_1Solver.html</anchorfile>
      <anchor>afed8390a758da6367f4f28c1602df807</anchor>
      <arglist>(const std::string solverName, const Teuchos::RCP&lt; Matrix &gt; A)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Amesos::FunctionMap&lt; Superlumt, Scalar &gt;</name>
    <filename>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</filename>
    <templarg></templarg>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gssvx</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</anchorfile>
      <anchor>a878288e7256c0243625c2b3550ddfb53</anchor>
      <arglist>(SLUMT::superlumt_options_t *options, SLUMT::SuperMatrix *A, int *perm_c, int *perm_r, int *etree, SLUMT::equed_t *equed, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *R, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *C, SLUMT::SuperMatrix *L, SLUMT::SuperMatrix *U, void *work, int lwork, SLUMT::SuperMatrix *B, SLUMT::SuperMatrix *X, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *recip_pivot_growth, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *rcond, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *ferr, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *berr, SLUMT::superlu_memusage_t *mem_usage, SLUMT::Gstat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gstrs</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</anchorfile>
      <anchor>a1db42f9b033a4977ee5dde68e740f8e7</anchor>
      <arglist>(SLUMT::trans_t trans, SLUMT::SuperMatrix *L, SLUMT::SuperMatrix *U, int *perm_r, int *perm_c, SLUMT::SuperMatrix *B, SLUMT::Gstat_t *Gstat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gstrf</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</anchorfile>
      <anchor>a8ce194a9ca2fb0e0254d0f27516ab565</anchor>
      <arglist>(SLUMT::superlumt_options_t *options, SLUMT::SuperMatrix *A, int *perm_r, SLUMT::SuperMatrix *L, SLUMT::SuperMatrix *U, SLUMT::Gstat_t *stat, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_CompCol_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</anchorfile>
      <anchor>ab9319a70b3476203a728d9a4426dd02a</anchor>
      <arglist>(SLUMT::SuperMatrix *A, int numrows, int numcols, int nnz, typename TypeMap&lt; Superlumt, Scalar &gt;::type *nzval, int *rowind, int *colptr, SLUMT::Stype_t storage_t, SLUMT::Dtype_t data_t, SLUMT::Mtype_t mat_t)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>create_Dense_Matrix</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</anchorfile>
      <anchor>a9f89d29a4897a3a9dac42418aceb5898</anchor>
      <arglist>(SLUMT::SuperMatrix *X, int numrows, int numcols, typename TypeMap&lt; Superlumt, Scalar &gt;::type *x, int ldx, SLUMT::Stype_t storage_t, SLUMT::Dtype_t data_t, SLUMT::Mtype_t mat_t)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>gsequ</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</anchorfile>
      <anchor>abb1ed911ecbb96843720571eb7aa2441</anchor>
      <arglist>(SLUMT::SuperMatrix *A, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *r, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *c, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *rowcnd, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *colcnd, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *amax, int *info)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>laqgs</name>
      <anchorfile>structAmesos_1_1FunctionMap_3_01Superlumt_00_01Scalar_01_4.html</anchorfile>
      <anchor>ade5db113efe89e4f71d3c1493f6f1e54</anchor>
      <arglist>(SLUMT::SuperMatrix *A, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *r, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type *c, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type rowcnd, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type colcnd, typename TypeMap&lt; Superlumt, Scalar &gt;::magnitude_type amax, SLUMT::equed_t *equed)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::ConcreteMatrixAdapter&lt; Tpetra::CrsMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps &gt; &gt;</name>
    <filename>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <base>AbstractConcreteMatrixAdapter&lt; Tpetra::RowMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt;, Tpetra::CrsMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps &gt; &gt;</base>
    <member kind="typedef">
      <type>Tpetra::CrsMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps &gt;</type>
      <name>matrix_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a7b1a382f0dbb1e0aa138068a69e21e2a</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::scalar_t</type>
      <name>scalar_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a152c8c1c25a000693db287bf3f63b2bd</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::local_ordinal_t</type>
      <name>local_ordinal_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a31b3858647e9ea52d1567fa555a5782c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::global_ordinal_t</type>
      <name>global_ordinal_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a6a93c8f9fe29d1ae568004a327c451f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::node_t</type>
      <name>node_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>aebb20994cbd1b694488b8e1bd3b234db</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::global_size_t</type>
      <name>global_size_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a9105bbf0291c528505f9a3795415a9fc</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>LocalMatOps</type>
      <name>local_mat_ops_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a247011aa0c827302fb74d5ce06fc00c4</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ConcreteMatrixAdapter&lt; matrix_t &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>aac4d0910ae0c697cc5db42bc4d8ddaa9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ConcreteMatrixAdapter</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a4969e1b53d8c19c1cb3f526cb7fe669f</anchor>
      <arglist>(RCP&lt; matrix_t &gt; m)</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const MatrixAdapter&lt; matrix_t &gt; &gt;</type>
      <name>get_impl</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a881fc4befded1ee7b19f0d6597fee4e1</anchor>
      <arglist>(const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; map) const </arglist>
    </member>
    <member kind="typedef" protection="private">
      <type>AbstractConcreteMatrixAdapter&lt; Tpetra::RowMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt;, matrix_t &gt;</type>
      <name>super_t</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>a61364cd49f25aa375a652972f77a2552</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>MatrixAdapter&lt; Tpetra::RowMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt; &gt;</name>
      <anchorfile>classAmesos_1_1ConcreteMatrixAdapter_3_01Tpetra_1_1CrsMatrix_3_01Scalar_00_01LocalOrdinal_00_01Gf4979c56a0bdb61ece67db1bf2d7fa55.html</anchorfile>
      <anchor>aa937545c52ba529b494370bda6c35e49</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::MultiVecAdapter&lt; Tpetra::MultiVector&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt; &gt;</name>
    <filename>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <member kind="typedef">
      <type>Tpetra::MultiVector&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt;</type>
      <name>multivec_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>ab7036d111fa1354c5c96e51753b4f9d5</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Scalar</type>
      <name>scalar_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>ae0bd1754c4ced72c6986db459386261c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>LocalOrdinal</type>
      <name>local_ordinal_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a23866e0ea9ed81ed4bc9316530587ed8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>GlobalOrdinal</type>
      <name>global_ordinal_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a030dca814302cddb8112c7408e21e333</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Node</type>
      <name>node_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>afba6fc88b8df7e0f3338754cc13e3061</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Tpetra::global_size_t</type>
      <name>global_size_type</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>ada3e4cadf8ab76d13e262a4367aa424c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MultiVecAdapter</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a1f12d33a2ccdd0877f08ab39186ccf69</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MultiVecAdapter</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a50b9f031b4ffe3bfc17d11256c17c4ff</anchor>
      <arglist>(const MultiVecAdapter&lt; multivec_type &gt; &amp;adapter)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MultiVecAdapter</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a0d771bf67708360ac244af8fa6837e8c</anchor>
      <arglist>(const Teuchos::RCP&lt; multivec_type &gt; &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; multivec_type &gt;</type>
      <name>getAdaptee</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>aac3ab9947f7d3b040465daef9ccd96a4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MultiVecAdapter&lt; multivec_type &gt; &amp;</type>
      <name>scale</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a296f6f635ff3213f0ccdceee66a51594</anchor>
      <arglist>(const scalar_type alpha)</arglist>
    </member>
    <member kind="function">
      <type>MultiVecAdapter&lt; multivec_type &gt; &amp;</type>
      <name>update</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>af6ce15ceb96d0f26a5dbc1836ce9c88b</anchor>
      <arglist>(const scalar_type beta, const MultiVecAdapter&lt; multivec_type &gt; &amp;B, const scalar_type alpha)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isLocal</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a5a2191800c2440d04dc9b80b68534340</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Tpetra::Map&lt; local_ordinal_type, global_ordinal_type, node_type &gt; &gt; &amp;</type>
      <name>getMap</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>ac303e6e711fd55019a60dada723b13f0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Teuchos::RCP&lt; const Teuchos::Comm&lt; int &gt; &gt; &amp;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a621cb86a3fc97caa2d686f330bf1cbf3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalLength</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a4b59ba436480d11809dc9bd9de5f4f27</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumVectors</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a8ca52bef862018ccc97835b28dca9e9f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>global_size_type</type>
      <name>getGlobalLength</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a33a32ba84e7aaa0e33caf7589eecb194</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getGlobalNumVectors</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>aa7bf036d0f61dd401c4d6310848c3c1a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getStride</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a59a9604150aaf0cb99dd2e3165eb08ca</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isConstantStride</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>abb2653ac8293b92483530cfa1c3176d8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; const Tpetra::Vector&lt; scalar_type, local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>getVector</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a5083902bc213d340dcae9ac1e551e58f</anchor>
      <arglist>(size_t j) const </arglist>
    </member>
    <member kind="function">
      <type>Teuchos::RCP&lt; Tpetra::Vector&lt; scalar_type, local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>getVectorNonConst</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>ad1f7737101b8e384a73604b318589135</anchor>
      <arglist>(size_t j)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get1dCopy</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a43e0417ee91006a8173d15ddf067b073</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_type &gt; &amp;A, size_t lda, bool global_copy) const </arglist>
    </member>
    <member kind="function">
      <type>Teuchos::ArrayRCP&lt; scalar_type &gt;</type>
      <name>get1dViewNonConst</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a84b98f836b4be505ae855efd3e717352</anchor>
      <arglist>(bool local=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get2dCopy</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a777ed90af76be6df3a9af75ab5fbba78</anchor>
      <arglist>(Teuchos::ArrayView&lt; const Teuchos::ArrayView&lt; scalar_type &gt; &gt; A) const </arglist>
    </member>
    <member kind="function">
      <type>Teuchos::ArrayRCP&lt; Teuchos::ArrayRCP&lt; scalar_type &gt; &gt;</type>
      <name>get2dViewNonConst</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>aa98ef24e8e9c18250b1dc82cc44a7f5c</anchor>
      <arglist>(bool local=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>globalize</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>ab01ac4ce8a6571a694eccca504a3afe5</anchor>
      <arglist>(int root=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>globalize</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a9c8580f7b32f9d45b4cf5fca2bd2ab21</anchor>
      <arglist>(const Teuchos::ArrayView&lt; Value_t &gt; &amp;newVals, int root=0)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a48b50229ec681f6ecacd0c9c5b6eeeed</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a265a83b3cbe3377037851db84c3b15d5</anchor>
      <arglist>(Teuchos::FancyOStream &amp;os, const Teuchos::EVerbosityLevel verbLevel) const </arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const char *</type>
      <name>name</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a8ca8cc5a38fe566499782846eaf772c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>localize</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a0714342928c104056d86a7e8b406eed6</anchor>
      <arglist>(bool root) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Tpetra::MultiVector&lt; scalar_type, local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>mv_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>ad8f2c794f494044456eaf91d6e07c8c2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Tpetra::MultiVector&lt; scalar_type, local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>l_mv_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>aec5c7763f26be14a9dbb565a51acdbb1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Tpetra::MultiVector&lt; scalar_type, local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>l_l_mv_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a49cf65e22c6c0f8dd749dfb44aaa2f09</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Tpetra::Import&lt; local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>importer_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a0f27587798cc46697f6a0a2b66a19399</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; Tpetra::Export&lt; local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>exporter_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a8aba0b9104f71c479e444b32c42d60e5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; const Tpetra::Map&lt; local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>l_map_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a75f526cd1b9277393957427b10306579</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Teuchos::RCP&lt; const Tpetra::Map&lt; local_ordinal_type, global_ordinal_type, node_type &gt; &gt;</type>
      <name>o_map_</name>
      <anchorfile>classAmesos_1_1MultiVecAdapter_3_01Tpetra_1_1MultiVector_3_01Scalar_00_01LocalOrdinal_00_01GlobalOrdinal_00_01Node_01_4_01_4.html</anchorfile>
      <anchor>a6ef6004c5a4e08b8d03b03e6c65edd93</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Amesos::AbstractConcreteMatrixAdapter&lt; Tpetra::RowMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt;, DerivedMat &gt;</name>
    <filename>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <base>MatrixAdapter&lt; DerivedMat &gt;</base>
    <member kind="typedef">
      <type>Tpetra::RowMatrix&lt; Scalar, LocalOrdinal, GlobalOrdinal, Node &gt;</type>
      <name>matrix_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a60d8eabc29c12cc9413544e20a663053</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Scalar</type>
      <name>scalar_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a01091a7f22a60f08b44cb18f8da68559</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>LocalOrdinal</type>
      <name>local_ordinal_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a77c78ea19be050c401685ae3b27a804c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>GlobalOrdinal</type>
      <name>global_ordinal_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a151184b3a2f5a906971697deac16e391</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Node</type>
      <name>node_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>aa6aeced853d3234269b22dbf70b0a354</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>super_t::global_size_t</type>
      <name>global_size_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>abf7a417c6e40df87c1677de474d6907b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>AbstractConcreteMatrixAdapter&lt; matrix_t, DerivedMat &gt;</type>
      <name>type</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>ac2de87006c2613c4c0c7954d92314171</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Util::no_special_impl</type>
      <name>get_crs_spec</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>aaedf677389f092d1a5b3ed92e98ca598</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Util::no_special_impl</type>
      <name>get_ccs_spec</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>ada9cf38bda2a809ea95a2508c52b8195</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Util::row_access</type>
      <name>major_access</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>ab2b1e9a14289d118f4de04144152ce82</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>ConcreteMatrixAdapter&lt; DerivedMat &gt;</type>
      <name>adapter_t</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a3bfa635f15a7ee77ba78dc64f5822e9c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>AbstractConcreteMatrixAdapter</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a4d35a94ba18a68f78e5ebe992baca0e5</anchor>
      <arglist>(RCP&lt; matrix_t &gt; m)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getGlobalRowCopy_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a3cd0d7da520b2bb8cbf6ce0edfa20a0b</anchor>
      <arglist>(global_ordinal_t row, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getGlobalColCopy_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a5f3daa31a5004e04b589968f0de0b8bd</anchor>
      <arglist>(global_ordinal_t col, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const </arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a1a69b62f70c99487b1c6c3d903ede46f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a65281fe0d3f255c3d92b25e8f91af78d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMaxRowNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>afebde8bc3eecca42708a57822101a501</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMaxColNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a3734101273fa4a88ff8659225721d939</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getGlobalRowNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>aa391900e4586fee8c1370d46e1a3d9c6</anchor>
      <arglist>(global_ordinal_t row) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalRowNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a57587adb86365c8c9abe1f556f5437d8</anchor>
      <arglist>(local_ordinal_t row) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getGlobalColNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a4fe1179ce04db43dad17ce2520ee8a6d</anchor>
      <arglist>(global_ordinal_t col) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalColNNZ_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>aad9186ecdd22763f44daea143d330937</anchor>
      <arglist>(local_ordinal_t col) const </arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>getRowMap_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>add5404f4abe4d48dc3eea59de8f9776d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>getColMap_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a4652898e083b80b9b2e566fe24845c96</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a0352a08791fe2a88334d737a44329457</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isLocallyIndexed_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>aea49903ec1cdbd791b1ececf47f29a1a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isGloballyIndexed_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a4a70de2e4acb97caa3ade99ceed7b6ef</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const super_t &gt;</type>
      <name>get_impl</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>adc8d456018fad3474ae8e78f28acd6b7</anchor>
      <arglist>(const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; map) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a7da98dde7decdcd4006c4203193a1edb</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; rowmap, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCrs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a9b9abaec0efdf184ac502bcf0ef08b3c</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; colind, const Teuchos::ArrayView&lt; global_size_t &gt; rowptr, global_size_t &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a54a1fb9c3f9203c6a7ec4ba9c87dae25</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; colmap, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getCcs</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad85223e6476499e4ab135661bf8339aa</anchor>
      <arglist>(const Teuchos::ArrayView&lt; scalar_t &gt; nzval, const Teuchos::ArrayView&lt; global_ordinal_t &gt; rowind, const Teuchos::ArrayView&lt; global_size_t &gt; colptr, global_size_t &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Util::Arbitrary) const</arglist>
    </member>
    <member kind="function">
      <type>const RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>getComm</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aae644193c73b59c573720132be955431</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNumRows</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ac729469ca361ac261590f143624fa6b4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNumCols</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a2c1f4eed533cf858db94067149c70c0b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>global_size_t</type>
      <name>getGlobalNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aa144b89709fd20ed4627d5966997132e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumRows</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aa72c5a28011f01c5f09f1be50df5de2e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNumCols</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ada217a4504703e97a0ed2e165ec9b766</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getLocalNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ac1fcfa87de3617bc6fee93f884ca6bf7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t &gt; &gt;</type>
      <name>getRowMap</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a079a384379ea6ba42f82f254ec322ab2</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t &gt; &gt;</type>
      <name>getColMap</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad095ea1b4a87c53258bc95344540d369</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>description</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aac6aec881feab1fa3d30c4c3a7d3a2f1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>describe</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a0e6d8d17887b32364e946403a2617365</anchor>
      <arglist>(Teuchos::FancyOStream &amp;out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getGlobalRowCopy</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a82a04efccdb9fc4c677f6971b0f2a9ff</anchor>
      <arglist>(global_ordinal_t row, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>getGlobalColCopy</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a1930fc8e940fa82ae8933c8fbd61c58d</anchor>
      <arglist>(global_ordinal_t col, const Teuchos::ArrayView&lt; global_ordinal_t &gt; &amp;indices, const Teuchos::ArrayView&lt; scalar_t &gt; &amp;vals, size_t &amp;nnz) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getMaxRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aec3505d353a2b1bfac631222aa6e582d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getMaxColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>af4fca7c867756723a39d7251eb5e214e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getGlobalRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a8d5fb3fc6c2db1015f75399f880d8e1e</anchor>
      <arglist>(global_ordinal_t row) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getLocalRowNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a6286f7e449913ec698831cbe9f855ed0</anchor>
      <arglist>(local_ordinal_t row) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getGlobalColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a139320c33a3e6954e1263ac083c56fe9</anchor>
      <arglist>(global_ordinal_t col) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>size_t</type>
      <name>getLocalColNNZ</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a59b1d89868a0044cdb5d6442c505ba21</anchor>
      <arglist>(local_ordinal_t col) const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>isLocallyIndexed</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a37bec34548dee2f1bfb9ab42f5c872a6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>isGloballyIndexed</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a760657acc1d3e7abfcb35b689ed7c3a5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>RCP&lt; const type &gt;</type>
      <name>get</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ab2e42b619d1ecc2796bce408615bd15d</anchor>
      <arglist>(const Teuchos::Ptr&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt; map) const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const RCP&lt; const DerivedMat &gt;</type>
      <name>mat_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>aaa7c67b31929c0252ceb73e3a5fc8abf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>row_map_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a91680cc808fdc40b64a01ab7ca6ca554</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Tpetra::Map&lt; local_ordinal_t, global_ordinal_t, node_t &gt; &gt;</type>
      <name>col_map_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a06617c51632049ccc3dee8fe7cdcf4f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>RCP&lt; const Teuchos::Comm&lt; int &gt; &gt;</type>
      <name>comm_</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>a898d99aef23d15a31a92fc54ea9a5116</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef" protection="private">
      <type>MatrixAdapter&lt; DerivedMat &gt;</type>
      <name>super_t</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>aaec0928265733ba5daf74cb5c171f41d</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>MatrixAdapter&lt; DerivedMat &gt;</name>
      <anchorfile>classAmesos_1_1AbstractConcreteMatrixAdapter_3_01Tpetra_1_1RowMatrix_3_01Scalar_00_01LocalOrdinadbc14e7d1b871e5aa0964a7a23b8b63d.html</anchorfile>
      <anchor>a0a0a3629b22324189309040873178be5</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend class</type>
      <name>Util::get_cxs_helper</name>
      <anchorfile>classAmesos_1_1MatrixAdapter.html</anchorfile>
      <anchor>ad17a6004ba9322d05cb84067e15a3016</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Amesos::TypeMap</name>
    <filename>structAmesos_1_1TypeMap.html</filename>
    <templarg>ConcreteSolver</templarg>
    <templarg></templarg>
  </compound>
  <compound kind="struct">
    <name>Amesos::Util::has_special_impl</name>
    <filename>structAmesos_1_1Util_1_1has__special__impl.html</filename>
  </compound>
  <compound kind="struct">
    <name>Amesos::Util::row_access</name>
    <filename>structAmesos_1_1Util_1_1row__access.html</filename>
  </compound>
  <compound kind="struct">
    <name>Amesos::Util::col_access</name>
    <filename>structAmesos_1_1Util_1_1col__access.html</filename>
  </compound>
  <compound kind="struct">
    <name>Amesos::Util::get_cxs_helper</name>
    <filename>structAmesos_1_1Util_1_1get__cxs__helper.html</filename>
    <templarg>Matrix</templarg>
    <templarg>S</templarg>
    <templarg>GO</templarg>
    <templarg>GS</templarg>
    <templarg>Op</templarg>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>a5b13e00632fce03633bdff1db8db5cfa</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>a175f9f8006d3c767a225ca1881f4bbcf</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>af8f35e2282d630a52668a16d7136ee94</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; typename Matrix::local_ordinal_t, typename Matrix::global_ordinal_t, typename Matrix::node_t &gt; &gt; map, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Amesos::Util::get_ccs_helper</name>
    <filename>structAmesos_1_1Util_1_1get__ccs__helper.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <base>get_cxs_helper&lt; Matrix, S, GO, GS, get_ccs_func&lt; Matrix &gt; &gt;</base>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>a5b13e00632fce03633bdff1db8db5cfa</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>a175f9f8006d3c767a225ca1881f4bbcf</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>af8f35e2282d630a52668a16d7136ee94</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; typename Matrix::local_ordinal_t, typename Matrix::global_ordinal_t, typename Matrix::node_t &gt; &gt; map, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
    <docanchor file="structAmesos_1_1Util_1_1get__ccs__helper">get_ccs_helper_example</docanchor>
  </compound>
  <compound kind="struct">
    <name>Amesos::Util::get_crs_helper</name>
    <filename>structAmesos_1_1Util_1_1get__crs__helper.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <base>get_cxs_helper&lt; Matrix, S, GO, GS, get_crs_func&lt; Matrix &gt; &gt;</base>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>a5b13e00632fce03633bdff1db8db5cfa</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, EDistribution distribution, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>a175f9f8006d3c767a225ca1881f4bbcf</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>do_get</name>
      <anchorfile>structAmesos_1_1Util_1_1get__cxs__helper.html</anchorfile>
      <anchor>af8f35e2282d630a52668a16d7136ee94</anchor>
      <arglist>(const Teuchos::Ptr&lt; Matrix &gt; mat, const ArrayView&lt; S &gt; nzvals, const ArrayView&lt; GO &gt; indices, const ArrayView&lt; GS &gt; pointers, GS &amp;nnz, const Teuchos::Ptr&lt; const Tpetra::Map&lt; typename Matrix::local_ordinal_t, typename Matrix::global_ordinal_t, typename Matrix::node_t &gt; &gt; map, EStorage_Ordering ordering=Arbitrary)</arglist>
    </member>
  </compound>
  <compound kind="dir">
    <name>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</name>
    <path>/home/etbavie/dev/Trilinos/preCopyrightTrilinos/amesos2/src/</path>
    <filename>dir_a9f82e0d0a1fd3fb7ca2436fb62c836b.html</filename>
    <file>Amesos2.hpp</file>
    <file>Amesos2_AbstractConcreteMatrixAdapter.hpp</file>
    <file>Amesos2_ConcreteMatrixAdapter.hpp</file>
    <file>Amesos2_Control.cpp</file>
    <file>Amesos2_Control.hpp</file>
    <file>Amesos2_EpetraCrsMatrix_MatrixAdapter.cpp</file>
    <file>Amesos2_EpetraCrsMatrix_MatrixAdapter.hpp</file>
    <file>Amesos2_EpetraCrsMatrix_MatrixAdapter_decl.hpp</file>
    <file>Amesos2_EpetraCrsMatrix_MatrixAdapter_def.hpp</file>
    <file>Amesos2_EpetraMultiVecAdapter.cpp</file>
    <file>Amesos2_EpetraMultiVecAdapter.hpp</file>
    <file>Amesos2_EpetraMultiVecAdapter_decl.hpp</file>
    <file>Amesos2_EpetraMultiVecAdapter_def.hpp</file>
    <file>Amesos2_EpetraRowMatrix_AbstractMatrixAdapter.cpp</file>
    <file>Amesos2_EpetraRowMatrix_AbstractMatrixAdapter.hpp</file>
    <file>Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_decl.hpp</file>
    <file>Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_def.hpp</file>
    <file>Amesos2_Factory.cpp</file>
    <file>Amesos2_Factory.hpp</file>
    <file>Amesos2_Factory_decl.hpp</file>
    <file>Amesos2_Factory_def.hpp</file>
    <file>Amesos2_FunctionMap.hpp</file>
    <file>Amesos2_MatrixAdapter.cpp</file>
    <file>Amesos2_MatrixAdapter.hpp</file>
    <file>Amesos2_MatrixAdapter_decl.hpp</file>
    <file>Amesos2_MatrixAdapter_def.hpp</file>
    <file>Amesos2_MatrixHelper.hpp</file>
    <file>Amesos2_MatrixTraits.hpp</file>
    <file>Amesos2_MultiVecAdapter.cpp</file>
    <file>Amesos2_MultiVecAdapter.hpp</file>
    <file>Amesos2_MultiVecAdapter_decl.hpp</file>
    <file>Amesos2_Solver.cpp</file>
    <file>Amesos2_Solver.hpp</file>
    <file>Amesos2_Solver_decl.hpp</file>
    <file>Amesos2_SolverCore.cpp</file>
    <file>Amesos2_SolverCore.hpp</file>
    <file>Amesos2_SolverCore_decl.hpp</file>
    <file>Amesos2_SolverCore_def.hpp</file>
    <file>Amesos2_Status.cpp</file>
    <file>Amesos2_Status.hpp</file>
    <file>Amesos2_Superlu.cpp</file>
    <file>Amesos2_Superlu.hpp</file>
    <file>Amesos2_Superlu_decl.hpp</file>
    <file>Amesos2_Superlu_def.hpp</file>
    <file>Amesos2_Superlu_FunctionMap.hpp</file>
    <file>Amesos2_Superlu_MatrixHelper.hpp</file>
    <file>Amesos2_Superlu_TypeMap.hpp</file>
    <file>Amesos2_Superludist.cpp</file>
    <file>Amesos2_Superludist.hpp</file>
    <file>Amesos2_Superludist_decl.hpp</file>
    <file>Amesos2_Superludist_def.hpp</file>
    <file>Amesos2_Superludist_FunctionMap.hpp</file>
    <file>Amesos2_Superludist_MatrixHelper.hpp</file>
    <file>Amesos2_Superludist_TypeMap.hpp</file>
    <file>Amesos2_Superlumt.cpp</file>
    <file>Amesos2_Superlumt.hpp</file>
    <file>Amesos2_Superlumt_decl.hpp</file>
    <file>Amesos2_Superlumt_def.hpp</file>
    <file>Amesos2_Superlumt_FunctionMap.hpp</file>
    <file>Amesos2_Superlumt_MatrixHelper.hpp</file>
    <file>Amesos2_Superlumt_TypeMap.hpp</file>
    <file>Amesos2_Timers.hpp</file>
    <file>Amesos2_TpetraCrsMatrix_MatrixAdapter.cpp</file>
    <file>Amesos2_TpetraCrsMatrix_MatrixAdapter.hpp</file>
    <file>Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp</file>
    <file>Amesos2_TpetraCrsMatrix_MatrixAdapter_def.hpp</file>
    <file>Amesos2_TpetraMultiVecAdapter.cpp</file>
    <file>Amesos2_TpetraMultiVecAdapter.hpp</file>
    <file>Amesos2_TpetraMultiVecAdapter_decl.hpp</file>
    <file>Amesos2_TpetraMultiVecAdapter_def.hpp</file>
    <file>Amesos2_TpetraRowMatrix_AbstractMatrixAdapter.cpp</file>
    <file>Amesos2_TpetraRowMatrix_AbstractMatrixAdapter.hpp</file>
    <file>Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_decl.hpp</file>
    <file>Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_def.hpp</file>
    <file>Amesos2_TypeMap.hpp</file>
    <file>Amesos2_Util.cpp</file>
    <file>Amesos2_Util.hpp</file>
    <file>Amesos2_Util_decl.hpp</file>
    <file>Amesos2_Util_def.hpp</file>
    <file>Amesos2_Util_is_same.hpp</file>
    <file>Amesos2_Version.cpp</file>
    <file>Amesos2_Version.hpp</file>
  </compound>
</tagfile>
