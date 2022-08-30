BATCHED -- KokkosKernels batched functor-level interfaces
=========================================================

innerlu
-------
CodeCleanup-TODO: Move Decl file to dense/impl/KokkosBatched_InnerLU_Internal.hpp

applypivot
----------
.. doxygenstruct:: KokkosBatched::TeamVectorApplyPivot
    :members:

qr_withcolumnpivoting
---------------------
.. doxygenstruct:: KokkosBatched::TeamVectorQR_WithColumnPivoting
    :members:

addradial
---------
.. doxygenstruct:: KokkosBatched::SerialAddRadial
    :members:
.. doxygenstruct:: KokkosBatched::TeamAddRadial
    :members:

householder
-----------
.. doxygenstruct:: KokkosBatched::SerialHouseholder
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorHouseholder
    :members:

set
---
.. doxygenstruct:: KokkosBatched::SerialSet
    :members:
.. doxygenstruct:: KokkosBatched::TeamSet
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorSet
    :members:

scale
-----
.. doxygenstruct:: KokkosBatched::SerialScale
    :members:
.. doxygenstruct:: KokkosBatched::TeamScale
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorScale
    :members:

setidentity
-----------
.. doxygenstruct:: KokkosBatched::SerialSetIdentity
    :members:
.. doxygenstruct:: KokkosBatched::TeamSetIdentity
    :members:
.. doxygenstruct:: KokkosBatched::SetIdentity
    :members:

applyhouseholder
----------------
.. doxygenstruct:: KokkosBatched::SerialApplyHouseholder
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorApplyHouseholder
    :members:

innermultipledotproduct
-----------------------
CodeCleanup-TODO: Move Decl file to dense/impl/KokkosBatched_InnerMultipleDotProduct_Internal.hpp

lu
--
.. doxygenstruct:: KokkosBatched::SerialLU
    :members:
.. doxygenstruct:: KokkosBatched::TeamLU
    :members:
.. doxygenstruct:: KokkosBatched::LU
    :members:

solveutv
--------
.. doxygenstruct:: KokkosBatched::TeamVectorSolveUTV
    :members:

utv
---
.. doxygenstruct:: KokkosBatched::TeamVectorUTV
    :members:

inverselu
---------
CodeCleanup-TODO: Move Decl file to dense/impl/KokkosBatched_InverseLU_Internal.hpp

svd
---
.. doxygenstruct:: KokkosBatched::SerialSVD
    :members:

eigendecomposition
------------------
.. doxygenstruct:: KokkosBatched::SerialEigendecomposition
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorEigendecomposition
    :members:

trtri
-----
.. doxygenstruct:: KokkosBatched::SerialTrtri
    :members:

qr
--
.. doxygenstruct:: KokkosBatched::SerialQR
    :members:
.. doxygenstruct:: KokkosBatched::TeamQR
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorQR
    :members:
.. doxygenstruct:: KokkosBatched::QR
    :members:

trmm
----
.. doxygenstruct:: KokkosBatched::SerialTrmm
    :members:

trsm
----
.. doxygenstruct:: KokkosBatched::SerialTrsm
    :members:
.. doxygenstruct:: KokkosBatched::TeamTrsm
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorTrsm
    :members:
.. doxygenstruct:: KokkosBatched::Trsm
    :members:

innergemmfixa
-------------
CodeCleanup-TODO: Move Decl file to dense/impl/KokkosBatched_InnerGemmFixA_Internal.hpp

innergemmfixb
-------------
CodeCleanup-TODO: Move Decl file to dense/impl/KokkosBatched_InnerGemmFixB_Internal.hpp

innergemmfixc
-------------
CodeCleanup-TODO: Move Decl file to dense/impl/KokkosBatched_InnerGemmFixC_Internal.hpp

applyq
------
.. doxygenstruct:: KokkosBatched::SerialApplyQ
    :members:
.. doxygenstruct:: KokkosBatched::TeamApplyQ
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorApplyQ
    :members:
.. doxygenstruct:: KokkosBatched::ApplyQ
    :members:

copy
----
.. doxygenstruct:: KokkosBatched::SerialCopy
    :members:
.. doxygenstruct:: KokkosBatched::TeamCopy
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorCopy
    :members:
.. doxygenstruct:: KokkosBatched::Copy
    :members:

innertrsm
---------
CodeCleanup-TODO: Move Decl file to dense/impl/KokkosBatched_InnerTrsm_Internal.hpp

solvelu
-------
.. doxygenstruct:: KokkosBatched::SerialSolveLU
    :members:
.. doxygenstruct:: KokkosBatched::TeamSolveLU
    :members:
.. doxygenstruct:: KokkosBatched::SolveLU
    :members:

xpay
----
.. doxygenstruct:: KokkosBatched::SerialXpay
    :members:
.. doxygenstruct:: KokkosBatched::TeamXpay
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorXpay
    :members:

axpy
----
.. doxygenstruct:: KokkosBatched::SerialAxpy
    :members:
.. doxygenstruct:: KokkosBatched::TeamAxpy
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorAxpy
    :members:

gemv
----
.. doxygenstruct:: KokkosBatched::SerialGemv
    :members:
.. doxygenstruct:: KokkosBatched::TeamGemv
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorGemv
    :members:
.. doxygenstruct:: KokkosBatched::Gemv
    :members:

dot
---
.. doxygenstruct:: KokkosBatched::SerialDot
    :members:
.. doxygenstruct:: KokkosBatched::TeamDot
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorDot
    :members:

hadamardproduct
---------------
.. doxygenstruct:: KokkosBatched::SerialHadamardProduct
    :members:
.. doxygenstruct:: KokkosBatched::TeamHadamardProduct
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorHadamardProduct
    :members:
.. doxygenstruct:: KokkosBatched::HadamardProduct
    :members:

vector
------
CodeCleanup-TODO: Move Decl file to dense/impl/

trsv
----
.. doxygenstruct:: KokkosBatched::SerialTrsv
    :members:
.. doxygenstruct:: KokkosBatched::TeamTrsv
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorTrsv
    :members:
.. doxygenstruct:: KokkosBatched::Trsv
    :members:

gemm
----
.. doxygenstruct:: KokkosBatched::SerialGemm
    :members:
.. doxygenstruct:: KokkosBatched::TeamGemm
    :members:
.. doxygenstruct:: KokkosBatched::TeamVectorGemm
    :members:
.. doxygenstruct:: KokkosBatched::Gemm
    :members: