Karen:  will put Phipps' code into Zoltan2, get outline of test compiling
Put it all on a branch.
Erik:  will take over to make sure it works correctly.

typedef Tpetra::CrsMatrix<#> matrix_type;
typedef Tpetra::MultiVector<#> multivector_type;

Zoltan2::TpetraCrsColorer<matrix_type> colorer(J);
{

  if (Jacobian)
    if IanNotReady && !testZoltan2anyway
      ZoltanCrsColorer(PartialDistance2)
    else if IanReady || testZoltan2anyway
      Zoltan2CrsColorer(PartialDistance2)

  else if (Hessian)
    // Hessian is already symmetric; 
    // code currently still uses partial distance 2.
    // Should use Distance2 instead?
    if IanNotReady && !testZoltan2anyway
      ZoltanCrsColorer(Distance2)
    else if IanReady || testZoltan2anyway
      Zoltan2CrsColorer(Distance2)


}

main()
{

//KDD Tpetra::ZoltanCrsColorer<matrix_type> colorer(J);
Teuchos::ParameterList coloring_params;

coloring_params.set("MatrixType", "Jacobian");  // Valid values are Jacobian or Hessian

coloring_params.set("Symmetrize", true);  // We always do this in SPARC, although I don't remember why exactly

// Compute coloring
colorer.computeColoring(coloring_params);

// Check coloring
const bool valid_coloring = colorer.checkColoring();
if (!valid_coloring) { #}

// Compute seed matrix -- this is what the application wants
multivector_type V(#);  // Dense matrix of 0/1 indicating the compression via coloring
colorer.computeSeedMatrix(V);


// To test the result...
// Compute compression vector
multivector_type W(#);  // W is the compressed matrix
J->apply(W, V);

// Reconstruct matrix from compression vector
matrix_type Jp(#);
colorer.reconstructMatrix(W, Jp);

// Check J = Jp somehow
// KDD is there a way to do this comparison in Tpetra?

}

Test cases -- UserInputForTests can generate Galeri or read files:
-  tri-diagonal matrix -- can check the number of colors
-  galeri matrix

Want to test with fitted and non-fitted maps
Call regular and fitted versions of functions

Test both with and without Symmetrize -- test both to exercise both sets of 
callbacks in Zoltan

