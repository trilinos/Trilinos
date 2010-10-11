/*
  Direct translation of parts of Galeri matrix generator.
*/
#ifndef __MATRIX_FACTORY_HPP__
#define __MATRIX_FACTORY_HPP__

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "matrix_types.hpp"

#include <iostream>

template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
CreateCrsMatrix(const std::string &MatrixType,
                Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > Map,
                Teuchos::ParameterList& List)
{
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > returnMatrix;
  if (MatrixType == "Laplace1D") {

    GlobalOrdinal nx = List.get("nx", -1);
    if (nx == -1)
    {
      GlobalOrdinal n = Map->getGlobalNumElements();
      nx = (GlobalOrdinal)sqrt((Scalar)n);
      TEST_FOR_EXCEPTION(nx*nx != n, std::logic_error, "You need to specify nx.");
    }

    //return(TriDiag<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Map, nx, 2.0, -1.0, -1.0));
    returnMatrix = TriDiag<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Map, nx, 2.0, -1.0, -1.0);

  } else if (MatrixType == "Laplace2D") {

    GlobalOrdinal nx = List.get("nx", -1);
    GlobalOrdinal ny = List.get("ny", -1);
    if (nx == -1 || ny == -1)
    {
      GlobalOrdinal n = Map->getGlobalNumElements();
      nx = (GlobalOrdinal)sqrt((Scalar)n);
      ny = nx;
      TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
    }

    //return(Cross2D<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Map, nx, ny, 4.0, -1.0, -1.0, -1.0, -1.0));
    returnMatrix = Cross2D<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Map, nx, ny, 4.0, -1.0, -1.0, -1.0, -1.0);

  } else if (MatrixType == "Star2D") {

    GlobalOrdinal nx = List.get("nx", -1);
    GlobalOrdinal ny = List.get("ny", -1);

    Scalar a = List.get("a", 8.0);
    Scalar b = List.get("b", -1.0);
    Scalar c = List.get("c", -1.0);
    Scalar d = List.get("d", -1.0);
    Scalar e = List.get("e", -1.0);
    Scalar z1 = List.get("z1", -1.0);
    Scalar z2 = List.get("z2", -1.0);
    Scalar z3 = List.get("z3", -1.0);
    Scalar z4 = List.get("z4", -1.0);

    //return(Star2D(Map, nx, ny, a, b, c, d, e, z1, z2, z3, z4));
    returnMatrix = Star2D(Map, nx, ny, a, b, c, d, e, z1, z2, z3, z4);

  } else if (MatrixType == "BigStar2D") {

    GlobalOrdinal nx = List.get("nx", -1);
    GlobalOrdinal ny = List.get("ny", -1);

    Scalar a = List.get("a", 20.0);
    Scalar b = List.get("b", -8.0);
    Scalar c = List.get("c", -8.0);
    Scalar d = List.get("d", -8.0);
    Scalar e = List.get("e", -8.0);
    Scalar z1 = List.get("z1", 2.0);
    Scalar z2 = List.get("z2", 2.0);
    Scalar z3 = List.get("z3", 2.0);
    Scalar z4 = List.get("z4", 2.0);
    Scalar bb = List.get("bb", 1.0);
    Scalar cc = List.get("cc", 1.0);
    Scalar dd = List.get("dd", 1.0);
    Scalar ee = List.get("ee", 1.0);

    //return(BigStar2D(Map, nx, ny, a, b, c, d, e, z1, z2, z3, z4, bb, cc, dd, ee));
    returnMatrix = BigStar2D(Map, nx, ny, a, b, c, d, e, z1, z2, z3, z4, bb, cc, dd, ee);

  } else if (MatrixType == "Laplace3D") {

    GlobalOrdinal nx = List.get("nx", -1);
    GlobalOrdinal ny = List.get("ny", -1);
    GlobalOrdinal nz = List.get("nz", -1);
    if (nx == -1 || ny == -1 || nz == -1)
    {
      GlobalOrdinal n = Map->getGlobalNumElements();
      nx = (GlobalOrdinal)pow((Scalar)n, 0.33334);
      ny = nx; nz = nx;
      TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
    } 
    //return(Cross3D(Map, nx, ny, nz, 6.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0));
    returnMatrix = Cross3D(Map, nx, ny, nz, 6.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);

  } else if (MatrixType == "Brick3D") {

    GlobalOrdinal nx = List.get("nx", -1);
    GlobalOrdinal ny = List.get("ny", -1);
    GlobalOrdinal nz = List.get("nz", -1);
    if (nx == -1 || ny == -1 || nz == -1)
    {
      GlobalOrdinal n = Map->getGlobalNumElements();
      nx = (GlobalOrdinal)pow((Scalar)n, 0.33334);
      ny = nx; nz = nx;
      TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
    } 
    //return(Brick3D(Map, nx, ny, nz, 26.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0));
    returnMatrix = Brick3D(Map, nx, ny, nz, 26.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);

  } else {

     TEST_FOR_EXCEPTION(true,
                        std::logic_error,
                        "`MatrixType' has incorrect value (" << MatrixType
                        << ") in input to function CreateCrsMatrix()."
                        << "Check the documentation for a list of valid choices");
  } //if-else

  returnMatrix->setObjectLabel(MatrixType);
  return(returnMatrix);

} // CreateCrsMatrix()

#endif //ifndef __MATRIX_FACTORY_HPP__
