/*
  Direct translation of parts of Galeri to use Cthulhu rather than Epetra.
*/

#ifndef __MATRIX_TYPES_HPP__
#define  __MATRIX_TYPES_HPP__

#include "Cthulhu_Map.hpp"
#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_CrsMatrixFactory.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ArrayView.hpp"

/* prototypes */
template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(const GlobalOrdinal i,
                              const GlobalOrdinal nx, const GlobalOrdinal ny,
                              GlobalOrdinal & left, GlobalOrdinal & right, 
                              GlobalOrdinal & lower, GlobalOrdinal & upper);

template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                        GlobalOrdinal& left, GlobalOrdinal& right, GlobalOrdinal& lower, GlobalOrdinal& upper,
                        GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2);

template <typename GlobalOrdinal>
void GetNeighboursCartesian3d(const GlobalOrdinal i,
                              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                              GlobalOrdinal& left, GlobalOrdinal& right,
                              GlobalOrdinal& lower, GlobalOrdinal& upper,
                              GlobalOrdinal& below, GlobalOrdinal& above);

template <typename GlobalOrdinal, typename Scalar>
void Fill9PointStencil(const GlobalOrdinal center,
                       std::vector<Scalar>  &Values,
                       std::vector<GlobalOrdinal>  &Indices,
                       GlobalOrdinal &NumEntries,
                       const GlobalOrdinal nx, const GlobalOrdinal ny,
                       const GlobalOrdinal nz,
                       const Scalar b,  const Scalar c,  const Scalar d,  const Scalar e,
                       const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                       GlobalOrdinal left=-2, GlobalOrdinal right=-2,
                       GlobalOrdinal lower=-2, GlobalOrdinal upper=-2);

/* end of prototypes */

/* ****************************************************************************************************** *
 *    Laplace 1D
 * ****************************************************************************************************** */
template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node>
Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TriDiag(Teuchos::RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > Map,
        const GlobalOrdinal nx,
        const Scalar a, const Scalar b, const Scalar c)
{
  Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix;

  Matrix = Cthulhu::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Map,  3);

  LocalOrdinal NumMyElements = Map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = Map->getNodeElementList();

  GlobalOrdinal NumGlobalElements = Map->getGlobalNumElements();

  LocalOrdinal NumEntries;
  LocalOrdinal nnz=2;
  std::vector<Scalar> Values(nnz);
  std::vector<LocalOrdinal> Indices(nnz);

  for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i)
  {
    if (MyGlobalElements[i] == 0)
    {
      // off-diagonal for first row
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = c;
    }
    else if (MyGlobalElements[i] == NumGlobalElements - 1)
    {
      // off-diagonal for last row
      Indices[0] = NumGlobalElements - 2;
      NumEntries = 1;
      Values[0] = b;
    }
    else
    {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1] = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0] = c;
      NumEntries = 2;
    }

    // put the off-diagonal entries
    // Cthulhu wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
    Matrix->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    Matrix->insertGlobalValues(MyGlobalElements[i],
                               Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                               Teuchos::tuple<Scalar>(a) );
  } //for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i)

  Matrix->fillComplete();

  return(Matrix);
} //TriDiag

/* ****************************************************************************************************** *
 *    Laplace 2D
 * ****************************************************************************************************** */
template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node>
Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
Cross2D(Teuchos::RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > Map,
        const GlobalOrdinal nx, const GlobalOrdinal ny,
        const Scalar a, const Scalar b, const Scalar c, 
        const Scalar d, const Scalar e)
{
  Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix;

  Matrix = Cthulhu::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Map,  5);

  LocalOrdinal NumMyElements = Map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = Map->getNodeElementList();

  LocalOrdinal left, right, lower, upper;
  LocalOrdinal nnz=4;
  std::vector<Scalar> Values(nnz);
  std::vector<LocalOrdinal> Indices(nnz);

  //    e
  //  b a c
  //    d
  for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i) 
  {
    LocalOrdinal NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
                 left, right, lower, upper);

    if (left != -1) 
    {
      Indices[NumEntries] = left;
      Values[NumEntries] = b;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Indices[NumEntries] = right;
      Values[NumEntries] = c;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e;
      ++NumEntries;
    }
    // put the off-diagonal entries
    // Cthulhu wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
    Matrix->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    Matrix->insertGlobalValues(MyGlobalElements[i],
                               Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                               Teuchos::tuple<Scalar>(a) );
  }
  Matrix->fillComplete();

  return(Matrix);
} //Cross2D

/* ****************************************************************************************************** *
 *    Star2D
 * ****************************************************************************************************** */
template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node>
Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
Star2D(Teuchos::RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > Map,
        const GlobalOrdinal nx, const GlobalOrdinal ny,
        const Scalar a, const Scalar b, const Scalar c, 
        const Scalar d, const Scalar e,
        const Scalar z1, const Scalar z2,
        const Scalar z3, const Scalar z4)
{
  Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix;

  Matrix = Cthulhu::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Map,  9);

  LocalOrdinal NumMyElements = Map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = Map->getNodeElementList();

  GlobalOrdinal left, right, lower, upper;
  //Scalar Values[9];
  //GlobalOrdinal    Indices[9];
  std::vector<Scalar> Values(9);
  std::vector<GlobalOrdinal>    Indices(9);

  //  z3  e  z4
  //   b  a  c
  //  z1  d  z2
  for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i) 
  {
    GlobalOrdinal NumEntries = 0;
    /*
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
                 left, right, lower, upper);
    */
#define OLD_CODE
#ifndef OLD_CODE
    Fill9PointStencil<GlobalOrdinal,Scalar>(MyGlobalElements[i], Values, Indices, NumEntries, nx, ny,
                   b, c, d, e, z1, z2, z3, z4);
#else
    if (left != -1) 
    {
      Values[NumEntries] = b;
      Indices[NumEntries] = left;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Values[NumEntries] = c;
      Indices[NumEntries] = right;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Values[NumEntries] = d;
      Indices[NumEntries] = lower;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Values[NumEntries] = e;
      Indices[NumEntries] = upper;
      ++NumEntries;
    }
    if (left != -1 && lower != -1) 
    {
      Values[NumEntries] = z1;
      Indices[NumEntries] = lower - 1;
      ++NumEntries;
    }
    if (right != -1 && lower != -1) 
    {
      Values[NumEntries] = z2;
      Indices[NumEntries] = lower + 1;
      ++NumEntries;
    }
    if (left != -1 && upper != -1) 
    {
      Values[NumEntries] = z3;
      Indices[NumEntries] = upper - 1;
      ++NumEntries;
    }
    if (right != -1 && upper != -1) 
    {
      Values[NumEntries] = z4;
      Indices[NumEntries] = upper + 1;
      ++NumEntries;
    }
#endif
    
    Values[NumEntries] = a;
    Indices[NumEntries] = MyGlobalElements[i];
    ++NumEntries;

    Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
    Matrix->insertGlobalValues(MyGlobalElements[i], iv, av);
  }

  Matrix->fillComplete();

  return(Matrix);
} //Star2D

/* ****************************************************************************************************** *
 *    BigStar2D (2D Biharmonic operator, for example)
 * ****************************************************************************************************** */
template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node>
Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
BigStar2D(Teuchos::RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > Map,
          const GlobalOrdinal nx, const GlobalOrdinal ny,
          const Scalar a, const Scalar b, const Scalar c, 
          const Scalar d, const Scalar e,
          const Scalar z1, const Scalar z2,
          const Scalar z3, const Scalar z4,
          const Scalar bb, const Scalar cc, const Scalar dd, const Scalar ee)
{
  Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix;

  Matrix = Cthulhu::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Map,  13);

  LocalOrdinal NumMyElements = Map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = Map->getNodeElementList();

  GlobalOrdinal left, right, lower, upper;
  GlobalOrdinal left2, right2, lower2, upper2;
  Scalar Values[13];
  LocalOrdinal Indices[13];

  //        ee
  //    z3  e  z4
  // bb  b  a  c  cc
  //    z1  d  z2
  //        dd
  
  for (GlobalOrdinal i = 0 ; i < NumMyElements ; ++i) 
  {
    GlobalOrdinal NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, left, right, lower, upper,
                             left2, right2, lower2, upper2);

    if (left != -1) 
    {
      Values[NumEntries] = b;
      Indices[NumEntries] = left;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Values[NumEntries] = c;
      Indices[NumEntries] = right;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Values[NumEntries] = d;
      Indices[NumEntries] = lower;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Values[NumEntries] = e;
      Indices[NumEntries] = upper;
      ++NumEntries;
    }
    if (left != -1 && lower != -1) 
    {
      Values[NumEntries] = z1;
      Indices[NumEntries] = lower - 1;
      ++NumEntries;
    }
    if (right != -1 && lower != -1) 
    {
      Values[NumEntries] = z2;
      Indices[NumEntries] = lower + 1;
      ++NumEntries;
    }
    if (left != -1 && upper != -1) 
    {
      Values[NumEntries] = z3;
      Indices[NumEntries] = upper - 1;
      ++NumEntries;
    }
    if (right != -1 && upper != -1) 
    {
      Values[NumEntries] = z4;
      Indices[NumEntries] = upper + 1;
      ++NumEntries;
    }
    if (left2 != -1) 
    {
      Values[NumEntries] = bb;
      Indices[NumEntries] = left2;
      ++NumEntries;
    }
    if (right2 != -1) 
    {
      Values[NumEntries] = cc;
      Indices[NumEntries] = right2;
      ++NumEntries;
    }
    if (lower2 != -1) 
    {
      Values[NumEntries] = dd;
      Indices[NumEntries] = lower2;
      ++NumEntries;
    }
    if (upper2 != -1) 
    {
      Values[NumEntries] = ee;
      Indices[NumEntries] = upper2;
      ++NumEntries;
    }
    
    Values[NumEntries] = a;
    Indices[NumEntries] = MyGlobalElements[i];
    ++NumEntries;

    Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
    Matrix->insertGlobalValues(MyGlobalElements[i], iv, av);
  }

  Matrix->fillComplete();

  return(Matrix);
} //BigStar2D

/* ****************************************************************************************************** *
 *    Laplace3D
 * ****************************************************************************************************** */
template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node>
Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
Cross3D(Teuchos::RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > Map,
        const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
        const Scalar a, const Scalar b, const Scalar c, 
        const Scalar d, const Scalar e,
        const Scalar f, const Scalar g)
{
  Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix;

  Matrix = Cthulhu::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Map, 7);

  LocalOrdinal NumMyElements = Map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = Map->getNodeElementList();

  GlobalOrdinal left, right, lower, upper, below, above;
  std::vector<Scalar> Values(6);
  std::vector<GlobalOrdinal> Indices(6);

  //    e
  //  b a c
  //    d
  // + f below and g above
  
  for (GlobalOrdinal i = 0 ; i < NumMyElements ; ++i) 
  {
    GlobalOrdinal NumEntries = 0;
    GetNeighboursCartesian3d(MyGlobalElements[i], nx, ny, nz,
                 left, right, lower, upper, below, above);

    if (left != -1) 
    {
      Indices[NumEntries] = left;
      Values[NumEntries] = b;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Indices[NumEntries] = right;
      Values[NumEntries] = c;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e;
      ++NumEntries;
    }
    if (below != -1) 
    {
      Indices[NumEntries] = below;
      Values[NumEntries] = f;
      ++NumEntries;
    }
    if (above != -1) 
    {
      Indices[NumEntries] = above;
      Values[NumEntries] = g;
      ++NumEntries;
    }
    // put the off-diagonal entries
    Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
    Matrix->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    Matrix->insertGlobalValues(MyGlobalElements[i],
                               Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                               Teuchos::tuple<Scalar>(a) );
  }
  Matrix->fillComplete();

  return(Matrix);
} //Cross3D

/* ****************************************************************************************************** *
 *    3D 27-point stencil
 * ****************************************************************************************************** */
template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node>
Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
Brick3D(Teuchos::RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> > Map,
        const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
        const Scalar a, const Scalar b, const Scalar c, 
        const Scalar d, const Scalar e,
        const Scalar f, const Scalar g)
{
  Teuchos::RCP< Cthulhu::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix;

  Matrix = Cthulhu::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Map, 27);

  LocalOrdinal NumMyElements = Map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = Map->getNodeElementList();

  GlobalOrdinal left, right, lower, upper, below, above;
  std::vector<Scalar> Values(26);
  std::vector<GlobalOrdinal> Indices(26);

  // upper plane
  //   e  d  e
  //   d  b  d
  //   e  d  e

  // middle plane
  //   c  b  c
  //   b  a  b
  //   c  b  c

  // lower plane
  //   e  d  e
  //   d  b  d
  //   e  d  e

  for (GlobalOrdinal i = 0 ; i < NumMyElements ; ++i) 
  {
    GlobalOrdinal NumEntries = 0;
    // First process the middle plane and figure out the
    // centers of the upper and lower planes
    GetNeighboursCartesian3d(MyGlobalElements[i], nx, ny, nz,
                 left, right, lower, upper, below, above);

#ifdef DEBUG_STENCIL
    std::cout << "middle plane, center " << MyGlobalElements[i]  << std::endl;
#endif
    Fill9PointStencil(MyGlobalElements[i],Values, Indices,
                 NumEntries, nx, ny, nz,
                 b, b, b, b, c, c, c, c,
                 left, right, lower, upper);

    //process lower plane, centered on "below"
    if (below != -1) 
    {
      Indices[NumEntries] = below;
      Values[NumEntries] = b;
      ++NumEntries;
#ifdef DEBUG_STENCIL
    std::cout << "lower plane, center " << below  << std::endl;
#endif
      Fill9PointStencil(below, Values, Indices, NumEntries, nx, ny, nz,
                   d, d, d, d, e, e, e, e);
    }
    //process upper plane, centered on "upper"
    if (above != -1) 
    {
      Indices[NumEntries] = above;
      Values[NumEntries] = b;
      ++NumEntries;
#ifdef DEBUG_STENCIL
    std::cout << "upper plane, center " << above  << std::endl;
#endif
      Fill9PointStencil(above, Values, Indices, NumEntries, nx, ny, nz,
                   d, d, d, d, e, e, e, e);
    }


    // put the off-diagonal entries
    Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
    Matrix->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    Matrix->insertGlobalValues(MyGlobalElements[i],
                               Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                               Teuchos::tuple<Scalar>(a) );
  } //for (GlobalOrdinal i = 0 ; i < NumMyElements ; ++i) 
  Matrix->fillComplete();

  return(Matrix);
} //Brick3D

/* ****************************************************************************************************** *
 *    Utilities
 * ****************************************************************************************************** */

/* GetNeighboursCartesian2d */
template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                              GlobalOrdinal & left, GlobalOrdinal & right,
                              GlobalOrdinal & lower, GlobalOrdinal & upper)
{
  GlobalOrdinal ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)      left = -1;
  else              left = i - 1;
  if (ix == nx - 1) right = -1;
  else              right = i + 1;
  if (iy == 0)      lower = -1;
  else              lower = i - nx;
  if (iy == ny - 1) upper = -1;
  else              upper = i + nx;
}

/* GetNeighboursCartesian2d */
template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                        GlobalOrdinal& left, GlobalOrdinal& right, GlobalOrdinal& lower, GlobalOrdinal& upper,
                        GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2)
{
  GlobalOrdinal ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)      left = -1;
  else              left = i - 1;
  if (ix == nx - 1) right = -1;
  else              right = i + 1;
  if (iy == 0)      lower = -1;
  else              lower = i - nx;
  if (iy == ny - 1) upper = -1;
  else              upper = i + nx;

  if (ix <= 1)      left2 = -1;
  else              left2 = i - 2;
  if (ix >= nx - 2) right2 = -1;
  else              right2 = i + 2;
  if (iy <= 1)      lower2 = -1;
  else              lower2 = i - 2 * nx;
  if (iy >= ny - 2) upper2 = -1;
  else              upper2 = i + 2 * nx;
}

/* GetNeighboursCartesian3d */
template <typename GlobalOrdinal>
void GetNeighboursCartesian3d(const GlobalOrdinal i,
                              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                              GlobalOrdinal& left, GlobalOrdinal& right,
                              GlobalOrdinal& lower, GlobalOrdinal& upper,
                              GlobalOrdinal& below, GlobalOrdinal& above)
{
  GlobalOrdinal ixy, iz;
  ixy = i % (nx * ny);

  iz = (i - ixy) / (nx * ny); 

  if (iz == 0)      below = -1;
  else              below = i - nx * ny;
  if (iz == nz - 1) above = -1;
  else              above = i + nx * ny;

  GetNeighboursCartesian2d(ixy, nx, ny, left, right, lower, upper);

  if (left  != -1) left  += iz * (nx * ny);
  if (right != -1) right += iz * (nx * ny);
  if (lower != -1) lower += iz * (nx * ny);
  if (upper != -1) upper += iz * (nx * ny);
}

/* Fill9PointStencil */
template <typename GlobalOrdinal, typename Scalar>
void Fill9PointStencil(const GlobalOrdinal center,
                       std::vector<Scalar>  &Values,
                       std::vector<GlobalOrdinal>  &Indices,
                       GlobalOrdinal &NumEntries,
                       const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                       const Scalar b,  const Scalar c,  const Scalar d,  const Scalar e,
                       const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                       GlobalOrdinal left, GlobalOrdinal right,
                       GlobalOrdinal lower, GlobalOrdinal upper)
{
  //  z3  e  z4
  //   b  .  c
  //  z1  d  z2
  GlobalOrdinal below,above;
#ifdef DEBUG_STENCIL
  std::cout << "   ";
#endif
  if (left == -2) {
    GetNeighboursCartesian3d(center, nx, ny, nz, left, right, lower, upper, below, above);
  }

  if (left != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << left << " ";
#endif
    Values[NumEntries] = b;
    Indices[NumEntries] = left;
    ++NumEntries;
  }
  if (right != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << right << " ";
#endif
    Values[NumEntries] = c;
    Indices[NumEntries] = right;
    ++NumEntries;
  }
  if (lower != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << lower << " ";
#endif
    Values[NumEntries] = d;
    Indices[NumEntries] = lower;
    ++NumEntries;
  }
  if (upper != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << upper << " ";
#endif
    Values[NumEntries] = e;
    Indices[NumEntries] = upper;
    ++NumEntries;
  }
  if (left != -1 && lower != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << lower-1 << " ";
#endif
    Values[NumEntries] = z1;
    Indices[NumEntries] = lower - 1;
    ++NumEntries;
  }
  if (right != -1 && lower != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << lower+1 << " ";
#endif
    Values[NumEntries] = z2;
    Indices[NumEntries] = lower + 1;
    ++NumEntries;
  }
  if (left != -1 && upper != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << upper-1 << " ";
#endif
    Values[NumEntries] = z2;
    Values[NumEntries] = z3;
    Indices[NumEntries] = upper - 1;
    ++NumEntries;
  }
  if (right != -1 && upper != -1) 
  {
#ifdef DEBUG_STENCIL
    std::cout << upper+1 << " ";
#endif
    Values[NumEntries] = z4;
    Indices[NumEntries] = upper + 1;
    ++NumEntries;
  }
#ifdef DEBUG_STENCIL
  std::cout << std::endl;
#endif
} //Fill9PointStencil()

#endif //ifndef __MATRIX_TYPES_HPP__
