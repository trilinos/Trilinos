#ifndef INC_MemoryMgmt_h
#define INC_MemoryMgmt_h

namespace GAASP {

//-----------------------------------------------------------------------------
//	MemoryMgmt.h
//	This is temporary until arrays are made "smart".
//-----------------------------------------------------------------------------

template <class T> T* CreateArray (int const dim)
{
    return new T [dim];		// not initialized!
}

template <class T> T** CreateArray (int const numRows, int const numCols)
{
    T** p = CreateArray<T*> (numRows);
    T** pp = p;
    for ( int i = 0; i < numRows; ++i, ++pp )
	(*pp) = CreateArray<T> (numCols);
    return p;
}

template <class T> void DeleteArray ( T* & array1D )
{
    delete [] array1D;
    array1D = 0;
}

template <class T> void DeleteArray ( T** & array2D, int const numRows )
{
    if ( array2D )
    {
	for ( int i = 0; i < numRows; ++i )
	    delete [] array2D[ i ];
	delete [] array2D;
	array2D = 0;
    }
}

} // namespace GAASP

#endif // INC_MemoryMgmt_h
