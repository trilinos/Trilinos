#ifndef _IFPACK_CrsGRAPH_H_
#define _IFPACK_CrsGRAPH_H_

class Ifpack_CrsGraph
{
public:
    virtual ~Ifpack_CrsGraph() {}

    virtual int NumRows() const = 0;
    virtual int NumCols() const = 0;
    virtual int IndexBase() const = 0;
    virtual int NumIndices(int Row) const = 0;
    virtual int * ExtractRowCopy(int Row, int LenOfIndices, int & NumIndices, int *& Indices) const = 0;

};

#endif // _IFPACK_CrsGRAPH_H_
