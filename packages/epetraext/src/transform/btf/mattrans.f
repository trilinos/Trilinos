        subroutine mattrans (m, n, ja, ia, jao, iao)
        integer ia(m+1), iao(n+1), ja(*), jao(*)
c-------------------------------------------------------------------
c transpose a matrix stored in a, ja, ia format.
c ---------------
c input arguments:
c m     = row dimension of A.
c n     = column dimension of A.
c ja    = integer array of size nnz containing the column positions
c         of the corresponding elements in a.
c ia    = integer of size n+1. ia(k) contains the position in a, ja of
c         the beginning of the k-th row.
c output arguments:
c jao   = integer array of size nnz containing the column indices.
c iao   = integer array of size m+1 containing the "ia" index array of
c         the transposed matrix.
c--------------------------------------------------------------------
c
c  count the number of elements in every column of a and row of ao
c
        do 1 i=1, n+1
 1           iao(i) = 0
        do 3 i=1, m
                k1 = ia(i)
                k2 = ia(i+1) -1
                do 2 k=k1, k2
                        j = ja(k)+1
                        iao(j) = iao(j)+1
 2              continue
 3      continue
c find addresses of new first elements..
        iao(1) = 1
        do 4 i=1, n
 4      iao(i+1) = iao(i) + iao(i+1)
c now do the actual copying.
        do 6 i=1, m
                k1 = ia(i)
                k2 = ia(i+1)-1
                do 62 k=k1,k2
                        j = ja(k)
                        next = iao(j)
                        jao(next) = i
                        iao(j) = next+1
 62             continue
 6      continue
c reshift iao
        do 7 i = n, 1, -1
 7         iao(i+1) = iao(i)
        iao(1) = 1
c--------------- end of mattrans ---------------------------------
        end
