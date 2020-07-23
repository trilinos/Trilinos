#ifndef __TACHO_BLAS_TEAM_HPP__
#define __TACHO_BLAS_TEAM_HPP__

/// \file  Tacho_Blas_TEAM.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

namespace Tacho {

    template<typename T>
    struct BlasTeam {
      struct Impl {
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void set(MemberType &member,
                 int m, 
                 const T alpha, 
                 /* */ T *__restrict__ a, int as0) {
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member,m),[&](const int &i) {
              a[i*as0] = alpha;
            });
        }
        
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void scale(MemberType &member,
                   int m, 
                   const T alpha, 
                   /* */ T *__restrict__ a, int as0) {
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member,m),[&](const int &i) {
              a[i*as0] *= alpha;
            });
        }
        
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void set(MemberType &member,
                 int m, int n, 
                 const T alpha, 
                 /* */ T *__restrict__ a, int as0, int as1) {
          if (as0 == 1 || as0 < as1) 
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,m),[&](const int &i) {
                    a[i*as0+j*as1] = alpha;
                  });
              });
          else 
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,m),[&](const int &i) {
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,n),[&](const int &j) {
                    a[i*as0+j*as1] = alpha;
                  });
              });
        }
        
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void scale(MemberType &member,
                   int m, int n, 
                   const T alpha, 
                   /* */ T *__restrict__ a, int as0, int as1) {
          if (as0 == 1 || as0 < as1) 
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,m),[&](const int &i) {
                    a[i*as0+j*as1] *= alpha;
                  });
              });
          else 
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member,m),[&](const int &i) {
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,n),[&](const int &j) {
                    a[i*as0+j*as1] *= alpha;
                  });
              });
        }
        
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void set_upper(MemberType &member,
                       int m, int n, int offset,
                       const T alpha, 
                       /* */ T *__restrict__ a, int as0, int as1) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,j+1-offset),[&](const int &i) {
                  a[i*as0+j*as1] = alpha;
                });
            });
        }
        
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void scale_upper(MemberType &member,
                         int m, int n, int offset, 
                         const T alpha, 
                         /* */ T *__restrict__ a, int as0, int as1) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,j+1-offset),[&](const int &i) {
                  a[i*as0+j*as1] *= alpha;
                });
            });
        }

        
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void set_lower(MemberType &member,
                       int m, int n, int offset,
                       const T alpha, 
                       /* */ T *__restrict__ a, int as0, int as1) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
              const int jj = j + offset;
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,n-j-offset),[&](const int &i) {
                  a[(i+jj)*as0+j*as1] = alpha;
                });
            });
        }
        
        template<typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void scale_lower(MemberType &member,
                         int m, int n, int offset, 
                         const T alpha, 
                         /* */ T *__restrict__ a, int as0, int as1) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
              const int jj = j + offset;
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,n-j-offset),[&](const int &i) {
                  a[(i+jj)*as0+j*as1] *= alpha;
                });
            });
        }

        template<typename ConjType, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void gemv(MemberType &member, const ConjType &cj,
                  const int m, const int n, 
                  const T alpha, 
                  const T *__restrict__ A, const int as0, const int as1,
                  const T *__restrict__ x, const int xs0,
                  const T beta,
                  /* */ T *__restrict__ y, const int ys0) {          
          const T one(1), zero(0);

          if (beta == zero) set  (member, m, zero, y, ys0);
          if (beta != one ) scale(member, m, beta, y, ys0);
          
          if (alpha != zero) {
            if (m <=0 || n <=0) return;

            member.team_barrier();            
            {
              if (as0 == 1) { 
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,m),[&](const int &i) {
                    T t(0);
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
                        t += cj(A[i*as0+j*as1])*x[j*xs0];
                      });                  
                    Kokkos::atomic_add(&y[i*ys0], alpha*t);
                  });
              } else {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(member,m),[&](const int &i) {
                    T t(0);
                    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,n),[&](const int &j) {
                        t += cj(A[i*as0+j*as1])*x[j*xs0];
                      });                  
                    Kokkos::atomic_add(&y[i*ys0], alpha*t);
                  });
              }
            }
          }
        }

        template<typename ConjType, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void trsv_upper(MemberType &member, const ConjType &cjA,
                        const char diag, 
                        const int m, 
                        const T *__restrict__ A, const int as0, const int as1,
                        /* */ T *__restrict__ b, const int bs0) {          
          if (m <= 0) return;
          
          const bool use_unit_diag = diag == 'U'|| diag == 'u';
          T *__restrict__ b0 = b;
          for (int p=(m-1);p>=0;--p) {
            const int iend = p;
            
            const T *__restrict__ a01   = A+p*as1;
            /**/  T *__restrict__ beta1 = b+p*bs0;
            
            /// make sure the previous iteration update is done
            member.team_barrier();
            T local_beta1 = *beta1;
            if (!use_unit_diag) {
              const T alpha11 = cjA(A[p*as0+p*as1]);
              local_beta1 /= alpha11;
              /// before modifying beta1 we need make sure
              /// that every local_beta1 has the previous beta1 value
              member.team_barrier();
              Kokkos::single(Kokkos::PerTeam(member), [&]() {
                  *beta1 = local_beta1;
                });
            }
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member,iend),[&](const int &i) {
                b0[i*bs0] -= cjA(a01[i*as0]) * local_beta1;
              });
          }
        } 


        template<typename ConjType, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void trsv_lower(MemberType &member, const ConjType &cjA,
                        const char diag, 
                        const int m, 
                        const T *__restrict__ A, const int as0, const int as1,
                        /* */ T *__restrict__ b, const int bs0) {          
          if (m <= 0) return;

          const bool use_unit_diag = diag == 'U'|| diag == 'u';
          //T *__restrict__ b0 = b;
          for (int p=0;p<m;++p) {
            const int iend = m-p-1;
            
            const T
              *__restrict__ a21   = iend ? A+(p+1)*as0+p*as1 : NULL;
            
            T
              *__restrict__ beta1 =        b+p*bs0,
              *__restrict__ b2    = iend ? beta1+bs0 : NULL;
            
            /// make sure that the previous iteration update is done
            member.team_barrier();
            T local_beta1 = *beta1;
            if (!use_unit_diag) {
              const T alpha11 = A[p*as0+p*as1];
              local_beta1 /= alpha11;
              /// before modifying beta1 we need make sure
              /// that every local_beta1 has the previous beta1 value
              member.team_barrier();
              Kokkos::single(Kokkos::PerTeam(member), [&]() {              
                  *beta1 = local_beta1;
                });
            }
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member,iend),[&](const int &i) {
                b2[i*bs0] -= a21[i*as0] * local_beta1;
              });
          }
        }

       
        template<typename ConjTypeA, typename ConjTypeB, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void gemm(MemberType &member, const ConjTypeA &cjA, const ConjTypeB &cjB,
                  const int m, const int n, const int k,
                  const T alpha, 
                  const T *__restrict__ A, const int as0, const int as1,
                  const T *__restrict__ B, const int bs0, const int bs1,
                  const T beta,
                  /* */ T *__restrict__ C, const int cs0, const int cs1) {
          const T one(1), zero(0);
          
          if      (beta == zero) set  (member, m, n, zero, C, cs0, cs1);
          else if (beta != one ) scale(member, m, n, beta, C, cs0, cs1);
          
          if (alpha != zero) {
            if (m <= 0 || n <= 0 || k <= 0) return;
            
            member.team_barrier();
            {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,m),[&](const int &i) {
                      const T 
                        *__restrict__ pA = A+i*as0, 
                        *__restrict__ pB = B+j*bs1;
                      T c(0);
                      for (int p=0;p<k;++p)
                        c += cjA(pA[p*as1])*cjB(pB[p*bs0]);
                      C[i*cs0+j*cs1] += alpha*c;
                    });
                });
            }
          } 
        }


        template<typename ConjTypeA, typename ConjTypeB, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void herk_upper(MemberType &member, const ConjTypeA &cjA, const ConjTypeB &cjB,  
                        const int n, const int k,
                        const T alpha, 
                        const T *__restrict__ A, const int as0, const int as1,
                        const T beta,
                        /* */ T *__restrict__ C, const int cs0, const int cs1) {
          const T one(1), zero(0);
          
          if      (beta == zero) set_upper  (member, n, n, 0, zero, C, cs0, cs1);
          else if (beta != one ) scale_upper(member, n, n, 0, beta, C, cs0, cs1);
          
          if (alpha != zero) {
            if (n <= 0 || k <= 0) return;
            
            member.team_barrier();
            {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
                  const T *__restrict__ pA = A+j*as0;
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,j+1),[&](const int &i) {
                      const T *__restrict__ pB = A+i*as0;
                      T c(0);
                      for (int p=0;p<k;++p) 
                        c += cjA(pA[p*as1])*cjB(pB[p*as1]);
                      C[i*cs0+j*cs1] += alpha*c;
                    });
                });
            }
          } 
        }


        template<typename ConjTypeA, typename ConjTypeB, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void herk_lower(MemberType &member, const ConjTypeA &cjA, const ConjTypeB &cjB,
                        const int n, const int k,
                        const T alpha, 
                        const T *__restrict__ A, const int as0, const int as1,
                        const T beta,
                        /* */ T *__restrict__ C, const int cs0, const int cs1) {
          const T one(1), zero(0);
          
          if      (beta == zero) set_lower  (member, n, n, 0, zero, C, cs0, cs1);
          else if (beta != one ) scale_lower(member, n, n, 0, beta, C, cs0, cs1);
          
          if (alpha != zero) {
            if (n <= 0 || k <= 0) return;
            
            member.team_barrier();
            {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,n),[&](const int &j) {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,n-j),[&](const int &i) {
                      const int ii = i+j;
                      const T 
                        *__restrict__ pA = A+j*as0,
                        *__restrict__ pB = A+ii*as0;
                      T c(0);
                      for (int p=0;p<k;++p)
                        c += cjA(pA[p*as1])*cjB(pB[p*as1]);
                      C[ii*cs0+j*cs1] += alpha*c;
                    });
                });
            }
          } 
        }


        template<typename ConjType, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void trsm_left_lower(MemberType &member, const ConjType &cjA, 
                             const char diag, 
                             const int m, const int n, 
                             const T alpha, 
                             const T *__restrict__ A, const int as0, const int as1,
                             /* */ T *__restrict__ B, const int bs0, const int bs1) {          
          const T one(1), zero(0);

          if (alpha == zero)   set  (member, m, n, zero,  B, bs0, bs1);
          else {
            if (alpha != one)  scale(member, m, n, alpha, B, bs0, bs1);
            if (m <= 0 || n <= 0) return;

            const bool use_unit_diag = diag == 'U'|| diag == 'u';
            for (int p=0;p<m;++p) {
              const int iend = m-p-1, jend = n;

              const T
                *__restrict__ a21 = iend ? A+(p+1)*as0+p*as1 : NULL;
              
              T
                *__restrict__ b1t =        B+p*bs0,
                *__restrict__ B2  = iend ? B+(p+1)*bs0 : NULL;
              
              member.team_barrier();
              if (!use_unit_diag) {
                const T alpha11 = cjA(A[p*as0+p*as1]);
                Kokkos::parallel_for(Kokkos::TeamVectorRange(member,jend),[&](const int &j) {
                    b1t[j*bs1] /= alpha11;
                  });
                member.team_barrier();
              }

              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,jend),[&](const int &j) {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,iend),[&](const int &i) {
                      B2[i*bs0+j*bs1] -= cjA(a21[i*as0]) * b1t[j*bs1];
                    });
                });
            }
          }
        } 
        

        template<typename ConjType, typename MemberType>
        static 
        KOKKOS_INLINE_FUNCTION
        void trsm_left_upper(MemberType &member, const ConjType &cjA, 
                             const char diag, 
                             const int m, const int n, 
                             const T alpha, 
                             const T *__restrict__ A, const int as0, const int as1,
                             /* */ T *__restrict__ B, const int bs0, const int bs1) {          
          const T one(1.0), zero(0.0);

          // note that parallel range is different ( m*n vs m-1*n);
          if (alpha == zero)  set  (member, m, n, zero,  B, bs0, bs1);
          else {
            if (alpha != one) scale(member, m, n, alpha, B, bs0, bs1);
            if (m <= 0 || n <= 0) return;

            const bool use_unit_diag = diag == 'U'|| diag == 'u';
            T *__restrict__ B0 = B;
            for (int p=(m-1);p>=0;--p) {
              const int iend = p, jend = n;

              const T *__restrict__ a01 = A+p*as1;
              /**/  T *__restrict__ b1t = B+p*bs0;

              member.team_barrier();
              if (!use_unit_diag) {
                const T alpha11 = cjA(A[p*as0+p*as1]);
                Kokkos::parallel_for(Kokkos::TeamVectorRange(member,jend),[&](const int &j) {
                    b1t[j*bs1] /= alpha11;
                  });
                member.team_barrier();
              }

              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,jend),[&](const int &j) {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,iend),[&](const int &i) {
                      B0[i*bs0+j*bs1] -= cjA(a01[i*as0]) * b1t[j*bs1];
                    });
                });
            }
          }
        }
        
      };
      
      template<typename MemberType>
      static 
      KOKKOS_INLINE_FUNCTION
      void gemv(MemberType &member,
                const char trans, 
                const int m, const int n, 
                const T alpha, 
                const T *__restrict__ a, const int lda,
                const T *__restrict__ x, const int xs,
                const T beta,
                /* */ T *__restrict__ y, const int ys) {
        switch (trans) {
        case 'N':
        case 'n': {
          const NoConjugate cj;
          Impl::gemv(member, cj,
                     m, n, 
                     alpha, 
                     a, 1, lda, 
                     x, xs,
                     beta,
                     y, ys);
          break;
        }
        case 'T':
        case 't': {
          const NoConjugate cj;
          Impl::gemv(member, cj,
                     n, m, 
                     alpha, 
                     a, lda, 1, 
                     x, xs,
                     beta,
                     y, ys);
          break;
        }
        case 'C':
        case 'c': {
          const Conjugate cj;
          Impl::gemv(member, cj,
                     n, m, 
                     alpha, 
                     a, lda, 1, 
                     x, xs,
                     beta,
                     y, ys);
          break;
        }
        default:
          Kokkos::abort("Invalid trans character");
        }
      }


      template<typename MemberType>
      static 
      KOKKOS_INLINE_FUNCTION
      void trsv(MemberType &member,
                const char uplo, const char trans, const char diag, 
                const int m, 
                const T *__restrict__ a, const int lda,
                /* */ T *__restrict__ b, const int bs) {
        if (uplo == 'U' || uplo == 'u') {
          switch (trans) {
          case 'N':
          case 'n': {
            NoConjugate cjA;
            Impl::trsv_upper(member, cjA, diag,
                             m, 
                             a, 1, lda,
                             b, bs);
            break;
          }
          case 'T':
          case 't': {
            NoConjugate cjA;
            Impl::trsv_lower(member, cjA, diag,
                             m, 
                             a, lda, 1,
                             b, bs);        
            break;
          }
          case 'C':
          case 'c': {
            Conjugate cjA;
            Impl::trsv_lower(member, cjA, diag,
                             m, 
                             a, lda, 1,
                             b, bs);        
            break;        
          }
          default:
            Kokkos::abort("trans is not valid");
          }
        } else if (uplo == 'L' || uplo == 'l') {
          switch (trans) {
          case 'N':
          case 'n': {
            NoConjugate cjA;
            Impl::trsv_lower(member, cjA, diag,
                             m, 
                             a, 1, lda,
                             b, bs);
            break;
          }
          case 'T':
          case 't': {
            NoConjugate cjA;
            Impl::trsv_upper(member, cjA, diag,
                             m, 
                             a, lda, 1,
                             b, bs);        
            break;
          }
          case 'C':
          case 'c': {
            Conjugate cjA;
            Impl::trsv_upper(member, cjA, diag,
                             m, 
                             a, lda, 1,
                             b, bs);        
            break;        
          }
          default:
            Kokkos::abort("trans is not valid");
          }
        }
      }
      
      template<typename MemberType>
      static 
      KOKKOS_INLINE_FUNCTION
      void gemm(MemberType &member, 
                const char transa, const char transb, 
                const int m, const int n, const int k,
                const T alpha, 
                const T *__restrict__ a, int lda,
                const T *__restrict__ b, int ldb,
                const T beta,
                /* */ T *__restrict__ c, int ldc) {

        if (transa == 'N' || transa == 'n') {
          const NoConjugate cjA;
          switch (transb) {
          case 'N':
          case 'n': {
            const NoConjugate cjB;
            Impl::gemm(member, cjA, cjB,
                       m, n, k,
                       alpha,
                       a, 1, lda,
                       b, 1, ldb,
                       beta,
                       c, 1, ldc);
            break;
          }
          case 'T':
          case 't': {
            const NoConjugate cjB;
            Impl::gemm(member, cjA, cjB, 
                       m, n, k,
                       alpha,
                       a, 1, lda,
                       b, ldb, 1,
                       beta,
                       c, 1, ldc);
            break;
          }
          case 'C':
          case 'c': {
            const Conjugate cjB;
            Impl::gemm(member, cjA, cjB,
                       m, n, k,
                       alpha,
                       a, 1, lda,
                       b, ldb, 1,
                       beta,
                       c, 1, ldc);
            break;
          }
          default:
            Kokkos::abort("transa is no trans but transb is not valid");
          }
        } else if (transa == 'T' || transa == 't') {
          const NoConjugate cjA;          
          switch (transb) {
          case 'N':
          case 'n': {
            const NoConjugate cjB;          
            Impl::gemm(member, cjA, cjB,
                       m, n, k,
                       alpha,
                       a, lda, 1,
                       b, 1, ldb,
                       beta,
                       c, 1, ldc);
            break;
          }
          case 'T':
          case 't': {
            const NoConjugate cjB;          
            Impl::gemm(member, cjA, cjB,
                       m, n, k,
                       alpha,
                       a, lda, 1,
                       b, ldb, 1,
                       beta,
                       c, 1, ldc);
            break;
          }
          case 'C':
          case 'c': {
            const Conjugate cjB;          
            Impl::gemm(member, cjA, cjB,
                       m, n, k,
                       alpha,
                       a, lda, 1,
                       b, ldb, 1,
                       beta,
                       c, 1, ldc);
            break;
          }
          default:
            Kokkos::abort("transa is trans but transb is not valid");
          }
        } else if (transa == 'C' || transa == 'c') {
          const Conjugate cjA;          
          switch (transb) {
          case 'N':
          case 'n': {
            const NoConjugate cjB;          
            Impl::gemm(member, cjA, cjB,
                       m, n, k,
                       alpha,
                       a, lda, 1,
                       b, 1, ldb,
                       beta,
                       c, 1, ldc);
            break;
          }
          case 'T':
          case 't': {
            const NoConjugate cjB;          
            Impl::gemm(member, cjA, cjB, 
                       m, n, k,
                       alpha,
                       a, lda, 1,
                       b, ldb, 1,
                       beta,
                       c, 1, ldc);
            break;
          }
          case 'C':
          case 'c': {
            const Conjugate cjB;          
            Impl::gemm(member, cjA, cjB,
                       m, n, k,
                       alpha,
                       a, lda, 1,
                       b, ldb, 1,
                       beta,
                       c, 1, ldc);
            break;
          }
          default:
            Kokkos::abort("transa is conj trans but transb is not valid");
          }
        } else { 
          Kokkos::abort("transa is not valid");          
        }
      }

      
      template<typename MemberType>
      static 
      KOKKOS_INLINE_FUNCTION
      void herk(MemberType &member, 
                const char uplo, const char trans, 
                const int n, const int k,
                const T alpha, 
                const T *__restrict__ a, const int lda,
                const T beta,
                /* */ T *__restrict__ c, const int ldc) {
        if (uplo  == 'U' || uplo  == 'u')
          switch (trans) {
          case 'N':
          case 'n': {
            const NoConjugate cjA;
            const Conjugate cjB;
            Impl::herk_upper(member, cjA, cjB, 
                             n, k, 
                             alpha, 
                             a, 1, lda,
                             beta, 
                             c, 1, ldc);
            break;
          }
          case 'C':
          case 'c': {
            const NoConjugate cjA;
            const Conjugate cjB;
            Impl::herk_upper(member, cjA, cjB,
                             n, k,
                             alpha, 
                             a, lda, 1,
                             beta, 
                             c, 1, ldc);
            break;
          }
          default:
            Kokkos::abort("trans is not valid");                    
          }
        else if (uplo  == 'L' || uplo  == 'l')
          switch (trans) {
          case 'N':
          case 'n': {
            const NoConjugate cjA;
            const Conjugate cjB;
            Impl::herk_lower(member, cjA, cjB,
                             n, k, 
                             alpha, 
                             a, 1, lda,
                             beta, 
                             c, 1, ldc);
            break;
          }
          case 'C':
          case 'c': {
            const NoConjugate cjA;
            const Conjugate cjB;
            Impl::herk_lower(member, cjA, cjB, 
                             n, k, 
                             alpha, 
                             a, lda, 1,
                             beta, 
                             c, 1, ldc);
            break;
          }
          default:
            Kokkos::abort("trans is not valid");                    
          }
        else
          Kokkos::abort("uplo is not valid");                              
      }

      
      template<typename MemberType>
      static 
      KOKKOS_INLINE_FUNCTION
      void trsm(MemberType &member, 
                const char side, const char uplo, const char trans, const char diag,
                const int m, const int n, 
                const T alpha, 
                const T *__restrict__ a, const int lda,
                /* */ T *__restrict__ b, const int ldb) {
        ///
        /// side left 
        ///
        if (side == 'L' || side == 'l') {
          if (uplo == 'U' || uplo == 'u') {
            switch (trans) {
            case 'N':
            case 'n': {
              NoConjugate cjA;
              Impl::trsm_left_upper(member, cjA, 
                                    diag, 
                                    m, n, 
                                    alpha, 
                                    a, 1, lda,
                                    b, 1, ldb);
              break;
            }
            case 'T':
            case 't': {
              NoConjugate cjA;
              Impl::trsm_left_lower(member, cjA, 
                                    diag, 
                                    m, n, 
                                    alpha, 
                                    a, lda, 1,
                                    b, 1, ldb);
              break;
            }
            case 'C':
            case 'c': {
              Conjugate cjA;
              Impl::trsm_left_lower(member, cjA, 
                                    diag, 
                                    m, n, 
                                    alpha, 
                                    a, lda, 1,
                                    b, 1, ldb);
              break;        
            }
            default:
              Kokkos::abort("trans is not valid");
            }
          } else if (uplo == 'L' || uplo == 'l') {
            switch (trans) {
            case 'N':
            case 'n': {
              NoConjugate cjA;
              Impl::trsm_left_lower(member, cjA, 
                                    diag, 
                                    m, n, 
                                    alpha, 
                                    a, 1, lda,
                                    b, 1, ldb);
              break;
            }
            case 'T':
            case 't': {
              NoConjugate cjA;
              Impl::trsm_left_upper(member, cjA, 
                                    diag, 
                                    m, n, 
                                    alpha, 
                                    a, lda, 1,
                                    b, 1, ldb);
              break;
            }
            case 'C':
            case 'c': {
              Conjugate cjA;
              Impl::trsm_left_upper(member, cjA, 
                                    diag, 
                                    m, n, 
                                    alpha, 
                                    a, lda, 1,
                                    b, 1, ldb);
              break;        
            }
            default:
              Kokkos::abort("trans is not valid");
            }
          }
        } 
        ///
        /// side right 
        ///
        else if (side == 'R' || side == 'r') {
          Kokkos::abort("right side is not implemented");
        }
      }
      
    };
    
}

#endif
