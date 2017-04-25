!     PARAMETERS:
!
!     LDA: MAXIMUM NUMBER OF DA-VECTORS;    CAN BE CHANGED QUITE ARBITRARILY
!     LST: LENGTH OF MAIN STORAGE STACK;    CAN BE CHANGED QUITE ARBITRARILY
!     LEA: MAXIMUM NUMBER OF MONOMIALS;     CAN BE INCREASED FOR LARGE NO, NV ( (no+nv) over no )
!     LIA: DIMENSION OF IA1, IA2;           CAN BE INCREASED FOR LARGE NO, NV ( (no+1)^nv/2 )
!     LNO: MAXIMUM ORDER;                   CAN BE INCREASED TO ABOUT 1000
!     LNV: MAXIMUM NUMBER OF VARIABLES;     CAN BE INCREASED TO ABOUT 1000
!
!-----------------------------------------------------------------------------1

      integer       lda, lea, lia, lno, lnv
      integer       lst

!      parameter (lda=70000, lst=300000000, lea=171000, lia=66000,       &
!      parameter (lda=70000, lst=900000000, lea=171000, lia=66000,       &
!     &           lno=15, lnv=7)
!     Laptop
!      parameter (lda=30000, lst=100000000, lea=3005, lia=8000,          &
!     &           lno=5, lnv=10)
!      parameter (lda=100000, lst=700000000, lea=11440, lia=11000,        &
!     &           lno=10, lnv=7)
!      parameter (lda=30000, lst=900000000, lea=81000, lia=370000000,     &
!     &           lno=5, lnv=22)
!      parameter (lda=100000, lst=2147000000, lea=500000, lia=80000,         &
!     &           lno=16, lnv=7)
!      parameter (lda=100000, lst=1000000000, lea=500000, lia=80000,         &
!     &           lno=16, lnv=7)
!      parameter (lda=200000, lst=2000000000, lea=500000, lia=80000,         &
!      parameter (lda=30000, lst=200000000, lea=500000, lia=80000,         &
      parameter (lda=80000, lst=700000000, lea=500000, lia=80000,         &
     &           lno=16, lnv=7)


      integer       nda, ndamaxi
      common /fordes/   nda, ndamaxi

      double precision  cc, eps, epsmac
      common /da/       cc(lst), eps, epsmac

      integer       i1, i2, ie1, ie2, ieo
      integer       ia1, ia2, ifi, idano
      integer       idanv, idalm, idall
      integer       idapo
      integer       nst, nomax, nvmax, nmmax, nocut, lfi
      common /dai/      i1(lst), i2(lst), ie1(lea), ie2(lea), ieo(lea), &
     &                  ia1(0:lia), ia2(0:lia), ifi(lea), idano(lda),   &
     &                  idanv(lda), idapo(lda), idalm(lda), idall(lda), &
     &                  nst, nomax, nvmax, nmmax, nocut, lfi

      double precision  facint
      common /factor/   facint(0:lno)
!-----------------------------------------------------------------------------9
