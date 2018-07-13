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

      integer lda, lea, lia, lno, lnv
!     Fortran-77 restriction: integer*4 for array index.
!     integer*4: 2^(4*8-1) - 1 = 2147483647.
      integer lst

!     Max lst for Fortran array integer*4 index.
!      parameter (lda=100000, lst=1700000000, lea=500000, lia=80000,      &
!     &           lno=11, lnv=7)

!     For Laptop, NO = 7.
!      parameter (lda=100000, lst=400000000, lea=500000, lia=80000,         &
!     &           lno=10, lnv=7)
!     For Laptop, NO = 9.
!      parameter (lda=100000, lst=700000000, lea=500000, lia=80000,         &
!     &           lno=10, lnv=7)
!     For Workstation.
!      parameter (lda=100000, lst=900000000, lea=500000, lia=80000,         &
!     &           lno=10, lnv=7)
!     For Cluster, NO = 11. Increase virtual memory for SGE by:
!       qsub -l mem_free=50G,h_vmem=50G -q...
      parameter (lda=100000, lst=2000000000, lea=500000, lia=80000,        &
     &           lno=11, lnv=7)


      integer         nda, ndamaxi
      common /fordes/ nda, ndamaxi

      double precision cc, eps, epsmac
      common /da/      cc(lst), eps, epsmac

      integer      i1, i2, ie1, ie2, ieo
      integer      ia1, ia2, ifi, idano
      integer      idanv, idalm, idall
      integer      idapo
      integer      nst, nomax, nvmax, nmmax, nocut, lfi
      common /dai/ i1(lst), i2(lst), ie1(lea), ie2(lea), ieo(lea),         &
     &             ia1(0:lia), ia2(0:lia), ifi(lea), idano(lda),           &
     &             idanv(lda), idapo(lda), idalm(lda), idall(lda),         &
     &             nst, nomax, nvmax, nmmax, nocut, lfi

      double precision facint
      common /factor/  facint(0:lno)
!-----------------------------------------------------------------------------9
