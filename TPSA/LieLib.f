      module lielib
      contains

      subroutine lieinit(no1,nv1,nd1,ndpt1,iref1,nis)                   &
     &           bind(C, name="lieinit_")
      use iso_c_binding
      implicit none
      integer(C_LONG) no1,nv1,nd1,ndpt1,iref1,nis

      integer i,ndc1,ndim
      double precision ang,ra,st
!! Lieinit initializes AD Package and Lielib
      parameter (ndim=3)
      dimension st(ndim),ang(ndim),ra(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer iref
      common /resfile/iref
      integer itu
      common /tunedef/itu
      integer lienot
      common /dano/lienot
      integer idpr
      common /printing/ idpr
      double precision xintex
      common /integratedex/ xintex(0:20)
      integer idao,is,iscrri
      double precision rs
      common/dascr/is(100),rs(100),iscrri(100),idao
      integer nplane
      double precision epsplane,xplane
      common /choice/ xplane(ndim),epsplane,nplane(ndim)
!+CA DASCR
      call daexter
      do i=1,ndim
        nplane(i)=2*i-1
        ang(i)=0.d0
        ra(i)=0.d0
        st(i)=1.d0
      enddo
      no=no1
      nv=nv1
      nd=nd1
      nd2=2*nd1
      do i=1,100
        is(i)=0
      enddo
      idpr=1
      write(*, 100) 'lieinit: no = ', no, ', nv = ', nv, ', nd = ',     &
     &      nd, ', nd2 = ', nd2, ', ndpt = ', ndpt, ', idpr = ', idpr
 100  format(6(a, i0))
      call daini(no,nv,0)
      if(nis.gt.0)call etallnom(is,nis,'$$IS      ')
      if(ndpt1.eq.0) then
        ndpt=0
        ndt=0
        ndc1=0
      else
        ndpt=ndpt1
        ndc1=1
        if(ndpt.eq.nd2) then
          ndt=nd2-1
        else
          ndt=nd2
          if(ndpt.ne.nd2-1) then
            write(6,*) ' LETHAL ERROR IN LIEINIT', ndpt, nd2-1
            stop
          endif
        endif
      endif
      ndc=ndc1
      ndc2=2*ndc1
      iref=0
      call initpert(st,ang,ra)
      iref=iref1
      if(iref1.eq.0) then
        itu=0
      else
        itu=1
      endif
      if(iref1.eq.0) iref=-1

      if(idpr.eq.1) write(6, 200) ' NO = ', no,' IN DA-CALCULATIONS '
 200  format(a, i0, a)

      do i=0,20
        xintex(i)=0.d0
      enddo
      xintex(          0)=       1.000000000000000
      xintex(          1)=  5.000000000000000e-001
      xintex(          2)=  8.333333333333334e-002
      xintex(          3)=  0.000000000000000e+000
      xintex(          4)= -1.388888888888898e-003
      xintex(          5)=  0.000000000000000e+000
      xintex(          6)=  3.306878306878064e-005
      xintex(          7)= 0.d0
      xintex(          8)= -8.267195767165669e-007
      xintex(          9)=  0.d0
      xintex(         10)=  4.592886537931051e-008

      return
      end subroutine

      subroutine flowpara(ifl,jtu)
      implicit none
      integer iflow,jtune
      common /vecflow/ iflow,jtune
      integer ifl,jtu
      iflow=ifl
      jtune=jtu
      return
      end subroutine

      subroutine pertpeek(st,ang,ra)
      implicit none
      integer i,ndim,ndim2,nreso,ntt
      double precision ang,ra,st
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
      dimension st(ndim),ang(ndim),ra(ndim)
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer mx,nres
      common /reson/mx(ndim,nreso),nres
      integer iref
      common /resfile/iref
      do i=1,nd
        st(i)=sta(i)
        ang(i)=angle(i)
        ra(i)=rad(i)
      enddo
      return
      end subroutine

      subroutine inputres(mx1,nres1)
      implicit none
      integer i,j,ndim,ndim2,nreso,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
      integer mx1(ndim,nreso),nres1
      integer mx,nres
      common /reson/mx(ndim,nreso),nres

      nres=nres1
      do i=1,nreso
        do j=1,ndim
          mx(j,i)=0
        enddo
      enddo

      do i=1,nres
        do j=1,ndim
          mx(j,i)=mx1(j,i)
        enddo
      enddo
      return
      end subroutine

      subroutine respoke(mres,nre,ire)
      implicit none
      integer i,ire,j,ndim,ndim2,nre,nreso,ntt
      double precision ang,ra,st
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
      integer mres(ndim,nreso)
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer mx,nres
      common /reson/mx(ndim,nreso),nres
      integer iref
      common /resfile/iref
      dimension ang(ndim),ra(ndim),st(ndim)
      iref=ire
      nres=nre
      do j=1,nreso
        do i=1,nd
          mx(i,j)=mres(i,j)
        enddo
      enddo
      call initpert(st,ang,ra)
      return
      end subroutine

      subroutine liepeek(iia,icoast)
      implicit none
      integer ndim,ndim2,nreso,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer iia(*),icoast(*)

      iia(1)=no
      iia(2)=nv
      iia(3)=nd
      iia(4)=nd2
      icoast(1)=ndc
      icoast(2)=ndc2
      icoast(3)=ndt
      icoast(4)=ndpt

      return
      end subroutine

      subroutine lienot(not)
      implicit none
      integer no,not

      call danot(not)
      no=not

      return
      end subroutine

      subroutine etallnom1(x,nom)
      implicit none
      integer n,nd2
! CREATES A AD-VARIABLE WHICH CAN BE DESTROYED BY DADAL
! allocates vector of n polynomials and give it the name NOM=A10
      integer x,i1(4),i2(4)
      character*10 nom
      x=0
      call daallno1(x,nom)
      return
      end subroutine

      subroutine etallnom(x,n,nom)
      implicit none
      integer i,n,nd2
! CREATES A AD-VARIABLE WHICH CAN BE DESTROYED BY DADAL
! allocates vector of n polynomials and give it the name NOM=A10
      integer x(*),i1(4),i2(4)
      character*10 nom
      do i=1,iabs(n)
        x(i)=0
      enddo
      call daallno(x,iabs(n),nom)
      if(n.lt.0) then
        call liepeek(i1,i2)
        nd2=i1(4)
        do i=nd2+1,-n
          call davar(x(i),0.d0,i)
        enddo
      endif
      return
      end subroutine

      subroutine etall(x,n)
      implicit none
      integer i,n,nd2
! allocates vector of n polynomials
      integer x(*),i1(4),i2(4)
      do i=1,iabs(n)
        x(i)=0
      enddo
      call daallno(x,iabs(n),'ETALL     ')
      if(n.lt.0) then
        call liepeek(i1,i2)
        nd2=i1(4)
        do i=nd2+1,-n
          call davar(x(i),0.d0,i)
        enddo
      endif
      return
      end subroutine

      subroutine etall1(x)
      implicit none
      integer x
      call daallno1(x,'ETALL     ')
      return
      end subroutine

      subroutine etppulnv(x,xi,xff)
      implicit none
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer  x(*)
      double precision xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nv
        xii(i)=xi(i)
      enddo
      do i=nv+1,ntt
        xii(i)=0.d0
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nv
        xff(i)=xf(i)
      enddo

      return
      end subroutine

      subroutine etmtree(y,x) bind(C, name="etmtree_")
      use iso_c_binding
      implicit none
      integer(C_LONG) y(*), x(*)

      integer i,ie,iv,ndim,ndim2,nt,ntt
! ROUTINES USING THE MAP IN AD-FORM
! Note, no must be set to nomax
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie(ntt),iv(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      nt=nv-nd2
      if(nt.gt.0) then
        call etallnom(ie,nt,'IE        ')
        do i=nd2+1,nv
          call davar(ie(i-nd2),0.d0,i)
        enddo
        do i=nd2+1,nv
          iv(i)=ie(i-nd2)
        enddo
      endif
      do i=1,nd2
        iv(i)=y(i)
      enddo
      call mtree(iv,nv,x,nv)
      if(nt.gt.0) then
        call dadal(ie,nt)
      endif
      return
      end subroutine

      subroutine etppush(x,xi) bind(C, name="etppush_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*)
      real(C_DOUBLE) xi(*)

      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      double precision xf(ntt),xii(ntt)

      do i=1,nd2
        xii(i)=xi(i)
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nd2
        xi(i)=xf(i)
      enddo

      return
      end subroutine

      subroutine etppush2(x,xi,xff) bind(C, name="etppush2_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*)
      real(C_DOUBLE) xi(*), xff(*)

      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      double precision xf(ntt),xii(ntt)

      do i=1,nd2
        xii(i)=xi(i)
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nd2
        xff(i)=xf(i)
      enddo

      return
      end subroutine

      subroutine ppushlnv(x,xi,xff,nd1)
      implicit none
      integer i,nd1,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer  x(*)
      double precision xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nd1
        xii(i)=xi(i)
      enddo
      do i=nd1+1,ntt
        xii(i)=0.d0
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nd1
        xff(i)=xf(i)
      enddo

      return
      end subroutine

      subroutine etcct(x,y,z) bind(C, name="etcct_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*),y(*),z(*)

      integer i,ie,iv,ndim,ndim2,nt,ntt
!  Z=XoY
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie(ntt),iv(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      nt=nv-nd2
      if(nt.gt.0) then
        call etallnom(ie,nt,'IE        ')
        do i=nd2+1,nv
          call davar(ie(i-nd2),0.d0,i)
        enddo
        do i=nd2+1,nv
          iv(i)=ie(i-nd2)
        enddo
      endif
      do i=1,nd2
        iv(i)=y(i)
      enddo
      call dacct(x,nd2,iv,nv,z,nd2)
      if(nt.gt.0) then
        call dadal(ie,nt)
      endif
      return
      end subroutine

      subroutine trx(h,rh,y) bind(C, name="trx_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h,rh
      integer(C_LONG) y(*)
      integer i,ie,iv,ndim,ndim2,nt,ntt
!  :RH: = Y :H: Y^-1 =  :HoY:
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie(ntt),iv(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
!
      integer h_buf(1), rh_buf(1)

      nt=nv-nd2
      if(nt.gt.0) then
        call etallnom(ie,nt,'IE        ')
        do i=nd2+1,nv
          call davar(ie(i-nd2),0.d0,i)
        enddo
        do i=nd2+1,nv
          iv(i)=ie(i-nd2)
        enddo
      endif
      do i=1,nd2
        iv(i)=y(i)
      enddo
      h_buf(1) = h
      rh_buf(1) = rh
      call dacct(h_buf(1),1,iv,nv,rh_buf(1),1)
      rh = rh_buf(1)
      if(nt.gt.0) then
        call dadal(ie,nt)
      endif
      return
      end subroutine

      subroutine trxflo(h,rh,y)
      implicit none
      integer j,k,ndim,ndim2,ntt
!  *RH* = Y *H* Y^-1  CHANGE OF A VECTOR FLOW OPERATOR
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer h(*),rh(*),y(*)
      integer yi(ndim2),ht(ndim2),b1,b2
!
!
      call etallnom(yi,nd2  ,'YI        ')
      call etallnom(ht,nd2  ,'HT        ')
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')

      call etinv(y,yi)
!----- HT= H o Y
      call etcct(h,y,ht)
!----
      call daclrd(rh)
      do j=1,nd2
        do k=1,nd2
          call dader(k,yi(j),b1)
          call trx(b1,b2,y)
          call damul(b2,ht(k),b1)
          call daadd(b1,rh(j),b2)
          call dacop(b2,rh(j))
        enddo
      enddo

      call dadal1(b2)
      call dadal1(b1)
      call dadal(ht,nd2)
      call dadal(yi,nd2)
      return
      end subroutine

      subroutine simil(a,x,ai,y)
      implicit none
      integer ndim,ndim2,ntt
!  Y= AoXoAI
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer  x(*),y(*),a(*),ai(*)

      integer w(ndim2),v(ndim2)
!
      call etallnom(w,nd2  ,'W         ')
      call etallnom(v,nd2  ,'V         ')

      call etcct(a,x,w)
      call etcct(w,ai,v)

      call dacopd(v,y)

      call dadal(v,nd2)
      call dadal(w,nd2)
      return
      end subroutine

      subroutine etini(x)
      implicit none
      integer i,ndim,ndim2,ntt
!  X=IDENTITY
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer x(*)
!*DAEXT(NO,NV) X(NDIM2)
      do i=1,nd2
        call davar(x(i),0.d0,i)
      enddo
      return
      end subroutine

      subroutine etinv(x,y) bind(C, name="etinv_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*), y(*)

      integer i,ie1,ie2,iv1,iv2,ndim,ndim2,nt,ntt
! Y=X^-1
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie1(ntt),ie2(ntt),iv1(ntt),iv2(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      nt=nv-nd2
      if(nt.gt.0) then
        do i=1,nt
          ie1(i)=0
          ie2(i)=0
        enddo
        call etallnom(ie1,nt,'IE1       ')
        call etallnom(ie2,nt,'IE2       ')
        do i=nd2+1,nv
          call davar(ie1(i-nd2),0.d0,i)
        enddo
        do i=nd2+1,nv
          iv1(i)=ie1(i-nd2)
          iv2(i)=ie2(i-nd2)
        enddo
      endif
      do i=1,nd2
        iv1(i)=x(i)
        iv2(i)=y(i)
      enddo

      call dainv(iv1,nv,iv2,nv)
      if(nt.gt.0) then
        call dadal(ie2,nt)
        call dadal(ie1,nt)
      endif
      return
      end subroutine

      subroutine etpin(x,y,jj) bind(C, name="etpin_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*), y(*), jj(*)

      integer i,ie1,ie2,iv1,iv2,ndim,ndim2,nt,ntt
!  Y=PARTIAL INVERSION OF X SEE BERZ'S PACKAGE
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie1(ntt),ie2(ntt),iv1(ntt),iv2(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      nt=nv-nd2
      if(nt.gt.0) then
        do i=1,nt
          ie1(i)=0
          ie2(i)=0
        enddo
        call etallnom(ie1,nt,'IE1       ')
        call etallnom(ie2,nt,'IE2       ')
        do i=nd2+1,nv
          call davar(ie1(i-nd2),0.d0,i)
        enddo
        do i=nd2+1,nv
          iv1(i)=ie1(i-nd2)
          iv2(i)=ie2(i-nd2)
        enddo
      endif
      do i=1,nd2
        iv1(i)=x(i)
        iv2(i)=y(i)
      enddo

      call dapin(iv1,nv,iv2,nv,jj)
      if(nt.gt.0) then
        call dadal(ie2,nt)
        call dadal(ie1,nt)
      endif
      return
      end subroutine

      subroutine dapek0(v,x,jj)
      implicit none
      integer i,jj,ndim2,ntt
      double precision x
!- MORE EXTENSIONS OF BASIC BERZ'S PACKAGE
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*),jd(ntt)
      dimension x(*)
      do i=1,ntt
        jd(i)=0
      enddo
      do i=1,jj
        call dapek(v(i),jd,x(i))
      enddo
      return
      end subroutine

      subroutine dapok0(v,x,jj)
      implicit none
      integer i,jj,ndim2,ntt
      double precision x
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*),jd(ntt)
      dimension x(*)
      do i=1,ntt
        jd(i)=0
      enddo
      do i=1,jj
        call dapok(v(i),jd,x(i))
      enddo
      return
      end subroutine

      subroutine dapokzer(v,jj)
      implicit none
      integer i,jj,ndim2,ntt
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*),jd(ntt)
      do i=1,ntt
        jd(i)=0
      enddo
      do i=1,jj
        call dapok(v(i),jd,0.d0)
      enddo
      return
      end subroutine

      subroutine davar0(v,x,jj)
      implicit none
      integer i,jj,ndim2,ntt
      double precision x
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*)
      dimension x(*)
      do i=1,jj
        call davar(v(i),x(i),i)
      enddo
      return
      end subroutine

      subroutine comcfu(b,f1,f2,c) bind(C, name="comcfu_")
      use iso_c_binding
      implicit none
      integer(C_LONG) b(*), c(*)

      abstract interface
        function f1(j) bind(C)
          import :: c_long, c_double
          real(c_long), intent(in) :: j(*)
          real(c_double) :: f1
        end function
      end interface

      abstract interface
        function f2(j) bind(C)
          import :: c_long, c_double
          real(c_long), intent(in) :: j(*)
          real(c_double) :: f2
        end function
      end interface

! Complex dacfu
      integer t(4)
      call etall(t,4)

      call dacfu(b(1),f1,t(1))
      call dacfu(b(1),f2,t(2))
      call dacfu(b(2),f1,t(3))
      call dacfu(b(2),f2,t(4))

      call dasub(t(1),t(4),c(1))
      call daadd(t(2),t(3),c(2))
      call dadal(t,4)
      return
      end subroutine

      subroutine take(h,m,ht) bind(C, name="take_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h, ht, m

      integer i,ndim,ntt
      double precision r
!  HT= H_M  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
      parameter (ndim=3)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer j(ntt)

      integer b1,b2,b3
!
!
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')
      call etallnom1(b3,'B3        ')

      if(no.ge.2) then
        if(m.eq.0) then
          do i=1,ntt
            j(i)=0
          enddo
          call dapek(h,j,r)
          call dacon(ht,r)
        else
          call danot(m)
          call dacop(h,b1)
          call danot(m-1)
          call dacop(b1,b2)
          call danot(no)
          call dasub(b1,b2,b3)
          call dacop(b3,ht)
        endif
      else
        do i=1,ntt
          j(i)=0
        enddo
        if(m.eq.0) then
          call dapek(h,j,r)
          call dacon(ht,r)
        elseif(m.eq.1)  then
          do i=1,nv
            j(i)=1
            call dapek(h,j,r)
            call dapok(b3,j,r)
            j(i)=0
          enddo
          call dacop(b3,ht)
        else
          call daclr(ht)
        endif
      endif

      call dadal1(b3)
      call dadal1(b2)
      call dadal1(b1)
      return
      end subroutine

      subroutine taked(h,m,ht) bind(C, name="taked_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h(*), m, ht(*)

      integer i,ndim2,ntt
!  \VEC{HT}= \VEC{H_M}  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer j(ntt)

      integer b1,b2,x(ndim2)
!
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')
      call etallnom(x,nd2  ,'X         ')


      do i=1,ntt
        j(i)=0
      enddo

      do   i=1,nd2
        call take(h(i),m,ht(i))
      enddo
      call dadal(x,nd2)
      call dadal1(b2)
      call dadal1(b1)
      return
      end subroutine

      subroutine daclrd(h)
      implicit none
      integer i,ndim2,ntt
! clear a map : a vector of nd2 polynomials
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer h(*)
      do i=1,nd2
        call daclr(h(i))
      enddo
      return
      end subroutine

      subroutine dacopd(h,ht)
      implicit none
      integer i,ndim2,ntt
!    H goes into HT  (nd2 array)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer h(*),ht(*)
      do i=1,nd2
        call dacop(h(i),ht(i))
      enddo
      return
      end subroutine

      subroutine dacmud(h,sca,ht)
      implicit none
      integer i,ndim2,ntt
      double precision sca
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer h(*),ht(*)
      do i=1,nd2
        call dacmu(h(i),sca,ht(i))
      enddo
      return
      end subroutine

      subroutine dalind(h,rh,ht,rt,hr)
      implicit none
      integer i,ndim2
      double precision rh,rt
      parameter (ndim2=6)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer h(*),ht(*),hr(*)

      integer b(ndim2)
!
      call etallnom(b,nd2  ,'B         ')

      do i=1,nd2
        call dalin(h(i),rh,ht(i),rt,b(i))
      enddo
      call dacopd(b,hr)
      call dadal(b,nd2)
      return
      end subroutine

      subroutine prresflo(h,eps,mfile) bind(C, name="prresflo_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h(*), mfile
      real(C_DOUBLE) eps

      integer i,ndim2,ntt
      double precision deps,filtres
!  print a map   in resonance basis for human consumption (useless)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer b(ndim2),c(ndim2)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ifilt
      common /filtr/ ifilt
      external filtres
      call etall(b,nd2)
      call etall(c,nd2)
      call dacopd(h,c)
      do i=1,nd2
        ifilt=(-1)**i
        call  dacfu(c(i),filtres,h(i))
      enddo

      deps=-1.d0
      call daeps(deps)
      call daeps(eps)

      call dacopd(c,h)
      call daeps(deps)
      call  dadal(c,nd2)
      call  dadal(b,nd2)
      return
      end subroutine

      real(C_DOUBLE) function filtres(j) bind(C, name="filtres_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ic,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      integer J(NTT)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer ifilt
      common /filtr/ ifilt
      filtres=1.d0
      ic=0
      do i=1,(nd2-ndc2)
        ic=ic+j(i)*(-1)**(i+1)
      enddo
      ic=ic+ifilt
      if(ic.lt.0) filtres=0.d0
      if(ic.eq.0.and.ifilt.eq.1) then
        filtres=0.0d0
      endif
      return
      end function

      subroutine daflo(h,x,y) bind(C, name="daflo_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h(*), x, y

      integer i,ndim,ndim2,ntt
! LIE EXPONENT ROUTINES WITH FLOW OPERATORS

!     \VEC{H}.GRAD X =Y
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer b1,b2,b3
!
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')
      call etallnom1(b3,'B3        ')

      call daclr(b1)
      call daclr(b2)
      do i=1,nd2
        call dader(i,x,b2)
        call damul(b2,h(i),b3)
        call daadd(b3,b1,b2)
        call dacop(b2,b1)
      enddo
      call dacop(b1,y)
      call dadal1(b3)
      call dadal1(b2)
      call dadal1(b1)
      return
      end subroutine

      subroutine daflod(h,x,y) bind(C, name="daflod_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h(*), x(*), y(*)

      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer b1(ndim2),b2(ndim2)
!
      call etall(b1,nd2)
      call etall(b2,nd2)

      call dacopd(h,b1)
      call dacopd(x,b2)

      do i=1,nd2
        call daflo(b1,b2(i),y(i))
      enddo

      call dadal(b1,nd2)
      call dadal(b2,nd2)
      return
      end subroutine

      subroutine intd(v,h,sca) bind(C, name="intd_")
      use iso_c_binding
      implicit none
      integer(C_LONG) v(*), h
      real(C_DOUBLE) sca

      integer i,ndim,ndim2,ntt
      double precision dlie
! IF SCA=-1.D0
!     \VEC{V}.GRAD   = J GRAD H . GRAD = :H:

! IF SCA=1.D0
!     \VEC{V}.GRAD  = GRAD H . GRAD
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      external dlie

      integer b1,b2,b3,b4,x(ndim2)
!
!
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')
      call etallnom1(b3,'B3        ')
      call etallnom1(b4,'B4        ')
      call etallnom(x,nd2  ,'X         ')

      call daclr(b4)
      call daclr(h)
      call etini(x)
      do i=1,nd
        call dacfu(v(2*i-1),dlie,b3)
        call dacfu(v(2*i),dlie,b1)
        call damul(b1,x(2*i-1),b2)
        call damul(b3,x(2*i),b1)
        call dalin(b2,1.d0,b1,sca,b3)
        call daadd(b3,b4,b2)
        call dacop(b2,b4)
      enddo
      call dacop(b4,h)
      call dadal(x,nd2)
      call dadal1(b4)
      call dadal1(b3)
      call dadal1(b2)
      call dadal1(b1)
      return
      end subroutine

      subroutine difd(h1,v,sca) bind(C, name="difd_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h1, v(*)
      real(C_DOUBLE)  sca

      integer i,ndim,ndim2,ntt
! INVERSE OF INTD ROUTINE
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer b1,h
!     J.B. 01/02/23: Changed etall1 to etallnom1.
      call etallnom1(b1, 'B1        ')
      call etallnom1(h, 'H         ')
      call dacop(h1,h)
      do i=1,nd
        call dader(2*i-1,h,v(2*i))
        call dader(2*i,h,b1)
        call dacmu(b1,sca,v(2*i-1))
      enddo
      call dadal1(h)
      call dadal1(b1)
      return
      end subroutine

      subroutine expflo(h,x,y,eps,nrmax)
      implicit none
      integer i,ndim,ndim2,nrmax,ntt
      double precision coe,eps,r,rbefore
! DOES EXP( \VEC{H} ) X = Y
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer idpr
      common /printing/ idpr
      integer h(*),x,y
      integer b1,b2,b3,b4
      logical more
!
!
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')
      call etallnom1(b3,'B3        ')
      call etallnom1(b4,'B4        ')

      call dacop(x,b4)
      call dacop(x,b1)
      more=.true.
      rbefore=1.d30
      do i=1,nrmax
        coe=1.d0/dble(i)
        call dacmu(b1,coe,b2)
        call daflo(h,b2,b1)
        call daadd(b4,b1,b3)
        call daabs(b1,r)
        if(more) then
          if(r.gt.eps) then
            rbefore=r
            goto 100
          else
            rbefore=r
            more=.false.
          endif
        else
          if(r.ge.rbefore) then
            call dacop(b3,y)
            call dadal1(b4)
            call dadal1(b3)
            call dadal1(b2)
            call dadal1(b1)
            return
          endif
          rbefore=r
        endif
100     continue
        call dacop(b3,b4)
      enddo
      if(idpr.ge.0) then
        write(6,*) ' NORM ',eps,' NEVER REACHED IN EXPFLO '
        write(6,*) 'NEW IDPR '
        read(5,*) idpr
      endif
      call dacop(b3,y)
      call dadal1(b4)
      call dadal1(b3)
      call dadal1(b2)
      call dadal1(b1)
      return
      end subroutine

      subroutine expflod(h,x,w,eps,nrmax)
      implicit none
      integer j,ndim,ndim2,nrmax,ntt
      double precision eps
! DOES EXP( \VEC{H} ) \VEC{X} = \VEC{Y}
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer x(*),w(*),h(*)
      integer b0,v(ndim2)
!
!
      call etallnom1(b0,'B0        ')
      call etallnom(v,nd2  ,'V         ')

      call dacopd(x,v)
      do j=1,nd2
        call expflo(h,v(j),b0,eps,nrmax)
        call dacop(b0,v(j))
      enddo
      call dacopd(v,w)
      call dadal(v,nd2)
      call dadal1(b0)
      return
      end subroutine

      subroutine facflo(h,x,w,nrmin,nrmax,sca,ifac)
      implicit none
      integer i,ifac,ndim,ndim2,nmax,nrmax,nrmin,ntt
      double precision eps,sca
! IFAC=1
! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX ) X= Y
! IFAC=-1
! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) X= Y
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer x,w,h(*)
      integer bm(ndim2),b0(ndim2),v
!
      call etallnom(bm,nd2  ,'BM        ')
      call etallnom(b0,nd2  ,'B0        ')
      call etallnom1(v,'V         ')

      call dacop(x,v)

      eps=-1.d0
      call daeps(eps)
      nmax=100
!
! IFAC =1 ---> V = EXP(:SCA*H(NRMAX):)...EXP(:SCA*H(NRMIN):)X
      if(ifac.eq.1) then
        do i=nrmax,nrmin,-1
          call taked(h,i,b0)
          call dacmud(b0,sca,bm)

          call expflo(bm,v,b0(1),eps,nmax)
          call dacop(b0(1),v)
        enddo
      else
! IFAC =-1 ---> V = EXP(:SCA*H(NRMIN):)...EXP(:SCA*H(NRMAX):)X
        do i=nrmin,nrmax
          call taked(h,i,b0)
          call dacmud(b0,sca,bm)

          call expflo(bm,v,b0(1),eps,nmax)
          call dacop(b0(1),v)
        enddo
      endif
      call dacop(v,w)
      call dadal1(v)
      call dadal(b0,nd2)
      call dadal(bm,nd2)
      return
      end subroutine

      subroutine facflod(h,x,w,nrmin,nrmax,sca,ifac)
      implicit none
      integer i,ifac,ndim,ndim2,nrmax,nrmin,ntt
      double precision sca
! IFAC=1
! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX )  \VEC{X}= \VEC{Y}
! IFAC=-1
! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) \VEC{X}= \VEC{Y}
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer x(*),w(*),h(*)

      do i=1,nd2
        call facflo(h,x(i),w(i),nrmin,nrmax,sca,ifac)
      enddo

      return
      end subroutine

      subroutine fexpo(h,x,w,nrmin,nrmax,sca,ifac)                      &
     &           bind(C, name="fexpo_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h, x(*), w(*), nrmin, nrmax, ifac
      real(C_DOUBLE) sca

      integer ndim,ndim2,nrma,nrmi,ntt
!   WRAPPED ROUTINES FOR THE OPERATOR  \VEC{H}=:H:
! WRAPPING FACFLOD
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      integer v(ndim2)

      nrmi=nrmin-1
      nrma=nrmax-1
      call etall(v,nd2)
      call difd(h,v,-1.d0)
      call facflod(v,x,w,nrmi,nrma,sca,ifac)

      call dadal(v,nd2)

      return
      end subroutine

      subroutine etcom(x,y,h)
      implicit none
      integer i,j,ndim,ndim2,ntt
! ETCOM TAKES THE BRACKET OF TWO VECTOR FIELDS.
      parameter (ndim2=6)
      parameter (ndim=3)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer h(*),x(*),y(*),t1,t2,t3(ndim2)

      call etall1(t1)
      call etall1(t2)
      call etall(t3,nd2)

      do j=1,nd2
        do i=1,nd2

          call dader(i,x(j),t1)
          call dader(i,y(j),t2)
          call damul(x(i),t2,t2)
          call damul(y(i),t1,t1)
          call dalin(t2,1.d0,t1,-1.d0,t1)
          call daadd(t1,t3(j),t3(j))

        enddo
      enddo

      call dacopd(t3,h)

      call dadal1(t1)
      call dadal1(t2)
      call dadal(t3,nd2)
      return
      end subroutine

      subroutine etpoi(x,y,h) bind(C, name="etpoi_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x, y, h

      integer i,ndim,ndim2,ntt
! ETPOI TAKES THE POISSON BRACKET OF TWO FUNCTIONS
      parameter (ndim2=6)
      parameter (ndim=3)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer t1,t2,t3

      call etall1(t1)
      call etall1(t2)
      call etall1(t3)

      do i=1,nd

        call dader(2*i-1,x,t1)
        call dader(2*i,y,t2)
        call damul(t1,t2,t1)

        call dalin(t1,1.d0,t3,1.d0,t3)
        call dader(2*i-1,y,t1)
        call dader(2*i,x,t2)
        call damul(t1,t2,t1)

        call dalin(t1,-1.d0,t3,1.d0,t3)

      enddo

      call dacop(t3,h)

      call dadal1(t1)
      call dadal1(t2)
      call dadal1(t3)
      return
      end subroutine

      subroutine exp1d(h,x,y,eps,non) bind(C, name="exp1d_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h, x, y, non
      real(C_DOUBLE) eps

      integer ndim,ndim2,ntt
! WRAPPING EXPFLO
      parameter (ndim2=6)
      parameter (ndim=3)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer idpr
      common /printing/ idpr

      integer v(ndim2)

      call etall(v,nd2)
      call difd(h,v,-1.d0)
      call expflo(v,x,y,eps,non)

      call dadal(v,nd2)

      return
      end subroutine

      subroutine expnd2(h,x,w,eps,nrmax) bind(C, name="expnd2_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h, x(*), w(*), nrmax
      real(C_DOUBLE) eps

      integer j,ndim,ndim2,ntt
! WRAPPING EXPFLOD USING EXP1D
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      integer b0,v(ndim2)
!
!
      call etallnom1(b0,'B0        ')
      call etallnom(v,nd2  ,'V         ')

      call dacopd(x,v)
      do j=1,nd2
        call exp1d(h,v(j),b0,eps,nrmax)
        call dacop(b0,v(j))
      enddo
      call dacopd(v,w)
      call dadal(v,nd2)
      call dadal1(b0)
      return
      end subroutine

      subroutine flofacg(xy,h,epsone) bind(C, name="flofacg_")
      use iso_c_binding
      implicit none
      integer(C_LONG) xy(*), h(*), epsone

      integer i,k,kk,ndim,ndim2,nrmax,ntt
      double precision eps,r,xn,xnbefore,xnorm,xnorm1,xx
! GENERAL ONE EXPONENT FACTORIZATION
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer idpr
      common /printing/ idpr
      logical more
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      double precision xintex
      common /integratedex/ xintex(0:20)
      integer x(ndim2)
      integer v(ndim2),w(ndim2),t(ndim2), z(ndim2)
      integer jj(ntt)
      jj(1)=1
!
      call etallnom(v,nd2  ,'V         ')
      call etallnom(w,nd2  ,'W         ')
      call etallnom(t,nd2  ,'T         ')
      call etallnom(x,nd2  ,'Z         ')
      call etallnom(z,nd2  ,'Z         ')

      call etini(v)
      call daclrd(w)
      xnorm1=0.d0
      do i=1,nd2
        call daabs(xy(i),r)
        xnorm1=xnorm1+r
      enddo
      xnbefore=1.d36
      more=.false.
      eps=1.e-9
      nrmax=1000
      xn=10000.d0
      do k=1,nrmax
        call dacmud(h,-1.d0,t)
        call expflod(t,xy,x,eps,nrmax)
        call dalind(x,1.d0,v,-1.d0,t)
! write(20,*) "$$$$$$$$$$$$$$",k,"$$$$$$$$$$$$$$$$$$$$"
! call daprid(t,1,1,20)
        if(xn.lt.epsone) then
          if(idpr.ge.0) write(6,*) "xn quadratic",xn
          call daflod(t,t,w)
          call dalind(t,1.d0,w,-0.5d0,t)
          call dacopd(t,z)
          call dacopd(t,w)
!  second order in W
          call etcom(h,w,x)
          call etcom(x,w,x)
!  END OF  order in W

          do kk=1,10
            call etcom(h,w,w)
            call dalind(z,1.d0,w,xintex(kk),z)
          enddo
          call dacopd(z,t)
          xx=1.d0/12.d0
          call dalind(x,xx,h,1.d0,h)
        endif

        call dalind(t,1.d0,h,1.d0,h)
        xnorm=0.d0
        do i=1,nd2
          call daabs(t(i),r)
          xnorm=xnorm+r
        enddo
        xn=xnorm/xnorm1
        if(xn.ge.epsone.and.(idpr.ge.0)) write(6,*)" xn linear ",xn
        if(xn.lt.eps.or.more) then
          more=.true.
          if(xn.ge.xnbefore) goto 1000
          xnbefore=xn
        endif
      enddo
1000  if(idpr.ge.0) write(6,*) " iteration " , k
      call dadal(x,nd2)
      call dadal(w,nd2)
      call dadal(v,nd2)
      call dadal(t,nd2)
      call dadal(z,nd2)
      return
      end subroutine

      subroutine flofac(xy,x,h)
      implicit none
      integer k,ndim,ndim2,ntt
! GENERAL DRAGT-FINN FACTORIZATION
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer xy(*),x(*),h(*)
      integer v(ndim2),w(ndim2)
!
      call etallnom(v,nd2  ,'V         ')
      call etallnom(w,nd2  ,'W         ')

      call dacopd(xy,x)
      call dacopd(x,v)
      call daclrd(w)
      call danot(1)
      call etinv(v,w)
      call danot(no)
      call etcct(x,w,v)
      call danot(1)
      call dacopd(xy,x)
      call danot(no)
      call dacopd(v,w)
      call daclrd(h)
      do k=2,no
        call taked(w,k,v)
        call dalind(v,1.d0,h,1.d0,h)
        call facflod(h,w,v,k,k,-1.d0,-1)
        call dacopd(v,w)
      enddo
      call dadal(w,nd2)
      call dadal(v,nd2)
      return
      end subroutine

      subroutine liefact(xy,x,h) bind(C, name="liefact_")
      use iso_c_binding
      implicit none
      integer(C_LONG) xy(*), x(*), h

      integer ndim,ndim2,ntt
! SYMPLECTIC DRAGT-FINN FACTORIZATION WRAPPING FLOFAC
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      integer v(ndim2)

      call etall(v,nd2)

      call flofac(xy,x,v)
      call intd(v,h,-1.d0)
!
      call dadal(v,nd2)

      return
      end subroutine

      logical(C_BOOL) function mapnorm(x,ft,a2,a1,xy,h,nord)            &
     &                         bind(C, name="mapnorm_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*), a1(*), a2(*), ft, xy(*), h, nord

      integer isi,ndim,ndim2,ntt
!--NORMALIZATION ROUTINES OF LIELIB
!- WRAPPING MAPNORMF
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer hf(ndim2),ftf(ndim2)
      logical*1 mapnormf

      call etall(ftf,nd2)
      call etall(hf,nd2)
      isi=0
      mapnorm = mapnormf(x,ftf,a2,a1,xy,hf,nord,isi)
      call intd(hf,h,-1.d0)
      call intd(ftf,ft,-1.d0)
      call dadal(ftf,nd2)
      call dadal(hf,nd2)

      return
      end function

      subroutine gettura(psq,radsq) bind(C, name="gettura_")
      use iso_c_binding
      implicit none
      integer(C_LONG) ndim
      parameter (ndim=3)
      real(C_DOUBLE) psq(ndim),radsq(ndim)

      integer ik,ndim2,ntt
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      double precision ps,rads
      common /tunerad/ ps(ndim),rads(ndim)

      do ik=1,nd
        psq(ik)=ps(ik)
        radsq(ik)=rads(ik)
      enddo

      return
      end subroutine

      subroutine setidpr(idprint,nplan)
      implicit none
      integer idprint,ik,ndim,ndim2,nplan
      parameter (ndim=3)
      parameter (ndim2=6)
      dimension nplan(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer idpr
      common /printing/ idpr
      integer nplane
      double precision epsplane,xplane
      common /choice/ xplane(ndim),epsplane,nplane(ndim)

      do ik=1,nd
        nplane(ik)=nplan(ik)
      enddo
      idpr=idprint

      return
      end subroutine

      subroutine idprset(idprint) bind(C, name="idprset_")
      use iso_c_binding
      implicit none
      integer(C_LONG) idprint

      integer ndim,ndim2
      parameter (ndim=3)
      parameter (ndim2=6)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer idpr
      common /printing/ idpr
      integer nplane
      double precision epsplane,xplane
      common /choice/ xplane(ndim),epsplane,nplane(ndim)

      idpr=idprint

      return
      end subroutine

      logical(C_BOOL) function mapnormf(x,ft,a2,a1,xy,h,nord,isi)       &
     &                         bind(C, name="mapnormf_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*), a1(*), a2(*), ft(*), xy(*), h(*), nord, isi

      integer ij,ndim,ndim2,ntt
      double precision angle,p,rad,st,x2pi,x2pii
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension angle(ndim),st(ndim),p(ndim),rad(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer itu
      common /tunedef/itu
      integer idpr
      common /printing/ idpr
      integer iflow,jtune
      common /vecflow/ iflow,jtune
      double precision ps,rads
      common /tunerad/ ps(ndim),rads(ndim)
      integer a1i(ndim2),a2i(ndim2)
      logical*1 midbflo
!
      call etallnom(a1i,nd2  ,'A1I       ')
      call etallnom(a2i,nd2  ,'A2I       ')
!     frank/etienne
      do itu=1,ndim
        angle(itu)=0.d0
        p(itu)=0.d0
        st(itu)=0.d0
        rad(itu)=0.d0
        ps(itu)=0.d0
        rads(itu)=0.d0
      enddo
      jtune=isi
      x2pii=1.d0/datan(1.d0)/8.d0
      x2pi=datan(1.d0)*8.d0
      call dacopd(x,xy)
! go to fix point in the parameters + pt to order nord>=1
      call gofix(xy,a1,a1i,nord)
      call simil(a1i,xy,a1,xy)
! linear part
      mapnormf = midbflo(xy,a2,a2i,angle,rad,st)
      do ij=1,nd-ndc
        p(ij)=angle(ij)*(st(ij)*(x2pii-1.d0)+1.d0)
      enddo
      if(ndc.eq.1) p(nd)=angle(nd)
      if(idpr.ge.0) then
        write(6,*) 'tune    ',(p(ij),ij=1,nd)
        write(6,*) 'damping ', (rad(ij),ij=1,nd)
      endif
      do ij=1,nd       !  -ndc    Frank
        ps(ij)=p(ij)
        rads(ij)=rad(ij)
      enddo
      call initpert(st,angle,rad)
      call simil(a2i,xy,a2,xy)
      call dacopd(xy,a2i)
!        write(6,*) 'Entering orderflo'
      call orderflo(h,ft,xy,angle,rad)
      do ij=1,nd-ndc
        p(ij)=angle(ij)
        if(angle(ij).gt.x2pi/2.d0.and.st(ij).gt.0.d0.and.itu.eq.1)then
          p(ij)=angle(ij)-x2pi
          write(6,*) ij,' TH TUNE MODIFIED IN H2 TO ',p(ij)/x2pi
        endif
      enddo
      call h2pluflo(h,p,rad)
!      CALL TAKED(A2I,1,XY)
      call taked(a2i,1,a1i)
      call etcct(xy,a1i,xy)

      call dadal(a2i,nd2)
      call dadal(a1i,nd2)
      return
      end function

      subroutine gofix(xy,a1,a1i,nord) bind(C, name="gofix_")
      use iso_c_binding
      implicit none
      integer(C_LONG) xy(*), a1(*), a1i(*), nord

      integer i,ndim,ndim2,ntt
      double precision xic
! GETTING TO THE FIXED POINT AND CHANGING TIME APPROPRIATELY IN THE
! COASTING BEAM CASE

!****************************************************************
! X = A1 XY A1I WHERE X IS TO THE FIXED POINT TO ORDER NORD
! for ndpt not zero, works in all cases. (coasting beam: eigenvalue
!1 in Jordan form)
!****************************************************************
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt

      integer x(ndim2),w(ndim2),v(ndim2),rel(ndim2)
!
      call etallnom(x,nd2  ,  'X         ')
      call etallnom(w,nd2  ,  'W         ')
      call etallnom(v,nd2  ,  'V         ')
      call etallnom(rel,nd2  ,'REL       ')

! COMPUTATION OF A1 AND A1I USING DAINV
      call etini(rel)

      call danot(nord)

      call etini(v)

      do i=1,nd2-ndc2
        call dacop(xy(i),x(i))
        call dalin(x(i),1.d0,rel(i),-1.d0,v(i))
      enddo
      call etinv(v,w)
      call daclrd(x)
      if(ndc.eq.1) then
        call davar(x(ndpt),0.d0,ndpt)
      endif
      call etcct(w,x,v)
      if(ndc.eq.1) then
        call daclr(v(nd2))
        call daclr(v(nd2-ndc))
      endif
      call dalind(rel,1.d0,v,1.d0,a1)
      call dalind(rel,1.d0,v,-1.d0,a1i)

      if(ndpt.ne.0) then

!  CORRECTIONS
        call daclrd(w)
        call daclrd(v)
        call daclrd(x)

        do i=1,nd2-ndc2
          call dalin(a1(i),1.d0,rel(i),-1.d0,w(i))
        enddo

!      COMPUTE Deta/Ddelta
        call dacopd(w,a1)

        do i=1,nd2-ndc2
          call dader(ndpt,w(i),w(i))
        enddo
!      COMPUTE J*Deta/dDELTA

        do i=1,nd-ndc
          call dacmu(w(2*i),1.d0,v(2*i-1) )
          call dacmu(w(2*i-1),-1.d0,v(2*i) )
        enddo

        xic=(-1)**(ndt)

        do i=1,nd2-ndc2
          call damul(v(i),rel(i),x(1))
          call daadd(x(1),w(ndt),w(ndt))
          call dacop(a1(i),w(i))
        enddo
        call dacmu(w(ndt),xic,w(ndt))

        call expflod(w,rel,a1,1.d-7,10000)
! END OF  CORRECTIONS

        call etinv(a1,a1i)

      endif

      call danot(no)

      call dadal(rel,nd2)
      call dadal(v,nd2)
      call dadal(w,nd2)
      call dadal(x,nd2)
      return
      end subroutine

      double precision function transver(j)
      implicit none
      integer i,ic,ndim
! USED IN A DACFU CALL OF GOFIX
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      integer J(NTT)
      integer j(*)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt

      transver=1.d0
      ic=0
      do i=1,nd2-ndc2
        ic=ic+j(i)
      enddo
      if(ic.ne.1) transver=0.d0
      return
      end function

      subroutine orderflo(h,ft,x,ang,ra)
      implicit none
      integer k,ndim,ndim2,ntt
      double precision ang,ra
!-   NONLINEAR NORMALIZATION PIECE OF MAPNORMF
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer idpr
      common /printing/ idpr
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      dimension ang(ndim),ra(ndim)
      integer x(*),h(*),ft(*)
      integer w(ndim2),v(ndim2),rel(ndim2)
      integer roi(ndim2)
      integer b1(ndim2),b5(ndim2),b6(ndim2),b9(ndim2)
!
      call etallnom(w,nd2  ,'W         ')
      call etallnom(v,nd2  ,'V         ')
      call etallnom(rel,nd2  ,'REL       ')
      call etallnom(roi,nd2  ,'ROI       ')
      call etallnom(b1,nd2  ,'B1        ')
      call etallnom(b5,nd2  ,'B5        ')
      call etallnom(b6,nd2  ,'B6        ')
      call etallnom(b9,nd2  ,'B9        ')
      call rotiflo(roi,ang,ra)
      call etini(rel)
      call daclrd(h)
      call daclrd(ft)
      call etcct(x,roi,x)
      do k=2,no
! IF K>2 V = H(K)^-1 X(K)
        call facflod(h,x,v,2,k-1,-1.d0,-1)
! EXTRACTING K TH DEGREE OF V ----> W
        call taked(v,k,w)
!  write(16,*) "$$$$$$$$  K  $$$$$$$$$$", k
! W = EXP(B5) + ...
        call dacopd(w,b5)
!      CALL INTD(W,B5,-1.D0)
! B5 ON EXIT IS THE NEW CONTRIBUTION TO H
! B6 IS THE NEW CONTRIBUTION TO FT
        call nuanaflo(b5,b6)
        call dalind(b5,1.d0,h,1.d0,b1)
        call dacopd(b1,h)
! EXP(B9) = EXP( : ROTI B6 :)
        call trxflo(b6,b9,roi)

! V = EXP(-B6) REL
        call facflod(b6,rel,v,k,k,-1.d0,1)
! W = V o X
        call etcct(v,x,w)
        if(idpr.ge.0) then
          write(6,*) ' ORDERFLO K= ', k
        endif
! X = EXP(B9) W
        call facflod(b9,w,x,k,k,1.d0,1)
! B6 IS THE NEW CONTRIBUTION TO FT
        call dalind(b6,1.d0,ft,1.d0,b1)
        call dacopd(b1,ft)
      enddo
      call dadal(b9,nd2)
      call dadal(b6,nd2)
      call dadal(b5,nd2)
      call dadal(b1,nd2)
      call dadal(roi,nd2)
      call dadal(rel,nd2)
      call dadal(v,nd2)
      call dadal(w,nd2)
      return
      end subroutine

      subroutine nuanaflo(h,ft)
      implicit none
      integer i,ndim,ndim2,ntt
      double precision dfilt,filt,xgam,xgbm
! RESONANCE DENOMINATOR OPERATOR (1-R^-1)^-1
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer iflow,jtune
      common /vecflow/ iflow,jtune
      external xgam,xgbm,dfilt,filt
      integer h(*),ft(*),br(ndim2),bi(ndim2),c(ndim2),ci(ndim2)
      integer t1(2),t2(2)

      call etall(br,nd2)
      call etall(bi,nd2)
      call etall(c,nd2)
      call etall(ci,nd2)

      call ctorflo(h,br,bi)

! FILTERING RESONANCES AND TUNE SHIFTS
! ASSUMING REALITY I.E. B(2*I-1)=CMPCJG(B(2*I))

      do i=1,nd2
        iflow=i
        call dacfu(br(i),filt,c(i))
        call dacfu(bi(i),filt,ci(i))
      enddo
      call rtocflo(c,ci,h)

      do i=1,nd2

        iflow=i
        call dacfu(br(i),dfilt,br(i))
        call dacfu(bi(i),dfilt,bi(i))
      enddo
!  NOW WE MUST REORDER C AND CI TO SEPARATE THE REAL AND IMAGINARY PART
! THIS IS NOT NECESSARY WITH :H: OPERATORS

      do i=1,nd2
        t1(1)=br(i)
        t1(2)=bi(i)
        t2(1)=c(i)
        t2(2)=ci(i)
        iflow=i
        call comcfu(t1,xgam,xgbm,t2)
      enddo

      call rtocflo(c,ci,ft)

      call dadal(br,nd2)
      call dadal(bi,nd2)
      call dadal(c,nd2)
      call dadal(ci,nd2)

      return
      end subroutine

      real(C_DOUBLE) function xgam(j) bind(C, name="xgam_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ic,ij,ik,ndim,ndim2
      double precision ad,ans,as,ex,exh
! XGAM AND XGBM ARE THE EIGENVALUES OF THE OPERATOR NEWANAFLO
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
      integer iflow,jtune
      common /vecflow/ iflow,jtune
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
!      integer J(NTT),JJ(NDIM),JP(NDIM)
      integer jj(ndim),jp(ndim)
      xgam=0.d0
      ad=0.d0
      as=0.d0
      ic=0
      do i=1,nd-ndc
        ik=2*i-1
        ij=2*i
        jp(i)=j(ik)+j(ij)
        jj(i)=j(ik)-j(ij)
        if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
          jp(i)=jp(i)-1
        endif
        ic=ic+iabs(jj(i))
      enddo

      do i=1,nd-ndc
        ad=dsta(i)*dble(jj(i))*angle(i)-dble(jp(i))*rad(i)+ad
        as=sta(i)*dble(jj(i))*angle(i)+as
      enddo

      exh=dexp(ad/2.d0)
      ex=exh**2
      ans=4.d0*ex*(dsinh(ad/2.d0)**2+dsin(as/2.d0)**2)
      xgam=2.d0*(-exh*dsinh(ad/2.d0)+ex*dsin(as/2.d0)**2)/ans

      return
      end function

      real(C_DOUBLE) function xgbm(j) bind(C, name="xgbm_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ic,ij,ik,ndim,ndim2
      double precision ad,ans,as,ex,exh
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
      integer iflow,jtune
      common /vecflow/ iflow,jtune
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
!      integer J(NTT),JJ(NDIM),JP(NDIM)
      integer jj(ndim),jp(ndim)
      xgbm=0.d0
      ad=0.d0
      as=0.d0
      ic=0
      do i=1,nd-ndc
        ik=2*i-1
        ij=2*i
        jp(i)=j(ik)+j(ij)
        jj(i)=j(ik)-j(ij)
        if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
          jp(i)=jp(i)-1
        endif
        ic=ic+iabs(jj(i))
      enddo

      do i=1,nd-ndc
        ad=dsta(i)*dble(jj(i))*angle(i)-dble(jp(i))*rad(i)+ad
        as=sta(i)*dble(jj(i))*angle(i)+as
      enddo

      exh=dexp(ad/2.d0)
      ex=exh**2
      ans=4.d0*ex*(dsinh(ad/2.d0)**2+dsin(as/2.d0)**2)
      xgbm=dsin(as)*ex/ans

      return
      end function

      real(C_DOUBLE) function filt(j) bind(C, name="filt_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ic,ic1,ic2,ij,ik,ji,ndim,ndim2,nreso
!  PROJECTION FUNCTIONS ON THE KERNEL ANMD RANGE OF (1-R^-1)
!-  THE KERNEL OF (1-R^-1)
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
      parameter (nreso=20)
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer mx,nres
      common /reson/mx(ndim,nreso),nres
      integer iflow,jtune
      common /vecflow/ iflow,jtune
!      integer J(NTT),JJ(NDIM)
      integer jj(ndim)

      filt=1.d0

      ic=0
      do i=1,nd-ndc
        ik=2*i-1
        ij=2*i
        jj(i)=j(ik)-j(ij)
        if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
        endif
        ic=ic+iabs(jj(i))
      enddo

      if(ic.eq.0.and.jtune.eq.0) return

      do i=1,nres
        ic1=1
        ic2=1
        do ji=1,nd-ndc
          if(mx(ji,i).ne.jj(ji)) ic1=0
          if(mx(ji,i).ne.-jj(ji)) ic2=0
          if(ic1.eq.0.and.ic2.eq.0) goto 3
        enddo
        return
 3      continue
      enddo

      filt=0.d0
      return
      end function

      real(C_DOUBLE) function dfilt(j) bind(C, name="dfilt_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer ndim,ndim2,nreso
      double precision fil,filt
!-  THE RANGE OF (1-R^-1)^1
!- CALLS FILT AND EXCHANGES 1 INTO 0 AND 0 INTO 1.
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
      parameter (nreso=20)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer mx,nres
      common /reson/mx(ndim,nreso),nres
      external filt
!      integer J(NTT)

      fil=filt(j)
      if(fil.gt.0.5d0) then
        dfilt=0.d0
      else
        dfilt=1.d0
      endif
      return
      end function

      subroutine dhdjflo(h,t)
      implicit none
      integer i,ndim,ndim2,ntt
      double precision coe,x2pi
! CONVENIENT TUNE SHIFT FINDED FOR SYMPLECTIC CASE (NU,DL)(H)=T
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer h(*),t(*)

      integer b1(ndim2),b2(ndim2),bb1,bb2
!
      call etall(b1,nd2)
      call etall(b2,nd2)
      call etall1(bb1)
      call etall1(bb2)

      x2pi=datan(1.d0)*8.d0
      call ctorflo(h,b1,b2)
      coe=1.d0/x2pi

      do i=1,nd-ndc
        call datra(2*i,b2(2*i),bb1)
        call dacmu(bb1,coe,t(i+nd))
        call dacop(t(i+nd),bb1)
        call daclr(bb2)
        call rtoc(bb1,bb2,bb1)
        call dacop(bb1,t(i))
      enddo

      if(ndpt.ne.0) then
        call dacop(h(ndt),t(nd))
        call dacop(b1(ndt),t(nd2))
      endif

      call dadal1(bb2)
      call dadal1(bb1)
      call dadal(b2,nd2)
      call dadal(b1,nd2)
      return
      end subroutine

      subroutine dhdj(h,t) bind(C, name="dhdj_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h, t(*)

      integer i,ndim,ndim2,ntt
      double precision coe,x2pi
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt

      integer b1,b2,bb1,bb2
!
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')
      call etallnom1(bb1,'BB1       ')
      call etallnom1(bb2,'BB2       ')

      x2pi=datan(1.d0)*8.d0
      call ctor(h,b1,b2)
      coe=-2.d0/x2pi
      do i=1,nd-ndc
        call dader(2*i-1,b1,b2)
        call datra(2*i,b2,bb2)
        call dacmu(bb2,coe,t(i+nd))
        call dacop(t(i+nd),bb2)
        call daclr(b2)
        call rtoc(bb2,b2,bb1)
        call dacop(bb1,t(i))
      enddo

      if(ndpt.eq.nd2) then
        call dader(ndpt,h,t(nd))
        call dader(ndpt,b1,t(nd2))
        call dacmu(t(nd),-1.d0,t(nd))
        call dacmu(t(nd2),-1.d0,t(nd2))
      endif
      if(ndt.eq.nd2) then
        call dader(ndpt,h,t(nd))
        call dader(ndpt,b1,t(nd2))
      endif
      call dadal1(bb2)
      call dadal1(bb1)
      call dadal1(b2)
      call dadal1(b1)
      return
      end subroutine

      subroutine h2pluflo(h,ang,ra)
      implicit none
      integer i,j,ndim,ndim2,ntt
      double precision ang,r1,r2,ra,st
! POKES IN \VEC{H}  ANGLES AND DAMPING COEFFFICIENTS
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      dimension ang(ndim),st(ndim),ra(ndim),j(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer h(*)
!*DAEXT(NO,NV) H

      do i=1,nd
        st(i)=2.d0*sta(i)-1.d0
      enddo

      do i=1,ntt
        j(i)=0
      enddo

      do i=1,nd-ndc
        j(2*i-1)=1
        r1=-ang(i)
!-----
        call dapok(h(2*i),j,r1)

        r2=ra(i)
        call dapok(h(2*i-1),j,r2)
        j(2*i-1)=0

        j(2*i)=1
        r1=ang(i)*st(i)
        call dapok(h(2*i-1),j,r1)
        call dapok(h(2*i),j,r2)
        j(2*i)=0

      enddo

      if(ndpt.eq.nd2-1) then
        j(ndpt)=1
        call dapok(h(ndt),j,ang(nd))
      elseif(ndpt.eq.nd2) then
        j(ndpt)=1
        call dapok(h(ndt),j,-ang(nd))
      endif
      return
      end subroutine

      subroutine rotflo(ro,ang,ra)
      implicit none
      integer i,ndim,ndim2,ntt
      double precision ang,ch,co,ra,sh,si,sim,xx
! CREATES R AND R^-1 USING THE EXISTING ANGLES AND DAMPING
! COULD BE REPLACED BY A CALL H2PLUFLO FOLLOWED BY EXPFLOD
! CREATES R
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      dimension co(ndim),si(ndim),ang(ndim),ra(ndim)
      integer j(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer ro(*)
      call daclrd(ro)
      do i=1,nd-ndc
        xx=dexp(ra(i))
        if(ista(i).eq.0) then
          call hyper(ang(i),ch,sh)
          co(i)=ch*xx
          si(i)=-sh*xx
        else
          co(i)=dcos(ang(i))*xx
          si(i)=dsin(ang(i))*xx
        endif
      enddo
      do i=1,nd-ndc
        if(ista(i).eq.0)then
          sim=si(i)
        else
          sim=-si(i)
        endif
        j(2*i-1)=1
        call dapok(ro(2*i-1),j,co(i))
        call dapok(ro(2*i),j,sim)
        j(2*i-1)=0
        j(2*i)=1
        call dapok(ro(2*i),j,co(i))
        call dapok(ro(2*i-1),j,si(i))
        j(2*i)=0
      enddo

      if(ndc.eq.1) then
        j(ndt)=1
        call dapok(ro(ndt),j,1.d0)
        call dapok(ro(ndpt),j,0.d0)
        j(ndt)=0
        j(ndpt)=1
        call dapok(ro(ndt),j,ang(nd))
        call dapok(ro(ndpt),j,1.d0)
        j(ndpt)=0
      endif

      return
      end subroutine

      subroutine rotiflo(roi,ang,ra)
      implicit none
      integer i,ndim,ndim2,ntt
      double precision ang,ch,co,ra,sh,si,sim,simv,xx
! CREATES  R^-1
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      dimension co(ndim),si(ndim),ang(ndim),ra(ndim)
      integer j(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer roi(*)

      do i=1,10
        j(i)=0
      enddo

      call daclrd(roi)
      do i=1,nd-ndc
        xx=dexp(-ra(i))
        if(ista(i).eq.0) then
          call hyper(ang(i),ch,sh)
          co(i)=ch*xx
          si(i)=-sh*xx
        else
          co(i)=dcos(ang(i))*xx
          si(i)=dsin(ang(i))*xx
        endif
      enddo
      do i=1,nd-ndc
        if(ista(i).eq.0)then
          sim=si(i)
        else
          sim=-si(i)
        endif
        j(2*i-1)=1
        call dapok(roi(2*i-1),j,co(i))
        simv=-sim
        call dapok(roi(2*i),j,simv)
        j(2*i-1)=0
        j(2*i)=1
        simv=-si(i)
        call dapok(roi(2*i),j,co(i))
        call dapok(roi(2*i-1),j,simv)
        j(2*i)=0
      enddo

      if(ndc.eq.1) then
        j(ndt)=1
        call dapok(roi(ndt),j,1.d0)
        call dapok(roi(ndpt),j,0.d0)
        j(ndt)=0
        j(ndpt)=1
        call dapok(roi(ndt),j,-ang(nd))
        call dapok(roi(ndpt),j,1.d0)
        j(ndpt)=0
      endif

      return
      end subroutine

      subroutine hyper(a,ch,sh)
      implicit none
      double precision a,ch,sh,x,xi
!   USED IN ROTIFLO AND ROTFLO
      x=dexp(a)
      xi=1.d0/x
      ch=(x+xi)/2.d0
      sh=(x-xi)/2.d0
      return
      end subroutine

      subroutine ctor(c1,r2,i2) bind(C, name="ctor_")
      use iso_c_binding
      implicit none
      integer(C_LONG) c1, r2, i2

      integer ndim2,ntt
! CHANGES OF BASIS
!   C1------> R2+I R1
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer b1,b2,x(ndim2)
!
!
      call etallnom1(b1,'B1        ')
      call etallnom1(b2,'B2        ')
      call etallnom(x,nd2  ,'X         ')

      call ctoi(c1,b1)
      call etcjg(x)
      call trx(b1,b2,x)
      call dalin(b1,.5d0,b2,.5d0,r2)
      call dalin(b1,.5d0,b2,-.5d0,i2)
      call dadal(x,nd2)
      call dadal1(b2)
      call dadal1(b1)
      return
      end subroutine

      subroutine rtoc(r1,i1,c2) bind(C, name="rtoc_")
      use iso_c_binding
      implicit none
      integer(C_LONG) r1, i1, c2

      integer ndim2,ntt
!  INVERSE OF CTOR
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2

      integer b1
!
      call etallnom1(b1,'B1        ')

      call daadd(r1,i1,b1)
      call itoc(b1,c2)
      call dadal1(b1)
      return
      end subroutine

      subroutine ctorflo(c,dr,di)
      implicit none
      integer ndim,ndim2,ntt
! FLOW CTOR
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer dr(*),di(*),c(*)

      call ctord(c,dr,di)
      call resvec(dr,di,dr,di)

      return
      end subroutine

      subroutine rtocflo(dr,di,c)
      implicit none
      integer ndim,ndim2,ntt
! FLOW RTOC
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer dr(*),di(*),c(*),er(ndim2),ei(ndim2)

      call etall(er,nd2)
      call etall(ei,nd2)

      call reelflo(dr,di,er,ei)
      call rtocd(er,ei,c)

      call dadal(er,nd2)
      call dadal(ei,nd2)

      return
      end subroutine

      subroutine ctord(c,cr,ci)
      implicit none
      integer i,ndim,ndim2,ntt
! ROUTINES USED IN THE INTERMEDIATE STEPS OF CTORFLO AND RTOCFLO
! SAME AS CTOR  OVER ARRAYS CONTAINING ND2 COMPONENTS
! ROUTINE USEFUL IN INTERMEDIATE FLOW CHANGE OF BASIS
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer c(*),ci(*),cr(*)
      do i=1,nd2
        call ctor(c(i),cr(i),ci(i))
      enddo
      return
      end subroutine

      subroutine rtocd(cr,ci,c)
      implicit none
      integer i,ndim,ndim2,ntt
!  INVERSE OF CTORD
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer c(*),ci(*),cr(*)
      do i=1,nd2
        call rtoc(cr(i),ci(i),c(i))
      enddo
      return
      end subroutine

      subroutine resvec(cr,ci,dr,di)
      implicit none
      integer i,ndim,ndim2,ntt
! DOES THE SPINOR PART IN CTORFLO
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer dr(*),di(*),ci(*),cr(*),tr(2),ti(2)

      call etall(tr,2)
      call etall(ti,2)

      do i=1,nd-ndc
        if(ista(i).eq.1) then
          call dasub(cr(2*i-1),ci(2*i),tr(1))
          call daadd(ci(2*i-1),cr(2*i),ti(1))
          call daadd(cr(2*i-1),ci(2*i),tr(2))
          call dasub(ci(2*i-1),cr(2*i),ti(2))
          call dacop(tr(1),dr(2*i-1))
          call dacop(tr(2),dr(2*i))
          call dacop(ti(1),di(2*i-1))
          call dacop(ti(2),di(2*i))
        else
          call daadd(cr(2*i-1),cr(2*i),tr(1))
          call daadd(ci(2*i-1),ci(2*i),ti(1))
          call dasub(cr(2*i-1),cr(2*i),tr(2))
          call dasub(ci(2*i-1),ci(2*i),ti(2))
          call dacop(tr(1),dr(2*i-1))
          call dacop(tr(2),dr(2*i))
          call dacop(ti(1),di(2*i-1))
          call dacop(ti(2),di(2*i))
        endif
      enddo

      do i=nd2-ndc2+1,nd2
        call dacop(cr(i),dr(i))
        call dacop(ci(i),di(i))
      enddo

      call dadal(tr,2)
      call dadal(ti,2)
      return
      end subroutine

      subroutine reelflo(c,ci,f,fi)
      implicit none
      integer i,ndim,ndim2,ntt
! DOES THE SPINOR PART IN RTOCFLO
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer c(*),ci(*),f(*),fi(*),e(ndim2),ei(ndim2)

      call etall(e,nd2)
      call etall(ei,nd2)

      do i=1,nd-ndc
        call dalin(c(2*i-1),0.5d0,c(2*i),0.5d0,e(2*i-1))
        call dalin(ci(2*i-1),0.5d0,ci(2*i),0.5d0,ei(2*i-1))
        if(ista(i).eq.1) then
          call dalin(ci(2*i-1),0.5d0,ci(2*i),-0.5d0,e(2*i))
          call dalin(c(2*i-1),-0.5d0,c(2*i),0.5d0,ei(2*i))
        else
          call dalin(ci(2*i-1),0.5d0,ci(2*i),-0.5d0,ei(2*i))
          call dalin(c(2*i-1),0.5d0,c(2*i),-0.5d0,e(2*i))
        endif
      enddo

      do i=nd2-ndc2+1,nd2
        call dacop(c(i),e(i))
        call dacop(ci(i),ei(i))
      enddo

      call dacopd(e,f)
      call dacopd(ei,fi)

      call dadal(e,nd2)
      call dadal(ei,nd2)
      return
      end subroutine

      subroutine compcjg(cr,ci,dr,di)
      implicit none
      integer ndim,ndim2,ntt
! TAKES THE COMPLEX CONJUGATE IN RESONANCE BASIS OF A POLYNOMIAL
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer dr,di,ci,cr,x(ndim2)

      call etall(x,nd2)

      call etcjg(x)
      call trx(cr,dr,x)
      call trx(ci,di,x)
      call dacmu(di,-1.d0,di)

      call dadal(x,nd2)
      return
      end subroutine

      logical(C_BOOL) function midbflo(c,a2,a2i,q,a,st)                 &
     &                         bind(C, name="midbflo_")
      use iso_c_binding
      implicit none
      integer(C_LONG) ndim
      parameter (ndim=3)
      integer(C_LONG) c(*), a2(*), a2i(*)
      real(C_DOUBLE) q(ndim), a(ndim), st(ndim)

      integer i,j,ndim2,ntt
      double precision ch,cm,cr,r,sa,sai,shm,                           &
     &x2pi
! LINEAR EXACT NORMALIZATION USING EIGENVALUE PACKAGE OF NERI
      parameter (ntt=40)
      parameter (ndim2=6)
      integer jx(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      dimension cr(ndim2,ndim2)
      dimension sa(ndim2,ndim2),sai(ndim2,ndim2),cm(ndim2,ndim2)
      logical*1 mapflol

!*DAEXT(NO,NV) C(NDIM2),A2(NDIM2),A2I(NDIM2)
      x2pi=datan(1.d0)*8.d0

      do i=1,ntt
        jx(i)=0
      enddo

!     frank/etienne
      do i=1,ndim
        st(i)=0d0
        q(i)=0d0
        a(i)=0d0
      enddo
!     frank/etienne
      do i=1,ndim2
!     frank/etienne
        do j=1,ndim2
          sai(i,j)=0.d0
          sa(i,j)=0.d0
          cm(i,j)=0.d0
          cr(i,j)=0.d0
        enddo
      enddo

      do i=1,nd2
        do j=1,nd2
          jx(j)=1
          call  dapek(c(i),jx,r)
          jx(j)=0
          cm(i,j)=r
        enddo
      enddo

      midbflo = mapflol(sa,sai,cr,cm,st)
      do i=1,nd-ndc
        if(st(i)+0.001.gt.1.d0) then
          a(i)=dsqrt(cr(2*i-1,2*i-1)**2+cr(2*i-1,2*i)**2)
          q(i)=dacos(cr(2*i-1,2*i-1)/a(i))
          a(i)=dlog(a(i))
          if(cr(2*i-1,2*i).lt.0.d0) q(i)=x2pi-q(i)
        else
          a(i)=dsqrt(cr(2*i-1,2*i-1)**2-cr(2*i-1,2*i)**2)
          ch=cr(2*i-1,2*i-1)/a(i)
          shm=cr(2*i-1,2*i)/a(i)
!       CH=CH+DSQRT(CH**2-1.D0)
!       q(i)=DLOG(CH)
          q(i)=-dlog(ch+shm)
!       IF(cr(2*i-1,2*i).gt.0.d0) Q(I)=-Q(I)
          a(i)=dlog(a(i))
        endif
      enddo

      if(ndc.eq.0) then
        if(st(3)+0.001.gt.1.d0.and.nd.eq.3.and.q(nd).gt.0.5d0)          &
     &  q(3)=q(3)-x2pi
      else
        q(nd)=cr(ndt,ndpt)
      endif

      call daclrd(a2)
      call daclrd(a2i)

      do i=1,nd2
        do j=1,nd2
          jx(j)=1
          r=sa(i,j)
          if(r.ne.0.d0)call  dapok(a2(i),jx,r)
          jx(j)=1
          r=sai(i,j)
          if(r.ne.0.d0)call  dapok(a2i(i),jx,r)
          jx(j)=0
        enddo
      enddo

      return
      end function

      logical(C_BOOL) function mapflol(sa,sai,cr,cm,st)                 &
     &                         bind(C, name="mapflol_")
      use iso_c_binding
      implicit none
      integer(C_LONG) ndim, ndim2
      parameter (ndim=3, ndim2=6)
      real(C_DOUBLE) sa(ndim2,ndim2),sai(ndim2,ndim2), cr(ndim2,ndim2), &
     &               cm(ndim2,ndim2), st(ndim)

      integer i,ier,iunst,j,l,n,n1
      double precision ap,ax,                                           &
     &p,rd,rd1,ri,rr,s1,vi,vr,w,x,x2pi,xd,xj,xsu,xx
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
!---- FROM TRACKING CODE
      integer idpr
      common /printing/ idpr
      integer nplane
      double precision epsplane,xplane
      common /choice/ xplane(ndim),epsplane,nplane(ndim)
! ---------------------
      dimension xj(ndim2,ndim2),n(ndim),x(ndim)
      dimension rr(ndim2),ri(ndim2),xx(ndim)                            &
     &,w(ndim2,ndim2)
      dimension vr(ndim2,ndim2),vi(ndim2,ndim2),s1(ndim2,ndim2),p(ndim2)
      logical*1 eig6

      x2pi=datan(1.d0)*8.d0
      n1=0
!     frank/etienne
      do i=1,ndim2
        do j=1,ndim2
          cr(j,i)=cm(i,j)
          xj(i,j)=0.d0
          s1(i,j)=0.d0
        enddo
      enddo
!     frank/etienne
      do i=1,ndim
        n(i)=0
        xj(2*i-1,2*i)=1.d0
        xj(2*i,2*i-1)=-1.d0
      enddo
!     frank/etienne
      do i=1,ndim2
        do j=1,ndim2
          sai(i,j)=0.d0
          w(i,j)=cm(i,j)
        enddo
      enddo
      if(ndc.eq.1) then
        s1(nd2-ndc,nd2-ndc)=1.d0
        s1(nd2,nd2)=1.d0
        sai(nd2-ndc,nd2-ndc)=1.d0
        sai(nd2,nd2)=1.d0
      endif
      call mulnd2(xj,w)
      call mulnd2(cr,w)
      if(idpr.ge.0.or.idpr.eq.-102) then
        write(6,*)'Check of the symplectic condition on the linear part'
        xsu=0.d0
        do i=1,nd2
          write(6,'(6(2x,g23.16))') ( w(i,j), j = 1, nd2 )
          do j=1,nd2
            xsu=xsu+dabs(w(i,j))
          enddo
        enddo
        write(6,*)'deviation for symplecticity ',100.d0*(xsu-nd2)/xsu,  &
     &  ' %'
      endif
      mapflol = eig6(cr,rr,ri,vr,vi)
      if(idpr.ge.0) then
        write(6,*) '   '
        write(6,*) '       Index         Real Part  ',                  &
     &  '       ArcSin(Imaginary Part)/2/pi'
        write(6,*) '   '
        do i=1,nd-ndc
          rd1=dsqrt(rr(2*i-1)**2+ri(2*i-1)**2)
          rd=dsqrt(rr(2*i)**2+ri(2*i)**2)
          write(6,*) 2*i-1,rr(2*i-1),dasin(ri(2*i-1)/rd1)/x2pi
          write(6,*) 2*i,rr(2*i),dasin(ri(2*i)/rd)/x2pi
          write(6,*) ' alphas ', dlog(dsqrt(rd*rd1))
        enddo
        if ( idpr.ge. 0) then
          write(6,*)                                                    &
     &    ' select ',nd-ndc,                                            &
     &    ' eigenplanes (odd integers <0 real axis)'
          read(5,*) (n(i),i=1,nd-ndc)
        else
          n(1) = 1
          n(2) = 3
          n(3) = 5
        endif
      elseif(idpr.eq.-100) then
        do i=1,nd-ndc
          n(i)=nplane(i)
        enddo
      elseif(idpr.eq.-101.or.idpr.eq.-102) then
        do i=1,nd-ndc
          if(ri(2*i).ne.0.d0) then
            n(i)=2*i-1
          else
            n(i)=-2*i+1
          endif
        enddo
      else
        do i=1,nd-ndc
          n(i)=2*i-1
        enddo
      endif
      iunst=0
      do i=1,nd-ndc                  ! Frank NDC  kept
        if(n(i).lt.0) then
          n(i)=-n(i)
          st(i)=0.d0
          iunst=1
        else
          st(i)=1.d0
        endif
        x(i)=0.d0
        xx(i)=1.d0
        do j=1,nd-ndc
          x(i)=vr(2*j-1,n(i))*vi(2*j,n(i))-vr(2*j,n(i))*vi(2*j-1,n(i))+ &
     &    x(i)
        enddo
      enddo

      do i=1,nd-ndc
        if(x(i).lt.0.d0) xx(i)=-1.d0
        x(i)=dsqrt(dabs(x(i)))
      enddo
      do i=1,nd2-ndc2
        do j=1,nd-ndc
          if(st(j)+0.001.gt.1.d0) then
            sai(2*j-1,i)=vr(i,n(j))*xx(j)/x(j)
            sai(2*j,i)=vi(i,n(j))/x(j)
          else
            ax=vr(i,n(j))*xx(j)/x(j)
            ap=vi(i,n(j))/x(j)
            sai(2*j-1,i)=(ax+ap)/dsqrt(2.d0)
            sai(2*j,i)=(ap-ax)/dsqrt(2.d0)
          endif
        enddo
      enddo
      if(idpr.eq.-101.or.idpr.eq.-102) then
        call movearou(sai)
      endif
! adjust sa such that sa(1,2)=0 and sa(3,4)=0. (courant-snyder-edwards-teng
! phase advances)
      if(iunst.ne.1) then
        do i=1,nd-ndc
          p(i)=datan(-sai(2*i-1,2*i)/sai(2*i,2*i))
          s1(2*i-1,2*i-1)=dcos(p(i))
          s1(2*i,2*i)=dcos(p(i))
          s1(2*i-1,2*i)=dsin(p(i))
          s1(2*i,2*i-1)=-dsin(p(i))
        enddo
        call mulnd2(s1,sai)
! adjust sa to have sa(1,1)>0 and sa(3,3)>0 rotate by pi if necessary.
        do i=1,nd-ndc
          xd=1.d0
          if(sai(2*i-1,2*i-1).lt.0.d0) xd=-1.d0
          s1(2*i-1,2*i-1)=xd
          s1(2*i-1,2*i)=0.d0
          s1(2*i,2*i-1)=0.d0
          s1(2*i,2*i)=xd
        enddo
        call mulnd2(s1,sai)
! sa is now uniquely and unambigeously determined.
      endif
      do i=1,nd2
        do l=1,nd2
          sa(i,l)=sai(i,l)
        enddo
      enddo
      call matinv(sai,sa,nd2,6,ier)

      call mulnd2(sai,cm)
      do i=1,nd2
        do j=1,nd2
          cr(i,j)=sa(i,j)
        enddo
      enddo

      call mulnd2(cm,cr)

      return
      end function

      subroutine mulnd2(rt,r)
      implicit none
      integer i,ia,j,ndim,ndim2
      double precision r,rt,rtt
      parameter (ndim2=6)
      parameter (ndim=3)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      dimension rt(ndim2,ndim2),r(ndim2,ndim2),rtt(ndim2,ndim2)
      do i=1,nd2
        do j=1,nd2
          rtt(i,j)=0.d0
        enddo
      enddo
      do i=1,nd2
        do j=1,nd2
          do ia=1,nd2
            rtt(i,ia)=rt(i,j)*r(j,ia)+rtt(i,ia)
          enddo
        enddo
      enddo

      do i=1,nd2
        do j=1,nd2
          r(i,j)=rtt(i,j)
        enddo
      enddo
      return
      end subroutine

      subroutine movearou(rt)
      implicit none
      integer ipause, mypause
      integer i,ic,j,ndim,ndim2
      double precision rt,rto,s,xr,xrold,xy,xyz,xz,xzy,yz
      parameter (ndim2=6)
      parameter (ndim=3)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer idpr
      common /printing/ idpr
      dimension rt(ndim2,ndim2),rto(ndim2,ndim2)
      dimension xy(ndim2,ndim2),xz(ndim2,ndim2),yz(ndim2,ndim2)
      dimension xyz(ndim2,ndim2),xzy(ndim2,ndim2)
      dimension s(ndim2,ndim2)
      do i=1,nd2
        do j=1,nd2
          s(i,j)=0.d0
          s(i,i)=1.d0
          xy(i,j)=0.d0
          xz(i,j)=0.d0
          yz(i,j)=0.d0
          xyz(i,j)=0.d0
          xzy(i,j)=0.d0
        enddo
      enddo

      xy(1,3)=1.d0
      xy(3,1)=1.d0
      xy(2,4)=1.d0
      xy(4,2)=1.d0
      xy(5,5)=1.d0
      xy(6,6)=1.d0

      xz(1,5)=1.d0
      xz(5,1)=1.d0
      xz(2,6)=1.d0
      xz(6,2)=1.d0
      xz(3,3)=1.d0
      xz(4,4)=1.d0

      yz(3,5)=1.d0
      yz(5,3)=1.d0
      yz(4,6)=1.d0
      yz(6,4)=1.d0
      yz(1,1)=1.d0
      yz(2,2)=1.d0

      xyz(1,3)=1.d0
      xyz(3,5)=1.d0
      xyz(5,1)=1.d0
      xyz(2,4)=1.d0
      xyz(4,6)=1.d0
      xyz(6,2)=1.d0

      xzy(1,5)=1.d0
      xzy(5,3)=1.d0
      xzy(3,1)=1.d0
      xzy(2,6)=1.d0
      xzy(6,4)=1.d0
      xzy(4,2)=1.d0

      ic=0
      xrold=1000000000.d0
      call movemul(rt,s,rto,xr)
! write(6,*) xr,xrold
!  do i=1,6
!       write(6,'(6(1x,1pe12.5))') (RTO(i,j),j=1,6)
!  enddo
!  ipause=mypause(0)
      if(xr.lt.xrold) then
        xrold=xr
      endif

      if(nd.ge.2) then
        call movemul(rt,xy,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=1
        endif
      endif

      if(nd.eq.3) then
        call movemul(rt,xz,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=2
        endif
        call movemul(rt,yz,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=3
        endif
        call movemul(rt,xyz,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=4
        endif
        call movemul(rt,xzy,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=5
        endif
      endif

      if(ic.eq.0) then
        call movemul(rt,s,rto,xr)
        if(idpr.gt.-101) write(6,*) " no exchanged"
      elseif(ic.eq.1) then
        call movemul(rt,xy,rto,xr)
        if(idpr.gt.-101) write(6,*) " x-y exchanged"
      elseif(ic.eq.2) then
        call movemul(rt,xz,rto,xr)
        if(idpr.gt.-101) write(6,*) " x-z exchanged"
      elseif(ic.eq.3) then
        call movemul(rt,yz,rto,xr)
        if(idpr.gt.-101) write(6,*) " y-z exchanged"
      elseif(ic.eq.4) then
        call movemul(rt,xyz,rto,xr)
        if(idpr.gt.-101) write(6,*) " x-y-z permuted"
      elseif(ic.eq.5) then
        call movemul(rt,xzy,rto,xr)
        if(idpr.gt.-101) write(6,*) " x-z-y permuted"
      endif

      do i=1,nd2
        do j=1,nd2
          rt(i,j)=rto(i,j)
        enddo
      enddo

      return
      end subroutine

      subroutine movemul(rt,xy,rto,xr)
      implicit none
      integer i,j,k,ndim,ndim2
      double precision rt,rto,xr,xy
      parameter (ndim2=6)
      parameter (ndim=3)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      dimension rt(ndim2,ndim2)
      dimension xy(ndim2,ndim2),rto(ndim2,ndim2)

      do i=1,nd2
        do j=1,nd2
          rto(i,j)=0.d0
        enddo
      enddo

      do  i=1,nd2
        do  j=1,nd2
          do  k=1,nd2
            rto(i,k)=xy(i,j)*rt(j,k)+rto(i,k)
          enddo
        enddo
      enddo

      xr=0.d0
      do i=1,nd2
        do j=1,nd2
          xr=xr+dabs(rto(i,j))
        enddo
      enddo
      do i=1,nd
        xr=xr-dabs(rto(2*i-1,2*i-1))
        xr=xr-dabs(rto(2*i-1,2*i))
        xr=xr-dabs(rto(2*i,2*i))
        xr=xr-dabs(rto(2*i,2*i-1))
      enddo
      return
      end subroutine

      subroutine initpert(st,ang,ra)
      implicit none
      integer i,ndim,ndim2,nn,nreso
      double precision ang,ra,st
!   X-RATED
!- SETS UP ALL THE COMMON BLOCKS RELEVANT TO NORMAL FORM AND THE BASIS
!- CHANGES INSIDE  MAPNORMF
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (nreso=20)
      dimension st(ndim),ang(ndim),ra(ndim)
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer mx,nres
      common /reson/mx(ndim,nreso),nres
      integer iref
      common /resfile/iref

      if(iref.gt.0) then
        write(6,*) iref
        read(iref,*) nres
        if(nres.ge.nreso) then
          write(6,*) ' NRESO IN LIELIB TOO SMALL '
          stop999
        endif
      elseif(iref.eq.0) then
        nres=0
      endif
      if(nres.ne.0) write(6,*)' warning resonances left in the map'
      if(iref.gt.0) then
        do i=1,nres
          read(iref,*) (mx(nn,i),nn=1,nd-ndc)
        enddo
      endif
      do i=nres+1,nreso
        do nn=1,ndim
          mx(nn,i)=0
        enddo
      enddo
!      frank/Etienne
      do i=1,ndim
        angle(i)=0.d0
        rad(i)=0.d0
        sta(i)=0.d0
        dsta(i)=1.d0-sta(i)
        ista(i)=0
        idsta(i)=0
      enddo
      do i=1,nd        !  Frank          -ndc
        angle(i)=ang(i)
        rad(i)=ra(i)
        sta(i)=st(i)
        dsta(i)=1.d0-sta(i)
      enddo
      do i=1,nd
        ista(i)=idint(sta(i)+.01)
        idsta(i)=idint(dsta(i)+.01)
      enddo
      return
      end subroutine

      real(C_DOUBLE) function dlie(j) bind(C, name="dlie_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      integer J(NTT)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      dlie=0.d0
      do i=1,nd
        dlie=dble(j(2*i-1)+j(2*i))+dlie
      enddo
      dlie=dlie+1.d0
      dlie=1.d0/dlie
      return
      end function

      real(C_DOUBLE) function rext(j) bind(C, name="rext_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,lie,mo,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      lie=0
      do i=1,nd-ndc
        lie=ista(i)*j(2*i)+lie
      enddo
      mo=mod(lie,4)+1
      goto(11,12,13,14),mo
 11   rext = 1.d0
      return
 12   rext = -1.d0
      return
 13   rext = -1.d0
      return
 14   rext = 1.d0
      return
      end function

      subroutine cpart(h,ch) bind(C, name="cpart_")
      use iso_c_binding
      implicit none
      integer(C_LONG) h,ch
      integer ndim,ntt
      double precision rext
      parameter (ndim=3)
      parameter (ntt=40)
      external rext
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      call dacfu(h,rext,ch)
      return
      end subroutine

      subroutine ctoi(f1,f2) bind(C, name="ctoi_")
      use iso_c_binding
      implicit none
      integer(C_LONG) f1,f2
      integer ndim2,ntt
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer b1,x(ndim2)
!
!
      call etallnom1(b1,'B1        ')
      call etallnom(x,nd2  ,'X         ')

      call cpart(f1,b1)
      call etctr(x)
      call trx(b1,f2,x)
      call dadal(x,nd2)
      call dadal1(b1)
      return
      end subroutine

      subroutine itoc(f1,f2)
      implicit none
      integer ndim2,ntt
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer f1,f2
      integer b1,x(ndim2)
!
      call etallnom1(b1,'B1        ')
      call etallnom(x,nd2  ,'X         ')

      call etrtc(x)
      call trx(f1,b1,x)
      call cpart(b1,f2)
      call dadal(x,nd2)
      call dadal1(b1)
      return
      end subroutine

      subroutine etrtc(x) bind(C, name="etrtc_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*)
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt

      integer rel(ndim2)
!
!
      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        call daadd(rel(2*i-1),rel(2*i),x(2*i-1))
        call dasub(rel(2*i-1),rel(2*i),x(2*i))
      enddo
      call dadal(rel,nd2)
      return
      end subroutine

      subroutine etctr(x) bind(C, name="etctr_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*)
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer rel(ndim2)
!
!
      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        call dalin(rel(2*i-1),.5d0,rel(2*i),.5d0,x(2*i-1))
        call dalin(rel(2*i-1),.5d0,rel(2*i),-.5d0,x(2*i))
      enddo
      call dadal(rel,nd2)
      return
      end subroutine

      subroutine etcjg(x) bind(C, name="etcjg_")
      use iso_c_binding
      implicit none
      integer(C_LONG) x(*)
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt

      integer rel(ndim2)
!
!
      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        if(ista(i).eq.1) then
          call dacop(rel(2*i-1),x(2*i))
          call dacop(rel(2*i),x(2*i-1))
        else
          call dacop(rel(2*i-1),x(2*i-1))
          call dacop(rel(2*i),x(2*i))
        endif
      enddo
      call dadal(rel,nd2)
      return
      end subroutine

      logical(C_BOOL) function eig6(fm,reval,aieval,revec,aievec)       &
     &                         bind(C, name="eig6_")
      use iso_c_binding
      implicit none
      integer(C_LONG) ndim2
      parameter (ndim2=6)
      real(C_DOUBLE) fm(ndim2,ndim2), aieval(ndim2), reval(ndim2),      &
     &               revec(ndim2,ndim2), aievec(ndim2,ndim2)

      integer jet
!**************************************************************************

!  Diagonalization routines of NERI

!ccccccccccccccccc
!
!  this routine finds the eigenvalues and eigenvectors
!  of the full matrix fm.
!  the eigenvectors are normalized so that the real and
!  imaginary part of vectors 1, 3, and 5 have +1 antisymmetric
!  product:
!      revec1 J aivec1 = 1 ; revec3 J aivec3 = 1 ;
!      revec5 J aivec5 = 1.
!  the eigenvectors 2 ,4, and 6 have the opposite normalization.
!  written by F. Neri, Feb 26 1986.
!
      integer nn
      integer ilo,ihi,mdim,info
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      double precision aa(ndim2,ndim2)
      integer i,i1
      double precision ort(ndim2),vv(ndim2,ndim2)

      eig6 = .true.
!  copy matrix to temporary storage (the matrix aa is destroyed)
      do i=1,nd2-ndc2
        do i1=1,nd2-ndc2
          aa(i1,i) = fm(i1,i)
        enddo
      enddo
      ilo = 1
      ihi = nd2-ndc2
      mdim = ndim2
      nn = nd2-ndc2
!  compute eigenvalues and eigenvectors using double
!  precision Eispack routines:
      call ety(mdim,nn,ilo,ihi,aa,ort)
      call etyt(mdim,nn,ilo,ihi,aa,ort,vv)
      call ety2(mdim,nn,ilo,ihi,aa,reval,aieval,vv,info)
      if ( info .ne. 0 ) then
        write(6,*) '  ERROR IN EIG6'
        eig6 = .false.
        return
      endif
!      call neigv(vv,pbkt)
      do i=1,nd-ndc
        do jet=1,nd2-ndc2
          revec(jet,2*i-1)=vv(jet,2*i-1)
          revec(jet,2*i)=vv(jet,2*i-1)
          aievec(jet,2*i-1)=vv(jet,2*i)
          aievec(jet,2*i)=-vv(jet,2*i)
        enddo
      enddo
      do i=1,nd2-ndc2
        if(dabs(reval(i)**2+aieval(i)**2 -1.d0).gt.1.d-10) then
          write(6,*) ' EIG6: Eigenvalues off the unit circle!'
          eig6 = .false.
        endif
      enddo
      return
      end function

      subroutine ety(nm,n,low,igh,a,ort)
      implicit none
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision a(nm,n),ort(igh)
      double precision f,g,h,scale
!
!     this subroutine is a translation of the algol procedure orthes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     orthogonal similarity transformations.
!
!     on input-
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n,
!
!        a contains the input matrix.
!
!     on output-
!
!        a contains the hessenberg matrix.  information about
!          the orthogonal transformations used in the reduction
!          is stored in the remaining triangle under the
!          hessenberg matrix,
!
!        ort contains further information about the transformations.
!          only elements low through igh are used.
!
!     fortran routine by b. s. garbow
!     modified by filippo neri.
!
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do m = kp1, la
        h = 0.0
        ort(m) = 0.0
        scale = 0.0
!     ********** scale column (algol tol then not needed) **********
        do i = m, igh
          scale = scale + dabs(a(i,m-1))
        enddo
!
        if (scale .eq. 0.0) go to 180
        mp = m + igh
!     ********** for i=igh step -1 until m do -- **********
        do ii = m, igh
          i = mp - ii
          ort(i) = a(i,m-1) / scale
          h = h + ort(i) * ort(i)
        enddo
!
        g = -dsign(dsqrt(h),ort(m))
        h = h - ort(m) * g
        ort(m) = ort(m) - g
!     ********** form (i-(u*ut)/h) * a **********
        do j = m, n
          f = 0.0
!     ********** for i=igh step -1 until m do -- **********
          do ii = m, igh
            i = mp - ii
            f = f + ort(i) * a(i,j)
          enddo
!
          f = f / h
!
          do i = m, igh
            a(i,j) = a(i,j) - f * ort(i)
          enddo
!
        enddo
!     ********** form (i-(u*ut)/h)*a*(i-(u*ut)/h) **********
        do i = 1, igh
          f = 0.0
!     ********** for j=igh step -1 until m do -- **********
          do jj = m, igh
            j = mp - jj
            f = f + ort(j) * a(i,j)
          enddo
!
          f = f / h
!
          do j = m, igh
            a(i,j) = a(i,j) - f * ort(j)
          enddo
!
        enddo
!
        ort(m) = scale * ort(m)
        a(m,m-1) = scale * g
  180   continue
      enddo
!
  200 return
!     ********** last card of ety **********
      end subroutine

      subroutine etyt(nm,n,low,igh,a,ort,z)
      implicit none
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),ort(igh),z(nm,n)
      double precision g
!
!     this subroutine is a translation of the algol procedure ortrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the orthogonal similarity
!     transformations used in the reduction of a real general
!     matrix to upper hessenberg form by  ety.
!
!     on input-
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n,
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  orthes
!          in its strict lower triangle,
!
!          ort contains further information about the trans-
!          formations used in the reduction by  ety.
!          only elements low through igh are used.
!
!     on output-
!
!        z contains the transformation matrix produced in the
!          reduction by  ety,
!
!        ort has been altered.
!
!     fortran routine by b. s. garbow.
!     modified by f. neri.
!
!
!     ********** initialize z to identity matrix **********
      do i = 1, n
!
        do j = 1, n
          z(i,j) = 0.0
        enddo
!
        z(i,i) = 1.0
      enddo
!
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
!     ********** for mp=igh-1 step -1 until low+1 do -- **********
      do mm = 1, kl
        mp = igh - mm
        if (a(mp,mp-1) .eq. 0.0) go to 140
        mp1 = mp + 1
!
        do i = mp1, igh
          ort(i) = a(i,mp-1)
        enddo
!
        do j = mp, igh
          g = 0.0
!
          do i = mp, igh
            g = g + ort(i) * z(i,j)
          enddo
!     ********** divisor below is negative of h formed in orthes.
!                double division avoids possible underflow **********
          g = (g / ort(mp)) / a(mp,mp-1)
!
          do i = mp, igh
            z(i,j) = z(i,j) + g * ort(i)
          enddo
!
        enddo
!
  140   continue
      enddo
!
  200 return
!     ********** last card of etyt **********
      end subroutine

      subroutine ety2(nm,n,low,igh,h,wr,wi,z,ierr)
      implicit none
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,                      &
     &igh,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n),z(nm,n)
      double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,machep
      logical notlas
      double precision z3r,z3i
!
!
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
!
!     on input-
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n,
!
!        h contains the upper hessenberg matrix,
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output-
!
!        h has been destroyed,
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n,
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found,
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 200 iterations.
!
!     arithmetic is double precision. complex division
!     is simulated by routin etdiv.
!
!     fortran routine by b. s. garbow.
!     modified by f. neri.
!
!
!     ********** machep is a machine dependent parameter specifying
!                the relative precision of floating point arithmetic.
!
!                **********
      machep = 1.d-17
!     machep = r1mach(4)
!
      ierr = 0
      norm = 0.0
      k = 1
!     ********** store roots isolated by balanc
!                and compute matrix norm **********
      do i = 1, n
!
        do j = k, n
          norm = norm + dabs(h(i,j))
        enddo
!
        k = i
        if (i .ge. low .and. i .le. igh) go to 50
        wr(i) = h(i,i)
        wi(i) = 0.0
   50   continue
      enddo
!
      en = igh
      t = 0.0
!     ********** search for next eigenvalues **********
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
!     ********** look for single small sub-diagonal element
!                for l=en step -1 until low do -- **********
   70 do ll = low, en
        l = en + low - ll
        if (l .eq. low) go to 100
        s = dabs(h(l-1,l-1)) + dabs(h(l,l))
        if (s .eq. 0.0) s = norm
        if (dabs(h(l,l-1)) .le. machep * s) go to 100
      enddo
!     ********** form shift **********
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (its .eq. 200) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     ********** form exceptional shift **********
      t = t + x
!
      do i = low, en
        h(i,i) = h(i,i) - x
      enddo
!
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75 * s
      y = x
      w = -0.4375 * s * s
  130 its = its + 1
!     ********** look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- **********
      do mm = l, enm2
        m = enm2 + l - mm
        zz = h(m,m)
        r = x - zz
        s = y - zz
        p = (r * s - w) / h(m+1,m) + h(m,m+1)
        q = h(m+1,m+1) - zz - r - s
        r = h(m+2,m+1)
        s = dabs(p) + dabs(q) + dabs(r)
        p = p / s
        q = q / s
        r = r / s
        if (m .eq. l) go to 150
        if (dabs(h(m,m-1)) * (dabs(q) + dabs(r)) .le. machep * dabs(p)  &
     &  * (dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))) go to 150
      enddo
!
  150 mp2 = m + 2
!
      do i = mp2, en
        h(i,i-2) = 0.0
        if (i .eq. mp2) go to 160
        h(i,i-3) = 0.0
  160   continue
      enddo
!     ********** double qr step involving rows l to en and
!                columns m to en **********
      do k = m, na
        notlas = k .ne. na
        if (k .eq. m) go to 170
        p = h(k,k-1)
        q = h(k+1,k-1)
        r = 0.0
        if (notlas) r = h(k+2,k-1)
        x = dabs(p) + dabs(q) + dabs(r)
        if (x .eq. 0.0) go to 260
        p = p / x
        q = q / x
        r = r / x
  170   s = dsign(dsqrt(p*p+q*q+r*r),p)
        if (k .eq. m) go to 180
        h(k,k-1) = -s * x
        go to 190
  180   if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190   p = p + s
        x = p / s
        y = q / s
        zz = r / s
        q = q / p
        r = r / p
!     ********** row modification **********
        do j = k, n
          p = h(k,j) + q * h(k+1,j)
          if (.not. notlas) go to 200
          p = p + r * h(k+2,j)
          h(k+2,j) = h(k+2,j) - p * zz
  200     h(k+1,j) = h(k+1,j) - p * y
          h(k,j) = h(k,j) - p * x
        enddo
!
        j = min0(en,k+3)
!     ********** column modification **********
        do i = 1, j
          p = x * h(i,k) + y * h(i,k+1)
          if (.not. notlas) go to 220
          p = p + zz * h(i,k+2)
          h(i,k+2) = h(i,k+2) - p * r
  220     h(i,k+1) = h(i,k+1) - p * q
          h(i,k) = h(i,k) - p
        enddo
!     ********** accumulate transformations **********
        do i = low, igh
          p = x * z(i,k) + y * z(i,k+1)
          if (.not. notlas) go to 240
          p = p + zz * z(i,k+2)
          z(i,k+2) = z(i,k+2) - p * r
  240     z(i,k+1) = z(i,k+1) - p * q
          z(i,k) = z(i,k) - p
        enddo
!
  260   continue
      enddo
!
      go to 70
!     ********** one root found **********
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0
      en = na
      go to 60
!     ********** two roots found **********
  280 p = (y - x) / 2.0
      q = p * p + w
      zz = dsqrt(dabs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0) go to 320
!     ********** real pair **********
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0) wr(en) = x - w / zz
      wi(na) = 0.0
      wi(en) = 0.0
      x = h(en,na)
      s = dabs(x) + dabs(zz)
      p = x / s
      q = zz / s
      r = dsqrt(p*p+q*q)
      p = p / r
      q = q / r
!     ********** row modification **********
      do j = na, n
        zz = h(na,j)
        h(na,j) = q * zz + p * h(en,j)
        h(en,j) = q * h(en,j) - p * zz
      enddo
!     ********** column modification **********
      do i = 1, en
        zz = h(i,na)
        h(i,na) = q * zz + p * h(i,en)
        h(i,en) = q * h(i,en) - p * zz
      enddo
!     ********** accumulate transformations **********
      do i = low, igh
        zz = z(i,na)
        z(i,na) = q * zz + p * z(i,en)
        z(i,en) = q * z(i,en) - p * zz
      enddo
!
      go to 330
!     ********** complex pair **********
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!     ********** all roots found.  backsubstitute to find
!                vectors of upper triangular form **********
  340 if (norm .eq. 0.0) go to 1001
!     ********** for en=n step -1 until 1 do -- **********
      do nn = 1, n
        en = n + 1 - nn
        p = wr(en)
        q = wi(en)
        na = en - 1
        if (q.lt.0) goto 710
        if (q.eq.0) goto 600
        if (q.gt.0) goto 800
!     ********** real vector **********
  600   m = en
        h(en,en) = 1.0
        if (na .eq. 0) go to 800
!     ********** for i=en-1 step -1 until 1 do -- **********
        do ii = 1, na
          i = en - ii
          w = h(i,i) - p
          r = h(i,en)
          if (m .gt. na) go to 620
!
          do j = m, na
            r = r + h(i,j) * h(j,en)
          enddo
!
  620     if (wi(i) .ge. 0.0) go to 630
          zz = w
          s = r
          go to 700
  630     m = i
          if (wi(i) .ne. 0.0) go to 640
          t = w
          if (w .eq. 0.0) t = machep * norm
          h(i,en) = -r / t
          go to 700
!     ********** solve real equations **********
  640     x = h(i,i+1)
          y = h(i+1,i)
          q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
          t = (x * s - zz * r) / q
          h(i,en) = t
          if (dabs(x) .le. dabs(zz)) go to 650
          h(i+1,en) = (-r - w * t) / x
          go to 700
  650     h(i+1,en) = (-s - y * t) / zz
  700     continue
        enddo
!     ********** end real vector **********
        go to 800
!     ********** complex vector **********
  710   m = na
!     ********** last vector component chosen imaginary so that
!                eigenvector matrix is triangular **********
        if (dabs(h(en,na)) .le. dabs(h(na,en))) go to 720
        h(na,na) = q / h(en,na)
        h(na,en) = -(h(en,en) - p) / h(en,na)
        go to 730
! 720    z3 = cmplx(0.0,-h(na,en)) / cmplx(h(na,na)-p,q)
!        h(na,na) = real(z3)
!        h(na,en) = aimag(z3)
  720   call etdiv(z3r,z3i,0.d0,-h(na,en),h(na,na)-p,q)
        h(na,na) = z3r
        h(na,en) = z3i
  730   h(en,na) = 0.0
        h(en,en) = 1.0
        enm2 = na - 1
        if (enm2 .eq. 0) go to 800
!     ********** for i=en-2 step -1 until 1 do -- **********
        do ii = 1, enm2
          i = na - ii
          w = h(i,i) - p
          ra = 0.0
          sa = h(i,en)
!
          do j = m, na
            ra = ra + h(i,j) * h(j,na)
            sa = sa + h(i,j) * h(j,en)
          enddo
!
          if (wi(i) .ge. 0.0) go to 770
          zz = w
          r = ra
          s = sa
          go to 790
  770     m = i
          if (wi(i) .ne. 0.0) go to 780
!           z3 = cmplx(-ra,-sa) / cmplx(w,q)
!           h(i,na) = real(z3)
!           h(i,en) = aimag(z3)
          call etdiv(z3r,z3i,-ra,-sa,w,q)
          h(i,na) = z3r
          h(i,en) = z3i
          go to 790
!     ********** solve complex equations **********
  780     x = h(i,i+1)
          y = h(i+1,i)
          vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
          vi = (wr(i) - p) * 2.0 * q
          if (vr .eq. 0.0 .and. vi .eq. 0.0) vr = machep * norm         &
     &    * (dabs(w) + dabs(q) + dabs(x) + dabs(y) + dabs(zz))
!           z3 = cmplx(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra) / cmplx(vr,vi)
!           h(i,na) = real(z3)
!           h(i,en) = aimag(z3)
          call etdiv(z3r,z3i,x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi)
          h(i,na) = z3r
          h(i,en) = z3i
          if (dabs(x) .le. dabs(zz) + dabs(q)) go to 785
          h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
          h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
          go to 790
! 785       z3 = cmplx(-r-y*h(i,na),-s-y*h(i,en)) / cmplx(zz,q)
!           h(i+1,na) = real(z3)
!           h(i+1,en) = aimag(z3)
  785     call etdiv(z3r,z3i,-r-y*h(i,na),-s-y*h(i,en),zz,q)
          h(i+1,na) = z3r
          h(i+1,en) = z3i
  790     continue
        enddo
!     ********** end complex vector **********
  800   continue
      enddo
!     ********** end back substitution.
!                vectors of isolated roots **********
      do i = 1, n
        if (i .ge. low .and. i .le. igh) go to 840
!
        do j = i, n
          z(i,j) = h(i,j)
        enddo
!
  840   continue
      enddo
!     ********** multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- **********
      do jj = low, n
        j = n + low - jj
        m = min0(j,igh)
!
        do i = low, igh
          zz = 0.0
!
          do k = low, m
            zz = zz + z(i,k) * h(k,j)
          enddo
!
          z(i,j) = zz
        enddo
      enddo
!
      go to 1001
!     ********** set error -- no convergence to an
!                eigenvalue after 200 iterations **********
 1000 ierr = en
 1001 return
!     ********** last card of ety2 **********
      end subroutine

      subroutine etdiv(a,b,c,d,e,f)
      implicit none
!   computes the complex division
!     a + ib = (c + id)/(e + if)
!  very slow, but tries to be as accurate as
!  possible by changing the order of the
!  operations, so to avoid under(over)flow
!  problems.
!  Written by F. Neri Feb. 12 1986
!
      double precision a,b,c,d,e,f
      double precision s,t
      double precision cc,dd,ee,ff
      double precision temp
      integer flip
      flip = 0
      cc = c
      dd = d
      ee = e
      ff = f
      if( dabs(f).ge.dabs(e) ) then
        ee = f
        ff = e
        cc = d
        dd = c
        flip = 1
      endif
      s = 1.d0/ee
      t = 1.d0/(ee+ ff*(ff*s))
      if ( dabs(ff) .ge. dabs(s) ) then
        temp = ff
        ff = s
        s = temp
      endif
      if( dabs(dd) .ge. dabs(s) ) then
        a = t*(cc + s*(dd*ff))
      else if ( dabs(dd) .ge. dabs(ff) ) then
        a = t*(cc + dd*(s*ff))
      else
        a = t*(cc + ff*(s*dd))
      endif
      if ( dabs(cc) .ge. dabs(s)) then
        b = t*(dd - s*(cc*ff))
      else if ( dabs(cc) .ge. dabs(ff)) then
        b = t*(dd - cc*(s*ff))
      else
        b = t*(dd - ff*(s*cc))
      endif
      if (flip.ne.0 ) then
        b = -b
      endif
      return
      end subroutine

      subroutine sympl3(m)
!**********************************************************
!
!    SYMPL3
!
!
!   On return ,the matrix m(*,*), supposed to be almost
!   symplectic on entry is made exactly symplectic by
!   using a non iterative, constructive method.
!
!**********************************************************
!
!  Written by F. Neri  Feb 7 1986
!
      implicit none
      integer n
      parameter ( n = 3 )
      integer kp,kq,lp,lq,jp,jq,i
      double precision m(2*n,2*n)
      double precision qq,pq,qp,pp
!
      do kp=2,2*n,2
        kq = kp-1
        do lp=2,kp-2,2
          lq = lp-1
          qq = 0.d0
          pq = 0.d0
          qp = 0.d0
          pp = 0.d0
          do jp=2,2*n,2
            jq = jp-1
            qq = qq + m(lq,jq)*m(kq,jp) - m(lq,jp)*m(kq,jq)
            pq = pq + m(lp,jq)*m(kq,jp) - m(lp,jp)*m(kq,jq)
            qp = qp + m(lq,jq)*m(kp,jp) - m(lq,jp)*m(kp,jq)
            pp = pp + m(lp,jq)*m(kp,jp) - m(lp,jp)*m(kp,jq)
          enddo
!         write(6,*) qq,pq,qp,pp
          do i=1,2*n
            m(kq,i) = m(kq,i) - qq*m(lp,i) + pq*m(lq,i)
            m(kp,i) = m(kp,i) - qp*m(lp,i) + pp*m(lq,i)
          enddo
        enddo
        qp = 0.d0
        do jp=2,2*n,2
          jq = jp-1
          qp = qp + m(kq,jq)*m(kp,jp) - m(kq,jp)*m(kp,jq)
        enddo
!       write(6,*) qp
        do i=1,2*n
          m(kp,i) = m(kp,i)/qp
        enddo
!
!  Maybe the following is a better idea ( uses sqrt and is slower )
!       sign = 1.d0
!       if ( qp.lt.0.d0 ) sign = -1.d0
!  OR, BETTER:
!       if ( qp.le.0.d0 ) then complain
!       qp = dabs(qp)
!       qp = dsqrt(qp)
!       do 600 i=1,2*n
!         m(kq,i) = m(kq,i)/qp
!         m(kp,i) = sign*m(kp,i)/qp
! 600   continue
      enddo
      return
      end subroutine

      logical*1 function averaged(f,a,flag,fave)
      implicit none
      integer isi,ndim,ndim2,nord,ntt
      double precision avepol
!      TAKES THE AVERAGE OF A FUNCTION F
!  FLAG TRUE A=ONE TURN MAP
!       FALSE A=A_SCRIPT
!
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer idpr
      common /printing/ idpr
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer f,fave,a(*)
      integer cosi,sine
      logical flag
      external avepol

      integer a1(ndim2),a2(ndim2),xy(ndim2),hf(ndim2),ftf(ndim2)

      logical*1 mapnormf

      if(.not.flag) then
        call etall1(cosi)
        call etall1(sine)
        call trx(f,f,a)
        call ctor(f,cosi,sine)
        call dacfu(cosi,avepol,fave)
        call dadal1(cosi)
        call dadal1(sine)
      else

        call etall1(cosi)
        call etall1(sine)
        call etall(ftf,nd2)
        call etall(hf,nd2)
        call etall(a2,nd2)
        call etall(a1,nd2)
        call etall(xy,nd2)

        isi=0
        nord=1
        averaged = mapnormf(a,ftf,a2,a1,xy,hf,nord,isi)
        nord=no
        call etcct(a1,a2,xy)
        call facflod(ftf,xy,a1,2,nord,1.d0,-1)
        call trx(f,f,a1)
        call ctor(f,cosi,sine)
        call dacfu(cosi,avepol,fave)

        call dadal1(cosi)
        call dadal1(sine)
        call dadal(ftf,nd2)
        call dadal(hf,nd2)
        call dadal(a2,nd2)
        call dadal(a1,nd2)
        call dadal(xy,nd2)

      endif

      return
      end function

      real(C_DOUBLE) function avepol(j) bind(C, name="avepol_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      integer J(NTT)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      avepol=1.d0
      do i=1,(nd-ndc)
        if(j(2*i).ne.j(2*i-1)) then
          avepol=0.d0
          return
        endif
      enddo

      return
      end function

      logical function couplean(map1,tune,map2,oneturn)
      implicit none
      integer i,ndim,ndim2,no1,nord,ntt
      double precision crazy,tpi
!  map1 ascript a1 not there
!  tune 2 or 3 tunes

!   map2 ascript with a couple parameter in nv
!  oneturn map created with tunes and map2

      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer map1(*),oneturn(*),map2(*),ftf,hf
      integer xy(ndim2),m1(ndim2),m2(ndim2),a2(ndim2),a1(ndim2)
      integer cs,h

      double precision killnonl,planar,psq(ndim),radsq(ndim)
      double precision tune(ndim)
      external killnonl,planar

      logical*1 mapnorm

      call etall1(ftf)
      call etall1(hf)
      call etall(a1,nd2)
      call etall(a2,nd2)
      call etall(m1,nd2)
      call etall(m2,nd2)
      call etall(xy,nd2)
      call etall1(cs)
      call etall1(h)

!     map1 is an a-script, the last nv entry should be empty
!  this a-script should around the fixed point to all orders
!     one order is lost because I use PB-field

      tpi=datan(1.d0)*8.d0
      do i=1,nd2
        call dacfu(map1(i),killnonl,m1(i))
      enddo

      call etini(xy)
      call daclr(cs)

      do i=1,nd-ndc
        call dasqr(xy(2*i-1),a2(2*i-1))
        call dasqr(xy(2*i),a2(2*i))
        call daadd(a2(2*i-1),a2(2*i),ftf)
        crazy=-tune(i)*tpi/2.d0
        call dacmu(ftf,crazy,ftf)
        call daadd(ftf,cs,cs)
      enddo

      call etinv(m1,m2)
      call trx(cs,h,m2)

      call dacfu(h,planar,cs)
      call dasub(h,cs,h)
      call davar(a2(1),0.d0,nv)

      call damul(a2(1),h,h)
      call daadd(cs,h,h)
      call expnd2(h,xy,xy,1.d-9,1000)

      call dacopd(xy,oneturn)

      nord=1
      couplean = mapnorm(xy,ftf,a2,a1,m2,hf,nord)

      call gettura(psq,radsq)
      write(6,*) (psq(i),i=1,nd)

      call etini(xy)
      no1=no
      call fexpo(ftf,xy,xy,3,no1,1.d0,-1)
      call etcct(a2,xy,map2)
      call etcct(a1,map2,map2)

      call dadal1(ftf)
      call dadal1(hf)
      call dadal(a1,nd2)
      call dadal(a2,nd2)
      call dadal(m1,nd2)
      call dadal(m2,nd2)
      call dadal(xy,nd2)
      call dadal1(cs)
      call dadal1(h)

      return
      end function

      real(C_DOUBLE) function planar(j) bind(C, name="planar_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      integer J(NTT)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
      planar=0.d0
      do i=1,(nd-ndc)
        if(j(2*i).eq.j(2*i-1)) then
          planar=1.d0
          return
        endif
        if(j(2*i).eq.2) then
          planar=1.d0
          return
        endif
        if(j(2*i-1).eq.2) then
          planar=1.d0
          return
        endif
      enddo

      return
      end function

      real(C_DOUBLE) function killnonl(j) bind(C, name="killnonl_")
      use iso_c_binding
      implicit none
      integer(C_LONG) j(*)

      integer i,ic,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      integer J(NTT)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt

      killnonl=1.d0

      ic=0
      do i=1,nd2-ndc2
        ic=ic+j(i)
      enddo
      if(ic.gt.1) killnonl=0.d0
      if(j(nv).ne.0) killnonl=0.d0

      return
      end function

      subroutine fexpo1(h,x,w,nrmin,nrmax,sca,ifac)
      implicit none
      integer ifac,ndim,ndim2,nrma,nrmax,nrmi,nrmin,ntt
      double precision sca
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer x,w,h

      integer v(ndim2)

      nrmi=nrmin-1
      nrma=nrmax-1
      call etall(v,nd2)
      call difd(h,v,-1.d0)
      call facflo(v,x,w,nrmi,nrma,sca,ifac)
      call dadal(v,nd2)

      return
      end subroutine

      subroutine etcctpar(x,ix,xj,z)
      implicit none
      integer i,ie,ix,ndim,ndim2,ntt
      double precision xj
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension xj(*)
      dimension ie(ntt)
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
      integer  x(*),z(*)

      call etallnom(ie,nv,'IE        ')
      do i=1,nd2
        call davar(ie(i),0.d0,i)
      enddo
      do  i=nd2+1,nv
        call dacon(ie(i),xj(i-nd2))
      enddo

      call dacct(x,ix,ie,nv,z,ix)

      call dadal(ie,nv)
      return
      end subroutine

      end module lielib
