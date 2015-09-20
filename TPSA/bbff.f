      subroutine bbff(sepx,sepy,sigxx,sigyy,bbfx,bbfy,bbgx,bbgy)
C-----------------------------------------------------------------------
C      
C   Hirata's 6d beam-beam from BBC
C   DA VERSION for SIXTRACK courtesy Peter Leunissen
C   January 1999
C      
C-----------------------------------------------------------------------
C**********************************************************************
C BBFF gives transverse (f_x and f_y) and longitudinal(g_x and g_y)
C beam-beam kicks except for the kinematical term (nr_e/\gamma)
C SIGXX is \Sigma
C**********************************************************************
      implicit none
      integer idaa
      double precision dare,hundred,sqrpi2
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele = 1000,nblo = 400,nper = 16, nelb = 140,nblz =
     &15000,nzfz = 300000,mmul = 11)
      parameter(nran = 280000,ncom = 100,mran = 500,mpa = 6,nrco = 5,
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 15)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 160)
      double precision c180e0,c1e1,c1e12,c1e13,c1e15,c1e16,c1e2,c1e3,
     &c1e4,c1e6,c1m1,c1m10,c1m12,c1m13,c1m15,c1m18,c1m2,c1m21,c1m24,
     &c1m3,c1m36,c1m38,c1m6,c1m9,c2e3,c4e3,crade,four,half,one,pieni,
     &pmae,pmap,three,two,zero
      parameter(pieni = 1d-38)
      parameter(zero = 0.0d0,half = 0.5d0,one = 1.0d0)
      parameter(two = 2.0d0,three = 3.0d0,four = 4.0d0)
      parameter(c1e1 = 1.0d1,c1e2 = 1.0d2,c1m2 = 1.0d-2)
      parameter(c1e3 = 1.0d3,c2e3 = 2.0d3,c4e3 = 4.0d3,c1e4 = 1.0d4)
      parameter(c1e12 = 1.0d12,c1e13 = 1.0d13,c1e15 = 1.0d15,c1e16 =
     &1.0d16)
      parameter(c180e0 = 180.0d0,c1e6 = 1.0d6)
      parameter(c1m1 = 1.0d-1,c1m3 = 1.0d-3,c1m6 = 1.0d-6,c1m9 = 1.0d-9)
      parameter(c1m10 = 1.0d-10,c1m12 = 1.0d-12,c1m13 = 1.0d-13)
      parameter(c1m15 = 1.0d-15)
      parameter(c1m18 = 1.0d-18,c1m21 = 1.0d-21,c1m24 = 1.0d-24)
      parameter(c1m36 = 1.0d-36,c1m38 = 1.0d-38)
      parameter(pmap = 938.271998d0,pmae = .510998902d0)
      parameter(crade = 2.817940285d-15)
      parameter(sqrpi2 = 3.544907701811032d0,hundred = 100d0)
      integer icorr,idial,idptr,imod1,imod2,inorm,ipar,namp,ncor,nctype,
     &ndimf,nmom,nmom1,nmom2,nord,nord1,nordf,nsix,nvar,nvar2,nvarf
      double precision dpmax,preda,weig1,weig2
      character*16 coel
      common/dial/preda,idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
      common/norf/nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
      common/tcorr/icorr,nctype,namp,nmom,nmom1,nmom2,weig1,weig2,dpmax,
     &coel(10)
      save
C-----------------------------------------------------------------------
*FOX  B D ;
*FOX  D V DA EXT SEPX NORD NVAR ; D V DA EXT SEPY NORD NVAR ;
*FOX  D V DA EXT SIGXX NORD NVAR ; D V DA EXT SIGYY NORD NVAR ;
*FOX  D V DA EXT BBFX NORD NVAR ; D V DA EXT BBFY NORD NVAR ;
*FOX  D V DA EXT BBGX NORD NVAR ; D V DA EXT BBGY NORD NVAR ;
*FOX  D V DA INT SIGXY NORD NVAR ; D V DA INT EXPFAC NORD NVAR ;
*FOX  D V DA INT X NORD NVAR ; D V DA INT FAC NORD NVAR ;
*FOX  D V DA INT FAC2 NORD NVAR ; D V DA INT CONST NORD NVAR ;
*FOX  D V DA INT ARG1X NORD NVAR ; D V DA INT ARG1Y NORD NVAR ;
*FOX  D V DA INT ARG2X NORD NVAR ; D V DA INT ARG2Y NORD NVAR ;
*FOX  D V DA INT WX1 NORD NVAR ; D V DA INT WY1 NORD NVAR ;
*FOX  D V DA INT WX2 NORD NVAR ; D V DA INT WY2 NORD NVAR ;
*FOX  D V DA INT COMFAC NORD NVAR ; D V DA INT COMFAC2 NORD NVAR ;
*FOX  D V RE INT ZERO ; D V RE INT HALF ; D V RE INT ONE ;
*FOX  D V RE INT TWO ; D V RE INT FOUR ; D V RE INT HUNDRED ;
*FOX  D V RE INT SQRPI2 ;
*FOX  D F RE DARE 1 ;
*FOX  E D ;
*FOX{
      INTEGER SEPX    
      INTEGER SEPY    
      INTEGER SIGXX   
      INTEGER SIGYY   
      INTEGER BBFX    
      INTEGER BBFY    
      INTEGER BBGX    
      INTEGER BBGY    
      INTEGER SIGXY   
      INTEGER EXPFAC  
      INTEGER X       
      INTEGER FAC     
      INTEGER FAC2    
      INTEGER CONST   
      INTEGER ARG1X   
      INTEGER ARG1Y   
      INTEGER ARG2X   
      INTEGER ARG2Y   
      INTEGER WX1     
      INTEGER WY1     
      INTEGER WX2     
      INTEGER WY2     
      INTEGER COMFAC  
      INTEGER COMFAC2 
      INTEGER ISCRDA, ISCRRI,IDAO
      DOUBLE PRECISION RSCRRI
      COMMON/DASCR/ISCRDA(100),RSCRRI(100),ISCRRI(100),IDAO
      if(1.eq.1) then                                                 
         CALL DAKEY('FOX V2.1')
         CALL DAALL(SIGXY   ,1,'SIGXY     ',NORD,NVAR)
         CALL DAALL(EXPFAC  ,1,'EXPFAC    ',NORD,NVAR)
         CALL DAALL(X       ,1,'X         ',NORD,NVAR)
         CALL DAALL(FAC     ,1,'FAC       ',NORD,NVAR)
         CALL DAALL(FAC2    ,1,'FAC2      ',NORD,NVAR)
         CALL DAALL(CONST   ,1,'CONST     ',NORD,NVAR)
         CALL DAALL(ARG1X   ,1,'ARG1X     ',NORD,NVAR)
         CALL DAALL(ARG1Y   ,1,'ARG1Y     ',NORD,NVAR)
         CALL DAALL(ARG2X   ,1,'ARG2X     ',NORD,NVAR)
         CALL DAALL(ARG2Y   ,1,'ARG2Y     ',NORD,NVAR)
         CALL DAALL(WX1     ,1,'WX1       ',NORD,NVAR)
         CALL DAALL(WY1     ,1,'WY1       ',NORD,NVAR)
         CALL DAALL(WX2     ,1,'WX2       ',NORD,NVAR)
         CALL DAALL(WY2     ,1,'WY2       ',NORD,NVAR)
         CALL DAALL(COMFAC  ,1,'COMFAC    ',NORD,NVAR)
         CALL DAALL(COMFAC2 ,1,'COMFAC2   ',NORD,NVAR)
      ENDIF
      IDAA = IDAO
*FOX}
C-----------------------------------------------------------------------
      if(dare(sigxx).eq.dare(sigyy)) then
*FOX    X=SEPX*SEPX+SEPY*SEPY ;                                         *FOX
      CALL DAMUL(SEPX       ,SEPX       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(SEPY       ,SEPY       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),X          )                          
        if(abs(dare(sigxx)+dare(sigyy)).gt.pieni) then
*FOX      CONST=X/(SIGXX+SIGYY) ;                                       *FOX
      CALL DAADD(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DADIV(X          ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DACOP(ISCRDA(  2+IDAA),CONST      )                          
        else
*FOX      CONST=ZERO ;                                                  *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(CONST      ,RSCRRI(100))                               
        endif
*FOX    EXPFAC=EXP(-CONST) ;                                            *FOX
      CALL DACMU(CONST      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAFUN('EXP ',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DACOP(ISCRDA(  2+IDAA),EXPFAC     )                          
        if(abs(dare(x)).gt.pieni) then
*FOX      BBFX=TWO*SEPX*(ONE-EXPFAC)/X ;                                *FOX
      CALL DASUC(EXPFAC     ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(SEPX       ,ONE*TWO        ,ISCRDA(  2+IDAA))          
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DADIV(ISCRDA(  3+IDAA),X          ,ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),BBFX       )                          
*FOX      BBFY=TWO*SEPY*(ONE-EXPFAC)/X ;                                *FOX
      CALL DASUC(EXPFAC     ,ONE*ONE        ,ISCRDA(  1+IDAA))          
      CALL DACMU(SEPY       ,ONE*TWO        ,ISCRDA(  2+IDAA))          
      CALL DAMUL(ISCRDA(  2+IDAA),ISCRDA(  1+IDAA),ISCRDA(  3+IDAA))    
      CALL DADIV(ISCRDA(  3+IDAA),X          ,ISCRDA(  4+IDAA))         
      CALL DACOP(ISCRDA(  4+IDAA),BBFY       )                          
*FOX      COMFAC=-SEPX*BBFX+SEPY*BBFY ;                                 *FOX
      CALL DACMU(SEPX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DAMUL(ISCRDA(  1+IDAA),BBFX       ,ISCRDA(  2+IDAA))         
      CALL DAMUL(SEPY       ,BBFY       ,ISCRDA(  3+IDAA))              
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))    
      CALL DACOP(ISCRDA(  4+IDAA),COMFAC     )                          
          if(dare(sigxx).lt.zero) then
*FOX        SIGXX=-SIGXX ;                                              *FOX
      CALL DACMU(SIGXX      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),SIGXX      )                          
          endif
          if(dare(sigyy).lt.zero) then
*FOX        SIGYY=-SIGYY ;                                              *FOX
      CALL DACMU(SIGYY      ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),SIGYY      )                          
          endif
*FOX      COMFAC2=(SIGXX+SIGYY)*(SIGXX+SIGYY) ;                         *FOX
      CALL DAADD(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DAADD(SIGXX      ,SIGYY      ,ISCRDA(  2+IDAA))              
      CALL DAMUL(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),COMFAC2    )                          
*FOX      BBGX=(COMFAC+FOUR*SEPX*SEPX*CONST/X*EXPFAC)/(TWO*X) ;         *FOX
      CALL DACMU(SEPX       ,ONE*FOUR       ,ISCRDA(  1+IDAA))          
      CALL DAMUL(ISCRDA(  1+IDAA),SEPX       ,ISCRDA(  2+IDAA))         
      CALL DAMUL(ISCRDA(  2+IDAA),CONST      ,ISCRDA(  3+IDAA))         
      CALL DADIV(ISCRDA(  3+IDAA),X          ,ISCRDA(  4+IDAA))         
      CALL DAMUL(ISCRDA(  4+IDAA),EXPFAC     ,ISCRDA(  5+IDAA))         
      CALL DACMU(X          ,ONE*TWO        ,ISCRDA(  6+IDAA))          
      CALL DAADD(COMFAC     ,ISCRDA(  5+IDAA),ISCRDA(  7+IDAA))         
      CALL DADIV(ISCRDA(  7+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DACOP(ISCRDA(  8+IDAA),BBGX       )                          
*FOX      BBGY=(-COMFAC+FOUR*SEPY*SEPY*CONST/X*EXPFAC)/(TWO*X) ;        *FOX
      CALL DACMU(COMFAC     ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(SEPY       ,ONE*FOUR       ,ISCRDA(  2+IDAA))          
      CALL DAMUL(ISCRDA(  2+IDAA),SEPY       ,ISCRDA(  3+IDAA))         
      CALL DAMUL(ISCRDA(  3+IDAA),CONST      ,ISCRDA(  4+IDAA))         
      CALL DADIV(ISCRDA(  4+IDAA),X          ,ISCRDA(  5+IDAA))         
      CALL DAMUL(ISCRDA(  5+IDAA),EXPFAC     ,ISCRDA(  6+IDAA))         
      CALL DACMU(X          ,ONE*TWO        ,ISCRDA(  7+IDAA))          
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  6+IDAA),ISCRDA(  8+IDAA))    
      CALL DADIV(ISCRDA(  8+IDAA),ISCRDA(  7+IDAA),ISCRDA(  9+IDAA))    
      CALL DACOP(ISCRDA(  9+IDAA),BBGY       )                          
        else
*FOX      BBFX=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBFX       ,RSCRRI(100))                               
*FOX      BBFY=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBFY       ,RSCRRI(100))                               
*FOX      BBGX=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBGX       ,RSCRRI(100))                               
*FOX      BBGY=ZERO ;                                                   *FOX
      RSCRRI(100) = ZERO                                                
      CALL DACON(BBGY       ,RSCRRI(100))                               
        endif
      else
*FOX    X=SEPX*SEPX/SIGXX+SEPY*SEPY/SIGYY ;                             *FOX
      CALL DAMUL(SEPX       ,SEPX       ,ISCRDA(  1+IDAA))              
      CALL DADIV(ISCRDA(  1+IDAA),SIGXX      ,ISCRDA(  2+IDAA))         
      CALL DAMUL(SEPY       ,SEPY       ,ISCRDA(  3+IDAA))              
      CALL DADIV(ISCRDA(  3+IDAA),SIGYY      ,ISCRDA(  4+IDAA))         
      CALL DAADD(ISCRDA(  2+IDAA),ISCRDA(  4+IDAA),ISCRDA(  5+IDAA))    
      CALL DACOP(ISCRDA(  5+IDAA),X          )                          
*FOX    FAC2=TWO*(SIGXX-SIGYY) ;                                        *FOX
      CALL DASUB(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*TWO        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FAC2       )                          
        if(dare(sigxx).lt.dare(sigyy)) then
*FOX      FAC2=TWO*(SIGYY-SIGXX) ;                                      *FOX
      CALL DASUB(SIGYY      ,SIGXX      ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*TWO        ,ISCRDA(  2+IDAA))     
      CALL DACOP(ISCRDA(  2+IDAA),FAC2       )                          
        endif
*FOX    FAC=SQRT(FAC2) ;                                                *FOX
      CALL DAFUN('SQRT',FAC2       ,ISCRDA(  1+IDAA))                   
      CALL DACOP(ISCRDA(  1+IDAA),FAC        )                          
*FOX    CONST=SQRPI2/FAC ;                                              *FOX
      CALL DADIC(FAC        ,ONE*SQRPI2     ,ISCRDA(  1+IDAA))          
      CALL DACOP(ISCRDA(  1+IDAA),CONST      )                          
*FOX    SIGXY=SQRT(SIGXX/SIGYY) ;                                       *FOX
      CALL DADIV(SIGXX      ,SIGYY      ,ISCRDA(  1+IDAA))              
      CALL DAFUN('SQRT',ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))              
      CALL DACOP(ISCRDA(  2+IDAA),SIGXY      )                          
*FOX    ARG1X=(SEPX/FAC) ;                                              *FOX
      CALL DADIV(SEPX       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG1X      )                          
        if(dare(sepx).lt.zero) then
*FOX      ARG1X=-(SEPX/FAC) ;                                           *FOX
      CALL DADIV(SEPX       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DACOP(ISCRDA(  2+IDAA),ARG1X      )                          
        endif
*FOX    ARG1Y=(SEPY/FAC) ;                                              *FOX
      CALL DADIV(SEPY       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG1Y      )                          
        if(dare(sepy).lt.zero) then
*FOX      ARG1Y=-(SEPY/FAC) ;                                           *FOX
      CALL DADIV(SEPY       ,FAC        ,ISCRDA(  1+IDAA))              
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DACOP(ISCRDA(  2+IDAA),ARG1Y      )                          
        endif
        call errff(arg1x,arg1y,wy1,wx1)
        if(dare(x).lt.hundred) then
*FOX      EXPFAC=EXP(-X*HALF) ;                                         *FOX
      CALL DACMU(X          ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACMU(ISCRDA(  1+IDAA),ONE*HALF       ,ISCRDA(  2+IDAA))     
      CALL DAFUN('EXP ',ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))              
      CALL DACOP(ISCRDA(  3+IDAA),EXPFAC     )                          
*FOX      ARG2X=ARG1X/SIGXY ;                                           *FOX
      CALL DADIV(ARG1X      ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG2X      )                          
*FOX      ARG2Y=ARG1Y*SIGXY ;                                           *FOX
      CALL DAMUL(ARG1Y      ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),ARG2Y      )                          
          call errff(arg2x,arg2y,wy2,wx2)
*FOX      BBFX=CONST*(WX1-EXPFAC*WX2) ;                                 *FOX
      CALL DAMUL(EXPFAC     ,WX2        ,ISCRDA(  1+IDAA))              
      CALL DASUB(WX1        ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DAMUL(CONST      ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),BBFX       )                          
*FOX      BBFY=CONST*(WY1-EXPFAC*WY2) ;                                 *FOX
      CALL DAMUL(EXPFAC     ,WY2        ,ISCRDA(  1+IDAA))              
      CALL DASUB(WY1        ,ISCRDA(  1+IDAA),ISCRDA(  2+IDAA))         
      CALL DAMUL(CONST      ,ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),BBFY       )                          
          if(dare(sepx).lt.zero) then
*FOX        BBFX=-BBFX ;                                                *FOX
      CALL DACMU(BBFX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFX       )                          
          endif
          if(dare(sepy).lt.zero) then
*FOX        BBFY=-BBFY ;                                                *FOX
      CALL DACMU(BBFY       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFY       )                          
          endif
*FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;                                  *FOX
      CALL DAMUL(SEPX       ,BBFX       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(SEPY       ,BBFY       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),COMFAC     )                          
*FOX      BBGX=-(COMFAC+TWO*(EXPFAC/SIGXY-ONE))/FAC2 ;                  *FOX
      CALL DADIV(EXPFAC     ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACSU(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*TWO        ,ISCRDA(  3+IDAA))     
      CALL DAADD(COMFAC     ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DACMU(ISCRDA(  4+IDAA),ONE*(-ONE       ),ISCRDA(  5+IDAA))   
      CALL DADIV(ISCRDA(  5+IDAA),FAC2       ,ISCRDA(  6+IDAA))         
      CALL DACOP(ISCRDA(  6+IDAA),BBGX       )                          
*FOX      BBGY= (COMFAC+TWO*(EXPFAC*SIGXY-ONE))/FAC2 ;                  *FOX
      CALL DAMUL(EXPFAC     ,SIGXY      ,ISCRDA(  1+IDAA))              
      CALL DACSU(ISCRDA(  1+IDAA),ONE*ONE        ,ISCRDA(  2+IDAA))     
      CALL DACMU(ISCRDA(  2+IDAA),ONE*TWO        ,ISCRDA(  3+IDAA))     
      CALL DAADD(COMFAC     ,ISCRDA(  3+IDAA),ISCRDA(  4+IDAA))         
      CALL DADIV(ISCRDA(  4+IDAA),FAC2       ,ISCRDA(  5+IDAA))         
      CALL DACOP(ISCRDA(  5+IDAA),BBGY       )                          
        else
*FOX      BBFX=CONST*WX1 ;                                              *FOX
      CALL DAMUL(CONST      ,WX1        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),BBFX       )                          
*FOX      BBFY=CONST*WY1 ;                                              *FOX
      CALL DAMUL(CONST      ,WY1        ,ISCRDA(  1+IDAA))              
      CALL DACOP(ISCRDA(  1+IDAA),BBFY       )                          
          if(dare(sepx).lt.zero) then
*FOX        BBFX=-BBFX ;                                                *FOX
      CALL DACMU(BBFX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFX       )                          
          endif
          if(dare(sepy).lt.zero) then
*FOX        BBFY=-BBFY ;                                                *FOX
      CALL DACMU(BBFY       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBFY       )                          
          endif
*FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;                                  *FOX
      CALL DAMUL(SEPX       ,BBFX       ,ISCRDA(  1+IDAA))              
      CALL DAMUL(SEPY       ,BBFY       ,ISCRDA(  2+IDAA))              
      CALL DAADD(ISCRDA(  1+IDAA),ISCRDA(  2+IDAA),ISCRDA(  3+IDAA))    
      CALL DACOP(ISCRDA(  3+IDAA),COMFAC     )                          
*FOX      BBGX=-(COMFAC-TWO)/FAC2 ;                                     *FOX
      CALL DACSU(COMFAC     ,ONE*TWO        ,ISCRDA(  1+IDAA))          
      CALL DACMU(ISCRDA(  1+IDAA),ONE*(-ONE       ),ISCRDA(  2+IDAA))   
      CALL DADIV(ISCRDA(  2+IDAA),FAC2       ,ISCRDA(  3+IDAA))         
      CALL DACOP(ISCRDA(  3+IDAA),BBGX       )                          
*FOX      BBGY= -BBGX ;                                                 *FOX
      CALL DACMU(BBGX       ,ONE*(-ONE       ),ISCRDA(  1+IDAA))        
      CALL DACOP(ISCRDA(  1+IDAA),BBGY       )                          
        endif
      endif
        CALL DADAL(COMFAC2 ,1)                                                  
        CALL DADAL(COMFAC  ,1)                                                  
        CALL DADAL(WY2     ,1)                                                  
        CALL DADAL(WX2     ,1)                                                  
        CALL DADAL(WY1     ,1)                                                  
        CALL DADAL(WX1     ,1)                                                  
        CALL DADAL(ARG2Y   ,1)                                                  
        CALL DADAL(ARG2X   ,1)                                                  
        CALL DADAL(ARG1Y   ,1)                                                  
        CALL DADAL(ARG1X   ,1)                                                  
        CALL DADAL(CONST   ,1)                                                  
        CALL DADAL(FAC2    ,1)                                                  
        CALL DADAL(FAC     ,1)                                                  
        CALL DADAL(X       ,1)                                                  
        CALL DADAL(EXPFAC  ,1)                                                  
        CALL DADAL(SIGXY   ,1)                                                  
C     DADAL AUTOMATIC INCLUSION
      return
      end
