      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     TEST1 linear EOS+corr initial boundary 28 Mei 2020
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN, IK   
      DIMENSION YA(10), EK(4,10), Y(10)
      DOUBLE PRECISION LAMBDA, LMDCC

c---------------------------------------------------------------------
      CL=1.0D0
      YDY=-1.15
C      XC=5.0D-5   
      XC=1.0D-3
C---------------------------------------------------------------------      


c      OPEN (unit=2,STATUS='unknown',FILE='CatatanB145SC.dat')
       OPEN (unit=3,STATUS='unknown',FILE='radmass.dat')
       OPEN (unit=1,STATUS='unknown',FILE='profil.dat')

C     IM = NUMBER OF EQUATIONS Y(I)=Pressure, Y(2)=NS Mass and Y(3)=E density

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      IM=5
      IN=IM-1
      
      IK=0

       DO 10 IL=600,1,-1
       FIXEDIL=300
C        IL=FIXEDIL
       PC=1.D0*IL
C        WRITE(*,*) IL
 
       YA(10)=CL
       YA(9)=YDY
       EDC=FED(PC)
       DEDPC=DEDP(PC)

       PCC=PC-2.D0/3.D0*PI*GS*XC*XC*(PC+EDC)*(3.D0*PC+EDC)
       MCC=4.D0*PI*XC*XC*XC*EDC/(3.D0*MSS)
       ALPCC=-2.D0/3.D0*PI*GS*XC*XC*(3.D0*PC+EDC)
       
       EDENTC=EDC
       PRESSTC=PCC
       ALPHAPC=(GS/(XC*XC))
     &   *(MCC*MSS+4.D0*PI*XC*XC*XC*PRESSTC)
       PRESSPC=-(EDC+PCC)*ALPHAPC
       MASSTPC=4.D0*PI*XC*XC*EDENTC
       CNST2=(EDC+PCC)*(1.D0-DEDPC)
       SIGMAP=-YDY*(2.D0*GS*MSS*MCC/XC*XC)*PCC
     &       +YDY*(2.D0*GS*MSS*MCC/XC)*PRESSPC
     &       +YDY*(2.D0*GS*PCC/XC)*MASSTPC
       LMDCC=1.D0/3.D0*(-2.D0*CNST2*GS*PI*XC*XC*
     &      (3.D0*PCC*PCC+4.D0*PC*EDC+EDC*EDC)) 
     &       +2.D0*(EDC+PCC)*SIGMAP
     
C        PRINT*, SIGMAP
C        STOP
       IK=0
       
       
 29   Y(1)=PCC
      Y(2)=MCC
      Y(3)=ALPCC
      Y(4)=LMDCC
      
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      
      P0=Y(1)
 
      Y(5)=FED(P0)
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D0
      NS=32
      !XL=30.0D3

      H=PU/NS
c     XP should be larger than XC in order to avoid the unphysical behavior
c      near center!!    
      XP=1.0D0
      HH=H/(2.0D0)

      IF (XP.LT.XC) THEN
            WRITE(*,*) "XP=",XP," is NOT larger than XC=",XC
            STOP
      END IF

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH

C    COMPUTE K1, L1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K2, L2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
            
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      


         YA(5)=FED(P0)

         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K3, L3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  
            YA(5)=FED(P0)          

         XA=XM


         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K4, L4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)
         
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

 

            YA(5)=FED(P0)
     

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)
         
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  


          Y(5)=FED(P0)
          
       END DO
       
       XA=XP
       PRESS=Y(1)
       MASST=Y(2)
       ALPHA=Y(3)
       LAMBDA=Y(4)
       EDEN=Y(5)

      
      IF (IL.EQ.FIXEDIL .AND. IK.EQ.1) THEN 
            WRITE(1,*)(XP/1.D3),PRESS,MASST,ALPHA,LAMBDA,EDEN
      ENDIF
     

       PS=Y(1)
       PMIN=1.0D-8

      IF (PS .GT. PMIN  ) GOTO 28
      
      
C       WRITE(*,*)ALPCC,LMDCC,ALPHA,LAMBDA
      IF (ABS(ALPHA).GT.1.D-3) THEN
C             PRINT *,"PRINT 1"
            ALPCC=ALPCC-ALPHA
            GOTO 29
      ENDIF
      IF (ABS(LAMBDA).GT.1.D-2) THEN
C             PRINT *,"PRINT 2"
            LMDCC=LMDCC-LAMBDA
            GOTO 29
      ENDIF
      IF (IK.EQ.0) THEN
C             PRINT *,"PRINT 3"
            IK=IK+1
C             PRINT *,"IK=",IK
            GOTO 29
      ENDIF
      IF (IK.EQ.1) THEN 
C             PRINT *,"PRINT 4"
C             WRITE(*,*)ALPCC,LMDCC,ALPHA,LAMBDA
            WRITE(*,*)PC,(EDC/1.D3),(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP 
            WRITE(3,*)PC,(EDC/1.D3),(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP 
      ENDIF

c      WRITE(2,*)IL,(XP/1.D3), (Y(I),I=1,IM)
C         WRITE(*,*)PC,(EDC/1.D3),(XP/1.D3),Y(2),
C      &  2.0D0*GS*Y(2)*MSS/XP  
 
   

  10   CONTINUE
      
      
      END

      SUBROUTINE FUNCT(EK,J,YA,XA,H)
C     *********************************************************
C     DEFINES TOV EQUATIONS 
C
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(4,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
c--------------------------------------------------------------

      CT=YA(10)
      EDEN=YA(5)
      PRESS=YA(1)
      MASST=YA(2)
      ALPHA=YA(3)
      LAMBDA=YA(4)
      
      YDY=YA(9)
      SIGMA=YDY*(2.D0*GS*MSS*MASST/XA)*PRESS
      
      EDENT=EDEN+2.D0*PI*GS*CT
     &        *(LAMBDA-1.5D0*(EDEN+PRESS)**2
     &        +2.D0*SIGMA*(EDEN+PRESS))

      PRESST=PRESS-2.D0*PI*GS*CT
     &        *(LAMBDA+0.5D0*(EDEN+PRESS)**2
     &        -2.D0*SIGMA*(EDEN+PRESS))
      
      ALPHAP=(GS*MASST*MSS/(XA*XA))
     &   /(1.D0-2.D0*GS*MASST*MSS/XA)
     &   *(1.D0+4.D0*PI*XA*XA*XA/(MASST*MSS)*PRESST)
      
      PRESSP=-(EDEN+PRESS)*ALPHAP
     &       -2.D0*SIGMA/XA
     
      MASSTP=4.D0*PI*XA*XA*EDENT
      
      SIGMAP=-YDY*(2.D0*GS*MSS*MASST/XA*XA)*PRESS
     &       +YDY*(2.D0*GS*MSS*MASST/XA)*PRESSP
     &       +YDY*(2.D0*GS*PRESS/XA)*MASSTP
      
      LAMBDAP=PRESSP*((EDEN+PRESS)
     &       *(1.D0-DEDP(PRESS))
     &       +2.D0*SIGMA*(DEDP(PRESS)-1.D0))
     &       +8.D0*SIGMA*(EDEN+PRESS-SIGMA)/XA
     &       +2.D0*SIGMAP*(EDEN+PRESS)
  
      EK(J,1)=H*PRESSP
      
      EK(J,2)=H/MSS*MASSTP
      
      EK(J,3)=H*ALPHAP
      
      EK(J,4)=H*LAMBDAP
     
C       WRITE(*,*)CT,EDEN,PRESS,MASST,ALPHA,LAMBDA,SIGMA
C       WRITE(*,*)EK(J,1),EK(J,2),EK(J,3),EK(J,4)
C       STOP

      RETURN
      END

 
C-------------------------------------------------
C     MIT BAG B=145 MeV^4
c-------------------------------------------------
  
      
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      HC  = 197.327D0
      B=(145.D0)**4  ! PC until 800
      !B=(185.D0)**4  ! PC until 1600
      BMF=B/(HC*HC*HC)
      FED= 3.D0*P0+4.D0*BMF
   
       RETURN
       END
c-----------------------------------------------------------------------  
      
C-----------------------------------------------------------------------
C     DE/DP AS A FUNCTION OF (P0) for NS (ANTO'S VERSION)
C-----------------------------------------------------------------------
      FUNCTION DEDP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL FED
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
   
   
      DEDP = (FED(x4)-8.D0*FED(x2)+8.D0*FED(x1)-FED(x3))/(12.D0*h)
c      WRITE(*,*)xa,DEDP
      RETURN
      END
 
c---------------------------------------------------------------------    

C-----------------------------------------------------------------------
C     DS/DP AS A FUNCTION OF (P0) for NS (ANTO'S VERSION)
C-----------------------------------------------------------------------
      FUNCTION DSDP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL SIG
      h  = 1.D-5
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      
      DSDP=(SIG(x4)-8.D0*SIG(x2)+8.D0*SIG(x1)-SIG(x3))/(12.D0*h)
      
      RETURN
      END
 
c---------------------------------------------------------------------    

      FUNCTION SIG(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
      SIG = 0.D0
   
       RETURN
       END
