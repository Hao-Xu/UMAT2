      SUBROUTINE UMAT(STRESS,STRAN,DSTRAN,
     X    OMEGA,EPSID,fd0,
     X    PROPS,NTENS)
C
      include 'sord.inc'
C    
      DIMENSION STRESS(*),STRAN(*),DSTRAN(*),
     X OMEGA(NTENS),EPSID(NTENS),
     X PROPS(*)
      DIMENSION DSIG(NTENS),AMATDOM(NTENS,NTENS),
     3 DG_DY(NTENS),DF_DSIG(NTENS),TEMP(NTENS)
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
C
      PARAMETER (FTOL=1.D-6,PI=3.1415926526D0,ITMAX=250)
      PARAMETER (tiny=1.D-3)
C
C     MATERIAL PARAMETERS
C
      E0 = PROPS(1)
      ENU0 = PROPS(2)
      A1 = PROPS(3)
      A2 = PROPS(4)
      A3 = PROPS(5)
      A4 = PROPS(6)
      C0 = PROPS(7)
      C1 = PROPS(8)
      ALPHA = PROPS(9)
      OMEGA_11 = PROPS(10)
      OMEGA_22 = PROPS(11)
      OMEGA_33 = PROPS(12)
      OMEGA_12 = PROPS(13)
      OMEGA_23 = PROPS(14)
      OMEGA_13 = PROPS(15)
C
C====================================================================== 
C     FD = fd0
C
C      IF (TIME(2).EQ.0) THEN
C        OMEGA(1,1) = OMEGA_11
C        OMEGA(2,2) = OMEGA_22
C        OMEGA(3,3) = OMEGA_33
C      END IF
C
C --- CALCULATION OF THE EFFECTIVE ELASTIC STIFFNESS
C
      DO I=1,NTENS
        IF(I<=3)THEN
          STRAN(I)=STRAN(I)
          DSTRAN(I)=DSTRAN(I)
        ELSE
          STRAN(I)=STRAN(I)/2.D0
          DSTRAN(I)=DSTRAN(I)/2.D0
        ENDIF
      ENDDO
      dstrain=0
      DO I=1,NTENS
         dstrain=dstrain+DSTRAN(I)*DSTRAN(I)
      ENDDO
      dstrain=sqrt(dstrain)
      CALL DMAT(AMATDOM,OMEGA,NTENS,
     X          E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA)
C
C     ELASTICAL TRIAL 
C
       DO I=1,NTENS
          DSIG(I)=0.D0
         DO J=1,NTENS
           IF(J<=3)THEN
             DSIG(I)=DSIG(I)+AMATDOM(I,J)*DSTRAN(J)
           ELSE
             DSIG(I)=DSIG(I)+2.D0*AMATDOM(I,J)*DSTRAN(J)
           ENDIF
         END DO
       END DO

         do i=1,NTENS
           STRESS(i) = STRESS(i) + DSIG(i)
         enddo
      
C    
      IOPTFD = 1
      CALL FDDP(FD1,STRESS,OMEGA,NTENS,
     X E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA,IOPTFD)
C
C     DAMAGE OR NOT?
C
      IF (FD1.GT.FTOL) THEN
C
C     DMAGE OCCURS
C
        INC=0
        
        FDT=FD1
      DO WHILE (FDT.GT.0.D0 .AND.((FDT/FD1).GT.FTOL.AND.INC.LT.ITMAX))  
c--------------------------------------------------------------------
c       iopt = 0   no update of AMATDOM
c            = 1      update of AMATDOM
c-------------------------------------------------------------------
        iopt = 0
        if(INC.eq.0) iopt=0
        CALL FD_LAM(OMEGA,STRESS,NTENS,AMATDOM,
     X              E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA,
     X              H0,HP,DG_DY,DF_DSIG,TEMP,iopt)
      
        XL=FDT/(H0-HP)
       if (dstrain.gt.tiny)     XL=XL/1.5D0      

       DO I=1,NTENS
         OMEGA(I)=OMEGA(I)+XL*DG_DY(I)
         EPSID(I)=EPSID(I)+XL*DF_DSIG(I)
         STRESS(I)=STRESS(I)-XL*TEMP(I)
       ENDDO
       
c     write(*,"(A10,4e12.3)") 'TEMP = ',TEMP(1),TEMP(2),TEMP(3),XL
C
      IOPTFD = 2
         CALL FDDP(FDT,STRESS,OMEGA,NTENS,
     X          E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA,IOPTFD)
         INC=INC+1
      END DO
      IF(INC.GE.ITMAX) write (7,*) 'no convernece'
      END IF 
C
C
C
C --- UPDATING STATE VARIABLES
C
      IOPTFD = 3
      CALL FDDP(FD2,STRESS,OMEGA,NTENS,
     X E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA,IOPTFD)
C
      fd0=FD2
c     END IF 
      RETURN
      END
C
C
CC ====================================================================
C                              SUBROUTINES 
C
C
      SUBROUTINE DMAT(AMATDOM,OMEGA,NTENS,
     X                E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA)
CC ====================================================================
C
C                              D M A T 
C
C ====================================================================
      include 'sord.inc'
C
      DIMENSION AMATDOM(NTENS,NTENS),OMEGA(NTENS)
      DIMENSION AMATS(NTENS,NTENS)
     
C
      DATA ONE,TWO,HALF / 1.0D0,2.0D0,0.5D0 /
C
C      INITIALIZING TENSORS AND MATRIX
C
C       call dclear(AMATS,NTENS*NTENS)
C       call dclear(AMATDOM,NTENS*NTENS)
       DO I=1,NTENS
         DO J=1,NTENS
           AMATS(I,J)=0.0D0
           AMATDOM(I,J)=0.0D0
         END DO
       END DO
C
       TROMEGA = OMEGA(1)+OMEGA(2)+OMEGA(3)
       B1 = (1.D0+ENU0)/E0/2.D0
       B2 = ENU0/E0
C
       COE1=2.D0
       COE2=4.D0
C
       AMATS(1,1)=2.D0*B1-B2+2.D0*A1*TROMEGA+2.D0*(A2+A3)*OMEGA(1)
     X +2.D0*A4*TROMEGA
       AMATS(1,2)=-B2+2.D0*A1*TROMEGA+A3*(OMEGA(1)+OMEGA(2))
       AMATS(1,3)=-B2+2.D0*A1*TROMEGA+A3*(OMEGA(3)+OMEGA(1))
       AMATS(1,4)=COE1*(A2*OMEGA(4)+A3*OMEGA(4))
C
       AMATS(2,1)=AMATS(1,2)
       AMATS(2,2)=2.D0*B1-B2+2.D0*A1*TROMEGA+2.D0*(A2+A3)*OMEGA(2)
     X +2.D0*A4*TROMEGA
       AMATS(2,3)=-B2+2.D0*A1*TROMEGA+A3*(OMEGA(3)+OMEGA(2))
       AMATS(2,4)=COE1*(A2*OMEGA(4)+A3*OMEGA(4))
C
       AMATS(3,1)=AMATS(1,3)
       AMATS(3,2)=AMATS(2,3)
       AMATS(3,3)=2.D0*B1-B2+2.D0*A1*TROMEGA+2.D0*(A2+A3)*OMEGA(3)
     X +2.D0*A4*TROMEGA
       AMATS(3,4)=COE1*(A3*OMEGA(4))
C
       AMATS(4,1)=AMATS(1,4)
       AMATS(4,2)=AMATS(2,4)       
       AMATS(4,3)=AMATS(3,4)       
       AMATS(4,4)=COE2*(B1+0.5D0*A2*(OMEGA(1)+OMEGA(2))+A4*TROMEGA)
C 
       
       IF (NTENS.EQ.6) THEN
         AMATS(1,5)=COE1*(A3*OMEGA(5))
         AMATS(1,6)=COE1*(A2*OMEGA(6)+A3*OMEGA(6))
C
         AMATS(2,5)=COE1*(A2*OMEGA(5)+A3*OMEGA(5))
         AMATS(2,6)=COE1*(A3*OMEGA(6))
C
         AMATS(3,5)=COE1*(A2*OMEGA(5)+A3*OMEGA(5))
         AMATS(3,6)=COE1*(A2*OMEGA(6)+A3*OMEGA(6))
C
         AMATS(4,5)=COE2*0.5D0*A2*OMEGA(6)
         AMATS(4,6)=COE2*0.5D0*A2*OMEGA(5)
C
         AMATS(5,1)=AMATS(1,5)
         AMATS(5,2)=AMATS(2,5)       
         AMATS(5,3)=AMATS(3,5)       
         AMATS(5,4)=AMATS(4,5)
         AMATS(5,5)=COE2*(B1+0.5D0*A2*(OMEGA(3)+OMEGA(2))+A4*TROMEGA)
         AMATS(5,6)=COE2*0.5D0*A2*OMEGA(4)
C
         AMATS(6,1)=AMATS(1,6)
         AMATS(6,2)=AMATS(2,6)       
         AMATS(6,3)=AMATS(3,6)       
         AMATS(6,4)=AMATS(4,6)
         AMATS(6,5)=AMATS(5,6)
         AMATS(6,6)=COE2*(B1+0.5D0*A2*(OMEGA(3)+OMEGA(1))+A4*TROMEGA)

       ENDIF
C
        CALL INVERSE(AMATS,NTENS,NTENS,AMATDOM)
c       CALL invert_matrix(AMATDOM,AMATS,NTENS)
C
      RETURN
      END
C
      SUBROUTINE MAT1_MAT2(NTENS,VECTOR,TENSOR,FACT)
C ====================================================================
C 
C                MAT1_MAT2 : VECTOR TO TENOR  
C
C ====================================================================
      include 'sord.inc'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,3
        DO J=1,3
          TENSOR(I,J)=0.0D0
        END DO
      END DO
C
      TENSOR(1 , 1) = VECTOR( 1 )
      TENSOR(2 , 2) = VECTOR( 2 )
      TENSOR(3 , 3) = VECTOR( 3 )
      TENSOR(1 , 2) = VECTOR( 4 )*FACT
      TENSOR(2 , 1) = TENSOR(1 , 2)
      
      IF (NTENS.EQ.6) THEN  
        TENSOR(2 , 3) = VECTOR( 5 )*FACT
        TENSOR(3 , 2) = TENSOR(2 , 3)
        TENSOR(1 , 3) = VECTOR( 6 )*FACT
        TENSOR(3 , 1) = TENSOR(1 , 3)
      ELSEIF (NTENS.EQ.4) THEN
        TENSOR(2 , 3) = 0.D0
        TENSOR(3 , 2) = 0.D0
        TENSOR(1 , 3) = 0.D0
        TENSOR(3 , 1) = 0.D0
      ENDIF
  
C
      RETURN
      END
C
      SUBROUTINE FDDP(FD,SIGMA,OMEGA,NTENS,
     X E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA,IOPTFD)
C=========================================================================
C
C       Damage Yield Function (Pseudo-Drucker-Pager)
C
C=========================================================================
      include 'sord.inc'
C
      DIMENSION SIGMA(NTENS),OMEGA(NTENS)
      DIMENSION SIGSIG(NTENS),P_1(NTENS,NTENS),P_1Y(NTENS),SIJ(NTENS),
     X E1(NTENS),YD1(NTENS)
C
      DATA ONE,TWO,HALF / 1.0D0,2.0D0,0.5D0 /
C
      DO I=1,NTENS
          IF(I<=3)THEN
            E1(I)=1.D0
          ELSE
            E1(I)=0.D0
          ENDIF
      END DO
C
C
      TRSIGMA = SIGMA(1)+SIGMA(2)+SIGMA(3)
      TROMEGA = OMEGA(1)+OMEGA(2)+OMEGA(3)
c     write(6,*)' TRSIGMA = ', TRSIGMA
c     write(6,*)' ',OMEGA(1), ' ',OMEGA(2),' ',OMEGA(3)
C
      CALL Aik_Bkj(SIGMA,SIGMA,SIGSIG,NTENS)
C
      TRSIGSIG =SIGSIG(1)+SIGSIG(2)+SIGSIG(3)
C
      DO I=1,NTENS
          YD1(I)= A1*TRSIGMA**2.D0*E1(I)+A2*SIGSIG(I)+
     1        A3*TRSIGMA*SIGMA(I)+A4*TRSIGSIG*E1(I)
      END DO
C
      CALL MATP_1(P_1,SIGMA,NTENS)
C
      DO I=1,NTENS
            P_1Y(I)=0.D0
        DO J=1,NTENS
          IF (J<=3) THEN
            P_1Y(I)=P_1Y(I)+P_1(I,J)*YD1(J)
          ELSE
            P_1Y(I)=P_1Y(I)+2.D0*P_1(I,J)*YD1(J)
          ENDIF  
        ENDDO
      ENDDO
C
      TRY = P_1Y(1)+P_1Y(2)+P_1Y(3)
C
      DO I= 1,NTENS
          SIJ(I) = P_1Y(I)-1.D0/3.D0*TRY*E1(I)
      END DO
C
      SS=0.
      DO I=1,NTENS
          IF (I<=3) THEN
            SS=SS+SIJ(I)*SIJ(I)
          ELSE
            SS=SS+2.D0*SIJ(I)*SIJ(I)
          ENDIF  
      ENDDO
      FD=ALPHA*TRY-C0-C1*TROMEGA

c     IF(FD>0.)THEN
c       write (7,*) 'GT_DSID:  Stresses are in tension'
c       write(*,*)'Nelm = ',Nelm
c       write(*,"(A10,I8)")'IOPT = ',IOPTFD
c       write(*,"(A10,4e12.3)")'SIGMA = ',SIGMA(1),SIGMA(2),SIGMA(3),
c    X  SIGMA(4)
c       write(*,*)'P1Y = '
c       write(*,"(4e12.3)")P_1Y(1),P_1Y(2),P_1Y(3),P_1Y(4)
c       write(*,*)'TRY = ',TRY
c       call pend('test terminated')
c     ENDIF  
      FD = SQRT(HALF*SS)+ALPHA*TRY-C0-C1*TROMEGA   !THE SIGN BEFORE ALPHA IS "+" DUE TO MECHANICAL CONVENTION
c     write(6,*) FD, ' ', SS, " ", TRY
C
      RETURN
      END
C
      SUBROUTINE MATP_1(P_1,SIGMA,NTENS)
C======================================================================
C
C                          MATP_1	  
C
C=====================================================================
      include 'sord.inc'
C
      DIMENSION P_1(NTENS,NTENS),SIGMA(NTENS)
      DIMENSION AN1(3),AN2(3),AN3(3),S(3),IV1(3),SIGMA1(3,3)
      DIMENSION ANAN1(NTENS),ANAN2(NTENS),ANAN3(NTENS)
      REAL*8 EIGVALR(3),EIGVALI(3),EIGVEC(3,3),FV1(3),SIGMA0(3,3)
C
      PARAMETER (FTOL=1.D-6,ONE=1.D0,N3=3)

C

      CALL MAT1_MAT2(NTENS,SIGMA,SIGMA1,ONE)
       DO I=1,3
         AN1(I)=0.D0
         AN2(I)=0.D0
         AN3(I)=0.D0
         S(I)=0.D0
         EIGVALR(I)=0.D0
         EIGVALI(I)=0.D0
         FV1(I)=0.D0
         IV1(I)=0.D0
         DO J=1,3
           SIGMA0(I,J)=SIGMA1(I,J)
           EIGVEC(I,J)=0.D0
         END DO
       END DO 
C
C
      DO I=1,NTENS
        DO J=1,NTENS
           P_1(I,J)=0.D0
        ENDDO
      ENDDO
C
      CALL GEIGEN(N3,SIGMA0,EIGVALR,EIGVALI,EIGVEC,FV1,IV1)
C 
C
      DO I=1,3
        IF (ABS(EIGVALR(I)).LT.FTOL) EIGVALR(I)=0.D0
      END DO
C
      DO I=1,3
        AN1(I)=EIGVEC(I,1)
        AN2(I)=EIGVEC(I,2)
        AN3(I)=EIGVEC(I,3)
      END DO
C
      DO I=1,3 
        IF (EIGVALR(I).GE.0.D0) THEN
          S(I)=1.D0
        ELSE
          S(I)=-1.D0
        END IF
      END DO
C	  
      DO I=1,3
         ANAN1(I)=AN1(I)*AN1(I)
         ANAN2(I)=AN2(I)*AN2(I)
         ANAN3(I)=AN3(I)*AN3(I)
      ENDDO
C
      ANAN1(4)=AN1(1)*AN1(2)
      ANAN2(4)=AN2(1)*AN2(2)
      ANAN3(4)=AN3(1)*AN3(2)
C 
      IF(NTENS.EQ.6)THEN
        ANAN1(5)=AN1(2)*AN1(3)
        ANAN2(5)=AN2(2)*AN2(3)
        ANAN3(5)=AN3(2)*AN3(3)
C
        ANAN1(6)=AN1(3)*AN1(1)
        ANAN2(6)=AN2(3)*AN2(1)
        ANAN3(6)=AN3(3)*AN3(1)
      ENDIF
        
       DO I=1,NTENS
         DO J=1,NTENS
           P_1(I,J)=S(1)*ANAN1(I)*ANAN1(J)+
     1              S(2)*ANAN2(I)*ANAN2(J)+
     2              S(3)*ANAN3(I)*ANAN3(J)
         END DO
       END DO
     
      RETURN
      END
C
C
      SUBROUTINE Aik_Bkj(A,B,C,NTENS)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      include 'sord.inc'
C
      DIMENSION A(NTENS),B(NTENS),C(NTENS)
C
      IF (NTENS.EQ.4) THEN
        C(1)=A(1)*B(1)+A(4)*B(4)
        C(2)=A(4)*B(4)+A(2)*B(2)
        C(3)=A(3)*B(3)
        C(4)=A(1)*B(4)+A(4)*B(2)
      ELSEIF(NTENS.EQ.6) THEN
        C(1)=A(1)*B(1)+A(4)*B(4)+A(6)*B(6)
        C(2)=A(4)*B(4)+A(2)*B(2)+A(5)*B(5)
        C(3)=A(6)*B(6)+A(5)*B(5)+A(3)*B(3)
        C(4)=A(1)*B(4)+A(4)*B(2)+A(6)*B(5)
        C(5)=A(4)*B(6)+A(2)*B(5)+A(5)*B(3)
        C(6)=A(1)*B(6)+A(4)*B(5)+A(6)*B(3)
      ENDIF 
C
      RETURN
      END
C
C
      SUBROUTINE MATP_2(P_2,SIGMA,NTENS)
C====================================================================================================
C
C                          MATP_2	  
C
C=====================================================================================================
C
      include 'sord.inc'
C
      DIMENSION P_2(NTENS,NTENS),SIGMA(NTENS)
      DIMENSION ANAN1(NTENS),ANAN2(NTENS),ANAN3(NTENS)
      DIMENSION AN1(3),AN2(3),AN3(3),S(3),IV1(3)
      REAL*8 EIGVALR(3),EIGVALI(3),EIGVEC(3,3),FV1(3),SIGMA0(3,3)
C
      PARAMETER (FTOL=1.D-6,ONE=1.,N3=3)
      CALL MAT1_MAT2(NTENS,SIGMA,SIGMA0,ONE)
       DO I=1,3
         AN1(I)=0.D0
         AN2(I)=0.D0
         AN3(I)=0.D0
         S(I)=0.D0
         EIGVALR(I)=0.D0
         EIGVALI(I)=0.D0
         FV1(I)=0.D0
         IV1(I)=0.D0
         DO J=1,3
           EIGVEC(I,J)=0.D0
         END DO
       END DO 
C
       DO I=1,NTENS
         DO J=1,NTENS
           P_2(I,J)=0.D0
         ENDDO
       ENDDO
C
      CALL GEIGEN(N3,SIGMA0,EIGVALR,EIGVALI,EIGVEC,FV1,IV1)
C
C      DO I=1,3
C        IF (ABS(EIGVALR(I)).LT.(1.D-10)) EIGVALR(I)=0.D0
C      END DO
      DO I=1,3
        AN1(I)=EIGVEC(I,1)
        AN2(I)=EIGVEC(I,2)
        AN3(I)=EIGVEC(I,3)
      END DO
      DO I=1,3 
        RES=EIGVALR(I)-MIN(EIGVALR(1),EIGVALR(2),EIGVALR(3))
        write (6,*) 'RES=',RES
        IF (RES.GT.0.0D0) THEN     !CHECK
          S(I)=1.0D0       
        ELSE
          S(I)=0.0D0
        END IF
        IF (abs(RES)<FTOL) S(I)=0.0D0
      END DO
        write(6,*)'strss=',EIGVALR(1),' ',EIGVALR(2),' ',EIGVALR(3)
C	  
      DO I=1,3
         ANAN1(I)=AN1(I)*AN1(I)
         ANAN2(I)=AN2(I)*AN2(I)
         ANAN3(I)=AN3(I)*AN3(I)
      ENDDO
C
      ANAN1(4)=AN1(1)*AN1(2)
      ANAN2(4)=AN2(1)*AN2(2)
      ANAN3(4)=AN3(1)*AN3(2)
C 
      IF(NTENS.EQ.6)THEN
        ANAN1(5)=AN1(2)*AN1(3)
        ANAN2(5)=AN2(2)*AN2(3)
        ANAN3(5)=AN3(2)*AN3(3)
C
        ANAN1(6)=AN1(3)*AN1(1)
        ANAN2(6)=AN2(3)*AN2(1)
        ANAN3(6)=AN3(3)*AN3(1)
      ENDIF
        
       DO I=1,NTENS
         DO J=1,NTENS
           P_2(I,J)=S(1)*ANAN1(I)*ANAN1(J)+
     1              S(2)*ANAN2(I)*ANAN2(J)+
     2              S(3)*ANAN3(I)*ANAN3(J)
         END DO
       END DO
c     
      RETURN
      END
C
C
      SUBROUTINE DY_DSIGFUN(SIGMA,DY_DSIG,NTENS,
     X E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA)
C=======================================================================
C                                                                      *
C                           DY_DSIGFUN                                 *
C                                                                      *
C=======================================================================
      include 'sord.inc'
C
      DIMENSION SIGMA(NTENS),DY_DSIG(NTENS,NTENS)
C
C
      TRSIG=SIGMA(1)+SIGMA(2)+SIGMA(3)
C
        DY_DSIG(1,1)=2.*A1*TRSIG+2.*A2*SIGMA(1)+A3*(SIGMA(1)+TRSIG)
     X   +2.*A4*SIGMA(1)
        DY_DSIG(1,2)=2.*A1*TRSIG+A3*SIGMA(1)+2.*A4*SIGMA(2)
        DY_DSIG(1,3)=2.*A1*TRSIG+A3*SIGMA(1)+2.*A4*SIGMA(3)
        DY_DSIG(1,4)=A2*SIGMA(4)+2.*A4*SIGMA(4)
C
        DY_DSIG(2,1)=2.*A1*TRSIG+A3*SIGMA(2)+2.*A4*SIGMA(1)
        DY_DSIG(2,2)=2.*A1*TRSIG+2.*A2*SIGMA(2)+A3*(SIGMA(2)+TRSIG)
     X   +2.*A4*SIGMA(2)
        DY_DSIG(2,3)=2.*A1*TRSIG+A3*SIGMA(2)+2.*A4*SIGMA(3)
        DY_DSIG(2,4)=A2*SIGMA(4)+2.*A4*SIGMA(4)
C
        DY_DSIG(3,1)=2.*A1*TRSIG+A3*SIGMA(3)+2.*A4*SIGMA(1)
        DY_DSIG(3,2)=2.*A1*TRSIG+A3*SIGMA(3)+2.*A4*SIGMA(2)
        DY_DSIG(3,3)=2.*A1*TRSIG+2.*A2*SIGMA(3)+A3*(SIGMA(3)+TRSIG)
     X   +2.*A4*SIGMA(3)
        DY_DSIG(3,4)=2.*A4*SIGMA(4)
C
        DY_DSIG(4,1)=A2*SIGMA(4)+A3*SIGMA(4)
        DY_DSIG(4,2)=A2*SIGMA(4)+A3*SIGMA(4)
        DY_DSIG(4,3)=A3*SIGMA(4)
        DY_DSIG(4,4)=0.5*A2*(SIGMA(1)+SIGMA(2))+0.5*A3*TRSIG
C
        IF(NTENS.EQ.6)THEN
        DY_DSIG(1,5)=2.*A4*SIGMA(5)
        DY_DSIG(1,6)=A2*SIGMA(6)+2.*A4*SIGMA(6)
C
        DY_DSIG(2,5)=A2*SIGMA(5)+2.*A4*SIGMA(5)
        DY_DSIG(2,6)=2.*A4*SIGMA(6)
C
        DY_DSIG(3,5)=A2*SIGMA(5)+2.*A4*SIGMA(5)
        DY_DSIG(3,6)=A2*SIGMA(6)+2.*A4*SIGMA(6)
C
        DY_DSIG(4,5)=0.5*A2*SIGMA(6)
        DY_DSIG(4,6)=0.5*A2*SIGMA(5)
C
        DY_DSIG(5,1)=A3*SIGMA(5)
        DY_DSIG(5,2)=A2*SIGMA(5)+A3*SIGMA(5)
        DY_DSIG(5,3)=A2*SIGMA(5)+A3*SIGMA(5)
        DY_DSIG(5,4)=DY_DSIG(4,5)
        DY_DSIG(5,5)=0.5*A2*(SIGMA(2)+SIGMA(3))+0.5*A3*TRSIG
        DY_DSIG(5,6)=0.5*A2*SIGMA(4)
C
        DY_DSIG(6,1)=A2*SIGMA(6)+A3*SIGMA(6)
        DY_DSIG(6,2)=A3*SIGMA(6)
        DY_DSIG(6,3)=A2*SIGMA(6)+A3*SIGMA(6)
        DY_DSIG(6,4)=DY_DSIG(4,6)
        DY_DSIG(6,5)=DY_DSIG(5,6)
        DY_DSIG(6,6)=0.5*A2*(SIGMA(3)+SIGMA(1))+0.5*A3*TRSIG
        ENDIF

      RETURN
      END
C
C
C
      SUBROUTINE FD_LAM(OMEGA,SIGMA,NTENS,AMATDOM,
     X                 E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA,
     X                 H0,HP,DG_DY,DF_DSIG,TEMP,iopt)
C==================================================================
C
C             ITERATION TO SOLVE FOR LAMBDA
C
C===================================================================
      include 'sord.inc'
      DIMENSION OMEGA(NTENS),SIGMA(NTENS),AMATDOM(NTENS,NTENS),
     x DG_DY(NTENS),DF_DSIG(NTENS),TEMP(NTENS)
      DIMENSION E(NTENS),YD1(NTENS),P_1(NTENS,NTENS),
     1 DY_DSIG(NTENS,NTENS),SIGSIG(NTENS),
     2 F1IJ(NTENS),F2IJ(NTENS),P1YD1(NTENS),F2P2(NTENS),EP1(NTENS),
     3 EEP1(NTENS,NTENS),DF_DY(NTENS),DS_DOMEGA(NTENS,NTENS,NTENS),
     4 P_3(NTENS,NTENS),DF_DOME(NTENS),F1P3(NTENS),
     5 P_2(NTENS,NTENS),ZEROS(NTENS),TEMP1(NTENS),TEMP2(NTENS,NTENS)
C
C
      DATA ONE,TWO,HALF / 1.0D0,2.0D0,0.5D0 /
C
      DO I=1,NTENS
          IF(I<=3)THEN
            E(I)=1.
          ELSE
            E(I)=0.
          ENDIF
      END DO
C
      TRSIGMA = SIGMA(1)+SIGMA(2)+SIGMA(3)
      TROMEGA = OMEGA(1)+OMEGA(2)+OMEGA(3)
C
      CALL Aik_Bkj(SIGMA,SIGMA,SIGSIG,NTENS)
C
      TRSIGSIG =SIGSIG(1)+SIGSIG(2)+SIGSIG(3)
C
      DO I=1,NTENS
          ZEROS(I)=0.
      END DO
      
      DO I=1,NTENS
          YD1(I)= A1*TRSIGMA**2.*E(I)+A2*SIGSIG(I)+
     1        A3*TRSIGMA*SIGMA(I)+A4*TRSIGSIG*E(I)
      ENDDO
C
      CALL DY_DSIGFUN(SIGMA,DY_DSIG,NTENS,
     X E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA)
      CALL MATP_1(P_1,SIGMA,NTENS)
C
      CALL MATP_2(P_2,SIGMA,NTENS)
      DO I=1,NTENS
        P1YD1(I)=0.
        DO J=1,NTENS
          IF (J<=3) THEN
            P1YD1(I)=P1YD1(I)+P_1(I,J)*YD1(J)
          ELSE
            P1YD1(I)=P1YD1(I)+2.*P_1(I,J)*YD1(J)
          ENDIF  
        ENDDO
      ENDDO
      P1YD1E=P1YD1(1)+P1YD1(2)+P1YD1(3)
c     write(6,*) ' P1YD1E = ',P1YD1E
c     write(6,*) ' YD1= ',(YD1(i),i=1,3)
C
      DO I=1,NTENS
        F2IJ(I)=0.
        DO J=1,NTENS
          IF (J<=3) THEN
            F2IJ(I)=F2IJ(I)+P_2(I,J)*YD1(J)
          ELSE
            F2IJ(I)=F2IJ(I)+2.*P_2(I,J)*YD1(J)
          ENDIF  
        ENDDO
      ENDDO
C
      DO I=1,NTENS
        F2P2(I)=0.
        DO J=1,NTENS
          IF (J<=3) THEN
            F2P2(I)=F2P2(I)+P_2(J,I)*F2IJ(J)
          ELSE
            F2P2(I)=F2P2(I)+2.*P_2(J,I)*F2IJ(J)
          ENDIF  
        ENDDO
      ENDDO
C      
      F2F2=0.
      DO I=1,NTENS
        IF (I<=3) THEN
          F2F2=F2F2+F2IJ(I)*F2IJ(I)
        ELSE
          F2F2=F2F2+2.*F2IJ(I)*F2IJ(I)
        ENDIF  
      ENDDO
C
      DO I=1,NTENS
        EP1(I)=0.
        DO J=1,NTENS
          IF (J<=3) THEN
            EP1(I)=EP1(I)+P_1(J,I)*E(J)
          ELSE
            EP1(I)=EP1(I)+2.*P_1(J,I)*E(J)
          ENDIF  
        ENDDO
      ENDDO
C     
      DO I=1,NTENS
        DO J=1,NTENS
          EEP1(I,J)=E(I)*EP1(J)
        ENDDO
      ENDDO
C
      DO I=1,NTENS
        DO J=1,NTENS
          P_3(I,J)=P_1(I,J)-1./3.*EEP1(I,J)
        END DO
      END DO
C
      DO I=1,NTENS
          F1IJ(I)=P1YD1(I)-1./3.*P1YD1E*E(I)
          DF_DOME(I)=-C1*E(I)
          IF(F2F2.EQ.0.) THEN
              DG_DY(I)=0.
          ELSE  
             DG_DY(I)=F2P2(I)/SQRT(2.0*F2F2)
          END IF
      END DO
c     write(6,*) ' f2ij= ', (F2IJ(i),i=1,3)
C      
      DO I=1,NTENS
        F1P3(I)=0.
        DO J=1,NTENS
          IF (J<=3) THEN
            F1P3(I)=F1P3(I)+F1IJ(J)*P_3(J,I)
          ELSE
            F1P3(I)=F1P3(I)+2.*F1IJ(J)*P_3(J,I)
          ENDIF  
        ENDDO
      ENDDO
C      
      F1F1=0.
      DO I=1,NTENS
        IF (I<=3) THEN
          F1F1=F1F1+F1IJ(I)*F1IJ(I)
        ELSE
          F1F1=F1F1+2.*F1IJ(I)*F1IJ(I)
        ENDIF  
      ENDDO
C
      DO I=1,NTENS
        DF_DY(I)=F1P3(I)/SQRT(2.*F1F1)+ALPHA*EP1(I)
      END DO
C     
      DO I=1,NTENS
        DF_DSIG(I)=0.D0
        DO J=1,NTENS
          IF (J<=3) THEN
            DF_DSIG(I)=DF_DSIG(I)+DF_DY(J)*DY_DSIG(J,I)
          ELSE
            DF_DSIG(I)=DF_DSIG(I)+2.D0*DF_DY(J)*DY_DSIG(J,I)
          ENDIF  
        ENDDO
      ENDDO
C
      
      HP=0.D0
      DO I=1,NTENS
        IF(I<=3)THEN
        HP=HP+DF_DOME(I)*DG_DY(I)
        ELSE
        HP=HP+2.D0*DF_DOME(I)*DG_DY(I)
        ENDIF 
      ENDDO
c
      H0=0.D0
c
      CALL AMATDSDOMEGA(DS_DOMEGA,NTENS,
     X                 E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA)
c
      DO I=1,NTENS
        DO J=1,NTENS
          TEMP2(I,J)=0.D0
          DO K=1,NTENS
            IF (K<=3) THEN
            TEMP2(I,J)=TEMP2(I,J)+SIGMA(K)*DS_DOMEGA(K,I,J)
            ELSE
            TEMP2(I,J)=TEMP2(I,J)+2.D0*SIGMA(K)*DS_DOMEGA(K,I,J)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      DO I=1,NTENS
        TEMP1(I)=0.D0
        DO J=1,NTENS
         IF (J<=3)THEN
           TEMP1(I)=TEMP1(I)+TEMP2(I,J)*DG_DY(J)
         ELSE
           TEMP1(I)=TEMP1(I)+2.0D0*TEMP2(I,J)*DG_DY(J)
         ENDIF
        ENDDO
      ENDDO
C
      DO I=1,NTENS
         TEMP1(I)=TEMP1(I)+DF_DSIG(I)
      ENDDO
     
      if(iopt.eq.1) then
      CALL DMAT(AMATDOM,OMEGA,NTENS,
     X          E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA)
      endif
      DO I=1,NTENS
        TEMP(I)=0.D0
        DO J=1,NTENS
          IF(J<=3)THEN
          TEMP(I)=TEMP(I)+AMATDOM(I,J)*TEMP1(J)
          ELSE
          TEMP(I)=TEMP(I)+2.D0*AMATDOM(I,J)*TEMP1(J)
          ENDIF 
        ENDDO
      ENDDO
      DO I=1,NTENS
        IF(I<=3)THEN
          H0=H0+DF_DSIG(I)*TEMP(I)
        ELSE
          H0=H0+2.D0*DF_DSIG(I)*TEMP(I)
        ENDIF 
      ENDDO
      RETURN
      END
C
      SUBROUTINE AMATDSDOMEGA(DS_DO,NTENS,
     X                 E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA)
      include 'sord.inc'
      DIMENSION DS_DO(NTENS,NTENS,NTENS)


C ------- page 1 -------
      DS_DO(1,1,1)=2.*A1+2.*A2+2.*A3+2.*A4
      DS_DO(1,2,1)=2.*A1+A3
      DS_DO(1,3,1)=2.*A1+A3  
      DS_DO(1,4,1)=0.
C
      DS_DO(2,1,1)=DS_DO(1,2,1)
      DS_DO(2,2,1)=2.*A1+2.*A4
      DS_DO(2,3,1)=2.*A1
      DS_DO(2,4,1)=0.
C 
      DS_DO(3,1,1)=DS_DO(1,3,1)
      DS_DO(3,2,1)=DS_DO(2,3,1)
      DS_DO(3,3,1)=2.*A1+2.*A4
      DS_DO(3,4,1)=0.
C
      DS_DO(4,1,1)=DS_DO(1,4,1)
      DS_DO(4,2,1)=DS_DO(2,4,1)
      DS_DO(4,3,1)=DS_DO(3,4,1)
      DS_DO(4,4,1)=0.5*A2+A4
C
      IF(NTENS.EQ.6)THEN
      DS_DO(1,5,1)=0.
      DS_DO(1,6,1)=0.
      DS_DO(2,5,1)=0.
      DS_DO(2,6,1)=0.
      DS_DO(3,5,1)=0.
      DS_DO(3,6,1)=0.
      DS_DO(4,5,1)=0.
      DS_DO(4,6,1)=0.
      DS_DO(5,1,1)=DS_DO(1,5,1)
      DS_DO(5,2,1)=DS_DO(2,5,1)
      DS_DO(5,3,1)=DS_DO(3,5,1)
      DS_DO(5,4,1)=DS_DO(4,5,1)
      DS_DO(5,5,1)=A4
      DS_DO(5,6,1)=0.
C
      DS_DO(6,1,1)=DS_DO(1,6,1)
      DS_DO(6,2,1)=DS_DO(2,6,1)
      DS_DO(6,3,1)=DS_DO(3,6,1)
      DS_DO(6,4,1)=DS_DO(4,6,1)
      DS_DO(6,5,1)=DS_DO(5,6,1)
      DS_DO(6,6,1)=0.5*A2+A4
      ENDIF
C
C ------- page 2 -------
      DS_DO(1,1,2)=2.*A1+2.*A4
      DS_DO(1,2,2)=2.*A1+A3
      DS_DO(1,3,2)=2.*A1  
      DS_DO(1,4,2)=0.
C
      DS_DO(2,1,2)=DS_DO(1,2,2)
      DS_DO(2,2,2)=2.*A1+2.*A2+2.*A3+2.*A4
      DS_DO(2,3,2)=2.*A1+A3
      DS_DO(2,4,2)=0.
C 
      DS_DO(3,1,2)=DS_DO(1,3,2)
      DS_DO(3,2,2)=DS_DO(2,3,2)
      DS_DO(3,3,2)=2.*A1+2.*A4
      DS_DO(3,4,2)=0.
C
      DS_DO(4,1,2)=DS_DO(1,4,2)
      DS_DO(4,2,2)=DS_DO(2,4,2)
      DS_DO(4,3,2)=DS_DO(3,4,2)
      DS_DO(4,4,2)=0.5*A2+A4
C
      IF(NTENS.EQ.6)THEN
      DS_DO(1,5,2)=0.
      DS_DO(1,6,2)=0.
      DS_DO(2,5,2)=0.
      DS_DO(2,6,2)=0.
      DS_DO(3,5,2)=0.
      DS_DO(3,6,2)=0.
      DS_DO(4,5,2)=0.
      DS_DO(4,6,2)=0.
      DS_DO(5,1,2)=DS_DO(1,5,2)
      DS_DO(5,2,2)=DS_DO(2,5,2)
      DS_DO(5,3,2)=DS_DO(3,5,2)
      DS_DO(5,4,2)=DS_DO(4,5,2)
      DS_DO(5,5,2)=0.5*A2+A4
      DS_DO(5,6,2)=0.
C
      DS_DO(6,1,2)=DS_DO(1,6,2)
      DS_DO(6,2,2)=DS_DO(2,6,2)
      DS_DO(6,3,2)=DS_DO(3,6,2)
      DS_DO(6,4,2)=DS_DO(4,6,2)
      DS_DO(6,5,2)=DS_DO(5,6,2)
      DS_DO(6,6,2)=0.5*A2+A4
      ENDIF
C
C ------- page 3 -------
      DS_DO(1,1,3)=2.*A1+2.*A4
      DS_DO(1,2,3)=2.*A1
      DS_DO(1,3,3)=2.*A1+A3  
      DS_DO(1,4,3)=0.
C
      DS_DO(2,1,3)=DS_DO(1,2,3)
      DS_DO(2,2,3)=2.*A1+2.*A4
      DS_DO(2,3,3)=2.*A1+A3
      DS_DO(2,4,3)=0.
C 
      DS_DO(3,1,3)=DS_DO(1,3,3)
      DS_DO(3,2,3)=DS_DO(2,3,3)
      DS_DO(3,3,3)=2.*A1+2.*A2+2.*A3+2.*A4
      DS_DO(3,4,3)=0.
C
      DS_DO(4,1,3)=DS_DO(1,4,3)
      DS_DO(4,2,3)=DS_DO(2,4,3)
      DS_DO(4,3,3)=DS_DO(3,4,3)
      DS_DO(4,4,3)=A4
C
      IF(NTENS.EQ.6)THEN
      DS_DO(1,5,3)=0.
      DS_DO(1,6,3)=0.
      DS_DO(2,5,3)=0.
      DS_DO(2,6,3)=0.
      DS_DO(3,5,3)=0.
      DS_DO(3,6,3)=0.
      DS_DO(4,5,3)=0.
      DS_DO(4,6,3)=0.
      DS_DO(5,1,3)=DS_DO(1,5,3)
      DS_DO(5,2,3)=DS_DO(2,5,3)
      DS_DO(5,3,3)=DS_DO(3,5,3)
      DS_DO(5,4,3)=DS_DO(4,5,3)
      DS_DO(5,5,3)=0.5*A2+A4
      DS_DO(5,6,3)=0.
C
      DS_DO(6,1,3)=DS_DO(1,6,3)
      DS_DO(6,2,3)=DS_DO(2,6,3)
      DS_DO(6,3,3)=DS_DO(3,6,3)
      DS_DO(6,4,3)=DS_DO(4,6,3)
      DS_DO(6,5,3)=DS_DO(5,6,3)
      DS_DO(6,6,3)=0.5*A2+A4
      ENDIF
C
C ------- page 4 -------
      DS_DO(1,1,4)=0.
      DS_DO(1,2,4)=0.
      DS_DO(1,3,4)=0.  
      DS_DO(1,4,4)=0.5*A2+0.5*A3
C
      DS_DO(2,1,4)=DS_DO(1,2,4)
      DS_DO(2,2,4)=0.
      DS_DO(2,3,4)=0.
      DS_DO(2,4,4)=0.5*A2+0.5*A3
C 
      DS_DO(3,1,4)=DS_DO(1,3,4)
      DS_DO(3,2,4)=DS_DO(2,3,4)
      DS_DO(3,3,4)=0.
      DS_DO(3,4,4)=0.5*A3
C
      DS_DO(4,1,4)=DS_DO(1,4,4)
      DS_DO(4,2,4)=DS_DO(2,4,4)
      DS_DO(4,3,4)=DS_DO(3,4,4)
      DS_DO(4,4,4)=0.
C
      IF(NTENS.EQ.6)THEN
      DS_DO(1,5,4)=0.
      DS_DO(1,6,4)=0.
      DS_DO(2,5,4)=0.
      DS_DO(2,6,4)=0.
      DS_DO(3,5,4)=0.
      DS_DO(3,6,4)=0.
      DS_DO(4,5,4)=0.
      DS_DO(4,6,4)=0.
      DS_DO(5,1,4)=DS_DO(1,5,4)
      DS_DO(5,2,4)=DS_DO(2,5,4)
      DS_DO(5,3,4)=DS_DO(3,5,4)
      DS_DO(5,4,4)=DS_DO(4,5,4)
      DS_DO(5,5,4)=0.
      DS_DO(5,6,4)=0.25*A2
C
      DS_DO(6,1,4)=DS_DO(1,6,4)
      DS_DO(6,2,4)=DS_DO(2,6,4)
      DS_DO(6,3,4)=DS_DO(3,6,4)
      DS_DO(6,4,4)=DS_DO(4,6,4)
      DS_DO(6,5,4)=DS_DO(5,6,4)
      DS_DO(6,6,4)=0.
      ENDIF
C
      IF(NTENS.EQ.6)THEN
C ------- page 5 -------
      DS_DO(1,1,5)=0.
      DS_DO(1,2,5)=0.
      DS_DO(1,3,5)=0.  
      DS_DO(1,4,5)=0.
      DS_DO(1,5,5)=0.5*A3
      DS_DO(1,6,5)=0.
C
      DS_DO(2,1,5)=DS_DO(1,2,5)
      DS_DO(2,2,5)=0.
      DS_DO(2,3,5)=0.
      DS_DO(2,4,5)=0.
      DS_DO(2,5,5)=0.5*A2+0.5*A3
      DS_DO(2,6,5)=0.
C 
      DS_DO(3,1,5)=DS_DO(1,3,5)
      DS_DO(3,2,5)=DS_DO(2,3,5)
      DS_DO(3,3,5)=0.
      DS_DO(3,4,5)=0.
      DS_DO(3,5,5)=0.5*A2+0.5*A3
      DS_DO(3,6,5)=0.
C
      DS_DO(4,1,5)=DS_DO(1,4,5)
      DS_DO(4,2,5)=DS_DO(2,4,5)
      DS_DO(4,3,5)=DS_DO(3,4,5)
      DS_DO(4,4,5)=0.
      DS_DO(4,5,5)=0.
      DS_DO(4,6,5)=0.25*A2
C
      DS_DO(5,1,5)=DS_DO(1,5,5)
      DS_DO(5,2,5)=DS_DO(2,5,5)
      DS_DO(5,3,5)=DS_DO(3,5,5)
      DS_DO(5,4,5)=DS_DO(4,5,5)
      DS_DO(5,5,5)=0.
      DS_DO(5,6,5)=0.
C
      DS_DO(6,1,5)=DS_DO(1,6,5)
      DS_DO(6,2,5)=DS_DO(2,6,5)
      DS_DO(6,3,5)=DS_DO(3,6,5)
      DS_DO(6,4,5)=DS_DO(4,6,5)
      DS_DO(6,5,5)=DS_DO(5,6,5)
      DS_DO(6,6,5)=0.
C
C ------- page 6 -------
      DS_DO(1,1,6)=0.
      DS_DO(1,2,6)=0.
      DS_DO(1,3,6)=0.  
      DS_DO(1,4,6)=0.
      DS_DO(1,5,6)=0.
      DS_DO(1,6,6)=0.5*A2+0.5*A3
C
      DS_DO(2,1,6)=DS_DO(1,2,6)
      DS_DO(2,2,6)=0.
      DS_DO(2,3,6)=0.
      DS_DO(2,4,6)=0.
      DS_DO(2,5,6)=0.
      DS_DO(2,6,6)=0.5*A3
C 
      DS_DO(3,1,6)=DS_DO(1,3,6)
      DS_DO(3,2,6)=DS_DO(2,3,6)
      DS_DO(3,3,6)=0.
      DS_DO(3,4,6)=0.
      DS_DO(3,5,6)=0.
      DS_DO(3,6,6)=0.5*A2+0.5*A3
C
      DS_DO(4,1,6)=DS_DO(1,4,6)
      DS_DO(4,2,6)=DS_DO(2,4,6)
      DS_DO(4,3,6)=DS_DO(3,4,6)
      DS_DO(4,4,6)=0.
      DS_DO(4,5,6)=0.25*A2
      DS_DO(4,6,6)=0.
C
      DS_DO(5,1,6)=DS_DO(1,5,6)
      DS_DO(5,2,6)=DS_DO(2,5,6)
      DS_DO(5,3,6)=DS_DO(3,5,6)
      DS_DO(5,4,6)=DS_DO(4,5,6)
      DS_DO(5,5,6)=0.
      DS_DO(5,6,6)=0.
C
      DS_DO(6,1,6)=DS_DO(1,6,6)
      DS_DO(6,2,6)=DS_DO(2,6,6)
      DS_DO(6,3,6)=DS_DO(3,6,6)
      DS_DO(6,4,6)=DS_DO(4,6,6)
      DS_DO(6,5,6)=DS_DO(5,6,6)
      DS_DO(6,6,6)=0.
      ENDIF

      RETURN
      END

