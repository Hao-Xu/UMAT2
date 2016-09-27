      SUBROUTINE GEIGEN( N,A,WR,WI,Z,FV1,IV1)
C=======================================================================
C   25/01/77            MEMBER NAME  GEIGEN   (SO)          FORTRAN
C   Input
C	N - dimension
C	A - input array
C   Output
C	A   - array of eigenvectors
C       WR  - array of eigenvalues (real)
C	WI  - array of eigenvalues (imagin)
C	Z   - working array
C	FV1 - working array
C	IV1 - working array
C======================================================================
CCCCCCCOMMON/ KONTR/KTRW(80)
      include 'sord.inc'
C
      INTEGER*4 IV1(N)
      REAL*8 A(N,N),WR(N),WI(N),Z(N,N),FV1(N)
CCC
CCC   SUBROUTINES FROM EISPECK-PACKAGE
      CALL BALANC(N,N,A, IS1,IS2,FV1)
      CALL ELMHES(N,N,IS1,IS2,A,IV1)
      CALL GEIELTRAN(N,N,IS1,IS2,A,IV1,Z)
      CALL HQR2(N,N,IS1,IS2,A,WR,WI,Z,IERR)
CG    IF(IERR.NE.0)  CALL WRITEI('IERR',1,IERR,1)
      CALL BALBAK(N,N,IS1,IS2,FV1,N,Z)
C
CCC---
CG    IF(KTRW(4).EQ.1)   CALL WRITED('WR  ',1,WR,N)
CG    IF(KTRW(4).EQ.1)   CALL WRITED('WI  ',1,WI,N)
      RETURN
      END
C
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
C===============================================================================
C===============================================================================
      include 'sord.inc'
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL*8 A(NM,N),SCALE(N)
      REAL*8 C,F,G,R,S,B2,RADIX
      REAL*8 DABS
      LOGICAL NOCONV
C
C     ***** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C          THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION
C     ******
      RADIX=16
C
      B2=RADIX*RADIX
      K=1
      L=N
      GO TO 100
C     ***** IN-LINE PROCEDURE FOR ROW AND
C          COLUMN EXCHANGE*****
20    SCALE(M)=J
      IF(J.EQ.M) GO TO 50
C
       DO 30 I=1,L
      F=A(I,J)
      A(I,J)=A(I,M)
      A(I,M)=F
30    CONTINUE
C
       DO  40 I=K,N
      F=A(J,I)
      A(J,I)=A(M,I)
      A(M,I)=F
40    CONTINUE
C
50    GO TO (80,130),IEXC
C     *****SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C          AND PUSH THEM DOWN*****
80    IF(L.EQ.1) GO TO 280
      L=L-1
C     *****FOR J=L STEP -1 UNTIL 1 DO -- *****
100   DO 120 JJ=1,L
      J=L+1-JJ
C
       DO 110 I=1,L
      IF(I.EQ.J) GO TO 110
      IF(A(J,I).NE.0.0) GO TO 120
110   CONTINUE
C
      M=L
      IEXC=1
      GO TO 20
120   CONTINUE
C
      GO TO 140
C     *****SEARCH FOR COLUMS ISOLATING AN EIGENVALUE
C          AND PUSH THEM LEFT*****
130   K=K+1
C
140   DO 170 J=K,L
C
       DO 150 I=K,L
      IF(I.EQ.J) GO TO 150
      IF(A(I,J).NE.0.0) GO TO 170
150   CONTINUE
C
      M=K
      IEXC=2
      GO TO 20
170   CONTINUE
C     *****NOW BALANCE THE SUBMATRIX IN ROWS K TO L*****
       DO 180 I=K,L
180   SCALE(I)=1.0
C     *****ITERATIVE LOOP FOR NORM REDUCTION*****
190   NOCONV=.FALSE.
C
       DO 270 I=K,L
      C=0.0
      R=0.0
C
       DO 200 J=K,L
      IF(J.EQ.I) GO TO 200
      C=C+DABS(A(J,I))
      R=R+DABS(A(I,J))
200   CONTINUE
C
      G=R/RADIX
      F=1.0
      S=C+R
210   IF(C.GE.G) GO TO 220
      F=F*RADIX
      C=C*B2
      GO TO 210
220   G=R*RADIX
230   IF(C.LT.G) GO TO 240
      F=F/RADIX
      C=C/B2
      GO TO 230
C     *****NOW BALANCE*****
240   IF((C+R)/F.GE.0.95*S) GO TO 270
      G=1.0/F
      SCALE(I)=SCALE(I)*F
      NOCONV=.TRUE.
C
       DO 250 J=K,N
250   A(I,J)=A(I,J)*G
C
       DO 260 J=1,L
260   A(J,I)=A(J,I)*F
C
270   CONTINUE
C
      IF(NOCONV) GO TO 190
C
280   LOW=K
      IGH=L
      RETURN
      END
C
C   02/03/76            MEMBER NAME  ELMHES   (QUELLE)      FORTRAN
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
C=================================================================================
C=================================================================================
      include 'sord.inc'
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,NM1,MP1
      REAL*8 A(NM,N)
      REAL*8 X,Y
      REAL*8 DABS
      INTEGER INT(IGH)
C
      LA=IGH-1
      KP1=LOW+1
      IF(LA.LT.KP1) GO TO 200
C
       DO 180 M=KP1,LA
      MM1=M-1
      X=0.0
      I=M
C
       DO 100 J=M,IGH
      IF(DABS(A(J,MM1)).LE.DABS(X)) GO TO 100
      X=A(J,MM1)
      I=J
100   CONTINUE
C
      INT(M)=I
      IF(I.EQ.M) GO TO 130
C     *****INTERCHANGE ROWS AND COLUMNS OF A*****
       DO 110 J=MM1,N
      Y=A(I,J)
      A(I,J)=A(M,J)
      A(M,J)=Y
110   CONTINUE
C
       DO 120 J=1,IGH
      Y=A(J,I)
      A(J,I)=A(J,M)
      A(J,M)=Y
120   CONTINUE
C     *****END INTERCHANGE*****
130   IF(X.EQ.0.0) GO TO 180
      MP1=M+1
C
       DO 160 I=MP1,IGH
      Y=A(I,MM1)
      IF(Y.EQ.0.0) GO TO 160
      Y=Y/X
      A(I,MM1)=Y
C
       DO 140 J=M,N
140   A(I,J)=A(I,J)-Y*A(M,J)
C
       DO 150 J=1,IGH
150   A(J,M)=A(J,M)+Y*A(J,I)
C
160   CONTINUE
C
180   CONTINUE
C
200   RETURN
      END
C   02/03/76            MEMBER NAME  ELTRAN   (QUELLE)      FORTRAN
      SUBROUTINE GEIELTRAN(NM,N,LOW,IGH,A,INT,Z)
C======================================================================================
C-======================================================================================
      include 'sord.inc'
      INTEGER I,J,N,KL,NM,MP,NM1,IGH,LOW,MP1
      REAL*8 A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
C
C     *****INITIALIZE Z TO IDENTITY MATRIX*****
       DO 80 I=1,N
C
       DO 60 J=1,N
60    Z(I,J)=0.0
C
      Z(I,I)=1.0
80    CONTINUE
C
      KL=IGH-LOW-1
      IF(KL.LT.1) GO TO 200
C     *****FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- *****
       DO 140 MM=1,KL
      MP=IGH-MM
      MP1=MP+1
C
       DO 100 I=MP1,IGH
100   Z(I,MP)=A(I,MP-1)
C
      I=INT(MP)
      IF(I.EQ.MP) GO TO 140
C
       DO 130 J=MP,IGH
      Z(MP,J)=Z(I,J)
      Z(I,J)=0.0
130   CONTINUE
C
      Z(I,MP)=1.0
140   CONTINUE
C
200   RETURN
      END
C   02/03/76            MEMBER NAME  HQR2     (QUELLE)      FORTRAN
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C=====================================================================================
C=====================================================================================
      include 'sord.inc'
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITS,LOW,MP2,ENM2,IERR
      REAL*8 H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL*8 P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
      REAL*8 DSQRT,DABS,DSIGN
      INTEGER MIN0
      LOGICAL NOTLAS
      COMPLEX*16 Z3
      COMPLEX*16 DCMPLX
      REAL*8 T3(2)
      EQUIVALENCE(Z3,T3(1))
C
C     *****MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C          THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C     *****
      MACHEP= 1.0D-16
C
      IERR=0
C     *****STORE ROOTS ISOLATED BY BALANC*****
       DO 50 I=1,N
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 50
      WR(I)=H(I,I)
      WI(I)=0.0
50    CONTINUE
C
      EN=IGH
      T=0.0
C     *****SEARCH FOR NEXT EIGENVALUES*****
60    IF(EN.LT.LOW) GO TO 340
      ITS=0
      NA=EN-1
      ENM2=NA-1
C     *****LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C          FOR L=EN STEP -1 UNTIL LOW DO -- *****
70    DO 80 LL=LOW,EN
      L=EN+LOW-LL
      IF(L.EQ.LOW) GO TO 100
      IF(DABS(H(L,L-1)).LE.MACHEP*(DABS(H(L-1,L-1))
     X    +DABS(H(L,L)))) GO TO 100
80    CONTINUE
C     *****FORM SHIFT*****
100   X=H(EN,EN)
      IF(L.EQ.EN) GO TO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF(L.EQ.NA) GO TO 280
      IF(ITS.EQ.30) GO TO 1000
      IF(ITS.NE.10.AND.ITS.NE.20) GO TO 130
C     *****FORM EXCEPTIONAL SHIFT*****
      T=T+X
C
       DO 120 I=LOW,EN
120   H(I,I)=H(I,I)-X
C
      S=DABS(H(EN,NA))+DABS(H(NA,ENM2))
      X=0.75*S
      Y=X
      W=-0.4375*S*S
130   ITS=ITS+1
C     LOOK FOR TWO CONSECUTIVE SMALL
C          SUB-DIAGONAL ELEMENTS.
C          FOR M=EN-2 STEP -1 UNTIL L DO -- *****
       DO 140 MM=L,ENM2
      M=ENM2+L-MM
      ZZ=H(M,M)
      R=X-ZZ
      S=Y-ZZ
      P=(R*S-W)/H(M+1,M)+H(M,M+1)
      Q=H(M+1,M+1)-ZZ-R-S
      R=H(M+2,M+1)
      S=DABS(P)+DABS(Q)+DABS(R)
      P=P/S
      Q=Q/S
      R=R/S
      IF(M.EQ.L) GO TO 150
      IF(DABS(H(M,M-1))*(DABS(Q)+DABS(R)).LE.MACHEP*DABS(P)
     X  *(DABS(H(M-1,M-1))+DABS(ZZ)+DABS(H(M+1,M+1)))) GO TO 150
140   CONTINUE
C
150   MP2=M+2
C
       DO 160 I=MP2,EN
      H(I,I-2)=0.0
      IF(I.EQ.MP2) GO TO 160
      H(I,I-3)=0.0
160   CONTINUE
C     *****DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C          COLUMNS M TO EN*****
       DO 260 K=M,NA
      NOTLAS=K.NE.NA
      IF(K.EQ.M) GO TO 170
      P=H(K,K-1)
      Q=H(K+1,K-1)
      R=0.0
      IF(NOTLAS) R=H(K+2,K-1)
      X=DABS(P)+DABS(Q)+DABS(R)
      IF(X.EQ.0.0) GO TO 260
      P=P/X
      Q=Q/X
      R=R/X
170   S=DSIGN(DSQRT(P*P+Q*Q+R*R),P)
      IF(K.EQ.M) GO TO 180
      H(K,K-1)=-S*X
      GO TO 190
180   IF(L.NE.M) H(K,K-1)=-H(K,K-1)
190   P=P+S
      X=P/S
      Y=Q/S
      ZZ=R/S
      Q=Q/P
      R=R/P
C     *****ROW MODIFICATION*****
       DO 210 J=K,N
      P=H(K,J)+Q*H(K+1,J)
      IF(.NOT.NOTLAS) GO TO 200
      P=P+R*H(K+2,J)
      H(K+2,J)=H(K+2,J)-P*ZZ
200   H(K+1,J)=H(K+1,J)-P*Y
      H(K,J)=H(K,J)-P*X
210   CONTINUE
C
      J=MIN0(EN,K+3)
C     *****COLUMN MODIFICATION*****
       DO 230 I=1,J
      P=X*H(I,K)+Y*H(I,K+1)
      IF(.NOT.NOTLAS) GO TO 220
      P=P+ZZ*H(I,K+2)
      H(I,K+2)=H(I,K+2)-P*R
220   H(I,K+1)=H(I,K+1)-P*Q
      H(I,K)=H(I,K)-P
230   CONTINUE
C     *****ACCUMULATE TRANSFORMATIONS*****
       DO 250 I=LOW,IGH
      P=X*Z(I,K)+Y*Z(I,K+1)
      IF(.NOT.NOTLAS) GO TO 240
      P=P+ZZ*Z(I,K+2)
      Z(I,K+2)=Z(I,K+2)-P*R
240   Z(I,K+1)=Z(I,K+1)-P*Q
      Z(I,K)=Z(I,K)-P
250   CONTINUE
C
260   CONTINUE
C
      GO TO 70
C     *****ONE ROOT FOUND*****
270   H(EN,EN)=X+T
      WR(EN)=H(EN,EN)
      WI(EN)=0.0
      EN=NA
      GO TO 60
C     *****TWO ROOTS FOUND*****
280   P=(Y-X)/2.0
      Q=P*P+W
      ZZ=DSQRT(DABS(Q))
      H(EN,EN)=X+T
      X=H(EN,EN)
      H(NA,NA)=Y+T
      IF(Q.LT.0.0) GO TO 320
C     *****REAL PAIR*****
      ZZ=P+DSIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF(ZZ.NE.0.0) WR(EN)=X-W/ZZ
      WI(NA)=0.0
      WI(EN)=0.0
      X=H(EN,NA)
      R=DSQRT(X*X+ZZ*ZZ)
      P=X/R
      Q=ZZ/R
C     *****ROW MODIFICATION*****
       DO 290 J=NA,N
      ZZ=H(NA,J)
      H(NA,J)=Q*ZZ+P*H(EN,J)
      H(EN,J)=Q*H(EN,J)-P*ZZ
290   CONTINUE
C     *****COLUMN MODIFICATION*****
       DO 300 I=1,EN
      ZZ=H(I,NA)
      H(I,NA)=Q*ZZ+P*H(I,EN)
      H(I,EN)=Q*H(I,EN)-P*ZZ
300   CONTINUE
C     *****ACCUMULATE TRANSFORMATIONS*****
       DO 310 I=LOW,IGH
      ZZ=Z(I,NA)
      Z(I,NA)=Q*ZZ+P*Z(I,EN)
      Z(I,EN)=Q*Z(I,EN)-P*ZZ
310   CONTINUE
C
      GO TO 330
C     *****COMPLEX PAIR*****
320   WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
330   EN=ENM2
      GO TO 60
C     *****ALL ROOTS FOUND.BACKSUBSTITUTE TO FIND
C          VECTORS OF UPPER TRIANGULAR FORM*****
340   NORM=0.0
      K=1
C
       DO 360 I=1,N
C
       DO 350 J=K,N
350   NORM=NORM+DABS(H(I,J))
C
      K=I
360   CONTINUE
C
      IF(NORM.EQ.0.0) GO TO 1001
C     *****FOR EN=N STEP -1 UNTIL 1 DO -- *****
       DO 800 NN=1,N
      EN=N+1-NN
      P=WR(EN)
      Q=WI(EN)
      NA=EN-1
      IF(Q) 710,600,800
C     *****REAL VECTOR*****
600   M=EN
      H(EN,EN)=1.0
      IF(NA.EQ.0) GO TO 800
C     *****FOR I=EN-1 STEP -1 UNTIL 1 DO --*****
       DO 700 II=1,NA
      I=EN-II
      W=H(I,I)-P
      R=H(I,EN)
      IF(M.GT.NA) GO TO 620
C
       DO 610 J=M,NA
610   R=R+H(I,J)*H(J,EN)
C
620   IF(WI(I).GE.0.0) GO TO 630
      ZZ=W
      S=R
      GO TO 700
630   M=I
      IF(WI(I).NE.0.0) GO TO 640
      T=W
      IF(W.EQ.0.0) T=MACHEP*NORM
      H(I,EN)=-R/T
      GO TO 700
C************* SOLVE REAL EQUATIONS *******************
640   X=H(I,I+1)
      Y=H(I+1,I)
      Q=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
      T=(X*S-ZZ*R)/Q
      H(I,EN)=T
      IF(DABS(X).LE.DABS(ZZ)) GOTO 650
      H(I+1,EN)=(-R-W*T)/X
      GOTO 700
650   H(I+1,EN)=(-S-Y*T)/ZZ
700   CONTINUE
C     *****END REAL VECTOR*****
      GO TO 800
C     *****COMPLEX VECTOR*****
710   M=NA
C     *****LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C          EIGENVECTOR MATRIX IS TRIANGULAR*****
      IF(DABS(H(EN,NA)).LE.DABS(H(NA,EN))) GO TO 720
      H(NA,NA)=Q/H(EN,NA)
      H(NA,EN)=-(H(EN,EN)-P)/H(EN,NA)
      GO TO 730
720   Z3=DCMPLX(0.D+0,-H(NA,EN))/DCMPLX(H(NA,NA)-P,Q)
      H(NA,NA)=T3(1)
      H(NA,EN)=T3(2)
730   H(EN,NA)=0.0
      H(EN,EN)=1.0
      ENM2=NA-1
      IF(ENM2.EQ.0) GO TO 800
C
       DO 790 II=1,ENM2
      I=NA-II
      W=H(I,I)-P
      RA=0.0
      SA=H(I,EN)
C
       DO 760 J=M,NA
      RA=RA+H(I,J)*H(J,NA)
      SA=SA+H(I,J)*H(J,EN)
760   CONTINUE
C
      IF(WI(I).GE.0.0) GO TO 770
      ZZ=W
      R=RA
      S=SA
      GO TO 790
770   M=I
      IF(WI(I).NE.0.0) GO TO 780
      Z3=DCMPLX(-RA,-SA)/DCMPLX(W,Q)
      H(I,NA)=T3(1)
      H(I,EN)=T3(2)
      GO TO 790
C     *****SOLVE COMPLEX EQUATIONS*****
780   X=H(I,I+1)
      Y=H(I+1,I)
      VR=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
      VI=(WR(I)-P)*2.0*Q
      IF(VR.EQ.0.0.AND.VI.EQ.0.0) VR=MACHEP*NORM
     X  *(DABS(W)+DABS(Q)+DABS(X)+DABS(Y)+DABS(ZZ))
      Z3=DCMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/DCMPLX(VR,VI)
      H(I,NA)=T3(1)
      H(I,EN)=T3(2)
      IF(DABS(X).LE.DABS(ZZ)+DABS(Q)) GO TO 785
      H(I+1,NA)=(-RA-W*H(I,NA)+Q*H(I,EN))/X
      H(I+1,EN)=(-SA-W*H(I,EN)-Q*H(I,NA))/X
      GO TO 790
785   Z3=DCMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN))/DCMPLX(ZZ,Q)
      H(I+1,NA)=T3(1)
      H(I+1,EN)=T3(2)
790   CONTINUE
C     *****END COMPLEX VECTOR*****
800   CONTINUE
C     *****END BACK SUBSTITUTION.
C          VECTORS OF ISOLATED ROOTS*****
       DO 840 I=1,N
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 840
C
       DO 820 J=1,N
820   Z(I,J)=H(I,J)
C
840   CONTINUE
C     *****MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C          VECTORS OF ORIGINAL FULL MATRIX.
C          FOR J=N STEP -1 UNTIL LOW DO --*****
       DO 880 JJ=LOW,N
      J=N+LOW-JJ
      M=MIN0(J,IGH)
C
       DO 880 I=LOW,IGH
      ZZ=0.0
C
       DO 860 K=LOW,M
860   ZZ=ZZ+Z(I,K)*H(K,J)
C
      Z(I,J)=ZZ
880   CONTINUE
C
      GO TO 1001
C     *****SET ERROR -- NO CONVERGENCE TO AN
C          EIGENVALUE AFTER 30 ITERATIONS*****
1000  IERR=EN
1001  RETURN
      END
C   02/03/76            MEMBER NAME  BALBAK   (QUELLE)      FORTRAN
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
C========================================================================
C========================================================================
      include 'sord.inc'
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL*8 SCALE(N),Z(NM,M)
      REAL*8 S,SUM
C
      IF(IGH.EQ.LOW) GO TO 120
C
       DO 110 I=LOW,IGH
      S=SCALE(I)
C     *****LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C          IF THE FOREGOING STATEMENT IS REPLACED BY
C          S=1.0/SCALE(I).*****
       DO 100 J=1,M
100   Z(I,J)=Z(I,J)*S
C
110   CONTINUE
C     *****FOR I=LOW-1 STEP -1 UNTIL 1,
C          IGH+1 STEP 1 UNTIL N DO -- *****
120   DO 140 II=1,N
      I=II
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 140
      IF(I.LT.LOW) I=LOW-II
      K=SCALE(I)
      IF(K.EQ.I) GO TO 140
C
       DO 130 J=1,M
      S=Z(I,J)
      Z(I,J)=Z(K,J)
      Z(K,J)=S
130   CONTINUE
C
140   CONTINUE
C
       DO  170 J=1,M
      SUM=0.0
       DO  150 I=1,N
      SUM=SUM+Z(I,J)*Z(I,J)
150   CONTINUE
      SUM=DSQRT(SUM)
       DO  160 I=1,N
      Z(I,J)=Z(I,J)/SUM
160   CONTINUE
170   CONTINUE
C
      RETURN
      END

