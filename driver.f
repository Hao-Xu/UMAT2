      PROGRAM MAIN
      INCLUDE 'sord.inc'
 
C      REAL*16
      PARAMETER (NTENS=6,NPROPS=12)
      DIMENSION STRESS(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),OMEGA(NTENS),EPSID(NTENS),
     3 PROPS(NPROPS)

      DO I=1,NTENS
         DSTRAN(I) = 0.
         STRAN(I)  = 0.
         STRESS(I) = 0.
         OMEGA(I)  = 0.
         EPSID(I)  = 0.
      ENDDO
      DSTRAN(3) = 0.00002
      PROPS(1)  = 6.8E10
      PROPS(2)  = 0.21
      PROPS(3)  = 1.26E-13
      PROPS(4)  = 3.94E-11
      PROPS(5)  = -1.26E-12
      PROPS(6)  = 2.51E-12
      PROPS(7)  = 0.E6
      PROPS(8)  = 2.2E6
      PROPS(9)  = 0.231
      PROPS(10) = 0.
      PROPS(11) = 0.
      PROPS(12) = 0.

      ISTEP = 1000

      WRITE(*,*) '========================================'
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '=         ABAQUS UMAT DRIVER           ='
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '========================================'

      open(unit=666, file='result.out',status='unknown')
      WRITE(666,6660) 
      DO INC=1,ISTEP
c        WRITE(6,*) '*** Loading increment ', INC, ' ***'
c         WRITE(*,*) 'STRESS = ',(STRESS(I),I=1,NTENS)
c         WRITE(*,*) 'CMNAME = ', CMNAME
c         WRITE(*,*) 'PROPS  = ',(PROPS(I),I=1,NPROPS) 
      CALL UMAT(STRESS,STRAN,DSTRAN,
     2 OMEGA,EPSID,fd0,
     3 PROPS,NTENS)
      DO I=1,NTENS
         STRAN(I)=STRAN(I)+DSTRAN(I)
      enddo
         WRITE(666,'(I5,2X,10(E10.4,2X))') INC, (STRAN(I),I=1,3),
     1 (STRESS(I),I=1,3),(OMEGA(I),I=1,3),fd0 
      END DO
      close(666)
6660  format ('# STEP   STRAN(1)    STRAN(2)    STRAN(3)   STRESS(1)   '
     1 'STRESS(2)   STRESS(1)    OMEGA(1)    OMEGA(2)    OMEGA(3)'
     2 '   fd0')
      END

