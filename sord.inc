

C********************************************************************
C If real is used use nb = 1 zero = 0.0E+0.0 (CRAY CONVEX)
C If real*8 is used use nb = 2 zero = 0.0D+0.0
C If Complex*16 is used use nb = 4 zero = (0.0,0.0)
C If Complex is used use nb = 2 zero = (0.0,0.0) (CRAY CONVEX)
C********************************************************************






         implicit real*8 (a-h,o-z)
c implicit double precision (a-h,o-z)
         integer ipr
         parameter (ipr = 2)
c parameter (zero = 0.0d0)

