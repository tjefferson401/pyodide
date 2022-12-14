* -*- fortran -*-
      SUBROUTINE sBICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, 
     $                    INFO,NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      real   RESID
      INTEGER            NDX1, NDX2
      real   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      real   X( * ), B( * ), WORK( LDW,* )
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW .gt. = max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          .gt.   0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*           .ls.   0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N  .ls.  0
*                   -2: LDW  .ls.  N
*                   -3: Maximum number of iterations ITER  .ls. = 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO  .ls.  BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA  .ls.  BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in func GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*
*     .. Parameters ..
      real             ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, P, PHAT, V, S, SHAT, T, MAXIT,
     $                   NEED1, NEED2
      real   TOL,
     $     RHOTOL, OMEGATOL, sGETBREAK, 
     $     sNRM2
      real   ALPHA, BETA, RHO, RHO1, OMEGA, TMPVAL,
     $     sdot
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Funcs ..
      EXTERNAL           sGETBREAK, sAXPY, sCOPY, 
     $      sdot, sNRM2, sSCAL
*     ..
*     .. Intrinsic Funcs ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      P    = 3
      V    = 4
      T    = 5
      PHAT = 6
      SHAT = 7
      S    = 1
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((T - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((T - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set parameter tolerances.
*
      RHOTOL = sGETBREAK()
      OMEGATOL = sGETBREAK()
*
*     Set initial residual.
*
      CALL sCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( sNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL sMATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using RTLD[] as temp. storage.
*********CALL sCOPY(N, X, 1, WORK(1,RTLD), 1)
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
      ENDIF
*
*****************
 2    CONTINUE
*****************
*
      IF ( sNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30

      CALL sCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform BiConjugate Gradient Stabilized iteration.
*
      ITER = ITER + 1
*     
      RHO = sdot( N, WORK(1,RTLD), 1, WORK(1,R), 1 )
      IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*     
*        Compute vector P.
*
      IF ( ITER.GT.1 ) THEN
         BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
         CALL sAXPY( N, -OMEGA, WORK(1,V), 1, WORK(1,P), 1 )
         CALL sSCAL( N, BETA, WORK(1,P), 1 )
         TMPVAL = ONE
         CALL sAXPY( N, TMPVAL, WORK(1,R), 1, WORK(1,P), 1 )
      ELSE
         CALL sCOPY( N, WORK(1,R), 1, WORK(1,P), 1 )
      ENDIF
*
*        Compute direction adjusting vector PHAT and scalar ALPHA.
*
*********CALL PSOLVE( WORK(1,PHAT), WORK(1,P) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((P    - 1) * LDW) + 1
*     Prepare for return & return
      RLBL = 3
      IJOB = 2
      RETURN
*
*****************
 3    CONTINUE
*****************
*
*********CALL MATVEC( ONE, WORK(1,PHAT), ZERO, WORK(1,V) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((V    - 1) * LDW) + 1
*        Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 4
      IJOB = 1
      RETURN
*
*****************
 4    CONTINUE
*****************
*
      TMPVAL = sdot( N, WORK(1,RTLD), 1, WORK(1,V), 1 )
      IF (TMPVAL.EQ.0) THEN
*        Breakdown
         INFO = -11
         GO TO 20
      ENDIF
      ALPHA = RHO / TMPVAL
*
*        Early check for tolerance.
*
      CALL sAXPY( N, -ALPHA, WORK(1,V), 1, WORK(1,R), 1 )
      CALL sCOPY( N, WORK(1,R), 1, WORK(1,S), 1 )
      IF ( sNRM2( N, WORK(1,S), 1 ).LE.TOL ) THEN
         CALL sAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
         RESID = sNRM2( N, WORK(1,S), 1 )
         GO TO 30
      ELSE
*
*           Compute stabilizer vector SHAT and scalar OMEGA.
*
************CALL PSOLVE( WORK(1,SHAT), WORK(1,S) )
*
         NDX1 = ((SHAT - 1) * LDW) + 1
         NDX2 = ((S    - 1) * LDW) + 1
*     Prepare for return & return
         RLBL = 5
         IJOB = 2
         RETURN
      ENDIF
*
*****************
 5    CONTINUE
*****************
*
************CALL MATVEC( ONE, WORK(1,SHAT), ZERO, WORK(1,T) )
*
      NDX1 = ((SHAT - 1) * LDW) + 1
      NDX2 = ((T    - 1) * LDW) + 1
*           Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 6
      IJOB = 1
      RETURN
*
*****************
 6    CONTINUE
*****************
*
      OMEGA = sdot( N, WORK(1,T), 1, WORK(1,S), 1 ) / 
     $     sdot( N, WORK(1,T), 1, WORK(1,T), 1 )
*
*           Compute new solution approximation vector X.
*
      CALL sAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
      CALL sAXPY( N, OMEGA, WORK(1,SHAT), 1, X, 1 )
*     
*     Compute residual R, check for tolerance.
*
      CALL sAXPY( N, -OMEGA, WORK(1,T), 1, WORK(1,R), 1 )
*
************RESID = DNRM2( N, WORK(1,R), 1 ) / BNRM2 
************IF ( RESID.LE.TOL  ) GO TO 30
*
      NDX1 = NEED1
      NDX2 = NEED2
*     Prepare for resumption & return
      RLBL = 7
      IJOB = 4
      RETURN
*
*****************
 7    CONTINUE
*****************
      IF( INFO.EQ.1 ) GO TO 30
*     
      IF ( ITER.EQ.MAXIT ) THEN
         INFO = 1
         GO TO 20
      ENDIF
*
      IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         GO TO 25
      ELSE
         RHO1 = RHO
         GO TO 10
      ENDIF
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         INFO = -11
      ENDIF
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of BICGSTABREVCOM
*
      END
*     END SUBROUTINE sBICGSTABREVCOM


      SUBROUTINE dBICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, 
     $                    INFO,NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      double precision   RESID
      INTEGER            NDX1, NDX2
      double precision   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      double precision   X( * ), B( * ), WORK( LDW,* )
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW .gt. = max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          .gt.   0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*           .ls.   0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N  .ls.  0
*                   -2: LDW  .ls.  N
*                   -3: Maximum number of iterations ITER  .ls. = 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO  .ls.  BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA  .ls.  BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in func GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*
*     .. Parameters ..
      double precision             ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, P, PHAT, V, S, SHAT, T, MAXIT,
     $                   NEED1, NEED2
      double precision   TOL,
     $     RHOTOL, OMEGATOL, dGETBREAK, 
     $     dNRM2
      double precision   ALPHA, BETA, RHO, RHO1, OMEGA, TMPVAL,
     $     ddot
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Funcs ..
      EXTERNAL           dGETBREAK, dAXPY, dCOPY, 
     $      ddot, dNRM2, dSCAL
*     ..
*     .. Intrinsic Funcs ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      P    = 3
      V    = 4
      T    = 5
      PHAT = 6
      SHAT = 7
      S    = 1
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((T - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((T - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set parameter tolerances.
*
      RHOTOL = dGETBREAK()
      OMEGATOL = dGETBREAK()
*
*     Set initial residual.
*
      CALL dCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( dNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL dMATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using RTLD[] as temp. storage.
*********CALL dCOPY(N, X, 1, WORK(1,RTLD), 1)
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
      ENDIF
*
*****************
 2    CONTINUE
*****************
*
      IF ( dNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30

      CALL dCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform BiConjugate Gradient Stabilized iteration.
*
      ITER = ITER + 1
*     
      RHO = ddot( N, WORK(1,RTLD), 1, WORK(1,R), 1 )
      IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*     
*        Compute vector P.
*
      IF ( ITER.GT.1 ) THEN
         BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
         CALL dAXPY( N, -OMEGA, WORK(1,V), 1, WORK(1,P), 1 )
         CALL dSCAL( N, BETA, WORK(1,P), 1 )
         TMPVAL = ONE
         CALL dAXPY( N, TMPVAL, WORK(1,R), 1, WORK(1,P), 1 )
      ELSE
         CALL dCOPY( N, WORK(1,R), 1, WORK(1,P), 1 )
      ENDIF
*
*        Compute direction adjusting vector PHAT and scalar ALPHA.
*
*********CALL PSOLVE( WORK(1,PHAT), WORK(1,P) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((P    - 1) * LDW) + 1
*     Prepare for return & return
      RLBL = 3
      IJOB = 2
      RETURN
*
*****************
 3    CONTINUE
*****************
*
*********CALL MATVEC( ONE, WORK(1,PHAT), ZERO, WORK(1,V) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((V    - 1) * LDW) + 1
*        Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 4
      IJOB = 1
      RETURN
*
*****************
 4    CONTINUE
*****************
*
      TMPVAL = ddot( N, WORK(1,RTLD), 1, WORK(1,V), 1 )
      IF (TMPVAL.EQ.0) THEN
*        Breakdown
         INFO = -11
         GO TO 20
      ENDIF
      ALPHA = RHO / TMPVAL
*
*        Early check for tolerance.
*
      CALL dAXPY( N, -ALPHA, WORK(1,V), 1, WORK(1,R), 1 )
      CALL dCOPY( N, WORK(1,R), 1, WORK(1,S), 1 )
      IF ( dNRM2( N, WORK(1,S), 1 ).LE.TOL ) THEN
         CALL dAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
         RESID = dNRM2( N, WORK(1,S), 1 )
         GO TO 30
      ELSE
*
*           Compute stabilizer vector SHAT and scalar OMEGA.
*
************CALL PSOLVE( WORK(1,SHAT), WORK(1,S) )
*
         NDX1 = ((SHAT - 1) * LDW) + 1
         NDX2 = ((S    - 1) * LDW) + 1
*     Prepare for return & return
         RLBL = 5
         IJOB = 2
         RETURN
      ENDIF
*
*****************
 5    CONTINUE
*****************
*
************CALL MATVEC( ONE, WORK(1,SHAT), ZERO, WORK(1,T) )
*
      NDX1 = ((SHAT - 1) * LDW) + 1
      NDX2 = ((T    - 1) * LDW) + 1
*           Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 6
      IJOB = 1
      RETURN
*
*****************
 6    CONTINUE
*****************
*
      OMEGA = ddot( N, WORK(1,T), 1, WORK(1,S), 1 ) / 
     $     ddot( N, WORK(1,T), 1, WORK(1,T), 1 )
*
*           Compute new solution approximation vector X.
*
      CALL dAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
      CALL dAXPY( N, OMEGA, WORK(1,SHAT), 1, X, 1 )
*     
*     Compute residual R, check for tolerance.
*
      CALL dAXPY( N, -OMEGA, WORK(1,T), 1, WORK(1,R), 1 )
*
************RESID = DNRM2( N, WORK(1,R), 1 ) / BNRM2 
************IF ( RESID.LE.TOL  ) GO TO 30
*
      NDX1 = NEED1
      NDX2 = NEED2
*     Prepare for resumption & return
      RLBL = 7
      IJOB = 4
      RETURN
*
*****************
 7    CONTINUE
*****************
      IF( INFO.EQ.1 ) GO TO 30
*     
      IF ( ITER.EQ.MAXIT ) THEN
         INFO = 1
         GO TO 20
      ENDIF
*
      IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         GO TO 25
      ELSE
         RHO1 = RHO
         GO TO 10
      ENDIF
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         INFO = -11
      ENDIF
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of BICGSTABREVCOM
*
      END
*     END SUBROUTINE dBICGSTABREVCOM


      SUBROUTINE cBICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, 
     $                    INFO,NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      real   RESID
      INTEGER            NDX1, NDX2
      complex   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      complex   X( * ), B( * ), WORK( LDW,* )
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW .gt. = max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          .gt.   0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*           .ls.   0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N  .ls.  0
*                   -2: LDW  .ls.  N
*                   -3: Maximum number of iterations ITER  .ls. = 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO  .ls.  BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA  .ls.  BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in func GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*
*     .. Parameters ..
      real             ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, P, PHAT, V, S, SHAT, T, MAXIT,
     $                   NEED1, NEED2
      real   TOL,
     $     RHOTOL, OMEGATOL, sGETBREAK, 
     $     scNRM2
      complex   ALPHA, BETA, RHO, RHO1, OMEGA, TMPVAL,
     $     wcdotc
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Funcs ..
      EXTERNAL           sGETBREAK, cAXPY, cCOPY, 
     $      wcdotc, scNRM2, cSCAL
*     ..
*     .. Intrinsic Funcs ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      P    = 3
      V    = 4
      T    = 5
      PHAT = 6
      SHAT = 7
      S    = 1
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((T - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((T - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set parameter tolerances.
*
      RHOTOL = sGETBREAK()
      OMEGATOL = sGETBREAK()
*
*     Set initial residual.
*
      CALL cCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( scNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL cMATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using RTLD[] as temp. storage.
*********CALL cCOPY(N, X, 1, WORK(1,RTLD), 1)
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
      ENDIF
*
*****************
 2    CONTINUE
*****************
*
      IF ( scNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30

      CALL cCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform BiConjugate Gradient Stabilized iteration.
*
      ITER = ITER + 1
*     
      RHO = wcdotc( N, WORK(1,RTLD), 1, WORK(1,R), 1 )
      IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*     
*        Compute vector P.
*
      IF ( ITER.GT.1 ) THEN
         BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
         CALL cAXPY( N, -OMEGA, WORK(1,V), 1, WORK(1,P), 1 )
         CALL cSCAL( N, BETA, WORK(1,P), 1 )
         TMPVAL = ONE
         CALL cAXPY( N, TMPVAL, WORK(1,R), 1, WORK(1,P), 1 )
      ELSE
         CALL cCOPY( N, WORK(1,R), 1, WORK(1,P), 1 )
      ENDIF
*
*        Compute direction adjusting vector PHAT and scalar ALPHA.
*
*********CALL PSOLVE( WORK(1,PHAT), WORK(1,P) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((P    - 1) * LDW) + 1
*     Prepare for return & return
      RLBL = 3
      IJOB = 2
      RETURN
*
*****************
 3    CONTINUE
*****************
*
*********CALL MATVEC( ONE, WORK(1,PHAT), ZERO, WORK(1,V) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((V    - 1) * LDW) + 1
*        Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 4
      IJOB = 1
      RETURN
*
*****************
 4    CONTINUE
*****************
*
      TMPVAL = wcdotc( N, WORK(1,RTLD), 1, WORK(1,V), 1 )
      IF (TMPVAL.EQ.0) THEN
*        Breakdown
         INFO = -11
         GO TO 20
      ENDIF
      ALPHA = RHO / TMPVAL
*
*        Early check for tolerance.
*
      CALL cAXPY( N, -ALPHA, WORK(1,V), 1, WORK(1,R), 1 )
      CALL cCOPY( N, WORK(1,R), 1, WORK(1,S), 1 )
      IF ( scNRM2( N, WORK(1,S), 1 ).LE.TOL ) THEN
         CALL cAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
         RESID = scNRM2( N, WORK(1,S), 1 )
         GO TO 30
      ELSE
*
*           Compute stabilizer vector SHAT and scalar OMEGA.
*
************CALL PSOLVE( WORK(1,SHAT), WORK(1,S) )
*
         NDX1 = ((SHAT - 1) * LDW) + 1
         NDX2 = ((S    - 1) * LDW) + 1
*     Prepare for return & return
         RLBL = 5
         IJOB = 2
         RETURN
      ENDIF
*
*****************
 5    CONTINUE
*****************
*
************CALL MATVEC( ONE, WORK(1,SHAT), ZERO, WORK(1,T) )
*
      NDX1 = ((SHAT - 1) * LDW) + 1
      NDX2 = ((T    - 1) * LDW) + 1
*           Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 6
      IJOB = 1
      RETURN
*
*****************
 6    CONTINUE
*****************
*
      OMEGA = wcdotc( N, WORK(1,T), 1, WORK(1,S), 1 ) / 
     $     wcdotc( N, WORK(1,T), 1, WORK(1,T), 1 )
*
*           Compute new solution approximation vector X.
*
      CALL cAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
      CALL cAXPY( N, OMEGA, WORK(1,SHAT), 1, X, 1 )
*     
*     Compute residual R, check for tolerance.
*
      CALL cAXPY( N, -OMEGA, WORK(1,T), 1, WORK(1,R), 1 )
*
************RESID = DNRM2( N, WORK(1,R), 1 ) / BNRM2 
************IF ( RESID.LE.TOL  ) GO TO 30
*
      NDX1 = NEED1
      NDX2 = NEED2
*     Prepare for resumption & return
      RLBL = 7
      IJOB = 4
      RETURN
*
*****************
 7    CONTINUE
*****************
      IF( INFO.EQ.1 ) GO TO 30
*     
      IF ( ITER.EQ.MAXIT ) THEN
         INFO = 1
         GO TO 20
      ENDIF
*
      IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         GO TO 25
      ELSE
         RHO1 = RHO
         GO TO 10
      ENDIF
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         INFO = -11
      ENDIF
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of BICGSTABREVCOM
*
      END
*     END SUBROUTINE cBICGSTABREVCOM


      SUBROUTINE zBICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, 
     $                    INFO,NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      double precision   RESID
      INTEGER            NDX1, NDX2
      double complex   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      double complex   X( * ), B( * ), WORK( LDW,* )
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW .gt. = max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          .gt.   0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*           .ls.   0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N  .ls.  0
*                   -2: LDW  .ls.  N
*                   -3: Maximum number of iterations ITER  .ls. = 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO  .ls.  BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA  .ls.  BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in func GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*
*     .. Parameters ..
      double precision             ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, P, PHAT, V, S, SHAT, T, MAXIT,
     $                   NEED1, NEED2
      double precision   TOL,
     $     RHOTOL, OMEGATOL, dGETBREAK, 
     $     dzNRM2
      double complex   ALPHA, BETA, RHO, RHO1, OMEGA, TMPVAL,
     $     wzdotc
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Funcs ..
      EXTERNAL           dGETBREAK, zAXPY, zCOPY, 
     $      wzdotc, dzNRM2, zSCAL
*     ..
*     .. Intrinsic Funcs ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      P    = 3
      V    = 4
      T    = 5
      PHAT = 6
      SHAT = 7
      S    = 1
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((T - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((T - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set parameter tolerances.
*
      RHOTOL = dGETBREAK()
      OMEGATOL = dGETBREAK()
*
*     Set initial residual.
*
      CALL zCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( dzNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL zMATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using RTLD[] as temp. storage.
*********CALL zCOPY(N, X, 1, WORK(1,RTLD), 1)
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
      ENDIF
*
*****************
 2    CONTINUE
*****************
*
      IF ( dzNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30

      CALL zCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform BiConjugate Gradient Stabilized iteration.
*
      ITER = ITER + 1
*     
      RHO = wzdotc( N, WORK(1,RTLD), 1, WORK(1,R), 1 )
      IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*     
*        Compute vector P.
*
      IF ( ITER.GT.1 ) THEN
         BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
         CALL zAXPY( N, -OMEGA, WORK(1,V), 1, WORK(1,P), 1 )
         CALL zSCAL( N, BETA, WORK(1,P), 1 )
         TMPVAL = ONE
         CALL zAXPY( N, TMPVAL, WORK(1,R), 1, WORK(1,P), 1 )
      ELSE
         CALL zCOPY( N, WORK(1,R), 1, WORK(1,P), 1 )
      ENDIF
*
*        Compute direction adjusting vector PHAT and scalar ALPHA.
*
*********CALL PSOLVE( WORK(1,PHAT), WORK(1,P) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((P    - 1) * LDW) + 1
*     Prepare for return & return
      RLBL = 3
      IJOB = 2
      RETURN
*
*****************
 3    CONTINUE
*****************
*
*********CALL MATVEC( ONE, WORK(1,PHAT), ZERO, WORK(1,V) )
*
      NDX1 = ((PHAT - 1) * LDW) + 1
      NDX2 = ((V    - 1) * LDW) + 1
*        Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 4
      IJOB = 1
      RETURN
*
*****************
 4    CONTINUE
*****************
*
      TMPVAL = wzdotc( N, WORK(1,RTLD), 1, WORK(1,V), 1 )
      IF (TMPVAL.EQ.0) THEN
*        Breakdown
         INFO = -11
         GO TO 20
      ENDIF
      ALPHA = RHO / TMPVAL
*
*        Early check for tolerance.
*
      CALL zAXPY( N, -ALPHA, WORK(1,V), 1, WORK(1,R), 1 )
      CALL zCOPY( N, WORK(1,R), 1, WORK(1,S), 1 )
      IF ( dzNRM2( N, WORK(1,S), 1 ).LE.TOL ) THEN
         CALL zAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
         RESID = dzNRM2( N, WORK(1,S), 1 )
         GO TO 30
      ELSE
*
*           Compute stabilizer vector SHAT and scalar OMEGA.
*
************CALL PSOLVE( WORK(1,SHAT), WORK(1,S) )
*
         NDX1 = ((SHAT - 1) * LDW) + 1
         NDX2 = ((S    - 1) * LDW) + 1
*     Prepare for return & return
         RLBL = 5
         IJOB = 2
         RETURN
      ENDIF
*
*****************
 5    CONTINUE
*****************
*
************CALL MATVEC( ONE, WORK(1,SHAT), ZERO, WORK(1,T) )
*
      NDX1 = ((SHAT - 1) * LDW) + 1
      NDX2 = ((T    - 1) * LDW) + 1
*           Prepare for return & return
      SCLR1 = ONE
      SCLR2 = ZERO
      RLBL = 6
      IJOB = 1
      RETURN
*
*****************
 6    CONTINUE
*****************
*
      OMEGA = wzdotc( N, WORK(1,T), 1, WORK(1,S), 1 ) / 
     $     wzdotc( N, WORK(1,T), 1, WORK(1,T), 1 )
*
*           Compute new solution approximation vector X.
*
      CALL zAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
      CALL zAXPY( N, OMEGA, WORK(1,SHAT), 1, X, 1 )
*     
*     Compute residual R, check for tolerance.
*
      CALL zAXPY( N, -OMEGA, WORK(1,T), 1, WORK(1,R), 1 )
*
************RESID = DNRM2( N, WORK(1,R), 1 ) / BNRM2 
************IF ( RESID.LE.TOL  ) GO TO 30
*
      NDX1 = NEED1
      NDX2 = NEED2
*     Prepare for resumption & return
      RLBL = 7
      IJOB = 4
      RETURN
*
*****************
 7    CONTINUE
*****************
      IF( INFO.EQ.1 ) GO TO 30
*     
      IF ( ITER.EQ.MAXIT ) THEN
         INFO = 1
         GO TO 20
      ENDIF
*
      IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         GO TO 25
      ELSE
         RHO1 = RHO
         GO TO 10
      ENDIF
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         INFO = -11
      ENDIF
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of BICGSTABREVCOM
*
      END
*     END SUBROUTINE zBICGSTABREVCOM


