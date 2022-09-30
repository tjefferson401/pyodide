*> \brief \b CGEMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CGEMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgemqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgemqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgemqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
*                          C, LDC, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*       ..
*       .. Array Arguments ..
*       COMPLEX   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGEMQRT overwrites the general complex M-by-N matrix C with
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q C            C Q
*> TRANS = 'C':    Q**H C            C Q**H
*>
*> where Q is a complex orthogonal matrix defined as the product of K
*> elementary reflectors:
*>
*>       Q = H(1) H(2) . . . H(K) = I - V T V**H
*>
*> generated using the compact WY representation as returned by CGEQRT.
*>
*> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left;
*>          = 'R': apply Q or Q**H from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'C':  Conjugate transpose, apply Q**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CGEQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CGEQRT in the first K columns of its array argument A.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CGEQRT, stored as a NB-by-N matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q C, Q**H C, C Q**H or C Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array. The dimension of WORK is
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexGEcomputational
*
*  =====================================================================
      SUBROUTINE CGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
     $                   C, LDC, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*     ..
*     .. Array Arguments ..
      COMPLEX   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, LDWORK, KF, Q
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CLARFB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'C' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF( LEFT ) THEN
         LDWORK = MAX( 1, N )
         Q = M
      ELSE IF ( RIGHT ) THEN
         LDWORK = MAX( 1, M )
         Q = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.Q ) THEN
         INFO = -5
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0)) THEN
         INFO = -6
      ELSE IF( LDV.LT.MAX( 1, Q ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL CLARFB( 'L', 'C', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL CLARFB( 'R', 'N', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL CLARFB( 'L', 'N', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL CLARFB( 'R', 'C', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      END IF
*
      RETURN
*
*     End of CGEMQRT
*
      END
*> \brief \b CGEQRFP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CGEQRFP + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrfp.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrfp.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrfp.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGEQR2P computes a QR factorization of a complex M-by-N matrix A:
*>
*>    A = Q * ( R ),
*>            ( 0 )
*>
*> where:
*>
*>    Q is a M-by-M orthogonal matrix;
*>    R is an upper-triangular N-by-N matrix with nonnegative diagonal
*>    entries;
*>    0 is a (M-N)-by-N zero matrix, if M > N.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if m >= n). The diagonal entries of R
*>          are real and nonnegative; the elements below the diagonal,
*>          with the array TAU, represent the unitary matrix Q as a
*>          product of min(m,n) elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**H
*>
*>  where tau is a complex scalar, and v is a complex vector with
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*>  and tau in TAU(i).
*>
*> See Lapack Working Note 203 for details
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEQR2P, CLARFB, CLARFT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQRFP', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'CGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'CGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL CGEQR2P( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL CLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H**H to A(i:m,i+ib:n) from the left
*
               CALL CLARFB( 'Left', 'Conjugate transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL CGEQR2P( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of CGEQRFP
*
      END
*> \brief \b CGEQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CGEQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDT, M, N, NB
*       ..
*       .. Array Arguments ..
*       COMPLEX A( LDA, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGEQRT computes a blocked QR factorization of a complex M-by-N matrix A
*> using the compact WY representation of Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if M >= N); the elements below the diagonal
*>          are the columns of V.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,MIN(M,N))
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See below
*>          for further details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix V stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal. For example, if M=5 and N=3, the matrix V is
*>
*>               V = (  1       )
*>                   ( v1  1    )
*>                   ( v1 v2  1 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  where the vi's represent the vectors which define H(i), which are returned
*>  in the matrix A.  The 1's along the diagonal of V are not stored in A.
*>
*>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-K matrix T as
*>
*>               T = (T1 T2 ... TB).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDT, M, N, NB
*     ..
*     .. Array Arguments ..
      COMPLEX A( LDA, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, IINFO, K
      LOGICAL    USE_RECURSIVE_QR
      PARAMETER( USE_RECURSIVE_QR=.TRUE. )
*     ..
*     .. External Subroutines ..
      EXTERNAL   CGEQRT2, CGEQRT3, CLARFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NB.LT.1 .OR. ( NB.GT.MIN(M,N) .AND. MIN(M,N).GT.0 ) )THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN
*
*     Blocked loop of length K
*
      DO I = 1, K,  NB
         IB = MIN( K-I+1, NB )
*
*     Compute the QR factorization of the current block A(I:M,I:I+IB-1)
*
         IF( USE_RECURSIVE_QR ) THEN
            CALL CGEQRT3( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         ELSE
            CALL CGEQRT2( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         END IF
         IF( I+IB.LE.N ) THEN
*
*     Update by applying H**H to A(I:M,I+IB:N) from the left
*
            CALL CLARFB( 'L', 'C', 'F', 'C', M-I+1, N-I-IB+1, IB,
     $                   A( I, I ), LDA, T( 1, I ), LDT,
     $                   A( I, I+IB ), LDA, WORK , N-I-IB+1 )
         END IF
      END DO
      RETURN
*
*     End of CGEQRT
*
      END
*> \brief \b CLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CLAHQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
*                          IHIZ, Z, LDZ, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    CLAHQR is an auxiliary routine called by CHSEQR to update the
*>    eigenvalues and Schur decomposition already computed by CHSEQR, by
*>    dealing with the Hessenberg submatrix in rows and columns ILO to
*>    IHI.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          = .TRUE. : the full Schur form T is required;
*>          = .FALSE.: only eigenvalues are required.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE. : the matrix of Schur vectors Z is required;
*>          = .FALSE.: Schur vectors are not required.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix H.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          It is assumed that H is already upper triangular in rows and
*>          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
*>          CLAHQR works primarily with the Hessenberg submatrix in rows
*>          and columns ILO to IHI, but applies transformations to all of
*>          H if WANTT is .TRUE..
*>          1 <= ILO <= max(1,IHI); IHI <= N.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX array, dimension (LDH,N)
*>          On entry, the upper Hessenberg matrix H.
*>          On exit, if INFO is zero and if WANTT is .TRUE., then H
*>          is upper triangular in rows and columns ILO:IHI.  If INFO
*>          is zero and if WANTT is .FALSE., then the contents of H
*>          are unspecified on exit.  The output state of H in case
*>          INF is positive is below under the description of INFO.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          The leading dimension of the array H. LDH >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX array, dimension (N)
*>          The computed eigenvalues ILO to IHI are stored in the
*>          corresponding elements of W. If WANTT is .TRUE., the
*>          eigenvalues are stored in the same order as on the diagonal
*>          of the Schur form returned in H, with W(i) = H(i,i).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>          Specify the rows of Z to which transformations must be
*>          applied if WANTZ is .TRUE..
*>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX array, dimension (LDZ,N)
*>          If WANTZ is .TRUE., on entry Z must contain the current
*>          matrix Z of transformations accumulated by CHSEQR, and on
*>          exit Z has been updated; transformations are applied only to
*>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*>          If WANTZ is .FALSE., Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           = 0:  successful exit
*>           > 0:  if INFO = i, CLAHQR failed to compute all the
*>                  eigenvalues ILO to IHI in a total of 30 iterations
*>                  per eigenvalue; elements i+1:ihi of W contain
*>                  those eigenvalues which have been successfully
*>                  computed.
*>
*>                  If INFO > 0 and WANTT is .FALSE., then on exit,
*>                  the remaining unconverged eigenvalues are the
*>                  eigenvalues of the upper Hessenberg matrix
*>                  rows and columns ILO through INFO of the final,
*>                  output value of H.
*>
*>                  If INFO > 0 and WANTT is .TRUE., then on exit
*>          (*)       (initial value of H)*U  = U*(final value of H)
*>                  where U is an orthogonal matrix.    The final
*>                  value of H is upper Hessenberg and triangular in
*>                  rows and columns INFO+1 through IHI.
*>
*>                  If INFO > 0 and WANTZ is .TRUE., then on exit
*>                      (final value of Z)  = (initial value of Z)*U
*>                  where U is the orthogonal matrix in (*)
*>                  (regardless of the value of WANTT.)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexOTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>     02-96 Based on modifications by
*>     David Day, Sandia National Laboratory, USA
*>
*>     12-04 Further modifications by
*>     Ralph Byers, University of Kansas, USA
*>     This is a modified version of CLAHQR from LAPACK version 3.0.
*>     It is (1) more robust against overflow and underflow and
*>     (2) adopts the more conservative Ahues & Tisseur stopping
*>     criterion (LAWN 122, 1997).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
*     ..
*
*  =========================================================
*
*     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
     $                   ONE = ( 1.0e0, 0.0e0 ) )
      REAL               RZERO, RONE, HALF
      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0, HALF = 0.5e0 )
      REAL               DAT1
      PARAMETER          ( DAT1 = 3.0e0 / 4.0e0 )
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 10 )
*     ..
*     .. Local Scalars ..
      COMPLEX            CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U,
     $                   V2, X, Y
      REAL               AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX,
     $                   SAFMIN, SMLNUM, SX, T2, TST, ULP
      INTEGER            I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M,
     $                   NH, NZ, KDEFL
*     ..
*     .. Local Arrays ..
      COMPLEX            V( 2 )
*     ..
*     .. External Functions ..
      COMPLEX            CLADIV
      REAL               SLAMCH
      EXTERNAL           CLADIV, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CLARFG, CSCAL, SLABAD
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CONJG, MAX, MIN, REAL, SQRT
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         W( ILO ) = H( ILO, ILO )
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 )
     $   H( IHI, IHI-2 ) = ZERO
*     ==== ensure that subdiagonal entries are real ====
      IF( WANTT ) THEN
         JLO = 1
         JHI = N
      ELSE
         JLO = ILO
         JHI = IHI
      END IF
      DO 20 I = ILO + 1, IHI
         IF( AIMAG( H( I, I-1 ) ).NE.RZERO ) THEN
*           ==== The following redundant normalization
*           .    avoids problems with both gradual and
*           .    sudden underflow in ABS(H(I,I-1)) ====
            SC = H( I, I-1 ) / CABS1( H( I, I-1 ) )
            SC = CONJG( SC ) / ABS( SC )
            H( I, I-1 ) = ABS( H( I, I-1 ) )
            CALL CSCAL( JHI-I+1, SC, H( I, I ), LDH )
            CALL CSCAL( MIN( JHI, I+1 )-JLO+1, CONJG( SC ), H( JLO, I ),
     $                  1 )
            IF( WANTZ )
     $         CALL CSCAL( IHIZ-ILOZ+1, CONJG( SC ), Z( ILOZ, I ), 1 )
         END IF
   20 CONTINUE
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL SLABAD( SAFMIN, SAFMAX )
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( NH ) / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITMAX is the total number of QR iterations allowed.
*
      ITMAX = 30 * MAX( 10, NH )
*
*     KDEFL counts the number of iterations since a deflation
*
      KDEFL = 0
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   30 CONTINUE
      IF( I.LT.ILO )
     $   GO TO 150
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      L = ILO
      DO 130 ITS = 0, ITMAX
*
*        Look for a single small subdiagonal element.
*
         DO 40 K = I, L + 1, -1
            IF( CABS1( H( K, K-1 ) ).LE.SMLNUM )
     $         GO TO 50
            TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO )
     $            TST = TST + ABS( REAL( H( K-1, K-2 ) ) )
               IF( K+1.LE.IHI )
     $            TST = TST + ABS( REAL( H( K+1, K ) ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some examples.  ====
            IF( ABS( REAL( H( K, K-1 ) ) ).LE.ULP*TST ) THEN
               AB = MAX( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               BA = MIN( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               AA = MAX( CABS1( H( K, K ) ),
     $              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( CABS1( H( K, K ) ),
     $              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM,
     $             ULP*( BB*( AA / S ) ) ) )GO TO 50
            END IF
   40    CONTINUE
   50    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 has split off.
*
         IF( L.GE.I )
     $      GO TO 140
         KDEFL = KDEFL + 1
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( MOD(KDEFL,2*KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = DAT1*ABS( REAL( H( I, I-1 ) ) )
            T = S + H( I, I )
         ELSE IF( MOD(KDEFL,KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = DAT1*ABS( REAL( H( L+1, L ) ) )
            T = S + H( L, L )
         ELSE
*
*           Wilkinson's shift.
*
            T = H( I, I )
            U = SQRT( H( I-1, I ) )*SQRT( H( I, I-1 ) )
            S = CABS1( U )
            IF( S.NE.RZERO ) THEN
               X = HALF*( H( I-1, I-1 )-T )
               SX = CABS1( X )
               S = MAX( S, CABS1( X ) )
               Y = S*SQRT( ( X / S )**2+( U / S )**2 )
               IF( SX.GT.RZERO ) THEN
                  IF( REAL( X / SX )*REAL( Y )+AIMAG( X / SX )*
     $                AIMAG( Y ).LT.RZERO )Y = -Y
               END IF
               T = T - U*CLADIV( U, ( X+Y ) )
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 60 M = I - 1, L + 1, -1
*
*           Determine the effect of starting the single-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.
*
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H11S = H11 - T
            H21 = REAL( H( M+1, M ) )
            S = CABS1( H11S ) + ABS( H21 )
            H11S = H11S / S
            H21 = H21 / S
            V( 1 ) = H11S
            V( 2 ) = H21
            H10 = REAL( H( M, M-1 ) )
            IF( ABS( H10 )*ABS( H21 ).LE.ULP*
     $          ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) )
     $          GO TO 70
   60    CONTINUE
         H11 = H( L, L )
         H22 = H( L+1, L+1 )
         H11S = H11 - T
         H21 = REAL( H( L+1, L ) )
         S = CABS1( H11S ) + ABS( H21 )
         H11S = H11S / S
         H21 = H21 / S
         V( 1 ) = H11S
         V( 2 ) = H21
   70    CONTINUE
*
*        Single-shift QR step
*
         DO 120 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix.
*
*           V(2) is always real before the call to CLARFG, and hence
*           after the call T2 ( = T1*V(2) ) is also real.
*
            IF( K.GT.M )
     $         CALL CCOPY( 2, H( K, K-1 ), 1, V, 1 )
            CALL CLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
            END IF
            V2 = V( 2 )
            T2 = REAL( T1*V2 )
*
*           Apply G from the left to transform the rows of the matrix
*           in columns K to I2.
*
            DO 80 J = K, I2
               SUM = CONJG( T1 )*H( K, J ) + T2*H( K+1, J )
               H( K, J ) = H( K, J ) - SUM
               H( K+1, J ) = H( K+1, J ) - SUM*V2
   80       CONTINUE
*
*           Apply G from the right to transform the columns of the
*           matrix in rows I1 to min(K+2,I).
*
            DO 90 J = I1, MIN( K+2, I )
               SUM = T1*H( J, K ) + T2*H( J, K+1 )
               H( J, K ) = H( J, K ) - SUM
               H( J, K+1 ) = H( J, K+1 ) - SUM*CONJG( V2 )
   90       CONTINUE
*
            IF( WANTZ ) THEN
*
*              Accumulate transformations in the matrix Z
*
               DO 100 J = ILOZ, IHIZ
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                  Z( J, K ) = Z( J, K ) - SUM
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*CONJG( V2 )
  100          CONTINUE
            END IF
*
            IF( K.EQ.M .AND. M.GT.L ) THEN
*
*              If the QR step was started at row M > L because two
*              consecutive small subdiagonals were found, then extra
*              scaling must be performed to ensure that H(M,M-1) remains
*              real.
*
               TEMP = ONE - T1
               TEMP = TEMP / ABS( TEMP )
               H( M+1, M ) = H( M+1, M )*CONJG( TEMP )
               IF( M+2.LE.I )
     $            H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
               DO 110 J = M, I
                  IF( J.NE.M+1 ) THEN
                     IF( I2.GT.J )
     $                  CALL CSCAL( I2-J, TEMP, H( J, J+1 ), LDH )
                     CALL CSCAL( J-I1, CONJG( TEMP ), H( I1, J ), 1 )
                     IF( WANTZ ) THEN
                        CALL CSCAL( NZ, CONJG( TEMP ), Z( ILOZ, J ), 1 )
                     END IF
                  END IF
  110          CONTINUE
            END IF
  120    CONTINUE
*
*        Ensure that H(I,I-1) is real.
*
         TEMP = H( I, I-1 )
         IF( AIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = ABS( TEMP )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I )
     $         CALL CSCAL( I2-I, CONJG( TEMP ), H( I, I+1 ), LDH )
            CALL CSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( WANTZ ) THEN
               CALL CSCAL( NZ, TEMP, Z( ILOZ, I ), 1 )
            END IF
         END IF
*
  130 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  140 CONTINUE
*
*     H(I,I-1) is negligible: one eigenvalue has converged.
*
      W( I ) = H( I, I )
*     reset deflation counter
      KDEFL = 0
*
*     return to start of the main loop with new value of I.
*
      I = L - 1
      GO TO 30
*
  150 CONTINUE
      RETURN
*
*     End of CLAHQR
*
      END
*> \brief \b CSYCONV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CSYCONV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX            A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CSYCONV convert A given by TRF into L and D and vice-versa.
*> Get Non-diag elements of D (returned in workspace) and
*> apply or reverse permutation done in TRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by CSYTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by CSYTRF.
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is COMPLEX array, dimension (N)
*>          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1
*>          or 2-by-2 block diagonal matrix D in LDLT.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexSYcomputational
*
*  =====================================================================
      SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = (0.0E+0,0.0E+0) )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, J
      COMPLEX            TEMP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSYCONV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*      A is UPPER
*
*      Convert A (A is upper)
*
*        Convert VALUE
*
         IF ( CONVERT ) THEN
            I=N
            E(1)=ZERO
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I-1,I)
                  E(I-1)=ZERO
                  A(I-1,I)=ZERO
                  I=I-1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I-1
            END DO
*
*        Convert PERMUTATIONS
*
         I=N
         DO WHILE ( I .GE. 1 )
            IF( IPIV(I) .GT. 0) THEN
               IP=IPIV(I)
               IF( I .LT. N) THEN
                  DO 12 J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
 12            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
               IF( I .LT. N) THEN
             DO 13 J= I+1,N
                 TEMP=A(IP,J)
                 A(IP,J)=A(I-1,J)
                 A(I-1,J)=TEMP
 13            CONTINUE
                ENDIF
                I=I-1
           ENDIF
           I=I-1
        END DO

         ELSE
*
*      Revert A (A is upper)
*
*
*        Revert PERMUTATIONS
*
            I=1
            DO WHILE ( I .LE. N )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF( I .LT. N) THEN
                  DO J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
                  END DO
                  ENDIF
               ELSE
                 IP=-IPIV(I)
                 I=I+1
                 IF( I .LT. N) THEN
                    DO J= I+1,N
                       TEMP=A(IP,J)
                       A(IP,J)=A(I-1,J)
                       A(I-1,J)=TEMP
                    END DO
                 ENDIF
               ENDIF
               I=I+1
            END DO
*
*        Revert VALUE
*
            I=N
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I-1,I)=E(I)
                  I=I-1
               ENDIF
               I=I-1
            END DO
         END IF
      ELSE
*
*      A is LOWER
*
         IF ( CONVERT ) THEN
*
*      Convert A (A is lower)
*
*
*        Convert VALUE
*
            I=1
            E(N)=ZERO
            DO WHILE ( I .LE. N )
               IF( I.LT.N .AND. IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I+1,I)
                  E(I+1)=ZERO
                  A(I+1,I)=ZERO
                  I=I+1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I+1
            END DO
*
*        Convert PERMUTATIONS
*
         I=1
         DO WHILE ( I .LE. N )
            IF( IPIV(I) .GT. 0 ) THEN
               IP=IPIV(I)
               IF (I .GT. 1) THEN
               DO 22 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I,J)
                 A(I,J)=TEMP
 22            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
              IF (I .GT. 1) THEN
              DO 23 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I+1,J)
                 A(I+1,J)=TEMP
 23           CONTINUE
              ENDIF
              I=I+1
           ENDIF
           I=I+1
        END DO
         ELSE
*
*      Revert A (A is lower)
*
*
*        Revert PERMUTATIONS
*
            I=N
            DO WHILE ( I .GE. 1 )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I,J)
                        A(I,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ELSE
                  IP=-IPIV(I)
                  I=I-1
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I+1,J)
                        A(I+1,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ENDIF
               I=I-1
            END DO
*
*        Revert VALUE
*
            I=1
            DO WHILE ( I .LE. N-1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I+1,I)=E(I)
                  I=I+1
               ENDIF
               I=I+1
            END DO
         END IF
      END IF

      RETURN
*
*     End of CSYCONV
*
      END
*> \brief \b CSYCONVF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CSYCONVF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconvf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconvf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconvf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX            A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> CSYCONVF converts the factorization output format used in
*> CSYTRF provided on entry in parameter A into the factorization
*> output format used in CSYTRF_RK (or CSYTRF_BK) that is stored
*> on exit in parameters A and E. It also converts in place details of
*> the intechanges stored in IPIV from the format used in CSYTRF into
*> the format used in CSYTRF_RK (or CSYTRF_BK).
*>
*> If parameter WAY = 'R':
*> CSYCONVF performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in CSYTRF_RK
*> (or CSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in CSYTRF that is stored
*> on exit in parameter A. It also converts in place details of
*> the intechanges stored in IPIV from the format used in CSYTRF_RK
*> (or CSYTRF_BK) into the format used in CSYTRF.
*>
*> CSYCONVF can also convert in Hermitian matrix case, i.e. between
*> formats used in CHETRF and CHETRF_RK (or CHETRF_BK).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          CSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          CSYTRF_RK or CSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          CSYTRF_RK or CSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          CSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is COMPLEX array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in,out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>
*>          1) If WAY ='C':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in CSYTRF.
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in CSYTRF_RK
*>          ( or CSYTRF_BK).
*>
*>          1) If WAY ='R':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in CSYTRF_RK
*>          ( or CSYTRF_BK).
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in CSYTRF.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexSYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE CSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           CSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSYCONVF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL CSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL CSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i-1 and IPIV(i-1),
*                 so this should be recorded in two consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I-1 )
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where k increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL CSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL CSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i+1 and IPIV(i+1),
*                 so this should be recorded in consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I+1 )
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of CSYCONVF
*
      END
*> \brief \b CSYCONVF_ROOK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CSYCONVF_ROOK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconvf_rook.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconvf_rook.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconvf_rook.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX            A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> CSYCONVF_ROOK converts the factorization output format used in
*> CSYTRF_ROOK provided on entry in parameter A into the factorization
*> output format used in CSYTRF_RK (or CSYTRF_BK) that is stored
*> on exit in parameters A and E. IPIV format for CSYTRF_ROOK and
*> CSYTRF_RK (or CSYTRF_BK) is the same and is not converted.
*>
*> If parameter WAY = 'R':
*> CSYCONVF_ROOK performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in CSYTRF_RK
*> (or CSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in CSYTRF_ROOK that is stored
*> on exit in parameter A. IPIV format for CSYTRF_ROOK and
*> CSYTRF_RK (or CSYTRF_BK) is the same and is not converted.
*>
*> CSYCONVF_ROOK can also convert in Hermitian matrix case, i.e. between
*> formats used in CHETRF_ROOK and CHETRF_RK (or CHETRF_BK).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          CSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          CSYTRF_RK or CSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          CSYTRF_RK or CSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          CSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is COMPLEX array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          On entry, details of the interchanges and the block
*>          structure of D as determined:
*>          1) by CSYTRF_ROOK, if WAY ='C';
*>          2) by CSYTRF_RK (or CSYTRF_BK), if WAY ='R'.
*>          The IPIV format is the same for all these routines.
*>
*>          On exit, is not changed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexSYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE CSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           CSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, IP2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSYCONVF_ROOK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1)
*                 in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                     IF( IP2.NE.(I-1) ) THEN
                        CALL CSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP2, I+1 ), LDA )
                     END IF
                  END IF
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i)
*                 in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP2.NE.(I-1) ) THEN
                        CALL CSWAP( N-I, A( IP2, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1)
*                 in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                     IF( IP2.NE.(I+1) ) THEN
                        CALL CSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP2, 1 ), LDA )
                     END IF
                  END IF
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i)
*                 in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP2.NE.(I+1) ) THEN
                        CALL CSWAP( I-1, A( IP2, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL CSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of CSYCONVF_ROOK
*
      END
*> \brief \b CTPMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CTPMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpmqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpmqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpmqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
*                           A, LDA, B, LDB, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*       ..
*       .. Array Arguments ..
*       COMPLEX   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
*      $          WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CTPMQRT applies a complex orthogonal matrix Q obtained from a
*> "triangular-pentagonal" complex block reflector H to a general
*> complex matrix C, which consists of two blocks A and B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left;
*>          = 'R': apply Q or Q**H from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'C':  Conjugate transpose, apply Q**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The order of the trapezoidal part of V.
*>          K >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CTPQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CTPQRT in B.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDV >= max(1,M);
*>          if SIDE = 'R', LDV >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CTPQRT, stored as a NB-by-K matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension
*>          (LDA,N) if SIDE = 'L' or
*>          (LDA,K) if SIDE = 'R'
*>          On entry, the K-by-N or M-by-K matrix A.
*>          On exit, A is overwritten by the corresponding block of
*>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDC >= max(1,K);
*>          If SIDE = 'R', LDC >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,N)
*>          On entry, the M-by-N matrix B.
*>          On exit, B is overwritten by the corresponding block of
*>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.
*>          LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array. The dimension of WORK is
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The columns of the pentagonal matrix V contain the elementary reflectors
*>  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
*>  trapezoidal block V2:
*>
*>        V = [V1]
*>            [V2].
*>
*>  The size of the trapezoidal block V2 is determined by the parameter L,
*>  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L
*>  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular;
*>  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
*>
*>  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K.
*>                      [B]
*>
*>  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K.
*>
*>  The complex orthogonal matrix Q is formed from V and T.
*>
*>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.
*>
*>  If TRANS='C' and SIDE='L', C is on exit replaced with Q**H * C.
*>
*>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.
*>
*>  If TRANS='C' and SIDE='R', C is on exit replaced with C * Q**H.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
     $                    A, LDA, B, LDB, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*     ..
*     .. Array Arguments ..
      COMPLEX   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
     $          WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, MB, LB, KF, LDAQ, LDVQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTPRFB, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'C' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF ( LEFT ) THEN
         LDVQ = MAX( 1, M )
         LDAQ = MAX( 1, K )
      ELSE IF ( RIGHT ) THEN
         LDVQ = MAX( 1, N )
         LDAQ = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      ELSE IF( L.LT.0 .OR. L.GT.K ) THEN
         INFO = -6
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0) ) THEN
         INFO = -7
      ELSE IF( LDV.LT.LDVQ ) THEN
         INFO = -9
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -11
      ELSE IF( LDA.LT.LDAQ ) THEN
         INFO = -13
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTPMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL CTPRFB( 'L', 'C', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL CTPRFB( 'R', 'N', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL CTPRFB( 'L', 'N', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL CTPRFB( 'R', 'C', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      END IF
*
      RETURN
*
*     End of CTPMQRT
*
      END
*> \brief \b CTPQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CTPQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*       ..
*       .. Array Arguments ..
*       COMPLEX A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CTPQRT computes a blocked QR factorization of a complex
*> "triangular-pentagonal" matrix C, which is composed of a
*> triangular block A and pentagonal block B, using the compact
*> WY representation for Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B.
*>          M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B, and the order of the
*>          triangular matrix A.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The number of rows of the upper trapezoidal part of B.
*>          MIN(M,N) >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  N >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the upper triangular N-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the upper triangular matrix R.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,N)
*>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
*>          are rectangular, and the last L rows are upper trapezoidal.
*>          On exit, B contains the pentagonal matrix V.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,N)
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complexOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The input matrix C is a (N+M)-by-N matrix
*>
*>               C = [ A ]
*>                   [ B ]
*>
*>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
*>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
*>  upper trapezoidal matrix B2:
*>
*>               B = [ B1 ]  <- (M-L)-by-N rectangular
*>                   [ B2 ]  <-     L-by-N upper trapezoidal.
*>
*>  The upper trapezoidal matrix B2 consists of the first L rows of a
*>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
*>  B is rectangular M-by-N; if M=L=N, B is upper triangular.
*>
*>  The matrix W stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal (of A) in the (N+M)-by-N input matrix C
*>
*>               C = [ A ]  <- upper triangular N-by-N
*>                   [ B ]  <- M-by-N pentagonal
*>
*>  so that W can be represented as
*>
*>               W = [ I ]  <- identity, N-by-N
*>                   [ V ]  <- M-by-N, same form as B.
*>
*>  Thus, all of information needed for W is contained on exit in B, which
*>  we call V above.  Note that V has the same form as B; that is,
*>
*>               V = [ V1 ] <- (M-L)-by-N rectangular
*>                   [ V2 ] <-     L-by-N upper trapezoidal.
*>
*>  The columns of V represent the vectors which define the H(i)'s.
*>
*>  The number of blocks is B = ceiling(N/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-N matrix T as
*>
*>               T = [T1 T2 ... TB].
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*     ..
*     .. Array Arguments ..
      COMPLEX A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, LB, MB, IINFO
*     ..
*     .. External Subroutines ..
      EXTERNAL   CTPQRT2, CTPRFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 .OR. (L.GT.MIN(M,N) .AND. MIN(M,N).GE.0)) THEN
         INFO = -3
      ELSE IF( NB.LT.1 .OR. (NB.GT.N .AND. N.GT.0)) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTPQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      DO I = 1, N, NB
*
*     Compute the QR factorization of the current block
*
         IB = MIN( N-I+1, NB )
         MB = MIN( M-L+I+IB-1, M )
         IF( I.GE.L ) THEN
            LB = 0
         ELSE
            LB = MB-M+L-I+1
         END IF
*
         CALL CTPQRT2( MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB,
     $                 T(1, I ), LDT, IINFO )
*
*     Update by applying H**H to B(:,I+IB:N) from the left
*
         IF( I+IB.LE.N ) THEN
            CALL CTPRFB( 'L', 'C', 'F', 'C', MB, N-I-IB+1, IB, LB,
     $                    B( 1, I ), LDB, T( 1, I ), LDT,
     $                    A( I, I+IB ), LDA, B( 1, I+IB ), LDB,
     $                    WORK, IB )
         END IF
      END DO
      RETURN
*
*     End of CTPQRT
*
      END
*> \brief \b DGEMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGEMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgemqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgemqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgemqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
*                          C, LDC, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMQRT overwrites the general real M-by-N matrix C with
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q C            C Q
*> TRANS = 'T':   Q**T C            C Q**T
*>
*> where Q is a real orthogonal matrix defined as the product of K
*> elementary reflectors:
*>
*>       Q = H(1) H(2) . . . H(K) = I - V T V**T
*>
*> generated using the compact WY representation as returned by DGEQRT.
*>
*> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**T from the Left;
*>          = 'R': apply Q or Q**T from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'C':  Transpose, apply Q**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CGEQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CGEQRT in the first K columns of its array argument A.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CGEQRT, stored as a NB-by-N matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array. The dimension of
*>          WORK is N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      SUBROUTINE DGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
     $                   C, LDC, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, LDWORK, KF, Q
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DLARFB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'T' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF( LEFT ) THEN
         LDWORK = MAX( 1, N )
         Q = M
      ELSE IF ( RIGHT ) THEN
         LDWORK = MAX( 1, M )
         Q = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.Q ) THEN
         INFO = -5
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0)) THEN
         INFO = -6
      ELSE IF( LDV.LT.MAX( 1, Q ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL DLARFB( 'L', 'T', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL DLARFB( 'R', 'N', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL DLARFB( 'L', 'N', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL DLARFB( 'R', 'T', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      END IF
*
      RETURN
*
*     End of DGEMQRT
*
      END
*> \brief \b DGEQRFP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGEQRFP + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrfp.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrfp.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrfp.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEQR2P computes a QR factorization of a real M-by-N matrix A:
*>
*>    A = Q * ( R ),
*>            ( 0 )
*>
*> where:
*>
*>    Q is a M-by-M orthogonal matrix;
*>    R is an upper-triangular N-by-N matrix with nonnegative diagonal
*>    entries;
*>    0 is a (M-N)-by-N zero matrix, if M > N.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if m >= n). The diagonal entries of R
*>          are nonnegative; the elements below the diagonal,
*>          with the array TAU, represent the orthogonal matrix Q as a
*>          product of min(m,n) elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*>  and tau in TAU(i).
*>
*> See Lapack Working Note 203 for details
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQR2P, DLARFB, DLARFT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRFP', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL DGEQR2P( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H**T to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL DGEQR2P( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGEQRFP
*
      END
*> \brief \b DGEQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGEQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDT, M, N, NB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A( LDA, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEQRT computes a blocked QR factorization of a real M-by-N matrix A
*> using the compact WY representation of Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if M >= N); the elements below the diagonal
*>          are the columns of V.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,MIN(M,N))
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See below
*>          for further details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix V stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal. For example, if M=5 and N=3, the matrix V is
*>
*>               V = (  1       )
*>                   ( v1  1    )
*>                   ( v1 v2  1 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  where the vi's represent the vectors which define H(i), which are returned
*>  in the matrix A.  The 1's along the diagonal of V are not stored in A.
*>
*>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-K matrix T as
*>
*>               T = (T1 T2 ... TB).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDT, M, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A( LDA, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, IINFO, K
      LOGICAL    USE_RECURSIVE_QR
      PARAMETER( USE_RECURSIVE_QR=.TRUE. )
*     ..
*     .. External Subroutines ..
      EXTERNAL   DGEQRT2, DGEQRT3, DLARFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NB.LT.1 .OR. ( NB.GT.MIN(M,N) .AND. MIN(M,N).GT.0 ) )THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN
*
*     Blocked loop of length K
*
      DO I = 1, K,  NB
         IB = MIN( K-I+1, NB )
*
*     Compute the QR factorization of the current block A(I:M,I:I+IB-1)
*
         IF( USE_RECURSIVE_QR ) THEN
            CALL DGEQRT3( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         ELSE
            CALL DGEQRT2( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         END IF
         IF( I+IB.LE.N ) THEN
*
*     Update by applying H**T to A(I:M,I+IB:N) from the left
*
            CALL DLARFB( 'L', 'T', 'F', 'C', M-I+1, N-I-IB+1, IB,
     $                   A( I, I ), LDA, T( 1, I ), LDT,
     $                   A( I, I+IB ), LDA, WORK , N-I-IB+1 )
         END IF
      END DO
      RETURN
*
*     End of DGEQRT
*
      END
*> \brief \b DLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAHQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
*                          ILOZ, IHIZ, Z, LDZ, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DLAHQR is an auxiliary routine called by DHSEQR to update the
*>    eigenvalues and Schur decomposition already computed by DHSEQR, by
*>    dealing with the Hessenberg submatrix in rows and columns ILO to
*>    IHI.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          = .TRUE. : the full Schur form T is required;
*>          = .FALSE.: only eigenvalues are required.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE. : the matrix of Schur vectors Z is required;
*>          = .FALSE.: Schur vectors are not required.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix H.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          It is assumed that H is already upper quasi-triangular in
*>          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
*>          ILO = 1). DLAHQR works primarily with the Hessenberg
*>          submatrix in rows and columns ILO to IHI, but applies
*>          transformations to all of H if WANTT is .TRUE..
*>          1 <= ILO <= max(1,IHI); IHI <= N.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is DOUBLE PRECISION array, dimension (LDH,N)
*>          On entry, the upper Hessenberg matrix H.
*>          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
*>          quasi-triangular in rows and columns ILO:IHI, with any
*>          2-by-2 diagonal blocks in standard form. If INFO is zero
*>          and WANTT is .FALSE., the contents of H are unspecified on
*>          exit.  The output state of H if INFO is nonzero is given
*>          below under the description of INFO.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          The leading dimension of the array H. LDH >= max(1,N).
*> \endverbatim
*>
*> \param[out] WR
*> \verbatim
*>          WR is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] WI
*> \verbatim
*>          WI is DOUBLE PRECISION array, dimension (N)
*>          The real and imaginary parts, respectively, of the computed
*>          eigenvalues ILO to IHI are stored in the corresponding
*>          elements of WR and WI. If two eigenvalues are computed as a
*>          complex conjugate pair, they are stored in consecutive
*>          elements of WR and WI, say the i-th and (i+1)th, with
*>          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
*>          eigenvalues are stored in the same order as on the diagonal
*>          of the Schur form returned in H, with WR(i) = H(i,i), and, if
*>          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
*>          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>          Specify the rows of Z to which transformations must be
*>          applied if WANTZ is .TRUE..
*>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
*>          If WANTZ is .TRUE., on entry Z must contain the current
*>          matrix Z of transformations accumulated by DHSEQR, and on
*>          exit Z has been updated; transformations are applied only to
*>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*>          If WANTZ is .FALSE., Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           = 0:  successful exit
*>           > 0:  If INFO = i, DLAHQR failed to compute all the
*>                  eigenvalues ILO to IHI in a total of 30 iterations
*>                  per eigenvalue; elements i+1:ihi of WR and WI
*>                  contain those eigenvalues which have been
*>                  successfully computed.
*>
*>                  If INFO > 0 and WANTT is .FALSE., then on exit,
*>                  the remaining unconverged eigenvalues are the
*>                  eigenvalues of the upper Hessenberg matrix rows
*>                  and columns ILO through INFO of the final, output
*>                  value of H.
*>
*>                  If INFO > 0 and WANTT is .TRUE., then on exit
*>          (*)       (initial value of H)*U  = U*(final value of H)
*>                  where U is an orthogonal matrix.    The final
*>                  value of H is upper Hessenberg and triangular in
*>                  rows and columns INFO+1 through IHI.
*>
*>                  If INFO > 0 and WANTZ is .TRUE., then on exit
*>                      (final value of Z)  = (initial value of Z)*U
*>                  where U is the orthogonal matrix in (*)
*>                  (regardless of the value of WANTT.)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     02-96 Based on modifications by
*>     David Day, Sandia National Laboratory, USA
*>
*>     12-04 Further modifications by
*>     Ralph Byers, University of Kansas, USA
*>     This is a modified version of DLAHQR from LAPACK version 3.0.
*>     It is (1) more robust against overflow and underflow and
*>     (2) adopts the more conservative Ahues & Tisseur stopping
*>     criterion (LAWN 122, 1997).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
*     ..
*
*  =========================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0 )
      DOUBLE PRECISION   DAT1, DAT2
      PARAMETER          ( DAT1 = 3.0d0 / 4.0d0, DAT2 = -0.4375d0 )
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 10 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S,
     $                   H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX,
     $                   SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST,
     $                   ULP, V2, V3
      INTEGER            I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ,
     $                   KDEFL 
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLABAD, DLANV2, DLARFG, DROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 )
     $   H( IHI, IHI-2 ) = ZERO
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( NH ) / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITMAX is the total number of QR iterations allowed.
*
      ITMAX = 30 * MAX( 10, NH )
*
*     KDEFL counts the number of iterations since a deflation
*
      KDEFL = 0
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   20 CONTINUE
      L = ILO
      IF( I.LT.ILO )
     $   GO TO 160
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 140 ITS = 0, ITMAX
*
*        Look for a single small subdiagonal element.
*
         DO 30 K = I, L + 1, -1
            IF( ABS( H( K, K-1 ) ).LE.SMLNUM )
     $         GO TO 40
            TST = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO )
     $            TST = TST + ABS( H( K-1, K-2 ) )
               IF( K+1.LE.IHI )
     $            TST = TST + ABS( H( K+1, K ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some cases.  ====
            IF( ABS( H( K, K-1 ) ).LE.ULP*TST ) THEN
               AB = MAX( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               BA = MIN( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               AA = MAX( ABS( H( K, K ) ),
     $              ABS( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( ABS( H( K, K ) ),
     $              ABS( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM,
     $             ULP*( BB*( AA / S ) ) ) )GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 or 2 has split off.
*
         IF( L.GE.I-1 )
     $      GO TO 150
         KDEFL = KDEFL + 1
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( MOD(KDEFL,2*KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
            H11 = DAT1*S + H( I, I )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE IF( MOD(KDEFL,KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = ABS( H( L+1, L ) ) + ABS( H( L+2, L+1 ) )
            H11 = DAT1*S + H( L, L )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE
*
*           Prepare to use Francis' double shift
*           (i.e. 2nd degree generalized Rayleigh quotient)
*
            H11 = H( I-1, I-1 )
            H21 = H( I, I-1 )
            H12 = H( I-1, I )
            H22 = H( I, I )
         END IF
         S = ABS( H11 ) + ABS( H12 ) + ABS( H21 ) + ABS( H22 )
         IF( S.EQ.ZERO ) THEN
            RT1R = ZERO
            RT1I = ZERO
            RT2R = ZERO
            RT2I = ZERO
         ELSE
            H11 = H11 / S
            H21 = H21 / S
            H12 = H12 / S
            H22 = H22 / S
            TR = ( H11+H22 ) / TWO
            DET = ( H11-TR )*( H22-TR ) - H12*H21
            RTDISC = SQRT( ABS( DET ) )
            IF( DET.GE.ZERO ) THEN
*
*              ==== complex conjugate shifts ====
*
               RT1R = TR*S
               RT2R = RT1R
               RT1I = RTDISC*S
               RT2I = -RT1I
            ELSE
*
*              ==== real shifts (use only one of them)  ====
*
               RT1R = TR + RTDISC
               RT2R = TR - RTDISC
               IF( ABS( RT1R-H22 ).LE.ABS( RT2R-H22 ) ) THEN
                  RT1R = RT1R*S
                  RT2R = RT1R
               ELSE
                  RT2R = RT2R*S
                  RT1R = RT2R
               END IF
               RT1I = ZERO
               RT2I = ZERO
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 50 M = I - 2, L, -1
*           Determine the effect of starting the double-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.  (The following uses scaling to avoid
*           overflows and most underflows.)
*
            H21S = H( M+1, M )
            S = ABS( H( M, M )-RT2R ) + ABS( RT2I ) + ABS( H21S )
            H21S = H( M+1, M ) / S
            V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )*
     $               ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S )
            V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R )
            V( 3 ) = H21S*H( M+2, M+1 )
            S = ABS( V( 1 ) ) + ABS( V( 2 ) ) + ABS( V( 3 ) )
            V( 1 ) = V( 1 ) / S
            V( 2 ) = V( 2 ) / S
            V( 3 ) = V( 3 ) / S
            IF( M.EQ.L )
     $         GO TO 60
            IF( ABS( H( M, M-1 ) )*( ABS( V( 2 ) )+ABS( V( 3 ) ) ).LE.
     $          ULP*ABS( V( 1 ) )*( ABS( H( M-1, M-1 ) )+ABS( H( M,
     $          M ) )+ABS( H( M+1, M+1 ) ) ) )GO TO 60
   50    CONTINUE
   60    CONTINUE
*
*        Double-shift QR step
*
         DO 130 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix. NR is the order of G.
*
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M )
     $         CALL DCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 )
     $            H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
*               ==== Use the following instead of
*               .    H( K, K-1 ) = -H( K, K-1 ) to
*               .    avoid a bug when v(2) and v(3)
*               .    underflow. ====
               H( K, K-1 ) = H( K, K-1 )*( ONE-T1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 70 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   70          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 80 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   80          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 90 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   90             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 100 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
  100          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 110 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  110          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 120 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  120             CONTINUE
               END IF
            END IF
  130    CONTINUE
*
  140 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  150 CONTINUE
*
      IF( L.EQ.I ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*
         WR( I ) = H( I, I )
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
*        Transform the 2-by-2 submatrix to standard Schur form,
*        and compute and store the eigenvalues.
*
         CALL DLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ),
     $                H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ),
     $                CS, SN )
*
         IF( WANTT ) THEN
*
*           Apply the transformation to the rest of H.
*
            IF( I2.GT.I )
     $         CALL DROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH,
     $                    CS, SN )
            CALL DROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
         END IF
         IF( WANTZ ) THEN
*
*           Apply the transformation to Z.
*
            CALL DROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
         END IF
      END IF
*     reset deflation counter
      KDEFL = 0
*
*     return to start of the main loop with new value of I.
*
      I = L - 1
      GO TO 20
*
  160 CONTINUE
      RETURN
*
*     End of DLAHQR
*
      END
*> \brief \b DSYCONV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DSYCONV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyconv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyconv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyconv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYCONV convert A given by TRF into L and D and vice-versa.
*> Get Non-diag elements of D (returned in workspace) and
*> apply or reverse permutation done in TRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by DSYTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by DSYTRF.
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N)
*>          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1
*>          or 2-by-2 block diagonal matrix D in LDLT.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleSYcomputational
*
*  =====================================================================
      SUBROUTINE DSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, J
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYCONV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*      A is UPPER
*
*      Convert A (A is upper)
*
*        Convert VALUE
*
         IF ( CONVERT ) THEN
            I=N
            E(1)=ZERO
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I-1,I)
                  E(I-1)=ZERO
                  A(I-1,I)=ZERO
                  I=I-1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I-1
            END DO
*
*        Convert PERMUTATIONS
*
         I=N
         DO WHILE ( I .GE. 1 )
            IF( IPIV(I) .GT. 0) THEN
               IP=IPIV(I)
               IF( I .LT. N) THEN
                  DO 12 J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
 12            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
               IF( I .LT. N) THEN
             DO 13 J= I+1,N
                 TEMP=A(IP,J)
                 A(IP,J)=A(I-1,J)
                 A(I-1,J)=TEMP
 13            CONTINUE
                ENDIF
                I=I-1
           ENDIF
           I=I-1
        END DO

         ELSE
*
*      Revert A (A is upper)
*
*
*        Revert PERMUTATIONS
*
            I=1
            DO WHILE ( I .LE. N )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF( I .LT. N) THEN
                  DO J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
                  END DO
                  ENDIF
               ELSE
                 IP=-IPIV(I)
                 I=I+1
                 IF( I .LT. N) THEN
                    DO J= I+1,N
                       TEMP=A(IP,J)
                       A(IP,J)=A(I-1,J)
                       A(I-1,J)=TEMP
                    END DO
                 ENDIF
               ENDIF
               I=I+1
            END DO
*
*        Revert VALUE
*
            I=N
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I-1,I)=E(I)
                  I=I-1
               ENDIF
               I=I-1
            END DO
         END IF
      ELSE
*
*      A is LOWER
*
         IF ( CONVERT ) THEN
*
*      Convert A (A is lower)
*
*
*        Convert VALUE
*
            I=1
            E(N)=ZERO
            DO WHILE ( I .LE. N )
               IF( I.LT.N .AND. IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I+1,I)
                  E(I+1)=ZERO
                  A(I+1,I)=ZERO
                  I=I+1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I+1
            END DO
*
*        Convert PERMUTATIONS
*
         I=1
         DO WHILE ( I .LE. N )
            IF( IPIV(I) .GT. 0 ) THEN
               IP=IPIV(I)
               IF (I .GT. 1) THEN
               DO 22 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I,J)
                 A(I,J)=TEMP
 22            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
              IF (I .GT. 1) THEN
              DO 23 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I+1,J)
                 A(I+1,J)=TEMP
 23           CONTINUE
              ENDIF
              I=I+1
           ENDIF
           I=I+1
        END DO
         ELSE
*
*      Revert A (A is lower)
*
*
*        Revert PERMUTATIONS
*
            I=N
            DO WHILE ( I .GE. 1 )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I,J)
                        A(I,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ELSE
                  IP=-IPIV(I)
                  I=I-1
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I+1,J)
                        A(I+1,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ENDIF
               I=I-1
            END DO
*
*        Revert VALUE
*
            I=1
            DO WHILE ( I .LE. N-1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I+1,I)=E(I)
                  I=I+1
               ENDIF
               I=I+1
            END DO
         END IF
      END IF

      RETURN
*
*     End of DSYCONV
*
      END
*> \brief \b DSYCONVF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DSYCONVF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyconvf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyconvf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyconvf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> DSYCONVF converts the factorization output format used in
*> DSYTRF provided on entry in parameter A into the factorization
*> output format used in DSYTRF_RK (or DSYTRF_BK) that is stored
*> on exit in parameters A and E. It also converts in place details of
*> the intechanges stored in IPIV from the format used in DSYTRF into
*> the format used in DSYTRF_RK (or DSYTRF_BK).
*>
*> If parameter WAY = 'R':
*> DSYCONVF performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in DSYTRF_RK
*> (or DSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in DSYTRF that is stored
*> on exit in parameter A. It also converts in place details of
*> the intechanges stored in IPIV from the format used in DSYTRF_RK
*> (or DSYTRF_BK) into the format used in DSYTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          DSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          DSYTRF_RK or DSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          DSYTRF_RK or DSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          DSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in,out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>
*>          1) If WAY ='C':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in DSYTRF.
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in DSYTRF_RK
*>          ( or DSYTRF_BK).
*>
*>          1) If WAY ='R':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in DSYTRF_RK
*>          ( or DSYTRF_BK).
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in DSYTRF.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleSYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE DSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           DSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYCONVF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL DSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL DSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i-1 and IPIV(i-1),
*                 so this should be recorded in two consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I-1 )
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where k increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL DSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL DSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i+1 and IPIV(i+1),
*                 so this should be recorded in consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I+1 )
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of DSYCONVF
*
      END
*> \brief \b DSYCONVF_ROOK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DSYCONVF_ROOK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyconvf_rook.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyconvf_rook.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyconvf_rook.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> DSYCONVF_ROOK converts the factorization output format used in
*> DSYTRF_ROOK provided on entry in parameter A into the factorization
*> output format used in DSYTRF_RK (or DSYTRF_BK) that is stored
*> on exit in parameters A and E. IPIV format for DSYTRF_ROOK and
*> DSYTRF_RK (or DSYTRF_BK) is the same and is not converted.
*>
*> If parameter WAY = 'R':
*> DSYCONVF_ROOK performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in DSYTRF_RK
*> (or DSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in DSYTRF_ROOK that is stored
*> on exit in parameter A. IPIV format for DSYTRF_ROOK and
*> DSYTRF_RK (or DSYTRF_BK) is the same and is not converted.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          DSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          DSYTRF_RK or DSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          DSYTRF_RK or DSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          DSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          On entry, details of the interchanges and the block
*>          structure of D as determined:
*>          1) by DSYTRF_ROOK, if WAY ='C';
*>          2) by DSYTRF_RK (or DSYTRF_BK), if WAY ='R'.
*>          The IPIV format is the same for all these routines.
*>
*>          On exit, is not changed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleSYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE DSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           DSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, IP2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYCONVF_ROOK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1)
*                 in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                     IF( IP2.NE.(I-1) ) THEN
                        CALL DSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP2, I+1 ), LDA )
                     END IF
                  END IF
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i)
*                 in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP2.NE.(I-1) ) THEN
                        CALL DSWAP( N-I, A( IP2, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1)
*                 in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                     IF( IP2.NE.(I+1) ) THEN
                        CALL DSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP2, 1 ), LDA )
                     END IF
                  END IF
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i)
*                 in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP2.NE.(I+1) ) THEN
                        CALL DSWAP( I-1, A( IP2, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL DSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of DSYCONVF_ROOK
*
      END
*> \brief \b DTPMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DTPMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpmqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpmqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpmqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
*                           A, LDA, B, LDB, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   V( LDV, * ), A( LDA, * ), B( LDB, * ),
*      $                   T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTPMQRT applies a real orthogonal matrix Q obtained from a
*> "triangular-pentagonal" real block reflector H to a general
*> real matrix C, which consists of two blocks A and B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**T from the Left;
*>          = 'R': apply Q or Q**T from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'T':  Transpose, apply Q**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The order of the trapezoidal part of V.
*>          K >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CTPQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CTPQRT in B.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDV >= max(1,M);
*>          if SIDE = 'R', LDV >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CTPQRT, stored as a NB-by-K matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension
*>          (LDA,N) if SIDE = 'L' or
*>          (LDA,K) if SIDE = 'R'
*>          On entry, the K-by-N or M-by-K matrix A.
*>          On exit, A is overwritten by the corresponding block of
*>          Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDC >= max(1,K);
*>          If SIDE = 'R', LDC >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,N)
*>          On entry, the M-by-N matrix B.
*>          On exit, B is overwritten by the corresponding block of
*>          Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.
*>          LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array. The dimension of WORK is
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The columns of the pentagonal matrix V contain the elementary reflectors
*>  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
*>  trapezoidal block V2:
*>
*>        V = [V1]
*>            [V2].
*>
*>  The size of the trapezoidal block V2 is determined by the parameter L,
*>  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L
*>  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular;
*>  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
*>
*>  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K.
*>                      [B]
*>
*>  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K.
*>
*>  The real orthogonal matrix Q is formed from V and T.
*>
*>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.
*>
*>  If TRANS='T' and SIDE='L', C is on exit replaced with Q**T * C.
*>
*>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.
*>
*>  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q**T.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
     $                    A, LDA, B, LDB, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   V( LDV, * ), A( LDA, * ), B( LDB, * ),
     $                   T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, MB, LB, KF, LDAQ, LDVQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTPRFB, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'T' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF ( LEFT ) THEN
         LDVQ = MAX( 1, M )
         LDAQ = MAX( 1, K )
      ELSE IF ( RIGHT ) THEN
         LDVQ = MAX( 1, N )
         LDAQ = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      ELSE IF( L.LT.0 .OR. L.GT.K ) THEN
         INFO = -6
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0) ) THEN
         INFO = -7
      ELSE IF( LDV.LT.LDVQ ) THEN
         INFO = -9
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -11
      ELSE IF( LDA.LT.LDAQ ) THEN
         INFO = -13
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTPMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL DTPRFB( 'L', 'T', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL DTPRFB( 'R', 'N', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL DTPRFB( 'L', 'N', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL DTPRFB( 'R', 'T', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      END IF
*
      RETURN
*
*     End of DTPMQRT
*
      END
*> \brief \b DTPQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DTPQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTPQRT computes a blocked QR factorization of a real
*> "triangular-pentagonal" matrix C, which is composed of a
*> triangular block A and pentagonal block B, using the compact
*> WY representation for Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B.
*>          M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B, and the order of the
*>          triangular matrix A.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The number of rows of the upper trapezoidal part of B.
*>          MIN(M,N) >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  N >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the upper triangular N-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the upper triangular matrix R.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,N)
*>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
*>          are rectangular, and the last L rows are upper trapezoidal.
*>          On exit, B contains the pentagonal matrix V.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,N)
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The input matrix C is a (N+M)-by-N matrix
*>
*>               C = [ A ]
*>                   [ B ]
*>
*>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
*>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
*>  upper trapezoidal matrix B2:
*>
*>               B = [ B1 ]  <- (M-L)-by-N rectangular
*>                   [ B2 ]  <-     L-by-N upper trapezoidal.
*>
*>  The upper trapezoidal matrix B2 consists of the first L rows of a
*>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
*>  B is rectangular M-by-N; if M=L=N, B is upper triangular.
*>
*>  The matrix W stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal (of A) in the (N+M)-by-N input matrix C
*>
*>               C = [ A ]  <- upper triangular N-by-N
*>                   [ B ]  <- M-by-N pentagonal
*>
*>  so that W can be represented as
*>
*>               W = [ I ]  <- identity, N-by-N
*>                   [ V ]  <- M-by-N, same form as B.
*>
*>  Thus, all of information needed for W is contained on exit in B, which
*>  we call V above.  Note that V has the same form as B; that is,
*>
*>               V = [ V1 ] <- (M-L)-by-N rectangular
*>                   [ V2 ] <-     L-by-N upper trapezoidal.
*>
*>  The columns of V represent the vectors which define the H(i)'s.
*>
*>  The number of blocks is B = ceiling(N/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-N matrix T as
*>
*>               T = [T1 T2 ... TB].
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, LB, MB, IINFO
*     ..
*     .. External Subroutines ..
      EXTERNAL   DTPQRT2, DTPRFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 .OR. (L.GT.MIN(M,N) .AND. MIN(M,N).GE.0)) THEN
         INFO = -3
      ELSE IF( NB.LT.1 .OR. (NB.GT.N .AND. N.GT.0)) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTPQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      DO I = 1, N, NB
*
*     Compute the QR factorization of the current block
*
         IB = MIN( N-I+1, NB )
         MB = MIN( M-L+I+IB-1, M )
         IF( I.GE.L ) THEN
            LB = 0
         ELSE
            LB = MB-M+L-I+1
         END IF
*
         CALL DTPQRT2( MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB,
     $                 T(1, I ), LDT, IINFO )
*
*     Update by applying H**T to B(:,I+IB:N) from the left
*
         IF( I+IB.LE.N ) THEN
            CALL DTPRFB( 'L', 'T', 'F', 'C', MB, N-I-IB+1, IB, LB,
     $                    B( 1, I ), LDB, T( 1, I ), LDT,
     $                    A( I, I+IB ), LDA, B( 1, I+IB ), LDB,
     $                    WORK, IB )
         END IF
      END DO
      RETURN
*
*     End of DTPQRT
*
      END
*> \brief \b SGEMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGEMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgemqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgemqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgemqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
*                          C, LDC, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*       ..
*       .. Array Arguments ..
*       REAL   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGEMQRT overwrites the general real M-by-N matrix C with
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q C            C Q
*> TRANS = 'T':   Q**T C            C Q**T
*>
*> where Q is a real orthogonal matrix defined as the product of K
*> elementary reflectors:
*>
*>       Q = H(1) H(2) . . . H(K) = I - V T V**T
*>
*> generated using the compact WY representation as returned by SGEQRT.
*>
*> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**T from the Left;
*>          = 'R': apply Q or Q**T from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'T':  Transpose, apply Q**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CGEQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is REAL array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CGEQRT in the first K columns of its array argument A.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is REAL array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CGEQRT, stored as a NB-by-N matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is REAL array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array. The dimension of WORK is
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realGEcomputational
*
*  =====================================================================
      SUBROUTINE SGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
     $                   C, LDC, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*     ..
*     .. Array Arguments ..
      REAL   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, LDWORK, KF, Q
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, SLARFB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'T' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF( LEFT ) THEN
         LDWORK = MAX( 1, N )
         Q = M
      ELSE IF ( RIGHT ) THEN
         LDWORK = MAX( 1, M )
         Q = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.Q ) THEN
         INFO = -5
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0)) THEN
         INFO = -6
      ELSE IF( LDV.LT.MAX( 1, Q ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL SLARFB( 'L', 'T', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL SLARFB( 'R', 'N', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL SLARFB( 'L', 'N', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL SLARFB( 'R', 'T', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      END IF
*
      RETURN
*
*     End of SGEMQRT
*
      END
*> \brief \b SGEQRFP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGEQRFP + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrfp.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrfp.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrfp.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGEQR2P computes a QR factorization of a real M-by-N matrix A:
*>
*>    A = Q * ( R ),
*>            ( 0 )
*>
*> where:
*>
*>    Q is a M-by-M orthogonal matrix;
*>    R is an upper-triangular N-by-N matrix with nonnegative diagonal
*>    entries;
*>    0 is a (M-N)-by-N zero matrix, if M > N.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if m >= n). The diagonal entries of R
*>          are nonnegative; the elements below the diagonal,
*>          with the array TAU, represent the orthogonal matrix Q as a
*>          product of min(m,n) elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is REAL array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*>  and tau in TAU(i).
*>
*> See Lapack Working Note 203 for details
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQR2P, SLARFB, SLARFT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQRFP', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'SGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL SGEQR2P( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H**T to A(i:m,i+ib:n) from the left
*
               CALL SLARFB( 'Left', 'Transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL SGEQR2P( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of SGEQRFP
*
      END
*> \brief \b SGEQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGEQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDT, M, N, NB
*       ..
*       .. Array Arguments ..
*       REAL A( LDA, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGEQRT computes a blocked QR factorization of a real M-by-N matrix A
*> using the compact WY representation of Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if M >= N); the elements below the diagonal
*>          are the columns of V.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is REAL array, dimension (LDT,MIN(M,N))
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See below
*>          for further details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix V stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal. For example, if M=5 and N=3, the matrix V is
*>
*>               V = (  1       )
*>                   ( v1  1    )
*>                   ( v1 v2  1 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  where the vi's represent the vectors which define H(i), which are returned
*>  in the matrix A.  The 1's along the diagonal of V are not stored in A.
*>
*>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-K matrix T as
*>
*>               T = (T1 T2 ... TB).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDT, M, N, NB
*     ..
*     .. Array Arguments ..
      REAL A( LDA, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, IINFO, K
      LOGICAL    USE_RECURSIVE_QR
      PARAMETER( USE_RECURSIVE_QR=.TRUE. )
*     ..
*     .. External Subroutines ..
      EXTERNAL   SGEQRT2, SGEQRT3, SLARFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NB.LT.1 .OR. ( NB.GT.MIN(M,N) .AND. MIN(M,N).GT.0 ) )THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN
*
*     Blocked loop of length K
*
      DO I = 1, K,  NB
         IB = MIN( K-I+1, NB )
*
*     Compute the QR factorization of the current block A(I:M,I:I+IB-1)
*
         IF( USE_RECURSIVE_QR ) THEN
            CALL SGEQRT3( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         ELSE
            CALL SGEQRT2( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         END IF
         IF( I+IB.LE.N ) THEN
*
*     Update by applying H**T to A(I:M,I+IB:N) from the left
*
            CALL SLARFB( 'L', 'T', 'F', 'C', M-I+1, N-I-IB+1, IB,
     $                   A( I, I ), LDA, T( 1, I ), LDT,
     $                   A( I, I+IB ), LDA, WORK , N-I-IB+1 )
         END IF
      END DO
      RETURN
*
*     End of SGEQRT
*
      END
*> \brief \b SLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SLAHQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
*                          ILOZ, IHIZ, Z, LDZ, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       REAL               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    SLAHQR is an auxiliary routine called by SHSEQR to update the
*>    eigenvalues and Schur decomposition already computed by SHSEQR, by
*>    dealing with the Hessenberg submatrix in rows and columns ILO to
*>    IHI.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          = .TRUE. : the full Schur form T is required;
*>          = .FALSE.: only eigenvalues are required.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE. : the matrix of Schur vectors Z is required;
*>          = .FALSE.: Schur vectors are not required.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix H.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          It is assumed that H is already upper quasi-triangular in
*>          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
*>          ILO = 1). SLAHQR works primarily with the Hessenberg
*>          submatrix in rows and columns ILO to IHI, but applies
*>          transformations to all of H if WANTT is .TRUE..
*>          1 <= ILO <= max(1,IHI); IHI <= N.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is REAL array, dimension (LDH,N)
*>          On entry, the upper Hessenberg matrix H.
*>          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
*>          quasi-triangular in rows and columns ILO:IHI, with any
*>          2-by-2 diagonal blocks in standard form. If INFO is zero
*>          and WANTT is .FALSE., the contents of H are unspecified on
*>          exit.  The output state of H if INFO is nonzero is given
*>          below under the description of INFO.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          The leading dimension of the array H. LDH >= max(1,N).
*> \endverbatim
*>
*> \param[out] WR
*> \verbatim
*>          WR is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] WI
*> \verbatim
*>          WI is REAL array, dimension (N)
*>          The real and imaginary parts, respectively, of the computed
*>          eigenvalues ILO to IHI are stored in the corresponding
*>          elements of WR and WI. If two eigenvalues are computed as a
*>          complex conjugate pair, they are stored in consecutive
*>          elements of WR and WI, say the i-th and (i+1)th, with
*>          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
*>          eigenvalues are stored in the same order as on the diagonal
*>          of the Schur form returned in H, with WR(i) = H(i,i), and, if
*>          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
*>          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>          Specify the rows of Z to which transformations must be
*>          applied if WANTZ is .TRUE..
*>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ,N)
*>          If WANTZ is .TRUE., on entry Z must contain the current
*>          matrix Z of transformations accumulated by SHSEQR, and on
*>          exit Z has been updated; transformations are applied only to
*>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*>          If WANTZ is .FALSE., Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           = 0:   successful exit
*>           > 0:   If INFO = i, SLAHQR failed to compute all the
*>                  eigenvalues ILO to IHI in a total of 30 iterations
*>                  per eigenvalue; elements i+1:ihi of WR and WI
*>                  contain those eigenvalues which have been
*>                  successfully computed.
*>
*>                  If INFO > 0 and WANTT is .FALSE., then on exit,
*>                  the remaining unconverged eigenvalues are the
*>                  eigenvalues of the upper Hessenberg matrix rows
*>                  and columns ILO through INFO of the final, output
*>                  value of H.
*>
*>                  If INFO > 0 and WANTT is .TRUE., then on exit
*>          (*)       (initial value of H)*U  = U*(final value of H)
*>                  where U is an orthogonal matrix.    The final
*>                  value of H is upper Hessenberg and triangular in
*>                  rows and columns INFO+1 through IHI.
*>
*>                  If INFO > 0 and WANTZ is .TRUE., then on exit
*>                      (final value of Z)  = (initial value of Z)*U
*>                  where U is the orthogonal matrix in (*)
*>                  (regardless of the value of WANTT.)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     02-96 Based on modifications by
*>     David Day, Sandia National Laboratory, USA
*>
*>     12-04 Further modifications by
*>     Ralph Byers, University of Kansas, USA
*>     This is a modified version of SLAHQR from LAPACK version 3.0.
*>     It is (1) more robust against overflow and underflow and
*>     (2) adopts the more conservative Ahues & Tisseur stopping
*>     criterion (LAWN 122, 1997).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      REAL               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
*     ..
*
*  =========================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0e0, ONE = 1.0e0, TWO = 2.0e0 )
      REAL               DAT1, DAT2
      PARAMETER          ( DAT1 = 3.0e0 / 4.0e0, DAT2 = -0.4375e0 )
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 10 )
*     ..
*     .. Local Scalars ..
      REAL               AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S,
     $                   H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX,
     $                   SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST,
     $                   ULP, V2, V3
      INTEGER            I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ,
     $                   KDEFL
*     ..
*     .. Local Arrays ..
      REAL               V( 3 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLABAD, SLANV2, SLARFG, SROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 )
     $   H( IHI, IHI-2 ) = ZERO
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL SLABAD( SAFMIN, SAFMAX )
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( NH ) / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITMAX is the total number of QR iterations allowed.
*
      ITMAX = 30 * MAX( 10, NH )
*
*     KDEFL counts the number of iterations since a deflation
*
      KDEFL = 0
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   20 CONTINUE
      L = ILO
      IF( I.LT.ILO )
     $   GO TO 160
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 140 ITS = 0, ITMAX
*
*        Look for a single small subdiagonal element.
*
         DO 30 K = I, L + 1, -1
            IF( ABS( H( K, K-1 ) ).LE.SMLNUM )
     $         GO TO 40
            TST = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO )
     $            TST = TST + ABS( H( K-1, K-2 ) )
               IF( K+1.LE.IHI )
     $            TST = TST + ABS( H( K+1, K ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some cases.  ====
            IF( ABS( H( K, K-1 ) ).LE.ULP*TST ) THEN
               AB = MAX( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               BA = MIN( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               AA = MAX( ABS( H( K, K ) ),
     $              ABS( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( ABS( H( K, K ) ),
     $              ABS( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM,
     $             ULP*( BB*( AA / S ) ) ) )GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 or 2 has split off.
*
         IF( L.GE.I-1 )
     $      GO TO 150
         KDEFL = KDEFL + 1
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( MOD(KDEFL,2*KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
            H11 = DAT1*S + H( I, I )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE IF( MOD(KDEFL,KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = ABS( H( L+1, L ) ) + ABS( H( L+2, L+1 ) )
            H11 = DAT1*S + H( L, L )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE
*
*           Prepare to use Francis' double shift
*           (i.e. 2nd degree generalized Rayleigh quotient)
*
            H11 = H( I-1, I-1 )
            H21 = H( I, I-1 )
            H12 = H( I-1, I )
            H22 = H( I, I )
         END IF
         S = ABS( H11 ) + ABS( H12 ) + ABS( H21 ) + ABS( H22 )
         IF( S.EQ.ZERO ) THEN
            RT1R = ZERO
            RT1I = ZERO
            RT2R = ZERO
            RT2I = ZERO
         ELSE
            H11 = H11 / S
            H21 = H21 / S
            H12 = H12 / S
            H22 = H22 / S
            TR = ( H11+H22 ) / TWO
            DET = ( H11-TR )*( H22-TR ) - H12*H21
            RTDISC = SQRT( ABS( DET ) )
            IF( DET.GE.ZERO ) THEN
*
*              ==== complex conjugate shifts ====
*
               RT1R = TR*S
               RT2R = RT1R
               RT1I = RTDISC*S
               RT2I = -RT1I
            ELSE
*
*              ==== real shifts (use only one of them)  ====
*
               RT1R = TR + RTDISC
               RT2R = TR - RTDISC
               IF( ABS( RT1R-H22 ).LE.ABS( RT2R-H22 ) ) THEN
                  RT1R = RT1R*S
                  RT2R = RT1R
               ELSE
                  RT2R = RT2R*S
                  RT1R = RT2R
               END IF
               RT1I = ZERO
               RT2I = ZERO
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 50 M = I - 2, L, -1
*           Determine the effect of starting the double-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.  (The following uses scaling to avoid
*           overflows and most underflows.)
*
            H21S = H( M+1, M )
            S = ABS( H( M, M )-RT2R ) + ABS( RT2I ) + ABS( H21S )
            H21S = H( M+1, M ) / S
            V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )*
     $               ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S )
            V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R )
            V( 3 ) = H21S*H( M+2, M+1 )
            S = ABS( V( 1 ) ) + ABS( V( 2 ) ) + ABS( V( 3 ) )
            V( 1 ) = V( 1 ) / S
            V( 2 ) = V( 2 ) / S
            V( 3 ) = V( 3 ) / S
            IF( M.EQ.L )
     $         GO TO 60
            IF( ABS( H( M, M-1 ) )*( ABS( V( 2 ) )+ABS( V( 3 ) ) ).LE.
     $          ULP*ABS( V( 1 ) )*( ABS( H( M-1, M-1 ) )+ABS( H( M,
     $          M ) )+ABS( H( M+1, M+1 ) ) ) )GO TO 60
   50    CONTINUE
   60    CONTINUE
*
*        Double-shift QR step
*
         DO 130 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix. NR is the order of G.
*
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M )
     $         CALL SCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL SLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 )
     $            H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
*               ==== Use the following instead of
*               .    H( K, K-1 ) = -H( K, K-1 ) to
*               .    avoid a bug when v(2) and v(3)
*               .    underflow. ====
               H( K, K-1 ) = H( K, K-1 )*( ONE-T1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 70 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   70          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 80 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   80          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 90 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   90             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 100 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
  100          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 110 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  110          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 120 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  120             CONTINUE
               END IF
            END IF
  130    CONTINUE
*
  140 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  150 CONTINUE
*
      IF( L.EQ.I ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*
         WR( I ) = H( I, I )
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
*        Transform the 2-by-2 submatrix to standard Schur form,
*        and compute and store the eigenvalues.
*
         CALL SLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ),
     $                H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ),
     $                CS, SN )
*
         IF( WANTT ) THEN
*
*           Apply the transformation to the rest of H.
*
            IF( I2.GT.I )
     $         CALL SROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH,
     $                    CS, SN )
            CALL SROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
         END IF
         IF( WANTZ ) THEN
*
*           Apply the transformation to Z.
*
            CALL SROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
         END IF
      END IF
*     reset deflation counter
      KDEFL = 0
*
*     return to start of the main loop with new value of I.
*
      I = L - 1
      GO TO 20
*
  160 CONTINUE
      RETURN
*
*     End of SLAHQR
*
      END
*> \brief \b SSYCONV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SSYCONV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SSYCONV convert A given by TRF into L and D and vice-versa.
*> Get Non-diag elements of D (returned in workspace) and
*> apply or reverse permutation done in TRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by SSYTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by SSYTRF.
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is REAL array, dimension (N)
*>          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1
*>          or 2-by-2 block diagonal matrix D in LDLT.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realSYcomputational
*
*  =====================================================================
      SUBROUTINE SSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, J
      REAL               TEMP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYCONV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*      A is UPPER
*
*      Convert A (A is upper)
*
*        Convert VALUE
*
         IF ( CONVERT ) THEN
            I=N
            E(1)=ZERO
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I-1,I)
                  E(I-1)=ZERO
                  A(I-1,I)=ZERO
                  I=I-1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I-1
            END DO
*
*        Convert PERMUTATIONS
*
         I=N
         DO WHILE ( I .GE. 1 )
            IF( IPIV(I) .GT. 0) THEN
               IP=IPIV(I)
               IF( I .LT. N) THEN
                  DO 12 J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
 12            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
               IF( I .LT. N) THEN
             DO 13 J= I+1,N
                 TEMP=A(IP,J)
                 A(IP,J)=A(I-1,J)
                 A(I-1,J)=TEMP
 13            CONTINUE
                ENDIF
                I=I-1
           ENDIF
           I=I-1
        END DO

         ELSE
*
*      Revert A (A is upper)
*
*
*        Revert PERMUTATIONS
*
            I=1
            DO WHILE ( I .LE. N )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF( I .LT. N) THEN
                  DO J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
                  END DO
                  ENDIF
               ELSE
                 IP=-IPIV(I)
                 I=I+1
                 IF( I .LT. N) THEN
                    DO J= I+1,N
                       TEMP=A(IP,J)
                       A(IP,J)=A(I-1,J)
                       A(I-1,J)=TEMP
                    END DO
                 ENDIF
               ENDIF
               I=I+1
            END DO
*
*        Revert VALUE
*
            I=N
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I-1,I)=E(I)
                  I=I-1
               ENDIF
               I=I-1
            END DO
         END IF
      ELSE
*
*      A is LOWER
*
         IF ( CONVERT ) THEN
*
*      Convert A (A is lower)
*
*
*        Convert VALUE
*
            I=1
            E(N)=ZERO
            DO WHILE ( I .LE. N )
               IF( I.LT.N .AND. IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I+1,I)
                  E(I+1)=ZERO
                  A(I+1,I)=ZERO
                  I=I+1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I+1
            END DO
*
*        Convert PERMUTATIONS
*
         I=1
         DO WHILE ( I .LE. N )
            IF( IPIV(I) .GT. 0 ) THEN
               IP=IPIV(I)
               IF (I .GT. 1) THEN
               DO 22 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I,J)
                 A(I,J)=TEMP
 22            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
              IF (I .GT. 1) THEN
              DO 23 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I+1,J)
                 A(I+1,J)=TEMP
 23           CONTINUE
              ENDIF
              I=I+1
           ENDIF
           I=I+1
        END DO
         ELSE
*
*      Revert A (A is lower)
*
*
*        Revert PERMUTATIONS
*
            I=N
            DO WHILE ( I .GE. 1 )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I,J)
                        A(I,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ELSE
                  IP=-IPIV(I)
                  I=I-1
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I+1,J)
                        A(I+1,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ENDIF
               I=I-1
            END DO
*
*        Revert VALUE
*
            I=1
            DO WHILE ( I .LE. N-1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I+1,I)=E(I)
                  I=I+1
               ENDIF
               I=I+1
            END DO
         END IF
      END IF

      RETURN
*
*     End of SSYCONV
*
      END
*> \brief \b SSYCONVF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SSYCONVF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconvf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconvf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconvf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> SSYCONVF converts the factorization output format used in
*> SSYTRF provided on entry in parameter A into the factorization
*> output format used in SSYTRF_RK (or SSYTRF_BK) that is stored
*> on exit in parameters A and E. It also converts in place details of
*> the intechanges stored in IPIV from the format used in SSYTRF into
*> the format used in SSYTRF_RK (or SSYTRF_BK).
*>
*> If parameter WAY = 'R':
*> SSYCONVF performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in SSYTRF_RK
*> (or SSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in SSYTRF that is stored
*> on exit in parameter A. It also converts in place details of
*> the intechanges stored in IPIV from the format used in SSYTRF_RK
*> (or SSYTRF_BK) into the format used in SSYTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          SSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          SSYTRF_RK or SSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          SSYTRF_RK or SSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          SSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is REAL array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in,out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>
*>          1) If WAY ='C':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in SSYTRF.
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in SSYTRF_RK
*>          ( or SSYTRF_BK).
*>
*>          1) If WAY ='R':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in SSYTRF_RK
*>          ( or SSYTRF_BK).
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in SSYTRF.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup singleSYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE SSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           SSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYCONVF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL SSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL SSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i-1 and IPIV(i-1),
*                 so this should be recorded in two consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I-1 )
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where k increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL SSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL SSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i+1 and IPIV(i+1),
*                 so this should be recorded in consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I+1 )
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of SSYCONVF
*
      END
*> \brief \b SSYCONVF_ROOK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SSYCONVF_ROOK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconvf_rook.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconvf_rook.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconvf_rook.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> SSYCONVF_ROOK converts the factorization output format used in
*> SSYTRF_ROOK provided on entry in parameter A into the factorization
*> output format used in SSYTRF_RK (or SSYTRF_BK) that is stored
*> on exit in parameters A and E. IPIV format for SSYTRF_ROOK and
*> SSYTRF_RK (or SSYTRF_BK) is the same and is not converted.
*>
*> If parameter WAY = 'R':
*> SSYCONVF_ROOK performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in SSYTRF_RK
*> (or SSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in SSYTRF_ROOK that is stored
*> on exit in parameter A. IPIV format for SSYTRF_ROOK and
*> SSYTRF_RK (or SSYTRF_BK) is the same and is not converted.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          SSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          SSYTRF_RK or SSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          SSYTRF_RK or SSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          SSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is REAL array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          On entry, details of the interchanges and the block
*>          structure of D as determined:
*>          1) by SSYTRF_ROOK, if WAY ='C';
*>          2) by SSYTRF_RK (or SSYTRF_BK), if WAY ='R'.
*>          The IPIV format is the same for all these routines.
*>
*>          On exit, is not changed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup singleSYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE SSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           SSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, IP2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYCONVF_ROOK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1)
*                 in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                     IF( IP2.NE.(I-1) ) THEN
                        CALL SSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP2, I+1 ), LDA )
                     END IF
                  END IF
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i)
*                 in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP2.NE.(I-1) ) THEN
                        CALL SSWAP( N-I, A( IP2, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1)
*                 in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                     IF( IP2.NE.(I+1) ) THEN
                        CALL SSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP2, 1 ), LDA )
                     END IF
                  END IF
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i)
*                 in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP2.NE.(I+1) ) THEN
                        CALL SSWAP( I-1, A( IP2, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL SSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of SSYCONVF_ROOK
*
      END
*> \brief \b STPMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download STPMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpmqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpmqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpmqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
*                           A, LDA, B, LDB, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*       ..
*       .. Array Arguments ..
*       REAL   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
*      $          WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> STPMQRT applies a real orthogonal matrix Q obtained from a
*> "triangular-pentagonal" real block reflector H to a general
*> real matrix C, which consists of two blocks A and B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q^T from the Left;
*>          = 'R': apply Q or Q^T from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'T':  Transpose, apply Q^T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The order of the trapezoidal part of V.
*>          K >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CTPQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is REAL array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CTPQRT in B.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDV >= max(1,M);
*>          if SIDE = 'R', LDV >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is REAL array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CTPQRT, stored as a NB-by-K matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension
*>          (LDA,N) if SIDE = 'L' or
*>          (LDA,K) if SIDE = 'R'
*>          On entry, the K-by-N or M-by-K matrix A.
*>          On exit, A is overwritten by the corresponding block of
*>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDC >= max(1,K);
*>          If SIDE = 'R', LDC >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,N)
*>          On entry, the M-by-N matrix B.
*>          On exit, B is overwritten by the corresponding block of
*>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.
*>          LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array. The dimension of WORK is
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The columns of the pentagonal matrix V contain the elementary reflectors
*>  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
*>  trapezoidal block V2:
*>
*>        V = [V1]
*>            [V2].
*>
*>  The size of the trapezoidal block V2 is determined by the parameter L,
*>  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L
*>  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular;
*>  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
*>
*>  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K.
*>                      [B]
*>
*>  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K.
*>
*>  The real orthogonal matrix Q is formed from V and T.
*>
*>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.
*>
*>  If TRANS='T' and SIDE='L', C is on exit replaced with Q^T * C.
*>
*>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.
*>
*>  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q^T.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
     $                    A, LDA, B, LDB, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*     ..
*     .. Array Arguments ..
      REAL   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
     $          WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, MB, LB, KF, LDAQ, LDVQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           STPRFB, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'T' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF ( LEFT ) THEN
         LDVQ = MAX( 1, M )
         LDAQ = MAX( 1, K )
      ELSE IF ( RIGHT ) THEN
         LDVQ = MAX( 1, N )
         LDAQ = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      ELSE IF( L.LT.0 .OR. L.GT.K ) THEN
         INFO = -6
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0) ) THEN
         INFO = -7
      ELSE IF( LDV.LT.LDVQ ) THEN
         INFO = -9
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -11
      ELSE IF( LDA.LT.LDAQ ) THEN
         INFO = -13
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STPMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL STPRFB( 'L', 'T', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL STPRFB( 'R', 'N', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL STPRFB( 'L', 'N', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL STPRFB( 'R', 'T', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      END IF
*
      RETURN
*
*     End of STPMQRT
*
      END
*> \brief \b STPQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download STPQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE STPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*       ..
*       .. Array Arguments ..
*       REAL A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> STPQRT computes a blocked QR factorization of a real
*> "triangular-pentagonal" matrix C, which is composed of a
*> triangular block A and pentagonal block B, using the compact
*> WY representation for Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B.
*>          M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B, and the order of the
*>          triangular matrix A.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The number of rows of the upper trapezoidal part of B.
*>          MIN(M,N) >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  N >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the upper triangular N-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the upper triangular matrix R.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,N)
*>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
*>          are rectangular, and the last L rows are upper trapezoidal.
*>          On exit, B contains the pentagonal matrix V.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is REAL array, dimension (LDT,N)
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The input matrix C is a (N+M)-by-N matrix
*>
*>               C = [ A ]
*>                   [ B ]
*>
*>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
*>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
*>  upper trapezoidal matrix B2:
*>
*>               B = [ B1 ]  <- (M-L)-by-N rectangular
*>                   [ B2 ]  <-     L-by-N upper trapezoidal.
*>
*>  The upper trapezoidal matrix B2 consists of the first L rows of a
*>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
*>  B is rectangular M-by-N; if M=L=N, B is upper triangular.
*>
*>  The matrix W stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal (of A) in the (N+M)-by-N input matrix C
*>
*>               C = [ A ]  <- upper triangular N-by-N
*>                   [ B ]  <- M-by-N pentagonal
*>
*>  so that W can be represented as
*>
*>               W = [ I ]  <- identity, N-by-N
*>                   [ V ]  <- M-by-N, same form as B.
*>
*>  Thus, all of information needed for W is contained on exit in B, which
*>  we call V above.  Note that V has the same form as B; that is,
*>
*>               V = [ V1 ] <- (M-L)-by-N rectangular
*>                   [ V2 ] <-     L-by-N upper trapezoidal.
*>
*>  The columns of V represent the vectors which define the H(i)'s.
*>
*>  The number of blocks is B = ceiling(N/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-N matrix T as
*>
*>               T = [T1 T2 ... TB].
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE STPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*     ..
*     .. Array Arguments ..
      REAL A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, LB, MB, IINFO
*     ..
*     .. External Subroutines ..
      EXTERNAL   STPQRT2, STPRFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 .OR. (L.GT.MIN(M,N) .AND. MIN(M,N).GE.0)) THEN
         INFO = -3
      ELSE IF( NB.LT.1 .OR. (NB.GT.N .AND. N.GT.0)) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STPQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      DO I = 1, N, NB
*
*     Compute the QR factorization of the current block
*
         IB = MIN( N-I+1, NB )
         MB = MIN( M-L+I+IB-1, M )
         IF( I.GE.L ) THEN
            LB = 0
         ELSE
            LB = MB-M+L-I+1
         END IF
*
         CALL STPQRT2( MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB,
     $                 T(1, I ), LDT, IINFO )
*
*     Update by applying H^H to B(:,I+IB:N) from the left
*
         IF( I+IB.LE.N ) THEN
            CALL STPRFB( 'L', 'T', 'F', 'C', MB, N-I-IB+1, IB, LB,
     $                    B( 1, I ), LDB, T( 1, I ), LDT,
     $                    A( I, I+IB ), LDA, B( 1, I+IB ), LDB,
     $                    WORK, IB )
         END IF
      END DO
      RETURN
*
*     End of STPQRT
*
      END
*> \brief \b ZGEMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgemqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgemqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgemqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
*                          C, LDC, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*       ..
*       .. Array Arguments ..
*       COMPLEX*16   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEMQRT overwrites the general complex M-by-N matrix C with
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q C            C Q
*> TRANS = 'C':    Q**H C            C Q**H
*>
*> where Q is a complex orthogonal matrix defined as the product of K
*> elementary reflectors:
*>
*>       Q = H(1) H(2) . . . H(K) = I - V T V**H
*>
*> generated using the compact WY representation as returned by ZGEQRT.
*>
*> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left;
*>          = 'R': apply Q or Q**H from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'C':  Conjugate transpose, apply Q**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CGEQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CGEQRT in the first K columns of its array argument A.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CGEQRT, stored as a NB-by-N matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q C, Q**H C, C Q**H or C Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array. The dimension of WORK is
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      SUBROUTINE ZGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
     $                   C, LDC, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
*     ..
*     .. Array Arguments ..
      COMPLEX*16   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, LDWORK, KF, Q
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'C' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF( LEFT ) THEN
         LDWORK = MAX( 1, N )
         Q = M
      ELSE IF ( RIGHT ) THEN
         LDWORK = MAX( 1, M )
         Q = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.Q ) THEN
         INFO = -5
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0)) THEN
         INFO = -6
      ELSE IF( LDV.LT.MAX( 1, Q ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL ZLARFB( 'L', 'C', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            CALL ZLARFB( 'R', 'N', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL ZLARFB( 'L', 'N', 'F', 'C', M-I+1, N, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( I, 1 ), LDC, WORK, LDWORK )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            CALL ZLARFB( 'R', 'C', 'F', 'C', M, N-I+1, IB,
     $                   V( I, I ), LDV, T( 1, I ), LDT,
     $                   C( 1, I ), LDC, WORK, LDWORK )
         END DO
*
      END IF
*
      RETURN
*
*     End of ZGEMQRT
*
      END
*> \brief \b ZGEQRFP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEQRFP + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrfp.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrfp.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrfp.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEQR2P computes a QR factorization of a complex M-by-N matrix A:
*>
*>    A = Q * ( R ),
*>            ( 0 )
*>
*> where:
*>
*>    Q is a M-by-M orthogonal matrix;
*>    R is an upper-triangular N-by-N matrix with nonnegative diagonal
*>    entries;
*>    0 is a (M-N)-by-N zero matrix, if M > N.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if m >= n). The diagonal entries of R
*>          are real and nonnegative; The elements below the diagonal,
*>          with the array TAU, represent the unitary matrix Q as a
*>          product of min(m,n) elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16GEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**H
*>
*>  where tau is a complex scalar, and v is a complex vector with
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*>  and tau in TAU(i).
*>
*> See Lapack Working Note 203 for details
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEQR2P, ZLARFB, ZLARFT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQRFP', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'ZGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL ZGEQR2P( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL ZLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H**H to A(i:m,i+ib:n) from the left
*
               CALL ZLARFB( 'Left', 'Conjugate transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL ZGEQR2P( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGEQRFP
*
      END
*> \brief \b ZGEQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDT, M, N, NB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEQRT computes a blocked QR factorization of a complex M-by-N matrix A
*> using the compact WY representation of Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if M >= N); the elements below the diagonal
*>          are the columns of V.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,MIN(M,N))
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See below
*>          for further details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16GEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix V stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal. For example, if M=5 and N=3, the matrix V is
*>
*>               V = (  1       )
*>                   ( v1  1    )
*>                   ( v1 v2  1 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  where the vi's represent the vectors which define H(i), which are returned
*>  in the matrix A.  The 1's along the diagonal of V are not stored in A.
*>
*>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-K matrix T as
*>
*>               T = (T1 T2 ... TB).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDT, M, N, NB
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, IINFO, K
      LOGICAL    USE_RECURSIVE_QR
      PARAMETER( USE_RECURSIVE_QR=.TRUE. )
*     ..
*     .. External Subroutines ..
      EXTERNAL   ZGEQRT2, ZGEQRT3, ZLARFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NB.LT.1 .OR. ( NB.GT.MIN(M,N) .AND. MIN(M,N).GT.0 ) )THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) RETURN
*
*     Blocked loop of length K
*
      DO I = 1, K,  NB
         IB = MIN( K-I+1, NB )
*
*     Compute the QR factorization of the current block A(I:M,I:I+IB-1)
*
         IF( USE_RECURSIVE_QR ) THEN
            CALL ZGEQRT3( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         ELSE
            CALL ZGEQRT2( M-I+1, IB, A(I,I), LDA, T(1,I), LDT, IINFO )
         END IF
         IF( I+IB.LE.N ) THEN
*
*     Update by applying H**H to A(I:M,I+IB:N) from the left
*
            CALL ZLARFB( 'L', 'C', 'F', 'C', M-I+1, N-I-IB+1, IB,
     $                   A( I, I ), LDA, T( 1, I ), LDT,
     $                   A( I, I+IB ), LDA, WORK , N-I-IB+1 )
         END IF
      END DO
      RETURN
*
*     End of ZGEQRT
*
      END
*> \brief \b ZLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAHQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
*                          IHIZ, Z, LDZ, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAHQR is an auxiliary routine called by CHSEQR to update the
*>    eigenvalues and Schur decomposition already computed by CHSEQR, by
*>    dealing with the Hessenberg submatrix in rows and columns ILO to
*>    IHI.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          = .TRUE. : the full Schur form T is required;
*>          = .FALSE.: only eigenvalues are required.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE. : the matrix of Schur vectors Z is required;
*>          = .FALSE.: Schur vectors are not required.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix H.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          It is assumed that H is already upper triangular in rows and
*>          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
*>          ZLAHQR works primarily with the Hessenberg submatrix in rows
*>          and columns ILO to IHI, but applies transformations to all of
*>          H if WANTT is .TRUE..
*>          1 <= ILO <= max(1,IHI); IHI <= N.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>          On entry, the upper Hessenberg matrix H.
*>          On exit, if INFO is zero and if WANTT is .TRUE., then H
*>          is upper triangular in rows and columns ILO:IHI.  If INFO
*>          is zero and if WANTT is .FALSE., then the contents of H
*>          are unspecified on exit.  The output state of H in case
*>          INF is positive is below under the description of INFO.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          The leading dimension of the array H. LDH >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (N)
*>          The computed eigenvalues ILO to IHI are stored in the
*>          corresponding elements of W. If WANTT is .TRUE., the
*>          eigenvalues are stored in the same order as on the diagonal
*>          of the Schur form returned in H, with W(i) = H(i,i).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>          Specify the rows of Z to which transformations must be
*>          applied if WANTZ is .TRUE..
*>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>          If WANTZ is .TRUE., on entry Z must contain the current
*>          matrix Z of transformations accumulated by CHSEQR, and on
*>          exit Z has been updated; transformations are applied only to
*>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*>          If WANTZ is .FALSE., Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           = 0:   successful exit
*>           > 0:   if INFO = i, ZLAHQR failed to compute all the
*>                  eigenvalues ILO to IHI in a total of 30 iterations
*>                  per eigenvalue; elements i+1:ihi of W contain
*>                  those eigenvalues which have been successfully
*>                  computed.
*>
*>                  If INFO > 0 and WANTT is .FALSE., then on exit,
*>                  the remaining unconverged eigenvalues are the
*>                  eigenvalues of the upper Hessenberg matrix
*>                  rows and columns ILO through INFO of the final,
*>                  output value of H.
*>
*>                  If INFO > 0 and WANTT is .TRUE., then on exit
*>          (*)       (initial value of H)*U  = U*(final value of H)
*>                  where U is an orthogonal matrix.    The final
*>                  value of H is upper Hessenberg and triangular in
*>                  rows and columns INFO+1 through IHI.
*>
*>                  If INFO > 0 and WANTZ is .TRUE., then on exit
*>                      (final value of Z)  = (initial value of Z)*U
*>                  where U is the orthogonal matrix in (*)
*>                  (regardless of the value of WANTT.)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>     02-96 Based on modifications by
*>     David Day, Sandia National Laboratory, USA
*>
*>     12-04 Further modifications by
*>     Ralph Byers, University of Kansas, USA
*>     This is a modified version of ZLAHQR from LAPACK version 3.0.
*>     It is (1) more robust against overflow and underflow and
*>     (2) adopts the more conservative Ahues & Tisseur stopping
*>     criterion (LAWN 122, 1997).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
*     ..
*
*  =========================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ),
     $                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE, HALF
      PARAMETER          ( RZERO = 0.0d0, RONE = 1.0d0, HALF = 0.5d0 )
      DOUBLE PRECISION   DAT1
      PARAMETER          ( DAT1 = 3.0d0 / 4.0d0 )
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 10 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U,
     $                   V2, X, Y
      DOUBLE PRECISION   AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX,
     $                   SAFMIN, SMLNUM, SX, T2, TST, ULP
      INTEGER            I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M,
     $                   NH, NZ, KDEFL
*     ..
*     .. Local Arrays ..
      COMPLEX*16         V( 2 )
*     ..
*     .. External Functions ..
      COMPLEX*16         ZLADIV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           ZLADIV, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, ZCOPY, ZLARFG, ZSCAL
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SQRT
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         W( ILO ) = H( ILO, ILO )
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 )
     $   H( IHI, IHI-2 ) = ZERO
*     ==== ensure that subdiagonal entries are real ====
      IF( WANTT ) THEN
         JLO = 1
         JHI = N
      ELSE
         JLO = ILO
         JHI = IHI
      END IF
      DO 20 I = ILO + 1, IHI
         IF( DIMAG( H( I, I-1 ) ).NE.RZERO ) THEN
*           ==== The following redundant normalization
*           .    avoids problems with both gradual and
*           .    sudden underflow in ABS(H(I,I-1)) ====
            SC = H( I, I-1 ) / CABS1( H( I, I-1 ) )
            SC = DCONJG( SC ) / ABS( SC )
            H( I, I-1 ) = ABS( H( I, I-1 ) )
            CALL ZSCAL( JHI-I+1, SC, H( I, I ), LDH )
            CALL ZSCAL( MIN( JHI, I+1 )-JLO+1, DCONJG( SC ),
     $                  H( JLO, I ), 1 )
            IF( WANTZ )
     $         CALL ZSCAL( IHIZ-ILOZ+1, DCONJG( SC ), Z( ILOZ, I ), 1 )
         END IF
   20 CONTINUE
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( NH ) / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITMAX is the total number of QR iterations allowed.
*
      ITMAX = 30 * MAX( 10, NH )
*
*     KDEFL counts the number of iterations since a deflation
*
      KDEFL = 0
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   30 CONTINUE
      IF( I.LT.ILO )
     $   GO TO 150
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      L = ILO
      DO 130 ITS = 0, ITMAX
*
*        Look for a single small subdiagonal element.
*
         DO 40 K = I, L + 1, -1
            IF( CABS1( H( K, K-1 ) ).LE.SMLNUM )
     $         GO TO 50
            TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO )
     $            TST = TST + ABS( DBLE( H( K-1, K-2 ) ) )
               IF( K+1.LE.IHI )
     $            TST = TST + ABS( DBLE( H( K+1, K ) ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some examples.  ====
            IF( ABS( DBLE( H( K, K-1 ) ) ).LE.ULP*TST ) THEN
               AB = MAX( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               BA = MIN( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
               AA = MAX( CABS1( H( K, K ) ),
     $              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( CABS1( H( K, K ) ),
     $              CABS1( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM,
     $             ULP*( BB*( AA / S ) ) ) )GO TO 50
            END IF
   40    CONTINUE
   50    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 has split off.
*
         IF( L.GE.I )
     $      GO TO 140
         KDEFL = KDEFL + 1
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( MOD(KDEFL,2*KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = DAT1*ABS( DBLE( H( I, I-1 ) ) )
            T = S + H( I, I )
         ELSE IF( MOD(KDEFL,KEXSH).EQ.0 ) THEN
*
*           Exceptional shift.
*
            S = DAT1*ABS( DBLE( H( L+1, L ) ) )
            T = S + H( L, L )
         ELSE
*
*           Wilkinson's shift.
*
            T = H( I, I )
            U = SQRT( H( I-1, I ) )*SQRT( H( I, I-1 ) )
            S = CABS1( U )
            IF( S.NE.RZERO ) THEN
               X = HALF*( H( I-1, I-1 )-T )
               SX = CABS1( X )
               S = MAX( S, CABS1( X ) )
               Y = S*SQRT( ( X / S )**2+( U / S )**2 )
               IF( SX.GT.RZERO ) THEN
                  IF( DBLE( X / SX )*DBLE( Y )+DIMAG( X / SX )*
     $                DIMAG( Y ).LT.RZERO )Y = -Y
               END IF
               T = T - U*ZLADIV( U, ( X+Y ) )
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 60 M = I - 1, L + 1, -1
*
*           Determine the effect of starting the single-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.
*
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H11S = H11 - T
            H21 = DBLE( H( M+1, M ) )
            S = CABS1( H11S ) + ABS( H21 )
            H11S = H11S / S
            H21 = H21 / S
            V( 1 ) = H11S
            V( 2 ) = H21
            H10 = DBLE( H( M, M-1 ) )
            IF( ABS( H10 )*ABS( H21 ).LE.ULP*
     $          ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) )
     $          GO TO 70
   60    CONTINUE
         H11 = H( L, L )
         H22 = H( L+1, L+1 )
         H11S = H11 - T
         H21 = DBLE( H( L+1, L ) )
         S = CABS1( H11S ) + ABS( H21 )
         H11S = H11S / S
         H21 = H21 / S
         V( 1 ) = H11S
         V( 2 ) = H21
   70    CONTINUE
*
*        Single-shift QR step
*
         DO 120 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix.
*
*           V(2) is always real before the call to ZLARFG, and hence
*           after the call T2 ( = T1*V(2) ) is also real.
*
            IF( K.GT.M )
     $         CALL ZCOPY( 2, H( K, K-1 ), 1, V, 1 )
            CALL ZLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
            END IF
            V2 = V( 2 )
            T2 = DBLE( T1*V2 )
*
*           Apply G from the left to transform the rows of the matrix
*           in columns K to I2.
*
            DO 80 J = K, I2
               SUM = DCONJG( T1 )*H( K, J ) + T2*H( K+1, J )
               H( K, J ) = H( K, J ) - SUM
               H( K+1, J ) = H( K+1, J ) - SUM*V2
   80       CONTINUE
*
*           Apply G from the right to transform the columns of the
*           matrix in rows I1 to min(K+2,I).
*
            DO 90 J = I1, MIN( K+2, I )
               SUM = T1*H( J, K ) + T2*H( J, K+1 )
               H( J, K ) = H( J, K ) - SUM
               H( J, K+1 ) = H( J, K+1 ) - SUM*DCONJG( V2 )
   90       CONTINUE
*
            IF( WANTZ ) THEN
*
*              Accumulate transformations in the matrix Z
*
               DO 100 J = ILOZ, IHIZ
                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                  Z( J, K ) = Z( J, K ) - SUM
                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*DCONJG( V2 )
  100          CONTINUE
            END IF
*
            IF( K.EQ.M .AND. M.GT.L ) THEN
*
*              If the QR step was started at row M > L because two
*              consecutive small subdiagonals were found, then extra
*              scaling must be performed to ensure that H(M,M-1) remains
*              real.
*
               TEMP = ONE - T1
               TEMP = TEMP / ABS( TEMP )
               H( M+1, M ) = H( M+1, M )*DCONJG( TEMP )
               IF( M+2.LE.I )
     $            H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
               DO 110 J = M, I
                  IF( J.NE.M+1 ) THEN
                     IF( I2.GT.J )
     $                  CALL ZSCAL( I2-J, TEMP, H( J, J+1 ), LDH )
                     CALL ZSCAL( J-I1, DCONJG( TEMP ), H( I1, J ), 1 )
                     IF( WANTZ ) THEN
                        CALL ZSCAL( NZ, DCONJG( TEMP ), Z( ILOZ, J ),
     $                              1 )
                     END IF
                  END IF
  110          CONTINUE
            END IF
  120    CONTINUE
*
*        Ensure that H(I,I-1) is real.
*
         TEMP = H( I, I-1 )
         IF( DIMAG( TEMP ).NE.RZERO ) THEN
            RTEMP = ABS( TEMP )
            H( I, I-1 ) = RTEMP
            TEMP = TEMP / RTEMP
            IF( I2.GT.I )
     $         CALL ZSCAL( I2-I, DCONJG( TEMP ), H( I, I+1 ), LDH )
            CALL ZSCAL( I-I1, TEMP, H( I1, I ), 1 )
            IF( WANTZ ) THEN
               CALL ZSCAL( NZ, TEMP, Z( ILOZ, I ), 1 )
            END IF
         END IF
*
  130 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  140 CONTINUE
*
*     H(I,I-1) is negligible: one eigenvalue has converged.
*
      W( I ) = H( I, I )
*     reset deflation counter
      KDEFL = 0
*
*     return to start of the main loop with new value of I.
*
      I = L - 1
      GO TO 30
*
  150 CONTINUE
      RETURN
*
*     End of ZLAHQR
*
      END
*> \brief \b ZSYCONV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZSYCONV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYCONV converts A given by ZHETRF into L and D or vice-versa.
*> Get nondiagonal elements of D (returned in workspace) and
*> apply or reverse permutation done in TRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by ZSYTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by ZSYTRF.
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is COMPLEX*16 array, dimension (N)
*>          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1
*>          or 2-by-2 block diagonal matrix D in LDLT.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16SYcomputational
*
*  =====================================================================
      SUBROUTINE ZSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = (0.0D+0,0.0D+0) )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, J
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYCONV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*           Convert VALUE
*
            I=N
            E(1)=ZERO
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I-1,I)
                  E(I-1)=ZERO
                  A(I-1,I)=ZERO
                  I=I-1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I-1
            END DO
*
*           Convert PERMUTATIONS
*
            I=N
            DO WHILE ( I .GE. 1 )
               IF( IPIV(I) .GT. 0) THEN
                  IP=IPIV(I)
                  IF( I .LT. N) THEN
                     DO 12 J= I+1,N
                       TEMP=A(IP,J)
                       A(IP,J)=A(I,J)
                       A(I,J)=TEMP
 12                  CONTINUE
                  ENDIF
               ELSE
                  IP=-IPIV(I)
                  IF( I .LT. N) THEN
                     DO 13 J= I+1,N
                        TEMP=A(IP,J)
                        A(IP,J)=A(I-1,J)
                        A(I-1,J)=TEMP
 13                  CONTINUE
                  ENDIF
                  I=I-1
               ENDIF
               I=I-1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*           Revert PERMUTATIONS
*
            I=1
            DO WHILE ( I .LE. N )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF( I .LT. N) THEN
                  DO J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
                  END DO
                  ENDIF
               ELSE
                 IP=-IPIV(I)
                 I=I+1
                 IF( I .LT. N) THEN
                    DO J= I+1,N
                       TEMP=A(IP,J)
                       A(IP,J)=A(I-1,J)
                       A(I-1,J)=TEMP
                    END DO
                 ENDIF
               ENDIF
               I=I+1
            END DO
*
*           Revert VALUE
*
            I=N
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I-1,I)=E(I)
                  I=I-1
               ENDIF
               I=I-1
            END DO
         END IF
*
      ELSE
*
*        A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*           Convert VALUE
*
            I=1
            E(N)=ZERO
            DO WHILE ( I .LE. N )
               IF( I.LT.N .AND. IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I+1,I)
                  E(I+1)=ZERO
                  A(I+1,I)=ZERO
                  I=I+1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I+1
            END DO
*
*           Convert PERMUTATIONS
*
            I=1
            DO WHILE ( I .LE. N )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF (I .GT. 1) THEN
                     DO 22 J= 1,I-1
                        TEMP=A(IP,J)
                        A(IP,J)=A(I,J)
                        A(I,J)=TEMP
 22                  CONTINUE
                  ENDIF
               ELSE
                  IP=-IPIV(I)
                  IF (I .GT. 1) THEN
                     DO 23 J= 1,I-1
                        TEMP=A(IP,J)
                        A(IP,J)=A(I+1,J)
                        A(I+1,J)=TEMP
 23                  CONTINUE
                  ENDIF
                  I=I+1
               ENDIF
               I=I+1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*           Revert PERMUTATIONS
*
            I=N
            DO WHILE ( I .GE. 1 )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I,J)
                        A(I,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ELSE
                  IP=-IPIV(I)
                  I=I-1
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I+1,J)
                        A(I+1,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ENDIF
               I=I-1
            END DO
*
*           Revert VALUE
*
            I=1
            DO WHILE ( I .LE. N-1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I+1,I)=E(I)
                  I=I+1
               ENDIF
               I=I+1
            END DO
         END IF
      END IF
*
      RETURN
*
*     End of ZSYCONV
*
      END
*> \brief \b ZSYCONVF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZSYCONVF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconvf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconvf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconvf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> ZSYCONVF converts the factorization output format used in
*> ZSYTRF provided on entry in parameter A into the factorization
*> output format used in ZSYTRF_RK (or ZSYTRF_BK) that is stored
*> on exit in parameters A and E. It also converts in place details of
*> the intechanges stored in IPIV from the format used in ZSYTRF into
*> the format used in ZSYTRF_RK (or ZSYTRF_BK).
*>
*> If parameter WAY = 'R':
*> ZSYCONVF performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in ZSYTRF_RK
*> (or ZSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in ZSYTRF that is stored
*> on exit in parameter A. It also converts in place details of
*> the intechanges stored in IPIV from the format used in ZSYTRF_RK
*> (or ZSYTRF_BK) into the format used in ZSYTRF.
*>
*> ZSYCONVF can also convert in Hermitian matrix case, i.e. between
*> formats used in ZHETRF and ZHETRF_RK (or ZHETRF_BK).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          ZSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          ZSYTRF_RK or ZSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          ZSYTRF_RK or ZSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          ZSYTRF:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is COMPLEX*16 array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in,out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>
*>          1) If WAY ='C':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in ZSYTRF.
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in ZSYTRF_RK
*>          ( or ZSYTRF_BK).
*>
*>          1) If WAY ='R':
*>          On entry, details of the interchanges and the block
*>          structure of D in the format used in ZSYTRF_RK
*>          ( or ZSYTRF_BK).
*>          On exit, details of the interchanges and the block
*>          structure of D in the format used in ZSYTRF.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16SYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE ZSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           ZSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYCONVF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL ZSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.(I-1) ) THEN
                        CALL ZSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i-1 and IPIV(i-1),
*                 so this should be recorded in two consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I-1 )
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where k increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL ZSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is no interchnge of rows i and and IPIV(i),
*                 so this should be reflected in IPIV format for
*                 *SYTRF_RK ( or *SYTRF_BK)
*
                  IPIV( I ) = I
*
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS and IPIV
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.(I+1) ) THEN
                        CALL ZSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                  END IF
*
*                 Convert IPIV
*                 There is one interchange of rows i+1 and IPIV(i+1),
*                 so this should be recorded in consecutive entries
*                 in IPIV format for *SYTRF
*
                  IPIV( I ) = IPIV( I+1 )
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of ZSYCONVF
*
      END
*> \brief \b ZSYCONVF_ROOK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZSYCONVF_ROOK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconvf_rook.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconvf_rook.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconvf_rook.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> If parameter WAY = 'C':
*> ZSYCONVF_ROOK converts the factorization output format used in
*> ZSYTRF_ROOK provided on entry in parameter A into the factorization
*> output format used in ZSYTRF_RK (or ZSYTRF_BK) that is stored
*> on exit in parameters A and E. IPIV format for ZSYTRF_ROOK and
*> ZSYTRF_RK (or ZSYTRF_BK) is the same and is not converted.
*>
*> If parameter WAY = 'R':
*> ZSYCONVF_ROOK performs the conversion in reverse direction, i.e.
*> converts the factorization output format used in ZSYTRF_RK
*> (or ZSYTRF_BK) provided on entry in parameters A and E into
*> the factorization output format used in ZSYTRF_ROOK that is stored
*> on exit in parameter A. IPIV format for ZSYTRF_ROOK and
*> ZSYTRF_RK (or ZSYTRF_BK) is the same and is not converted.
*>
*> ZSYCONVF_ROOK can also convert in Hermitian matrix case, i.e. between
*> formats used in ZHETRF_ROOK and ZHETRF_RK (or ZHETRF_BK).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix A.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, contains factorization details in format used in
*>          ZSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          ZSYTRF_RK or ZSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains factorization details in format used in
*>          ZSYTRF_RK or ZSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                are stored on exit in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, contains factorization details in format used in
*>          ZSYTRF_ROOK:
*>            a) all elements of the symmetric block diagonal
*>               matrix D on the diagonal of A and on superdiagonal
*>               (or subdiagonal) of A, and
*>            b) If UPLO = 'U': multipliers used to obtain factor U
*>               in the superdiagonal part of A.
*>               If UPLO = 'L': multipliers used to obtain factor L
*>               in the superdiagonal part of A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is COMPLEX*16 array, dimension (N)
*>
*>          1) If WAY ='C':
*>
*>          On entry, just a workspace.
*>
*>          On exit, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
*>
*>          2) If WAY = 'R':
*>
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          On exit, is not changed
*> \endverbatim
*.
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          On entry, details of the interchanges and the block
*>          structure of D as determined:
*>          1) by ZSYTRF_ROOK, if WAY ='C';
*>          2) by ZSYTRF_RK (or ZSYTRF_BK), if WAY ='R'.
*>          The IPIV format is the same for all these routines.
*>
*>          On exit, is not changed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16SYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*  =====================================================================
      SUBROUTINE ZSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           ZSWAP, XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, IP2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYCONVF_ROOK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Begin A is UPPER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is upper)
*
*
*           Convert VALUE
*
*           Assign superdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = N
            E( 1 ) = ZERO
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  E( I ) = A( I-1, I )
                  E( I-1 ) = ZERO
                  A( I-1, I ) = ZERO
                  I = I - 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I - 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1)
*                 in A(1:i,N-i:N)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( N-I, A( I, I+1 ), LDA,
     $                              A( IP, I+1 ), LDA )
                     END IF
                     IF( IP2.NE.(I-1) ) THEN
                        CALL ZSWAP( N-I, A( I-1, I+1 ), LDA,
     $                              A( IP2, I+1 ), LDA )
                     END IF
                  END IF
                  I = I - 1
*
               END IF
               I = I - 1
            END DO
*
         ELSE
*
*           Revert A (A is upper)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of upper part of A
*           in reverse factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
*
                  IP = IPIV( I )
                  IF( I.LT.N ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i)
*                 in A(1:i,N-i:N)
*
                  I = I + 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I-1 )
                  IF( I.LT.N ) THEN
                     IF( IP2.NE.(I-1) ) THEN
                        CALL ZSWAP( N-I, A( IP2, I+1 ), LDA,
     $                              A( I-1, I+1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( N-I, A( IP, I+1 ), LDA,
     $                              A( I, I+1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I + 1
            END DO
*
*           Revert VALUE
*           Assign superdiagonal entries of D from array E to
*           superdiagonal entries of A.
*
            I = N
            DO WHILE ( I.GT.1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I-1, I ) = E( I )
                  I = I - 1
               END IF
               I = I - 1
            END DO
*
*        End A is UPPER
*
         END IF
*
      ELSE
*
*        Begin A is LOWER
*
         IF ( CONVERT ) THEN
*
*           Convert A (A is lower)
*
*
*           Convert VALUE
*           Assign subdiagonal entries of D to array E and zero out
*           corresponding entries in input storage A
*
            I = 1
            E( N ) = ZERO
            DO WHILE ( I.LE.N )
               IF( I.LT.N .AND. IPIV(I).LT.0 ) THEN
                  E( I ) = A( I+1, I )
                  E( I+1 ) = ZERO
                  A( I+1, I ) = ZERO
                  I = I + 1
               ELSE
                  E( I ) = ZERO
               END IF
               I = I + 1
            END DO
*
*           Convert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in factorization order where i increases from 1 to N
*
            I = 1
            DO WHILE ( I.LE.N )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1)
*                 in A(i:N,1:i-1)
*
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( I-1, A( I, 1 ), LDA,
     $                              A( IP, 1 ), LDA )
                     END IF
                     IF( IP2.NE.(I+1) ) THEN
                        CALL ZSWAP( I-1, A( I+1, 1 ), LDA,
     $                              A( IP2, 1 ), LDA )
                     END IF
                  END IF
                  I = I + 1
*
               END IF
               I = I + 1
            END DO
*
         ELSE
*
*           Revert A (A is lower)
*
*
*           Revert PERMUTATIONS
*
*           Apply permutations to submatrices of lower part of A
*           in reverse factorization order where i decreases from N to 1
*
            I = N
            DO WHILE ( I.GE.1 )
               IF( IPIV( I ).GT.0 ) THEN
*
*                 1-by-1 pivot interchange
*
*                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
*
                  IP = IPIV( I )
                  IF ( I.GT.1 ) THEN
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               ELSE
*
*                 2-by-2 pivot interchange
*
*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i)
*                 in A(i:N,1:i-1)
*
                  I = I - 1
                  IP = -IPIV( I )
                  IP2 = -IPIV( I+1 )
                  IF ( I.GT.1 ) THEN
                     IF( IP2.NE.(I+1) ) THEN
                        CALL ZSWAP( I-1, A( IP2, 1 ), LDA,
     $                              A( I+1, 1 ), LDA )
                     END IF
                     IF( IP.NE.I ) THEN
                        CALL ZSWAP( I-1, A( IP, 1 ), LDA,
     $                              A( I, 1 ), LDA )
                     END IF
                  END IF
*
               END IF
               I = I - 1
            END DO
*
*           Revert VALUE
*           Assign subdiagonal entries of D from array E to
*           subgiagonal entries of A.
*
            I = 1
            DO WHILE ( I.LE.N-1 )
               IF( IPIV( I ).LT.0 ) THEN
                  A( I + 1, I ) = E( I )
                  I = I + 1
               END IF
               I = I + 1
            END DO
*
         END IF
*
*        End A is LOWER
*
      END IF

      RETURN
*
*     End of ZSYCONVF_ROOK
*
      END
*> \brief \b ZTPMQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZTPMQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpmqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpmqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpmqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
*                           A, LDA, B, LDB, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER SIDE, TRANS
*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*       ..
*       .. Array Arguments ..
*       COMPLEX*16   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
*      $          WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTPMQRT applies a complex orthogonal matrix Q obtained from a
*> "triangular-pentagonal" complex block reflector H to a general
*> complex matrix C, which consists of two blocks A and B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left;
*>          = 'R': apply Q or Q**H from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'C':  Conjugate transpose, apply Q**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The order of the trapezoidal part of V.
*>          K >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size used for the storage of T.  K >= NB >= 1.
*>          This must be the same value of NB used to generate T
*>          in CTPQRT.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          CTPQRT in B.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If SIDE = 'L', LDV >= max(1,M);
*>          if SIDE = 'R', LDV >= max(1,N).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,K)
*>          The upper triangular factors of the block reflectors
*>          as returned by CTPQRT, stored as a NB-by-K matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension
*>          (LDA,N) if SIDE = 'L' or
*>          (LDA,K) if SIDE = 'R'
*>          On entry, the K-by-N or M-by-K matrix A.
*>          On exit, A is overwritten by the corresponding block of
*>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDC >= max(1,K);
*>          If SIDE = 'R', LDC >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On entry, the M-by-N matrix B.
*>          On exit, B is overwritten by the corresponding block of
*>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.
*>          LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array. The dimension of WORK is
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16OTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The columns of the pentagonal matrix V contain the elementary reflectors
*>  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
*>  trapezoidal block V2:
*>
*>        V = [V1]
*>            [V2].
*>
*>  The size of the trapezoidal block V2 is determined by the parameter L,
*>  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L
*>  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular;
*>  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
*>
*>  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K.
*>                      [B]
*>
*>  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K.
*>
*>  The complex orthogonal matrix Q is formed from V and T.
*>
*>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.
*>
*>  If TRANS='C' and SIDE='L', C is on exit replaced with Q**H * C.
*>
*>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.
*>
*>  If TRANS='C' and SIDE='R', C is on exit replaced with C * Q**H.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
     $                    A, LDA, B, LDB, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER SIDE, TRANS
      INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
*     ..
*     .. Array Arguments ..
      COMPLEX*16   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
     $          WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN
      INTEGER            I, IB, MB, LB, KF, LDAQ, LDVQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZTPRFB, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     .. Test the input arguments ..
*
      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'C' )
      NOTRAN = LSAME( TRANS, 'N' )
*
      IF ( LEFT ) THEN
         LDVQ = MAX( 1, M )
         LDAQ = MAX( 1, K )
      ELSE IF ( RIGHT ) THEN
         LDVQ = MAX( 1, N )
         LDAQ = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      ELSE IF( L.LT.0 .OR. L.GT.K ) THEN
         INFO = -6
      ELSE IF( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0) ) THEN
         INFO = -7
      ELSE IF( LDV.LT.LDVQ ) THEN
         INFO = -9
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -11
      ELSE IF( LDA.LT.LDAQ ) THEN
         INFO = -13
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTPMQRT', -INFO )
         RETURN
      END IF
*
*     .. Quick return if possible ..
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
*
      IF( LEFT .AND. TRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL ZTPRFB( 'L', 'C', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. NOTRAN ) THEN
*
         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL ZTPRFB( 'R', 'N', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      ELSE IF( LEFT .AND. NOTRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-M+L-I+1
            END IF
            CALL ZTPRFB( 'L', 'N', 'F', 'C', MB, N, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( I, 1 ), LDA, B, LDB, WORK, IB )
         END DO
*
      ELSE IF( RIGHT .AND. TRAN ) THEN
*
         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            IF( I.GE.L ) THEN
               LB = 0
            ELSE
               LB = MB-N+L-I+1
            END IF
            CALL ZTPRFB( 'R', 'C', 'F', 'C', M, MB, IB, LB,
     $                   V( 1, I ), LDV, T( 1, I ), LDT,
     $                   A( 1, I ), LDA, B, LDB, WORK, M )
         END DO
*
      END IF
*
      RETURN
*
*     End of ZTPMQRT
*
      END
*> \brief \b ZTPQRT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZTPQRT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpqrt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpqrt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpqrt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTPQRT computes a blocked QR factorization of a complex
*> "triangular-pentagonal" matrix C, which is composed of a
*> triangular block A and pentagonal block B, using the compact
*> WY representation for Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix B.
*>          M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix B, and the order of the
*>          triangular matrix A.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>          The number of rows of the upper trapezoidal part of B.
*>          MIN(M,N) >= L >= 0.  See Further Details.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The block size to be used in the blocked QR.  N >= NB >= 1.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the upper triangular N-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the upper triangular matrix R.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
*>          are rectangular, and the last L rows are upper trapezoidal.
*>          On exit, B contains the pentagonal matrix V.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,N)
*>          The upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See Further Details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (NB*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16OTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The input matrix C is a (N+M)-by-N matrix
*>
*>               C = [ A ]
*>                   [ B ]
*>
*>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
*>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
*>  upper trapezoidal matrix B2:
*>
*>               B = [ B1 ]  <- (M-L)-by-N rectangular
*>                   [ B2 ]  <-     L-by-N upper trapezoidal.
*>
*>  The upper trapezoidal matrix B2 consists of the first L rows of a
*>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
*>  B is rectangular M-by-N; if M=L=N, B is upper triangular.
*>
*>  The matrix W stores the elementary reflectors H(i) in the i-th column
*>  below the diagonal (of A) in the (N+M)-by-N input matrix C
*>
*>               C = [ A ]  <- upper triangular N-by-N
*>                   [ B ]  <- M-by-N pentagonal
*>
*>  so that W can be represented as
*>
*>               W = [ I ]  <- identity, N-by-N
*>                   [ V ]  <- M-by-N, same form as B.
*>
*>  Thus, all of information needed for W is contained on exit in B, which
*>  we call V above.  Note that V has the same form as B; that is,
*>
*>               V = [ V1 ] <- (M-L)-by-N rectangular
*>                   [ V2 ] <-     L-by-N upper trapezoidal.
*>
*>  The columns of V represent the vectors which define the H(i)'s.
*>
*>  The number of blocks is B = ceiling(N/NB), where each
*>  block is of order NB except for the last block, which is of order
*>  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB
*>  for the last block) T's are stored in the NB-by-N matrix T as
*>
*>               T = [T1 T2 ... TB].
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER    I, IB, LB, MB, IINFO
*     ..
*     .. External Subroutines ..
      EXTERNAL   ZTPQRT2, ZTPRFB, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 .OR. (L.GT.MIN(M,N) .AND. MIN(M,N).GE.0)) THEN
         INFO = -3
      ELSE IF( NB.LT.1 .OR. (NB.GT.N .AND. N.GT.0)) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDT.LT.NB ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTPQRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      DO I = 1, N, NB
*
*     Compute the QR factorization of the current block
*
         IB = MIN( N-I+1, NB )
         MB = MIN( M-L+I+IB-1, M )
         IF( I.GE.L ) THEN
            LB = 0
         ELSE
            LB = MB-M+L-I+1
         END IF
*
         CALL ZTPQRT2( MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB,
     $                 T(1, I ), LDT, IINFO )
*
*     Update by applying H**H to B(:,I+IB:N) from the left
*
         IF( I+IB.LE.N ) THEN
            CALL ZTPRFB( 'L', 'C', 'F', 'C', MB, N-I-IB+1, IB, LB,
     $                    B( 1, I ), LDB, T( 1, I ), LDT,
     $                    A( I, I+IB ), LDA, B( 1, I+IB ), LDB,
     $                    WORK, IB )
         END IF
      END DO
      RETURN
*
*     End of ZTPQRT
*
      END
