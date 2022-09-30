/* lapack_extras.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b CGEMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgemqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgemqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgemqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEMQRT overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q C            C Q */
/* > TRANS = 'C':    Q**H C            C Q**H */
/* > */
/* > where Q is a complex orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* >       Q = H(1) H(2) . . . H(K) = I - V T V**H */
/* > */
/* > generated using the compact WY representation as returned by CGEQRT. */
/* > */
/* > Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Conjugate transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CGEQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGEQRT in the first K columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CGEQRT, stored as a NB-by-N matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q C, Q**H C, C Q**H or C Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgemqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *nb, complex *v, integer *ldv, complex *t, 
	integer *ldt, complex *c__, integer *ldc, complex *work, integer *
	info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, q, ib, kf;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int clarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, complex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static integer ldwork;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldwork = max(1,*n);
	q = *m;
    } else if (right) {
	ldwork = max(1,*m);
	q = *n;
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > q) {
	*info = -5;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -6;
    } else if (*ldv < max(1,q)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    } else if (*ldc < max(1,*m)) {
	*info = -12;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGEMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *m - i__ + 1;
	    clarfb_("L", "C", "F", "C", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *n - i__ + 1;
	    clarfb_("R", "N", "F", "C", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *m - i__ + 1;
	    clarfb_("L", "N", "F", "C", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *n - i__ + 1;
	    clarfb_("R", "C", "F", "C", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    }

    return 0;

/*     End of CGEMQRT */

} /* cgemqrt_ */

/* > \brief \b CGEQRFP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrfp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrfp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrfp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEQR2P computes a QR factorization of a complex M-by-N matrix A: */
/* > */
/* >    A = Q * ( R ), */
/* >            ( 0 ) */
/* > */
/* > where: */
/* > */
/* >    Q is a M-by-M orthogonal matrix; */
/* >    R is an upper-triangular N-by-N matrix with nonnegative diagonal */
/* >    entries; */
/* >    0 is a (M-N)-by-N zero matrix, if M > N. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n). The diagonal entries of R */
/* >          are real and nonnegative; the elements below the diagonal, */
/* >          with the array TAU, represent the unitary matrix Q as a */
/* >          product of min(m,n) elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,N). */
/* >          For optimum performance LWORK >= N*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > */
/* > See Lapack Working Note 203 for details */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgeqrfp_(integer *m, integer *n, complex *a, integer *
	lda, complex *tau, complex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int clarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, complex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), clarft_(char *, char *, integer *, 
	    integer *, complex *, integer *, complex *, complex *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int cgeqr2p_(integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
    lwkopt = *n * nb;
    work[1].r = (real) lwkopt, work[1].i = 0.f;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGEQRFP", &i__1, (ftnlen)7);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	work[1].r = 1.f, work[1].i = 0.f;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "CGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "CGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

	    i__3 = *m - i__ + 1;
	    cgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *m - i__ + 1;
		clarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**H to A(i:m,i+ib:n) from the left */

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		clarfb_("Left", "Conjugate transpose", "Forward", "Columnwise"
			, &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &
			work[1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, 
			&work[ib + 1], &ldwork, (ftnlen)4, (ftnlen)19, (
			ftnlen)7, (ftnlen)10);
	    }
/* L10: */
	}
    } else {
	i__ = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	cgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		1], &iinfo);
    }

    work[1].r = (real) iws, work[1].i = 0.f;
    return 0;

/*     End of CGEQRFP */

} /* cgeqrfp_ */

/* > \brief \b CGEQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDT, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEQRT computes a blocked QR factorization of a complex M-by-N matrix A */
/* > using the compact WY representation of Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if M >= N); the elements below the diagonal */
/* >          are the columns of V. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,MIN(M,N)) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix V stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* >               V = (  1       ) */
/* >                   ( v1  1    ) */
/* >                   ( v1 v2  1 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  where the vi's represent the vectors which define H(i), which are returned */
/* >  in the matrix A.  The 1's along the diagonal of V are not stored in A. */
/* > */
/* >  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-K matrix T as */
/* > */
/* >               T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgeqrt_(integer *m, integer *n, integer *nb, complex *a, 
	integer *lda, complex *t, integer *ldt, complex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ib, iinfo;
    extern /* Subroutine */ int clarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, complex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), 
	    cgeqrt2_(integer *, integer *, complex *, integer *, complex *, 
	    integer *, integer *), cgeqrt3_(integer *, integer *, complex *, 
	    integer *, complex *, integer *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nb < 1 || *nb > min(*m,*n) && min(*m,*n) > 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*ldt < *nb) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGEQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	return 0;
    }

/*     Blocked loop of length K */

    i__1 = k;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	i__3 = k - i__ + 1;
	ib = min(i__3,*nb);

/*     Compute the QR factorization of the current block A(I:M,I:I+IB-1) */

	if (TRUE_) {
	    i__3 = *m - i__ + 1;
	    cgeqrt3_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	} else {
	    i__3 = *m - i__ + 1;
	    cgeqrt2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	}
	if (i__ + ib <= *n) {

/*     Update by applying H**H to A(I:M,I+IB:N) from the left */

	    i__3 = *m - i__ + 1;
	    i__4 = *n - i__ - ib + 1;
	    i__5 = *n - i__ - ib + 1;
	    clarfb_("L", "C", "F", "C", &i__3, &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + 
		    ib) * a_dim1], lda, &work[1], &i__5, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
	}
    }
    return 0;

/*     End of CGEQRT */

} /* cgeqrt_ */

/* > \brief \b CLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using th
e double-shift/single-shift QR algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAHQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, */
/*                          IHIZ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLAHQR is an auxiliary routine called by CHSEQR to update the */
/* >    eigenvalues and Schur decomposition already computed by CHSEQR, by */
/* >    dealing with the Hessenberg submatrix in rows and columns ILO to */
/* >    IHI. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >          = .TRUE. : the full Schur form T is required; */
/* >          = .FALSE.: only eigenvalues are required. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          = .TRUE. : the matrix of Schur vectors Z is required; */
/* >          = .FALSE.: Schur vectors are not required. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          It is assumed that H is already upper triangular in rows and */
/* >          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1). */
/* >          CLAHQR works primarily with the Hessenberg submatrix in rows */
/* >          and columns ILO to IHI, but applies transformations to all of */
/* >          H if WANTT is .TRUE.. */
/* >          1 <= ILO <= max(1,IHI); IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is COMPLEX array, dimension (LDH,N) */
/* >          On entry, the upper Hessenberg matrix H. */
/* >          On exit, if INFO is zero and if WANTT is .TRUE., then H */
/* >          is upper triangular in rows and columns ILO:IHI.  If INFO */
/* >          is zero and if WANTT is .FALSE., then the contents of H */
/* >          are unspecified on exit.  The output state of H in case */
/* >          INF is positive is below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX array, dimension (N) */
/* >          The computed eigenvalues ILO to IHI are stored in the */
/* >          corresponding elements of W. If WANTT is .TRUE., the */
/* >          eigenvalues are stored in the same order as on the diagonal */
/* >          of the Schur form returned in H, with W(i) = H(i,i). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* >          ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* >          IHIZ is INTEGER */
/* >          Specify the rows of Z to which transformations must be */
/* >          applied if WANTZ is .TRUE.. */
/* >          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ,N) */
/* >          If WANTZ is .TRUE., on entry Z must contain the current */
/* >          matrix Z of transformations accumulated by CHSEQR, and on */
/* >          exit Z has been updated; transformations are applied only to */
/* >          the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
/* >          If WANTZ is .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:  successful exit */
/* >           > 0:  if INFO = i, CLAHQR failed to compute all the */
/* >                  eigenvalues ILO to IHI in a total of 30 iterations */
/* >                  per eigenvalue; elements i+1:ihi of W contain */
/* >                  those eigenvalues which have been successfully */
/* >                  computed. */
/* > */
/* >                  If INFO > 0 and WANTT is .FALSE., then on exit, */
/* >                  the remaining unconverged eigenvalues are the */
/* >                  eigenvalues of the upper Hessenberg matrix */
/* >                  rows and columns ILO through INFO of the final, */
/* >                  output value of H. */
/* > */
/* >                  If INFO > 0 and WANTT is .TRUE., then on exit */
/* >          (*)       (initial value of H)*U  = U*(final value of H) */
/* >                  where U is an orthogonal matrix.    The final */
/* >                  value of H is upper Hessenberg and triangular in */
/* >                  rows and columns INFO+1 through IHI. */
/* > */
/* >                  If INFO > 0 and WANTZ is .TRUE., then on exit */
/* >                      (final value of Z)  = (initial value of Z)*U */
/* >                  where U is the orthogonal matrix in (*) */
/* >                  (regardless of the value of WANTT.) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >     02-96 Based on modifications by */
/* >     David Day, Sandia National Laboratory, USA */
/* > */
/* >     12-04 Further modifications by */
/* >     Ralph Byers, University of Kansas, USA */
/* >     This is a modified version of CLAHQR from LAPACK version 3.0. */
/* >     It is (1) more robust against overflow and underflow and */
/* >     (2) adopts the more conservative Ahues & Tisseur stopping */
/* >     criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *w, 
	integer *iloz, integer *ihiz, complex *z__, integer *ldz, integer *
	info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7;

    /* Builtin functions */
    double r_imag(complex *);
    void r_cnjg(complex *, complex *);
    double c_abs(complex *);
    void c_sqrt(complex *, complex *), pow_ci(complex *, complex *, integer *)
	    ;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real s;
    static complex t, u, v[2], x, y;
    static integer i1, i2;
    static complex t1;
    static real t2;
    static complex v2;
    static real aa, ab, ba, bb, h10;
    static complex h11;
    static real h21;
    static complex h22, sc;
    static integer nh, nz;
    static real sx;
    static integer jhi;
    static complex h11s;
    static integer jlo, its;
    static real ulp;
    static complex sum;
    static real tst;
    static complex temp;
    static integer kdefl;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), ccopy_(integer *, complex *, integer *, complex *, 
	    integer *);
    static integer itmax;
    static real rtemp;
    extern /* Subroutine */ int slabad_(real *, real *), clarfg_(integer *, 
	    complex *, complex *, integer *, complex *);
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);
    extern doublereal slamch_(char *, ftnlen);
    static real safmin, safmax, smlnum;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ========================================================= */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }
    if (*ilo == *ihi) {
	i__1 = *ilo;
	i__2 = *ilo + *ilo * h_dim1;
	w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;
	return 0;
    }

/*     ==== clear out the trash ==== */
    i__1 = *ihi - 3;
    for (j = *ilo; j <= i__1; ++j) {
	i__2 = j + 2 + j * h_dim1;
	h__[i__2].r = 0.f, h__[i__2].i = 0.f;
	i__2 = j + 3 + j * h_dim1;
	h__[i__2].r = 0.f, h__[i__2].i = 0.f;
/* L10: */
    }
    if (*ilo <= *ihi - 2) {
	i__1 = *ihi + (*ihi - 2) * h_dim1;
	h__[i__1].r = 0.f, h__[i__1].i = 0.f;
    }
/*     ==== ensure that subdiagonal entries are real ==== */
    if (*wantt) {
	jlo = 1;
	jhi = *n;
    } else {
	jlo = *ilo;
	jhi = *ihi;
    }
    i__1 = *ihi;
    for (i__ = *ilo + 1; i__ <= i__1; ++i__) {
	if (r_imag(&h__[i__ + (i__ - 1) * h_dim1]) != 0.f) {
/*           ==== The following redundant normalization */
/*           .    avoids problems with both gradual and */
/*           .    sudden underflow in ABS(H(I,I-1)) ==== */
	    i__2 = i__ + (i__ - 1) * h_dim1;
	    i__3 = i__ + (i__ - 1) * h_dim1;
	    r__3 = (r__1 = h__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&h__[i__ 
		    + (i__ - 1) * h_dim1]), dabs(r__2));
	    q__1.r = h__[i__2].r / r__3, q__1.i = h__[i__2].i / r__3;
	    sc.r = q__1.r, sc.i = q__1.i;
	    r_cnjg(&q__2, &sc);
	    r__1 = c_abs(&sc);
	    q__1.r = q__2.r / r__1, q__1.i = q__2.i / r__1;
	    sc.r = q__1.r, sc.i = q__1.i;
	    i__2 = i__ + (i__ - 1) * h_dim1;
	    r__1 = c_abs(&h__[i__ + (i__ - 1) * h_dim1]);
	    h__[i__2].r = r__1, h__[i__2].i = 0.f;
	    i__2 = jhi - i__ + 1;
	    cscal_(&i__2, &sc, &h__[i__ + i__ * h_dim1], ldh);
/* Computing MIN */
	    i__3 = jhi, i__4 = i__ + 1;
	    i__2 = min(i__3,i__4) - jlo + 1;
	    r_cnjg(&q__1, &sc);
	    cscal_(&i__2, &q__1, &h__[jlo + i__ * h_dim1], &c__1);
	    if (*wantz) {
		i__2 = *ihiz - *iloz + 1;
		r_cnjg(&q__1, &sc);
		cscal_(&i__2, &q__1, &z__[*iloz + i__ * z_dim1], &c__1);
	    }
	}
/* L20: */
    }

    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion. */

    safmin = slamch_("SAFE MINIMUM", (ftnlen)12);
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    ulp = slamch_("PRECISION", (ftnlen)9);
    smlnum = safmin * ((real) nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

    if (*wantt) {
	i1 = 1;
	i2 = *n;
    }

/*     ITMAX is the total number of QR iterations allowed. */

    itmax = max(10,nh) * 30;

/*     KDEFL counts the number of iterations since a deflation */

    kdefl = 0;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or */
/*     H(L,L-1) is negligible so that the matrix splits. */

    i__ = *ihi;
L30:
    if (i__ < *ilo) {
	goto L150;
    }

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

    l = *ilo;
    i__1 = itmax;
    for (its = 0; its <= i__1; ++its) {

/*        Look for a single small subdiagonal element. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {
	    i__3 = k + (k - 1) * h_dim1;
	    if ((r__1 = h__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&h__[k + (k 
		    - 1) * h_dim1]), dabs(r__2)) <= smlnum) {
		goto L50;
	    }
	    i__3 = k - 1 + (k - 1) * h_dim1;
	    i__4 = k + k * h_dim1;
	    tst = (r__1 = h__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&h__[k - 
		    1 + (k - 1) * h_dim1]), dabs(r__2)) + ((r__3 = h__[i__4]
		    .r, dabs(r__3)) + (r__4 = r_imag(&h__[k + k * h_dim1]), 
		    dabs(r__4)));
	    if (tst == 0.f) {
		if (k - 2 >= *ilo) {
		    i__3 = k - 1 + (k - 2) * h_dim1;
		    tst += (r__1 = h__[i__3].r, dabs(r__1));
		}
		if (k + 1 <= *ihi) {
		    i__3 = k + 1 + k * h_dim1;
		    tst += (r__1 = h__[i__3].r, dabs(r__1));
		}
	    }
/*           ==== The following is a conservative small subdiagonal */
/*           .    deflation criterion due to Ahues & Tisseur (LAWN 122, */
/*           .    1997). It has better mathematical foundation and */
/*           .    improves accuracy in some examples.  ==== */
	    i__3 = k + (k - 1) * h_dim1;
	    if ((r__1 = h__[i__3].r, dabs(r__1)) <= ulp * tst) {
/* Computing MAX */
		i__3 = k + (k - 1) * h_dim1;
		i__4 = k - 1 + k * h_dim1;
		r__5 = (r__1 = h__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&h__[
			k + (k - 1) * h_dim1]), dabs(r__2)), r__6 = (r__3 = 
			h__[i__4].r, dabs(r__3)) + (r__4 = r_imag(&h__[k - 1 
			+ k * h_dim1]), dabs(r__4));
		ab = dmax(r__5,r__6);
/* Computing MIN */
		i__3 = k + (k - 1) * h_dim1;
		i__4 = k - 1 + k * h_dim1;
		r__5 = (r__1 = h__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&h__[
			k + (k - 1) * h_dim1]), dabs(r__2)), r__6 = (r__3 = 
			h__[i__4].r, dabs(r__3)) + (r__4 = r_imag(&h__[k - 1 
			+ k * h_dim1]), dabs(r__4));
		ba = dmin(r__5,r__6);
		i__3 = k - 1 + (k - 1) * h_dim1;
		i__4 = k + k * h_dim1;
		q__2.r = h__[i__3].r - h__[i__4].r, q__2.i = h__[i__3].i - 
			h__[i__4].i;
		q__1.r = q__2.r, q__1.i = q__2.i;
/* Computing MAX */
		i__5 = k + k * h_dim1;
		r__5 = (r__1 = h__[i__5].r, dabs(r__1)) + (r__2 = r_imag(&h__[
			k + k * h_dim1]), dabs(r__2)), r__6 = (r__3 = q__1.r, 
			dabs(r__3)) + (r__4 = r_imag(&q__1), dabs(r__4));
		aa = dmax(r__5,r__6);
		i__3 = k - 1 + (k - 1) * h_dim1;
		i__4 = k + k * h_dim1;
		q__2.r = h__[i__3].r - h__[i__4].r, q__2.i = h__[i__3].i - 
			h__[i__4].i;
		q__1.r = q__2.r, q__1.i = q__2.i;
/* Computing MIN */
		i__5 = k + k * h_dim1;
		r__5 = (r__1 = h__[i__5].r, dabs(r__1)) + (r__2 = r_imag(&h__[
			k + k * h_dim1]), dabs(r__2)), r__6 = (r__3 = q__1.r, 
			dabs(r__3)) + (r__4 = r_imag(&q__1), dabs(r__4));
		bb = dmin(r__5,r__6);
		s = aa + ab;
/* Computing MAX */
		r__1 = smlnum, r__2 = ulp * (bb * (aa / s));
		if (ba * (ab / s) <= dmax(r__1,r__2)) {
		    goto L50;
		}
	    }
/* L40: */
	}
L50:
	l = k;
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

	    i__2 = l + (l - 1) * h_dim1;
	    h__[i__2].r = 0.f, h__[i__2].i = 0.f;
	}

/*        Exit from loop if a submatrix of order 1 has split off. */

	if (l >= i__) {
	    goto L140;
	}
	++kdefl;

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

	if (! (*wantt)) {
	    i1 = l;
	    i2 = i__;
	}

	if (kdefl % 20 == 0) {

/*           Exceptional shift. */

	    i__2 = i__ + (i__ - 1) * h_dim1;
	    s = (r__1 = h__[i__2].r, dabs(r__1)) * .75f;
	    i__2 = i__ + i__ * h_dim1;
	    q__1.r = s + h__[i__2].r, q__1.i = h__[i__2].i;
	    t.r = q__1.r, t.i = q__1.i;
	} else if (kdefl % 10 == 0) {

/*           Exceptional shift. */

	    i__2 = l + 1 + l * h_dim1;
	    s = (r__1 = h__[i__2].r, dabs(r__1)) * .75f;
	    i__2 = l + l * h_dim1;
	    q__1.r = s + h__[i__2].r, q__1.i = h__[i__2].i;
	    t.r = q__1.r, t.i = q__1.i;
	} else {

/*           Wilkinson's shift. */

	    i__2 = i__ + i__ * h_dim1;
	    t.r = h__[i__2].r, t.i = h__[i__2].i;
	    c_sqrt(&q__2, &h__[i__ - 1 + i__ * h_dim1]);
	    c_sqrt(&q__3, &h__[i__ + (i__ - 1) * h_dim1]);
	    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * 
		    q__3.i + q__2.i * q__3.r;
	    u.r = q__1.r, u.i = q__1.i;
	    s = (r__1 = u.r, dabs(r__1)) + (r__2 = r_imag(&u), dabs(r__2));
	    if (s != 0.f) {
		i__2 = i__ - 1 + (i__ - 1) * h_dim1;
		q__2.r = h__[i__2].r - t.r, q__2.i = h__[i__2].i - t.i;
		q__1.r = q__2.r * .5f, q__1.i = q__2.i * .5f;
		x.r = q__1.r, x.i = q__1.i;
		sx = (r__1 = x.r, dabs(r__1)) + (r__2 = r_imag(&x), dabs(r__2)
			);
/* Computing MAX */
		r__3 = s, r__4 = (r__1 = x.r, dabs(r__1)) + (r__2 = r_imag(&x)
			, dabs(r__2));
		s = dmax(r__3,r__4);
		q__5.r = x.r / s, q__5.i = x.i / s;
		pow_ci(&q__4, &q__5, &c__2);
		q__7.r = u.r / s, q__7.i = u.i / s;
		pow_ci(&q__6, &q__7, &c__2);
		q__3.r = q__4.r + q__6.r, q__3.i = q__4.i + q__6.i;
		c_sqrt(&q__2, &q__3);
		q__1.r = s * q__2.r, q__1.i = s * q__2.i;
		y.r = q__1.r, y.i = q__1.i;
		if (sx > 0.f) {
		    q__1.r = x.r / sx, q__1.i = x.i / sx;
		    q__2.r = x.r / sx, q__2.i = x.i / sx;
		    if (q__1.r * y.r + r_imag(&q__2) * r_imag(&y) < 0.f) {
			q__3.r = -y.r, q__3.i = -y.i;
			y.r = q__3.r, y.i = q__3.i;
		    }
		}
		q__4.r = x.r + y.r, q__4.i = x.i + y.i;
		cladiv_(&q__3, &u, &q__4);
		q__2.r = u.r * q__3.r - u.i * q__3.i, q__2.i = u.r * q__3.i + 
			u.i * q__3.r;
		q__1.r = t.r - q__2.r, q__1.i = t.i - q__2.i;
		t.r = q__1.r, t.i = q__1.i;
	    }
	}

/*        Look for two consecutive small subdiagonal elements. */

	i__2 = l + 1;
	for (m = i__ - 1; m >= i__2; --m) {

/*           Determine the effect of starting the single-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible. */

	    i__3 = m + m * h_dim1;
	    h11.r = h__[i__3].r, h11.i = h__[i__3].i;
	    i__3 = m + 1 + (m + 1) * h_dim1;
	    h22.r = h__[i__3].r, h22.i = h__[i__3].i;
	    q__1.r = h11.r - t.r, q__1.i = h11.i - t.i;
	    h11s.r = q__1.r, h11s.i = q__1.i;
	    i__3 = m + 1 + m * h_dim1;
	    h21 = h__[i__3].r;
	    s = (r__1 = h11s.r, dabs(r__1)) + (r__2 = r_imag(&h11s), dabs(
		    r__2)) + dabs(h21);
	    q__1.r = h11s.r / s, q__1.i = h11s.i / s;
	    h11s.r = q__1.r, h11s.i = q__1.i;
	    h21 /= s;
	    v[0].r = h11s.r, v[0].i = h11s.i;
	    v[1].r = h21, v[1].i = 0.f;
	    i__3 = m + (m - 1) * h_dim1;
	    h10 = h__[i__3].r;
	    if (dabs(h10) * dabs(h21) <= ulp * (((r__1 = h11s.r, dabs(r__1)) 
		    + (r__2 = r_imag(&h11s), dabs(r__2))) * ((r__3 = h11.r, 
		    dabs(r__3)) + (r__4 = r_imag(&h11), dabs(r__4)) + ((r__5 =
		     h22.r, dabs(r__5)) + (r__6 = r_imag(&h22), dabs(r__6)))))
		    ) {
		goto L70;
	    }
/* L60: */
	}
	i__2 = l + l * h_dim1;
	h11.r = h__[i__2].r, h11.i = h__[i__2].i;
	i__2 = l + 1 + (l + 1) * h_dim1;
	h22.r = h__[i__2].r, h22.i = h__[i__2].i;
	q__1.r = h11.r - t.r, q__1.i = h11.i - t.i;
	h11s.r = q__1.r, h11s.i = q__1.i;
	i__2 = l + 1 + l * h_dim1;
	h21 = h__[i__2].r;
	s = (r__1 = h11s.r, dabs(r__1)) + (r__2 = r_imag(&h11s), dabs(r__2)) 
		+ dabs(h21);
	q__1.r = h11s.r / s, q__1.i = h11s.i / s;
	h11s.r = q__1.r, h11s.i = q__1.i;
	h21 /= s;
	v[0].r = h11s.r, v[0].i = h11s.i;
	v[1].r = h21, v[1].i = 0.f;
L70:

/*        Single-shift QR step */

	i__2 = i__ - 1;
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. */

/*           V(2) is always real before the call to CLARFG, and hence */
/*           after the call T2 ( = T1*V(2) ) is also real. */

	    if (k > m) {
		ccopy_(&c__2, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
	    }
	    clarfg_(&c__2, v, &v[1], &c__1, &t1);
	    if (k > m) {
		i__3 = k + (k - 1) * h_dim1;
		h__[i__3].r = v[0].r, h__[i__3].i = v[0].i;
		i__3 = k + 1 + (k - 1) * h_dim1;
		h__[i__3].r = 0.f, h__[i__3].i = 0.f;
	    }
	    v2.r = v[1].r, v2.i = v[1].i;
	    q__1.r = t1.r * v2.r - t1.i * v2.i, q__1.i = t1.r * v2.i + t1.i * 
		    v2.r;
	    t2 = q__1.r;

/*           Apply G from the left to transform the rows of the matrix */
/*           in columns K to I2. */

	    i__3 = i2;
	    for (j = k; j <= i__3; ++j) {
		r_cnjg(&q__3, &t1);
		i__4 = k + j * h_dim1;
		q__2.r = q__3.r * h__[i__4].r - q__3.i * h__[i__4].i, q__2.i =
			 q__3.r * h__[i__4].i + q__3.i * h__[i__4].r;
		i__5 = k + 1 + j * h_dim1;
		q__4.r = t2 * h__[i__5].r, q__4.i = t2 * h__[i__5].i;
		q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		sum.r = q__1.r, sum.i = q__1.i;
		i__4 = k + j * h_dim1;
		i__5 = k + j * h_dim1;
		q__1.r = h__[i__5].r - sum.r, q__1.i = h__[i__5].i - sum.i;
		h__[i__4].r = q__1.r, h__[i__4].i = q__1.i;
		i__4 = k + 1 + j * h_dim1;
		i__5 = k + 1 + j * h_dim1;
		q__2.r = sum.r * v2.r - sum.i * v2.i, q__2.i = sum.r * v2.i + 
			sum.i * v2.r;
		q__1.r = h__[i__5].r - q__2.r, q__1.i = h__[i__5].i - q__2.i;
		h__[i__4].r = q__1.r, h__[i__4].i = q__1.i;
/* L80: */
	    }

/*           Apply G from the right to transform the columns of the */
/*           matrix in rows I1 to min(K+2,I). */

/* Computing MIN */
	    i__4 = k + 2;
	    i__3 = min(i__4,i__);
	    for (j = i1; j <= i__3; ++j) {
		i__4 = j + k * h_dim1;
		q__2.r = t1.r * h__[i__4].r - t1.i * h__[i__4].i, q__2.i = 
			t1.r * h__[i__4].i + t1.i * h__[i__4].r;
		i__5 = j + (k + 1) * h_dim1;
		q__3.r = t2 * h__[i__5].r, q__3.i = t2 * h__[i__5].i;
		q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		sum.r = q__1.r, sum.i = q__1.i;
		i__4 = j + k * h_dim1;
		i__5 = j + k * h_dim1;
		q__1.r = h__[i__5].r - sum.r, q__1.i = h__[i__5].i - sum.i;
		h__[i__4].r = q__1.r, h__[i__4].i = q__1.i;
		i__4 = j + (k + 1) * h_dim1;
		i__5 = j + (k + 1) * h_dim1;
		r_cnjg(&q__3, &v2);
		q__2.r = sum.r * q__3.r - sum.i * q__3.i, q__2.i = sum.r * 
			q__3.i + sum.i * q__3.r;
		q__1.r = h__[i__5].r - q__2.r, q__1.i = h__[i__5].i - q__2.i;
		h__[i__4].r = q__1.r, h__[i__4].i = q__1.i;
/* L90: */
	    }

	    if (*wantz) {

/*              Accumulate transformations in the matrix Z */

		i__3 = *ihiz;
		for (j = *iloz; j <= i__3; ++j) {
		    i__4 = j + k * z_dim1;
		    q__2.r = t1.r * z__[i__4].r - t1.i * z__[i__4].i, q__2.i =
			     t1.r * z__[i__4].i + t1.i * z__[i__4].r;
		    i__5 = j + (k + 1) * z_dim1;
		    q__3.r = t2 * z__[i__5].r, q__3.i = t2 * z__[i__5].i;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    sum.r = q__1.r, sum.i = q__1.i;
		    i__4 = j + k * z_dim1;
		    i__5 = j + k * z_dim1;
		    q__1.r = z__[i__5].r - sum.r, q__1.i = z__[i__5].i - 
			    sum.i;
		    z__[i__4].r = q__1.r, z__[i__4].i = q__1.i;
		    i__4 = j + (k + 1) * z_dim1;
		    i__5 = j + (k + 1) * z_dim1;
		    r_cnjg(&q__3, &v2);
		    q__2.r = sum.r * q__3.r - sum.i * q__3.i, q__2.i = sum.r *
			     q__3.i + sum.i * q__3.r;
		    q__1.r = z__[i__5].r - q__2.r, q__1.i = z__[i__5].i - 
			    q__2.i;
		    z__[i__4].r = q__1.r, z__[i__4].i = q__1.i;
/* L100: */
		}
	    }

	    if (k == m && m > l) {

/*              If the QR step was started at row M > L because two */
/*              consecutive small subdiagonals were found, then extra */
/*              scaling must be performed to ensure that H(M,M-1) remains */
/*              real. */

		q__1.r = 1.f - t1.r, q__1.i = 0.f - t1.i;
		temp.r = q__1.r, temp.i = q__1.i;
		r__1 = c_abs(&temp);
		q__1.r = temp.r / r__1, q__1.i = temp.i / r__1;
		temp.r = q__1.r, temp.i = q__1.i;
		i__3 = m + 1 + m * h_dim1;
		i__4 = m + 1 + m * h_dim1;
		r_cnjg(&q__2, &temp);
		q__1.r = h__[i__4].r * q__2.r - h__[i__4].i * q__2.i, q__1.i =
			 h__[i__4].r * q__2.i + h__[i__4].i * q__2.r;
		h__[i__3].r = q__1.r, h__[i__3].i = q__1.i;
		if (m + 2 <= i__) {
		    i__3 = m + 2 + (m + 1) * h_dim1;
		    i__4 = m + 2 + (m + 1) * h_dim1;
		    q__1.r = h__[i__4].r * temp.r - h__[i__4].i * temp.i, 
			    q__1.i = h__[i__4].r * temp.i + h__[i__4].i * 
			    temp.r;
		    h__[i__3].r = q__1.r, h__[i__3].i = q__1.i;
		}
		i__3 = i__;
		for (j = m; j <= i__3; ++j) {
		    if (j != m + 1) {
			if (i2 > j) {
			    i__4 = i2 - j;
			    cscal_(&i__4, &temp, &h__[j + (j + 1) * h_dim1], 
				    ldh);
			}
			i__4 = j - i1;
			r_cnjg(&q__1, &temp);
			cscal_(&i__4, &q__1, &h__[i1 + j * h_dim1], &c__1);
			if (*wantz) {
			    r_cnjg(&q__1, &temp);
			    cscal_(&nz, &q__1, &z__[*iloz + j * z_dim1], &
				    c__1);
			}
		    }
/* L110: */
		}
	    }
/* L120: */
	}

/*        Ensure that H(I,I-1) is real. */

	i__2 = i__ + (i__ - 1) * h_dim1;
	temp.r = h__[i__2].r, temp.i = h__[i__2].i;
	if (r_imag(&temp) != 0.f) {
	    rtemp = c_abs(&temp);
	    i__2 = i__ + (i__ - 1) * h_dim1;
	    h__[i__2].r = rtemp, h__[i__2].i = 0.f;
	    q__1.r = temp.r / rtemp, q__1.i = temp.i / rtemp;
	    temp.r = q__1.r, temp.i = q__1.i;
	    if (i2 > i__) {
		i__2 = i2 - i__;
		r_cnjg(&q__1, &temp);
		cscal_(&i__2, &q__1, &h__[i__ + (i__ + 1) * h_dim1], ldh);
	    }
	    i__2 = i__ - i1;
	    cscal_(&i__2, &temp, &h__[i1 + i__ * h_dim1], &c__1);
	    if (*wantz) {
		cscal_(&nz, &temp, &z__[*iloz + i__ * z_dim1], &c__1);
	    }
	}

/* L130: */
    }

/*     Failure to converge in remaining number of iterations */

    *info = i__;
    return 0;

L140:

/*     H(I,I-1) is negligible: one eigenvalue has converged. */

    i__1 = i__;
    i__2 = i__ + i__ * h_dim1;
    w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;
/*     reset deflation counter */
    kdefl = 0;

/*     return to start of the main loop with new value of I. */

    i__ = l - 1;
    goto L30;

L150:
    return 0;

/*     End of CLAHQR */

} /* clahqr_ */

/* > \brief \b CSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYCONV convert A given by TRF into L and D and vice-versa. */
/* > Get Non-diag elements of D (returned in workspace) and */
/* > apply or reverse permutation done in TRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N) */
/* >          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1 */
/* >          or 2-by-2 block diagonal matrix D in LDLT. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csyconv_(char *uplo, char *way, integer *n, complex *a, 
	integer *lda, integer *ipiv, complex *e, integer *info, ftnlen 
	uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ip;
    static complex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --e;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CSYCONV", &i__1, (ftnlen)7);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

	if (convert) {
	    i__ = *n;
	    e[1].r = 0.f, e[1].i = 0.f;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ - 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ - 1;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		    i__1 = i__ - 1 + i__ * a_dim1;
		    a[i__1].r = 0.f, a[i__1].i = 0.f;
		    --i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		}
		--i__;
	    }

/*        Convert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L12: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ - 1 + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ - 1 + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L13: */
			}
		    }
		    --i__;
		}
		--i__;
	    }
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    ++i__;
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ - 1 + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ - 1 + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		}
		++i__;
	    }

/*        Revert VALUE */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ - 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    --i__;
		}
		--i__;
	    }
	}
    } else {

/*      A is LOWER */

	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

	    i__ = 1;
	    i__1 = *n;
	    e[i__1].r = 0.f, e[i__1].i = 0.f;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ + 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ + 1;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		    i__1 = i__ + 1 + i__ * a_dim1;
		    a[i__1].r = 0.f, a[i__1].i = 0.f;
		    ++i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		}
		++i__;
	    }

/*        Convert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L22: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + 1 + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + 1 + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L23: */
			}
		    }
		    ++i__;
		}
		++i__;
	    }
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = i__ + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = i__ + j * a_dim1;
			    i__3 = ip + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = ip + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    --i__;
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = i__ + 1 + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = i__ + 1 + j * a_dim1;
			    i__3 = ip + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = ip + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		}
		--i__;
	    }

/*        Revert VALUE */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ + 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    ++i__;
		}
		++i__;
	    }
	}
    }
    return 0;

/*     End of CSYCONV */

} /* csyconv_ */

/* > \brief \b CSYCONVF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYCONVF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv
f.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv
f.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv
f.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > CSYCONVF converts the factorization output format used in */
/* > CSYTRF provided on entry in parameter A into the factorization */
/* > output format used in CSYTRF_RK (or CSYTRF_BK) that is stored */
/* > on exit in parameters A and E. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in CSYTRF into */
/* > the format used in CSYTRF_RK (or CSYTRF_BK). */
/* > */
/* > If parameter WAY = 'R': */
/* > CSYCONVF performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in CSYTRF_RK */
/* > (or CSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in CSYTRF that is stored */
/* > on exit in parameter A. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in CSYTRF_RK */
/* > (or CSYTRF_BK) into the format used in CSYTRF. */
/* > */
/* > CSYCONVF can also convert in Hermitian matrix case, i.e. between */
/* > formats used in CHETRF and CHETRF_RK (or CHETRF_BK). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          CSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          CSYTRF_RK or CSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          CSYTRF_RK or CSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          CSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in CSYTRF. */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in CSYTRF_RK */
/* >          ( or CSYTRF_BK). */
/* > */
/* >          1) If WAY ='R': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in CSYTRF_RK */
/* >          ( or CSYTRF_BK). */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int csyconvf_(char *uplo, char *way, integer *n, complex *a, 
	integer *lda, complex *e, integer *ipiv, integer *info, ftnlen 
	uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ip;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CSYCONVF", &i__1, (ftnlen)8);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1].r = 0.f, e[1].i = 0.f;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ - 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ - 1;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		    i__1 = i__ - 1 + i__ * a_dim1;
		    a[i__1].r = 0.f, a[i__1].i = 0.f;
		    --i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i-1 and IPIV(i-1), */
/*                 so this should be recorded in two consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ - 1];

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ - 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    i__1 = *n;
	    e[i__1].r = 0.f, e[i__1].i = 0.f;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ + 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ + 1;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		    i__1 = i__ + 1 + i__ * a_dim1;
		    a[i__1].r = 0.f, a[i__1].i = 0.f;
		    ++i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where k increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i+1 and IPIV(i+1), */
/*                 so this should be recorded in consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ + 1];

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ + 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of CSYCONVF */

} /* csyconvf_ */

/* > \brief \b CSYCONVF_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYCONVF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv
f_rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv
f_rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv
f_rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > CSYCONVF_ROOK converts the factorization output format used in */
/* > CSYTRF_ROOK provided on entry in parameter A into the factorization */
/* > output format used in CSYTRF_RK (or CSYTRF_BK) that is stored */
/* > on exit in parameters A and E. IPIV format for CSYTRF_ROOK and */
/* > CSYTRF_RK (or CSYTRF_BK) is the same and is not converted. */
/* > */
/* > If parameter WAY = 'R': */
/* > CSYCONVF_ROOK performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in CSYTRF_RK */
/* > (or CSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in CSYTRF_ROOK that is stored */
/* > on exit in parameter A. IPIV format for CSYTRF_ROOK and */
/* > CSYTRF_RK (or CSYTRF_BK) is the same and is not converted. */
/* > */
/* > CSYCONVF_ROOK can also convert in Hermitian matrix case, i.e. between */
/* > formats used in CHETRF_ROOK and CHETRF_RK (or CHETRF_BK). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          CSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          CSYTRF_RK or CSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          CSYTRF_RK or CSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          CSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D as determined: */
/* >          1) by CSYTRF_ROOK, if WAY ='C'; */
/* >          2) by CSYTRF_RK (or CSYTRF_BK), if WAY ='R'. */
/* >          The IPIV format is the same for all these routines. */
/* > */
/* >          On exit, is not changed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int csyconvf_rook__(char *uplo, char *way, integer *n, 
	complex *a, integer *lda, complex *e, integer *ipiv, integer *info, 
	ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ip, ip2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CSYCONVF_ROOK", &i__1, (ftnlen)13);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1].r = 0.f, e[1].i = 0.f;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ - 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ - 1;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		    i__1 = i__ - 1 + i__ * a_dim1;
		    a[i__1].r = 0.f, a[i__1].i = 0.f;
		    --i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1) */
/*                 in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip2 + (i__ + 1) * a_dim1], lda);
			}
		    }
		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i) */
/*                 in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[ip2 + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = *n - i__;
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ - 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    i__1 = *n;
	    e[i__1].r = 0.f, e[i__1].i = 0.f;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ + 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ + 1;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		    i__1 = i__ + 1 + i__ * a_dim1;
		    a[i__1].r = 0.f, a[i__1].i = 0.f;
		    ++i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0.f, e[i__1].i = 0.f;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1) */
/*                 in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip2 + 
				    a_dim1], lda);
			}
		    }
		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i) */
/*                 in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[ip2 + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = i__ - 1;
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ + 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of CSYCONVF_ROOK */

} /* csyconvf_rook__ */

/* > \brief \b CTPMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpmqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpmqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpmqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/*      $          WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPMQRT applies a complex orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" complex block reflector H to a general */
/* > complex matrix C, which consists of two blocks A and B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Conjugate transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The order of the trapezoidal part of V. */
/* >          K >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CTPQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CTPQRT in B.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDV >= max(1,M); */
/* >          if SIDE = 'R', LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CTPQRT, stored as a NB-by-K matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension */
/* >          (LDA,N) if SIDE = 'L' or */
/* >          (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDC >= max(1,K); */
/* >          If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. */
/* >          LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The columns of the pentagonal matrix V contain the elementary reflectors */
/* >  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a */
/* >  trapezoidal block V2: */
/* > */
/* >        V = [V1] */
/* >            [V2]. */
/* > */
/* >  The size of the trapezoidal block V2 is determined by the parameter L, */
/* >  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L */
/* >  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular; */
/* >  if L=0, there is no trapezoidal block, hence V = V1 is rectangular. */
/* > */
/* >  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K. */
/* >                      [B] */
/* > */
/* >  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K. */
/* > */
/* >  The complex orthogonal matrix Q is formed from V and T. */
/* > */
/* >  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C. */
/* > */
/* >  If TRANS='C' and SIDE='L', C is on exit replaced with Q**H * C. */
/* > */
/* >  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* >  If TRANS='C' and SIDE='R', C is on exit replaced with C * Q**H. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctpmqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *l, integer *nb, complex *v, integer *ldv, 
	complex *t, integer *ldt, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *work, integer *info, ftnlen side_len, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, mb, kf, ldaq;
    static logical left, tran;
    static integer ldvq;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), ctprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, complex *, integer *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, complex *, integer *);
    static logical notran;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldvq = max(1,*m);
	ldaq = max(1,*k);
    } else if (right) {
	ldvq = max(1,*n);
	ldaq = max(1,*m);
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0) {
	*info = -5;
    } else if (*l < 0 || *l > *k) {
	*info = -6;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -7;
    } else if (*ldv < ldvq) {
	*info = -9;
    } else if (*ldt < *nb) {
	*info = -11;
    } else if (*lda < ldaq) {
	*info = -13;
    } else if (*ldb < max(1,*m)) {
	*info = -15;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTPMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *m - *l + i__ + ib - 1;
	    mb = min(i__3,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    ctprfb_("L", "C", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *n - *l + i__ + ib - 1;
	    mb = min(i__3,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    ctprfb_("R", "N", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *m - *l + i__ + ib - 1;
	    mb = min(i__2,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    ctprfb_("L", "N", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *n - *l + i__ + ib - 1;
	    mb = min(i__2,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    ctprfb_("R", "C", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    }

    return 0;

/*     End of CTPMQRT */

} /* ctpmqrt_ */

/* > \brief \b CTPQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTPQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPQRT computes a blocked QR factorization of a complex */
/* > "triangular-pentagonal" matrix C, which is composed of a */
/* > triangular block A and pentagonal block B, using the compact */
/* > WY representation for Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B, and the order of the */
/* >          triangular matrix A. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of rows of the upper trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the upper triangular N-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          On entry, the pentagonal M-by-N matrix B.  The first M-L rows */
/* >          are rectangular, and the last L rows are upper trapezoidal. */
/* >          On exit, B contains the pentagonal matrix V.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,N) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a (N+M)-by-N matrix */
/* > */
/* >               C = [ A ] */
/* >                   [ B ] */
/* > */
/* >  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* >  upper trapezoidal matrix B2: */
/* > */
/* >               B = [ B1 ]  <- (M-L)-by-N rectangular */
/* >                   [ B2 ]  <-     L-by-N upper trapezoidal. */
/* > */
/* >  The upper trapezoidal matrix B2 consists of the first L rows of a */
/* >  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is upper triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* >               C = [ A ]  <- upper triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* > */
/* >               W = [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* > */
/* >               V = [ V1 ] <- (M-L)-by-N rectangular */
/* >                   [ V2 ] <-     L-by-N upper trapezoidal. */
/* > */
/* >  The columns of V represent the vectors which define the H(i)'s. */
/* > */
/* >  The number of blocks is B = ceiling(N/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* >               T = [T1 T2 ... TB]. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctpqrt_(integer *m, integer *n, integer *l, integer *nb, 
	complex *a, integer *lda, complex *b, integer *ldb, complex *t, 
	integer *ldt, complex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, ib, lb, mb, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), ctprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, complex *, integer *, complex *, integer *, complex *, 
	    integer *, complex *, integer *, complex *, integer *), ctpqrt2_(integer *, integer *, integer *,
	     complex *, integer *, complex *, integer *, complex *, integer *,
	     integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*l < 0 || *l > min(*m,*n) && min(*m,*n) >= 0) {
	*info = -3;
    } else if (*nb < 1 || *nb > *n && *n > 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*m)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTPQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = *n;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*     Compute the QR factorization of the current block */

/* Computing MIN */
	i__3 = *n - i__ + 1;
	ib = min(i__3,*nb);
/* Computing MIN */
	i__3 = *m - *l + i__ + ib - 1;
	mb = min(i__3,*m);
	if (i__ >= *l) {
	    lb = 0;
	} else {
	    lb = mb - *m + *l - i__ + 1;
	}

	ctpqrt2_(&mb, &ib, &lb, &a[i__ + i__ * a_dim1], lda, &b[i__ * b_dim1 
		+ 1], ldb, &t[i__ * t_dim1 + 1], ldt, &iinfo);

/*     Update by applying H**H to B(:,I+IB:N) from the left */

	if (i__ + ib <= *n) {
	    i__3 = *n - i__ - ib + 1;
	    ctprfb_("L", "C", "F", "C", &mb, &i__3, &ib, &lb, &b[i__ * b_dim1 
		    + 1], ldb, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + ib) 
		    * a_dim1], lda, &b[(i__ + ib) * b_dim1 + 1], ldb, &work[1]
		    , &ib);
	}
    }
    return 0;

/*     End of CTPQRT */

} /* ctpqrt_ */

/* > \brief \b DGEMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgemqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgemqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgemqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEMQRT overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q C            C Q */
/* > TRANS = 'T':   Q**T C            C Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* >       Q = H(1) H(2) . . . H(K) = I - V T V**T */
/* > */
/* > generated using the compact WY representation as returned by DGEQRT. */
/* > */
/* > Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Transpose, apply Q**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CGEQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGEQRT in the first K columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CGEQRT, stored as a NB-by-N matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array. The dimension of */
/* >          WORK is N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgemqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *nb, doublereal *v, integer *ldv, doublereal *t, 
	integer *ldt, doublereal *c__, integer *ldc, doublereal *work, 
	integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, q, ib, kf;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static logical notran;
    static integer ldwork;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldwork = max(1,*n);
	q = *m;
    } else if (right) {
	ldwork = max(1,*m);
	q = *n;
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > q) {
	*info = -5;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -6;
    } else if (*ldv < max(1,q)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    } else if (*ldc < max(1,*m)) {
	*info = -12;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *m - i__ + 1;
	    dlarfb_("L", "T", "F", "C", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *n - i__ + 1;
	    dlarfb_("R", "N", "F", "C", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *m - i__ + 1;
	    dlarfb_("L", "N", "F", "C", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *n - i__ + 1;
	    dlarfb_("R", "T", "F", "C", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    }

    return 0;

/*     End of DGEMQRT */

} /* dgemqrt_ */

/* > \brief \b DGEQRFP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrfp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrfp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrfp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQR2P computes a QR factorization of a real M-by-N matrix A: */
/* > */
/* >    A = Q * ( R ), */
/* >            ( 0 ) */
/* > */
/* > where: */
/* > */
/* >    Q is a M-by-M orthogonal matrix; */
/* >    R is an upper-triangular N-by-N matrix with nonnegative diagonal */
/* >    entries; */
/* >    0 is a (M-N)-by-N zero matrix, if M > N. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n). The diagonal entries of R */
/* >          are nonnegative; the elements below the diagonal, */
/* >          with the array TAU, represent the orthogonal matrix Q as a */
/* >          product of min(m,n) elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,N). */
/* >          For optimum performance LWORK >= N*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > */
/* > See Lapack Working Note 203 for details */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgeqrfp_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarft_(char *, char *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), xerbla_(char *, integer 
	    *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int dgeqr2p_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
    lwkopt = *n * nb;
    work[1] = (doublereal) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQRFP", &i__1, (ftnlen)7);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	work[1] = 1.;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

	    i__3 = *m - i__ + 1;
	    dgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *m - i__ + 1;
		dlarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**T to A(i:m,i+ib:n) from the left */

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib 
			+ 1], &ldwork, (ftnlen)4, (ftnlen)9, (ftnlen)7, (
			ftnlen)10);
	    }
/* L10: */
	}
    } else {
	i__ = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	dgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		1], &iinfo);
    }

    work[1] = (doublereal) iws;
    return 0;

/*     End of DGEQRFP */

} /* dgeqrfp_ */

/* > \brief \b DGEQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDT, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQRT computes a blocked QR factorization of a real M-by-N matrix A */
/* > using the compact WY representation of Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if M >= N); the elements below the diagonal */
/* >          are the columns of V. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,MIN(M,N)) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix V stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* >               V = (  1       ) */
/* >                   ( v1  1    ) */
/* >                   ( v1 v2  1 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  where the vi's represent the vectors which define H(i), which are returned */
/* >  in the matrix A.  The 1's along the diagonal of V are not stored in A. */
/* > */
/* >  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-K matrix T as */
/* > */
/* >               T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgeqrt_(integer *m, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *t, integer *ldt, doublereal *work, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ib, iinfo;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dgeqrt2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), dgeqrt3_(integer *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nb < 1 || *nb > min(*m,*n) && min(*m,*n) > 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*ldt < *nb) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	return 0;
    }

/*     Blocked loop of length K */

    i__1 = k;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	i__3 = k - i__ + 1;
	ib = min(i__3,*nb);

/*     Compute the QR factorization of the current block A(I:M,I:I+IB-1) */

	if (TRUE_) {
	    i__3 = *m - i__ + 1;
	    dgeqrt3_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	} else {
	    i__3 = *m - i__ + 1;
	    dgeqrt2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	}
	if (i__ + ib <= *n) {

/*     Update by applying H**T to A(I:M,I+IB:N) from the left */

	    i__3 = *m - i__ + 1;
	    i__4 = *n - i__ - ib + 1;
	    i__5 = *n - i__ - ib + 1;
	    dlarfb_("L", "T", "F", "C", &i__3, &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + 
		    ib) * a_dim1], lda, &work[1], &i__5, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
	}
    }
    return 0;

/*     End of DGEQRT */

} /* dgeqrt_ */

/* > \brief \b DLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using th
e double-shift/single-shift QR algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAHQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, */
/*                          ILOZ, IHIZ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLAHQR is an auxiliary routine called by DHSEQR to update the */
/* >    eigenvalues and Schur decomposition already computed by DHSEQR, by */
/* >    dealing with the Hessenberg submatrix in rows and columns ILO to */
/* >    IHI. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >          = .TRUE. : the full Schur form T is required; */
/* >          = .FALSE.: only eigenvalues are required. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          = .TRUE. : the matrix of Schur vectors Z is required; */
/* >          = .FALSE.: Schur vectors are not required. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          It is assumed that H is already upper quasi-triangular in */
/* >          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless */
/* >          ILO = 1). DLAHQR works primarily with the Hessenberg */
/* >          submatrix in rows and columns ILO to IHI, but applies */
/* >          transformations to all of H if WANTT is .TRUE.. */
/* >          1 <= ILO <= max(1,IHI); IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION array, dimension (LDH,N) */
/* >          On entry, the upper Hessenberg matrix H. */
/* >          On exit, if INFO is zero and if WANTT is .TRUE., H is upper */
/* >          quasi-triangular in rows and columns ILO:IHI, with any */
/* >          2-by-2 diagonal blocks in standard form. If INFO is zero */
/* >          and WANTT is .FALSE., the contents of H are unspecified on */
/* >          exit.  The output state of H if INFO is nonzero is given */
/* >          below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WR */
/* > \verbatim */
/* >          WR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is DOUBLE PRECISION array, dimension (N) */
/* >          The real and imaginary parts, respectively, of the computed */
/* >          eigenvalues ILO to IHI are stored in the corresponding */
/* >          elements of WR and WI. If two eigenvalues are computed as a */
/* >          complex conjugate pair, they are stored in consecutive */
/* >          elements of WR and WI, say the i-th and (i+1)th, with */
/* >          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the */
/* >          eigenvalues are stored in the same order as on the diagonal */
/* >          of the Schur form returned in H, with WR(i) = H(i,i), and, if */
/* >          H(i:i+1,i:i+1) is a 2-by-2 diagonal block, */
/* >          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* >          ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* >          IHIZ is INTEGER */
/* >          Specify the rows of Z to which transformations must be */
/* >          applied if WANTZ is .TRUE.. */
/* >          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ,N) */
/* >          If WANTZ is .TRUE., on entry Z must contain the current */
/* >          matrix Z of transformations accumulated by DHSEQR, and on */
/* >          exit Z has been updated; transformations are applied only to */
/* >          the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
/* >          If WANTZ is .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:  successful exit */
/* >           > 0:  If INFO = i, DLAHQR failed to compute all the */
/* >                  eigenvalues ILO to IHI in a total of 30 iterations */
/* >                  per eigenvalue; elements i+1:ihi of WR and WI */
/* >                  contain those eigenvalues which have been */
/* >                  successfully computed. */
/* > */
/* >                  If INFO > 0 and WANTT is .FALSE., then on exit, */
/* >                  the remaining unconverged eigenvalues are the */
/* >                  eigenvalues of the upper Hessenberg matrix rows */
/* >                  and columns ILO through INFO of the final, output */
/* >                  value of H. */
/* > */
/* >                  If INFO > 0 and WANTT is .TRUE., then on exit */
/* >          (*)       (initial value of H)*U  = U*(final value of H) */
/* >                  where U is an orthogonal matrix.    The final */
/* >                  value of H is upper Hessenberg and triangular in */
/* >                  rows and columns INFO+1 through IHI. */
/* > */
/* >                  If INFO > 0 and WANTZ is .TRUE., then on exit */
/* >                      (final value of Z)  = (initial value of Z)*U */
/* >                  where U is the orthogonal matrix in (*) */
/* >                  (regardless of the value of WANTT.) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     02-96 Based on modifications by */
/* >     David Day, Sandia National Laboratory, USA */
/* > */
/* >     12-04 Further modifications by */
/* >     Ralph Byers, University of Kansas, USA */
/* >     This is a modified version of DLAHQR from LAPACK version 3.0. */
/* >     It is (1) more robust against overflow and underflow and */
/* >     (2) adopts the more conservative Ahues & Tisseur stopping */
/* >     criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal 
	*wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__, 
	integer *ldz, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal s, v[3];
    static integer i1, i2;
    static doublereal t1, t2, t3, v2, v3, aa, ab, ba, bb, h11, h12, h21, h22, 
	    cs;
    static integer nh;
    static doublereal sn;
    static integer nr;
    static doublereal tr;
    static integer nz;
    static doublereal det, h21s;
    static integer its;
    static doublereal ulp, sum, tst, rt1i, rt2i, rt1r, rt2r;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer kdefl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer itmax;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal safmin, safmax, rtdisc, smlnum;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ========================================================= */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }
    if (*ilo == *ihi) {
	wr[*ilo] = h__[*ilo + *ilo * h_dim1];
	wi[*ilo] = 0.;
	return 0;
    }

/*     ==== clear out the trash ==== */
    i__1 = *ihi - 3;
    for (j = *ilo; j <= i__1; ++j) {
	h__[j + 2 + j * h_dim1] = 0.;
	h__[j + 3 + j * h_dim1] = 0.;
/* L10: */
    }
    if (*ilo <= *ihi - 2) {
	h__[*ihi + (*ihi - 2) * h_dim1] = 0.;
    }

    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion. */

    safmin = dlamch_("SAFE MINIMUM", (ftnlen)12);
    safmax = 1. / safmin;
    dlabad_(&safmin, &safmax);
    ulp = dlamch_("PRECISION", (ftnlen)9);
    smlnum = safmin * ((doublereal) nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

    if (*wantt) {
	i1 = 1;
	i2 = *n;
    }

/*     ITMAX is the total number of QR iterations allowed. */

    itmax = max(10,nh) * 30;

/*     KDEFL counts the number of iterations since a deflation */

    kdefl = 0;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
/*     H(L,L-1) is negligible so that the matrix splits. */

    i__ = *ihi;
L20:
    l = *ilo;
    if (i__ < *ilo) {
	goto L160;
    }

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 or 2 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

    i__1 = itmax;
    for (its = 0; its <= i__1; ++its) {

/*        Look for a single small subdiagonal element. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {
	    if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= smlnum) {
		goto L40;
	    }
	    tst = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs(d__1)) + (d__2 = 
		    h__[k + k * h_dim1], abs(d__2));
	    if (tst == 0.) {
		if (k - 2 >= *ilo) {
		    tst += (d__1 = h__[k - 1 + (k - 2) * h_dim1], abs(d__1));
		}
		if (k + 1 <= *ihi) {
		    tst += (d__1 = h__[k + 1 + k * h_dim1], abs(d__1));
		}
	    }
/*           ==== The following is a conservative small subdiagonal */
/*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122, */
/*           .    1997). It has better mathematical foundation and */
/*           .    improves accuracy in some cases.  ==== */
	    if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= ulp * tst) {
/* Computing MAX */
		d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)), d__4 = (
			d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
		ab = max(d__3,d__4);
/* Computing MIN */
		d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)), d__4 = (
			d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
		ba = min(d__3,d__4);
/* Computing MAX */
		d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)), d__4 = (d__2 =
			 h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], 
			abs(d__2));
		aa = max(d__3,d__4);
/* Computing MIN */
		d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)), d__4 = (d__2 =
			 h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], 
			abs(d__2));
		bb = min(d__3,d__4);
		s = aa + ab;
/* Computing MAX */
		d__1 = smlnum, d__2 = ulp * (bb * (aa / s));
		if (ba * (ab / s) <= max(d__1,d__2)) {
		    goto L40;
		}
	    }
/* L30: */
	}
L40:
	l = k;
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

	    h__[l + (l - 1) * h_dim1] = 0.;
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

	if (l >= i__ - 1) {
	    goto L150;
	}
	++kdefl;

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

	if (! (*wantt)) {
	    i1 = l;
	    i2 = i__;
	}

	if (kdefl % 20 == 0) {

/*           Exceptional shift. */

	    s = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs(d__1)) + (d__2 = 
		    h__[i__ - 1 + (i__ - 2) * h_dim1], abs(d__2));
	    h11 = s * .75 + h__[i__ + i__ * h_dim1];
	    h12 = s * -.4375;
	    h21 = s;
	    h22 = h11;
	} else if (kdefl % 10 == 0) {

/*           Exceptional shift. */

	    s = (d__1 = h__[l + 1 + l * h_dim1], abs(d__1)) + (d__2 = h__[l + 
		    2 + (l + 1) * h_dim1], abs(d__2));
	    h11 = s * .75 + h__[l + l * h_dim1];
	    h12 = s * -.4375;
	    h21 = s;
	    h22 = h11;
	} else {

/*           Prepare to use Francis' double shift */
/*           (i.e. 2nd degree generalized Rayleigh quotient) */

	    h11 = h__[i__ - 1 + (i__ - 1) * h_dim1];
	    h21 = h__[i__ + (i__ - 1) * h_dim1];
	    h12 = h__[i__ - 1 + i__ * h_dim1];
	    h22 = h__[i__ + i__ * h_dim1];
	}
	s = abs(h11) + abs(h12) + abs(h21) + abs(h22);
	if (s == 0.) {
	    rt1r = 0.;
	    rt1i = 0.;
	    rt2r = 0.;
	    rt2i = 0.;
	} else {
	    h11 /= s;
	    h21 /= s;
	    h12 /= s;
	    h22 /= s;
	    tr = (h11 + h22) / 2.;
	    det = (h11 - tr) * (h22 - tr) - h12 * h21;
	    rtdisc = sqrt((abs(det)));
	    if (det >= 0.) {

/*              ==== complex conjugate shifts ==== */

		rt1r = tr * s;
		rt2r = rt1r;
		rt1i = rtdisc * s;
		rt2i = -rt1i;
	    } else {

/*              ==== real shifts (use only one of them)  ==== */

		rt1r = tr + rtdisc;
		rt2r = tr - rtdisc;
		if ((d__1 = rt1r - h22, abs(d__1)) <= (d__2 = rt2r - h22, abs(
			d__2))) {
		    rt1r *= s;
		    rt2r = rt1r;
		} else {
		    rt2r *= s;
		    rt1r = rt2r;
		}
		rt1i = 0.;
		rt2i = 0.;
	    }
	}

/*        Look for two consecutive small subdiagonal elements. */

	i__2 = l;
	for (m = i__ - 2; m >= i__2; --m) {
/*           Determine the effect of starting the double-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible.  (The following uses scaling to avoid */
/*           overflows and most underflows.) */

	    h21s = h__[m + 1 + m * h_dim1];
	    s = (d__1 = h__[m + m * h_dim1] - rt2r, abs(d__1)) + abs(rt2i) + 
		    abs(h21s);
	    h21s = h__[m + 1 + m * h_dim1] / s;
	    v[0] = h21s * h__[m + (m + 1) * h_dim1] + (h__[m + m * h_dim1] - 
		    rt1r) * ((h__[m + m * h_dim1] - rt2r) / s) - rt1i * (rt2i 
		    / s);
	    v[1] = h21s * (h__[m + m * h_dim1] + h__[m + 1 + (m + 1) * h_dim1]
		     - rt1r - rt2r);
	    v[2] = h21s * h__[m + 2 + (m + 1) * h_dim1];
	    s = abs(v[0]) + abs(v[1]) + abs(v[2]);
	    v[0] /= s;
	    v[1] /= s;
	    v[2] /= s;
	    if (m == l) {
		goto L60;
	    }
	    if ((d__1 = h__[m + (m - 1) * h_dim1], abs(d__1)) * (abs(v[1]) + 
		    abs(v[2])) <= ulp * abs(v[0]) * ((d__2 = h__[m - 1 + (m - 
		    1) * h_dim1], abs(d__2)) + (d__3 = h__[m + m * h_dim1], 
		    abs(d__3)) + (d__4 = h__[m + 1 + (m + 1) * h_dim1], abs(
		    d__4)))) {
		goto L60;
	    }
/* L50: */
	}
L60:

/*        Double-shift QR step */

	i__2 = i__ - 1;
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. NR is the order of G. */

/* Computing MIN */
	    i__3 = 3, i__4 = i__ - k + 1;
	    nr = min(i__3,i__4);
	    if (k > m) {
		dcopy_(&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
	    }
	    dlarfg_(&nr, v, &v[1], &c__1, &t1);
	    if (k > m) {
		h__[k + (k - 1) * h_dim1] = v[0];
		h__[k + 1 + (k - 1) * h_dim1] = 0.;
		if (k < i__ - 1) {
		    h__[k + 2 + (k - 1) * h_dim1] = 0.;
		}
	    } else if (m > l) {
/*               ==== Use the following instead of */
/*               .    H( K, K-1 ) = -H( K, K-1 ) to */
/*               .    avoid a bug when v(2) and v(3) */
/*               .    underflow. ==== */
		h__[k + (k - 1) * h_dim1] *= 1. - t1;
	    }
	    v2 = v[1];
	    t2 = t1 * v2;
	    if (nr == 3) {
		v3 = v[2];
		t3 = t1 * v3;

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

		i__3 = i2;
		for (j = k; j <= i__3; ++j) {
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1] 
			    + v3 * h__[k + 2 + j * h_dim1];
		    h__[k + j * h_dim1] -= sum * t1;
		    h__[k + 1 + j * h_dim1] -= sum * t2;
		    h__[k + 2 + j * h_dim1] -= sum * t3;
/* L70: */
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

/* Computing MIN */
		i__4 = k + 3;
		i__3 = min(i__4,i__);
		for (j = i1; j <= i__3; ++j) {
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			     + v3 * h__[j + (k + 2) * h_dim1];
		    h__[j + k * h_dim1] -= sum * t1;
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
		    h__[j + (k + 2) * h_dim1] -= sum * t3;
/* L80: */
		}

		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

		    i__3 = *ihiz;
		    for (j = *iloz; j <= i__3; ++j) {
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1] + v3 * z__[j + (k + 2) * z_dim1];
			z__[j + k * z_dim1] -= sum * t1;
			z__[j + (k + 1) * z_dim1] -= sum * t2;
			z__[j + (k + 2) * z_dim1] -= sum * t3;
/* L90: */
		    }
		}
	    } else if (nr == 2) {

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

		i__3 = i2;
		for (j = k; j <= i__3; ++j) {
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1];
		    h__[k + j * h_dim1] -= sum * t1;
		    h__[k + 1 + j * h_dim1] -= sum * t2;
/* L100: */
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

		i__3 = i__;
		for (j = i1; j <= i__3; ++j) {
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			    ;
		    h__[j + k * h_dim1] -= sum * t1;
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
/* L110: */
		}

		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

		    i__3 = *ihiz;
		    for (j = *iloz; j <= i__3; ++j) {
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1];
			z__[j + k * z_dim1] -= sum * t1;
			z__[j + (k + 1) * z_dim1] -= sum * t2;
/* L120: */
		    }
		}
	    }
/* L130: */
	}

/* L140: */
    }

/*     Failure to converge in remaining number of iterations */

    *info = i__;
    return 0;

L150:

    if (l == i__) {

/*        H(I,I-1) is negligible: one eigenvalue has converged. */

	wr[i__] = h__[i__ + i__ * h_dim1];
	wi[i__] = 0.;
    } else if (l == i__ - 1) {

/*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged. */

/*        Transform the 2-by-2 submatrix to standard Schur form, */
/*        and compute and store the eigenvalues. */

	dlanv2_(&h__[i__ - 1 + (i__ - 1) * h_dim1], &h__[i__ - 1 + i__ * 
		h_dim1], &h__[i__ + (i__ - 1) * h_dim1], &h__[i__ + i__ * 
		h_dim1], &wr[i__ - 1], &wi[i__ - 1], &wr[i__], &wi[i__], &cs, 
		&sn);

	if (*wantt) {

/*           Apply the transformation to the rest of H. */

	    if (i2 > i__) {
		i__1 = i2 - i__;
		drot_(&i__1, &h__[i__ - 1 + (i__ + 1) * h_dim1], ldh, &h__[
			i__ + (i__ + 1) * h_dim1], ldh, &cs, &sn);
	    }
	    i__1 = i__ - i1 - 1;
	    drot_(&i__1, &h__[i1 + (i__ - 1) * h_dim1], &c__1, &h__[i1 + i__ *
		     h_dim1], &c__1, &cs, &sn);
	}
	if (*wantz) {

/*           Apply the transformation to Z. */

	    drot_(&nz, &z__[*iloz + (i__ - 1) * z_dim1], &c__1, &z__[*iloz + 
		    i__ * z_dim1], &c__1, &cs, &sn);
	}
    }
/*     reset deflation counter */
    kdefl = 0;

/*     return to start of the main loop with new value of I. */

    i__ = l - 1;
    goto L20;

L160:
    return 0;

/*     End of DLAHQR */

} /* dlahqr_ */

/* > \brief \b DSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYCONV convert A given by TRF into L and D and vice-versa. */
/* > Get Non-diag elements of D (returned in workspace) and */
/* > apply or reverse permutation done in TRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1 */
/* >          or 2-by-2 block diagonal matrix D in LDLT. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsyconv_(char *uplo, char *way, integer *n, doublereal *
	a, integer *lda, integer *ipiv, doublereal *e, integer *info, ftnlen 
	uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, j, ip;
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --e;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYCONV", &i__1, (ftnlen)7);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

	if (convert) {
	    i__ = *n;
	    e[1] = 0.;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
		    e[i__ - 1] = 0.;
		    a[i__ - 1 + i__ * a_dim1] = 0.;
		    --i__;
		} else {
		    e[i__] = 0.;
		}
		--i__;
	    }

/*        Convert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = temp;
/* L12: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
			    a[i__ - 1 + j * a_dim1] = temp;
/* L13: */
			}
		    }
		    --i__;
		}
		--i__;
	    }
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = temp;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    ++i__;
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
			    a[i__ - 1 + j * a_dim1] = temp;
			}
		    }
		}
		++i__;
	    }

/*        Revert VALUE */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
		    --i__;
		}
		--i__;
	    }
	}
    } else {

/*      A is LOWER */

	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

	    i__ = 1;
	    e[*n] = 0.;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
		    e[i__ + 1] = 0.;
		    a[i__ + 1 + i__ * a_dim1] = 0.;
		    ++i__;
		} else {
		    e[i__] = 0.;
		}
		++i__;
	    }

/*        Convert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = temp;
/* L22: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
			    a[i__ + 1 + j * a_dim1] = temp;
/* L23: */
			}
		    }
		    ++i__;
		}
		++i__;
	    }
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = temp;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    --i__;
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[i__ + 1 + j * a_dim1];
			    a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = temp;
			}
		    }
		}
		--i__;
	    }

/*        Revert VALUE */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
		    ++i__;
		}
		++i__;
	    }
	}
    }
    return 0;

/*     End of DSYCONV */

} /* dsyconv_ */

/* > \brief \b DSYCONVF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYCONVF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyconv
f.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyconv
f.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyconv
f.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > DSYCONVF converts the factorization output format used in */
/* > DSYTRF provided on entry in parameter A into the factorization */
/* > output format used in DSYTRF_RK (or DSYTRF_BK) that is stored */
/* > on exit in parameters A and E. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in DSYTRF into */
/* > the format used in DSYTRF_RK (or DSYTRF_BK). */
/* > */
/* > If parameter WAY = 'R': */
/* > DSYCONVF performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in DSYTRF_RK */
/* > (or DSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in DSYTRF that is stored */
/* > on exit in parameter A. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in DSYTRF_RK */
/* > (or DSYTRF_BK) into the format used in DSYTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          DSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          DSYTRF_RK or DSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          DSYTRF_RK or DSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          DSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in DSYTRF. */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in DSYTRF_RK */
/* >          ( or DSYTRF_BK). */
/* > */
/* >          1) If WAY ='R': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in DSYTRF_RK */
/* >          ( or DSYTRF_BK). */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int dsyconvf_(char *uplo, char *way, integer *n, doublereal *
	a, integer *lda, doublereal *e, integer *ipiv, integer *info, ftnlen 
	uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, ip;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYCONVF", &i__1, (ftnlen)8);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1] = 0.;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
		    e[i__ - 1] = 0.;
		    a[i__ - 1 + i__ * a_dim1] = 0.;
		    --i__;
		} else {
		    e[i__] = 0.;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i-1 and IPIV(i-1), */
/*                 so this should be recorded in two consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ - 1];

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    e[*n] = 0.;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
		    e[i__ + 1] = 0.;
		    a[i__ + 1 + i__ * a_dim1] = 0.;
		    ++i__;
		} else {
		    e[i__] = 0.;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where k increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i+1 and IPIV(i+1), */
/*                 so this should be recorded in consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ + 1];

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of DSYCONVF */

} /* dsyconvf_ */

/* > \brief \b DSYCONVF_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYCONVF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyconv
f_rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyconv
f_rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyconv
f_rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > DSYCONVF_ROOK converts the factorization output format used in */
/* > DSYTRF_ROOK provided on entry in parameter A into the factorization */
/* > output format used in DSYTRF_RK (or DSYTRF_BK) that is stored */
/* > on exit in parameters A and E. IPIV format for DSYTRF_ROOK and */
/* > DSYTRF_RK (or DSYTRF_BK) is the same and is not converted. */
/* > */
/* > If parameter WAY = 'R': */
/* > DSYCONVF_ROOK performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in DSYTRF_RK */
/* > (or DSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in DSYTRF_ROOK that is stored */
/* > on exit in parameter A. IPIV format for DSYTRF_ROOK and */
/* > DSYTRF_RK (or DSYTRF_BK) is the same and is not converted. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          DSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          DSYTRF_RK or DSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          DSYTRF_RK or DSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          DSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D as determined: */
/* >          1) by DSYTRF_ROOK, if WAY ='C'; */
/* >          2) by DSYTRF_RK (or DSYTRF_BK), if WAY ='R'. */
/* >          The IPIV format is the same for all these routines. */
/* > */
/* >          On exit, is not changed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int dsyconvf_rook__(char *uplo, char *way, integer *n, 
	doublereal *a, integer *lda, doublereal *e, integer *ipiv, integer *
	info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, ip, ip2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYCONVF_ROOK", &i__1, (ftnlen)13);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1] = 0.;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
		    e[i__ - 1] = 0.;
		    a[i__ - 1 + i__ * a_dim1] = 0.;
		    --i__;
		} else {
		    e[i__] = 0.;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1) */
/*                 in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip2 + (i__ + 1) * a_dim1], lda);
			}
		    }
		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i) */
/*                 in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[ip2 + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = *n - i__;
			    dswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    e[*n] = 0.;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
		    e[i__ + 1] = 0.;
		    a[i__ + 1 + i__ * a_dim1] = 0.;
		    ++i__;
		} else {
		    e[i__] = 0.;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1) */
/*                 in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip2 + 
				    a_dim1], lda);
			}
		    }
		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i) */
/*                 in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[ip2 + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = i__ - 1;
			    dswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of DSYCONVF_ROOK */

} /* dsyconvf_rook__ */

/* > \brief \b DTPMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpmqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpmqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpmqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   V( LDV, * ), A( LDA, * ), B( LDB, * ), */
/*      $                   T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPMQRT applies a real orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" real block reflector H to a general */
/* > real matrix C, which consists of two blocks A and B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The order of the trapezoidal part of V. */
/* >          K >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CTPQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CTPQRT in B.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDV >= max(1,M); */
/* >          if SIDE = 'R', LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CTPQRT, stored as a NB-by-K matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension */
/* >          (LDA,N) if SIDE = 'L' or */
/* >          (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDC >= max(1,K); */
/* >          If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. */
/* >          LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The columns of the pentagonal matrix V contain the elementary reflectors */
/* >  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a */
/* >  trapezoidal block V2: */
/* > */
/* >        V = [V1] */
/* >            [V2]. */
/* > */
/* >  The size of the trapezoidal block V2 is determined by the parameter L, */
/* >  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L */
/* >  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular; */
/* >  if L=0, there is no trapezoidal block, hence V = V1 is rectangular. */
/* > */
/* >  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K. */
/* >                      [B] */
/* > */
/* >  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K. */
/* > */
/* >  The real orthogonal matrix Q is formed from V and T. */
/* > */
/* >  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C. */
/* > */
/* >  If TRANS='T' and SIDE='L', C is on exit replaced with Q**T * C. */
/* > */
/* >  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* >  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q**T. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtpmqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *l, integer *nb, doublereal *v, integer *ldv, 
	doublereal *t, integer *ldt, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *work, integer *info, ftnlen side_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, mb, kf, ldaq;
    static logical left, tran;
    static integer ldvq;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), 
	
	dtprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical notran;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldvq = max(1,*m);
	ldaq = max(1,*k);
    } else if (right) {
	ldvq = max(1,*n);
	ldaq = max(1,*m);
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0) {
	*info = -5;
    } else if (*l < 0 || *l > *k) {
	*info = -6;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -7;
    } else if (*ldv < ldvq) {
	*info = -9;
    } else if (*ldt < *nb) {
	*info = -11;
    } else if (*lda < ldaq) {
	*info = -13;
    } else if (*ldb < max(1,*m)) {
	*info = -15;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTPMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *m - *l + i__ + ib - 1;
	    mb = min(i__3,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    dtprfb_("L", "T", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *n - *l + i__ + ib - 1;
	    mb = min(i__3,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    dtprfb_("R", "N", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *m - *l + i__ + ib - 1;
	    mb = min(i__2,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    dtprfb_("L", "N", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *n - *l + i__ + ib - 1;
	    mb = min(i__2,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    dtprfb_("R", "T", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    }

    return 0;

/*     End of DTPMQRT */

} /* dtpmqrt_ */

/* > \brief \b DTPQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPQRT computes a blocked QR factorization of a real */
/* > "triangular-pentagonal" matrix C, which is composed of a */
/* > triangular block A and pentagonal block B, using the compact */
/* > WY representation for Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B, and the order of the */
/* >          triangular matrix A. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of rows of the upper trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the upper triangular N-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          On entry, the pentagonal M-by-N matrix B.  The first M-L rows */
/* >          are rectangular, and the last L rows are upper trapezoidal. */
/* >          On exit, B contains the pentagonal matrix V.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a (N+M)-by-N matrix */
/* > */
/* >               C = [ A ] */
/* >                   [ B ] */
/* > */
/* >  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* >  upper trapezoidal matrix B2: */
/* > */
/* >               B = [ B1 ]  <- (M-L)-by-N rectangular */
/* >                   [ B2 ]  <-     L-by-N upper trapezoidal. */
/* > */
/* >  The upper trapezoidal matrix B2 consists of the first L rows of a */
/* >  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is upper triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* >               C = [ A ]  <- upper triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* > */
/* >               W = [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* > */
/* >               V = [ V1 ] <- (M-L)-by-N rectangular */
/* >                   [ V2 ] <-     L-by-N upper trapezoidal. */
/* > */
/* >  The columns of V represent the vectors which define the H(i)'s. */
/* > */
/* >  The number of blocks is B = ceiling(N/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* >               T = [T1 T2 ... TB]. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtpqrt_(integer *m, integer *n, integer *l, integer *nb, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	t, integer *ldt, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, ib, lb, mb, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), 
	
dtprfb_(char *, char *, char *, char *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, integer *), 
		
		dtpqrt2_(integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*l < 0 || *l > min(*m,*n) && min(*m,*n) >= 0) {
	*info = -3;
    } else if (*nb < 1 || *nb > *n && *n > 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*m)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTPQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = *n;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*     Compute the QR factorization of the current block */

/* Computing MIN */
	i__3 = *n - i__ + 1;
	ib = min(i__3,*nb);
/* Computing MIN */
	i__3 = *m - *l + i__ + ib - 1;
	mb = min(i__3,*m);
	if (i__ >= *l) {
	    lb = 0;
	} else {
	    lb = mb - *m + *l - i__ + 1;
	}

	dtpqrt2_(&mb, &ib, &lb, &a[i__ + i__ * a_dim1], lda, &b[i__ * b_dim1 
		+ 1], ldb, &t[i__ * t_dim1 + 1], ldt, &iinfo);

/*     Update by applying H**T to B(:,I+IB:N) from the left */

	if (i__ + ib <= *n) {
	    i__3 = *n - i__ - ib + 1;
	    dtprfb_("L", "T", "F", "C", &mb, &i__3, &ib, &lb, &b[i__ * b_dim1 
		    + 1], ldb, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + ib) 
		    * a_dim1], lda, &b[(i__ + ib) * b_dim1 + 1], ldb, &work[1]
		    , &ib);
	}}
    return 0;

/*     End of DTPQRT */

} /* dtpqrt_ */

/* > \brief \b SGEMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgemqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgemqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgemqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEMQRT overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q C            C Q */
/* > TRANS = 'T':   Q**T C            C Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* >       Q = H(1) H(2) . . . H(K) = I - V T V**T */
/* > */
/* > generated using the compact WY representation as returned by SGEQRT. */
/* > */
/* > Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CGEQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is REAL array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGEQRT in the first K columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CGEQRT, stored as a NB-by-N matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgemqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *nb, real *v, integer *ldv, real *t, integer *
	ldt, real *c__, integer *ldc, real *work, integer *info, ftnlen 
	side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, q, ib, kf;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int slarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static integer ldwork;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldwork = max(1,*n);
	q = *m;
    } else if (right) {
	ldwork = max(1,*m);
	q = *n;
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > q) {
	*info = -5;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -6;
    } else if (*ldv < max(1,q)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    } else if (*ldc < max(1,*m)) {
	*info = -12;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *m - i__ + 1;
	    slarfb_("L", "T", "F", "C", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *n - i__ + 1;
	    slarfb_("R", "N", "F", "C", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *m - i__ + 1;
	    slarfb_("L", "N", "F", "C", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *n - i__ + 1;
	    slarfb_("R", "T", "F", "C", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    }

    return 0;

/*     End of SGEMQRT */

} /* sgemqrt_ */

/* > \brief \b SGEQRFP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrfp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrfp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrfp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEQR2P computes a QR factorization of a real M-by-N matrix A: */
/* > */
/* >    A = Q * ( R ), */
/* >            ( 0 ) */
/* > */
/* > where: */
/* > */
/* >    Q is a M-by-M orthogonal matrix; */
/* >    R is an upper-triangular N-by-N matrix with nonnegative diagonal */
/* >    entries; */
/* >    0 is a (M-N)-by-N zero matrix, if M > N. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n). The diagonal entries of R */
/* >          are nonnegative; the elements below the diagonal, */
/* >          with the array TAU, represent the orthogonal matrix Q as a */
/* >          product of min(m,n) elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,N). */
/* >          For optimum performance LWORK >= N*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > */
/* > See Lapack Working Note 203 for details */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgeqrfp_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int slarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarft_(char *, char *, integer *, integer *, 
	    real *, integer *, real *, real *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sgeqr2p_(integer *, integer *, real *, 
	    integer *, real *, real *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
    lwkopt = *n * nb;
    work[1] = (real) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEQRFP", &i__1, (ftnlen)7);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	work[1] = 1.f;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "SGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "SGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

	    i__3 = *m - i__ + 1;
	    sgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *m - i__ + 1;
		slarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**T to A(i:m,i+ib:n) from the left */

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		slarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib 
			+ 1], &ldwork, (ftnlen)4, (ftnlen)9, (ftnlen)7, (
			ftnlen)10);
	    }
/* L10: */
	}
    } else {
	i__ = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	sgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		1], &iinfo);
    }

    work[1] = (real) iws;
    return 0;

/*     End of SGEQRFP */

} /* sgeqrfp_ */

/* > \brief \b SGEQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDT, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEQRT computes a blocked QR factorization of a real M-by-N matrix A */
/* > using the compact WY representation of Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if M >= N); the elements below the diagonal */
/* >          are the columns of V. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,MIN(M,N)) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix V stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* >               V = (  1       ) */
/* >                   ( v1  1    ) */
/* >                   ( v1 v2  1 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  where the vi's represent the vectors which define H(i), which are returned */
/* >  in the matrix A.  The 1's along the diagonal of V are not stored in A. */
/* > */
/* >  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-K matrix T as */
/* > */
/* >               T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgeqrt_(integer *m, integer *n, integer *nb, real *a, 
	integer *lda, real *t, integer *ldt, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ib, iinfo;
    extern /* Subroutine */ int slarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), sgeqrt2_(
	    integer *, integer *, real *, integer *, real *, integer *, 
	    integer *), sgeqrt3_(integer *, integer *, real *, integer *, 
	    real *, integer *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nb < 1 || *nb > min(*m,*n) && min(*m,*n) > 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*ldt < *nb) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	return 0;
    }

/*     Blocked loop of length K */

    i__1 = k;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	i__3 = k - i__ + 1;
	ib = min(i__3,*nb);

/*     Compute the QR factorization of the current block A(I:M,I:I+IB-1) */

	if (TRUE_) {
	    i__3 = *m - i__ + 1;
	    sgeqrt3_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	} else {
	    i__3 = *m - i__ + 1;
	    sgeqrt2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	}
	if (i__ + ib <= *n) {

/*     Update by applying H**T to A(I:M,I+IB:N) from the left */

	    i__3 = *m - i__ + 1;
	    i__4 = *n - i__ - ib + 1;
	    i__5 = *n - i__ - ib + 1;
	    slarfb_("L", "T", "F", "C", &i__3, &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + 
		    ib) * a_dim1], lda, &work[1], &i__5, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
	}
    }
    return 0;

/*     End of SGEQRT */

} /* sgeqrt_ */

/* > \brief \b SLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using th
e double-shift/single-shift QR algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAHQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, */
/*                          ILOZ, IHIZ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLAHQR is an auxiliary routine called by SHSEQR to update the */
/* >    eigenvalues and Schur decomposition already computed by SHSEQR, by */
/* >    dealing with the Hessenberg submatrix in rows and columns ILO to */
/* >    IHI. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >          = .TRUE. : the full Schur form T is required; */
/* >          = .FALSE.: only eigenvalues are required. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          = .TRUE. : the matrix of Schur vectors Z is required; */
/* >          = .FALSE.: Schur vectors are not required. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          It is assumed that H is already upper quasi-triangular in */
/* >          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless */
/* >          ILO = 1). SLAHQR works primarily with the Hessenberg */
/* >          submatrix in rows and columns ILO to IHI, but applies */
/* >          transformations to all of H if WANTT is .TRUE.. */
/* >          1 <= ILO <= max(1,IHI); IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is REAL array, dimension (LDH,N) */
/* >          On entry, the upper Hessenberg matrix H. */
/* >          On exit, if INFO is zero and if WANTT is .TRUE., H is upper */
/* >          quasi-triangular in rows and columns ILO:IHI, with any */
/* >          2-by-2 diagonal blocks in standard form. If INFO is zero */
/* >          and WANTT is .FALSE., the contents of H are unspecified on */
/* >          exit.  The output state of H if INFO is nonzero is given */
/* >          below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WR */
/* > \verbatim */
/* >          WR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is REAL array, dimension (N) */
/* >          The real and imaginary parts, respectively, of the computed */
/* >          eigenvalues ILO to IHI are stored in the corresponding */
/* >          elements of WR and WI. If two eigenvalues are computed as a */
/* >          complex conjugate pair, they are stored in consecutive */
/* >          elements of WR and WI, say the i-th and (i+1)th, with */
/* >          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the */
/* >          eigenvalues are stored in the same order as on the diagonal */
/* >          of the Schur form returned in H, with WR(i) = H(i,i), and, if */
/* >          H(i:i+1,i:i+1) is a 2-by-2 diagonal block, */
/* >          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* >          ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* >          IHIZ is INTEGER */
/* >          Specify the rows of Z to which transformations must be */
/* >          applied if WANTZ is .TRUE.. */
/* >          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ,N) */
/* >          If WANTZ is .TRUE., on entry Z must contain the current */
/* >          matrix Z of transformations accumulated by SHSEQR, and on */
/* >          exit Z has been updated; transformations are applied only to */
/* >          the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
/* >          If WANTZ is .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:   successful exit */
/* >           > 0:   If INFO = i, SLAHQR failed to compute all the */
/* >                  eigenvalues ILO to IHI in a total of 30 iterations */
/* >                  per eigenvalue; elements i+1:ihi of WR and WI */
/* >                  contain those eigenvalues which have been */
/* >                  successfully computed. */
/* > */
/* >                  If INFO > 0 and WANTT is .FALSE., then on exit, */
/* >                  the remaining unconverged eigenvalues are the */
/* >                  eigenvalues of the upper Hessenberg matrix rows */
/* >                  and columns ILO through INFO of the final, output */
/* >                  value of H. */
/* > */
/* >                  If INFO > 0 and WANTT is .TRUE., then on exit */
/* >          (*)       (initial value of H)*U  = U*(final value of H) */
/* >                  where U is an orthogonal matrix.    The final */
/* >                  value of H is upper Hessenberg and triangular in */
/* >                  rows and columns INFO+1 through IHI. */
/* > */
/* >                  If INFO > 0 and WANTZ is .TRUE., then on exit */
/* >                      (final value of Z)  = (initial value of Z)*U */
/* >                  where U is the orthogonal matrix in (*) */
/* >                  (regardless of the value of WANTT.) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup realOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     02-96 Based on modifications by */
/* >     David Day, Sandia National Laboratory, USA */
/* > */
/* >     12-04 Further modifications by */
/* >     Ralph Byers, University of Kansas, USA */
/* >     This is a modified version of SLAHQR from LAPACK version 3.0. */
/* >     It is (1) more robust against overflow and underflow and */
/* >     (2) adopts the more conservative Ahues & Tisseur stopping */
/* >     criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, real *h__, integer *ldh, real *wr, real *
	wi, integer *iloz, integer *ihiz, real *z__, integer *ldz, integer *
	info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m;
    static real s, v[3];
    static integer i1, i2;
    static real t1, t2, t3, v2, v3, aa, ab, ba, bb, h11, h12, h21, h22, cs;
    static integer nh;
    static real sn;
    static integer nr;
    static real tr;
    static integer nz;
    static real det, h21s;
    static integer its;
    static real ulp, sum, tst, rt1i, rt2i, rt1r, rt2r;
    extern /* Subroutine */ int srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *);
    static integer kdefl, itmax;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), slanv2_(real *, real *, real *, real *, real *, real *
	    , real *, real *, real *, real *), slabad_(real *, real *);
    extern doublereal slamch_(char *, ftnlen);
    static real safmin;
    extern /* Subroutine */ int slarfg_(integer *, real *, real *, integer *, 
	    real *);
    static real safmax, rtdisc, smlnum;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ========================================================= */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }
    if (*ilo == *ihi) {
	wr[*ilo] = h__[*ilo + *ilo * h_dim1];
	wi[*ilo] = 0.f;
	return 0;
    }

/*     ==== clear out the trash ==== */
    i__1 = *ihi - 3;
    for (j = *ilo; j <= i__1; ++j) {
	h__[j + 2 + j * h_dim1] = 0.f;
	h__[j + 3 + j * h_dim1] = 0.f;
/* L10: */
    }
    if (*ilo <= *ihi - 2) {
	h__[*ihi + (*ihi - 2) * h_dim1] = 0.f;
    }

    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion. */

    safmin = slamch_("SAFE MINIMUM", (ftnlen)12);
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    ulp = slamch_("PRECISION", (ftnlen)9);
    smlnum = safmin * ((real) nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

    if (*wantt) {
	i1 = 1;
	i2 = *n;
    }

/*     ITMAX is the total number of QR iterations allowed. */

    itmax = max(10,nh) * 30;

/*     KDEFL counts the number of iterations since a deflation */

    kdefl = 0;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
/*     H(L,L-1) is negligible so that the matrix splits. */

    i__ = *ihi;
L20:
    l = *ilo;
    if (i__ < *ilo) {
	goto L160;
    }

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 or 2 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

    i__1 = itmax;
    for (its = 0; its <= i__1; ++its) {

/*        Look for a single small subdiagonal element. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {
	    if ((r__1 = h__[k + (k - 1) * h_dim1], dabs(r__1)) <= smlnum) {
		goto L40;
	    }
	    tst = (r__1 = h__[k - 1 + (k - 1) * h_dim1], dabs(r__1)) + (r__2 =
		     h__[k + k * h_dim1], dabs(r__2));
	    if (tst == 0.f) {
		if (k - 2 >= *ilo) {
		    tst += (r__1 = h__[k - 1 + (k - 2) * h_dim1], dabs(r__1));
		}
		if (k + 1 <= *ihi) {
		    tst += (r__1 = h__[k + 1 + k * h_dim1], dabs(r__1));
		}
	    }
/*           ==== The following is a conservative small subdiagonal */
/*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122, */
/*           .    1997). It has better mathematical foundation and */
/*           .    improves accuracy in some cases.  ==== */
	    if ((r__1 = h__[k + (k - 1) * h_dim1], dabs(r__1)) <= ulp * tst) {
/* Computing MAX */
		r__3 = (r__1 = h__[k + (k - 1) * h_dim1], dabs(r__1)), r__4 = 
			(r__2 = h__[k - 1 + k * h_dim1], dabs(r__2));
		ab = dmax(r__3,r__4);
/* Computing MIN */
		r__3 = (r__1 = h__[k + (k - 1) * h_dim1], dabs(r__1)), r__4 = 
			(r__2 = h__[k - 1 + k * h_dim1], dabs(r__2));
		ba = dmin(r__3,r__4);
/* Computing MAX */
		r__3 = (r__1 = h__[k + k * h_dim1], dabs(r__1)), r__4 = (r__2 
			= h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1],
			 dabs(r__2));
		aa = dmax(r__3,r__4);
/* Computing MIN */
		r__3 = (r__1 = h__[k + k * h_dim1], dabs(r__1)), r__4 = (r__2 
			= h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1],
			 dabs(r__2));
		bb = dmin(r__3,r__4);
		s = aa + ab;
/* Computing MAX */
		r__1 = smlnum, r__2 = ulp * (bb * (aa / s));
		if (ba * (ab / s) <= dmax(r__1,r__2)) {
		    goto L40;
		}
	    }
/* L30: */
	}
L40:
	l = k;
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

	    h__[l + (l - 1) * h_dim1] = 0.f;
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

	if (l >= i__ - 1) {
	    goto L150;
	}
	++kdefl;

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

	if (! (*wantt)) {
	    i1 = l;
	    i2 = i__;
	}

	if (kdefl % 20 == 0) {

/*           Exceptional shift. */

	    s = (r__1 = h__[i__ + (i__ - 1) * h_dim1], dabs(r__1)) + (r__2 = 
		    h__[i__ - 1 + (i__ - 2) * h_dim1], dabs(r__2));
	    h11 = s * .75f + h__[i__ + i__ * h_dim1];
	    h12 = s * -.4375f;
	    h21 = s;
	    h22 = h11;
	} else if (kdefl % 10 == 0) {

/*           Exceptional shift. */

	    s = (r__1 = h__[l + 1 + l * h_dim1], dabs(r__1)) + (r__2 = h__[l 
		    + 2 + (l + 1) * h_dim1], dabs(r__2));
	    h11 = s * .75f + h__[l + l * h_dim1];
	    h12 = s * -.4375f;
	    h21 = s;
	    h22 = h11;
	} else {

/*           Prepare to use Francis' double shift */
/*           (i.e. 2nd degree generalized Rayleigh quotient) */

	    h11 = h__[i__ - 1 + (i__ - 1) * h_dim1];
	    h21 = h__[i__ + (i__ - 1) * h_dim1];
	    h12 = h__[i__ - 1 + i__ * h_dim1];
	    h22 = h__[i__ + i__ * h_dim1];
	}
	s = dabs(h11) + dabs(h12) + dabs(h21) + dabs(h22);
	if (s == 0.f) {
	    rt1r = 0.f;
	    rt1i = 0.f;
	    rt2r = 0.f;
	    rt2i = 0.f;
	} else {
	    h11 /= s;
	    h21 /= s;
	    h12 /= s;
	    h22 /= s;
	    tr = (h11 + h22) / 2.f;
	    det = (h11 - tr) * (h22 - tr) - h12 * h21;
	    rtdisc = sqrt((dabs(det)));
	    if (det >= 0.f) {

/*              ==== complex conjugate shifts ==== */

		rt1r = tr * s;
		rt2r = rt1r;
		rt1i = rtdisc * s;
		rt2i = -rt1i;
	    } else {

/*              ==== real shifts (use only one of them)  ==== */

		rt1r = tr + rtdisc;
		rt2r = tr - rtdisc;
		if ((r__1 = rt1r - h22, dabs(r__1)) <= (r__2 = rt2r - h22, 
			dabs(r__2))) {
		    rt1r *= s;
		    rt2r = rt1r;
		} else {
		    rt2r *= s;
		    rt1r = rt2r;
		}
		rt1i = 0.f;
		rt2i = 0.f;
	    }
	}

/*        Look for two consecutive small subdiagonal elements. */

	i__2 = l;
	for (m = i__ - 2; m >= i__2; --m) {
/*           Determine the effect of starting the double-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible.  (The following uses scaling to avoid */
/*           overflows and most underflows.) */

	    h21s = h__[m + 1 + m * h_dim1];
	    s = (r__1 = h__[m + m * h_dim1] - rt2r, dabs(r__1)) + dabs(rt2i) 
		    + dabs(h21s);
	    h21s = h__[m + 1 + m * h_dim1] / s;
	    v[0] = h21s * h__[m + (m + 1) * h_dim1] + (h__[m + m * h_dim1] - 
		    rt1r) * ((h__[m + m * h_dim1] - rt2r) / s) - rt1i * (rt2i 
		    / s);
	    v[1] = h21s * (h__[m + m * h_dim1] + h__[m + 1 + (m + 1) * h_dim1]
		     - rt1r - rt2r);
	    v[2] = h21s * h__[m + 2 + (m + 1) * h_dim1];
	    s = dabs(v[0]) + dabs(v[1]) + dabs(v[2]);
	    v[0] /= s;
	    v[1] /= s;
	    v[2] /= s;
	    if (m == l) {
		goto L60;
	    }
	    if ((r__1 = h__[m + (m - 1) * h_dim1], dabs(r__1)) * (dabs(v[1]) 
		    + dabs(v[2])) <= ulp * dabs(v[0]) * ((r__2 = h__[m - 1 + (
		    m - 1) * h_dim1], dabs(r__2)) + (r__3 = h__[m + m * 
		    h_dim1], dabs(r__3)) + (r__4 = h__[m + 1 + (m + 1) * 
		    h_dim1], dabs(r__4)))) {
		goto L60;
	    }
/* L50: */
	}
L60:

/*        Double-shift QR step */

	i__2 = i__ - 1;
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. NR is the order of G. */

/* Computing MIN */
	    i__3 = 3, i__4 = i__ - k + 1;
	    nr = min(i__3,i__4);
	    if (k > m) {
		scopy_(&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
	    }
	    slarfg_(&nr, v, &v[1], &c__1, &t1);
	    if (k > m) {
		h__[k + (k - 1) * h_dim1] = v[0];
		h__[k + 1 + (k - 1) * h_dim1] = 0.f;
		if (k < i__ - 1) {
		    h__[k + 2 + (k - 1) * h_dim1] = 0.f;
		}
	    } else if (m > l) {
/*               ==== Use the following instead of */
/*               .    H( K, K-1 ) = -H( K, K-1 ) to */
/*               .    avoid a bug when v(2) and v(3) */
/*               .    underflow. ==== */
		h__[k + (k - 1) * h_dim1] *= 1.f - t1;
	    }
	    v2 = v[1];
	    t2 = t1 * v2;
	    if (nr == 3) {
		v3 = v[2];
		t3 = t1 * v3;

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

		i__3 = i2;
		for (j = k; j <= i__3; ++j) {
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1] 
			    + v3 * h__[k + 2 + j * h_dim1];
		    h__[k + j * h_dim1] -= sum * t1;
		    h__[k + 1 + j * h_dim1] -= sum * t2;
		    h__[k + 2 + j * h_dim1] -= sum * t3;
/* L70: */
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

/* Computing MIN */
		i__4 = k + 3;
		i__3 = min(i__4,i__);
		for (j = i1; j <= i__3; ++j) {
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			     + v3 * h__[j + (k + 2) * h_dim1];
		    h__[j + k * h_dim1] -= sum * t1;
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
		    h__[j + (k + 2) * h_dim1] -= sum * t3;
/* L80: */
		}

		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

		    i__3 = *ihiz;
		    for (j = *iloz; j <= i__3; ++j) {
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1] + v3 * z__[j + (k + 2) * z_dim1];
			z__[j + k * z_dim1] -= sum * t1;
			z__[j + (k + 1) * z_dim1] -= sum * t2;
			z__[j + (k + 2) * z_dim1] -= sum * t3;
/* L90: */
		    }
		}
	    } else if (nr == 2) {

/*              Apply G from the left to transform the rows of the matrix */
/*              in columns K to I2. */

		i__3 = i2;
		for (j = k; j <= i__3; ++j) {
		    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1];
		    h__[k + j * h_dim1] -= sum * t1;
		    h__[k + 1 + j * h_dim1] -= sum * t2;
/* L100: */
		}

/*              Apply G from the right to transform the columns of the */
/*              matrix in rows I1 to min(K+3,I). */

		i__3 = i__;
		for (j = i1; j <= i__3; ++j) {
		    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1]
			    ;
		    h__[j + k * h_dim1] -= sum * t1;
		    h__[j + (k + 1) * h_dim1] -= sum * t2;
/* L110: */
		}

		if (*wantz) {

/*                 Accumulate transformations in the matrix Z */

		    i__3 = *ihiz;
		    for (j = *iloz; j <= i__3; ++j) {
			sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * 
				z_dim1];
			z__[j + k * z_dim1] -= sum * t1;
			z__[j + (k + 1) * z_dim1] -= sum * t2;
/* L120: */
		    }
		}
	    }
/* L130: */
	}

/* L140: */
    }

/*     Failure to converge in remaining number of iterations */

    *info = i__;
    return 0;

L150:

    if (l == i__) {

/*        H(I,I-1) is negligible: one eigenvalue has converged. */

	wr[i__] = h__[i__ + i__ * h_dim1];
	wi[i__] = 0.f;
    } else if (l == i__ - 1) {

/*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged. */

/*        Transform the 2-by-2 submatrix to standard Schur form, */
/*        and compute and store the eigenvalues. */

	slanv2_(&h__[i__ - 1 + (i__ - 1) * h_dim1], &h__[i__ - 1 + i__ * 
		h_dim1], &h__[i__ + (i__ - 1) * h_dim1], &h__[i__ + i__ * 
		h_dim1], &wr[i__ - 1], &wi[i__ - 1], &wr[i__], &wi[i__], &cs, 
		&sn);

	if (*wantt) {

/*           Apply the transformation to the rest of H. */

	    if (i2 > i__) {
		i__1 = i2 - i__;
		srot_(&i__1, &h__[i__ - 1 + (i__ + 1) * h_dim1], ldh, &h__[
			i__ + (i__ + 1) * h_dim1], ldh, &cs, &sn);
	    }
	    i__1 = i__ - i1 - 1;
	    srot_(&i__1, &h__[i1 + (i__ - 1) * h_dim1], &c__1, &h__[i1 + i__ *
		     h_dim1], &c__1, &cs, &sn);
	}
	if (*wantz) {

/*           Apply the transformation to Z. */

	    srot_(&nz, &z__[*iloz + (i__ - 1) * z_dim1], &c__1, &z__[*iloz + 
		    i__ * z_dim1], &c__1, &cs, &sn);
	}
    }
/*     reset deflation counter */
    kdefl = 0;

/*     return to start of the main loop with new value of I. */

    i__ = l - 1;
    goto L20;

L160:
    return 0;

/*     End of SLAHQR */

} /* slahqr_ */

/* > \brief \b SSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYCONV convert A given by TRF into L and D and vice-versa. */
/* > Get Non-diag elements of D (returned in workspace) and */
/* > apply or reverse permutation done in TRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* >          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1 */
/* >          or 2-by-2 block diagonal matrix D in LDLT. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssyconv_(char *uplo, char *way, integer *n, real *a, 
	integer *lda, integer *ipiv, real *e, integer *info, ftnlen uplo_len, 
	ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, j, ip;
    static real temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --e;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYCONV", &i__1, (ftnlen)7);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*      A is UPPER */

/*      Convert A (A is upper) */

/*        Convert VALUE */

	if (convert) {
	    i__ = *n;
	    e[1] = 0.f;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
		    e[i__ - 1] = 0.f;
		    a[i__ - 1 + i__ * a_dim1] = 0.f;
		    --i__;
		} else {
		    e[i__] = 0.f;
		}
		--i__;
	    }

/*        Convert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = temp;
/* L12: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
			    a[i__ - 1 + j * a_dim1] = temp;
/* L13: */
			}
		    }
		    --i__;
		}
		--i__;
	    }
	} else {

/*      Revert A (A is upper) */


/*        Revert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = temp;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    ++i__;
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
			    a[i__ - 1 + j * a_dim1] = temp;
			}
		    }
		}
		++i__;
	    }

/*        Revert VALUE */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
		    --i__;
		}
		--i__;
	    }
	}
    } else {

/*      A is LOWER */

	if (convert) {

/*      Convert A (A is lower) */


/*        Convert VALUE */

	    i__ = 1;
	    e[*n] = 0.f;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
		    e[i__ + 1] = 0.f;
		    a[i__ + 1 + i__ * a_dim1] = 0.f;
		    ++i__;
		} else {
		    e[i__] = 0.f;
		}
		++i__;
	    }

/*        Convert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = temp;
/* L22: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
			    a[i__ + 1 + j * a_dim1] = temp;
/* L23: */
			}
		    }
		    ++i__;
		}
		++i__;
	    }
	} else {

/*      Revert A (A is lower) */


/*        Revert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = temp;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    --i__;
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    temp = a[i__ + 1 + j * a_dim1];
			    a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
			    a[ip + j * a_dim1] = temp;
			}
		    }
		}
		--i__;
	    }

/*        Revert VALUE */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
		    ++i__;
		}
		++i__;
	    }
	}
    }
    return 0;

/*     End of SSYCONV */

} /* ssyconv_ */

/* > \brief \b SSYCONVF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYCONVF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv
f.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv
f.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv
f.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > SSYCONVF converts the factorization output format used in */
/* > SSYTRF provided on entry in parameter A into the factorization */
/* > output format used in SSYTRF_RK (or SSYTRF_BK) that is stored */
/* > on exit in parameters A and E. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in SSYTRF into */
/* > the format used in SSYTRF_RK (or SSYTRF_BK). */
/* > */
/* > If parameter WAY = 'R': */
/* > SSYCONVF performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in SSYTRF_RK */
/* > (or SSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in SSYTRF that is stored */
/* > on exit in parameter A. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in SSYTRF_RK */
/* > (or SSYTRF_BK) into the format used in SSYTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          SSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          SSYTRF_RK or SSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          SSYTRF_RK or SSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          SSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF. */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF_RK */
/* >          ( or SSYTRF_BK). */
/* > */
/* >          1) If WAY ='R': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF_RK */
/* >          ( or SSYTRF_BK). */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup singleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int ssyconvf_(char *uplo, char *way, integer *n, real *a, 
	integer *lda, real *e, integer *ipiv, integer *info, ftnlen uplo_len, 
	ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, ip;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYCONVF", &i__1, (ftnlen)8);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1] = 0.f;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
		    e[i__ - 1] = 0.f;
		    a[i__ - 1 + i__ * a_dim1] = 0.f;
		    --i__;
		} else {
		    e[i__] = 0.f;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i-1 and IPIV(i-1), */
/*                 so this should be recorded in two consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ - 1];

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    e[*n] = 0.f;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
		    e[i__ + 1] = 0.f;
		    a[i__ + 1 + i__ * a_dim1] = 0.f;
		    ++i__;
		} else {
		    e[i__] = 0.f;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where k increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i+1 and IPIV(i+1), */
/*                 so this should be recorded in consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ + 1];

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of SSYCONVF */

} /* ssyconvf_ */

/* > \brief \b SSYCONVF_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYCONVF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv
f_rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv
f_rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv
f_rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > SSYCONVF_ROOK converts the factorization output format used in */
/* > SSYTRF_ROOK provided on entry in parameter A into the factorization */
/* > output format used in SSYTRF_RK (or SSYTRF_BK) that is stored */
/* > on exit in parameters A and E. IPIV format for SSYTRF_ROOK and */
/* > SSYTRF_RK (or SSYTRF_BK) is the same and is not converted. */
/* > */
/* > If parameter WAY = 'R': */
/* > SSYCONVF_ROOK performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in SSYTRF_RK */
/* > (or SSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in SSYTRF_ROOK that is stored */
/* > on exit in parameter A. IPIV format for SSYTRF_ROOK and */
/* > SSYTRF_RK (or SSYTRF_BK) is the same and is not converted. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          SSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          SSYTRF_RK or SSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          SSYTRF_RK or SSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          SSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D as determined: */
/* >          1) by SSYTRF_ROOK, if WAY ='C'; */
/* >          2) by SSYTRF_RK (or SSYTRF_BK), if WAY ='R'. */
/* >          The IPIV format is the same for all these routines. */
/* > */
/* >          On exit, is not changed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup singleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int ssyconvf_rook__(char *uplo, char *way, integer *n, real *
	a, integer *lda, real *e, integer *ipiv, integer *info, ftnlen 
	uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, ip, ip2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYCONVF_ROOK", &i__1, (ftnlen)13);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1] = 0.f;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
		    e[i__ - 1] = 0.f;
		    a[i__ - 1 + i__ * a_dim1] = 0.f;
		    --i__;
		} else {
		    e[i__] = 0.f;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1) */
/*                 in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip2 + (i__ + 1) * a_dim1], lda);
			}
		    }
		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i) */
/*                 in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[ip2 + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = *n - i__;
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    e[*n] = 0.f;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
		    e[i__ + 1] = 0.f;
		    a[i__ + 1 + i__ * a_dim1] = 0.f;
		    ++i__;
		} else {
		    e[i__] = 0.f;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1) */
/*                 in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip2 + 
				    a_dim1], lda);
			}
		    }
		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i) */
/*                 in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[ip2 + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = i__ - 1;
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of SSYCONVF_ROOK */

} /* ssyconvf_rook__ */

/* > \brief \b STPMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpmqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpmqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpmqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/*      $          WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPMQRT applies a real orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" real block reflector H to a general */
/* > real matrix C, which consists of two blocks A and B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q^T from the Left; */
/* >          = 'R': apply Q or Q^T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q^T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The order of the trapezoidal part of V. */
/* >          K >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CTPQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is REAL array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CTPQRT in B.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDV >= max(1,M); */
/* >          if SIDE = 'R', LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CTPQRT, stored as a NB-by-K matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension */
/* >          (LDA,N) if SIDE = 'L' or */
/* >          (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDC >= max(1,K); */
/* >          If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. */
/* >          LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The columns of the pentagonal matrix V contain the elementary reflectors */
/* >  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a */
/* >  trapezoidal block V2: */
/* > */
/* >        V = [V1] */
/* >            [V2]. */
/* > */
/* >  The size of the trapezoidal block V2 is determined by the parameter L, */
/* >  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L */
/* >  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular; */
/* >  if L=0, there is no trapezoidal block, hence V = V1 is rectangular. */
/* > */
/* >  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K. */
/* >                      [B] */
/* > */
/* >  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K. */
/* > */
/* >  The real orthogonal matrix Q is formed from V and T. */
/* > */
/* >  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C. */
/* > */
/* >  If TRANS='T' and SIDE='L', C is on exit replaced with Q^T * C. */
/* > */
/* >  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* >  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q^T. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stpmqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *l, integer *nb, real *v, integer *ldv, real *t, 
	integer *ldt, real *a, integer *lda, real *b, integer *ldb, real *
	work, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, mb, kf, ldaq;
    static logical left, tran;
    static integer ldvq;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), stprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, real *, integer *, real *, integer *, real *, integer *
	    , real *, integer *, real *, integer *);
    static logical notran;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldvq = max(1,*m);
	ldaq = max(1,*k);
    } else if (right) {
	ldvq = max(1,*n);
	ldaq = max(1,*m);
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0) {
	*info = -5;
    } else if (*l < 0 || *l > *k) {
	*info = -6;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -7;
    } else if (*ldv < ldvq) {
	*info = -9;
    } else if (*ldt < *nb) {
	*info = -11;
    } else if (*lda < ldaq) {
	*info = -13;
    } else if (*ldb < max(1,*m)) {
	*info = -15;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("STPMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *m - *l + i__ + ib - 1;
	    mb = min(i__3,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    stprfb_("L", "T", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *n - *l + i__ + ib - 1;
	    mb = min(i__3,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    stprfb_("R", "N", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *m - *l + i__ + ib - 1;
	    mb = min(i__2,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    stprfb_("L", "N", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *n - *l + i__ + ib - 1;
	    mb = min(i__2,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    stprfb_("R", "T", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    }

    return 0;

/*     End of STPMQRT */

} /* stpmqrt_ */

/* > \brief \b STPQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STPQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPQRT computes a blocked QR factorization of a real */
/* > "triangular-pentagonal" matrix C, which is composed of a */
/* > triangular block A and pentagonal block B, using the compact */
/* > WY representation for Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B, and the order of the */
/* >          triangular matrix A. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of rows of the upper trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the upper triangular N-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,N) */
/* >          On entry, the pentagonal M-by-N matrix B.  The first M-L rows */
/* >          are rectangular, and the last L rows are upper trapezoidal. */
/* >          On exit, B contains the pentagonal matrix V.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,N) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a (N+M)-by-N matrix */
/* > */
/* >               C = [ A ] */
/* >                   [ B ] */
/* > */
/* >  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* >  upper trapezoidal matrix B2: */
/* > */
/* >               B = [ B1 ]  <- (M-L)-by-N rectangular */
/* >                   [ B2 ]  <-     L-by-N upper trapezoidal. */
/* > */
/* >  The upper trapezoidal matrix B2 consists of the first L rows of a */
/* >  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is upper triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* >               C = [ A ]  <- upper triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* > */
/* >               W = [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* > */
/* >               V = [ V1 ] <- (M-L)-by-N rectangular */
/* >                   [ V2 ] <-     L-by-N upper trapezoidal. */
/* > */
/* >  The columns of V represent the vectors which define the H(i)'s. */
/* > */
/* >  The number of blocks is B = ceiling(N/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* >               T = [T1 T2 ... TB]. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stpqrt_(integer *m, integer *n, integer *l, integer *nb, 
	real *a, integer *lda, real *b, integer *ldb, real *t, integer *ldt, 
	real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, ib, lb, mb, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), stprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, real *, integer *, real *, integer *, real *, integer *
	    , real *, integer *, real *, integer *), stpqrt2_(integer *, integer *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*l < 0 || *l > min(*m,*n) && min(*m,*n) >= 0) {
	*info = -3;
    } else if (*nb < 1 || *nb > *n && *n > 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*m)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("STPQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = *n;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*     Compute the QR factorization of the current block */

/* Computing MIN */
	i__3 = *n - i__ + 1;
	ib = min(i__3,*nb);
/* Computing MIN */
	i__3 = *m - *l + i__ + ib - 1;
	mb = min(i__3,*m);
	if (i__ >= *l) {
	    lb = 0;
	} else {
	    lb = mb - *m + *l - i__ + 1;
	}

	stpqrt2_(&mb, &ib, &lb, &a[i__ + i__ * a_dim1], lda, &b[i__ * b_dim1 
		+ 1], ldb, &t[i__ * t_dim1 + 1], ldt, &iinfo);

/*     Update by applying H^H to B(:,I+IB:N) from the left */

	if (i__ + ib <= *n) {
	    i__3 = *n - i__ - ib + 1;
	    stprfb_("L", "T", "F", "C", &mb, &i__3, &ib, &lb, &b[i__ * b_dim1 
		    + 1], ldb, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + ib) 
		    * a_dim1], lda, &b[(i__ + ib) * b_dim1 + 1], ldb, &work[1]
		    , &ib);
	}
    }
    return 0;

/*     End of STPQRT */

} /* stpqrt_ */

/* > \brief \b ZGEMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgemqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgemqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgemqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEMQRT overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q C            C Q */
/* > TRANS = 'C':    Q**H C            C Q**H */
/* > */
/* > where Q is a complex orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* >       Q = H(1) H(2) . . . H(K) = I - V T V**H */
/* > */
/* > generated using the compact WY representation as returned by ZGEQRT. */
/* > */
/* > Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Conjugate transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CGEQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGEQRT in the first K columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CGEQRT, stored as a NB-by-N matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q C, Q**H C, C Q**H or C Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgemqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *nb, doublecomplex *v, integer *ldv, 
	doublecomplex *t, integer *ldt, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, q, ib, kf;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlarfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static logical notran;
    static integer ldwork;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldwork = max(1,*n);
	q = *m;
    } else if (right) {
	ldwork = max(1,*m);
	q = *n;
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > q) {
	*info = -5;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -6;
    } else if (*ldv < max(1,q)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    } else if (*ldc < max(1,*m)) {
	*info = -12;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGEMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *m - i__ + 1;
	    zlarfb_("L", "C", "F", "C", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
	    i__3 = *n - i__ + 1;
	    zlarfb_("R", "N", "F", "C", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *m - i__ + 1;
	    zlarfb_("L", "N", "F", "C", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    i__2 = *n - i__ + 1;
	    zlarfb_("R", "C", "F", "C", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
	}

    }

    return 0;

/*     End of ZGEMQRT */

} /* zgemqrt_ */

/* > \brief \b ZGEQRFP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrfp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrfp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrfp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQR2P computes a QR factorization of a complex M-by-N matrix A: */
/* > */
/* >    A = Q * ( R ), */
/* >            ( 0 ) */
/* > */
/* > where: */
/* > */
/* >    Q is a M-by-M orthogonal matrix; */
/* >    R is an upper-triangular N-by-N matrix with nonnegative diagonal */
/* >    entries; */
/* >    0 is a (M-N)-by-N zero matrix, if M > N. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n). The diagonal entries of R */
/* >          are real and nonnegative; The elements below the diagonal, */
/* >          with the array TAU, represent the unitary matrix Q as a */
/* >          product of min(m,n) elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,N). */
/* >          For optimum performance LWORK >= N*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16GEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > */
/* > See Lapack Working Note 203 for details */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqrfp_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int zlarft_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zgeqr2p_(integer *, integer *, doublecomplex *
	    , integer *, doublecomplex *, doublecomplex *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
    lwkopt = *n * nb;
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGEQRFP", &i__1, (ftnlen)7);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	work[1].r = 1., work[1].i = 0.;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

	    i__3 = *m - i__ + 1;
	    zgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *m - i__ + 1;
		zlarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**H to A(i:m,i+ib:n) from the left */

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		zlarfb_("Left", "Conjugate transpose", "Forward", "Columnwise"
			, &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &
			work[1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, 
			&work[ib + 1], &ldwork, (ftnlen)4, (ftnlen)19, (
			ftnlen)7, (ftnlen)10);
	    }
/* L10: */
	}
    } else {
	i__ = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	zgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		1], &iinfo);
    }

    work[1].r = (doublereal) iws, work[1].i = 0.;
    return 0;

/*     End of ZGEQRFP */

} /* zgeqrfp_ */

/* > \brief \b ZGEQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDT, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQRT computes a blocked QR factorization of a complex M-by-N matrix A */
/* > using the compact WY representation of Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if M >= N); the elements below the diagonal */
/* >          are the columns of V. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,MIN(M,N)) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16GEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix V stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* >               V = (  1       ) */
/* >                   ( v1  1    ) */
/* >                   ( v1 v2  1 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  where the vi's represent the vectors which define H(i), which are returned */
/* >  in the matrix A.  The 1's along the diagonal of V are not stored in A. */
/* > */
/* >  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-K matrix T as */
/* > */
/* >               T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqrt_(integer *m, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *t, integer *ldt, 
	doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ib, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlarfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), zgeqrt2_(integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *)
	    , zgeqrt3_(integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nb < 1 || *nb > min(*m,*n) && min(*m,*n) > 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*ldt < *nb) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGEQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	return 0;
    }

/*     Blocked loop of length K */

    i__1 = k;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	i__3 = k - i__ + 1;
	ib = min(i__3,*nb);

/*     Compute the QR factorization of the current block A(I:M,I:I+IB-1) */

	if (TRUE_) {
	    i__3 = *m - i__ + 1;
	    zgeqrt3_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	} else {
	    i__3 = *m - i__ + 1;
	    zgeqrt2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 
		    + 1], ldt, &iinfo);
	}
	if (i__ + ib <= *n) {

/*     Update by applying H**H to A(I:M,I+IB:N) from the left */

	    i__3 = *m - i__ + 1;
	    i__4 = *n - i__ - ib + 1;
	    i__5 = *n - i__ - ib + 1;
	    zlarfb_("L", "C", "F", "C", &i__3, &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + 
		    ib) * a_dim1], lda, &work[1], &i__5, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
	}
    }
    return 0;

/*     End of ZGEQRT */

} /* zgeqrt_ */

/* > \brief \b ZLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using th
e double-shift/single-shift QR algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAHQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, */
/*                          IHIZ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZLAHQR is an auxiliary routine called by CHSEQR to update the */
/* >    eigenvalues and Schur decomposition already computed by CHSEQR, by */
/* >    dealing with the Hessenberg submatrix in rows and columns ILO to */
/* >    IHI. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >          = .TRUE. : the full Schur form T is required; */
/* >          = .FALSE.: only eigenvalues are required. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          = .TRUE. : the matrix of Schur vectors Z is required; */
/* >          = .FALSE.: Schur vectors are not required. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          It is assumed that H is already upper triangular in rows and */
/* >          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1). */
/* >          ZLAHQR works primarily with the Hessenberg submatrix in rows */
/* >          and columns ILO to IHI, but applies transformations to all of */
/* >          H if WANTT is .TRUE.. */
/* >          1 <= ILO <= max(1,IHI); IHI <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is COMPLEX*16 array, dimension (LDH,N) */
/* >          On entry, the upper Hessenberg matrix H. */
/* >          On exit, if INFO is zero and if WANTT is .TRUE., then H */
/* >          is upper triangular in rows and columns ILO:IHI.  If INFO */
/* >          is zero and if WANTT is .FALSE., then the contents of H */
/* >          are unspecified on exit.  The output state of H in case */
/* >          INF is positive is below under the description of INFO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX*16 array, dimension (N) */
/* >          The computed eigenvalues ILO to IHI are stored in the */
/* >          corresponding elements of W. If WANTT is .TRUE., the */
/* >          eigenvalues are stored in the same order as on the diagonal */
/* >          of the Schur form returned in H, with W(i) = H(i,i). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* >          ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* >          IHIZ is INTEGER */
/* >          Specify the rows of Z to which transformations must be */
/* >          applied if WANTZ is .TRUE.. */
/* >          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ,N) */
/* >          If WANTZ is .TRUE., on entry Z must contain the current */
/* >          matrix Z of transformations accumulated by CHSEQR, and on */
/* >          exit Z has been updated; transformations are applied only to */
/* >          the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
/* >          If WANTZ is .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:   successful exit */
/* >           > 0:   if INFO = i, ZLAHQR failed to compute all the */
/* >                  eigenvalues ILO to IHI in a total of 30 iterations */
/* >                  per eigenvalue; elements i+1:ihi of W contain */
/* >                  those eigenvalues which have been successfully */
/* >                  computed. */
/* > */
/* >                  If INFO > 0 and WANTT is .FALSE., then on exit, */
/* >                  the remaining unconverged eigenvalues are the */
/* >                  eigenvalues of the upper Hessenberg matrix */
/* >                  rows and columns ILO through INFO of the final, */
/* >                  output value of H. */
/* > */
/* >                  If INFO > 0 and WANTT is .TRUE., then on exit */
/* >          (*)       (initial value of H)*U  = U*(final value of H) */
/* >                  where U is an orthogonal matrix.    The final */
/* >                  value of H is upper Hessenberg and triangular in */
/* >                  rows and columns INFO+1 through IHI. */
/* > */
/* >                  If INFO > 0 and WANTZ is .TRUE., then on exit */
/* >                      (final value of Z)  = (initial value of Z)*U */
/* >                  where U is the orthogonal matrix in (*) */
/* >                  (regardless of the value of WANTT.) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >     02-96 Based on modifications by */
/* >     David Day, Sandia National Laboratory, USA */
/* > */
/* >     12-04 Further modifications by */
/* >     Ralph Byers, University of Kansas, USA */
/* >     This is a modified version of ZLAHQR from LAPACK version 3.0. */
/* >     It is (1) more robust against overflow and underflow and */
/* >     (2) adopts the more conservative Ahues & Tisseur stopping */
/* >     criterion (LAWN 122, 1997). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    void z_sqrt(doublecomplex *, doublecomplex *), pow_zi(doublecomplex *, 
	    doublecomplex *, integer *);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal s;
    static doublecomplex t, u, v[2], x, y;
    static integer i1, i2;
    static doublecomplex t1;
    static doublereal t2;
    static doublecomplex v2;
    static doublereal aa, ab, ba, bb, h10;
    static doublecomplex h11;
    static doublereal h21;
    static doublecomplex h22, sc;
    static integer nh, nz;
    static doublereal sx;
    static integer jhi;
    static doublecomplex h11s;
    static integer jlo, its;
    static doublereal ulp;
    static doublecomplex sum;
    static doublereal tst;
    static doublecomplex temp;
    static integer kdefl;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer itmax;
    static doublereal rtemp;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin, safmax;
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ========================================================= */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }
    if (*ilo == *ihi) {
	i__1 = *ilo;
	i__2 = *ilo + *ilo * h_dim1;
	w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;
	return 0;
    }

/*     ==== clear out the trash ==== */
    i__1 = *ihi - 3;
    for (j = *ilo; j <= i__1; ++j) {
	i__2 = j + 2 + j * h_dim1;
	h__[i__2].r = 0., h__[i__2].i = 0.;
	i__2 = j + 3 + j * h_dim1;
	h__[i__2].r = 0., h__[i__2].i = 0.;
/* L10: */
    }
    if (*ilo <= *ihi - 2) {
	i__1 = *ihi + (*ihi - 2) * h_dim1;
	h__[i__1].r = 0., h__[i__1].i = 0.;
    }
/*     ==== ensure that subdiagonal entries are real ==== */
    if (*wantt) {
	jlo = 1;
	jhi = *n;
    } else {
	jlo = *ilo;
	jhi = *ihi;
    }
    i__1 = *ihi;
    for (i__ = *ilo + 1; i__ <= i__1; ++i__) {
	if (d_imag(&h__[i__ + (i__ - 1) * h_dim1]) != 0.) {
/*           ==== The following redundant normalization */
/*           .    avoids problems with both gradual and */
/*           .    sudden underflow in ABS(H(I,I-1)) ==== */
	    i__2 = i__ + (i__ - 1) * h_dim1;
	    i__3 = i__ + (i__ - 1) * h_dim1;
	    d__3 = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[i__ 
		    + (i__ - 1) * h_dim1]), abs(d__2));
	    z__1.r = h__[i__2].r / d__3, z__1.i = h__[i__2].i / d__3;
	    sc.r = z__1.r, sc.i = z__1.i;
	    d_cnjg(&z__2, &sc);
	    d__1 = z_abs(&sc);
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    sc.r = z__1.r, sc.i = z__1.i;
	    i__2 = i__ + (i__ - 1) * h_dim1;
	    d__1 = z_abs(&h__[i__ + (i__ - 1) * h_dim1]);
	    h__[i__2].r = d__1, h__[i__2].i = 0.;
	    i__2 = jhi - i__ + 1;
	    zscal_(&i__2, &sc, &h__[i__ + i__ * h_dim1], ldh);
/* Computing MIN */
	    i__3 = jhi, i__4 = i__ + 1;
	    i__2 = min(i__3,i__4) - jlo + 1;
	    d_cnjg(&z__1, &sc);
	    zscal_(&i__2, &z__1, &h__[jlo + i__ * h_dim1], &c__1);
	    if (*wantz) {
		i__2 = *ihiz - *iloz + 1;
		d_cnjg(&z__1, &sc);
		zscal_(&i__2, &z__1, &z__[*iloz + i__ * z_dim1], &c__1);
	    }
	}
/* L20: */
    }

    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion. */

    safmin = dlamch_("SAFE MINIMUM", (ftnlen)12);
    safmax = 1. / safmin;
    dlabad_(&safmin, &safmax);
    ulp = dlamch_("PRECISION", (ftnlen)9);
    smlnum = safmin * ((doublereal) nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

    if (*wantt) {
	i1 = 1;
	i2 = *n;
    }

/*     ITMAX is the total number of QR iterations allowed. */

    itmax = max(10,nh) * 30;

/*     KDEFL counts the number of iterations since a deflation */

    kdefl = 0;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or */
/*     H(L,L-1) is negligible so that the matrix splits. */

    i__ = *ihi;
L30:
    if (i__ < *ilo) {
	goto L150;
    }

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

    l = *ilo;
    i__1 = itmax;
    for (its = 0; its <= i__1; ++its) {

/*        Look for a single small subdiagonal element. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {
	    i__3 = k + (k - 1) * h_dim1;
	    if ((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[k + (k 
		    - 1) * h_dim1]), abs(d__2)) <= smlnum) {
		goto L50;
	    }
	    i__3 = k - 1 + (k - 1) * h_dim1;
	    i__4 = k + k * h_dim1;
	    tst = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[k - 1 
		    + (k - 1) * h_dim1]), abs(d__2)) + ((d__3 = h__[i__4].r, 
		    abs(d__3)) + (d__4 = d_imag(&h__[k + k * h_dim1]), abs(
		    d__4)));
	    if (tst == 0.) {
		if (k - 2 >= *ilo) {
		    i__3 = k - 1 + (k - 2) * h_dim1;
		    tst += (d__1 = h__[i__3].r, abs(d__1));
		}
		if (k + 1 <= *ihi) {
		    i__3 = k + 1 + k * h_dim1;
		    tst += (d__1 = h__[i__3].r, abs(d__1));
		}
	    }
/*           ==== The following is a conservative small subdiagonal */
/*           .    deflation criterion due to Ahues & Tisseur (LAWN 122, */
/*           .    1997). It has better mathematical foundation and */
/*           .    improves accuracy in some examples.  ==== */
	    i__3 = k + (k - 1) * h_dim1;
	    if ((d__1 = h__[i__3].r, abs(d__1)) <= ulp * tst) {
/* Computing MAX */
		i__3 = k + (k - 1) * h_dim1;
		i__4 = k - 1 + k * h_dim1;
		d__5 = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + (k - 1) * h_dim1]), abs(d__2)), d__6 = (d__3 = 
			h__[i__4].r, abs(d__3)) + (d__4 = d_imag(&h__[k - 1 + 
			k * h_dim1]), abs(d__4));
		ab = max(d__5,d__6);
/* Computing MIN */
		i__3 = k + (k - 1) * h_dim1;
		i__4 = k - 1 + k * h_dim1;
		d__5 = (d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + (k - 1) * h_dim1]), abs(d__2)), d__6 = (d__3 = 
			h__[i__4].r, abs(d__3)) + (d__4 = d_imag(&h__[k - 1 + 
			k * h_dim1]), abs(d__4));
		ba = min(d__5,d__6);
		i__3 = k - 1 + (k - 1) * h_dim1;
		i__4 = k + k * h_dim1;
		z__2.r = h__[i__3].r - h__[i__4].r, z__2.i = h__[i__3].i - 
			h__[i__4].i;
		z__1.r = z__2.r, z__1.i = z__2.i;
/* Computing MAX */
		i__5 = k + k * h_dim1;
		d__5 = (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + k * h_dim1]), abs(d__2)), d__6 = (d__3 = z__1.r, 
			abs(d__3)) + (d__4 = d_imag(&z__1), abs(d__4));
		aa = max(d__5,d__6);
		i__3 = k - 1 + (k - 1) * h_dim1;
		i__4 = k + k * h_dim1;
		z__2.r = h__[i__3].r - h__[i__4].r, z__2.i = h__[i__3].i - 
			h__[i__4].i;
		z__1.r = z__2.r, z__1.i = z__2.i;
/* Computing MIN */
		i__5 = k + k * h_dim1;
		d__5 = (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = d_imag(&h__[
			k + k * h_dim1]), abs(d__2)), d__6 = (d__3 = z__1.r, 
			abs(d__3)) + (d__4 = d_imag(&z__1), abs(d__4));
		bb = min(d__5,d__6);
		s = aa + ab;
/* Computing MAX */
		d__1 = smlnum, d__2 = ulp * (bb * (aa / s));
		if (ba * (ab / s) <= max(d__1,d__2)) {
		    goto L50;
		}
	    }
/* L40: */
	}
L50:
	l = k;
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

	    i__2 = l + (l - 1) * h_dim1;
	    h__[i__2].r = 0., h__[i__2].i = 0.;
	}

/*        Exit from loop if a submatrix of order 1 has split off. */

	if (l >= i__) {
	    goto L140;
	}
	++kdefl;

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

	if (! (*wantt)) {
	    i1 = l;
	    i2 = i__;
	}

	if (kdefl % 20 == 0) {

/*           Exceptional shift. */

	    i__2 = i__ + (i__ - 1) * h_dim1;
	    s = (d__1 = h__[i__2].r, abs(d__1)) * .75;
	    i__2 = i__ + i__ * h_dim1;
	    z__1.r = s + h__[i__2].r, z__1.i = h__[i__2].i;
	    t.r = z__1.r, t.i = z__1.i;
	} else if (kdefl % 10 == 0) {

/*           Exceptional shift. */

	    i__2 = l + 1 + l * h_dim1;
	    s = (d__1 = h__[i__2].r, abs(d__1)) * .75;
	    i__2 = l + l * h_dim1;
	    z__1.r = s + h__[i__2].r, z__1.i = h__[i__2].i;
	    t.r = z__1.r, t.i = z__1.i;
	} else {

/*           Wilkinson's shift. */

	    i__2 = i__ + i__ * h_dim1;
	    t.r = h__[i__2].r, t.i = h__[i__2].i;
	    z_sqrt(&z__2, &h__[i__ - 1 + i__ * h_dim1]);
	    z_sqrt(&z__3, &h__[i__ + (i__ - 1) * h_dim1]);
	    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
		    z__3.i + z__2.i * z__3.r;
	    u.r = z__1.r, u.i = z__1.i;
	    s = (d__1 = u.r, abs(d__1)) + (d__2 = d_imag(&u), abs(d__2));
	    if (s != 0.) {
		i__2 = i__ - 1 + (i__ - 1) * h_dim1;
		z__2.r = h__[i__2].r - t.r, z__2.i = h__[i__2].i - t.i;
		z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
		x.r = z__1.r, x.i = z__1.i;
		sx = (d__1 = x.r, abs(d__1)) + (d__2 = d_imag(&x), abs(d__2));
/* Computing MAX */
		d__3 = s, d__4 = (d__1 = x.r, abs(d__1)) + (d__2 = d_imag(&x),
			 abs(d__2));
		s = max(d__3,d__4);
		z__5.r = x.r / s, z__5.i = x.i / s;
		pow_zi(&z__4, &z__5, &c__2);
		z__7.r = u.r / s, z__7.i = u.i / s;
		pow_zi(&z__6, &z__7, &c__2);
		z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
		z_sqrt(&z__2, &z__3);
		z__1.r = s * z__2.r, z__1.i = s * z__2.i;
		y.r = z__1.r, y.i = z__1.i;
		if (sx > 0.) {
		    z__1.r = x.r / sx, z__1.i = x.i / sx;
		    z__2.r = x.r / sx, z__2.i = x.i / sx;
		    if (z__1.r * y.r + d_imag(&z__2) * d_imag(&y) < 0.) {
			z__3.r = -y.r, z__3.i = -y.i;
			y.r = z__3.r, y.i = z__3.i;
		    }
		}
		z__4.r = x.r + y.r, z__4.i = x.i + y.i;
		zladiv_(&z__3, &u, &z__4);
		z__2.r = u.r * z__3.r - u.i * z__3.i, z__2.i = u.r * z__3.i + 
			u.i * z__3.r;
		z__1.r = t.r - z__2.r, z__1.i = t.i - z__2.i;
		t.r = z__1.r, t.i = z__1.i;
	    }
	}

/*        Look for two consecutive small subdiagonal elements. */

	i__2 = l + 1;
	for (m = i__ - 1; m >= i__2; --m) {

/*           Determine the effect of starting the single-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible. */

	    i__3 = m + m * h_dim1;
	    h11.r = h__[i__3].r, h11.i = h__[i__3].i;
	    i__3 = m + 1 + (m + 1) * h_dim1;
	    h22.r = h__[i__3].r, h22.i = h__[i__3].i;
	    z__1.r = h11.r - t.r, z__1.i = h11.i - t.i;
	    h11s.r = z__1.r, h11s.i = z__1.i;
	    i__3 = m + 1 + m * h_dim1;
	    h21 = h__[i__3].r;
	    s = (d__1 = h11s.r, abs(d__1)) + (d__2 = d_imag(&h11s), abs(d__2))
		     + abs(h21);
	    z__1.r = h11s.r / s, z__1.i = h11s.i / s;
	    h11s.r = z__1.r, h11s.i = z__1.i;
	    h21 /= s;
	    v[0].r = h11s.r, v[0].i = h11s.i;
	    v[1].r = h21, v[1].i = 0.;
	    i__3 = m + (m - 1) * h_dim1;
	    h10 = h__[i__3].r;
	    if (abs(h10) * abs(h21) <= ulp * (((d__1 = h11s.r, abs(d__1)) + (
		    d__2 = d_imag(&h11s), abs(d__2))) * ((d__3 = h11.r, abs(
		    d__3)) + (d__4 = d_imag(&h11), abs(d__4)) + ((d__5 = 
		    h22.r, abs(d__5)) + (d__6 = d_imag(&h22), abs(d__6)))))) {
		goto L70;
	    }
/* L60: */
	}
	i__2 = l + l * h_dim1;
	h11.r = h__[i__2].r, h11.i = h__[i__2].i;
	i__2 = l + 1 + (l + 1) * h_dim1;
	h22.r = h__[i__2].r, h22.i = h__[i__2].i;
	z__1.r = h11.r - t.r, z__1.i = h11.i - t.i;
	h11s.r = z__1.r, h11s.i = z__1.i;
	i__2 = l + 1 + l * h_dim1;
	h21 = h__[i__2].r;
	s = (d__1 = h11s.r, abs(d__1)) + (d__2 = d_imag(&h11s), abs(d__2)) + 
		abs(h21);
	z__1.r = h11s.r / s, z__1.i = h11s.i / s;
	h11s.r = z__1.r, h11s.i = z__1.i;
	h21 /= s;
	v[0].r = h11s.r, v[0].i = h11s.i;
	v[1].r = h21, v[1].i = 0.;
L70:

/*        Single-shift QR step */

	i__2 = i__ - 1;
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. */

/*           V(2) is always real before the call to ZLARFG, and hence */
/*           after the call T2 ( = T1*V(2) ) is also real. */

	    if (k > m) {
		zcopy_(&c__2, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
	    }
	    zlarfg_(&c__2, v, &v[1], &c__1, &t1);
	    if (k > m) {
		i__3 = k + (k - 1) * h_dim1;
		h__[i__3].r = v[0].r, h__[i__3].i = v[0].i;
		i__3 = k + 1 + (k - 1) * h_dim1;
		h__[i__3].r = 0., h__[i__3].i = 0.;
	    }
	    v2.r = v[1].r, v2.i = v[1].i;
	    z__1.r = t1.r * v2.r - t1.i * v2.i, z__1.i = t1.r * v2.i + t1.i * 
		    v2.r;
	    t2 = z__1.r;

/*           Apply G from the left to transform the rows of the matrix */
/*           in columns K to I2. */

	    i__3 = i2;
	    for (j = k; j <= i__3; ++j) {
		d_cnjg(&z__3, &t1);
		i__4 = k + j * h_dim1;
		z__2.r = z__3.r * h__[i__4].r - z__3.i * h__[i__4].i, z__2.i =
			 z__3.r * h__[i__4].i + z__3.i * h__[i__4].r;
		i__5 = k + 1 + j * h_dim1;
		z__4.r = t2 * h__[i__5].r, z__4.i = t2 * h__[i__5].i;
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
		sum.r = z__1.r, sum.i = z__1.i;
		i__4 = k + j * h_dim1;
		i__5 = k + j * h_dim1;
		z__1.r = h__[i__5].r - sum.r, z__1.i = h__[i__5].i - sum.i;
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
		i__4 = k + 1 + j * h_dim1;
		i__5 = k + 1 + j * h_dim1;
		z__2.r = sum.r * v2.r - sum.i * v2.i, z__2.i = sum.r * v2.i + 
			sum.i * v2.r;
		z__1.r = h__[i__5].r - z__2.r, z__1.i = h__[i__5].i - z__2.i;
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
/* L80: */
	    }

/*           Apply G from the right to transform the columns of the */
/*           matrix in rows I1 to min(K+2,I). */

/* Computing MIN */
	    i__4 = k + 2;
	    i__3 = min(i__4,i__);
	    for (j = i1; j <= i__3; ++j) {
		i__4 = j + k * h_dim1;
		z__2.r = t1.r * h__[i__4].r - t1.i * h__[i__4].i, z__2.i = 
			t1.r * h__[i__4].i + t1.i * h__[i__4].r;
		i__5 = j + (k + 1) * h_dim1;
		z__3.r = t2 * h__[i__5].r, z__3.i = t2 * h__[i__5].i;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		sum.r = z__1.r, sum.i = z__1.i;
		i__4 = j + k * h_dim1;
		i__5 = j + k * h_dim1;
		z__1.r = h__[i__5].r - sum.r, z__1.i = h__[i__5].i - sum.i;
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
		i__4 = j + (k + 1) * h_dim1;
		i__5 = j + (k + 1) * h_dim1;
		d_cnjg(&z__3, &v2);
		z__2.r = sum.r * z__3.r - sum.i * z__3.i, z__2.i = sum.r * 
			z__3.i + sum.i * z__3.r;
		z__1.r = h__[i__5].r - z__2.r, z__1.i = h__[i__5].i - z__2.i;
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
/* L90: */
	    }

	    if (*wantz) {

/*              Accumulate transformations in the matrix Z */

		i__3 = *ihiz;
		for (j = *iloz; j <= i__3; ++j) {
		    i__4 = j + k * z_dim1;
		    z__2.r = t1.r * z__[i__4].r - t1.i * z__[i__4].i, z__2.i =
			     t1.r * z__[i__4].i + t1.i * z__[i__4].r;
		    i__5 = j + (k + 1) * z_dim1;
		    z__3.r = t2 * z__[i__5].r, z__3.i = t2 * z__[i__5].i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    sum.r = z__1.r, sum.i = z__1.i;
		    i__4 = j + k * z_dim1;
		    i__5 = j + k * z_dim1;
		    z__1.r = z__[i__5].r - sum.r, z__1.i = z__[i__5].i - 
			    sum.i;
		    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
		    i__4 = j + (k + 1) * z_dim1;
		    i__5 = j + (k + 1) * z_dim1;
		    d_cnjg(&z__3, &v2);
		    z__2.r = sum.r * z__3.r - sum.i * z__3.i, z__2.i = sum.r *
			     z__3.i + sum.i * z__3.r;
		    z__1.r = z__[i__5].r - z__2.r, z__1.i = z__[i__5].i - 
			    z__2.i;
		    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
/* L100: */
		}
	    }

	    if (k == m && m > l) {

/*              If the QR step was started at row M > L because two */
/*              consecutive small subdiagonals were found, then extra */
/*              scaling must be performed to ensure that H(M,M-1) remains */
/*              real. */

		z__1.r = 1. - t1.r, z__1.i = 0. - t1.i;
		temp.r = z__1.r, temp.i = z__1.i;
		d__1 = z_abs(&temp);
		z__1.r = temp.r / d__1, z__1.i = temp.i / d__1;
		temp.r = z__1.r, temp.i = z__1.i;
		i__3 = m + 1 + m * h_dim1;
		i__4 = m + 1 + m * h_dim1;
		d_cnjg(&z__2, &temp);
		z__1.r = h__[i__4].r * z__2.r - h__[i__4].i * z__2.i, z__1.i =
			 h__[i__4].r * z__2.i + h__[i__4].i * z__2.r;
		h__[i__3].r = z__1.r, h__[i__3].i = z__1.i;
		if (m + 2 <= i__) {
		    i__3 = m + 2 + (m + 1) * h_dim1;
		    i__4 = m + 2 + (m + 1) * h_dim1;
		    z__1.r = h__[i__4].r * temp.r - h__[i__4].i * temp.i, 
			    z__1.i = h__[i__4].r * temp.i + h__[i__4].i * 
			    temp.r;
		    h__[i__3].r = z__1.r, h__[i__3].i = z__1.i;
		}
		i__3 = i__;
		for (j = m; j <= i__3; ++j) {
		    if (j != m + 1) {
			if (i2 > j) {
			    i__4 = i2 - j;
			    zscal_(&i__4, &temp, &h__[j + (j + 1) * h_dim1], 
				    ldh);
			}
			i__4 = j - i1;
			d_cnjg(&z__1, &temp);
			zscal_(&i__4, &z__1, &h__[i1 + j * h_dim1], &c__1);
			if (*wantz) {
			    d_cnjg(&z__1, &temp);
			    zscal_(&nz, &z__1, &z__[*iloz + j * z_dim1], &
				    c__1);
			}
		    }
/* L110: */
		}
	    }
/* L120: */
	}

/*        Ensure that H(I,I-1) is real. */

	i__2 = i__ + (i__ - 1) * h_dim1;
	temp.r = h__[i__2].r, temp.i = h__[i__2].i;
	if (d_imag(&temp) != 0.) {
	    rtemp = z_abs(&temp);
	    i__2 = i__ + (i__ - 1) * h_dim1;
	    h__[i__2].r = rtemp, h__[i__2].i = 0.;
	    z__1.r = temp.r / rtemp, z__1.i = temp.i / rtemp;
	    temp.r = z__1.r, temp.i = z__1.i;
	    if (i2 > i__) {
		i__2 = i2 - i__;
		d_cnjg(&z__1, &temp);
		zscal_(&i__2, &z__1, &h__[i__ + (i__ + 1) * h_dim1], ldh);
	    }
	    i__2 = i__ - i1;
	    zscal_(&i__2, &temp, &h__[i1 + i__ * h_dim1], &c__1);
	    if (*wantz) {
		zscal_(&nz, &temp, &z__[*iloz + i__ * z_dim1], &c__1);
	    }
	}

/* L130: */
    }

/*     Failure to converge in remaining number of iterations */

    *info = i__;
    return 0;

L140:

/*     H(I,I-1) is negligible: one eigenvalue has converged. */

    i__1 = i__;
    i__2 = i__ + i__ * h_dim1;
    w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;
/*     reset deflation counter */
    kdefl = 0;

/*     return to start of the main loop with new value of I. */

    i__ = l - 1;
    goto L30;

L150:
    return 0;

/*     End of ZLAHQR */

} /* zlahqr_ */

/* > \brief \b ZSYCONV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconv
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconv
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconv
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYCONV converts A given by ZHETRF into L and D or vice-versa. */
/* > Get nondiagonal elements of D (returned in workspace) and */
/* > apply or reverse permutation done in TRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is COMPLEX*16 array, dimension (N) */
/* >          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1 */
/* >          or 2-by-2 block diagonal matrix D in LDLT. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsyconv_(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *e, 
	integer *info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ip;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --e;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZSYCONV", &i__1, (ftnlen)7);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */

/*           Convert VALUE */

	    i__ = *n;
	    e[1].r = 0., e[1].i = 0.;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ - 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ - 1;
		    e[i__1].r = 0., e[i__1].i = 0.;
		    i__1 = i__ - 1 + i__ * a_dim1;
		    a[i__1].r = 0., a[i__1].i = 0.;
		    --i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0., e[i__1].i = 0.;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L12: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ - 1 + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ - 1 + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L13: */
			}
		    }
		    --i__;
		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */

/*           Revert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    ++i__;
		    if (i__ < *n) {
			i__1 = *n;
			for (j = i__ + 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ - 1 + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ - 1 + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		}
		++i__;
	    }

/*           Revert VALUE */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ - 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    --i__;
		}
		--i__;
	    }
	}

    } else {

/*        A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */

/*           Convert VALUE */

	    i__ = 1;
	    i__1 = *n;
	    e[i__1].r = 0., e[i__1].i = 0.;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ + 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ + 1;
		    e[i__1].r = 0., e[i__1].i = 0.;
		    i__1 = i__ + 1 + i__ * a_dim1;
		    a[i__1].r = 0., a[i__1].i = 0.;
		    ++i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0., e[i__1].i = 0.;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L22: */
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = ip + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = ip + j * a_dim1;
			    i__3 = i__ + 1 + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = i__ + 1 + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L23: */
			}
		    }
		    ++i__;
		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */

/*           Revert PERMUTATIONS */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {
		    ip = ipiv[i__];
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = i__ + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = i__ + j * a_dim1;
			    i__3 = ip + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = ip + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		} else {
		    ip = -ipiv[i__];
		    --i__;
		    if (i__ > 1) {
			i__1 = i__ - 1;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = i__ + 1 + j * a_dim1;
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
			    i__2 = i__ + 1 + j * a_dim1;
			    i__3 = ip + j * a_dim1;
			    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
			    i__2 = ip + j * a_dim1;
			    a[i__2].r = temp.r, a[i__2].i = temp.i;
			}
		    }
		}
		--i__;
	    }

/*           Revert VALUE */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ + 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    ++i__;
		}
		++i__;
	    }
	}
    }

    return 0;

/*     End of ZSYCONV */

} /* zsyconv_ */

/* > \brief \b ZSYCONVF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYCONVF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconv
f.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconv
f.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconv
f.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > ZSYCONVF converts the factorization output format used in */
/* > ZSYTRF provided on entry in parameter A into the factorization */
/* > output format used in ZSYTRF_RK (or ZSYTRF_BK) that is stored */
/* > on exit in parameters A and E. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in ZSYTRF into */
/* > the format used in ZSYTRF_RK (or ZSYTRF_BK). */
/* > */
/* > If parameter WAY = 'R': */
/* > ZSYCONVF performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in ZSYTRF_RK */
/* > (or ZSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in ZSYTRF that is stored */
/* > on exit in parameter A. It also converts in place details of */
/* > the intechanges stored in IPIV from the format used in ZSYTRF_RK */
/* > (or ZSYTRF_BK) into the format used in ZSYTRF. */
/* > */
/* > ZSYCONVF can also convert in Hermitian matrix case, i.e. between */
/* > formats used in ZHETRF and ZHETRF_RK (or ZHETRF_BK). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          ZSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          ZSYTRF_RK or ZSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          ZSYTRF_RK or ZSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          ZSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is COMPLEX*16 array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in ZSYTRF. */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in ZSYTRF_RK */
/* >          ( or ZSYTRF_BK). */
/* > */
/* >          1) If WAY ='R': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in ZSYTRF_RK */
/* >          ( or ZSYTRF_BK). */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16SYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int zsyconvf_(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	integer *info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ip;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZSYCONVF", &i__1, (ftnlen)8);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1].r = 0., e[1].i = 0.;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ - 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ - 1;
		    e[i__1].r = 0., e[i__1].i = 0.;
		    i__1 = i__ - 1 + i__ * a_dim1;
		    a[i__1].r = 0., a[i__1].i = 0.;
		    --i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0., e[i__1].i = 0.;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__ - 1) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i-1 and IPIV(i-1), */
/*                 so this should be recorded in two consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ - 1];

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ - 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    i__1 = *n;
	    e[i__1].r = 0., e[i__1].i = 0.;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ + 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ + 1;
		    e[i__1].r = 0., e[i__1].i = 0.;
		    i__1 = i__ + 1 + i__ * a_dim1;
		    a[i__1].r = 0., a[i__1].i = 0.;
		    ++i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0., e[i__1].i = 0.;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where k increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

		    ipiv[i__] = i__;

		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__ + 1) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i+1 and IPIV(i+1), */
/*                 so this should be recorded in consecutive entries */
/*                 in IPIV format for *SYTRF */

		    ipiv[i__] = ipiv[i__ + 1];

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ + 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of ZSYCONVF */

} /* zsyconvf_ */

/* > \brief \b ZSYCONVF_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYCONVF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconv
f_rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconv
f_rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconv
f_rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > ZSYCONVF_ROOK converts the factorization output format used in */
/* > ZSYTRF_ROOK provided on entry in parameter A into the factorization */
/* > output format used in ZSYTRF_RK (or ZSYTRF_BK) that is stored */
/* > on exit in parameters A and E. IPIV format for ZSYTRF_ROOK and */
/* > ZSYTRF_RK (or ZSYTRF_BK) is the same and is not converted. */
/* > */
/* > If parameter WAY = 'R': */
/* > ZSYCONVF_ROOK performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in ZSYTRF_RK */
/* > (or ZSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in ZSYTRF_ROOK that is stored */
/* > on exit in parameter A. IPIV format for ZSYTRF_ROOK and */
/* > ZSYTRF_RK (or ZSYTRF_BK) is the same and is not converted. */
/* > */
/* > ZSYCONVF_ROOK can also convert in Hermitian matrix case, i.e. between */
/* > formats used in ZHETRF_ROOK and ZHETRF_RK (or ZHETRF_BK). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          ZSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          ZSYTRF_RK or ZSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          ZSYTRF_RK or ZSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          ZSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is COMPLEX*16 array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D as determined: */
/* >          1) by ZSYTRF_ROOK, if WAY ='C'; */
/* >          2) by ZSYTRF_RK (or ZSYTRF_BK), if WAY ='R'. */
/* >          The IPIV format is the same for all these routines. */
/* > */
/* >          On exit, is not changed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16SYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int zsyconvf_rook__(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	integer *info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ip, ip2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZSYCONVF_ROOK", &i__1, (ftnlen)13);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Begin A is UPPER */

	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = *n;
	    e[1].r = 0., e[1].i = 0.;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ - 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ - 1;
		    e[i__1].r = 0., e[i__1].i = 0.;
		    i__1 = i__ - 1 + i__ * a_dim1;
		    a[i__1].r = 0., a[i__1].i = 0.;
		    --i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0., e[i__1].i = 0.;
		}
		--i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1) */
/*                 in A(1:i,N-i:N) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
			}
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip2 + (i__ + 1) * a_dim1], lda);
			}
		    }
		    --i__;

		}
		--i__;
	    }

	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

		    ip = ipiv[i__];
		    if (i__ < *n) {
			if (ip != i__) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i) */
/*                 in A(1:i,N-i:N) */

		    ++i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ - 1];
		    if (i__ < *n) {
			if (ip2 != i__ - 1) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[ip2 + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = *n - i__;
			    zswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
			}
		    }

		}
		++i__;
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

	    i__ = *n;
	    while(i__ > 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ - 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    --i__;
		}
		--i__;
	    }

/*        End A is UPPER */

	}

    } else {

/*        Begin A is LOWER */

	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

	    i__ = 1;
	    i__1 = *n;
	    e[i__1].r = 0., e[i__1].i = 0.;
	    while(i__ <= *n) {
		if (i__ < *n && ipiv[i__] < 0) {
		    i__1 = i__;
		    i__2 = i__ + 1 + i__ * a_dim1;
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
		    i__1 = i__ + 1;
		    e[i__1].r = 0., e[i__1].i = 0.;
		    i__1 = i__ + 1 + i__ * a_dim1;
		    a[i__1].r = 0., a[i__1].i = 0.;
		    ++i__;
		} else {
		    i__1 = i__;
		    e[i__1].r = 0., e[i__1].i = 0.;
		}
		++i__;
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in factorization order where i increases from 1 to N */

	    i__ = 1;
	    while(i__ <= *n) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1) */
/*                 in A(i:N,1:i-1) */

		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
			}
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip2 + 
				    a_dim1], lda);
			}
		    }
		    ++i__;

		}
		++i__;
	    }

	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS */

/*           Apply permutations to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

	    i__ = *n;
	    while(i__ >= 1) {
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

		    ip = ipiv[i__];
		    if (i__ > 1) {
			if (ip != i__) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i) */
/*                 in A(i:N,1:i-1) */

		    --i__;
		    ip = -ipiv[i__];
		    ip2 = -ipiv[i__ + 1];
		    if (i__ > 1) {
			if (ip2 != i__ + 1) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[ip2 + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
			}
			if (ip != i__) {
			    i__1 = i__ - 1;
			    zswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
			}
		    }

		}
		--i__;
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

	    i__ = 1;
	    while(i__ <= *n - 1) {
		if (ipiv[i__] < 0) {
		    i__1 = i__ + 1 + i__ * a_dim1;
		    i__2 = i__;
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
		    ++i__;
		}
		++i__;
	    }

	}

/*        End A is LOWER */

    }
    return 0;

/*     End of ZSYCONVF_ROOK */

} /* zsyconvf_rook__ */

/* > \brief \b ZTPMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpmqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpmqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpmqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/*      $          WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPMQRT applies a complex orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" complex block reflector H to a general */
/* > complex matrix C, which consists of two blocks A and B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Conjugate transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The order of the trapezoidal part of V. */
/* >          K >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CTPQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CTPQRT in B.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDV >= max(1,M); */
/* >          if SIDE = 'R', LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CTPQRT, stored as a NB-by-K matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension */
/* >          (LDA,N) if SIDE = 'L' or */
/* >          (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDC >= max(1,K); */
/* >          If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. */
/* >          LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The columns of the pentagonal matrix V contain the elementary reflectors */
/* >  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a */
/* >  trapezoidal block V2: */
/* > */
/* >        V = [V1] */
/* >            [V2]. */
/* > */
/* >  The size of the trapezoidal block V2 is determined by the parameter L, */
/* >  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L */
/* >  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular; */
/* >  if L=0, there is no trapezoidal block, hence V = V1 is rectangular. */
/* > */
/* >  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K. */
/* >                      [B] */
/* > */
/* >  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K. */
/* > */
/* >  The complex orthogonal matrix Q is formed from V and T. */
/* > */
/* >  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C. */
/* > */
/* >  If TRANS='C' and SIDE='L', C is on exit replaced with Q**H * C. */
/* > */
/* >  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* >  If TRANS='C' and SIDE='R', C is on exit replaced with C * Q**H. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztpmqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *l, integer *nb, doublecomplex *v, integer *ldv, 
	doublecomplex *t, integer *ldt, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *work, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, mb, kf, ldaq;
    static logical left, tran;
    static integer ldvq;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    extern /* Subroutine */ int ztprfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

    if (left) {
	ldvq = max(1,*m);
	ldaq = max(1,*k);
    } else if (right) {
	ldvq = max(1,*n);
	ldaq = max(1,*m);
    }
    if (! left && ! right) {
	*info = -1;
    } else if (! tran && ! notran) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0) {
	*info = -5;
    } else if (*l < 0 || *l > *k) {
	*info = -6;
    } else if (*nb < 1 || *nb > *k && *k > 0) {
	*info = -7;
    } else if (*ldv < ldvq) {
	*info = -9;
    } else if (*ldt < *nb) {
	*info = -11;
    } else if (*lda < ldaq) {
	*info = -13;
    } else if (*ldb < max(1,*m)) {
	*info = -15;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZTPMQRT", &i__1, (ftnlen)7);
	return 0;
    }

/*     .. Quick return if possible .. */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && tran) {

	i__1 = *k;
	i__2 = *nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *m - *l + i__ + ib - 1;
	    mb = min(i__3,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    ztprfb_("L", "C", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && notran) {

	i__2 = *k;
	i__1 = *nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *k - i__ + 1;
	    ib = min(i__3,i__4);
/* Computing MIN */
	    i__3 = *n - *l + i__ + ib - 1;
	    mb = min(i__3,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    ztprfb_("R", "N", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    } else if (left && notran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *m - *l + i__ + ib - 1;
	    mb = min(i__2,*m);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *m + *l - i__ + 1;
	    }
	    ztprfb_("L", "N", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib);
	}

    } else if (right && tran) {

	kf = (*k - 1) / *nb * *nb + 1;
	i__1 = -(*nb);
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
/* Computing MIN */
	    i__2 = *n - *l + i__ + ib - 1;
	    mb = min(i__2,*n);
	    if (i__ >= *l) {
		lb = 0;
	    } else {
		lb = mb - *n + *l - i__ + 1;
	    }
	    ztprfb_("R", "C", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m);
	}

    }

    return 0;

/*     End of ZTPMQRT */

} /* ztpmqrt_ */

/* > \brief \b ZTPQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTPQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpqrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpqrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpqrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPQRT computes a blocked QR factorization of a complex */
/* > "triangular-pentagonal" matrix C, which is composed of a */
/* > triangular block A and pentagonal block B, using the compact */
/* > WY representation for Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B, and the order of the */
/* >          triangular matrix A. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of rows of the upper trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size to be used in the blocked QR.  N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the upper triangular N-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          On entry, the pentagonal M-by-N matrix B.  The first M-L rows */
/* >          are rectangular, and the last L rows are upper trapezoidal. */
/* >          On exit, B contains the pentagonal matrix V.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,N) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (NB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a (N+M)-by-N matrix */
/* > */
/* >               C = [ A ] */
/* >                   [ B ] */
/* > */
/* >  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* >  upper trapezoidal matrix B2: */
/* > */
/* >               B = [ B1 ]  <- (M-L)-by-N rectangular */
/* >                   [ B2 ]  <-     L-by-N upper trapezoidal. */
/* > */
/* >  The upper trapezoidal matrix B2 consists of the first L rows of a */
/* >  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is upper triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* >               C = [ A ]  <- upper triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* > */
/* >               W = [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* > */
/* >               V = [ V1 ] <- (M-L)-by-N rectangular */
/* >                   [ V2 ] <-     L-by-N upper trapezoidal. */
/* > */
/* >  The columns of V represent the vectors which define the H(i)'s. */
/* > */
/* >  The number of blocks is B = ceiling(N/NB), where each */
/* >  block is of order NB except for the last block, which is of order */
/* >  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block */
/* >  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB */
/* >  for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* >               T = [T1 T2 ... TB]. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztpqrt_(integer *m, integer *n, integer *l, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *t, integer *ldt, doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, ib, lb, mb, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), ztprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), 
	    ztpqrt2_(integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *);


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*l < 0 || *l > min(*m,*n) && min(*m,*n) >= 0) {
	*info = -3;
    } else if (*nb < 1 || *nb > *n && *n > 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*m)) {
	*info = -8;
    } else if (*ldt < *nb) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZTPQRT", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = *n;
    i__2 = *nb;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*     Compute the QR factorization of the current block */

/* Computing MIN */
	i__3 = *n - i__ + 1;
	ib = min(i__3,*nb);
/* Computing MIN */
	i__3 = *m - *l + i__ + ib - 1;
	mb = min(i__3,*m);
	if (i__ >= *l) {
	    lb = 0;
	} else {
	    lb = mb - *m + *l - i__ + 1;
	}

	ztpqrt2_(&mb, &ib, &lb, &a[i__ + i__ * a_dim1], lda, &b[i__ * b_dim1 
		+ 1], ldb, &t[i__ * t_dim1 + 1], ldt, &iinfo);

/*     Update by applying H**H to B(:,I+IB:N) from the left */

	if (i__ + ib <= *n) {
	    i__3 = *n - i__ - ib + 1;
	    ztprfb_("L", "C", "F", "C", &mb, &i__3, &ib, &lb, &b[i__ * b_dim1 
		    + 1], ldb, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + ib) 
		    * a_dim1], lda, &b[(i__ + ib) * b_dim1 + 1], ldb, &work[1]
		    , &ib);
	}
    }
    return 0;

/*     End of ZTPQRT */

} /* ztpqrt_ */


typedef struct { float real, imag; } npy_complex64;
typedef struct { double real, imag; } npy_complex128;
int sorcsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs, int *m, int *p, int *q, float *x11, int *ldx11, float *x12, int *ldx12, float *x21, int *ldx21, float *x22, int *ldx22, float *theta, float *u1, int *ldu1, float *u2, int *ldu2, float *v1t, int *ldv1t, float *v2t, int *ldv2t, float *work, int *lwork, int *iwork, int *info){ return 0; }
int dorcsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs, int *m, int *p, int *q, double *x11, int *ldx11, double *x12, int *ldx12, double *x21, int *ldx21, double *x22, int *ldx22, double *theta, double *u1, int *ldu1, double *u2, int *ldu2, double *v1t, int *ldv1t, double *v2t, int *ldv2t, double *work, int *lwork, int *iwork, int *info){ return 0; }
int zuncsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs, int *m, int *p, int *q, npy_complex128 *x11, int *ldx11, npy_complex128 *x12, int *ldx12, npy_complex128 *x21, int *ldx21, npy_complex128 *x22, int *ldx22, double *theta, npy_complex128 *u1, int *ldu1, npy_complex128 *u2, int *ldu2, npy_complex128 *v1t, int *ldv1t, npy_complex128 *v2t, int *ldv2t, npy_complex128 *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *info){ return 0; }
int cuncsd_(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs, int *m, int *p, int *q, npy_complex64 *x11, int *ldx11, npy_complex64 *x12, int *ldx12, npy_complex64 *x21, int *ldx21, npy_complex64 *x22, int *ldx22, float *theta, npy_complex64 *u1, int *ldu1, npy_complex64 *u2, int *ldu2, npy_complex64 *v1t, int *ldv1t, npy_complex64 *v2t, int *ldv2t, npy_complex64 *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *info){ return 0; }
