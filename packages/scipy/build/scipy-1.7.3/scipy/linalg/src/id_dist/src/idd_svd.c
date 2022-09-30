/* idd_svd.f -- translated by f2c (version 20160102).
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

/*       this file contains the following user-callable routines: */


/*       routine iddr_svd computes an approximation of specified rank */
/*       to a given matrix, in the usual SVD form U S V^T, */
/*       where U has orthonormal columns, V has orthonormal columns, */
/*       and S is diagonal. */

/*       routine iddp_svd computes an approximation of specified */
/*       precision to a given matrix, in the usual SVD form U S V^T, */
/*       where U has orthonormal columns, V has orthonormal columns, */
/*       and S is diagonal. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int iddr_svd__(integer *m, integer *n, doublereal *a, integer *krank, doublereal *u, doublereal *v, doublereal *s, integer * ier, doublereal *r__){
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, i__1, i__2;

    /* Local variables */
    static integer j, k, io;
    extern /* Subroutine */ int iddr_qrpiv__(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *);    static integer ldr, ldu;
    extern /* Subroutine */ int idd_qmatmat__(integer *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *), idd_transer__(integer *, integer *, doublereal *, doublereal *);    static integer iftranspose, info;
    static char jobz[1];
    static integer ldvt;
    extern /* Subroutine */ int idd_permuter__(integer *, integer *, integer *, integer *, doublereal *);    static integer lwork;
    extern /* Subroutine */ int idd_retriever__(integer *, integer *, doublereal *, integer *, doublereal *), dgesdd_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, ftnlen);

/*       constructs a rank-krank SVD  u diag(s) v^T  approximating a, */
/*       where u is an m x krank matrix whose columns are orthonormal, */
/*       v is an n x krank matrix whose columns are orthonormal, */
/*       and diag(s) is a diagonal krank x krank matrix whose entries */
/*       are all nonnegative. This routine combines a QR code */
/*       (which is based on plane/Householder reflections) */
/*       with the LAPACK routine dgesdd. */

/*       input: */
/*       m -- first dimension of a and u */
/*       n -- second dimension of a, and first dimension of v */
/*       a -- matrix to be SVD'd */
/*       krank -- desired rank of the approximation to a */

/*       output: */
/*       u -- left singular vectors of a corresponding */
/*            to the k greatest singular values of a */
/*       v -- right singular vectors of a corresponding */
/*            to the k greatest singular values of a */
/*       s -- k greatest singular values of a */
/*       ier -- 0 when the routine terminates successfully; */
/*              nonzero when the routine encounters an error */

/*       work: */
/*       r -- must be at least */
/*            (krank+2)*n+8*min(m,n)+15*krank**2+8*krank */
/*            real*8 elements long */

/*       _N.B._: This routine destroys a. Also, please beware that */
/*               the source code for this routine could be clearer. */



    /* Parameter adjustments */
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    --v;
    u_dim1 = *m;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --r__;

    /* Function Body */
    io = min(*m,*n) << 3;


    *ier = 0;


/*       Compute a pivoted QR decomposition of a. */

    iddr_qrpiv__(m, n, &a[a_offset], krank, &r__[1], &r__[io + 1]);


/*       Extract R from the QR decomposition. */

    idd_retriever__(m, n, &a[a_offset], krank, &r__[io + 1]);


/*       Rearrange R according to ind (which is stored in r). */

    idd_permuter__(krank, &r__[1], krank, n, &r__[io + 1]);


/*       Use LAPACK to SVD R, */
/*       storing the krank (krank x 1) left singular vectors */
/*       in r(io+krank*n+1 : io+krank*n+krank*krank). */

    *(unsigned char *)jobz = 'S';
    ldr = *krank;
/* Computing 2nd power */
    i__1 = *krank;
/* Computing 2nd power */
    i__2 = *krank;
    lwork = i__1 * i__1 * 3 + *n + (i__2 * i__2 << 2) + (*krank << 2) << 1;
    ldu = *krank;
    ldvt = *krank;

    dgesdd_(jobz, krank, n, &r__[io + 1], &ldr, &s[1], &r__[io + *krank * *n 
	    + 1], &ldu, &v[1], &ldvt, &r__[io + *krank * *n + *krank * *krank 
	    + 1], &lwork, &r__[1], &info, (ftnlen)1);

    if (info != 0) {
	*ier = info;
	return 0;
    }


/*       Multiply the U from R from the left by Q to obtain the U */
/*       for A. */

    i__1 = *krank;
    for (k = 1; k <= i__1; ++k) {

	i__2 = *krank;
	for (j = 1; j <= i__2; ++j) {
	    u[j + k * u_dim1] = r__[io + *krank * *n + j + *krank * (k - 1)];
	}

/* j */
	i__2 = *m;
	for (j = *krank + 1; j <= i__2; ++j) {
	    u[j + k * u_dim1] = 0.;
	}

/* j */
    }

/* k */
    iftranspose = 0;
    idd_qmatmat__(&iftranspose, m, n, &a[a_offset], krank, krank, &u[u_offset]
	    , &r__[1]);


/*       Transpose v to obtain r. */

    idd_transer__(krank, n, &v[1], &r__[1]);


/*       Copy r into v. */

    i__1 = *n * *krank;
    for (k = 1; k <= i__1; ++k) {
	v[k] = r__[k];
    }


/* k */
    return 0;
} /* iddr_svd__ */





/* Subroutine */ int iddp_svd__(integer *lw, doublereal *eps, integer *m, integer *n, doublereal *a, integer *krank, integer *iu, integer *iv, integer *is, doublereal *w, integer *ier){
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j, k, io, ls, lu, lv;
    extern /* Subroutine */ int iddp_qrpiv__(doublereal *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *);    static integer ldr, ldu, isi, ivi;
    extern /* Subroutine */ int idd_qmatmat__(integer *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *), idd_transer__(integer *, integer *, doublereal *, doublereal *);    static integer iftranspose, info;
    static char jobz[1];
    static integer ldvt;
    extern /* Subroutine */ int idd_permuter__(integer *, integer *, integer *, integer *, doublereal *);    static integer lwork;
    extern /* Subroutine */ int idd_retriever__(integer *, integer *, doublereal *, integer *, doublereal *), dgesdd_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, ftnlen);

/*       constructs a rank-krank SVD  U Sigma V^T  approximating a */
/*       to precision eps, where U is an m x krank matrix whose */
/*       columns are orthonormal, V is an n x krank matrix whose */
/*       columns are orthonormal, and Sigma is a diagonal krank x krank */
/*       matrix whose entries are all nonnegative. */
/*       The entries of U are stored in w, starting at w(iu); */
/*       the entries of V are stored in w, starting at w(iv). */
/*       The diagonal entries of Sigma are stored in w, */
/*       starting at w(is). This routine combines a QR code */
/*       (which is based on plane/Householder reflections) */
/*       with the LAPACK routine dgesdd. */

/*       input: */
/*       lw -- maximum usable length of w (in real*8 elements) */
/*       eps -- precision to which the SVD approximates a */
/*       m -- first dimension of a and u */
/*       n -- second dimension of a, and first dimension of v */
/*       a -- matrix to be SVD'd */

/*       output: */
/*       krank -- rank of the approximation to a */
/*       iu -- index in w of the first entry of the matrix */
/*             of orthonormal left singular vectors of a */
/*       iv -- index in w of the first entry of the matrix */
/*             of orthonormal right singular vectors of a */
/*       is -- index in w of the first entry of the array */
/*             of singular values of a */
/*       w -- array containing the singular values and singular vectors */
/*            of a; w doubles as a work array, and so must be at least */
/*            (krank+1)*(m+2*n+9)+8*min(m,n)+15*krank**2 */
/*            real*8 elements long, where krank is the rank */
/*            output by the present routine */
/*       ier -- 0 when the routine terminates successfully; */
/*              -1000 when lw is too small; */
/*              other nonzero values when dgesdd bombs */

/*       _N.B._: This routine destroys a. Also, please beware that */
/*               the source code for this routine could be clearer. */
/*               w must be at least */
/*               (krank+1)*(m+2*n+9)+8*min(m,n)+15*krank**2 */
/*               real*8 elements long, where krank is the rank */
/*               output by the present routine. */



    /* Parameter adjustments */
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;

    /* Function Body */
    io = min(*m,*n) << 3;


    *ier = 0;


/*       Compute a pivoted QR decomposition of a. */

    iddp_qrpiv__(eps, m, n, &a[a_offset], krank, &w[1], &w[io + 1]);


    if (*krank > 0) {


/*         Extract R from the QR decomposition. */

	idd_retriever__(m, n, &a[a_offset], krank, &w[io + 1]);


/*         Rearrange R according to ind (which is stored in w). */

	idd_permuter__(krank, &w[1], krank, n, &w[io + 1]);


/*         Use LAPACK to SVD R, */
/*         storing the krank (krank x 1) left singular vectors */
/*         in w(io+krank*n+1 : io+krank*n+krank*krank). */

	*(unsigned char *)jobz = 'S';
	ldr = *krank;
/* Computing 2nd power */
	i__1 = *krank;
/* Computing 2nd power */
	i__2 = *krank;
	lwork = i__1 * i__1 * 3 + *n + (i__2 * i__2 << 2) + (*krank << 2) << 
		1;
	ldu = *krank;
	ldvt = *krank;

	ivi = io + *krank * *n + *krank * *krank + lwork + 1;
	lv = *n * *krank;

	isi = ivi + lv;
	ls = *krank;

	if (*lw < isi + ls + *m * *krank - 1) {
	    *ier = -1000;
	    return 0;
	}

	dgesdd_(jobz, krank, n, &w[io + 1], &ldr, &w[isi], &w[io + *krank * *
		n + 1], &ldu, &w[ivi], &ldvt, &w[io + *krank * *n + *krank * *
		krank + 1], &lwork, &w[1], &info, (ftnlen)1);

	if (info != 0) {
	    *ier = info;
	    return 0;
	}


/*         Transpose w(ivi:ivi+lv-1) to obtain V. */

	*iv = 1;
	idd_transer__(krank, n, &w[ivi], &w[*iv]);


/*         Copy w(isi:isi+ls-1) into w(is:is+ls-1). */

	*is = *iv + lv;

	i__1 = ls;
	for (k = 1; k <= i__1; ++k) {
	    w[*is + k - 1] = w[isi + k - 1];
	}


/*         Multiply the U from R from the left by Q to obtain the U */
/*         for A. */

/* k */
	*iu = *is + ls;
	lu = *m * *krank;

	i__1 = *krank;
	for (k = 1; k <= i__1; ++k) {

	    i__2 = *krank;
	    for (j = 1; j <= i__2; ++j) {
		w[*iu - 1 + j + *krank * (k - 1)] = w[io + *krank * *n + j + *
			krank * (k - 1)];
	    }

/* j */
	}

/* k */
	for (k = *krank; k >= 1; --k) {

	    i__1 = *krank + 1;
	    for (j = *m; j >= i__1; --j) {
		w[*iu - 1 + j + *m * (k - 1)] = 0.;
	    }

/* j */
	    for (j = *krank; j >= 1; --j) {
		w[*iu - 1 + j + *m * (k - 1)] = w[*iu - 1 + j + *krank * (k - 
			1)];
	    }

/* j */
	}

/* k */
	iftranspose = 0;
	idd_qmatmat__(&iftranspose, m, n, &a[a_offset], krank, krank, &w[*iu],
		 &w[*iu + lu + 1]);


    }


/* krank .gt. 0 */
    return 0;
} /* iddp_svd__ */





/* Subroutine */ int idd_permuter__(integer *krank, integer *ind, integer *m, integer *n, doublereal *a){
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer j, k;
    static doublereal rswap;


/*       permutes the columns of a according to ind obtained */
/*       from routine iddr_qrpiv or iddp_qrpiv, assuming that */
/*       a = q r from iddr_qrpiv or iddp_qrpiv. */

/*       input: */
/*       krank -- rank specified to routine iddr_qrpiv */
/*                or obtained from routine iddp_qrpiv */
/*       ind -- indexing array obtained from routine iddr_qrpiv */
/*              or iddp_qrpiv */
/*       m -- first dimension of a */
/*       n -- second dimension of a */
/*       a -- matrix to be rearranged */

/*       output: */
/*       a -- rearranged matrix */



    /* Parameter adjustments */
    --ind;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    for (k = *krank; k >= 1; --k) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {

	    rswap = a[j + k * a_dim1];
	    a[j + k * a_dim1] = a[j + ind[k] * a_dim1];
	    a[j + ind[k] * a_dim1] = rswap;

	}
/* j */
    }


/* k */
    return 0;
} /* idd_permuter__ */





/* Subroutine */ int idd_retriever__(integer *m, integer *n, doublereal *a, integer *krank, doublereal *r__){
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    static integer j, k;


/*       extracts R in the QR decomposition specified by the output a */
/*       of the routine iddr_qrpiv or iddp_qrpiv */

/*       input: */
/*       m -- first dimension of a */
/*       n -- second dimension of a and r */
/*       a -- output of routine iddr_qrpiv or iddp_qrpiv */
/*       krank -- rank specified to routine iddr_qrpiv, */
/*                or output by routine iddp_qrpiv */

/*       output: */
/*       r -- triangular factor in the QR decomposition specified */
/*            by the output a of the routine iddr_qrpiv or iddp_qrpiv */



/*       Copy a into r and zero out the appropriate */
/*       Householder vectors that are stored in one triangle of a. */

    /* Parameter adjustments */
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    r_dim1 = *krank;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *krank;
	for (j = 1; j <= i__2; ++j) {
	    r__[j + k * r_dim1] = a[j + k * a_dim1];
	}
/* j */
    }

/* k */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (k < *krank) {
	    i__2 = *krank;
	    for (j = k + 1; j <= i__2; ++j) {
		r__[j + k * r_dim1] = 0.;
	    }
/* j */
	}
    }


/* k */
    return 0;
} /* idd_retriever__ */





/* Subroutine */ int idd_transer__(integer *m, integer *n, doublereal *a, doublereal *at){
    /* System generated locals */
    integer a_dim1, a_offset, at_dim1, at_offset, i__1, i__2;

    /* Local variables */
    static integer j, k;


/*       forms the transpose at of a. */

/*       input: */
/*       m -- first dimension of a and second dimension of at */
/*       n -- second dimension of a and first dimension of at */
/*       a -- matrix to be transposed */

/*       output: */
/*       at -- transpose of a */



    /* Parameter adjustments */
    at_dim1 = *n;
    at_offset = 1 + at_dim1;
    at -= at_offset;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    at[k + j * at_dim1] = a[j + k * a_dim1];
	}
/* j */
    }


/* k */
    return 0;
} /* idd_transer__ */

