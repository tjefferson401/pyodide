/* idz_svd.f -- translated by f2c (version 20160102).
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


/*       routine idzr_svd computes an approximation of specified rank */
/*       to a given matrix, in the usual SVD form U S V^*, */
/*       where U has orthonormal columns, V has orthonormal columns, */
/*       and S is diagonal. */

/*       routine idzp_svd computes an approximation of specified */
/*       precision to a given matrix, in the usual SVD form U S V^*, */
/*       where U has orthonormal columns, V has orthonormal columns, */
/*       and S is diagonal. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int idzr_svd__(integer *m, integer *n, doublecomplex *a, integer *krank, doublecomplex *u, doublecomplex *v, doublereal *s, integer *ier, doublecomplex *r__){
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j, k;
    extern /* Subroutine */ int idz_adjer__(integer *, integer *, doublecomplex *, doublecomplex *);    static integer ifadjoint, io;
    extern /* Subroutine */ int idzr_qrpiv__(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *);    static integer ldr, ldu;
    extern /* Subroutine */ int idz_qmatmat__(integer *, integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *);    static integer info;
    static char jobz[1];
    extern /* Subroutine */ int idz_permuter__(integer *, integer *, integer *, integer *, doublecomplex *);    static integer lwork;
    extern /* Subroutine */ int idz_retriever__(integer *, integer *, doublecomplex *, integer *, doublecomplex *);    static integer ldvadj;
    extern /* Subroutine */ int zgesdd_(char *, integer *, integer *, doublecomplex *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);

/*       constructs a rank-krank SVD  u diag(s) v^*  approximating a, */
/*       where u is an m x krank matrix whose columns are orthonormal, */
/*       v is an n x krank matrix whose columns are orthonormal, */
/*       and diag(s) is a diagonal krank x krank matrix whose entries */
/*       are all nonnegative. This routine combines a QR code */
/*       (which is based on plane/Householder reflections) */
/*       with the LAPACK routine zgesdd. */

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
/*            (krank+2)*n+8*min(m,n)+6*krank**2+8*krank */
/*            complex*16 elements long */

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

    idzr_qrpiv__(m, n, &a[a_offset], krank, &r__[1], &r__[io + 1]);


/*       Extract R from the QR decomposition. */

    idz_retriever__(m, n, &a[a_offset], krank, &r__[io + 1]);


/*       Rearrange R according to ind. */

    idz_permuter__(krank, &r__[1], krank, n, &r__[io + 1]);


/*       Use LAPACK to SVD r, */
/*       storing the krank (krank x 1) left singular vectors */
/*       in r(io+krank*n+1 : io+krank*n+krank*krank). */

    *(unsigned char *)jobz = 'S';
    ldr = *krank;
/* Computing 2nd power */
    i__1 = *krank;
    lwork = i__1 * i__1 + (*krank << 1) + *n << 1;
    ldu = *krank;
    ldvadj = *krank;

    zgesdd_(jobz, krank, n, &r__[io + 1], &ldr, &s[1], &r__[io + *krank * *n 
	    + 1], &ldu, &v[1], &ldvadj, &r__[io + *krank * *n + *krank * *
	    krank + 1], &lwork, &r__[io + *krank * *n + *krank * *krank + 
	    lwork + 1], &r__[1], &info, (ftnlen)1);

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
	    i__3 = j + k * u_dim1;
	    i__4 = io + *krank * *n + j + *krank * (k - 1);
	    u[i__3].r = r__[i__4].r, u[i__3].i = r__[i__4].i;
	}

/* j */
	i__2 = *m;
	for (j = *krank + 1; j <= i__2; ++j) {
	    i__3 = j + k * u_dim1;
	    u[i__3].r = 0., u[i__3].i = 0.;
	}

/* j */
    }

/* k */
    ifadjoint = 0;
    idz_qmatmat__(&ifadjoint, m, n, &a[a_offset], krank, krank, &u[u_offset], 
	    &r__[1]);


/*       Take the adjoint of v to obtain r. */

    idz_adjer__(krank, n, &v[1], &r__[1]);


/*       Copy r into v. */

    i__1 = *n * *krank;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	v[i__2].r = r__[i__3].r, v[i__2].i = r__[i__3].i;
    }


/* k */
    return 0;
} /* idzr_svd__ */





/* Subroutine */ int idzp_svd__(integer *lw, doublereal *eps, integer *m, integer *n, doublecomplex *a, integer *krank, integer *iu, integer * iv, integer *is, doublecomplex *w, integer *ier){
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j, k;
    extern /* Subroutine */ int idz_adjer__(integer *, integer *, doublecomplex *, doublecomplex *);    static integer ifadjoint, io, ls, lu, lv;
    extern /* Subroutine */ int idzp_qrpiv__(doublereal *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *);    static integer ldr, ldu, isi, ivi;
    extern /* Subroutine */ int idz_qmatmat__(integer *, integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *);    static integer info;
    static char jobz[1];
    extern /* Subroutine */ int idz_realcomp__(integer *, doublereal *, doublecomplex *), idz_permuter__(integer *, integer *, integer *, integer *, doublecomplex *);    static integer lwork;
    extern /* Subroutine */ int idz_retriever__(integer *, integer *, doublecomplex *, integer *, doublecomplex *);    static integer ldvadj;
    extern /* Subroutine */ int zgesdd_();


/*       constructs a rank-krank SVD  U Sigma V^*  approximating a */
/*       to precision eps, where U is an m x krank matrix whose */
/*       columns are orthonormal, V is an n x krank matrix whose */
/*       columns are orthonormal, and Sigma is a diagonal krank x krank */
/*       matrix whose entries are all nonnegative. */
/*       The entries of U are stored in w, starting at w(iu); */
/*       the entries of V are stored in w, starting at w(iv). */
/*       The diagonal entries of Sigma are stored in w, */
/*       starting at w(is). This routine combines a QR code */
/*       (which is based on plane/Householder reflections) */
/*       with the LAPACK routine zgesdd. */

/*       input: */
/*       lw -- maximum usable length of w (in complex*16 elements) */
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
/*             of singular values of a; the singular values are stored */
/*             as complex*16 numbers whose imaginary parts are zeros */
/*       w -- array containing the singular values and singular vectors */
/*            of a; w doubles as a work array, and so must be at least */
/*            (krank+1)*(m+2*n+9)+8*min(m,n)+6*krank**2 */
/*            complex*16 elements long, where krank is the rank */
/*            output by the present routine */
/*       ier -- 0 when the routine terminates successfully; */
/*              -1000 when lw is too small; */
/*              other nonzero values when zgesdd bombs */

/*       _N.B._: This routine destroys a. Also, please beware that */
/*               the source code for this routine could be clearer. */
/*               w must be at least */
/*               (krank+1)*(m+2*n+9)+8*min(m,n)+6*krank**2 */
/*               complex*16 elements long, where krank is the rank */
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

    idzp_qrpiv__(eps, m, n, &a[a_offset], krank, &w[1], &w[io + 1]);


    if (*krank > 0) {


/*         Extract R from the QR decomposition. */

	idz_retriever__(m, n, &a[a_offset], krank, &w[io + 1]);


/*         Rearrange R according to ind. */

	idz_permuter__(krank, &w[1], krank, n, &w[io + 1]);


/*         Use LAPACK to SVD R, */
/*         storing the krank (krank x 1) left singular vectors */
/*         in w(io+krank*n+1 : io+krank*n+krank*krank). */

	*(unsigned char *)jobz = 'S';
	ldr = *krank;
/* Computing 2nd power */
	i__1 = *krank;
	lwork = i__1 * i__1 + (*krank << 1) + *n << 1;
	ldu = *krank;
	ldvadj = *krank;

/* Computing 2nd power */
	i__1 = *krank;
	ivi = io + *krank * *n + *krank * *krank + lwork + i__1 * i__1 * 3 + (
		*krank << 2) + 1;
	lv = *n * *krank;

	isi = ivi + lv;
	ls = *krank;

	if (*lw < isi + ls + *m * *krank - 1) {
	    *ier = -1000;
	    return 0;
	}

	zgesdd_(jobz, krank, n, &w[io + 1], &ldr, &w[isi], &w[io + *krank * *
		n + 1], &ldu, &w[ivi], &ldvadj, &w[io + *krank * *n + *krank *
		 *krank + 1], &lwork, &w[io + *krank * *n + *krank * *krank + 
		lwork + 1], &w[1], &info, (ftnlen)1);

	if (info != 0) {
	    *ier = info;
	    return 0;
	}


/*         Take the adjoint of w(ivi:ivi+lv-1) to obtain V. */

	*iv = 1;
	idz_adjer__(krank, n, &w[ivi], &w[*iv]);


/*         Copy w(isi:isi+ls/2) into w(is:is+ls-1). */

	*is = *iv + lv;

	idz_realcomp__(&ls, &w[isi], &w[*is]);


/*         Multiply the U from R from the left by Q to obtain the U */
/*         for A. */

	*iu = *is + ls;
	lu = *m * *krank;

	i__1 = *krank;
	for (k = 1; k <= i__1; ++k) {

	    i__2 = *krank;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *iu - 1 + j + *krank * (k - 1);
		i__4 = io + *krank * *n + j + *krank * (k - 1);
		w[i__3].r = w[i__4].r, w[i__3].i = w[i__4].i;
	    }

/* j */
	}

/* k */
	for (k = *krank; k >= 1; --k) {

	    i__1 = *krank + 1;
	    for (j = *m; j >= i__1; --j) {
		i__2 = *iu - 1 + j + *m * (k - 1);
		w[i__2].r = 0., w[i__2].i = 0.;
	    }

/* j */
	    for (j = *krank; j >= 1; --j) {
		i__1 = *iu - 1 + j + *m * (k - 1);
		i__2 = *iu - 1 + j + *krank * (k - 1);
		w[i__1].r = w[i__2].r, w[i__1].i = w[i__2].i;
	    }

/* j */
	}

/* k */
	ifadjoint = 0;
	idz_qmatmat__(&ifadjoint, m, n, &a[a_offset], krank, krank, &w[*iu], &
		w[*iu + lu + 1]);


    }


/* krank .gt. 0 */
    return 0;
} /* idzp_svd__ */





/* Subroutine */ int idz_realcomp__(integer *n, doublereal *a, doublecomplex * b){
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k;


/*       copies the real*8 array a into the complex*16 array b. */

/*       input: */
/*       n -- length of a and b */
/*       a -- real*8 array to be copied into b */

/*       output: */
/*       b -- complex*16 copy of a */



    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	b[i__2].r = a[i__3], b[i__2].i = 0.;
    }


/* k */
    return 0;
} /* idz_realcomp__ */





/* Subroutine */ int idz_permuter__(integer *krank, integer *ind, integer *m, integer *n, doublecomplex *a){
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k;
    static doublecomplex cswap;


/*       permutes the columns of a according to ind obtained */
/*       from routine idzr_qrpiv or idzp_qrpiv, assuming that */
/*       a = q r from idzr_qrpiv or idzp_qrpiv. */

/*       input: */
/*       krank -- rank specified to routine idzr_qrpiv */
/*                or obtained from routine idzp_qrpiv */
/*       ind -- indexing array obtained from routine idzr_qrpiv */
/*              or idzp_qrpiv */
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

	    i__2 = j + k * a_dim1;
	    cswap.r = a[i__2].r, cswap.i = a[i__2].i;
	    i__2 = j + k * a_dim1;
	    i__3 = j + ind[k] * a_dim1;
	    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	    i__2 = j + ind[k] * a_dim1;
	    a[i__2].r = cswap.r, a[i__2].i = cswap.i;

	}
/* j */
    }


/* k */
    return 0;
} /* idz_permuter__ */





/* Subroutine */ int idz_retriever__(integer *m, integer *n, doublecomplex *a, integer *krank, doublecomplex *r__){
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j, k;


/*       extracts R in the QR decomposition specified by the output a */
/*       of the routine idzr_qrpiv or idzp_qrpiv */

/*       input: */
/*       m -- first dimension of a */
/*       n -- second dimension of a and r */
/*       a -- output of routine idzr_qrpiv or idzp_qrpiv */
/*       krank -- rank specified to routine idzr_qrpiv, */
/*                or output by routine idzp_qrpiv */

/*       output: */
/*       r -- triangular factor in the QR decomposition specified */
/*            by the output a of the routine idzr_qrpiv or idzp_qrpiv */



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
	    i__3 = j + k * r_dim1;
	    i__4 = j + k * a_dim1;
	    r__[i__3].r = a[i__4].r, r__[i__3].i = a[i__4].i;
	}
/* j */
    }

/* k */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (k < *krank) {
	    i__2 = *krank;
	    for (j = k + 1; j <= i__2; ++j) {
		i__3 = j + k * r_dim1;
		r__[i__3].r = 0., r__[i__3].i = 0.;
	    }
/* j */
	}
    }


/* k */
    return 0;
} /* idz_retriever__ */





/* Subroutine */ int idz_adjer__(integer *m, integer *n, doublecomplex *a, doublecomplex *aa){
    /* System generated locals */
    integer a_dim1, a_offset, aa_dim1, aa_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;


/*       forms the adjoint aa of a. */

/*       input: */
/*       m -- first dimension of a and second dimension of aa */
/*       n -- second dimension of a and first dimension of aa */
/*       a -- matrix whose adjoint is to be taken */

/*       output: */
/*       aa -- adjoint of a */



    /* Parameter adjustments */
    aa_dim1 = *n;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = k + j * aa_dim1;
	    d_cnjg(&z__1, &a[j + k * a_dim1]);
	    aa[i__3].r = z__1.r, aa[i__3].i = z__1.i;
	}
/* j */
    }


/* k */
    return 0;
} /* idz_adjer__ */

