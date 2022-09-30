/* idd_id2svd.f -- translated by f2c (version 20160102).
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


/*       routine idd_id2svd converts an approximation to a matrix */
/*       in the form of an ID to an approximation in the form of an SVD. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int idd_id2svd__(integer *m, integer *krank, doublereal *b, integer *n, integer *list, doublereal *proj, doublereal *u, doublereal *v, doublereal *s, integer *ier, doublereal *w){
    /* System generated locals */
    integer b_dim1, b_offset, proj_dim1, proj_offset, u_dim1, u_offset, 
	    v_dim1, v_offset, i__1;

    /* Local variables */
    static integer ip, ir, lp, it, lr, lt, lw;
    extern /* Subroutine */ int idd_id2svd0__(integer *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, integer *);    static integer ir2, ir3, lr2, lr3, iind, lind, iindt, lindt, iwork, lwork;


/*       converts an approximation to a matrix in the form of an ID */
/*       to an approximation in the form of an SVD. */

/*       input: */
/*       m -- first dimension of b */
/*       krank -- rank of the ID */
/*       b -- columns of the original matrix in the ID */
/*       list -- list of columns chosen from the original matrix */
/*               in the ID */
/*       n -- length of list and part of the second dimension of proj */
/*       proj -- projection coefficients in the ID */

/*       output: */
/*       u -- left singular vectors */
/*       v -- right singular vectors */
/*       s -- singular values */
/*       ier -- 0 when the routine terminates successfully; */
/*              nonzero otherwise */

/*       work: */
/*       w -- must be at least (krank+1)*(m+3*n)+26*krank**2 real*8 */
/*            elements long */

/*       _N.B._: This routine destroys b. */



/* Computing 2nd power */
    i__1 = *krank;
    /* Parameter adjustments */
    --s;
    u_dim1 = *m;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    b_dim1 = *m;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --w;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    proj_dim1 = *krank;
    proj_offset = 1 + proj_dim1;
    proj -= proj_offset;
    --list;

    /* Function Body */
    lw = 0;

    iwork = lw + 1;
/* Computing 2nd power */
    i__1 = *krank;
    lwork = i__1 * i__1 * 25;
    lw += lwork;

    ip = lw + 1;
    lp = *krank * *n;
    lw += lp;

    it = lw + 1;
    lt = *n * *krank;
    lw += lt;

    ir = lw + 1;
    lr = *krank * *n;
    lw += lr;

    ir2 = lw + 1;
    lr2 = *krank * *m;
    lw += lr2;

    ir3 = lw + 1;
    lr3 = *krank * *krank;
    lw += lr3;

    iind = lw + 1;
    lind = *n / 2 + 1;
    ++lw;

    iindt = lw + 1;
    lindt = *m / 2 + 1;
    ++lw;


    idd_id2svd0__(m, krank, &b[b_offset], n, &list[1], &proj[proj_offset], &u[
	    u_offset], &v[v_offset], &s[1], ier, &w[iwork], &w[ip], &w[it], &
	    w[ir], &w[ir2], &w[ir3], &w[iind], &w[iindt]);


    return 0;
} /* idd_id2svd__ */





/* Subroutine */ int idd_id2svd0__(integer *m, integer *krank, doublereal *b, integer *n, integer *list, doublereal *proj, doublereal *u, doublereal *v, doublereal *s, integer *ier, doublereal *work, doublereal *p, doublereal *t, doublereal *r__, doublereal *r2, doublereal *r3, integer *ind, integer *indt){
    /* System generated locals */
    integer b_dim1, b_offset, proj_dim1, proj_offset, p_dim1, p_offset, 
	    r_dim1, r_offset, r2_dim1, r2_offset, t_dim1, t_offset, r3_dim1, 
	    r3_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    static integer j, k;
    extern /* Subroutine */ int idd_rearr__(integer *, integer *, integer *, integer *, doublereal *), idd_rinqr__(integer *, integer *, doublereal *, integer *, doublereal *), iddr_qrpiv__(integer *, integer *, doublereal *, integer *, integer *, doublereal *);    static integer ldr, ldu;
    extern /* Subroutine */ int idd_qmatmat__(integer *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *);    static integer iftranspose, info;
    static char jobz[1];
    static integer ldvt;
    extern /* Subroutine */ int idd_reconint__(integer *, integer *, integer * , doublereal *, doublereal *), idd_mattrans__(integer *, integer * , doublereal *, doublereal *), idd_matmultt__(integer *, integer * , doublereal *, integer *, doublereal *, doublereal *);    static integer lwork;
    extern /* Subroutine */ int dgesdd_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, ftnlen);

/*       routine idd_id2svd serves as a memory wrapper */
/*       for the present routine (please see routine idd_id2svd */
/*       for further documentation). */





/* Computing 2nd power */
    i__1 = *krank;
    /* Parameter adjustments */
    --indt;
    r3_dim1 = *krank;
    r3_offset = 1 + r3_dim1;
    r3 -= r3_offset;
    r2_dim1 = *krank;
    r2_offset = 1 + r2_dim1;
    r2 -= r2_offset;
    --s;
    u_dim1 = *m;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    b_dim1 = *m;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --ind;
    r_dim1 = *krank;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    t_dim1 = *n;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    p_dim1 = *krank;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    proj_dim1 = *krank;
    proj_offset = 1 + proj_dim1;
    proj -= proj_offset;
    --list;
    --work;

    /* Function Body */
    *ier = 0;



/*       Construct the projection matrix p from the ID. */

    idd_reconint__(n, &list[1], krank, &proj[proj_offset], &p[p_offset]);



/*       Compute a pivoted QR decomposition of b. */

    iddr_qrpiv__(m, krank, &b[b_offset], krank, &ind[1], &r__[r_offset]);


/*       Extract r from the QR decomposition. */

    idd_rinqr__(m, krank, &b[b_offset], krank, &r__[r_offset]);


/*       Rearrange r according to ind. */

    idd_rearr__(krank, &ind[1], krank, krank, &r__[r_offset]);



/*       Transpose p to obtain t. */

    idd_mattrans__(krank, n, &p[p_offset], &t[t_offset]);


/*       Compute a pivoted QR decomposition of t. */

    iddr_qrpiv__(n, krank, &t[t_offset], krank, &indt[1], &r2[r2_offset]);


/*       Extract r2 from the QR decomposition. */

    idd_rinqr__(n, krank, &t[t_offset], krank, &r2[r2_offset]);


/*       Rearrange r2 according to indt. */

    idd_rearr__(krank, &indt[1], krank, krank, &r2[r2_offset]);



/*       Multiply r and r2^T to obtain r3. */

    idd_matmultt__(krank, krank, &r__[r_offset], krank, &r2[r2_offset], &r3[
	    r3_offset]);



/*       Use LAPACK to SVD r3. */

    *(unsigned char *)jobz = 'S';
    ldr = *krank;
/* Computing 2nd power */
    i__1 = *krank;
/* Computing 2nd power */
    i__2 = *krank;
    lwork = i__1 * i__1 * 25 - i__2 * i__2 - (*krank << 2);
    ldu = *krank;
    ldvt = *krank;

/* Computing 2nd power */
    i__1 = *krank;
/* Computing 2nd power */
    i__2 = *krank;
    dgesdd_(jobz, krank, krank, &r3[r3_offset], &ldr, &s[1], &work[1], &ldu, &
	    r__[r_offset], &ldvt, &work[i__1 * i__1 + (*krank << 2) + 1], &
	    lwork, &work[i__2 * i__2 + 1], &info, (ftnlen)1);

    if (info != 0) {
	*ier = info;
	return 0;
    }



/*       Multiply the u from r3 from the left by the q from b */
/*       to obtain the u for a. */

    i__1 = *krank;
    for (k = 1; k <= i__1; ++k) {

	i__2 = *krank;
	for (j = 1; j <= i__2; ++j) {
	    u[j + k * u_dim1] = work[j + *krank * (k - 1)];
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
    idd_qmatmat__(&iftranspose, m, krank, &b[b_offset], krank, krank, &u[
	    u_offset], &r2[r2_offset]);



/*       Transpose r to obtain r2. */

    idd_mattrans__(krank, krank, &r__[r_offset], &r2[r2_offset]);


/*       Multiply the v from r3 from the left by the q from p^T */
/*       to obtain the v for a. */

    i__1 = *krank;
    for (k = 1; k <= i__1; ++k) {

	i__2 = *krank;
	for (j = 1; j <= i__2; ++j) {
	    v[j + k * v_dim1] = r2[j + k * r2_dim1];
	}

/* j */
	i__2 = *n;
	for (j = *krank + 1; j <= i__2; ++j) {
	    v[j + k * v_dim1] = 0.;
	}

/* j */
    }

/* k */
    iftranspose = 0;
    idd_qmatmat__(&iftranspose, n, krank, &t[t_offset], krank, krank, &v[
	    v_offset], &r2[r2_offset]);


    return 0;
} /* idd_id2svd0__ */





/* Subroutine */ int idd_mattrans__(integer *m, integer *n, doublereal *a, doublereal *at){
    /* System generated locals */
    integer a_dim1, a_offset, at_dim1, at_offset, i__1, i__2;

    /* Local variables */
    static integer j, k;


/*       transposes a to obtain at. */

/*       input: */
/*       m -- first dimension of a, and second dimension of at */
/*       n -- second dimension of a, and first dimension of at */
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
} /* idd_mattrans__ */





/* Subroutine */ int idd_matmultt__(integer *l, integer *m, doublereal *a, integer *n, doublereal *b, doublereal *c__){
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum;


/*       multiplies a and b^T to obtain c. */

/*       input: */
/*       l -- first dimension of a and c */
/*       m -- second dimension of a and b */
/*       a -- leftmost matrix in the product c = a b^T */
/*       n -- first dimension of b and second dimension of c */
/*       b -- rightmost matrix in the product c = a b^T */

/*       output: */
/*       c -- product of a and b^T */



    /* Parameter adjustments */
    a_dim1 = *l;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *l;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {

	    sum = 0.;

	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
		sum += a[i__ + j * a_dim1] * b[k + j * b_dim1];
	    }

/* j */
	    c__[i__ + k * c_dim1] = sum;

	}
/* k */
    }


/* i */
    return 0;
} /* idd_matmultt__ */





/* Subroutine */ int idd_rearr__(integer *krank, integer *ind, integer *m, integer *n, doublereal *a){
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer j, k;
    static doublereal rswap;


/*       rearranges a according to ind obtained */
/*       from routines iddr_qrpiv or iddp_qrpiv, */
/*       assuming that a = q r, where q and r are from iddr_qrpiv */
/*       or iddp_qrpiv. */

/*       input: */
/*       krank -- rank obtained from routine iddp_qrpiv, */
/*                or provided to routine iddr_qrpiv */
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
} /* idd_rearr__ */





/* Subroutine */ int idd_rinqr__(integer *m, integer *n, doublereal *a, integer *krank, doublereal *r__){
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    static integer j, k;


/*       extracts R in the QR decomposition specified by the output a */
/*       of the routine iddr_qrpiv or iddp_qrpiv. */

/*       input: */
/*       m -- first dimension of a */
/*       n -- second dimension of a and r */
/*       a -- output of routine iddr_qrpiv or iddp_qrpiv */
/*       krank -- rank output by routine iddp_qrpiv (or specified */
/*                to routine iddr_qrpiv) */

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
} /* idd_rinqr__ */

