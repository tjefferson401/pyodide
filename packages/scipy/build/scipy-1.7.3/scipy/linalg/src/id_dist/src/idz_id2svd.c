/* idz_id2svd.f -- translated by f2c (version 20160102).
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


/*       routine idz_id2svd converts an approximation to a matrix */
/*       in the form of an ID to an approximation in the form of an SVD. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int idz_id2svd__(integer *m, integer *krank, doublecomplex * b, integer *n, integer *list, doublecomplex *proj, doublecomplex *u, doublecomplex *v, doublereal *s, integer *ier, doublecomplex *w){
    /* System generated locals */
    integer b_dim1, b_offset, proj_dim1, proj_offset, u_dim1, u_offset, 
	    v_dim1, v_offset, i__1;

    /* Local variables */
    static integer ip, ir, lp, it, lr, lt, lw;
    extern /* Subroutine */ int idz_id2svd0__(integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, doublereal *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, integer *);    static integer ir2, ir3, lr2, lr3, iind, lind, iindt, lindt, iwork, lwork;


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
/*       w -- must be at least (krank+1)*(m+3*n+10)+9*krank**2 */
/*            complex*16 elements long */

/*       _N.B._: This routine destroys b. */



/*       Allocate memory for idz_id2svd0. */

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
    lwork = (i__1 * i__1 << 3) + *krank * 10;
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
    lind = *n / 4 + 1;
    ++lw;

    iindt = lw + 1;
    lindt = *m / 4 + 1;
    ++lw;


    idz_id2svd0__(m, krank, &b[b_offset], n, &list[1], &proj[proj_offset], &u[
	    u_offset], &v[v_offset], &s[1], ier, &w[iwork], &w[ip], &w[it], &
	    w[ir], &w[ir2], &w[ir3], &w[iind], &w[iindt]);


    return 0;
} /* idz_id2svd__ */





/* Subroutine */ int idz_id2svd0__(integer *m, integer *krank, doublecomplex * b, integer *n, integer *list, doublecomplex *proj, doublecomplex *u, doublecomplex *v, doublereal *s, integer *ier, doublecomplex *work, doublecomplex *p, doublecomplex *t, doublecomplex *r__, doublecomplex *r2, doublecomplex *r3, integer *ind, integer *indt){
    /* System generated locals */
    integer b_dim1, b_offset, proj_dim1, proj_offset, p_dim1, p_offset, 
	    r_dim1, r_offset, r2_dim1, r2_offset, t_dim1, t_offset, r3_dim1, 
	    r3_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3, 
	    i__4;

    /* Local variables */
    static integer j, k, ifadjoint;
    extern /* Subroutine */ int idz_rearr__(integer *, integer *, integer *, integer *, doublecomplex *), idz_rinqr__(integer *, integer *, doublecomplex *, integer *, doublecomplex *), idz_matadj__( integer *, integer *, doublecomplex *, doublecomplex *), idzr_qrpiv__(integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *);    static integer ldr, ldu;
    extern /* Subroutine */ int idz_qmatmat__(integer *, integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, doublecomplex *);    static integer info;
    static char jobz[1];
    static integer ldvt;
    extern /* Subroutine */ int idz_reconint__(integer *, integer *, integer * , doublecomplex *, doublecomplex *), idz_matmulta__(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *);    static integer lwork;
    extern /* Subroutine */ int zgesdd_(char *, integer *, integer *, doublecomplex *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);

/*       routine idz_id2svd serves as a memory wrapper */
/*       for the present routine (please see routine idz_id2svd */
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
    --work;
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

    /* Function Body */
    *ier = 0;



/*       Construct the projection matrix p from the ID. */

    idz_reconint__(n, &list[1], krank, &proj[proj_offset], &p[p_offset]);



/*       Compute a pivoted QR decomposition of b. */

    idzr_qrpiv__(m, krank, &b[b_offset], krank, &ind[1], &r__[r_offset]);


/*       Extract r from the QR decomposition. */

    idz_rinqr__(m, krank, &b[b_offset], krank, &r__[r_offset]);


/*       Rearrange r according to ind. */

    idz_rearr__(krank, &ind[1], krank, krank, &r__[r_offset]);



/*       Take the adjoint of p to obtain t. */

    idz_matadj__(krank, n, &p[p_offset], &t[t_offset]);


/*       Compute a pivoted QR decomposition of t. */

    idzr_qrpiv__(n, krank, &t[t_offset], krank, &indt[1], &r2[r2_offset]);


/*       Extract r2 from the QR decomposition. */

    idz_rinqr__(n, krank, &t[t_offset], krank, &r2[r2_offset]);


/*       Rearrange r2 according to indt. */

    idz_rearr__(krank, &indt[1], krank, krank, &r2[r2_offset]);



/*       Multiply r and r2^* to obtain r3. */

    idz_matmulta__(krank, krank, &r__[r_offset], krank, &r2[r2_offset], &r3[
	    r3_offset]);



/*       Use LAPACK to SVD r3. */

    *(unsigned char *)jobz = 'S';
    ldr = *krank;
/* Computing 2nd power */
    i__1 = *krank;
/* Computing 2nd power */
    i__2 = *krank;
/* Computing 2nd power */
    i__3 = *krank;
    lwork = (i__1 * i__1 << 3) + *krank * 10 - (i__2 * i__2 + (*krank << 1) + 
	    i__3 * i__3 * 3 + (*krank << 2));
    ldu = *krank;
    ldvt = *krank;

/* Computing 2nd power */
    i__1 = *krank;
/* Computing 2nd power */
    i__2 = *krank;
/* Computing 2nd power */
    i__3 = *krank;
/* Computing 2nd power */
    i__4 = *krank;
    zgesdd_(jobz, krank, krank, &r3[r3_offset], &ldr, &s[1], &work[1], &ldu, &
	    r__[r_offset], &ldvt, &work[i__1 * i__1 + (*krank << 1) + i__2 * 
	    i__2 * 3 + (*krank << 2) + 1], &lwork, &work[i__3 * i__3 + (*
	    krank << 1) + 1], &work[i__4 * i__4 + 1], &info, (ftnlen)1);

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
	    i__3 = j + k * u_dim1;
	    i__4 = j + *krank * (k - 1);
	    u[i__3].r = work[i__4].r, u[i__3].i = work[i__4].i;
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
    idz_qmatmat__(&ifadjoint, m, krank, &b[b_offset], krank, krank, &u[
	    u_offset], &r2[r2_offset]);



/*       Take the adjoint of r to obtain r2. */

    idz_matadj__(krank, krank, &r__[r_offset], &r2[r2_offset]);


/*       Multiply the v from r3 from the left by the q from p^* */
/*       to obtain the v for a. */

    i__1 = *krank;
    for (k = 1; k <= i__1; ++k) {

	i__2 = *krank;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = j + k * v_dim1;
	    i__4 = j + k * r2_dim1;
	    v[i__3].r = r2[i__4].r, v[i__3].i = r2[i__4].i;
	}

/* j */
	i__2 = *n;
	for (j = *krank + 1; j <= i__2; ++j) {
	    i__3 = j + k * v_dim1;
	    v[i__3].r = 0., v[i__3].i = 0.;
	}

/* j */
    }

/* k */
    ifadjoint = 0;
    idz_qmatmat__(&ifadjoint, n, krank, &t[t_offset], krank, krank, &v[
	    v_offset], &r2[r2_offset]);


    return 0;
} /* idz_id2svd0__ */





/* Subroutine */ int idz_matadj__(integer *m, integer *n, doublecomplex *a, doublecomplex *aa){
    /* System generated locals */
    integer a_dim1, a_offset, aa_dim1, aa_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;


/*       Takes the adjoint of a to obtain aa. */

/*       input: */
/*       m -- first dimension of a, and second dimension of aa */
/*       n -- second dimension of a, and first dimension of aa */
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
} /* idz_matadj__ */





/* Subroutine */ int idz_matmulta__(integer *l, integer *m, doublecomplex *a, integer *n, doublecomplex *b, doublecomplex *c__){
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex sum;


/*       multiplies a and b^* to obtain c. */

/*       input: */
/*       l -- first dimension of a and c */
/*       m -- second dimension of a and b */
/*       a -- leftmost matrix in the product c = a b^* */
/*       n -- first dimension of b and second dimension of c */
/*       b -- rightmost matrix in the product c = a b^* */

/*       output: */
/*       c -- product of a and b^* */



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

	    sum.r = 0., sum.i = 0.;

	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
		i__4 = i__ + j * a_dim1;
		d_cnjg(&z__3, &b[k + j * b_dim1]);
		z__2.r = a[i__4].r * z__3.r - a[i__4].i * z__3.i, z__2.i = a[
			i__4].r * z__3.i + a[i__4].i * z__3.r;
		z__1.r = sum.r + z__2.r, z__1.i = sum.i + z__2.i;
		sum.r = z__1.r, sum.i = z__1.i;
	    }

/* j */
	    i__3 = i__ + k * c_dim1;
	    c__[i__3].r = sum.r, c__[i__3].i = sum.i;

	}
/* k */
    }


/* i */
    return 0;
} /* idz_matmulta__ */





/* Subroutine */ int idz_rearr__(integer *krank, integer *ind, integer *m, integer *n, doublecomplex *a){
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k;
    static doublecomplex cswap;


/*       rearranges a according to ind obtained */
/*       from routines idzr_qrpiv or idzp_qrpiv, */
/*       assuming that a = q r, where q and r are from idzr_qrpiv */
/*       or idzp_qrpiv. */

/*       input: */
/*       krank -- rank obtained from routine idzp_qrpiv, */
/*                or provided to routine idzr_qrpiv */
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
} /* idz_rearr__ */





/* Subroutine */ int idz_rinqr__(integer *m, integer *n, doublecomplex *a, integer *krank, doublecomplex *r__){
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j, k;


/*       extracts R in the QR decomposition specified by the output a */
/*       of the routine idzr_qrpiv or idzp_qrpiv. */

/*       input: */
/*       m -- first dimension of a */
/*       n -- second dimension of a and r */
/*       a -- output of routine idzr_qrpiv or idzp_qrpiv */
/*       krank -- rank output by routine idzp_qrpiv (or specified */
/*                to routine idzr_qrpiv) */

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
} /* idz_rinqr__ */

