/* iddp_asvd.f -- translated by f2c (version 20160102).
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


/*       routine iddp_asvd computes the SVD, to a specified precision, */
/*       of an arbitrary matrix. This routine is randomized. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int iddp_asvd__(integer *lw, doublereal *eps, integer *m, integer *n, doublereal *a, doublereal *winit, integer *krank, integer *iu, integer *iv, integer *is, doublereal *w, integer *ier){
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int iddp_aid__(doublereal *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal * );    static integer k;
    extern /* Subroutine */ int iddp_asvd0__(integer *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, doublereal *) ;    static integer ls, lu, lv, lw2, isi, iui, ivi, icol, lcol, iproj, ilist, 
	    lproj, llist, iwork, lwork;


/*       constructs a rank-krank SVD  U Sigma V^T  approximating a */
/*       to precision eps, where U is an m x krank matrix whose */
/*       columns are orthonormal, V is an n x krank matrix whose */
/*       columns are orthonormal, and Sigma is a diagonal krank x krank */
/*       matrix whose entries are all nonnegative. */
/*       The entries of U are stored in w, starting at w(iu); */
/*       the entries of V are stored in w, starting at w(iv). */
/*       The diagonal entries of Sigma are stored in w, */
/*       starting at w(is). This routine uses a randomized algorithm. */

/*       input: */
/*       lw -- maximum usable length (in real*8 elements) */
/*             of the array w */
/*       eps -- precision of the desired approximation */
/*       m -- number of rows in a */
/*       n -- number of columns in a */
/*       a -- matrix to be approximated; the present routine does not */
/*            alter a */
/*       winit -- initialization array that has been constructed */
/*                by routine idd_frmi */

/*       output: */
/*       krank -- rank of the SVD constructed */
/*       iu -- index in w of the first entry of the matrix */
/*             of orthonormal left singular vectors of a */
/*       iv -- index in w of the first entry of the matrix */
/*             of orthonormal right singular vectors of a */
/*       is -- index in w of the first entry of the array */
/*             of singular values of a */
/*       w -- array containing the singular values and singular vectors */
/*            of a; w doubles as a work array, and so must be at least */
/*            max( (krank+1)*(3*m+5*n+1)+25*krank**2, (2*n+1)*(n2+1) ) */
/*            real*8 elements long, where n2 is the greatest integer */
/*            less than or equal to m, such that n2 is */
/*            a positive integer power of two; krank is the rank output */
/*            by this routine */
/*       ier -- 0 when the routine terminates successfully; */
/*              -1000 when lw is too small; */
/*              other nonzero values when idd_id2svd bombs */

/*       _N.B._: w must be at least */
/*               max( (krank+1)*(3*m+5*n+1)+25*krank^2, (2*n+1)*(n2+1) ) */
/*               real*8 elements long, where n2 is the greatest integer */
/*               less than or equal to m, such that n2 is */
/*               a positive integer power of two; */
/*               krank is the rank output by this routine. */
/*               Also, the algorithm used by this routine is randomized. */



/*       Allocate memory in w. */

    /* Parameter adjustments */
    --winit;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;

    /* Function Body */
    lw2 = 0;

    ilist = lw2 + 1;
    llist = *n;
    lw2 += llist;

    iproj = lw2 + 1;


/*       ID a. */

    iddp_aid__(eps, m, n, &a[a_offset], &winit[1], krank, &w[ilist], &w[iproj]
	    );


    if (*krank > 0) {


/*         Allocate more memory in w. */

	lproj = *krank * (*n - *krank);
	lw2 += lproj;

	icol = lw2 + 1;
	lcol = *m * *krank;
	lw2 += lcol;

	iui = lw2 + 1;
	lu = *m * *krank;
	lw2 += lu;

	ivi = lw2 + 1;
	lv = *n * *krank;
	lw2 += lv;

	isi = lw2 + 1;
	ls = *krank;
	lw2 += ls;

	iwork = lw2 + 1;
/* Computing 2nd power */
	i__1 = *krank;
	lwork = (*krank + 1) * (*m + *n * 3) + i__1 * i__1 * 26;
	lw2 += lwork;


	if (*lw < lw2) {
	    *ier = -1000;
	    return 0;
	}


	iddp_asvd0__(m, n, &a[a_offset], krank, &w[ilist], &w[iproj], &w[iui],
		 &w[ivi], &w[isi], ier, &w[icol], &w[iwork]);
	if (*ier != 0) {
	    return 0;
	}


	*iu = 1;
	*iv = *iu + lu;
	*is = *iv + lv;


/*         Copy the singular values and singular vectors */
/*         into their proper locations. */

	i__1 = lu;
	for (k = 1; k <= i__1; ++k) {
	    w[*iu + k - 1] = w[iui + k - 1];
	}

/* k */
	i__1 = lv;
	for (k = 1; k <= i__1; ++k) {
	    w[*iv + k - 1] = w[ivi + k - 1];
	}

/* k */
	i__1 = ls;
	for (k = 1; k <= i__1; ++k) {
	    w[*is + k - 1] = w[isi + k - 1];
	}


/* k */
    }


/* krank .gt. 0 */
    return 0;
} /* iddp_asvd__ */





/* Subroutine */ int iddp_asvd0__(integer *m, integer *n, doublereal *a, integer *krank, integer *list, doublereal *proj, doublereal *u, doublereal *v, doublereal *s, integer *ier, doublereal *col, doublereal *work){
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, proj_dim1, 
	    proj_offset, col_dim1, col_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int idd_id2svd__(integer *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *), idd_copycols__(integer *, integer *, doublereal *, integer *, integer *, doublereal *);

/*       routine iddp_asvd serves as a memory wrapper */
/*       for the present routine (please see routine iddp_asvd */
/*       for further documentation). */



/*       Collect together the columns of a indexed by list into col. */

/* Computing 2nd power */
    i__1 = *krank;
    /* Parameter adjustments */
    --list;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    col_dim1 = *m;
    col_offset = 1 + col_dim1;
    col -= col_offset;
    --s;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *m;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    proj_dim1 = *krank;
    proj_offset = 1 + proj_dim1;
    proj -= proj_offset;

    /* Function Body */
    idd_copycols__(m, n, &a[a_offset], krank, &list[1], &col[col_offset]);


/*       Convert the ID to an SVD. */

    idd_id2svd__(m, krank, &col[col_offset], n, &list[1], &proj[proj_offset], 
	    &u[u_offset], &v[v_offset], &s[1], ier, &work[1]);


    return 0;
} /* iddp_asvd0__ */

