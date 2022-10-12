/* idz_qrpiv.f -- translated by f2c (version 20160102).
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


/*       routine idzp_qrpiv computes the pivoted QR decomposition */
/*       of a matrix via Householder transformations, */
/*       stopping at a specified precision of the decomposition. */

/*       routine idzr_qrpiv computes the pivoted QR decomposition */
/*       of a matrix via Householder transformations, */
/*       stopping at a specified rank of the decomposition. */

/*       routine idz_qmatvec applies to a single vector */
/*       the Q matrix (or its adjoint) in the QR decomposition */
/*       of a matrix, as described by the output of idzp_qrpiv or */
/*       idzr_qrpiv. If you're concerned about efficiency and want */
/*       to apply Q (or its adjoint) to multiple vectors, */
/*       use idz_qmatmat instead. */

/*       routine idz_qmatmat applies */
/*       to multiple vectors collected together */
/*       as a matrix the Q matrix (or its adjoint) */
/*       in the QR decomposition of a matrix, as described */
/*       by the output of idzp_qrpiv. If you don't want to provide */
/*       a work array and want to apply Q (or its adjoint) */
/*       to a single vector, use idz_qmatvec instead. */

/*       routine idz_qinqr reconstructs the Q matrix */
/*       in a QR decomposition from the data generated by idzp_qrpiv */
/*       or idzr_qrpiv. */

/*       routine idz_permmult multiplies together a bunch */
/*       of permutations. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int idz_permmult__(integer *m, integer *ind, integer *n, integer *indprod){
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, iswap;


/*       multiplies together the series of permutations in ind. */

/*       input: */
/*       m -- length of ind */
/*       ind(k) -- number of the slot with which to swap */
/*                 the k^th slot */
/*       n -- length of indprod and indprodinv */

/*       output: */
/*       indprod -- product of the permutations in ind, */
/*                  with the permutation swapping 1 and ind(1) */
/*                  taken leftmost in the product, */
/*                  that swapping 2 and ind(2) taken next leftmost, */
/*                  ..., that swapping krank and ind(krank) */
/*                  taken rightmost; indprod(k) is the number */
/*                  of the slot with which to swap the k^th slot */
/*                  in the product permutation */



    /* Parameter adjustments */
    --ind;
    --indprod;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	indprod[k] = k;
    }

/* k */
    for (k = *m; k >= 1; --k) {

/*         Swap indprod(k) and indprod(ind(k)). */

	iswap = indprod[k];
	indprod[k] = indprod[ind[k]];
	indprod[ind[k]] = iswap;

    }


/* k */
    return 0;
} /* idz_permmult__ */





/* Subroutine */ int idz_qinqr__(integer *m, integer *n, doublecomplex *a, integer *krank, doublecomplex *q){
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer ifrescal, j, k, mm;
    static doublereal scal;
    extern /* Subroutine */ int idz_houseapp__(integer *, doublecomplex *, doublecomplex *, integer *, doublereal *, doublecomplex *);

/*       constructs the matrix q from idzp_qrpiv or idzr_qrpiv */
/*       (see the routine idzp_qrpiv or idzr_qrpiv */
/*       for more information). */

/*       input: */
/*       m -- first dimension of a; also, right now, q is m x m */
/*       n -- second dimension of a */
/*       a -- matrix output by idzp_qrpiv or idzr_qrpiv */
/*            (and denoted the same there) */
/*       krank -- numerical rank output by idzp_qrpiv or idzr_qrpiv */
/*                (and denoted the same there) */

/*       output: */
/*       q -- unitary matrix implicitly specified by the data in a */
/*            from idzp_qrpiv or idzr_qrpiv */

/*       Note: */
/*       Right now, this routine simply multiplies */
/*       one after another the krank Householder matrices */
/*       in the full QR decomposition of a, */
/*       in order to obtain the complete m x m Q factor in the QR. */
/*       This routine should instead use the following */
/*       (more elaborate but more efficient) scheme */
/*       to construct a q dimensioned q(krank,m); this scheme */
/*       was introduced by Robert Schreiber and Charles Van Loan */
/*       in "A Storage-Efficient _WY_ Representation */
/*       for Products of Householder Transformations," */
/*       _SIAM Journal on Scientific and Statistical Computing_, */
/*       Vol. 10, No. 1, pp. 53-57, January, 1989: */

/*       Theorem 1. Suppose that Q = _1_ + YTY^* is */
/*       an m x m unitary matrix, */
/*       where Y is an m x k matrix */
/*       and T is a k x k upper triangular matrix. */
/*       Suppose also that P = _1_ - 2 v v^* is */
/*       a Householder matrix and Q_+ = QP, */
/*       where v is an m x 1 real vector, */
/*       normalized so that v^* v = 1. */
/*       Then, Q_+ = _1_ + Y_+ T_+ Y_+^*, */
/*       where Y_+ = (Y v) is the m x (k+1) matrix */
/*       formed by adjoining v to the right of Y, */
/*                 ( T   z ) */
/*       and T_+ = (       ) is */
/*                 ( 0  -2 ) */
/*       the (k+1) x (k+1) upper triangular matrix */
/*       formed by adjoining z to the right of T */
/*       and the vector (0 ... 0 -2) with k zeroes below (T z), */
/*       where z = -2 T Y^* v. */

/*       Now, suppose that A is a (rank-deficient) matrix */
/*       whose complete QR decomposition has */
/*       the blockwise partioned form */
/*           ( Q_11 Q_12 ) ( R_11 R_12 )   ( Q_11 ) */
/*       A = (           ) (           ) = (      ) (R_11 R_12). */
/*           ( Q_21 Q_22 ) (  0    0   )   ( Q_21 ) */
/*       Then, the only blocks of the orthogonal factor */
/*       in the above QR decomposition of A that matter are */
/*                                                        ( Q_11 ) */
/*       Q_11 and Q_21, _i.e._, only the block of columns (      ) */
/*                                                        ( Q_21 ) */
/*       interests us. */
/*       Suppose in addition that Q_11 is a k x k matrix, */
/*       Q_21 is an (m-k) x k matrix, and that */
/*       ( Q_11 Q_12 ) */
/*       (           ) = _1_ + YTY^*, as in Theorem 1 above. */
/*       ( Q_21 Q_22 ) */
/*       Then, Q_11 = _1_ + Y_1 T Y_1^* */
/*       and Q_21 = Y_2 T Y_1^*, */
/*       where Y_1 is the k x k matrix and Y_2 is the (m-k) x k matrix */
/*                   ( Y_1 ) */
/*       so that Y = (     ). */
/*                   ( Y_2 ) */

/*       So, you can calculate T and Y via the above recursions, */
/*       and then use these to compute the desired Q_11 and Q_21. */




/*       Zero all of the entries of q. */

    /* Parameter adjustments */
    q_dim1 = *m;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = j + k * q_dim1;
	    q[i__3].r = 0., q[i__3].i = 0.;
	}
/* j */
    }


/*       Place 1's along the diagonal of q. */

/* k */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * q_dim1;
	q[i__2].r = 1., q[i__2].i = 0.;
    }


/*       Apply the krank Householder transformations stored in a. */

/* k */
    for (k = *krank; k >= 1; --k) {
	i__1 = *m;
	for (j = k; j <= i__1; ++j) {
	    mm = *m - k + 1;
	    ifrescal = 1;
	    if (k < *m) {
		idz_houseapp__(&mm, &a[k + 1 + k * a_dim1], &q[k + j * q_dim1]
			, &ifrescal, &scal, &q[k + j * q_dim1]);
	    }
	}
/* j */
    }


/* k */
    return 0;
} /* idz_qinqr__ */





/* Subroutine */ int idz_qmatvec__(integer *ifadjoint, integer *m, integer *n, doublecomplex *a, integer *krank, doublecomplex *v){
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer ifrescal, k, mm;
    static doublereal scal;
    extern /* Subroutine */ int idz_houseapp__(integer *, doublecomplex *, doublecomplex *, integer *, doublereal *, doublecomplex *);

/*       applies to a single vector the Q matrix (or its adjoint) */
/*       which the routine idzp_qrpiv or idzr_qrpiv has stored */
/*       in a triangle of the matrix it produces (stored, incidentally, */
/*       as data for applying a bunch of Householder reflections). */
/*       Use the routine idz_qmatmat to apply the Q matrix */
/*       (or its adjoint) */
/*       to a bunch of vectors collected together as a matrix, */
/*       if you're concerned about efficiency. */

/*       input: */
/*       ifadjoint -- set to 0 for applying Q; */
/*                    set to 1 for applying the adjoint of Q */
/*       m -- first dimension of a and length of v */
/*       n -- second dimension of a */
/*       a -- data describing the qr decomposition of a matrix, */
/*            as produced by idzp_qrpiv or idzr_qrpiv */
/*       krank -- numerical rank */
/*       v -- vector to which Q (or its adjoint) is to be applied */

/*       output: */
/*       v -- vector to which Q (or its adjoint) has been applied */



    /* Parameter adjustments */
    --v;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ifrescal = 1;


    if (*ifadjoint == 0) {

	for (k = *krank; k >= 1; --k) {
	    mm = *m - k + 1;
	    if (k < *m) {
		idz_houseapp__(&mm, &a[k + 1 + k * a_dim1], &v[k], &ifrescal, 
			&scal, &v[k]);
	    }
	}

/* k */
    }


    if (*ifadjoint == 1) {

	i__1 = *krank;
	for (k = 1; k <= i__1; ++k) {
	    mm = *m - k + 1;
	    if (k < *m) {
		idz_houseapp__(&mm, &a[k + 1 + k * a_dim1], &v[k], &ifrescal, 
			&scal, &v[k]);
	    }
	}

/* k */
    }


    return 0;
} /* idz_qmatvec__ */





/* Subroutine */ int idz_qmatmat__(integer *ifadjoint, integer *m, integer *n, doublecomplex *a, integer *krank, integer *l, doublecomplex *b, doublereal *work){
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer ifrescal, j, k, mm;
    extern /* Subroutine */ int idz_houseapp__(integer *, doublecomplex *, doublecomplex *, integer *, doublereal *, doublecomplex *);

/*       applies to a bunch of vectors collected together as a matrix */
/*       the Q matrix (or its adjoint) which the routine idzp_qrpiv */
/*       or idzr_qrpiv has stored in a triangle of the matrix */
/*       it produces (stored, incidentally, as data */
/*       for applying a bunch of Householder reflections). */
/*       Use the routine idz_qmatvec to apply the Q matrix */
/*       (or its adjoint) */
/*       to a single vector, if you'd rather not provide a work array. */

/*       input: */
/*       ifadjoint -- set to 0 for applying Q; */
/*                    set to 1 for applying the adjoint of Q */
/*       m -- first dimension of both a and b */
/*       n -- second dimension of a */
/*       a -- data describing the qr decomposition of a matrix, */
/*            as produced by idzp_qrpiv or idzr_qrpiv */
/*       krank -- numerical rank */
/*       l -- second dimension of b */
/*       b -- matrix to which Q (or its adjoint) is to be applied */

/*       output: */
/*       b -- matrix to which Q (or its adjoint) has been applied */

/*       work: */
/*       work -- must be at least krank real*8 elements long */



    /* Parameter adjustments */
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    b_dim1 = *m;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    if (*ifadjoint == 0) {


/*         Handle the first iteration, j = 1, */
/*         calculating all scals (ifrescal = 1). */

	ifrescal = 1;

	j = 1;

	for (k = *krank; k >= 1; --k) {
	    if (k < *m) {
		mm = *m - k + 1;
		idz_houseapp__(&mm, &a[k + 1 + k * a_dim1], &b[k + j * b_dim1]
			, &ifrescal, &work[k], &b[k + j * b_dim1]);
	    }
	}


/* k */
	if (*l > 1) {

/*           Handle the other iterations, j > 1, */
/*           using the scals just computed (ifrescal = 0). */

	    ifrescal = 0;

	    i__1 = *l;
	    for (j = 2; j <= i__1; ++j) {

		for (k = *krank; k >= 1; --k) {
		    if (k < *m) {
			mm = *m - k + 1;
			idz_houseapp__(&mm, &a[k + 1 + k * a_dim1], &b[k + j *
				 b_dim1], &ifrescal, &work[k], &b[k + j * 
				b_dim1]);
		    }
		}

/* k */
	    }

/* j */
	}


/* j .gt. 1 */
    }


/* ifadjoint .eq. 0 */
    if (*ifadjoint == 1) {


/*         Handle the first iteration, j = 1, */
/*         calculating all scals (ifrescal = 1). */

	ifrescal = 1;

	j = 1;

	i__1 = *krank;
	for (k = 1; k <= i__1; ++k) {
	    if (k < *m) {
		mm = *m - k + 1;
		idz_houseapp__(&mm, &a[k + 1 + k * a_dim1], &b[k + j * b_dim1]
			, &ifrescal, &work[k], &b[k + j * b_dim1]);
	    }
	}


/* k */
	if (*l > 1) {

/*           Handle the other iterations, j > 1, */
/*           using the scals just computed (ifrescal = 0). */

	    ifrescal = 0;

	    i__1 = *l;
	    for (j = 2; j <= i__1; ++j) {

		i__2 = *krank;
		for (k = 1; k <= i__2; ++k) {
		    if (k < *m) {
			mm = *m - k + 1;
			idz_houseapp__(&mm, &a[k + 1 + k * a_dim1], &b[k + j *
				 b_dim1], &ifrescal, &work[k], &b[k + j * 
				b_dim1]);
		    }
		}

/* k */
	    }

/* j */
	}


/* j .gt. 1 */
    }


/* ifadjoint .eq. 1 */
    return 0;
} /* idz_qmatmat__ */





/* Subroutine */ int idzp_qrpiv__(doublereal *eps, integer *m, integer *n, doublecomplex *a, integer *krank, integer *ind, doublereal *ss){
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    static integer ifrescal, j, k;
    extern /* Subroutine */ int idz_house__(integer *, doublecomplex *, doublecomplex *, doublecomplex *, doublereal *);    static integer mm;
    static doublereal scal, feps;
    static integer kpiv;
    extern /* Subroutine */ int idz_houseapp__(integer *, doublecomplex *, doublecomplex *, integer *, doublereal *, doublecomplex *);    static doublecomplex cswap;
    static doublereal ssmax, rswap;
    static integer nupdate;
    static doublereal ssmaxin;


/*       computes the pivoted QR decomposition */
/*       of the matrix input into a, using Householder transformations, */
/*       _i.e._, transforms the matrix a from its input value in */
/*       to the matrix out with entry */

/*                               m */
/*       out(j,indprod(k))  =  Sigma  q(l,j) * in(l,k), */
/*                              l=1 */

/*       for all j = 1, ..., krank, and k = 1, ..., n, */

/*       where in = the a from before the routine runs, */
/*       out = the a from after the routine runs, */
/*       out(j,k) = 0 when j > k (so that out is triangular), */
/*       q(1:m,1), ..., q(1:m,krank) are orthonormal, */
/*       indprod is the product of the permutations given by ind, */
/*       (as computable via the routine permmult, */
/*       with the permutation swapping 1 and ind(1) taken leftmost */
/*       in the product, that swapping 2 and ind(2) taken next leftmost, */
/*       ..., that swapping krank and ind(krank) taken rightmost), */
/*       and with the matrix out satisfying */

/*                   krank */
/*       in(j,k)  =  Sigma  q(j,l) * out(l,indprod(k))  +  epsilon(j,k), */
/*                    l=1 */

/*       for all j = 1, ..., m, and k = 1, ..., n, */

/*       for some matrix epsilon such that */
/*       the root-sum-square of the entries of epsilon */
/*       <= the root-sum-square of the entries of in * eps. */
/*       Well, technically, this routine outputs the Householder vectors */
/*       (or, rather, their second through last entries) */
/*       in the part of a that is supposed to get zeroed, that is, */
/*       in a(j,k) with m >= j > k >= 1. */

/*       input: */
/*       eps -- relative precision of the resulting QR decomposition */
/*       m -- first dimension of a and q */
/*       n -- second dimension of a */
/*       a -- matrix whose QR decomposition gets computed */

/*       output: */
/*       a -- triangular (R) factor in the QR decompositon */
/*            of the matrix input into the same storage locations, */
/*            with the Householder vectors stored in the part of a */
/*            that would otherwise consist entirely of zeroes, that is, */
/*            in a(j,k) with m >= j > k >= 1 */
/*       krank -- numerical rank */
/*       ind(k) -- index of the k^th pivot vector; */
/*                 the following code segment will correctly rearrange */
/*                 the product b of q and the upper triangle of out */
/*                 so that b matches the input matrix in */
/*                 to relative precision eps: */

/*                 copy the non-rearranged product of q and out into b */
/*                 set k to krank */
/*                 [start of loop] */
/*                   swap b(1:m,k) and b(1:m,ind(k)) */
/*                   decrement k by 1 */
/*                 if k > 0, then go to [start of loop] */

/*       work: */
/*       ss -- must be at least n real*8 words long */

/*       _N.B._: This routine outputs the Householder vectors */
/*       (or, rather, their second through last entries) */
/*       in the part of a that is supposed to get zeroed, that is, */
/*       in a(j,k) with m >= j > k >= 1. */

/*       reference: */
/*       Golub and Van Loan, "Matrix Computations," 3rd edition, */
/*            Johns Hopkins University Press, 1996, Chapter 5. */



    /* Parameter adjustments */
    --ss;
    --ind;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    feps = 1e-17;


/*       Compute the sum of squares of the entries in each column of a, */
/*       the maximum of all such sums, and find the first pivot */
/*       (column with the greatest such sum). */

    ssmax = 0.;
    kpiv = 1;

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {

	ss[k] = 0.;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = k;
	    i__4 = k;
	    i__5 = j + k * a_dim1;
	    d_cnjg(&z__3, &a[j + k * a_dim1]);
	    z__2.r = a[i__5].r * z__3.r - a[i__5].i * z__3.i, z__2.i = a[i__5]
		    .r * z__3.i + a[i__5].i * z__3.r;
	    z__1.r = ss[i__4] + z__2.r, z__1.i = z__2.i;
	    ss[i__3] = z__1.r;
	}

/* j */
	if (ss[k] > ssmax) {
	    ssmax = ss[k];
	    kpiv = k;
	}

    }

/* k */
    ssmaxin = ssmax;

    nupdate = 0;


/*       While ssmax > eps**2*ssmaxin, krank < m, and krank < n, */
/*       do the following block of code, */
/*       which ends at the statement labeled 2000. */

    *krank = 0;
L1000:

/* Computing 2nd power */
    d__1 = *eps;
    if (ssmax <= d__1 * d__1 * ssmaxin || *krank >= *m || *krank >= *n) {
	goto L2000;
    }
    ++(*krank);


    mm = *m - *krank + 1;


/*         Perform the pivoting. */

    ind[*krank] = kpiv;

/*         Swap a(1:m,krank) and a(1:m,kpiv). */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + *krank * a_dim1;
	cswap.r = a[i__2].r, cswap.i = a[i__2].i;
	i__2 = j + *krank * a_dim1;
	i__3 = j + kpiv * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = j + kpiv * a_dim1;
	a[i__2].r = cswap.r, a[i__2].i = cswap.i;
    }

/*         Swap ss(krank) and ss(kpiv). */

/* j */
    rswap = ss[*krank];
    ss[*krank] = ss[kpiv];
    ss[kpiv] = rswap;


    if (*krank < *m) {


/*           Compute the data for the Householder transformation */
/*           which will zero a(krank+1,krank), ..., a(m,krank) */
/*           when applied to a, replacing a(krank,krank) */
/*           with the first entry of the result of the application */
/*           of the Householder matrix to a(krank:m,krank), */
/*           and storing entries 2 to mm of the Householder vector */
/*           in a(krank+1,krank), ..., a(m,krank) */
/*           (which otherwise would get zeroed upon application */
/*           of the Householder transformation). */

	idz_house__(&mm, &a[*krank + *krank * a_dim1], &a[*krank + *krank * 
		a_dim1], &a[*krank + 1 + *krank * a_dim1], &scal);
	ifrescal = 0;


/*           Apply the Householder transformation */
/*           to the lower right submatrix of a */
/*           with upper leftmost entry at position (krank,krank+1). */

	if (*krank < *n) {
	    i__1 = *n;
	    for (k = *krank + 1; k <= i__1; ++k) {
		idz_houseapp__(&mm, &a[*krank + 1 + *krank * a_dim1], &a[*
			krank + k * a_dim1], &ifrescal, &scal, &a[*krank + k *
			 a_dim1]);
	    }
/* k */
	}


/*           Update the sums-of-squares array ss. */

	i__1 = *n;
	for (k = *krank; k <= i__1; ++k) {
	    i__2 = k;
	    i__3 = k;
	    i__4 = *krank + k * a_dim1;
	    d_cnjg(&z__3, &a[*krank + k * a_dim1]);
	    z__2.r = a[i__4].r * z__3.r - a[i__4].i * z__3.i, z__2.i = a[i__4]
		    .r * z__3.i + a[i__4].i * z__3.r;
	    z__1.r = ss[i__3] - z__2.r, z__1.i = -z__2.i;
	    ss[i__2] = z__1.r;
	}


/*           Find the pivot (column with the greatest sum of squares */
/*           of its entries). */

/* k */
	ssmax = 0.;
	kpiv = *krank + 1;

	if (*krank < *n) {

	    i__1 = *n;
	    for (k = *krank + 1; k <= i__1; ++k) {

		if (ss[k] > ssmax) {
		    ssmax = ss[k];
		    kpiv = k;
		}

	    }

/* k */
	}


/*           Recompute the sums-of-squares and the pivot */
/*           when ssmax first falls below */
/*           sqrt((1000*feps)^2) * ssmaxin */
/*           and when ssmax first falls below */
/*           ((1000*feps)^2) * ssmaxin. */

/* krank .lt. n */
/* Computing 2nd power */
	d__1 = feps * 1000;
/* Computing 2nd power */
	d__2 = feps * 1000;
	if (ssmax < sqrt(d__1 * d__1) * ssmaxin && nupdate == 0 || ssmax < 
		d__2 * d__2 * ssmaxin && nupdate == 1) {

	    ++nupdate;

	    ssmax = 0.;
	    kpiv = *krank + 1;

	    if (*krank < *n) {

		i__1 = *n;
		for (k = *krank + 1; k <= i__1; ++k) {

		    ss[k] = 0.;
		    i__2 = *m;
		    for (j = *krank + 1; j <= i__2; ++j) {
			i__3 = k;
			i__4 = k;
			i__5 = j + k * a_dim1;
			d_cnjg(&z__3, &a[j + k * a_dim1]);
			z__2.r = a[i__5].r * z__3.r - a[i__5].i * z__3.i, 
				z__2.i = a[i__5].r * z__3.i + a[i__5].i * 
				z__3.r;
			z__1.r = ss[i__4] + z__2.r, z__1.i = z__2.i;
			ss[i__3] = z__1.r;
		    }

/* j */
		    if (ss[k] > ssmax) {
			ssmax = ss[k];
			kpiv = k;
		    }

		}

/* k */
	    }

/* krank .lt. n */
	}


    }


/* krank .lt. m */
    goto L1000;
L2000:


    return 0;
} /* idzp_qrpiv__ */





/* Subroutine */ int idzr_qrpiv__(integer *m, integer *n, doublecomplex *a, integer *krank, integer *ind, doublereal *ss){
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    static integer ifrescal, j, k;
    extern /* Subroutine */ int idz_house__(integer *, doublecomplex *, doublecomplex *, doublecomplex *, doublereal *);    static integer mm;
    static doublereal scal, feps;
    static integer kpiv, loop;
    extern /* Subroutine */ int idz_houseapp__(integer *, doublecomplex *, doublecomplex *, integer *, doublereal *, doublecomplex *);    static doublecomplex cswap;
    static doublereal ssmax;
    static integer loops;
    static doublereal rswap;
    static integer nupdate;
    static doublereal ssmaxin;


/*       computes the pivoted QR decomposition */
/*       of the matrix input into a, using Householder transformations, */
/*       _i.e._, transforms the matrix a from its input value in */
/*       to the matrix out with entry */

/*                               m */
/*       out(j,indprod(k))  =  Sigma  q(l,j) * in(l,k), */
/*                              l=1 */

/*       for all j = 1, ..., krank, and k = 1, ..., n, */

/*       where in = the a from before the routine runs, */
/*       out = the a from after the routine runs, */
/*       out(j,k) = 0 when j > k (so that out is triangular), */
/*       q(1:m,1), ..., q(1:m,krank) are orthonormal, */
/*       indprod is the product of the permutations given by ind, */
/*       (as computable via the routine permmult, */
/*       with the permutation swapping 1 and ind(1) taken leftmost */
/*       in the product, that swapping 2 and ind(2) taken next leftmost, */
/*       ..., that swapping krank and ind(krank) taken rightmost), */
/*       and with the matrix out satisfying */

/*                  min(m,n,krank) */
/*       in(j,k)  =     Sigma      q(j,l) * out(l,indprod(k)) */
/*                       l=1 */

/*                +  epsilon(j,k), */

/*       for all j = 1, ..., m, and k = 1, ..., n, */

/*       for some matrix epsilon whose norm is (hopefully) minimized */
/*       by the pivoting procedure. */
/*       Well, technically, this routine outputs the Householder vectors */
/*       (or, rather, their second through last entries) */
/*       in the part of a that is supposed to get zeroed, that is, */
/*       in a(j,k) with m >= j > k >= 1. */

/*       input: */
/*       m -- first dimension of a and q */
/*       n -- second dimension of a */
/*       a -- matrix whose QR decomposition gets computed */
/*       krank -- desired rank of the output matrix */
/*                (please note that if krank > m or krank > n, */
/*                then the rank of the output matrix will be */
/*                less than krank) */

/*       output: */
/*       a -- triangular (R) factor in the QR decompositon */
/*            of the matrix input into the same storage locations, */
/*            with the Householder vectors stored in the part of a */
/*            that would otherwise consist entirely of zeroes, that is, */
/*            in a(j,k) with m >= j > k >= 1 */
/*       ind(k) -- index of the k^th pivot vector; */
/*                 the following code segment will correctly rearrange */
/*                 the product b of q and the upper triangle of out */
/*                 so that b matches the input matrix in */
/*                 to relative precision eps: */

/*                 copy the non-rearranged product of q and out into b */
/*                 set k to krank */
/*                 [start of loop] */
/*                   swap b(1:m,k) and b(1:m,ind(k)) */
/*                   decrement k by 1 */
/*                 if k > 0, then go to [start of loop] */

/*       work: */
/*       ss -- must be at least n real*8 words long */

/*       _N.B._: This routine outputs the Householder vectors */
/*       (or, rather, their second through last entries) */
/*       in the part of a that is supposed to get zeroed, that is, */
/*       in a(j,k) with m >= j > k >= 1. */

/*       reference: */
/*       Golub and Van Loan, "Matrix Computations," 3rd edition, */
/*            Johns Hopkins University Press, 1996, Chapter 5. */



    /* Parameter adjustments */
    --ss;
    --ind;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    feps = 1e-17;


/*       Compute the sum of squares of the entries in each column of a, */
/*       the maximum of all such sums, and find the first pivot */
/*       (column with the greatest such sum). */

    ssmax = 0.;
    kpiv = 1;

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {

	ss[k] = 0.;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = k;
	    i__4 = k;
	    i__5 = j + k * a_dim1;
	    d_cnjg(&z__3, &a[j + k * a_dim1]);
	    z__2.r = a[i__5].r * z__3.r - a[i__5].i * z__3.i, z__2.i = a[i__5]
		    .r * z__3.i + a[i__5].i * z__3.r;
	    z__1.r = ss[i__4] + z__2.r, z__1.i = z__2.i;
	    ss[i__3] = z__1.r;
	}

/* j */
	if (ss[k] > ssmax) {
	    ssmax = ss[k];
	    kpiv = k;
	}

    }

/* k */
    ssmaxin = ssmax;

    nupdate = 0;


/*       Set loops = min(krank,m,n). */

    loops = *krank;
    if (*m < loops) {
	loops = *m;
    }
    if (*n < loops) {
	loops = *n;
    }

    i__1 = loops;
    for (loop = 1; loop <= i__1; ++loop) {


	mm = *m - loop + 1;


/*         Perform the pivoting. */

	ind[loop] = kpiv;

/*         Swap a(1:m,loop) and a(1:m,kpiv). */

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = j + loop * a_dim1;
	    cswap.r = a[i__3].r, cswap.i = a[i__3].i;
	    i__3 = j + loop * a_dim1;
	    i__4 = j + kpiv * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = j + kpiv * a_dim1;
	    a[i__3].r = cswap.r, a[i__3].i = cswap.i;
	}

/*         Swap ss(loop) and ss(kpiv). */

/* j */
	rswap = ss[loop];
	ss[loop] = ss[kpiv];
	ss[kpiv] = rswap;


	if (loop < *m) {


/*           Compute the data for the Householder transformation */
/*           which will zero a(loop+1,loop), ..., a(m,loop) */
/*           when applied to a, replacing a(loop,loop) */
/*           with the first entry of the result of the application */
/*           of the Householder matrix to a(loop:m,loop), */
/*           and storing entries 2 to mm of the Householder vector */
/*           in a(loop+1,loop), ..., a(m,loop) */
/*           (which otherwise would get zeroed upon application */
/*           of the Householder transformation). */

	    idz_house__(&mm, &a[loop + loop * a_dim1], &a[loop + loop * 
		    a_dim1], &a[loop + 1 + loop * a_dim1], &scal);
	    ifrescal = 0;


/*           Apply the Householder transformation */
/*           to the lower right submatrix of a */
/*           with upper leftmost entry at position (loop,loop+1). */

	    if (loop < *n) {
		i__2 = *n;
		for (k = loop + 1; k <= i__2; ++k) {
		    idz_houseapp__(&mm, &a[loop + 1 + loop * a_dim1], &a[loop 
			    + k * a_dim1], &ifrescal, &scal, &a[loop + k * 
			    a_dim1]);
		}
/* k */
	    }


/*           Update the sums-of-squares array ss. */

	    i__2 = *n;
	    for (k = loop; k <= i__2; ++k) {
		i__3 = k;
		i__4 = k;
		i__5 = loop + k * a_dim1;
		d_cnjg(&z__3, &a[loop + k * a_dim1]);
		z__2.r = a[i__5].r * z__3.r - a[i__5].i * z__3.i, z__2.i = a[
			i__5].r * z__3.i + a[i__5].i * z__3.r;
		z__1.r = ss[i__4] - z__2.r, z__1.i = -z__2.i;
		ss[i__3] = z__1.r;
	    }


/*           Find the pivot (column with the greatest sum of squares */
/*           of its entries). */

/* k */
	    ssmax = 0.;
	    kpiv = loop + 1;

	    if (loop < *n) {

		i__2 = *n;
		for (k = loop + 1; k <= i__2; ++k) {

		    if (ss[k] > ssmax) {
			ssmax = ss[k];
			kpiv = k;
		    }

		}

/* k */
	    }


/*           Recompute the sums-of-squares and the pivot */
/*           when ssmax first falls below */
/*           sqrt((1000*feps)^2) * ssmaxin */
/*           and when ssmax first falls below */
/*           ((1000*feps)^2) * ssmaxin. */

/* loop .lt. n */
/* Computing 2nd power */
	    d__1 = feps * 1000;
/* Computing 2nd power */
	    d__2 = feps * 1000;
	    if (ssmax < sqrt(d__1 * d__1) * ssmaxin && nupdate == 0 || ssmax <
		     d__2 * d__2 * ssmaxin && nupdate == 1) {

		++nupdate;

		ssmax = 0.;
		kpiv = loop + 1;

		if (loop < *n) {

		    i__2 = *n;
		    for (k = loop + 1; k <= i__2; ++k) {

			ss[k] = 0.;
			i__3 = *m;
			for (j = loop + 1; j <= i__3; ++j) {
			    i__4 = k;
			    i__5 = k;
			    i__6 = j + k * a_dim1;
			    d_cnjg(&z__3, &a[j + k * a_dim1]);
			    z__2.r = a[i__6].r * z__3.r - a[i__6].i * z__3.i, 
				    z__2.i = a[i__6].r * z__3.i + a[i__6].i * 
				    z__3.r;
			    z__1.r = ss[i__5] + z__2.r, z__1.i = z__2.i;
			    ss[i__4] = z__1.r;
			}

/* j */
			if (ss[k] > ssmax) {
			    ssmax = ss[k];
			    kpiv = k;
			}

		    }

/* k */
		}

/* loop .lt. n */
	    }


	}


/* loop .lt. m */
    }


/* loop */
    return 0;
} /* idzr_qrpiv__ */
