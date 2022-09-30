/* idz_frm.f -- translated by f2c (version 20160102).
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

/*       this file contains the following user-callable routines: */


/*       routine idz_frm transforms a vector via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */

/*       routine idz_sfrm transforms a vector into a vector */
/*       of specified length via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */

/*       routine idz_frmi initializes routine idz_frm. */

/*       routine idz_sfrmi initializes routine idz_sfrm. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int idz_frm__(integer *m, integer *n, doublecomplex *w, doublecomplex *x, doublecomplex *y){
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int idz_random_transf__(doublecomplex *, doublecomplex *, doublecomplex *);    static integer iw;
    extern /* Subroutine */ int idz_permute__(integer *, integer *, doublecomplex *, doublecomplex *), zfftf_(integer *, doublecomplex *, doublecomplex *), idz_subselect__(integer *, integer *, integer *, doublecomplex *, doublecomplex *);

/*       transforms x into y via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */
/*       In contrast to routine idz_sfrm, the present routine works best */
/*       when the length of the transformed vector is the integer n */
/*       output by routine idz_frmi, or when the length */
/*       is not specified, but instead determined a posteriori */
/*       using the output of the present routine. The transformed vector */
/*       output by the present routine is randomly permuted. */

/*       input: */
/*       m -- length of x */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m, as obtained */
/*            from the routine idz_frmi; n is the length of y */
/*       w -- initialization array constructed by routine idz_frmi */
/*       x -- vector to be transformed */

/*       output: */
/*       y -- transform of x */

/*       reference: */
/*       Halko, Martinsson, Tropp, "Finding structure with randomness: */
/*            probabilistic algorithms for constructing approximate */
/*            matrix decompositions," SIAM Review, 53 (2): 217-288, */
/*            2011. */



/*       Apply Rokhlin's random transformation to x, obtaining */
/*       w(16*m+71 : 17*m+70). */

    /* Parameter adjustments */
    --x;
    --w;
    --y;

    /* Function Body */
    i__1 = *m + 3 + *n;
    iw = (integer) w[i__1].r;
    idz_random_transf__(&x[1], &w[(*m << 4) + 71], &w[iw]);


/*       Subselect from  w(16*m+71 : 17*m+70)  to obtain y. */

    idz_subselect__(n, &w[3], m, &w[(*m << 4) + 71], &y[1]);


/*       Copy y into  w(16*m+71 : 16*m+n+70). */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = (*m << 4) + 70 + k;
	i__3 = k;
	w[i__2].r = y[i__3].r, w[i__2].i = y[i__3].i;
    }


/*       Fourier transform  w(16*m+71 : 16*m+n+70). */

/* k */
    zfftf_(n, &w[(*m << 4) + 71], &w[*m + 4 + *n]);


/*       Permute  w(16*m+71 : 16*m+n+70)  to obtain y. */

    idz_permute__(n, &w[*m + 3], &w[(*m << 4) + 71], &y[1]);


    return 0;
} /* idz_frm__ */





/* Subroutine */ int idz_sfrm__(integer *l, integer *m, integer *n, doublecomplex *w, doublecomplex *x, doublecomplex *y){
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int idz_sfft__(integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *), idz_random_transf__( doublecomplex *, doublecomplex *, doublecomplex *);    static integer iw;
    extern /* Subroutine */ int idz_subselect__(integer *, integer *, integer *, doublecomplex *, doublecomplex *);

/*       transforms x into y via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */
/*       In contrast to routine idz_frm, the present routine works best */
/*       when the length l of the transformed vector is known a priori. */

/*       input: */
/*       l -- length of y; l must be less than or equal to n */
/*       m -- length of x */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m, as obtained */
/*            from the routine idz_frmi */
/*       w -- initialization array constructed by routine idz_sfrmi */
/*       x -- vector to be transformed */

/*       output: */
/*       y -- transform of x */

/*       _N.B._: l must be less than or equal to n. */

/*       reference: */
/*       Halko, Martinsson, Tropp, "Finding structure with randomness: */
/*            probabilistic algorithms for constructing approximate */
/*            matrix decompositions," SIAM Review, 53 (2): 217-288, */
/*            2011. */



/*       Apply Rokhlin's random transformation to x, obtaining */
/*       w(19*m+71 : 20*m+70). */

    /* Parameter adjustments */
    --y;
    --x;
    --w;

    /* Function Body */
    i__1 = *m + 4 + *l;
    iw = (integer) w[i__1].r;
    idz_random_transf__(&x[1], &w[*m * 19 + 71], &w[iw]);


/*       Subselect from  w(19*m+71 : 20*m+70)  to obtain */
/*       w(20*m+71 : 20*m+n+70). */

    idz_subselect__(n, &w[4], m, &w[*m * 19 + 71], &w[*m * 20 + 71]);


/*       Fourier transform  w(20*m+71 : 20*m+n+70). */

    idz_sfft__(l, &w[*m + 4], n, &w[*m + 5 + *l], &w[*m * 20 + 71]);


/*       Copy the desired entries from  w(20*m+71 : 20*m+n+70) */
/*       to y. */

    idz_subselect__(l, &w[*m + 4], n, &w[*m * 20 + 71], &y[1]);


    return 0;
} /* idz_sfrm__ */





/* Subroutine */ int idz_permute__(integer *n, integer *ind, doublecomplex *x, doublecomplex *y){
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k;


/*       copy the entries of x into y, rearranged according */
/*       to the permutation specified by ind. */

/*       input: */
/*       n -- length of ind, x, and y */
/*       ind -- permutation of n objects */
/*       x -- vector to be permuted */

/*       output: */
/*       y -- permutation of x */



    /* Parameter adjustments */
    --y;
    --x;
    --ind;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = ind[k];
	y[i__2].r = x[i__3].r, y[i__2].i = x[i__3].i;
    }


/* k */
    return 0;
} /* idz_permute__ */





/* Subroutine */ int idz_subselect__(integer *n, integer *ind, integer *m, doublecomplex *x, doublecomplex *y){
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k;


/*       copies into y the entries of x indicated by ind. */

/*       input: */
/*       n -- number of entries of x to copy into y */
/*       ind -- indices of the entries in x to copy into y */
/*       m -- length of x */
/*       x -- vector whose entries are to be copied */

/*       output: */
/*       y -- collection of entries of x specified by ind */



    /* Parameter adjustments */
    --y;
    --ind;
    --x;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = ind[k];
	y[i__2].r = x[i__3].r, y[i__2].i = x[i__3].i;
    }


/* k */
    return 0;
} /* idz_subselect__ */





/* Subroutine */ int idz_frmi__(integer *m, integer *n, doublecomplex *w)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer l, ia, lw;
    extern /* Subroutine */ int id_randperm__(integer *, doublecomplex *);
    static integer keep;
    extern /* Subroutine */ int prinf_(char *, integer *, integer *, ftnlen), zffti_(integer *, doublecomplex *);    static integer nsteps;
    extern /* Subroutine */ int idz_poweroftwo__(integer *, integer *, integer *), idz_random_transf_init__(integer *, integer *, doublecomplex *, integer *);

/*       initializes data for the routine idz_frm. */

/*       input: */
/*       m -- length of the vector to be transformed */

/*       output: */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m */
/*       w -- initialization array to be used by routine idz_frm */


/*       glossary for the fully initialized w: */

/*       w(1) = m */
/*       w(2) = n */
/*       w(3:2+m) stores a permutation of m objects */
/*       w(3+m:2+m+n) stores a permutation of n objects */
/*       w(3+m+n) = address in w of the initialization array */
/*                  for idz_random_transf */
/*       w(4+m+n:int(w(3+m+n))-1) stores the initialization array */
/*                                for zfft */
/*       w(int(w(3+m+n)):16*m+70) stores the initialization array */
/*                                for idz_random_transf */


/*       _N.B._: n is an output of the present routine; */
/*               this routine changes n. */




/*       Find the greatest integer less than or equal to m */
/*       which is a power of two. */

    /* Parameter adjustments */
    --w;

    /* Function Body */
    idz_poweroftwo__(m, &l, n);


/*       Store m and n in w. */

    w[1].r = (doublereal) (*m), w[1].i = 0.;
    w[2].r = (doublereal) (*n), w[2].i = 0.;


/*       Store random permutations of m and n objects in w. */

    id_randperm__(m, &w[3]);
    id_randperm__(n, &w[*m + 3]);


/*       Store the address within w of the idz_random_transf_init */
/*       initialization data. */

    ia = *m + 4 + *n + (*n << 1) + 15;
    i__1 = *m + 3 + *n;
    w[i__1].r = (doublereal) ia, w[i__1].i = 0.;


/*       Store the initialization data for zfft in w. */

    zffti_(n, &w[*m + 4 + *n]);


/*       Store the initialization data for idz_random_transf_init in w. */

    nsteps = 3;
    idz_random_transf_init__(&nsteps, m, &w[ia], &keep);


/*       Calculate the total number of elements used in w. */

    lw = *m + 3 + *n + (*n << 1) + 15 + nsteps * 3 * *m + (*m << 1) + *m / 4 
	    + 50;

    if ((*m << 4) + 70 < lw) {
	prinf_("lw = *", &lw, &c__1, (ftnlen)6);
	i__1 = (*m << 4) + 70;
	prinf_("16m+70 = *", &i__1, &c__1, (ftnlen)10);
	s_stop("", (ftnlen)0);
    }


    return 0;
} /* idz_frmi__ */





/* Subroutine */ int idz_sfrmi__(integer *l, integer *m, integer *n, doublecomplex *w){
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int idz_sffti__(integer *, doublecomplex *, integer *, doublecomplex *);    static integer ia, lw;
    extern /* Subroutine */ int id_randperm__(integer *, doublecomplex *);
    static integer keep;
    extern /* Subroutine */ int prinf_(char *, integer *, integer *, ftnlen);
    static integer idummy, nsteps;
    extern /* Subroutine */ int idz_poweroftwo__(integer *, integer *, integer *), idz_random_transf_init__(integer *, integer *, doublecomplex *, integer *);

/*       initializes data for the routine idz_sfrm. */

/*       input: */
/*       l -- length of the transformed (output) vector */
/*       m -- length of the vector to be transformed */

/*       output: */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m */
/*       w -- initialization array to be used by routine idz_sfrm */


/*       glossary for the fully initialized w: */

/*       w(1) = m */
/*       w(2) = n */
/*       w(3) is unused */
/*       w(4:3+m) stores a permutation of m objects */
/*       w(4+m:3+m+l) stores the indices of the l outputs which idz_sfft */
/*                    calculates */
/*       w(4+m+l) = address in w of the initialization array */
/*                  for idz_random_transf */
/*       w(5+m+l:int(w(4+m+l))-1) stores the initialization array */
/*                                for idz_sfft */
/*       w(int(w(4+m+l)):19*m+70) stores the initialization array */
/*                                for idz_random_transf */


/*       _N.B._: n is an output of the present routine; */
/*               this routine changes n. */




/*       Find the greatest integer less than or equal to m */
/*       which is a power of two. */

    /* Parameter adjustments */
    --w;

    /* Function Body */
    idz_poweroftwo__(m, &idummy, n);


/*       Store m and n in w. */

    w[1].r = (doublereal) (*m), w[1].i = 0.;
    w[2].r = (doublereal) (*n), w[2].i = 0.;
    w[3].r = 0., w[3].i = 0.;


/*       Store random permutations of m and n objects in w. */

    id_randperm__(m, &w[4]);
    id_randperm__(n, &w[*m + 4]);


/*       Store the address within w of the idz_random_transf_init */
/*       initialization data. */

    ia = *m + 5 + *l + (*l << 1) + 15 + *n * 3;
    i__1 = *m + 4 + *l;
    w[i__1].r = (doublereal) ia, w[i__1].i = 0.;


/*       Store the initialization data for idz_sfft in w. */

    idz_sffti__(l, &w[*m + 4], n, &w[*m + 5 + *l]);


/*       Store the initialization data for idz_random_transf_init in w. */

    nsteps = 3;
    idz_random_transf_init__(&nsteps, m, &w[ia], &keep);


/*       Calculate the total number of elements used in w. */

    lw = *m + 4 + *l + (*l << 1) + 15 + *n * 3 + nsteps * 3 * *m + (*m << 1) 
	    + *m / 4 + 50;

    if (*m * 19 + 70 < lw) {
	prinf_("lw = *", &lw, &c__1, (ftnlen)6);
	i__1 = *m * 19 + 70;
	prinf_("19m+70 = *", &i__1, &c__1, (ftnlen)10);
	s_stop("", (ftnlen)0);
    }


    return 0;
} /* idz_sfrmi__ */





/* Subroutine */ int idz_poweroftwo__(integer *m, integer *l, integer *n)
{

/*       computes l = floor(log_2(m)) and n = 2**l. */

/*       input: */
/*       m -- integer whose log_2 is to be taken */

/*       output: */
/*       l -- floor(log_2(m)) */
/*       n -- 2**l */



    *l = 0;
    *n = 1;

L1000:
    ++(*l);
    *n <<= 1;
    if (*n <= *m) {
	goto L1000;
    }

    --(*l);
    *n /= 2;


    return 0;
} /* idz_poweroftwo__ */

