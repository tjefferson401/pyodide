/* idd_frm.f -- translated by f2c (version 20160102).
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


/*       routine idd_frm transforms a vector via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */

/*       routine idd_sfrm transforms a vector into a vector */
/*       of specified length via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */

/*       routine idd_frmi initializes routine idd_frm. */

/*       routine idd_sfrmi initializes routine idd_sfrm. */

/*       routine idd_pairsamps calculates the indices of the pairs */
/*       of integers to which the individual integers */
/*       in a specified set belong. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int idd_frm__(integer *m, integer *n, doublereal *w, doublereal *x, doublereal *y){
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int idd_random_transf__(doublereal *, doublereal * , doublereal *);    static integer iw;
    extern /* Subroutine */ int idd_permute__(integer *, integer *, doublereal *, doublereal *), dfftf_(integer *, doublereal *, doublereal *), idd_subselect__(integer *, integer *, integer *, doublereal *, doublereal *);

/*       transforms x into y via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */
/*       In contrast to routine idd_sfrm, the present routine works best */
/*       when the length of the transformed vector is the integer n */
/*       output by routine idd_frmi, or when the length */
/*       is not specified, but instead determined a posteriori */
/*       using the output of the present routine. The transformed vector */
/*       output by the present routine is randomly permuted. */

/*       input: */
/*       m -- length of x */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m, as obtained */
/*            from the routine idd_frmi; n is the length of y */
/*       w -- initialization array constructed by routine idd_frmi */
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
    iw = (integer) w[*m + 3 + *n];
    idd_random_transf__(&x[1], &w[(*m << 4) + 71], &w[iw]);


/*       Subselect from  w(16*m+71 : 17*m+70)  to obtain y. */

    idd_subselect__(n, &w[3], m, &w[(*m << 4) + 71], &y[1]);


/*       Copy y into  w(16*m+71 : 16*m+n+70). */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[(*m << 4) + 70 + k] = y[k];
    }


/*       Fourier transform  w(16*m+71 : 16*m+n+70). */

/* k */
    dfftf_(n, &w[(*m << 4) + 71], &w[*m + 4 + *n]);


/*       Permute  w(16*m+71 : 16*m+n+70)  to obtain y. */

    idd_permute__(n, &w[*m + 3], &w[(*m << 4) + 71], &y[1]);


    return 0;
} /* idd_frm__ */





/* Subroutine */ int idd_sfrm__(integer *l, integer *m, integer *n, doublereal *w, doublereal *x, doublereal *y){
    extern /* Subroutine */ int idd_sfft__(integer *, doublereal *, integer *, doublereal *, doublereal *);    static integer l2;
    extern /* Subroutine */ int idd_random_transf__(doublereal *, doublereal * , doublereal *);    static integer iw;
    extern /* Subroutine */ int idd_subselect__(integer *, integer *, integer *, doublereal *, doublereal *);

/*       transforms x into y via a composition */
/*       of Rokhlin's random transform, random subselection, and an FFT. */
/*       In contrast to routine idd_frm, the present routine works best */
/*       when the length l of the transformed vector is known a priori. */

/*       input: */
/*       l -- length of y; l must be less than or equal to n */
/*       m -- length of x */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m, as obtained */
/*            from the routine idd_sfrmi */
/*       w -- initialization array constructed by routine idd_sfrmi */
/*       x -- vector to be transformed */

/*       output: */
/*       y -- transform of x */

/*       _N.B._: l must be less than or equal to n. */

/*       reference: */
/*       Halko, Martinsson, Tropp, "Finding structure with randomness: */
/*            probabilistic algorithms for constructing approximate */
/*            matrix decompositions," SIAM Review, 53 (2): 217-288, */
/*            2011. */



/*       Retrieve the number of pairs of outputs to be calculated */
/*       via sfft. */

    /* Parameter adjustments */
    --y;
    --x;
    --w;

    /* Function Body */
    l2 = (integer) w[3];


/*       Apply Rokhlin's random transformation to x, obtaining */
/*       w(25*m+91 : 26*m+90). */

    iw = (integer) w[*m + 4 + *l + l2];
    idd_random_transf__(&x[1], &w[*m * 25 + 91], &w[iw]);


/*       Subselect from  w(25*m+91 : 26*m+90)  to obtain */
/*       w(26*m+91 : 26*m+n+90). */

    idd_subselect__(n, &w[4], m, &w[*m * 25 + 91], &w[*m * 26 + 91]);


/*       Fourier transform  w(26*m+91 : 26*m+n+90). */

    idd_sfft__(&l2, &w[*m + 4 + *l], n, &w[*m + 5 + *l + l2], &w[*m * 26 + 91]
	    );


/*       Copy the desired entries from  w(26*m+91 : 26*m+n+90) */
/*       to y. */

    idd_subselect__(l, &w[*m + 4], n, &w[*m * 26 + 91], &y[1]);


    return 0;
} /* idd_sfrm__ */





/* Subroutine */ int idd_pairsamps__(integer *n, integer *l, integer *ind, integer *l2, integer *ind2, integer *marker){
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;


/*       calculates the indices of the l2 pairs of integers */
/*       to which the l individual integers from ind belong. */
/*       The integers in ind may range from 1 to n. */

/*       input: */
/*       n -- upper bound on the integers in ind */
/*            (the number 1 must be a lower bound); */
/*            n must be even */
/*       l -- length of ind */
/*       ind -- integers selected from 1 to n */

/*       output: */
/*       l2 -- length of ind2 */
/*       ind2 -- indices in the range from 1 to n/2 of the pairs */
/*               of integers to which the entries of ind belong */

/*       work: */
/*       marker -- must be at least n/2 integer elements long */

/*       _N.B._: n must be even. */



/*       Unmark all pairs. */

    /* Parameter adjustments */
    --marker;
    --ind2;
    --ind;

    /* Function Body */
    i__1 = *n / 2;
    for (k = 1; k <= i__1; ++k) {
	marker[k] = 0;
    }


/*       Mark the required pairs. */

/* k */
    i__1 = *l;
    for (k = 1; k <= i__1; ++k) {
	++marker[(ind[k] + 1) / 2];
    }


/*       Record the required pairs in indpair. */

/* k */
    *l2 = 0;

    i__1 = *n / 2;
    for (k = 1; k <= i__1; ++k) {

	if (marker[k] != 0) {
	    ++(*l2);
	    ind2[*l2] = k;
	}

    }


/* k */
    return 0;
} /* idd_pairsamps__ */





/* Subroutine */ int idd_permute__(integer *n, integer *ind, doublereal *x, doublereal *y){
    /* System generated locals */
    integer i__1;

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
	y[k] = x[ind[k]];
    }


/* k */
    return 0;
} /* idd_permute__ */





/* Subroutine */ int idd_subselect__(integer *n, integer *ind, integer *m, doublereal *x, doublereal *y){
    /* System generated locals */
    integer i__1;

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
	y[k] = x[ind[k]];
    }


/* k */
    return 0;
} /* idd_subselect__ */





/* Subroutine */ int idd_frmi__(integer *m, integer *n, doublereal *w)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer l, ia, lw;
    extern /* Subroutine */ int id_randperm__(integer *, doublereal *);
    static integer keep;
    extern /* Subroutine */ int dffti_(integer *, doublereal *), prinf_(char * , integer *, integer *, ftnlen);    static integer nsteps;
    extern /* Subroutine */ int idd_poweroftwo__(integer *, integer *, integer *), idd_random_transf_init__(integer *, integer *, doublereal *, integer *);

/*       initializes data for the routine idd_frm. */

/*       input: */
/*       m -- length of the vector to be transformed */

/*       output: */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m */
/*       w -- initialization array to be used by routine idd_frm */


/*       glossary for the fully initialized w: */

/*       w(1) = m */
/*       w(2) = n */
/*       w(3:2+m) stores a permutation of m objects */
/*       w(3+m:2+m+n) stores a permutation of n objects */
/*       w(3+m+n) = address in w of the initialization array */
/*                  for idd_random_transf */
/*       w(4+m+n:int(w(3+m+n))-1) stores the initialization array */
/*                                for dfft */
/*       w(int(w(3+m+n)):16*m+70) stores the initialization array */
/*                                for idd_random_transf */


/*       _N.B._: n is an output of the present routine; */
/*               this routine changes n. */




/*       Find the greatest integer less than or equal to m */
/*       which is a power of two. */

    /* Parameter adjustments */
    --w;

    /* Function Body */
    idd_poweroftwo__(m, &l, n);


/*       Store m and n in w. */

    w[1] = (doublereal) (*m);
    w[2] = (doublereal) (*n);


/*       Store random permutations of m and n objects in w. */

    id_randperm__(m, &w[3]);
    id_randperm__(n, &w[*m + 3]);


/*       Store the address within w of the idd_random_transf_init */
/*       initialization data. */

    ia = *m + 4 + *n + (*n << 1) + 15;
    w[*m + 3 + *n] = (doublereal) ia;


/*       Store the initialization data for dfft in w. */

    dffti_(n, &w[*m + 4 + *n]);


/*       Store the initialization data for idd_random_transf_init in w. */

    nsteps = 3;
    idd_random_transf_init__(&nsteps, m, &w[ia], &keep);


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
} /* idd_frmi__ */





/* Subroutine */ int idd_sfrmi__(integer *l, integer *m, integer *n, doublereal *w){
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int idd_sffti__(integer *, doublereal *, integer * , doublereal *);    static integer l2, ia, lw;
    extern /* Subroutine */ int id_randperm__(integer *, doublereal *);
    static integer keep;
    extern /* Subroutine */ int idd_copyints__(integer *, integer *, integer *), prinf_(char *, integer *, integer *, ftnlen), idd_pairsamps__(integer *, integer *, integer *, integer *, integer *, integer *);    static integer idummy, nsteps;
    extern /* Subroutine */ int idd_poweroftwo__(integer *, integer *, integer *), idd_random_transf_init__(integer *, integer *, doublereal *, integer *);

/*       initializes data for the routine idd_sfrm. */

/*       input: */
/*       l -- length of the transformed (output) vector */
/*       m -- length of the vector to be transformed */

/*       output: */
/*       n -- greatest integer expressible as a positive integer power */
/*            of 2 that is less than or equal to m */
/*       w -- initialization array to be used by routine idd_sfrm */


/*       glossary for the fully initialized w: */

/*       w(1) = m */
/*       w(2) = n */
/*       w(3) = l2 */
/*       w(4:3+m) stores a permutation of m objects */
/*       w(4+m:3+m+l) stores the indices of the l outputs which idd_sfft */
/*                    calculates */
/*       w(4+m+l:3+m+l+l2) stores the indices of the l2 pairs of outputs */
/*                         which idd_sfft calculates */
/*       w(4+m+l+l2) = address in w of the initialization array */
/*                     for idd_random_transf */
/*       w(5+m+l+l2:int(w(4+m+l+l2))-1) stores the initialization array */
/*                                      for idd_sfft */
/*       w(int(w(4+m+l+l2)):25*m+90) stores the initialization array */
/*                                   for idd_random_transf */


/*       _N.B._: n is an output of the present routine; */
/*               this routine changes n. */




/*       Find the greatest integer less than or equal to m */
/*       which is a power of two. */

    /* Parameter adjustments */
    --w;

    /* Function Body */
    idd_poweroftwo__(m, &idummy, n);


/*       Store m and n in w. */

    w[1] = (doublereal) (*m);
    w[2] = (doublereal) (*n);


/*       Store random permutations of m and n objects in w. */

    id_randperm__(m, &w[4]);
    id_randperm__(n, &w[*m + 4]);


/*       Find the pairs of integers covering the integers in */
/*       w(4+m : 3+m+(l+1)/2). */

    idd_pairsamps__(n, l, (integer*)&w[*m + 4], &l2, (integer*)&w[*m + 4 + (*
	    l << 1)], (integer*)&w[*m + 4 + *l * 3]);
    w[3] = (doublereal) l2;
    idd_copyints__(&l2, &w[*m + 4 + (*l << 1)], &w[*m + 4 + *l]);


/*       Store the address within w of the idd_random_transf_init */
/*       initialization data. */

    ia = *m + 5 + *l + l2 + (l2 << 2) + 30 + (*n << 3);
    w[*m + 4 + *l + l2] = (doublereal) ia;


/*       Store the initialization data for idd_sfft in w. */

    idd_sffti__(&l2, &w[*m + 4 + *l], n, &w[*m + 5 + *l + l2]);


/*       Store the initialization data for idd_random_transf_init in w. */

    nsteps = 3;
    idd_random_transf_init__(&nsteps, m, &w[ia], &keep);


/*       Calculate the total number of elements used in w. */

    lw = *m + 4 + *l + l2 + (l2 << 2) + 30 + (*n << 3) + nsteps * 3 * *m + (*
	    m << 1) + *m / 4 + 50;

    if (*m * 25 + 90 < lw) {
	prinf_("lw = *", &lw, &c__1, (ftnlen)6);
	i__1 = *m * 25 + 90;
	prinf_("25m+90 = *", &i__1, &c__1, (ftnlen)10);
	s_stop("", (ftnlen)0);
    }


    return 0;
} /* idd_sfrmi__ */





/* Subroutine */ int idd_copyints__(integer *n, integer *ia, integer *ib)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;


/*       copies ia into ib. */

/*       input: */
/*       n -- length of ia and ib */
/*       ia -- array to be copied */

/*       output: */
/*       ib -- copy of ia */



    /* Parameter adjustments */
    --ib;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	ib[k] = ia[k];
    }


/* k */
    return 0;
} /* idd_copyints__ */





/* Subroutine */ int idd_poweroftwo__(integer *m, integer *l, integer *n)
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
} /* idd_poweroftwo__ */

