/* idd_sfft.f -- translated by f2c (version 20160102).
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


/*       routine idd_sffti initializes routine idd_sfft. */

/*       routine idd_sfft rapidly computes a subset of the entries */
/*       of the DFT of a vector, composed with permutation matrices */
/*       both on input and on output. */

/*       routine idd_ldiv finds the greatest integer less than or equal */
/*       to a specified integer, that is divisible by another (larger) */
/*       specified integer. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */




/* Subroutine */ int idd_ldiv__(integer *l, integer *n, integer *m)
{

/*       finds the greatest integer less than or equal to l */
/*       that divides n. */

/*       input: */
/*       l -- integer at least as great as m */
/*       n -- integer divisible by m */

/*       output: */
/*       m -- greatest integer less than or equal to l that divides n */



    *m = *l;

L1000:
    if (*m * (*n / *m) == *n) {
	goto L2000;
    }

    --(*m);
    goto L1000;

L2000:


    return 0;
} /* idd_ldiv__ */





/* Subroutine */ int idd_sffti__(integer *l, integer *ind, integer *n, doublecomplex *wsave){
    extern /* Subroutine */ int idd_sffti1__(integer *, integer *, doublereal *), idd_sffti2__(integer *, integer *, integer *, doublecomplex *);

/*       initializes wsave for using routine idd_sfft. */

/*       input: */
/*       l -- number of pairs of entries in the output of idd_sfft */
/*            to compute */
/*       ind -- indices of the pairs of entries in the output */
/*              of idd_sfft to compute; the indices must be chosen */
/*              in the range from 1 to n/2 */
/*       n -- length of the vector to be transformed */

/*       output: */
/*       wsave -- array needed by routine idd_sfft for processing */
/*                (the present routine does not use the last n elements */
/*                 of wsave, but routine idd_sfft does) */



    /* Parameter adjustments */
    --ind;
    --wsave;

    /* Function Body */
    if (*l == 1) {
	idd_sffti1__(&ind[1], n, &wsave[1]);
    }
    if (*l > 1) {
	idd_sffti2__(l, &ind[1], n, &wsave[1]);
    }


    return 0;
} /* idd_sffti__ */





/* Subroutine */ int idd_sffti1__(integer *ind, integer *n, doublereal *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), cos(doublereal), sin(
	    doublereal);

    /* Local variables */
    static integer k;
    static doublereal r1, fact, twopi;


/*       routine idd_sffti serves as a wrapper around */
/*       the present routine; please see routine idd_sffti */
/*       for documentation. */


    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    r1 = 1.;
    twopi = atan(r1) * 8;


    fact = 1 / sqrt(r1 * *n);


    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	wsave[k] = cos(twopi * (k - 1) * *ind / (r1 * *n)) * fact;
    }

/* k */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	wsave[*n + k] = -sin(twopi * (k - 1) * *ind / (r1 * *n)) * fact;
    }


/* k */
    return 0;
} /* idd_sffti1__ */





/* Subroutine */ int idd_sffti2__(integer *l, integer *ind, integer *n, doublecomplex *wsave){
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9, z__10,
	     z__11, z__12;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal);
    void z_exp(doublecomplex *, doublecomplex *);

    /* Local variables */
    extern /* Subroutine */ int idd_ldiv__(integer *, integer *, integer *);
    static integer i__, j, k, m;
    static doublereal r1;
    static doublecomplex ci;
    static integer ii;
    static doublereal fact;
    extern /* Subroutine */ int dffti_(integer *, doublecomplex *);
    static integer imodm, idivm;
    static doublereal twopi;
    static integer nblock;
    static doublecomplex twopii;


/*       routine idd_sffti serves as a wrapper around */
/*       the present routine; please see routine idd_sffti */
/*       for documentation. */


    /* Parameter adjustments */
    --ind;
    --wsave;

    /* Function Body */
    ci.r = 0.f, ci.i = 1.f;
    r1 = 1.;
    twopi = atan(r1) * 8;
    z__1.r = twopi * ci.r, z__1.i = twopi * ci.i;
    twopii.r = z__1.r, twopii.i = z__1.i;


/*       Determine the block lengths for the FFTs. */

    idd_ldiv__(l, n, &nblock);
    m = *n / nblock;


/*       Initialize wsave for using routine dfftf. */

    dffti_(&nblock, &wsave[1]);


/*       Calculate the coefficients in the linear combinations */
/*       needed for the direct portion of the calculation. */

    fact = 1 / sqrt(r1 * *n);

    ii = (*l << 1) + 15;

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {


	i__ = ind[j];


	if (i__ <= *n / 2 - m / 2) {

	    idivm = (i__ - 1) / m;
	    imodm = i__ - 1 - m * idivm;

	    i__2 = m;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = ii + m * (j - 1) + k;
		z__7.r = -twopii.r, z__7.i = -twopii.i;
		i__4 = k - 1;
		d__1 = (doublereal) i__4;
		z__6.r = d__1 * z__7.r, z__6.i = d__1 * z__7.i;
		d__2 = (doublereal) imodm;
		z__5.r = d__2 * z__6.r, z__5.i = d__2 * z__6.i;
		d__3 = r1 * m;
		z__4.r = z__5.r / d__3, z__4.i = z__5.i / d__3;
		z_exp(&z__3, &z__4);
		z__12.r = -twopii.r, z__12.i = -twopii.i;
		i__5 = k - 1;
		d__4 = (doublereal) i__5;
		z__11.r = d__4 * z__12.r, z__11.i = d__4 * z__12.i;
		i__6 = idivm + 1;
		d__5 = (doublereal) i__6;
		z__10.r = d__5 * z__11.r, z__10.i = d__5 * z__11.i;
		d__6 = r1 * *n;
		z__9.r = z__10.r / d__6, z__9.i = z__10.i / d__6;
		z_exp(&z__8, &z__9);
		z__2.r = z__3.r * z__8.r - z__3.i * z__8.i, z__2.i = z__3.r * 
			z__8.i + z__3.i * z__8.r;
		z__1.r = fact * z__2.r, z__1.i = fact * z__2.i;
		wsave[i__3].r = z__1.r, wsave[i__3].i = z__1.i;
	    }

/* k */
	}


/* i .le. n/2-m/2 */
	if (i__ > *n / 2 - m / 2) {

	    idivm = i__ / (m / 2);
	    imodm = i__ - m / 2 * idivm;

	    i__2 = m;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = ii + m * (j - 1) + k;
		z__6.r = -twopii.r, z__6.i = -twopii.i;
		i__4 = k - 1;
		d__1 = (doublereal) i__4;
		z__5.r = d__1 * z__6.r, z__5.i = d__1 * z__6.i;
		d__2 = (doublereal) imodm;
		z__4.r = d__2 * z__5.r, z__4.i = d__2 * z__5.i;
		d__3 = r1 * m;
		z__3.r = z__4.r / d__3, z__3.i = z__4.i / d__3;
		z_exp(&z__2, &z__3);
		z__1.r = fact * z__2.r, z__1.i = fact * z__2.i;
		wsave[i__3].r = z__1.r, wsave[i__3].i = z__1.i;
	    }

/* k */
	}


/* i .gt. n/2-m/2 */
    }


/* j */
    return 0;
} /* idd_sffti2__ */





/* Subroutine */ int idd_sfft__(integer *l, integer *ind, integer *n, doublecomplex *wsave, doublereal *v){
    extern /* Subroutine */ int idd_sfft1__(integer *, integer *, doublereal *, doublereal *), idd_sfft2__(integer *, integer *, integer *, doublereal *, doublecomplex *);

/*       computes a subset of the entries of the DFT of v, */
/*       composed with permutation matrices both on input and on output, */
/*       via a two-stage procedure (debugging code routine dfftf2 above */
/*       is supposed to calculate the full vector from which idd_sfft */
/*       returns a subset of the entries, when dfftf2 has */
/*       the same parameter nblock as in the present routine). */

/*       input: */
/*       l -- number of pairs of entries in the output to compute */
/*       ind -- indices of the pairs of entries in the output */
/*              to compute; the indices must be chosen */
/*              in the range from 1 to n/2 */
/*       n -- length of v; n must be a positive integer power of 2 */
/*       v -- vector to be transformed */
/*       wsave -- processing array initialized by routine idd_sffti */

/*       output: */
/*       v -- pairs of entries indexed by ind are given */
/*            their appropriately transformed values */

/*       _N.B._: n must be a positive integer power of 2. */

/*       references: */
/*       Sorensen and Burrus, "Efficient computation of the DFT with */
/*            only a subset of input or output points," */
/*            IEEE Transactions on Signal Processing, 41 (3): 1184-1200, */
/*            1993. */
/*       Woolfe, Liberty, Rokhlin, Tygert, "A fast randomized algorithm */
/*            for the approximation of matrices," Applied and */
/*            Computational Harmonic Analysis, 25 (3): 335-366, 2008; */
/*            Section 3.3. */



    /* Parameter adjustments */
    --ind;
    --v;
    --wsave;

    /* Function Body */
    if (*l == 1) {
	idd_sfft1__(&ind[1], n, &v[1], &wsave[1]);
    }
    if (*l > 1) {
	idd_sfft2__(l, &ind[1], n, &v[1], &wsave[1]);
    }


    return 0;
} /* idd_sfft__ */





/* Subroutine */ int idd_sfft1__(integer *ind, integer *n, doublereal *v, doublereal *wsave){
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal r1, fact, sumi, sumr, twopi;


/*       routine idd_sfft serves as a wrapper around */
/*       the present routine; please see routine idd_sfft */
/*       for documentation. */


    /* Parameter adjustments */
    --wsave;
    --v;

    /* Function Body */
    r1 = 1.;
    twopi = atan(r1) * 8;


    if (*ind < *n / 2) {


	sumr = 0.;

	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    sumr += wsave[k] * v[k];
	}


/* k */
	sumi = 0.;

	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    sumi += wsave[*n + k] * v[k];
	}


/* k */
    }


/* ind .lt. n/2 */
    if (*ind == *n / 2) {


	fact = 1 / sqrt(r1 * *n);


	sumr = 0.;

	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    sumr += v[k];
	}

/* k */
	sumr *= fact;


	sumi = 0.;

	i__1 = *n / 2;
	for (k = 1; k <= i__1; ++k) {
	    sumi += v[(k << 1) - 1];
	    sumi -= v[k * 2];
	}

/* k */
	sumi *= fact;


    }


/* ind .eq. n/2 */
    v[(*ind << 1) - 1] = sumr;
    v[*ind * 2] = sumi;


    return 0;
} /* idd_sfft1__ */





/* Subroutine */ int idd_sfft2__(integer *l, integer *ind, integer *n, doublereal *v, doublecomplex *wsave){
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int idd_ldiv__(integer *, integer *, integer *);
    static integer i__, j, k, m;
    static doublereal r1;
    static doublecomplex ci;
    static integer ii, iii;
    static doublecomplex sum;
    static doublereal fact, rsum;
    extern /* Subroutine */ int dfftf_(integer *, doublereal *, doublecomplex *);    static integer imodm, idivm;
    static doublereal twopi;
    static integer nblock;


/*       routine idd_sfft serves as a wrapper around */
/*       the present routine; please see routine idd_sfft */
/*       for documentation. */


    /* Parameter adjustments */
    --ind;
    --wsave;
    --v;

    /* Function Body */
    ci.r = 0.f, ci.i = 1.f;
    r1 = 1.;
    twopi = atan(r1) * 8;


/*       Determine the block lengths for the FFTs. */

    idd_ldiv__(l, n, &nblock);


    m = *n / nblock;


/*       FFT each block of length nblock of v. */

    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	dfftf_(&nblock, &v[nblock * (k - 1) + 1], &wsave[1]);
    }


/*       Transpose v to obtain wsave(2*l+15+2*n+1 : 2*l+15+3*n). */

/* k */
    iii = (*l << 1) + 15 + (*n << 1);

    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nblock / 2 - 1;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = iii + m * (j - 1) + k;
	    i__4 = nblock * (k - 1) + (j << 1);
	    i__5 = nblock * (k - 1) + (j << 1) + 1;
	    z__2.r = v[i__5] * ci.r, z__2.i = v[i__5] * ci.i;
	    z__1.r = v[i__4] + z__2.r, z__1.i = z__2.i;
	    wsave[i__3].r = z__1.r, wsave[i__3].i = z__1.i;
	}
/* j */
    }

/*       Handle the purely real frequency components separately. */

/* k */
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	i__2 = iii + m * (nblock / 2 - 1) + k;
	i__3 = nblock * (k - 1) + nblock;
	wsave[i__2].r = v[i__3], wsave[i__2].i = 0.;
	i__2 = iii + m * (nblock / 2) + k;
	i__3 = nblock * (k - 1) + 1;
	wsave[i__2].r = v[i__3], wsave[i__2].i = 0.;
    }


/*       Directly calculate the desired entries of v. */

/* k */
    ii = (*l << 1) + 15;

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {


	i__ = ind[j];


	if (i__ <= *n / 2 - m / 2) {

	    idivm = (i__ - 1) / m;
	    imodm = i__ - 1 - m * idivm;

	    sum.r = 0., sum.i = 0.;

	    i__2 = m;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = iii + m * idivm + k;
		i__4 = ii + m * (j - 1) + k;
		z__2.r = wsave[i__3].r * wsave[i__4].r - wsave[i__3].i * 
			wsave[i__4].i, z__2.i = wsave[i__3].r * wsave[i__4].i 
			+ wsave[i__3].i * wsave[i__4].r;
		z__1.r = sum.r + z__2.r, z__1.i = sum.i + z__2.i;
		sum.r = z__1.r, sum.i = z__1.i;
	    }

/* k */
	    i__2 = (i__ << 1) - 1;
	    v[i__2] = sum.r;
	    i__2 = i__ << 1;
	    z__2.r = -ci.r, z__2.i = -ci.i;
	    z__1.r = z__2.r * sum.r - z__2.i * sum.i, z__1.i = z__2.r * sum.i 
		    + z__2.i * sum.r;
	    v[i__2] = z__1.r;

	}


/* i .le. n/2-m/2 */
	if (i__ > *n / 2 - m / 2) {

	    if (i__ < *n / 2) {

		idivm = i__ / (m / 2);
		imodm = i__ - m / 2 * idivm;

		sum.r = 0., sum.i = 0.;

		i__2 = m;
		for (k = 1; k <= i__2; ++k) {
		    i__3 = iii + m * (nblock / 2) + k;
		    i__4 = ii + m * (j - 1) + k;
		    z__2.r = wsave[i__3].r * wsave[i__4].r - wsave[i__3].i * 
			    wsave[i__4].i, z__2.i = wsave[i__3].r * wsave[
			    i__4].i + wsave[i__3].i * wsave[i__4].r;
		    z__1.r = sum.r + z__2.r, z__1.i = sum.i + z__2.i;
		    sum.r = z__1.r, sum.i = z__1.i;
		}

/* k */
		i__2 = (i__ << 1) - 1;
		v[i__2] = sum.r;
		i__2 = i__ << 1;
		z__2.r = -ci.r, z__2.i = -ci.i;
		z__1.r = z__2.r * sum.r - z__2.i * sum.i, z__1.i = z__2.r * 
			sum.i + z__2.i * sum.r;
		v[i__2] = z__1.r;

	    }

	    if (i__ == *n / 2) {

		fact = 1 / sqrt(r1 * *n);


		rsum = 0.;

		i__2 = m;
		for (k = 1; k <= i__2; ++k) {
		    i__3 = iii + m * (nblock / 2) + k;
		    z__1.r = rsum + wsave[i__3].r, z__1.i = wsave[i__3].i;
		    rsum = z__1.r;
		}

/* k */
		v[*n - 1] = rsum * fact;


		rsum = 0.;

		i__2 = m / 2;
		for (k = 1; k <= i__2; ++k) {
		    i__3 = iii + m * (nblock / 2) + (k << 1) - 1;
		    z__1.r = rsum + wsave[i__3].r, z__1.i = wsave[i__3].i;
		    rsum = z__1.r;
		    i__3 = iii + m * (nblock / 2) + (k << 1);
		    z__1.r = rsum - wsave[i__3].r, z__1.i = -wsave[i__3].i;
		    rsum = z__1.r;
		}

/* k */
		v[*n] = rsum * fact;

	    }

	}


/* i .gt. n/2-m/2 */
    }


/* j */
    return 0;
} /* idd_sfft2__ */

