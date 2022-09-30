/* dcsrch.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int dcsrch_(doublereal *stp, doublereal *f, doublereal *g, 
	doublereal *ftol, doublereal *gtol, doublereal *xtol, char *task, 
	doublereal *stpmin, doublereal *stpmax, integer *isave, doublereal *
	dsave, ftnlen task_len)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal fm, gm, fx, fy, gx, gy, fxm, fym, gxm, gym, stx, sty;
    static integer stage;
    static doublereal finit, ginit, width, ftest, gtest, stmin, stmax, width1;
    static logical brackt;
    extern /* Subroutine */ int dcstep_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, doublereal *,
	     doublereal *);

/*     ********** */

/*     Subroutine dcsrch */

/*     This subroutine finds a step that satisfies a sufficient */
/*     decrease condition and a curvature condition. */

/*     Each call of the subroutine updates an interval with */
/*     endpoints stx and sty. The interval is initially chosen */
/*     so that it contains a minimizer of the modified function */

/*           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0). */

/*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
/*     interval is chosen so that it contains a minimizer of f. */

/*     The algorithm is designed to find a step that satisfies */
/*     the sufficient decrease condition */

/*           f(stp) <= f(0) + ftol*stp*f'(0), */

/*     and the curvature condition */

/*           abs(f'(stp)) <= gtol*abs(f'(0)). */

/*     If ftol is less than gtol and if, for example, the function */
/*     is bounded below, then there is always a step which satisfies */
/*     both conditions. */

/*     If no step can be found that satisfies both conditions, then */
/*     the algorithm stops with a warning. In this case stp only */
/*     satisfies the sufficient decrease condition. */

/*     A typical invocation of dcsrch has the following outline: */

/*     Evaluate the function at stp = 0.0d0; store in f. */
/*     Evaluate the gradient at stp = 0.0d0; store in g. */
/*     Choose a starting step stp. */

/*     task = 'START' */
/*  10 continue */
/*        call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax, */
/*    +               isave,dsave) */
/*        if (task .eq. 'FG') then */
/*           Evaluate the function and the gradient at stp */
/*           go to 10 */
/*           end if */

/*     NOTE: The user must not alter work arrays between calls. */

/*     The subroutine statement is */

/*       subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, */
/*                         task,isave,dsave) */
/*     where */

/*       stp is a double precision variable. */
/*         On entry stp is the current estimate of a satisfactory */
/*            step. On initial entry, a positive initial estimate */
/*            must be provided. */
/*         On exit stp is the current estimate of a satisfactory step */
/*            if task = 'FG'. If task = 'CONV' then stp satisfies */
/*            the sufficient decrease and curvature condition. */

/*       f is a double precision variable. */
/*         On initial entry f is the value of the function at 0. */
/*            On subsequent entries f is the value of the */
/*            function at stp. */
/*         On exit f is the value of the function at stp. */

/*       g is a double precision variable. */
/*         On initial entry g is the derivative of the function at 0. */
/*            On subsequent entries g is the derivative of the */
/*            function at stp. */
/*         On exit g is the derivative of the function at stp. */

/*       ftol is a double precision variable. */
/*         On entry ftol specifies a nonnegative tolerance for the */
/*            sufficient decrease condition. */
/*         On exit ftol is unchanged. */

/*       gtol is a double precision variable. */
/*         On entry gtol specifies a nonnegative tolerance for the */
/*            curvature condition. */
/*         On exit gtol is unchanged. */

/*       xtol is a double precision variable. */
/*         On entry xtol specifies a nonnegative relative tolerance */
/*            for an acceptable step. The subroutine exits with a */
/*            warning if the relative difference between sty and stx */
/*            is less than xtol. */
/*         On exit xtol is unchanged. */

/*       task is a character variable of length at least 60. */
/*         On initial entry task must be set to 'START'. */
/*         On exit task indicates the required action: */

/*            If task(1:2) = 'FG' then evaluate the function and */
/*            derivative at stp and call dcsrch again. */

/*            If task(1:4) = 'CONV' then the search is successful. */

/*            If task(1:4) = 'WARN' then the subroutine is not able */
/*            to satisfy the convergence conditions. The exit value of */
/*            stp contains the best point found during the search. */

/*            If task(1:5) = 'ERROR' then there is an error in the */
/*            input arguments. */

/*         On exit with convergence, a warning or an error, the */
/*            variable task contains additional information. */

/*       stpmin is a double precision variable. */
/*         On entry stpmin is a nonnegative lower bound for the step. */
/*         On exit stpmin is unchanged. */

/*       stpmax is a double precision variable. */
/*         On entry stpmax is a nonnegative upper bound for the step. */
/*         On exit stpmax is unchanged. */

/*       isave is an integer work array of dimension 2. */

/*       dsave is a double precision work array of dimension 13. */

/*     Subprograms called */

/*       MINPACK-2 ... dcstep */

/*     MINPACK-1 Project. June 1983. */
/*     Argonne National Laboratory. */
/*     Jorge J. More' and David J. Thuente. */

/*     MINPACK-2 Project. November 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick, Richard G. Carter, and Jorge J. More'. */

/*     ********** */
/*     Initialization block. */
    /* Parameter adjustments */
    --dsave;
    --isave;

    /* Function Body */
    if (s_cmp(task, "START", (ftnlen)5, (ftnlen)5) == 0) {
/*        Check the input arguments for errors. */
	if (*stp < *stpmin) {
	    s_copy(task, "ERROR: STP .LT. STPMIN", task_len, (ftnlen)22);
	}
	if (*stp > *stpmax) {
	    s_copy(task, "ERROR: STP .GT. STPMAX", task_len, (ftnlen)22);
	}
	if (*g >= 0.) {
	    s_copy(task, "ERROR: INITIAL G .GE. ZERO", task_len, (ftnlen)26);
	}
	if (*ftol < 0.) {
	    s_copy(task, "ERROR: FTOL .LT. ZERO", task_len, (ftnlen)21);
	}
	if (*gtol < 0.) {
	    s_copy(task, "ERROR: GTOL .LT. ZERO", task_len, (ftnlen)21);
	}
	if (*xtol < 0.) {
	    s_copy(task, "ERROR: XTOL .LT. ZERO", task_len, (ftnlen)21);
	}
	if (*stpmin < 0.) {
	    s_copy(task, "ERROR: STPMIN .LT. ZERO", task_len, (ftnlen)23);
	}
	if (*stpmax < *stpmin) {
	    s_copy(task, "ERROR: STPMAX .LT. STPMIN", task_len, (ftnlen)25);
	}
/*        Exit if there are errors on input. */
	if (s_cmp(task, "ERROR", (ftnlen)5, (ftnlen)5) == 0) {
	    return 0;
	}
/*        Initialize local variables. */
	brackt = FALSE_;
	stage = 1;
	finit = *f;
	ginit = *g;
	gtest = *ftol * ginit;
	width = *stpmax - *stpmin;
	width1 = width / .5;
/*        The variables stx, fx, gx contain the values of the step, */
/*        function, and derivative at the best step. */
/*        The variables sty, fy, gy contain the value of the step, */
/*        function, and derivative at sty. */
/*        The variables stp, f, g contain the values of the step, */
/*        function, and derivative at stp. */
	stx = 0.;
	fx = finit;
	gx = ginit;
	sty = 0.;
	fy = finit;
	gy = ginit;
	stmin = 0.;
	stmax = *stp + *stp * 4.;
	s_copy(task, "FG", task_len, (ftnlen)2);
	goto L10;
    } else {
/*        Restore local variables. */
	if (isave[1] == 1) {
	    brackt = TRUE_;
	} else {
	    brackt = FALSE_;
	}
	stage = isave[2];
	ginit = dsave[1];
	gtest = dsave[2];
	gx = dsave[3];
	gy = dsave[4];
	finit = dsave[5];
	fx = dsave[6];
	fy = dsave[7];
	stx = dsave[8];
	sty = dsave[9];
	stmin = dsave[10];
	stmax = dsave[11];
	width = dsave[12];
	width1 = dsave[13];
    }
/*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
/*     algorithm enters the second stage. */
    ftest = finit + *stp * gtest;
    if (stage == 1 && *f <= ftest && *g >= 0.) {
	stage = 2;
    }
/*     Test for warnings. */
    if (brackt && (*stp <= stmin || *stp >= stmax)) {
	s_copy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS", task_len, (
		ftnlen)41);
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
	s_copy(task, "WARNING: XTOL TEST SATISFIED", task_len, (ftnlen)28);
    }
    if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
	s_copy(task, "WARNING: STP = STPMAX", task_len, (ftnlen)21);
    }
    if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
	s_copy(task, "WARNING: STP = STPMIN", task_len, (ftnlen)21);
    }
/*     Test for convergence. */
    if (*f <= ftest && abs(*g) <= *gtol * (-ginit)) {
	s_copy(task, "CONVERGENCE", task_len, (ftnlen)11);
    }
/*     Test for termination. */
    if (s_cmp(task, "WARN", (ftnlen)4, (ftnlen)4) == 0 || s_cmp(task, "CONV", 
	    (ftnlen)4, (ftnlen)4) == 0) {
	goto L10;
    }
/*     A modified function is used to predict the step during the */
/*     first stage if a lower function value has been obtained but */
/*     the decrease is not sufficient. */
    if (stage == 1 && *f <= fx && *f > ftest) {
/*        Define the modified function and derivative values. */
	fm = *f - *stp * gtest;
	fxm = fx - stx * gtest;
	fym = fy - sty * gtest;
	gm = *g - gtest;
	gxm = gx - gtest;
	gym = gy - gtest;
/*        Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
		stmin, &stmax);
/*        Reset the function and derivative values for f. */
	fx = fxm + stx * gtest;
	fy = fym + sty * gtest;
	gx = gxm + gtest;
	gy = gym + gtest;
    } else {
/*       Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
		stmax);
    }
/*     Decide if a bisection step is needed. */
    if (brackt) {
	if ((d__1 = sty - stx, abs(d__1)) >= width1 * .66) {
	    *stp = stx + (sty - stx) * .5;
	}
	width1 = width;
	width = (d__1 = sty - stx, abs(d__1));
    }
/*     Set the minimum and maximum steps allowed for stp. */
    if (brackt) {
	stmin = min(stx,sty);
	stmax = max(stx,sty);
    } else {
	stmin = *stp + (*stp - stx) * 1.1;
	stmax = *stp + (*stp - stx) * 4.;
    }
/*     Force the step to be within the bounds stpmax and stpmin. */
    *stp = max(*stp,*stpmin);
    *stp = min(*stp,*stpmax);
/*     If further progress is not possible, let stp be the best */
/*     point obtained during the search. */
    if (brackt && (*stp <= stmin || *stp >= stmax) || brackt && stmax - stmin 
	    <= *xtol * stmax) {
	*stp = stx;
    }
/*     Obtain another function and derivative. */
    s_copy(task, "FG", task_len, (ftnlen)2);
L10:
/*     Save local variables. */
    if (brackt) {
	isave[1] = 1;
    } else {
	isave[1] = 0;
    }
    isave[2] = stage;
    dsave[1] = ginit;
    dsave[2] = gtest;
    dsave[3] = gx;
    dsave[4] = gy;
    dsave[5] = finit;
    dsave[6] = fx;
    dsave[7] = fy;
    dsave[8] = stx;
    dsave[9] = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;
    return 0;
} /* dcsrch_ */

