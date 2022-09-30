/* _arpack-f2pywrappers.f -- translated by f2c (version 20160102).
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

/* Common Block Declarations */

struct {
    integer logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps, 
	    msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, 
	    mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

struct {
    integer nopx, nbx, nrorth, nitref, nrstrt;
    real tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd, 
	    tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2, 
	    tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0, 
	    titref, trvec;
} timing_;

#define timing_1 timing_

/*     -*- fortran -*- */
/*     This file is autogenerated with f2py (version:1.21.4) */
/*     It contains Fortran 77 wrappers to fortran functions. */
/* Subroutine */ int f2pyinitdebug_(S_fp setupfunc)
{
    (*setupfunc)(&debug_1.logfil, &debug_1.ndigit, &debug_1.mgetv0, &
	    debug_1.msaupd, &debug_1.msaup2, &debug_1.msaitr, &debug_1.mseigt,
	     &debug_1.msapps, &debug_1.msgets, &debug_1.mseupd, &
	    debug_1.mnaupd, &debug_1.mnaup2, &debug_1.mnaitr, &debug_1.mneigh,
	     &debug_1.mnapps, &debug_1.mngets, &debug_1.mneupd, &
	    debug_1.mcaupd, &debug_1.mcaup2, &debug_1.mcaitr, &debug_1.mceigh,
	     &debug_1.mcapps, &debug_1.mcgets, &debug_1.mceupd);
    return 0;
} /* f2pyinitdebug_ */

/* Subroutine */ int f2pyinittiming_(S_fp setupfunc)
{
    (*setupfunc)(&timing_1.nopx, &timing_1.nbx, &timing_1.nrorth, &
	    timing_1.nitref, &timing_1.nrstrt, &timing_1.tsaupd, &
	    timing_1.tsaup2, &timing_1.tsaitr, &timing_1.tseigt, &
	    timing_1.tsgets, &timing_1.tsapps, &timing_1.tsconv, &
	    timing_1.tnaupd, &timing_1.tnaup2, &timing_1.tnaitr, &
	    timing_1.tneigh, &timing_1.tngets, &timing_1.tnapps, &
	    timing_1.tnconv, &timing_1.tcaupd, &timing_1.tcaup2, &
	    timing_1.tcaitr, &timing_1.tceigh, &timing_1.tcgets, &
	    timing_1.tcapps, &timing_1.tcconv, &timing_1.tmvopx, &
	    timing_1.tmvbx, &timing_1.tgetv0, &timing_1.titref, &
	    timing_1.trvec);
    return 0;
} /* f2pyinittiming_ */

