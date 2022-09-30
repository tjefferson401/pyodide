/* sstats.f -- translated by f2c (version 20160102).
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


/* \SCCS Information: @(#) */
/* FILE: stats.F   SID: 2.1   DATE OF SID: 4/19/96   RELEASE: 2 */
/*     %---------------------------------------------% */
/*     | Initialize statistic and timing information | */
/*     | for symmetric Arnoldi code.                 | */
/*     %---------------------------------------------% */
/* Subroutine */ int sstats_(void)
{
    static integer nbx, nopx;
    static real trvec, tmvbx, tgetv0, tsaup2;
    static integer nitref;
    static real titref, tseigt, tsaupd, tsaitr, tsgets, tsapps;
    static integer nrorth;
    static real tsconv;
    static integer nrstrt;
    static real tmvopx;

/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */
/*     %-----------------------% */
/*     | Executable Statements | */
/*     %-----------------------% */
/*     %--------------------------------% */
/*     | See stat.doc for documentation | */
/*     %--------------------------------% */

/* \SCCS Information: @(#) */
/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */

/*      real       t0, t1, t2, t3, t4, t5 */
/*      save       t0, t1, t2, t3, t4, t5 */

/*      integer    nopx, nbx, nrorth, nitref, nrstrt */
/*      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, */
/*     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, */
/*     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv, */
/*     &           tmvopx, tmvbx, tgetv0, titref, trvec */
/*      common /timing/ */
/*     &           nopx, nbx, nrorth, nitref, nrstrt, */
/*     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, */
/*     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, */
/*     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv, */
/*     &           tmvopx, tmvbx, tgetv0, titref, trvec */
    nopx = 0;
    nbx = 0;
    nrorth = 0;
    nitref = 0;
    nrstrt = 0;
    tsaupd = 0.f;
    tsaup2 = 0.f;
    tsaitr = 0.f;
    tseigt = 0.f;
    tsgets = 0.f;
    tsapps = 0.f;
    tsconv = 0.f;
    titref = 0.f;
    tgetv0 = 0.f;
    trvec = 0.f;
/*     %----------------------------------------------------% */
/*     | User time including reverse communication overhead | */
/*     %----------------------------------------------------% */
    tmvopx = 0.f;
    tmvbx = 0.f;
    return 0;

/*     End of sstats */

} /* sstats_ */

