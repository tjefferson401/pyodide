/* mvn-f2pywrappers.f -- translated by f2c (version 20160102).
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
    integer ivls;
} dkblck_;

#define dkblck_1 dkblck_

/*     -*- fortran -*- */
/*     This file is autogenerated with f2py (version:1.21.4) */
/*     It contains Fortran 77 wrappers to fortran functions. */
/* Subroutine */ int f2pyinitdkblck_(S_fp setupfunc)
{
    (*setupfunc)(&dkblck_1.ivls);
    return 0;
} /* f2pyinitdkblck_ */

