/* This file was generated by _generate_pyx.py. */
/* Do not edit this file directly. */

#ifndef SCIPY_LINALG_BLAS_FORTRAN_WRAPPERS_H
#define SCIPY_LINALG_BLAS_FORTRAN_WRAPPERS_H
#include "fortran_defs.h"
#include "numpy/arrayobject.h"

#ifdef __cplusplus
extern "C" {
#endif

int F_FUNC(caxpywrp, CAXPYWRP)(int *ret, int *n, npy_complex64 *ca, npy_complex64 *cx, int *incx, npy_complex64 *cy, int *incy);
int F_FUNC(ccopywrp, CCOPYWRP)(int *ret, int *n, npy_complex64 *cx, int *incx, npy_complex64 *cy, int *incy);
int F_FUNC(cdotcwrp, CDOTCWRP)(npy_complex64 *ret, int *n, npy_complex64 *cx, int *incx, npy_complex64 *cy, int *incy);
int F_FUNC(cdotuwrp, CDOTUWRP)(npy_complex64 *ret, int *n, npy_complex64 *cx, int *incx, npy_complex64 *cy, int *incy);
int F_FUNC(cgbmvwrp, CGBMVWRP)(int *ret, char *trans, int *m, int *n, int *kl, int *ku, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx, npy_complex64 *beta, npy_complex64 *y, int *incy);
int F_FUNC(cgemmwrp, CGEMMWRP)(int *ret, char *transa, char *transb, int *m, int *n, int *k, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *b, int *ldb, npy_complex64 *beta, npy_complex64 *c, int *ldc);
int F_FUNC(cgemvwrp, CGEMVWRP)(int *ret, char *trans, int *m, int *n, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx, npy_complex64 *beta, npy_complex64 *y, int *incy);
int F_FUNC(cgercwrp, CGERCWRP)(int *ret, int *m, int *n, npy_complex64 *alpha, npy_complex64 *x, int *incx, npy_complex64 *y, int *incy, npy_complex64 *a, int *lda);
int F_FUNC(cgeruwrp, CGERUWRP)(int *ret, int *m, int *n, npy_complex64 *alpha, npy_complex64 *x, int *incx, npy_complex64 *y, int *incy, npy_complex64 *a, int *lda);
int F_FUNC(chbmvwrp, CHBMVWRP)(int *ret, char *uplo, int *n, int *k, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx, npy_complex64 *beta, npy_complex64 *y, int *incy);
int F_FUNC(chemmwrp, CHEMMWRP)(int *ret, char *side, char *uplo, int *m, int *n, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *b, int *ldb, npy_complex64 *beta, npy_complex64 *c, int *ldc);
int F_FUNC(chemvwrp, CHEMVWRP)(int *ret, char *uplo, int *n, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx, npy_complex64 *beta, npy_complex64 *y, int *incy);
int F_FUNC(cherwrp, CHERWRP)(int *ret, char *uplo, int *n, float *alpha, npy_complex64 *x, int *incx, npy_complex64 *a, int *lda);
int F_FUNC(cher2wrp, CHER2WRP)(int *ret, char *uplo, int *n, npy_complex64 *alpha, npy_complex64 *x, int *incx, npy_complex64 *y, int *incy, npy_complex64 *a, int *lda);
int F_FUNC(cher2kwrp, CHER2KWRP)(int *ret, char *uplo, char *trans, int *n, int *k, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *b, int *ldb, float *beta, npy_complex64 *c, int *ldc);
int F_FUNC(cherkwrp, CHERKWRP)(int *ret, char *uplo, char *trans, int *n, int *k, float *alpha, npy_complex64 *a, int *lda, float *beta, npy_complex64 *c, int *ldc);
int F_FUNC(chpmvwrp, CHPMVWRP)(int *ret, char *uplo, int *n, npy_complex64 *alpha, npy_complex64 *ap, npy_complex64 *x, int *incx, npy_complex64 *beta, npy_complex64 *y, int *incy);
int F_FUNC(chprwrp, CHPRWRP)(int *ret, char *uplo, int *n, float *alpha, npy_complex64 *x, int *incx, npy_complex64 *ap);
int F_FUNC(chpr2wrp, CHPR2WRP)(int *ret, char *uplo, int *n, npy_complex64 *alpha, npy_complex64 *x, int *incx, npy_complex64 *y, int *incy, npy_complex64 *ap);
int F_FUNC(crotgwrp, CROTGWRP)(int *ret, npy_complex64 *ca, npy_complex64 *cb, float *c, npy_complex64 *s);
int F_FUNC(cscalwrp, CSCALWRP)(int *ret, int *n, npy_complex64 *ca, npy_complex64 *cx, int *incx);
int F_FUNC(csrotwrp, CSROTWRP)(int *ret, int *n, npy_complex64 *cx, int *incx, npy_complex64 *cy, int *incy, float *c, float *s);
int F_FUNC(csscalwrp, CSSCALWRP)(int *ret, int *n, float *sa, npy_complex64 *cx, int *incx);
int F_FUNC(cswapwrp, CSWAPWRP)(int *ret, int *n, npy_complex64 *cx, int *incx, npy_complex64 *cy, int *incy);
int F_FUNC(csymmwrp, CSYMMWRP)(int *ret, char *side, char *uplo, int *m, int *n, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *b, int *ldb, npy_complex64 *beta, npy_complex64 *c, int *ldc);
int F_FUNC(csyr2kwrp, CSYR2KWRP)(int *ret, char *uplo, char *trans, int *n, int *k, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *b, int *ldb, npy_complex64 *beta, npy_complex64 *c, int *ldc);
int F_FUNC(csyrkwrp, CSYRKWRP)(int *ret, char *uplo, char *trans, int *n, int *k, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *beta, npy_complex64 *c, int *ldc);
int F_FUNC(ctbmvwrp, CTBMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx);
int F_FUNC(ctbsvwrp, CTBSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx);
int F_FUNC(ctpmvwrp, CTPMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex64 *ap, npy_complex64 *x, int *incx);
int F_FUNC(ctpsvwrp, CTPSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex64 *ap, npy_complex64 *x, int *incx);
int F_FUNC(ctrmmwrp, CTRMMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *b, int *ldb);
int F_FUNC(ctrmvwrp, CTRMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx);
int F_FUNC(ctrsmwrp, CTRSMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, npy_complex64 *alpha, npy_complex64 *a, int *lda, npy_complex64 *b, int *ldb);
int F_FUNC(ctrsvwrp, CTRSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex64 *a, int *lda, npy_complex64 *x, int *incx);
int F_FUNC(dasumwrp, DASUMWRP)(double *ret, int *n, double *dx, int *incx);
int F_FUNC(daxpywrp, DAXPYWRP)(int *ret, int *n, double *da, double *dx, int *incx, double *dy, int *incy);
int F_FUNC(dcabs1wrp, DCABS1WRP)(double *ret, npy_complex128 *z);
int F_FUNC(dcopywrp, DCOPYWRP)(int *ret, int *n, double *dx, int *incx, double *dy, int *incy);
int F_FUNC(ddotwrp, DDOTWRP)(double *ret, int *n, double *dx, int *incx, double *dy, int *incy);
int F_FUNC(dgbmvwrp, DGBMVWRP)(int *ret, char *trans, int *m, int *n, int *kl, int *ku, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
int F_FUNC(dgemmwrp, DGEMMWRP)(int *ret, char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
int F_FUNC(dgemvwrp, DGEMVWRP)(int *ret, char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
int F_FUNC(dgerwrp, DGERWRP)(int *ret, int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);
int F_FUNC(dnrm2wrp, DNRM2WRP)(double *ret, int *n, double *x, int *incx);
int F_FUNC(drotwrp, DROTWRP)(int *ret, int *n, double *dx, int *incx, double *dy, int *incy, double *c, double *s);
int F_FUNC(drotgwrp, DROTGWRP)(int *ret, double *da, double *db, double *c, double *s);
int F_FUNC(drotmwrp, DROTMWRP)(int *ret, int *n, double *dx, int *incx, double *dy, int *incy, double *dparam);
int F_FUNC(drotmgwrp, DROTMGWRP)(int *ret, double *dd1, double *dd2, double *dx1, double *dy1, double *dparam);
int F_FUNC(dsbmvwrp, DSBMVWRP)(int *ret, char *uplo, int *n, int *k, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
int F_FUNC(dscalwrp, DSCALWRP)(int *ret, int *n, double *da, double *dx, int *incx);
int F_FUNC(dsdotwrp, DSDOTWRP)(double *ret, int *n, float *sx, int *incx, float *sy, int *incy);
int F_FUNC(dspmvwrp, DSPMVWRP)(int *ret, char *uplo, int *n, double *alpha, double *ap, double *x, int *incx, double *beta, double *y, int *incy);
int F_FUNC(dsprwrp, DSPRWRP)(int *ret, char *uplo, int *n, double *alpha, double *x, int *incx, double *ap);
int F_FUNC(dspr2wrp, DSPR2WRP)(int *ret, char *uplo, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *ap);
int F_FUNC(dswapwrp, DSWAPWRP)(int *ret, int *n, double *dx, int *incx, double *dy, int *incy);
int F_FUNC(dsymmwrp, DSYMMWRP)(int *ret, char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
int F_FUNC(dsymvwrp, DSYMVWRP)(int *ret, char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
int F_FUNC(dsyrwrp, DSYRWRP)(int *ret, char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda);
int F_FUNC(dsyr2wrp, DSYR2WRP)(int *ret, char *uplo, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);
int F_FUNC(dsyr2kwrp, DSYR2KWRP)(int *ret, char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
int F_FUNC(dsyrkwrp, DSYRKWRP)(int *ret, char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc);
int F_FUNC(dtbmvwrp, DTBMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, double *a, int *lda, double *x, int *incx);
int F_FUNC(dtbsvwrp, DTBSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, double *a, int *lda, double *x, int *incx);
int F_FUNC(dtpmvwrp, DTPMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, double *ap, double *x, int *incx);
int F_FUNC(dtpsvwrp, DTPSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, double *ap, double *x, int *incx);
int F_FUNC(dtrmmwrp, DTRMMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
int F_FUNC(dtrmvwrp, DTRMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, double *a, int *lda, double *x, int *incx);
int F_FUNC(dtrsmwrp, DTRSMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
int F_FUNC(dtrsvwrp, DTRSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, double *a, int *lda, double *x, int *incx);
int F_FUNC(dzasumwrp, DZASUMWRP)(double *ret, int *n, npy_complex128 *zx, int *incx);
int F_FUNC(dznrm2wrp, DZNRM2WRP)(double *ret, int *n, npy_complex128 *x, int *incx);
int F_FUNC(icamaxwrp, ICAMAXWRP)(int *ret, int *n, npy_complex64 *cx, int *incx);
int F_FUNC(idamaxwrp, IDAMAXWRP)(int *ret, int *n, double *dx, int *incx);
int F_FUNC(isamaxwrp, ISAMAXWRP)(int *ret, int *n, float *sx, int *incx);
int F_FUNC(izamaxwrp, IZAMAXWRP)(int *ret, int *n, npy_complex128 *zx, int *incx);
int F_FUNC(lsamewrp, LSAMEWRP)(int *ret, char *ca, char *cb);
int F_FUNC(sasumwrp, SASUMWRP)(float *ret, int *n, float *sx, int *incx);
int F_FUNC(saxpywrp, SAXPYWRP)(int *ret, int *n, float *sa, float *sx, int *incx, float *sy, int *incy);
int F_FUNC(scasumwrp, SCASUMWRP)(float *ret, int *n, npy_complex64 *cx, int *incx);
int F_FUNC(scnrm2wrp, SCNRM2WRP)(float *ret, int *n, npy_complex64 *x, int *incx);
int F_FUNC(scopywrp, SCOPYWRP)(int *ret, int *n, float *sx, int *incx, float *sy, int *incy);
int F_FUNC(sdotwrp, SDOTWRP)(float *ret, int *n, float *sx, int *incx, float *sy, int *incy);
int F_FUNC(sdsdotwrp, SDSDOTWRP)(float *ret, int *n, float *sb, float *sx, int *incx, float *sy, int *incy);
int F_FUNC(sgbmvwrp, SGBMVWRP)(int *ret, char *trans, int *m, int *n, int *kl, int *ku, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);
int F_FUNC(sgemmwrp, SGEMMWRP)(int *ret, char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
int F_FUNC(sgemvwrp, SGEMVWRP)(int *ret, char *trans, int *m, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);
int F_FUNC(sgerwrp, SGERWRP)(int *ret, int *m, int *n, float *alpha, float *x, int *incx, float *y, int *incy, float *a, int *lda);
int F_FUNC(snrm2wrp, SNRM2WRP)(float *ret, int *n, float *x, int *incx);
int F_FUNC(srotwrp, SROTWRP)(int *ret, int *n, float *sx, int *incx, float *sy, int *incy, float *c, float *s);
int F_FUNC(srotgwrp, SROTGWRP)(int *ret, float *sa, float *sb, float *c, float *s);
int F_FUNC(srotmwrp, SROTMWRP)(int *ret, int *n, float *sx, int *incx, float *sy, int *incy, float *sparam);
int F_FUNC(srotmgwrp, SROTMGWRP)(int *ret, float *sd1, float *sd2, float *sx1, float *sy1, float *sparam);
int F_FUNC(ssbmvwrp, SSBMVWRP)(int *ret, char *uplo, int *n, int *k, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);
int F_FUNC(sscalwrp, SSCALWRP)(int *ret, int *n, float *sa, float *sx, int *incx);
int F_FUNC(sspmvwrp, SSPMVWRP)(int *ret, char *uplo, int *n, float *alpha, float *ap, float *x, int *incx, float *beta, float *y, int *incy);
int F_FUNC(ssprwrp, SSPRWRP)(int *ret, char *uplo, int *n, float *alpha, float *x, int *incx, float *ap);
int F_FUNC(sspr2wrp, SSPR2WRP)(int *ret, char *uplo, int *n, float *alpha, float *x, int *incx, float *y, int *incy, float *ap);
int F_FUNC(sswapwrp, SSWAPWRP)(int *ret, int *n, float *sx, int *incx, float *sy, int *incy);
int F_FUNC(ssymmwrp, SSYMMWRP)(int *ret, char *side, char *uplo, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
int F_FUNC(ssymvwrp, SSYMVWRP)(int *ret, char *uplo, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);
int F_FUNC(ssyrwrp, SSYRWRP)(int *ret, char *uplo, int *n, float *alpha, float *x, int *incx, float *a, int *lda);
int F_FUNC(ssyr2wrp, SSYR2WRP)(int *ret, char *uplo, int *n, float *alpha, float *x, int *incx, float *y, int *incy, float *a, int *lda);
int F_FUNC(ssyr2kwrp, SSYR2KWRP)(int *ret, char *uplo, char *trans, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
int F_FUNC(ssyrkwrp, SSYRKWRP)(int *ret, char *uplo, char *trans, int *n, int *k, float *alpha, float *a, int *lda, float *beta, float *c, int *ldc);
int F_FUNC(stbmvwrp, STBMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, float *a, int *lda, float *x, int *incx);
int F_FUNC(stbsvwrp, STBSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, float *a, int *lda, float *x, int *incx);
int F_FUNC(stpmvwrp, STPMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, float *ap, float *x, int *incx);
int F_FUNC(stpsvwrp, STPSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, float *ap, float *x, int *incx);
int F_FUNC(strmmwrp, STRMMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb);
int F_FUNC(strmvwrp, STRMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, float *a, int *lda, float *x, int *incx);
int F_FUNC(strsmwrp, STRSMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb);
int F_FUNC(strsvwrp, STRSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, float *a, int *lda, float *x, int *incx);
int F_FUNC(zaxpywrp, ZAXPYWRP)(int *ret, int *n, npy_complex128 *za, npy_complex128 *zx, int *incx, npy_complex128 *zy, int *incy);
int F_FUNC(zcopywrp, ZCOPYWRP)(int *ret, int *n, npy_complex128 *zx, int *incx, npy_complex128 *zy, int *incy);
int F_FUNC(zdotcwrp, ZDOTCWRP)(npy_complex128 *ret, int *n, npy_complex128 *zx, int *incx, npy_complex128 *zy, int *incy);
int F_FUNC(zdotuwrp, ZDOTUWRP)(npy_complex128 *ret, int *n, npy_complex128 *zx, int *incx, npy_complex128 *zy, int *incy);
int F_FUNC(zdrotwrp, ZDROTWRP)(int *ret, int *n, npy_complex128 *cx, int *incx, npy_complex128 *cy, int *incy, double *c, double *s);
int F_FUNC(zdscalwrp, ZDSCALWRP)(int *ret, int *n, double *da, npy_complex128 *zx, int *incx);
int F_FUNC(zgbmvwrp, ZGBMVWRP)(int *ret, char *trans, int *m, int *n, int *kl, int *ku, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx, npy_complex128 *beta, npy_complex128 *y, int *incy);
int F_FUNC(zgemmwrp, ZGEMMWRP)(int *ret, char *transa, char *transb, int *m, int *n, int *k, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *b, int *ldb, npy_complex128 *beta, npy_complex128 *c, int *ldc);
int F_FUNC(zgemvwrp, ZGEMVWRP)(int *ret, char *trans, int *m, int *n, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx, npy_complex128 *beta, npy_complex128 *y, int *incy);
int F_FUNC(zgercwrp, ZGERCWRP)(int *ret, int *m, int *n, npy_complex128 *alpha, npy_complex128 *x, int *incx, npy_complex128 *y, int *incy, npy_complex128 *a, int *lda);
int F_FUNC(zgeruwrp, ZGERUWRP)(int *ret, int *m, int *n, npy_complex128 *alpha, npy_complex128 *x, int *incx, npy_complex128 *y, int *incy, npy_complex128 *a, int *lda);
int F_FUNC(zhbmvwrp, ZHBMVWRP)(int *ret, char *uplo, int *n, int *k, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx, npy_complex128 *beta, npy_complex128 *y, int *incy);
int F_FUNC(zhemmwrp, ZHEMMWRP)(int *ret, char *side, char *uplo, int *m, int *n, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *b, int *ldb, npy_complex128 *beta, npy_complex128 *c, int *ldc);
int F_FUNC(zhemvwrp, ZHEMVWRP)(int *ret, char *uplo, int *n, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx, npy_complex128 *beta, npy_complex128 *y, int *incy);
int F_FUNC(zherwrp, ZHERWRP)(int *ret, char *uplo, int *n, double *alpha, npy_complex128 *x, int *incx, npy_complex128 *a, int *lda);
int F_FUNC(zher2wrp, ZHER2WRP)(int *ret, char *uplo, int *n, npy_complex128 *alpha, npy_complex128 *x, int *incx, npy_complex128 *y, int *incy, npy_complex128 *a, int *lda);
int F_FUNC(zher2kwrp, ZHER2KWRP)(int *ret, char *uplo, char *trans, int *n, int *k, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *b, int *ldb, double *beta, npy_complex128 *c, int *ldc);
int F_FUNC(zherkwrp, ZHERKWRP)(int *ret, char *uplo, char *trans, int *n, int *k, double *alpha, npy_complex128 *a, int *lda, double *beta, npy_complex128 *c, int *ldc);
int F_FUNC(zhpmvwrp, ZHPMVWRP)(int *ret, char *uplo, int *n, npy_complex128 *alpha, npy_complex128 *ap, npy_complex128 *x, int *incx, npy_complex128 *beta, npy_complex128 *y, int *incy);
int F_FUNC(zhprwrp, ZHPRWRP)(int *ret, char *uplo, int *n, double *alpha, npy_complex128 *x, int *incx, npy_complex128 *ap);
int F_FUNC(zhpr2wrp, ZHPR2WRP)(int *ret, char *uplo, int *n, npy_complex128 *alpha, npy_complex128 *x, int *incx, npy_complex128 *y, int *incy, npy_complex128 *ap);
int F_FUNC(zrotgwrp, ZROTGWRP)(int *ret, npy_complex128 *ca, npy_complex128 *cb, double *c, npy_complex128 *s);
int F_FUNC(zscalwrp, ZSCALWRP)(int *ret, int *n, npy_complex128 *za, npy_complex128 *zx, int *incx);
int F_FUNC(zswapwrp, ZSWAPWRP)(int *ret, int *n, npy_complex128 *zx, int *incx, npy_complex128 *zy, int *incy);
int F_FUNC(zsymmwrp, ZSYMMWRP)(int *ret, char *side, char *uplo, int *m, int *n, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *b, int *ldb, npy_complex128 *beta, npy_complex128 *c, int *ldc);
int F_FUNC(zsyr2kwrp, ZSYR2KWRP)(int *ret, char *uplo, char *trans, int *n, int *k, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *b, int *ldb, npy_complex128 *beta, npy_complex128 *c, int *ldc);
int F_FUNC(zsyrkwrp, ZSYRKWRP)(int *ret, char *uplo, char *trans, int *n, int *k, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *beta, npy_complex128 *c, int *ldc);
int F_FUNC(ztbmvwrp, ZTBMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx);
int F_FUNC(ztbsvwrp, ZTBSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, int *k, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx);
int F_FUNC(ztpmvwrp, ZTPMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex128 *ap, npy_complex128 *x, int *incx);
int F_FUNC(ztpsvwrp, ZTPSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex128 *ap, npy_complex128 *x, int *incx);
int F_FUNC(ztrmmwrp, ZTRMMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *b, int *ldb);
int F_FUNC(ztrmvwrp, ZTRMVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx);
int F_FUNC(ztrsmwrp, ZTRSMWRP)(int *ret, char *side, char *uplo, char *transa, char *diag, int *m, int *n, npy_complex128 *alpha, npy_complex128 *a, int *lda, npy_complex128 *b, int *ldb);
int F_FUNC(ztrsvwrp, ZTRSVWRP)(int *ret, char *uplo, char *trans, char *diag, int *n, npy_complex128 *a, int *lda, npy_complex128 *x, int *incx);


#ifdef __cplusplus
}
#endif
#endif
