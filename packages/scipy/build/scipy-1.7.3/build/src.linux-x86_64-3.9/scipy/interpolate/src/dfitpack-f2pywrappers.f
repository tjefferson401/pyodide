C     -*- fortran -*-
C     This file is autogenerated with f2py (version:1.21.4)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapsplint (splintf2pywrap, t, n, c, k, a, b,
     & wrk)
      external splint
      integer n
      integer k
      real*8 a
      real*8 b
      real*8 t(n)
      real*8 c(n)
      real*8 wrk(n)
      real*8 splintf2pywrap, splint
      splintf2pywrap = splint(t, n, c, k, a, b, wrk)
      end


      subroutine f2pywrapdblint (dblintf2pywrap, tx, nx, ty, ny, c
     &, kx, ky, xb, xe, yb, ye, wrk)
      external dblint
      integer nx
      integer ny
      integer kx
      integer ky
      real*8 xb
      real*8 xe
      real*8 yb
      real*8 ye
      real*8 tx(nx)
      real*8 ty(ny)
      real*8 c((nx-kx-1)*(ny-ky-1))
      real*8 wrk(nx+ny-kx-ky-2)
      real*8 dblintf2pywrap, dblint
      dblintf2pywrap = dblint(tx, nx, ty, ny, c, kx, ky, xb, xe, y
     &b, ye, wrk)
      end


      subroutine f2pyinittypes(setupfunc)
      external setupfunc
      integer intvar
      common /types/ intvar
      call setupfunc(intvar)
      end


