      program plot2
!
!======================================================================
!
!    Copyright 2020  David J. Webb
!
!    This file is part of PS.F.
!
!    PS.F is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PS.F is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PS.F.  If not, see <https://www.gnu.org/licenses/>.
!
!======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  File: plot2.F
!
!  Simple compile:
!
!    gfortran -c ps.F
!    gfortran -o plot2 plot2.F ps.o
!
!=======================================================================
!
      use ps
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer, parameter :: isize = 400
!
      integer  :: i, j, k, ic(isize,isize), icol, isizeh
      real(DP) :: x, y, pi, dd, ra(64), ga(64), ba(64),
     &                          bx(isize+1),by(isize+1)
      complex(DP) :: ci, ca, z
!
!=======================================================================
!   ps variables
!=======================================================================
!
      real(DP) ::  xamin,xamax,yamin,yamax,xbmin,xbmax,ybmin,ybmax,
     &             xcmin,xcmax,ycmin,ycmax,xmmin,xmmax,ymmin,ymmax,
     &             xpmin,xpmax,ypmin,ypmax,dpcm
!
!=======================================================================
!   Generate Mandelbrot function
!=======================================================================
!
      pi = 4d0*atan(1d0)
      ci = (0d0,1d0)
      dd = 1d0/isize
      isizeh = isize/2

      ca = -0.909 -ci*0.275
      do j = 1,isize
      do i = 1,isize
        k = 0
        z = ca + (i-isizeh)*dd +(j-isizeh)*dd*ci
        k = 0
        do
         k = k + 1
         if(abs(z).gt.2d0.or.k.ge.64)exit
         z = z*z + ca
        enddo
        ic(i,j) = k
      enddo
      enddo

!
!=======================================================================
! initialise plotting constants for A4 horizontal
!=======================================================================
!
      call psinit("a4v","POST.eps")
      call query_pagesize(xpmin, xpmax, ypmin, ypmax)
!
!=======================================================================
!  set up plotting space
!=======================================================================
!
      xamin =    0d0
      xamax =  isize
      yamin =    0d0
      yamax =  isize
!
      xbmin = xamin
      xbmax = xamax
      ybmin = yamin
      ybmax = yamax
!
      dpcm  = 72d0/2.54         !  Convert cm to pixels
      xcmin =  5.0*dpcm         !  Position on page in pixels
      xcmax = 15.0*dpcm
      ycmin =  3.0*dpcm
      ycmax = 13.0*dpcm
!
      call ps_set_a_limits(xamin,xamax,yamin,yamax)
      call ps_set_c_limits(xcmin,xcmax,ycmin,ycmax)
!
!-----------------------------------------------------------------------
! Initialise spectrum and specify boundaries
!-----------------------------------------------------------------------
!
      do i=1,63
        x = i
        call col_steps1(x,1d0,64d0,64,ra(i),ga(i),ba(i),icol,1)
!        print *,i,ra(i),ga(i),ba(i),icol
      enddo
      ra(64) = 0d0 ; ga(64) = 0d0; ba(64) = 0d0
      do i=1,isize+1
        bx(i) = i-1
        by(i) = i-1
      enddo
!
!-----------------------------------------------------------------------
! Plot data
!-----------------------------------------------------------------------
!
      call contc(ic,isize,isize,isize,bx,by,64,ra,ga,ba)
!
!-----------------------------------------------------------------------
! draw plotting box
!-----------------------------------------------------------------------
!
      call setlinewidth(1.0d0)
      call ps_use_c_space
      call move2(xcmin,ycmax)
      call draw2(xcmin,ycmin)
      call draw2(xcmax,ycmin)
      call draw2(xcmax,ycmax)
      call closepath
      call stroke
!
!=======================================================================
!  set up second plotting space for ortho projection
!=======================================================================
!
      xamin =     0d0           !  Lat/Long limits of plot
      xamax =   100d0
      yamin =   -50d0
      yamax =    50d0
!
      xbmin = -1.2d0            ! Ortho globe has radius 1
      xbmax =  1.2d0
      ybmin = -1.2d0
      ybmax =  1.2d0
!
      dpcm  = 72d0/2.54         !  Convert cm to pixels
      xcmin =  5.0*dpcm         !  Position on page in pixels
      xcmax = 15.0*dpcm
      ycmin = 14.0*dpcm
      ycmax = 24.0*dpcm
!
      call ps_set_a_limits(xamin,xamax,yamin,yamax)
      call ps_set_b_limits(xbmin,xbmax,ybmin,ybmax)
      call ps_set_c_limits(xcmin,xcmax,ycmin,ycmax)

      call set_ortho(45d0,0d0,1)
      call setmap("ortho")

      do i=1,isize+1
        bx(i) = (i-1)*100d0/isize
        by(i) = (i-1)*100d0/isize - 50d0
      enddo
!
!-----------------------------------------------------------------------
!  Plot data.  With the mercator projection and this ortho projection
!              centred on the equator conc works without error
!              but when using a projection which curves horizontal
!              lines each individual cell should be drawn separately
!-----------------------------------------------------------------------
!
      call contc(ic,isize,isize,isize,bx,by,64,ra,ga,ba)
!
!-----------------------------------------------------------------------
!  plot grid
!-----------------------------------------------------------------------
!
      call setlinewidth(0.5d0)
      do j=-80,80,10
        y = j
        call newpath
        call move2(0d0,y)
        do i=1,360
          x = i
          call draw2(x,y)
        enddo
        call stroke
      enddo

      do i=0,350,10
        x = i
        call newpath
        call move2(x,-90d0)
        do j=-89,90
          y = j
          call draw2(x,y)
        enddo
        call stroke
      enddo
!
!-----------------------------------------------------------------------
! draw plotting box
!-----------------------------------------------------------------------
!
      call setlinewidth(1.0d0)
      call ps_use_c_space
      call move2(xcmin,ycmax)
      call draw2(xcmin,ycmin)
      call draw2(xcmax,ycmin)
      call draw2(xcmax,ycmax)
      call closepath
      call stroke
!
      call showpage
      call psclose

      stop

      end program plot2
