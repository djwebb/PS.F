      program plot1
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
!  File: plot1.F
!
!  Simple compile:
!
!    gfortran -c ps.F
!    gfortran -o plot1 plot1.F ps.o
!
!=======================================================================
!
      use ps
      implicit none
      integer,          parameter :: DP = kind(1.0d0)
!
      integer  :: i
      real(DP) :: x, y
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
! initialise plotting constants for A4 horizontal
!=======================================================================
!
      call psinit("a4h","POST.eps")
      call query_pagesize(xpmin, xpmax, ypmin, ypmax)
!
!=======================================================================
!  set up plotting space
!=======================================================================
!
      xamin =    0d0
      xamax =   10d0
      yamin =  -1.2d0
      yamax =   1.2d0
!
      xbmin = xamin
      xbmax = xamax
      ybmin = yamin
      ybmax = yamax
!
      dpcm  = 72d0/2.54         !  Convert cm to pixels
      xcmin =  3.0*dpcm         !  Position on page in pixels
      xcmax = 24.0*dpcm
      ycmin =  6.0*dpcm
      ycmax = 17.0*dpcm
!
      call ps_set_a_limits(xamin,xamax,yamin,yamax)
      call ps_set_c_limits(xcmin,xcmax,ycmin,ycmax)
!
!-----------------------------------------------------------------------
!  Annotate plot (Outline drawn at end)
!-----------------------------------------------------------------------
!
      call setfont(5,10)
      call setlinewidth(0.5d0)
!  Draw axes and tick marks
      call ticks(1,0d0,0.2d0,5,10)
      call ticks(2,0d0,0.1d0,5,5)
!  Position axis titles
      call setfont(5,12)
      call pstext(xcmax-50.0,ycmin-40.0,"Time")
      call pstext3(xcmin-30d0,ycmax-35d0,90d0,"Speed",1,0)
!
!  Draw zero line
!
      call setgrey(0.8d0)
      call move2(xamin,0d0)
      call draw2(xamax,0d0)
      call stroke
      call setblack
!
!-----------------------------------------------------------------------
! Plot data
!-----------------------------------------------------------------------
!
      call move2(0d0,0d0)
      do i= 1,100
        x = i*0.1d0
        y = sin(x)
        call draw2(x,y)
      enddo
!
!-----------------------------------------------------------------------
! draw plotting box
!-----------------------------------------------------------------------
!
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

      end program plot1
