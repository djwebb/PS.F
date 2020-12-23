!
!=======================================================================
!
!    File: ps.F
!    Copyright 2020  David J. Webb
!
!    This file is part of PS.F.
!    Released under licence GPL-3.0-or-later
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
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  Module version of my psplot routines
!
!  Compile to ps.o and ps.mod (psplot module)
!  i.e. "gfortran -c  ps.F"
!
!  Notes:  1.  In text use \ to escape characters ()\
!              i.e. \( \) \\  or for character codes \123
!======================================================================
!
      module ps
      public
      private         :: ps_plot, isnan2, pssubs_print_on_nan
      integer, parameter, private :: DP = kind(1.0d0)

      integer,private   :: iopost,ifontsize,imap0,imap1,imap2
      logical,private   :: lpsmask,lpvalid
      real(DP),private  :: xamin,xamax,yamin,yamax,
     &                     xbmin,xbmax,ybmin,ybmax,
     &                     xcmin,xcmax,ycmin,ycmax,
     &                     xmmin,xmmax,ymmin,ymmax,
     &                     xpmin,xpmax,ypmin,ypmax
      real(DP), private :: alpha = 0d0, beta = 0d0, gamma = 0d0
      real(DP), private :: orthoa(3,3), ortho_lon, ortho_lat, ortho_z

      type :: ps_node
        type(ps_node), pointer :: prev_node => null()
        integer   :: iopost,ifontsize,imap0,imap1,imap2
        logical   :: lpsmask,lpvalid
        real(DP)  :: xamin,xamax,yamin,yamax,
     &               xbmin,xbmax,ybmin,ybmax,
     &               xcmin,xcmax,ycmin,ycmax,
     &               xmmin,xmmax,ymmin,ymmax,
     &               xpmin,xpmax,ypmin,ypmax
        real(DP) :: alpha, beta, gamma
      end type ps_node

      type(ps_node), pointer :: curr_node=>null(), last_node =>null()

!
!======================================================================
!  Interfaces
!======================================================================
!
      interface psinit
        module procedure psinitf, psiniti, psinita
      end interface
      interface psplot
        module procedure ps_plot
      end interface

      contains

!
!======================================================================
!  Routines to open and close the file
!======================================================================
!

      subroutine psinita(page, filename)
      implicit none
      character(len=*), intent(in) :: page
      character(len=*), intent(in), optional :: filename
!
      integer :: idirn
      real(DP)  :: x1, y1, x2, y2, dpc=(72d0/2.54d0)
!
      x1 =  0.5*dpc ; y1 =  0.5*dpc
      if(page.eq."a4v")then
        y2 = 29.5*dpc ; x2 = 21.2*dpc ; idirn = 0
      elseif(page.eq."a4h")then
        x2 = 29.5*dpc ; y2 = 21.2*dpc ; idirn = 1
      else
        print *," Subroutine psinit/psinita error."
        print *," Page type: ",page
        print *,"           not recodnised."
        print *," Program stopping ... "
        stop
      endif
!
      if(present(filename))then
        call psinit(x1,y1,x2,y2,idirn,filename)
      else
        call psinit(x1,y1,x2,y2,idirn,"POST.eps")
      endif
      return
      end subroutine psinita
!
!======================================================================
!
      subroutine psiniti(rllx,rlly,rurx,rury,isize)
      implicit none
!
!      open output postscript file 'POST' and initialise
!
      integer, intent(in) :: isize
      real(DP), intent(in)  :: rllx,rlly,rurx,rury
      character*8 file1
!
      file1 = 'POST.eps'
      call psinitf(rllx,rlly,rurx,rury,isize,file1)
      return
      end subroutine
!
!======================================================================
!
      subroutine psinitf(rllx,rlly,rurx,rury,isize,file1)
      implicit none
!
!      postscript header information
!
!  If isize = 0, plot is a vertical A4 plot
!           = 1, plot is a horizontal A4 plot
!
      integer, intent(in) :: isize
      real(DP), intent(in)  :: rllx,rlly,rurx,rury
      character*(*), intent(in) :: file1
!
      real(DP) :: dpc, a4s, a4l
      integer  :: iend
!
!      integer :: iopost,ifontsize,imap0,imap1,imap2
!      logical lpsmask,lpvalid
!      common /psi/ iopost,ifontsize,imap0,imap1,imap2,lpsmask,lpvalid
!
      dpc = 72d0/2.54d0                    ! dots or points per centimeter
      a4s = 21d0*dpc                       ! A4 short side
      a4l = a4s*sqrt(2d0)                  ! A4 long side
      if(isize.lt.0.or.isize.gt.1)then
        print *,' routine psinit.  isize out of range '
        print *,' program stopping ...'
        stop
      endif
!
      iopost = 41
      lpsmask = .false.
      lpvalid = .false.
      ifontsize = 12
      imap0 = 1
      imap1 = 0
      imap2 = 0

      iend = len_trim(file1)
      if(iend.ne.0)then
        iend = index(file1,".",.true.)-1
      endif
      if(iend.eq.0)then
        print *," subroutine psinit : error"
        print *," output filename has zero length"
        stop
      endif

      xpmin = rllx
      xpmax = rurx
      ypmin = rlly
      ypmax = rury
!
!  gcc (f95) may overwrite a file leaving rubbish at the end.
!
      open(iopost,file=file1,access='sequential',
     &            form='formatted',status='unknown')
      rewind(iopost)
      write(iopost,11)file1(1:iend)
!
      if(isize.eq.0)then
        write(iopost,12)nint(rllx),nint(rlly),nint(rurx),nint(rury)
        write(iopost,13)
        write(iopost,14)ifontsize
      elseif(isize.eq.1)then
        write(iopost,12)nint(a4s-rury),nint(rllx),
     &                  nint(a4s-rlly),nint(rurx)
        write(iopost,13)
        write(iopost,14)ifontsize
        write(iopost,15)a4s,0
        write(iopost,16)
      endif
      write(iopost,17)
      write(iopost,171)
      write(iopost,18)
      write(iopost,19)
      write(iopost,20)
      write(iopost,21)
!
      call setlinejoin(0)
      call setlinecap(0)
      call setmiterlimit(2d0)
      return
!
  11  format('%!PS-Adobe-1.0',/,
     &       '%%Creator: David Webb',/,
     &       '%%Title: ps: ',a)
  12  format('%%BoundingBox:',4(1x,i0))
  13  format('%%Pages: 1',/,
     &       '%%EndComments')
  14  format('/Helvetica findfont',/,
     &       i2,' scalefont',/,
     &       'setfont')
  15  format(f0.3,1x,i1,' translate')
  16  format('90 rotate')
  17  format('/m {moveto} def')
  171 format('/r {rmoveto} def')
  18  format('/l {lineto} def')
  19  format('/n {newpath} def')
  20  format('/s {stroke} def')
  21  format('%%EndProlog',/,
     &       '%%Page: 1 1')
      end subroutine
!
!======================================================================
!
      subroutine showpage
      implicit none
!
!  Draw whole page
!
      if(lpvalid)call stroke
      write(iopost,11)
  11  format('showpage')
      return
      end subroutine
!
!======================================================================
!
      subroutine psclose
      implicit none
!
!   Close postscript file with trailer line
!
      write(iopost,11)
      close(iopost)
      iopost = -1
      return
  11  format("%%Trailer",/)
      end subroutine psclose
!
!======================================================================
!  Setting the a, b and c space limits
!======================================================================
!
!
!  set "a" limits - limits of plotting region on paper
!
      subroutine ps_set_a_limits(x1,x2,y1,y2)
      implicit none
      real(DP), intent(in) :: x1, x2, y1, y2
      xamin = x1 ; xamax = x2
      yamin = y1 ; yamax = y2
      if(xbmin.eq.0d0)xbmin = xamin
      if(xbmax.eq.0d0)xbmax = xamax
      if(ybmin.eq.0d0)ybmin = yamin
      if(ybmax.eq.0d0)ybmax = yamax
      return
      end subroutine
!
!======================================================================
!
!  set "b" limits - limits of plotting region after map transformation
!                 - usually equals a limits
!
      subroutine ps_set_b_limits(x1,x2,y1,y2)
      implicit none
      real(DP), intent(in) :: x1, x2, y1, y2
      xbmin = x1 ; xbmax = x2
      ybmin = y1 ; ybmax = y2
      return
      end subroutine
!
!======================================================================
!
!  Set "c" limits - limits of plotting region in on base map
!
      subroutine ps_set_c_limits(x1,x2,y1,y2)
      implicit none
      real(DP), intent(in) :: x1, x2, y1, y2
      xcmin = x1 ; xcmax = x2
      ycmin = y1 ; ycmax = y2
      return
      end subroutine
!
!======================================================================
!
      subroutine ps_use_a_space()
      imap0 = 1
      return
      end subroutine
!
!======================================================================
!
      subroutine ps_use_c_space()
      imap0 = 0
      return
      end subroutine
!
!======================================================================
!
      subroutine query_pagesize(xmin, xmax, ymin, ymax)
      implicit none
!
!  Subroutine to return current page size
!
      real(DP), intent(out) ::  xmin, xmax, ymin, ymax
!
      xmin = xpmin
      xmax = xpmax
      ymin = ypmin
      ymax = ypmax
!      return
      end subroutine query_pagesize
!
!======================================================================
!
      subroutine pstran(xa,ya,xc,yc)
      implicit none
!
!   routine to transfer from user coordinates to page (pixel)
!   coordinates
!
      real(DP), intent(in)  :: xa, ya
      real(DP), intent(out) :: xc, yc
!
      real(DP)  :: xb, yb
      integer, save :: icount=0
!
! 1.  Convert from user coordinates (longitude, latitude)
!             to linear base coordinates
!
      if(imap1.eq.0)then
        xb = xa
        yb = ya
      elseif(imap1.eq.1)then
        call transeu(xa,ya,xb,yb,1,1,alpha,beta,gamma,1)
        if(xb.lt.xbmin)then
          xc = xb+360d0
          if(xc-xbmax.lt.xbmin-xb)xb = xc
        elseif(xb.gt.xbmax)then
          xc = xb-360d0
          if(xbmin-xc.lt.xb-xbmax)xb = xc
        endif
!        icount = icount + 1
!        if(icount.eq.1)print *,'a-b-c',alpha,beta,gamma
!        if(icount.lt.30)then
!          print *,xa,ya,xb,yb
!        endif
      else
        print *,' error calling pstran - variable imap1 out of range'
        print *,' imap1 = ', imap1
        print *,' programme stopping ...'
      endif
!
      if(imap2.eq.1)then
        call mercator(xb,yb,1)
      elseif(imap2.eq.2)then
        call ortho(xb,yb,1)
      elseif(imap2.ne.0)then
        print *,' error calling pstran - variable imap2 out of range'
        print *,' imap2 = ', imap2
        print *,' programme stopping ...'
      endif
!
      xc = xcmin + (xb-xbmin)*(xcmax-xcmin)/(xbmax-xbmin)
      yc = ycmin + (yb-ybmin)*(ycmax-ycmin)/(ybmax-ybmin)
      return
      end subroutine
!
!======================================================================
!  Simple plot commands
!======================================================================
!
      subroutine newpath
      implicit none
!
      write(iopost,11)
  11  format("n")
      return
      end subroutine newpath
!
!======================================================================
!
      subroutine move2(x,y)
      implicit none
      real(DP), intent(in) :: x,y
!
      call ps_plot(x,y,0)
      return
      end subroutine move2
!
!======================================================================
!
      subroutine draw2(x,y)
      implicit none
      real(DP), intent(in) :: x,y
!
      call ps_plot(x,y,1)
      return
      end subroutine draw2
!
!======================================================================
!
      subroutine closepath
      implicit none
!
      write(iopost,11)
  11  format('closepath')
      return
      end subroutine closepath
!
!======================================================================
!
      subroutine endline
      implicit none
!
      call stroke
      return
      end subroutine endline
!
!======================================================================
!
      subroutine stroke
      implicit none

      if(lpvalid)then
        write(iopost,11)
        lpvalid = .false.
      endif
  11  format('s')
      return
      end subroutine stroke
!
!======================================================================
!
      subroutine fill
      implicit none

      if(lpvalid)then
        write(iopost,11)
        lpvalid = .false.
      endif
  11  format('fill')
      return
      end subroutine fill
!
!  Special Postscript commands
!

      subroutine arc(x, y, r, a1, a2)
      implicit none
!
!  postscript arc command
!
      real(DP), intent(in) :: x, y, r, a1, a2
      write(iopost,11)x,y,r,a1,a2
  11  format(5(f10.3,1x),"arc")
      return
      end subroutine arc
!
!======================================================================
!
      subroutine arcn(x, y, r, a1, a2)
      implicit none
!
!  postscript arcn command
!
      real(DP), intent(in) :: x, y, r, a1, a2
      write(iopost,11)x,y,r,a1,a2
  11  format(5(f10.3,1x),"arcn")
      return
      end subroutine arcn
!
!======================================================================
!
      subroutine arcto(x1, y1, x2, y2, r)
      implicit none
!
!  postscript arc command
!
      real(DP), intent(in) :: x1, y1, x2, y2, r
      real(DP) :: xx1, yy1, xx2, yy2, rr, ratio

      call pstran(x1,y1,xx1,yy1)
      call pstran(x2,y2,xx2,yy2)
      ratio = sqrt(((xx2-xx1)**2+(yy2-yy1)**2)
     &                /((x2-x1)**2+(y2-y1)**2))
      rr = r*ratio
      write(iopost,11)xx1, yy1, xx2, yy2, rr
      write(iopost,12)
  11  format(5(f0.3,1x),"arcto")
  12  format("4 {pop} repeat")        !  Remove returned values
      return
      end subroutine arcto
!
!======================================================================
!  Specify Colours
!======================================================================
!
      subroutine setrgb(r,g,b)
      implicit none
      real(DP), intent(in) :: r, g, b

      write(iopost,11)r,g,b
  11  format(3(f0.3,1x),'setrgbcolor')
      return
      end subroutine setrgb
!
!======================================================================
!
      subroutine setgray(gray)
      implicit none
!
!    white -> grey = 1d0
!    black -> grey = 0d0
!
      real(DP), intent(in) :: gray
!
      write(iopost,11)gray
  11  format(f0.3,' setgray')
      return
      end subroutine setgray
!
!======================================================================
!
      subroutine setgrey(gray)
      implicit none
      real(DP), intent(in) :: gray
!
      write(iopost,11)gray
  11  format(f0.3,' setgray')
      return
      end subroutine setgrey
!
!======================================================================
!
      subroutine setblack()
      call setgray(0d0)
      end subroutine setblack
!
!======================================================================
!
      subroutine setwhite()
      call setgray(1d0)
      end subroutine setwhite

!
!======================================================================
!  Specify Line Properties!
!======================================================================
!
      subroutine setlinewidth(w)
      implicit none
      real(DP), intent(in) :: w
!
!  The default width is 1.0
!
      write(iopost,11)w
  11  format(f0.3,' setlinewidth')
      return
      end subroutine setlinewidth
!
!======================================================================
! Routine to set a dashed line
!
      subroutine setdash(l1,l2)
      implicit none
      integer, intent(in) :: l1  ! length of 'black' segment in pixels
      integer, intent(in) :: l2  ! length of 'white' segment in pixels

      if(l1.eq.0)then
        write(iopost,11)
      else
        write(iopost,12)l1,l2
      endif
      return
  11  format("[] 0 setdash")
  12  format("[ ",i0," ",i0,"] 0 setdash")
      end subroutine setdash
!
!======================================================================
!
      subroutine setlinejoin(ijoin)
      implicit none
      integer, intent(in) :: ijoin
!
      write(iopost,11)ijoin
  11  format(i0,' setlinejoin')
      return
      end subroutine setlinejoin
!
!======================================================================
!
      subroutine setlinecap(icap)
      implicit none
      integer, intent(in) :: icap
!
      write(iopost,11)icap
  11  format(i0,' setlinecap')
      return
      end subroutine setlinecap
!
!======================================================================
!
      subroutine setmiterlimit(ratio)
      implicit none
!
      real(DP), intent(in) :: ratio
!
      write(iopost,11)ratio
  11  format(f0.3,' setmiterlimit')
      return
      end subroutine setmiterlimit
!
!======================================================================
!  Masks
!======================================================================
!
      subroutine setpsmask(x1,x2,y1,y2)
      implicit none
      real(DP) x1, x2, y1, y2
!
      call ps_set_mgrid(x1,x2,y1,y2)
      return
      end subroutine setpsmask
!
!======================================================================
!
      subroutine psmask(ii)
      implicit none
      integer :: ii

      if(ii.eq.0)then
         lpsmask = .false.
      elseif(ii.eq.1)then
         lpsmask = .true.
      else
        print *,' psmask called with invalid mask value:',ii
      endif
      return
      end subroutine psmask
!
!======================================================================
!
      subroutine ps_set_mgrid(x1,x2,y1,y2)
      implicit none
      real(DP), intent(in) :: x1, x2, y1, y2
      xmmin = x1 ; xmmax = x2
      ymmin = y1 ; ymmax = y2
      lpsmask = .false.
      return
      end subroutine ps_set_mgrid
!
!======================================================================
!  Plot text
!======================================================================
!
      subroutine pstext(xx,yy,string,irel)
      implicit none
      real(DP),intent(in) ::  xx, yy
      integer, optional, intent(in) :: irel
      character*(*),intent(in) :: string
!
      integer :: ii
!
!      integer :: iopost,ifontsize,imap0,imap1,imap2
!      logical lpsmask,lpvalid
!      common /psi/ iopost,ifontsize,imap0,imap1,imap2,lpsmask,lpvalid
!
      ii = max(lnblnk(string),1)
      if(present(irel))then
        write(iopost,10) xx,yy
      else
        write(iopost,11) xx,yy
      endif
      write(iopost,12) string(1:ii)
  10  format(2(f0.3,1x),'r')
  11  format(2(f0.3,1x),'m')
  12  format('(',a,') show')
      return
      end subroutine pstext
!
!======================================================================
!
      subroutine pstext2(xx,yy,string,iposn)
      implicit none
      integer :: iposn
      real(DP)     xx, yy
      character*(*) string
!
!  Position non-blank string relative to xx,yy
!  If iposn = 0, start at xx,yy
!             1, centre on xx,yy
!             2, end at xx,yy
!
      integer i1, i2, im
!      integer :: iopost,ifontsize,imap0,imap1,imap2
!      logical lpsmask,lpvalid
!      common /psi/ iopost,ifontsize,imap0,imap1,imap2,lpsmask,lpvalid
!
      i2 = max(lnblnk(string),1)
      i1 = 1
      im = len(string)
  10  if(string(i1:i1).eq.' '.and.i1.lt.im)then
        i1 = i1+1
        goto 10
      endif
!
      if(iposn.eq.0)then
        write(iopost,11) xx,yy
        write(iopost,12) string(i1:i2)
      elseif(iposn.eq.1)then
        write(iopost,21) string(i1:i2)
        write(iopost,22) xx,yy
      elseif(iposn.eq.2)then
        write(iopost,31) string(i1:i2)
        write(iopost,22) xx,yy
      else
        print *,' pstext2 error:  iposn out of range.'
      endif
      return
!
  11  format(2(f0.3,1x),'m')
  12  format('(',a,') show')
  21  format('(',a,') dup stringwidth pop 2 div neg')
  22  format(f0.3,' add ',/,f0.3,' m show ')
  31  format('(',a,') dup stringwidth pop neg')
      end subroutine pstext2
!
!======================================================================
!
      subroutine pstext3(xx,yy,angle,string,iposn,iback)
      implicit none
!
!  Position non-blank string relative to xx,yy with angle 'angle'
!  If iposn = 0, bottom left at xx,yy
!             1, bottom centre on xx,yy
!             2, bottom right at xx,yy
!             3, middle centre at xx,yy
!  If iback = 0,  plot text only
!             1, plot white background before plotting text
!             2, plot white background only
!
      integer :: iposn, iback
      real(DP)     xx, yy, angle
      character*(*) string
!
      integer i1, i2, im
!
      logical lgsave,langle,lstring
!      integer :: iopost,ifontsize,imap0,imap1,imap2
!      logical lpsmask,lpvalid
!      common /psi/ iopost,ifontsize,imap0,imap1,imap2,lpsmask,lpvalid
!
!      print *,'pstext3 - ',xx,yy,angle,iposn
!      print *,'pstect3   ',string(1:lnblnk(string))
      i2 = max(lnblnk(string),1)
      i1 = 1
      im = len(string)
      do while(string(i1:i1).eq.' '.and.i1.lt.i2)
        i1 = i1+1
      end do
!
! 1 move to position and put string on stack
!
      write(iopost,10)
      write(iopost,11) xx,yy,string(i1:i2)
  10  format('%% pstext3')
  11  format(2(f0.3,1x),'m ','(',a,')')
!
! 2 add (a) gsave, (b) angle, (c) stringlength and height if required
!
      langle = angle .ne. 0d0
      lstring = iback.ne.0 .or. iposn.ne.0
      lgsave = langle .or. lstring
      if(lgsave) write (iopost,12)
      if(langle) write (iopost,13) angle
      if(lstring) write (iopost,14) ifontsize*2./3.
  12  format('gsave')
  13  format(f0.3,' rotate')
  14  format('dup stringwidth pop ',f0.3)
!
!  3  move to start of string
!     if lstring: stack = string stringlength stringheight
!     else:       stack = string
!
      if(iposn.eq.1)then
        write(iopost,20)
  20  format('2 copy exch -.5 mul exch rmoveto')
      elseif(iposn.eq.2)then
        write(iopost,21)
  21  format('2 copy exch neg exch rmoveto')
      elseif(iposn.eq.3)then
        write(iopost,22)
  22  format('2 copy -.5 mul exch -.5 mul exch rmoveto')
      endif
!
!  4 plot background, remove stringlength and height
!
      if(iback.ne.0)then
        write(iopost,30)
        write(iopost,31)
        write(iopost,32)
        write(iopost,33)
      elseif(lstring)then
        write(iopost,34)
      endif
  30  format('gsave currentpoint translate .9 setgray -1 -1 moveto ')
  31  format('1 index 2 add 0 rlineto 0 1 index 2 add rlineto ')
  32  format('pop 2 add neg 0 rlineto closepath fill ')
  33  format('0 0 moveto grestore ')
  34  format('pop pop')
!
!  5  plot text
!
      if(iback.ne.2)then
        write(iopost,41)
  41  format('show')
      else
        write(iopost,42)
  42  format('pop')
      endif
!
!  6 restore
!
      if(lgsave)then
        write(iopost,51)
  51  format('grestore')
      endif
!
      return
      end subroutine pstext3
!
!======================================================================
!
      subroutine setfont (ifont, isize)
      implicit none
!
      integer :: ifont, isize
!
      character*32 fontname(13)
      data fontname /'Times-Roman'      ,'Times-Bold',
     &               'Times-Italic'     ,'Times-BoldItalic',
     &               'Helvetica'        ,'Helvetica-Bold',
     &               'Helvetica-Oblique','Helvetica-BoldOblique',
     &               'Courier-Roman'    ,'Courier-Bold',
     &               'Courier-Oblique'  ,'Courier-BoldOblique',
     &               'Symbol'/
!
      if(ifont.lt.1.or.ifont.gt.13)then
        print *,' routine setfont: font index number'//
     &                         ' out of range (0-13)'
        print *,' ifont = ',ifont
        return
      endif
      if(isize .lt.6.or.isize.gt.64)then
        print *,' routine setfont: font size out of range (6-64)'
        print *,' isize = ',isize
        return
      endif
!
      write(iopost,11)trim(fontname(ifont)),isize
  11  format('/',a,' findfont ',i0,' scalefont setfont')
      ifontsize = isize
      return
      end subroutine setfont
!
!======================================================================
!  Tick marks
!======================================================================
!
      subroutine ticks(iaxis,tbase,tincr,ntick,nannot,rfmt,ifmt,sign_a)
      implicit none
!
!  subroutine to draw tickmarks on a plot
!
!    iaxis  = 1 to annotate bottom horizontal axis
!    iaxis  = 2 to annotate left vertical axis
!    tbase  = base value for ticks
!    tincr  = increment value between ticks
!    ntick  = number of ticks between large ticks
!    nannot = number of ticks between annotations
!  Optional parameters
!    rfmt   = format for annotation     (to override default)
!    ifmt   = integer format for annotation
!    sign_a = change sign of annotation (may be used to make
!                                               depths positive)
!
      integer  :: iaxis, ntick, nannot
      real(DP) :: tbase, tincr
      integer,          optional :: sign_a
      character(len=*), optional :: rfmt, ifmt
!
      integer  :: i, imapp, imax, imin
      real(DP) :: x, y, dpc, xoff, yoff, xl, yl, tannot, xx, yy, dd
      character(len=8) :: string8
!
!
      dpc = 72d0/2.54d0
      yoff = 0.6d0
      xoff = 0.4d0
      xl = 0.15d0*dpc*abs(xbmax-xbmin)/(xcmax-xcmin)
      yl = 0.15d0*dpc*abs(ybmax-ybmin)/(ycmax-ycmin)
!      print *,xmax,xmin,xcmax,xcmin,xl
!      print *,ymax,ymin,ycmax,ycmin,xl
      if(tincr.eq.0d0)then
        print *,'error in calling routine ticks'
        print *,'  tincr equals zero'
        return
      endif
!
!  save value of imap0 and turn off mapping by ps_plot
!
      imapp = imap0
      imap0 = 0
      tannot = abs(tincr*nannot)
!
!  x axis
!
      if(iaxis.eq.1)then
        imin = nint((xamin-tbase)/tincr)
        imax = nint((xamax-tbase)/tincr)
        do i=imin,imax
          x = tbase + i*tincr
          if(x.ge.xamin.and.x.le.xamax)then
            call pstran(x,yamin,xx,yy)
            dd = -0.15*dpc
            call ps_plot(xx,yy,0)
            if(ntick.ne.0)then
              if(mod(i,ntick).eq.0)dd = dd+dd
            endif
            call ps_plot(xx,yy+dd,1)
            call stroke
            if(nannot.ne.0.and.mod(i,nannot).eq.0)then
              if(x.ne.0d0.and.present(sign_a))x = -x
              if(present(ifmt))then
                write(string8,ifmt)nint(x)
              elseif(present(rfmt))then
                write(string8,rfmt)x
              elseif(abs(x).lt.0.001d0*abs(tincr))then
                write(string8,'("0.0")')
              elseif(abs(x).gt.1d5.or.abs(x).lt.0.001d0)then
                write(string8,'(1pe8.1)')x
              elseif(tannot.ge.0.1d0)then
                write(string8,'(f8.1)')x
              elseif(tannot.ge.0.01d0)then
                write(string8,'(f8.2)')x
              else
                write(string8,'(f8.3)')x
              endif
              call pstext2(xx,yy-yoff*dpc,string8,1)
            endif
          endif
        end do
!
!  y axis
!
      elseif(iaxis.eq.2)then
        imin = nint((yamin-tbase)/tincr)
        imax = nint((yamax-tbase)/tincr)
        do i=imin,imax
          y = tbase + i*tincr
          if(y.ge.yamin.and.y.le.yamax)then
            call pstran(xamin,y,xx,yy)
            dd = -0.15*dpc
            call ps_plot(xx,yy,0)
            if(ntick.ne.0)then
              if(mod(i,ntick).eq.0)dd = dd+dd
            endif
            call ps_plot(xx+dd,yy,1)
            call stroke
            if(nannot.ne.0.and.mod(i,nannot).eq.0)then
              if(y.ne.0d0.and.present(sign_a))y = -y
              if(present(ifmt))then
                write(string8,ifmt)nint(y)
              elseif(present(rfmt))then
                write(string8,rfmt)y
              elseif(abs(y).lt.0.001d0*abs(tincr))then
                write(string8,'("0.0")')
              elseif(abs(y).gt.1d5)then
                write(string8,'(1pe8.1)')y
              elseif(tannot.ge.0.1d0)then
                write(string8,'(f8.1)')y
              elseif(tannot.ge.0.01d0)then
                write(string8,'(f8.2)')y
              else
                write(string8,'(f8.3)')y
              endif
              call pstext2(xx-xoff*dpc,yy-0.15*dpc,string8,2)
            endif
          endif
        end do
      endif
!
!  reset value of imap0
!
      imap0 = imapp
!
      return
      end subroutine ticks
!
!======================================================================
!
      subroutine ticks2(iaxis,tbase,tincr,ntick,nannot,rfmt,ifmt,
     &                  sign_a, tick_len, xoff1, yoff1)
      implicit none
!
!  subroutine to draw tickmarks on a plot
!
!    iaxis  = 1 to annotate bottom horizontal axis
!    iaxis  = 2 to annotate left vertical axis
!    tbase  = base value for ticks
!    tincr  = increment value between ticks
!    ntick  = number of ticks between large ticks
!    nannot = number of ticks between annotations
!  Optional parameters
!    rfmt   = format for annotation     (to override default)
!    ifmt   = integer format for annotation
!    sign_a = change sign of annotation (may be used to make
!                                               depths positive)
!    tick_len = length of ticks (default 0.15 cm)
!    xoff1     Unused ??
!    yoff1     Unused ??
!
      integer  :: iaxis, ntick, nannot
      real(DP) :: tbase, tincr
      integer,          optional :: sign_a
      character(len=*), optional :: rfmt, ifmt
      real(DP),         optional :: tick_len, xoff1, yoff1
!
      integer  :: i, imapp, imax, imin
      real(DP) :: x, y, dpc, xoff, yoff, xl, yl, tannot, xx, yy,
     &            tl, tl2
      character(len=8) :: string8
!
      dpc = 72d0/2.54d0
      xl = 0.15d0*dpc*abs(xbmax-xbmin)/(xcmax-xcmin)
      yl = 0.15d0*dpc*abs(ybmax-ybmin)/(ycmax-ycmin)
      tl = 0.15*dpc
      if(present(tick_len)) tl = abs(tick_len)*dpc
      tl2 = tl+tl
      xoff = tl2 + 0.1d0*dpc
      yoff = tl2 + 0.25d0*dpc
!      print *,xmax,xmin,xcmax,xcmin,xl
!      print *,ymax,ymin,ycmax,ycmin,xl
      if(tincr.eq.0d0)then
        print *,'error in calling routine ticks'
        print *,'  tincr equals zero'
        return
      endif
      if(xamin.ge.xamax)then
        print *," psplot routine ticks2 error:"
        print *,' a-space minimum greater or equal to maximum'
        return
      endif
!
!  save value of imap0 and turn off mapping by ps_plot
!
      imapp = imap0
      imap0 = 0
      tannot = abs(tincr*nannot)
!
!  x axis
!
      if(iaxis.eq.1)then
        imin = nint((xamin-tbase)/tincr)
        imax = nint((xamax-tbase)/tincr)
        do i=imin,imax
          x = tbase + i*tincr
          if(x.ge.xamin.and.x.le.xamax)then
            call pstran(x,yamin,xx,yy)
            if(ntick.ne.0.and.mod(i,max(ntick,1)).eq.0)then
              call ps_plot(xx,yy,0)
              call ps_plot(xx,yy-tl2,1)
              call stroke
            else
              call ps_plot(xx,yy,0)
              call ps_plot(xx,yy-tl,1)
              call stroke
            endif
            if(nannot.ne.0.and.mod(i,nannot).eq.0)then
              if(x.ne.0d0.and.present(sign_a))x = -x
              if(present(ifmt))then
                write(string8,ifmt)nint(x)
              elseif(present(rfmt))then
                write(string8,rfmt)x
              elseif(abs(x).gt.1d5.or.abs(x).lt.0.001d0)then
                write(string8,'(1pe8.1)')x
              elseif(tannot.ge.0.1d0)then
                write(string8,'(f8.1)')x
              elseif(tannot.ge.0.01d0)then
                write(string8,'(f8.2)')x
              else
                write(string8,'(f8.3)')x
              endif
              call pstext2(xx,yy-yoff,string8,1)
            endif
          endif
        end do
!
!  y axis
!
      elseif(iaxis.eq.2)then
        imin = nint((yamin-tbase)/tincr)
        imax = nint((yamax-tbase)/tincr)
        do i=imin,imax
          y = tbase + i*tincr
          if(y.ge.yamin.and.y.le.yamax)then
            call pstran(xamin,y,xx,yy)
            if(ntick.ne.0.and.mod(i,max(ntick,1)).eq.0)then
              call ps_plot(xx,yy,0)
              call ps_plot(xx-tl2,yy,1)
              call stroke
            else
              call ps_plot(xx,yy,0)
              call ps_plot(xx-tl,yy,1)
              call stroke
            endif
            if(nannot.ne.0.and.mod(i,nannot).eq.0)then
              if(y.ne.0d0.and.present(sign_a))y = -y
              if(present(ifmt))then
                write(string8,ifmt)nint(y)
              elseif(present(rfmt))then
                write(string8,rfmt)y
              elseif(abs(y).gt.1d5)then
                write(string8,'(1pe8.1)')y
              elseif(tannot.ge.0.1d0)then
                write(string8,'(f8.1)')y
              elseif(tannot.ge.0.01d0)then
                write(string8,'(f8.2)')y
              else
                write(string8,'(f8.3)')y
              endif
              call pstext2(xx-xoff,yy-0.1*dpc,string8,2)
            endif
          endif
        end do
      endif
!
!  secondry x axis
!
      if(iaxis.eq.3)then
        imin = nint((xamin-tbase)/tincr)
        imax = nint((xamax-tbase)/tincr)
        do i=imin,imax
          x = tbase + i*tincr
          if(x.ge.xamin.and.x.le.xamax)then
            call pstran(x,yamax,xx,yy)
            if(ntick.ne.0.and.mod(i,max(ntick,1)).eq.0)then
              call ps_plot(xx,yy,0)
              call ps_plot(xx,yy+tl2,1)
              call stroke
            else
              call ps_plot(xx,yy,0)
              call ps_plot(xx,yy+tl,1)
              call stroke
            endif
            if(nannot.ne.0.and.mod(i,nannot).eq.0)then
              if(x.ne.0d0.and.present(sign_a))x = -x
              if(present(ifmt))then
                write(string8,ifmt)nint(x)
              elseif(present(rfmt))then
                write(string8,rfmt)x
              elseif(abs(x).gt.1d5.or.abs(x).lt.0.001d0)then
                write(string8,'(1pe8.1)')x
              elseif(tannot.ge.0.1d0)then
                write(string8,'(f8.1)')x
              elseif(tannot.ge.0.01d0)then
                write(string8,'(f8.2)')x
              else
                write(string8,'(f8.3)')x
              endif
              call pstext2(xx,yy+yoff,string8,1)
            endif
          endif
        end do
!
!  secondry y axis
!
      elseif(iaxis.eq.4)then
        imin = nint((yamin-tbase)/tincr)
        imax = nint((yamax-tbase)/tincr)
        do i=imin,imax
          y = tbase + i*tincr
          if(y.ge.yamin.and.y.le.yamax)then
            call pstran(xamax,y,xx,yy)
            if(ntick.ne.0.and.mod(i,max(ntick,1)).eq.0)then
              call ps_plot(xx,yy,0)
              call ps_plot(xx+tl2,yy,1)
              call stroke
            else
              call ps_plot(xx,yy,0)
              call ps_plot(xx+tl,yy,1)
              call stroke
            endif
            if(nannot.ne.0.and.mod(i,nannot).eq.0)then
              if(y.ne.0d0.and.present(sign_a))y = -y
              if(present(ifmt))then
                write(string8,ifmt)nint(y)
              elseif(present(rfmt))then
                write(string8,rfmt)y
              elseif(abs(y).gt.1d5)then
                write(string8,'(1pe8.1)')y
              elseif(tannot.ge.0.1d0)then
                write(string8,'(f8.1)')y
              elseif(tannot.ge.0.01d0)then
                write(string8,'(f8.2)')y
              else
                write(string8,'(f8.3)')y
              endif
              call pstext2(xx-xoff,yy+0.15*dpc,string8,2)
            endif
          endif
        end do
      endif
!
!  reset value of imap0
!
      imap0 = imapp
!
      return
      end subroutine ticks2
!
!======================================================================
!
      subroutine ticks_y(iaxis, delta_year, delta_yeara,
     &                         delta_month, delta_montha, month_annot)
!
!  subroutine to draw tickmarks identifying years and months
!
!    iaxis        = axis to annotate 1 =bottom
!    delta_year   = interval between year ticks
!    delta_yeara  = interval between year annotation (0 = none)
!    depta_month  = interval between month ticks
!    delta_montha = interval between month annotation
!    month_annot  = Number of letters in month annotation (1 or 3)
!
      implicit none
      integer  :: iaxis, delta_year, delta_yeara,
     &            delta_month, delta_montha, month_annot
!
      integer  :: i, j, imapp, imax, imin
      real(DP) :: dpc, xl, yl, tl, tl2, xoff, yoff, x, xx, yy
      character(len=8) :: string8
      character(len=1) :: month1(12) = ["J","F","M","A","M","J",
     &                                  "J","A","S","O","N","D"]
      character(len=3) :: month3(12) = ["Jan","Feb","Mar","Apr",
     &                                  "May","Jun","Jul","Aug",
     &                                  "Sep","Oct","Nov","Dec"]
      integer  :: iyd(12) =
     &                [0,31,59,90,120,151,181,212,243,273,304,334]

!
      dpc = 72d0/2.54d0
      xl = 0.15d0*dpc*abs(xbmax-xbmin)/(xcmax-xcmin)
      yl = 0.15d0*dpc*abs(ybmax-ybmin)/(ycmax-ycmin)
      tl = 0.15*dpc
      tl2 = tl+tl
      xoff = tl2 + 0.1d0*dpc
      yoff = tl2 + 0.25d0*dpc

      if(delta_year.eq.0d0)then
        print *,'error in calling routine ticksy'
        print *,'  delta_year equals zero'
        return
      endif
      if(xamin.ge.xamax)then
        print *," psplot routine ticks2 error:"
        print *,' a-space minimum greater or equal to maximum'
        return
      endif
!
!  save value of imap0 and turn off mapping by ps_plot
!
      imapp = imap0
      imap0 = 0
!
!  x axis
!
      if(iaxis.eq.1)then
        imin = nint(xamin/delta_year)
        imax = nint(xamax/delta_year)
        do i=imin,imax
          x = i*delta_year
          print *,"AA",x,xamin,xamax, tl
          if(x.ge.xamin.and.x.le.xamax)then
            call pstran(x,yamin,xx,yy)
            print *,"  BB",x,yamin,xx,yy,yy-yl
            call ps_use_c_space
            call move2(xx,yy)
            call draw2(xx,yy-tl2)
            call ps_use_a_space
            call stroke
            if(delta_yeara.ne.0.and.mod(i,delta_yeara).eq.0)then
              write(string8,'(i4)')nint(x)
              call pstext2(xx,yy-1.5d0*yoff,string8,1)
            endif

!
!  Plot months
!
            if(delta_month.eq.1)then
              do j=1,11
                x = i*delta_year + iyd(j+1)/365d0
                if(x.lt.xamin.or.x.gt.xamax)exit
                call pstran(x,yamin,xx,yy)
                print *,"  BB",x,yamin,xx,yy,yy-yl
                call ps_use_c_space
                call move2(xx,yy)
                call draw2(xx,yy-tl)
                call ps_use_a_space
                call stroke
              enddo
              if(delta_montha.eq.1)then
                do j=1,12
                  x = i*delta_year + (j-0.5d0)/12d0
                  if(x.lt.xamin.or.x.gt.xamax)exit
                  call pstran(x,yamin,xx,yy)
                  if(month_annot.eq.3)then
                    call pstext2(xx,yy-yoff,month3(j),1)
                  elseif(month_annot.eq.1)then
                    call pstext2(xx,yy-yoff,month1(j),1)
                  endif
                enddo
              endif
            endif
          endif
        enddo
      endif
!
!  reset value of imap0
!
      imap0 = imapp
!
      return
      end subroutine ticks_y
!
!======================================================================
!
      subroutine axislimit(xx,xtick,nannot)
      implicit none
!
!  subroutine to choose tick and annotation limits for an axis.
!
!  input:  xx - a typical range to contain about 4 annotation marks
!               and 20 tick marks
!
!  output: xtick  - separation betweeen tick marks
!          nannot - number of tick marks per annotaion mark
!
!
      integer :: nannot
      real(DP)     xx, xtick
!
      integer :: n
      real(DP)     aa, x1, x2
!
      aa = log10(xx)
      n = aa
      if(aa.lt.0d0)n=n-1
      x1 = 10d0**n
      x2 = xx/x1
      if(x2.lt.1.5)then
        n=n-1
        x1 = x1/10d0
        x2 = x2*10d0
      endif
!
      if(x2.gt.6d0)then
        xtick = x1*0.5d0
        nannot = 4
      elseif(x2.gt.3d0)then
        xtick = x1*0.2d0
        nannot = 5
      else
        xtick = x1*0.1d0
        nannot = 5
      endif
      print 10,'axislimit ',xx,aa,x1,x2,xtick,n,nannot
  10  format(a,1p5e12.3,2i5)
      return
      end subroutine axislimit
!
!======================================================================
!  Define mapping
!======================================================================
!
      subroutine setmap(name)
      implicit none
!
!  routine to set mapping
!
      character*(*), intent(in) :: name
!
      imap1 = 0

      if(name(1:4).eq.'eqeq')then
        imap2 = 0
      elseif(name(1:8).eq.'mercator')then
        imap2 = 1
      elseif(name(1:5).eq.'ortho')then
        imap2 = 2
      else
        print *,' error calling setmap - name not recognised:'
        print *,  name(1:lnblnk(name))
        print *,' programme stopping ...'
      endif
!
      if(imap1.ne.0.or.imap2.ne.0)imap0=1
      return
      end subroutine setmap
!
!======================================================================
!
      subroutine setmap2(name1,name2)
      implicit none
!
!  routine to set double mapping
!  i.e. rotated lonlat grid to lonlat grid
!       lonlat grid to mercator
!
      character*(*),intent(in) :: name1,name2
!      integer :: iopost,ifontsize,imap0,imap1,imap2
!      logical lpsmask,lpvalid
!      common /psi/ iopost,ifontsize,imap0,imap1,imap2,lpsmask,lpvalid
!
      if(name1(1:4).eq.'eqeq')then
        imap1 = 0
      elseif(name1(1:5).eq.'euler')then
        imap1 = 1
      else
        print *,' error calling setmap2 - name1 not recognised:'
        print *,  name1(1:lnblnk(name1))
        print *,' programme stopping ...'
      endif
!
      if(name2(1:4).eq.'eqeq')then
        imap2 = 0
      elseif(name2(1:8).eq.'mercator')then
        imap2 = 1
      elseif(name2(1:5).eq.'ortho')then
        imap2 = 2
      else
        print *,' error calling setmap2 - name2 not recognised:'
        print *,  name2(1:lnblnk(name2))
        print *,' programme stopping ...'
      endif
!
      if(imap1.ne.0.or.imap2.ne.0)imap0=1
      return
      end subroutine setmap2
!
!======================================================================
!
      subroutine eqeq(x,y,i)
      implicit none
!
!-----------------------------------------------------------------------
!  eqeq mapping
!-----------------------------------------------------------------------
!
      real(DP)  :: x,y
      integer :: i
      return
      end subroutine eqeq
!
!======================================================================
!
      subroutine mercator(x,y,i)
      implicit none
!
!-----------------------------------------------------------------------
!  Mercator projection - angles in degrees
!-----------------------------------------------------------------------
!
      real(DP), intent(inout) :: x,y
      integer, intent(in)   :: i
!
      integer, save :: in = 0
      real(DP) :: y0
      real(DP), save :: pi,pid2,pid4,radian,degree
!      save in,pi,pid2,pid4,radian,degree
!
      if(in.eq.0)then
        in = 1
        pi = 4d0*atan(1d0)
        pid2 = pi/2d0
        pid4 = pi/4d0
        radian = pi/180d0
        degree = 180d0/pi
      endif
      y0 = y*radian
!
      if(i.ge.0)then
        if(abs(y0).gt.pid2)then
          y = sign(1d6,y0)
        elseif(y0.ge.0d0)then
          y = log(tan(pid4+y0*0.5d0))
        else
          y =-log(tan(pid4-y0*0.5d0))
        endif
      else
        if(y0.ge.0d0)then
          y = 2d0*atan(exp(y0))-pid2
        else
          y =-2d0*atan(exp(-y0))+pid2
        endif
      endif
      y = y*degree
      return
      end subroutine mercator
!
!======================================================================
!
      subroutine ortho(x,y,i)
!
!=======================================================================
!
!  Orthographic projection.  Convert latitude and longitude
!  coordinates into x and y coordinates on a flat plane through
!  the centre of the earth.
!
!  Note the earth radius is taken as equal to 1
!
!  Input
!    i  - +1 for forward transform, -1 for back transform
!    x  - longitude on forward transform (in degrees)
!    y  - latitude on forward transform  (in degrees)
!
!  Output
!    x  - x-coordinate on projected plane
!    y  - y coordinate on projected plane
!    ortho_z - z coordinate above projected plane.
!
!=======================================================================
!
      implicit none
      integer, intent(in)   :: i
      real(DP), intent(inout) :: x, y
!
      integer in
      real(DP)  aa(3), bb(3), xx, yy, radian, clon, slon, clat, slat
!      real(DP)  orthoa(3,3), ortho_lon, ortho_lat, ortho_z
!      common /orthoc/orthoa, ortho_lon, ortho_lat, ortho_z
      save    in, radian
      data    in/0/
!
      if(in.eq.0)then
        in = 1
        radian = atan(1d0)/45d0
      endif
      if(i.eq.1)then
        xx = x*radian
        yy = y*radian
        clon = cos(xx)
        slon = sin(xx)
        clat = cos(yy)
        slat = sin(yy)
!
        aa(1) = clat*slon
        aa(2) = slat
        aa(3) = clat*clon
!
        bb(1) = orthoa(1,1)*aa(1)+orthoa(1,2)*aa(2)+orthoa(1,3)*aa(3)
        bb(2) = orthoa(2,1)*aa(1)+orthoa(2,2)*aa(2)+orthoa(2,3)*aa(3)
        bb(3) = orthoa(3,1)*aa(1)+orthoa(3,2)*aa(2)+orthoa(3,3)*aa(3)
!
        x = bb(1)
        y = bb(2)
        ortho_z = bb(3)
      else
        print *
        print *,"  Subroutine ortho. Back transform not available"
        print *,"  Program stopping ..."
        print *
        stop
      endif
!
      return
      end subroutine ortho
!
!======================================================================
!
       subroutine set_ortho(rlon1,rlat1,ideg)
      implicit none
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!=======================================================================
!
!  Routine to initialise an orthographic projection of latitude
!  and longitude points.  This is one in which the projection lines are
!  parallel - thus it projects onto a plane through the Earth's centre
!  at right angles to the line of sight
!
!  This routine defines the line of sight as being along a line through
!  the point (rlong, rlat) and through the earths centre.  It
!  initialises a 3x3 matrix used to carry out the transformation.
!
!  Input:
!    rlon1  - longitude of reference point
!    rlat1  - latitude of reference point
!    ideg   - 1 if angles in degrees, 0 if in radians.
!
!=======================================================================
!
      integer ideg, in
      real(DP)  rlon1, rlat1
!
      real(DP)  radian, rlon, rlat, clon,slon,clat,slat
!      real(DP)  orthoa(3,3), ortho_lon, ortho_lat, ortho_z
      save    in,radian
!      common /orthoc/orthoa, ortho_lon, ortho_lat, ortho_z
      data    in/0/
!
      if(in.eq.0)then
        in = 1
        radian = atan(1d0)/45d0
      endif
      if(ideg.eq.1)then
        rlon = rlon1*radian
        rlat = rlat1*radian
      else
        rlon = rlon1
        rlat = rlat1
      endif
      ortho_lon = rlon
      ortho_lat = rlat
!
      clon = cos(rlon)
      slon = sin(rlon)
      clat = cos(rlat)
      slat = sin(rlat)
!
      orthoa(1,1) =  clon
      orthoa(1,2) =  0d0
      orthoa(1,3) = -slon
      orthoa(2,1) = -slat*slon
      orthoa(2,2) =  clat
      orthoa(2,3) = -slat*clon
      orthoa(3,1) =  clat*slon
      orthoa(3,2) =  slat
      orthoa(3,3) =  clat*clon
!
      return
      end subroutine set_ortho
!
!======================================================================
!
      subroutine seteu(alpha1,beta1,gamma1)
      implicit none
!
!  define euler angles
!
      real(DP), intent(in) :: alpha1,beta1,gamma1
!
      alpha = alpha1
      beta  = beta1
      gamma = gamma1
      return
      end subroutine seteu
!
!======================================================================
!
      subroutine gsave
      implicit none
!
      write(iopost,11)
  11  format('gsave')
      return
      end subroutine gsave
!
!======================================================================
!
      subroutine grestore
      implicit none

      write(iopost,12)
  12  format('grestore')
      return
      end subroutine grestore
!
!======================================================================
!
      subroutine grestoreall
      write(iopost,13)
  13  format('grestoreall')
      return
      end subroutine grestoreall
!
!======================================================================
!
      subroutine transeu(rlon1, rlat1, rlon2, rlat2, deg_in, deg_out,
     &                          alpha, beta,  gamma, deg_eu)
      implicit none
      real(DP)  rlon2, rlat2, rlon1, rlat1, alpha, beta, gamma
      integer deg_in, deg_out, deg_eu
!
!-----------------------------------------------------------------------
!23456789012345678901234567890123456789012345678901234567890123456789012
!
! PURPOSE To convert latitude and longitude in one spherical coordinate
!         system to latitude and longitude in a second system.  The
!         angles alpha, beta and gamma are the three Euler angles
!         defining the transform from the first system to the second.
!         The rotation is defined here is the one commonly used with
!         spherical harmonics:
!           1.  Rotate meridian by angle alpha around the pole.
!           2.  Rote pole down angle beta along the meridian.
!           3.  Rotate new meridian by angle gamma around the new pole.
!         Another commonly used definition, used in dynamics, is:
!           1.  Rotate meridian by angle alpha' around pole.
!           2.  Rote pole down angle beta', keeping the point where
!               the meridian crosses the equator stationary.  (This
!               direction defines the astronomical node).
!           3.  Rotate new meridian by angle gamma' around the new pole.
!          The angles are related by the equations:
!               alpha' = alpha + 90 degrees
!               beta'  = beta
!               gamma' = gamma - 90 degrees
!
!           The back transformation corresponding to (alpha,beta,gamma)
!           is (-gamma,-beta,-alpha).
!
! INPUT     rlat1  - model 1 latitude
!           rlon1  - model 1 longitude
!           deg_in - 1 if input values in degrees
!                    0 if input values in radians
!           deg_out- 1 if output values in degrees
!                    0 if output values in radians
!
!           alpha  - first Euler angle
!           beta   - second Euler angle (latitude)
!           gamma  - third Euler angle
!           deg_eu - 1 if Euler angles in degrees
!                    0 if Euler angles in radians
!
! OUTPUT    rlat2  - model 2 latitude
!           rlon2  - model 2 longitude
!
!  OCCAM    For the forward transform from the model 1 grid
!           to the model 2 grid, OCCAM uses:
!           alpha = 232 degrees
!           beta  =  90 degrees
!           gamma =   0 degrees
!
!  Last modified by djw, 15 June 2010
!-----------------------------------------------------------------------
!
      integer in
      real(DP)  a, b, c, p1, q1, p2, q2, s, t
      real(DP)  pi, twopi, radian, degree
      save in, pi,twopi,degree,radian
      data in/0/
!
!-----------------------------------------------------------------------
!  constants
!-----------------------------------------------------------------------
!
      if(in.eq.0)then
!        pi = 3.1415926535897932d0
        pi = 4d0*atan(1d0)
        twopi = 2d0*pi
        radian= pi/180d0
        degree= 180d0/pi
        in = 1
      endif
!
!-----------------------------------------------------------------------
!  copy across arguments and scale if necessary
!-----------------------------------------------------------------------
!
      if(deg_eu.ne.0)then
        a = alpha*radian
        b = beta*radian
        c = gamma*radian
      else
        a = alpha
        b = beta
        c = gamma
      endif
!
      if(deg_in.ne.0)then
        p1 = radian*(90d0-rlat1)
        q1 = radian*(rlon1-alpha)
      else
        p1 = 0.5d0*pi-rlat1
        q1 = rlon1-alpha
      endif
!
!  transform q1 to range -pi to +pi
!
      q1 = mod(q1,twopi)
      if(q1.gt.pi) q1=q1-twopi
      if(q1.le.-pi)q1=q1+twopi
!      print *
!      print 11,a,b,c
!      print 11,p1,q1
!  11  format(5x,5f12.6)
!
      p2 = acos(cos(b)*cos(p1)+sin(b)*sin(p1)*cos(q1))
!
      if(b.eq.0)then
        q2 = q1
      elseif(p2.eq.0d0)then
        q2 = 0d0
      else
        s = (p1+p2+b)*0.5d0
        t = sin(s-b)*sin(s-p2)/(sin(b)*sin(p2))
!  correct for small rounding errors near the limits
        if(t.lt.0d0)t = 0d0
        if(t.gt.1d0)t = 1d0
        q2 = 2d0*acos(sqrt(t))
        if(q1.lt.0d0)q2 = -q2
!
!        if(isnan2(q2))then
!          print *,"File transeu.F.  Error in routine transeu:"
!          print *," variable a = ",a
!          print *," variable b = ",b
!          print *," variable c = ",c
!          print *," variable p1 = ",q2
!          print *," variable q1 = ",q1
!          print *," variable s = ",s
!          print *," variable t = ",t
!          print *," variable p2 = ",p2
!          print *," variable q2 = ",q2
!        endif
!
!        print 11,s,t,q2
      endif
!
!-----------------------------------------------------------------------
!  rescale the output value then and return to a float
!  from the double values used in the conversion
!  program
!-----------------------------------------------------------------------
!
      if( deg_out.ne.0)then
        rlon2 = degree*(q2 - c)
        rlon2 = mod(rlon2,360d0)
        rlat2 = 90d0-degree*p2
      else
        rlon2 = q2 - c
        rlon2 = mod(rlon2,twopi)
        rlat2 = 0.5d0*pi-p2
      endif
      return
      end subroutine transeu
!
!======================================================================
!  Predefined spectra of colours
!======================================================================
!
      subroutine col_steps1(x,xmin,xmax,ncols,r,g,b,icol,ipal)
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  Routine to transform 'x' to red, green and blue colours from a
!  spectrum of 'ncols' colours  :  colour bar
!  Input:
!    x   -  input variable
!    xmin - minimum value (at start of first colour step)
!    xmax - maximum value (at end of last colour step)
!    ncols - number of colour steps
!  Output
!    r    - output red colour (range 0d0 to 1d0)
!    g    - output green colour (range 0d0 to 1d0)
!    b    - output blue colour (range 0d0 to 1d0)
!    icol - index of final colour in spectrum
!    ipal - (optional) pallet  1 = standard
!                              2 = red    - blue
!                              3 = yellow - green
!            negative values reverse colour sequence
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      implicit none
      integer, parameter :: DP = kind(1.0d0)
!
      real(DP),intent(in)  :: x, xmin, xmax
      integer, intent(in)  :: ncols
      real(DP),intent(out) :: r, g, b
      integer, intent(out) :: icol
      integer, intent(in), optional :: ipal
!
      real(DP) :: z
      integer  :: im, ip, kpal, i, j
      integer  :: ia(3)  = (/6, 5, 5/)
      real(DP) :: za(6,3), ra(6,3), ga(6,3), ba(6,3)
!
!  Rows 1, 2 and 3 of arrays corrrespond to ipal values 1, 2 and 3
!
      za = reshape((/0.00d0,0.10d0,0.34d0,0.67d0,0.90d0,1.00d0,
     &               0.00d0,0.15d0,0.50d0,0.85d0,1.00d0,0.00d0,
     &               0.00d0,0.25d0,0.50d0,0.75d0,1.00d0,0.00d0/),
     &               shape(za))
      ra = reshape((/0.80d0,0.30d0,0.00d0,1.00d0,1.00d0,0.60d0,
     &               0.00d0,0.00d0,0.99d0,1.00d0,0.60d0,0.00d0,
     &               0.00d0,0.00d0,0.80d0,1.00d0,0.80d0,0.00d0/),
     &               shape(ra))
      ga = reshape((/0.00d0,0.20d0,0.60d0,1.00d0,0.00d0,0.00d0,
     &               0.00d0,0.00d0,0.80d0,0.00d0,0.00d0,0.00d0,
     &               0.60d0,1.00d0,0.80d0,1.00d0,0.60d0,0.00d0/),
     &               shape(ga))
      ba = reshape((/1.00d0,1.00d0,0.30d0,0.00d0,0.00d0,0.00d0,
     &               0.60d0,1.00d0,0.99d0,0.00d0,0.00d0,0.00d0,
     &               0.70d0,0.00d0,0.80d0,0.00d0,0.00d0,0.00d0/),
     &               shape(ba))
!
      if(ncols.lt.3.or.xmin.eq.xmax)goto 999
      z = (x-xmin)/(xmax-xmin)
      if(z.lt.0d0.or.z.gt.1d0)goto 999
!
      kpal = 1
      if(present(ipal))then
        if(ipal.gt.3.or.ipal.lt.-3.or.ipal.eq.0)goto 999
        kpal = abs(ipal)
        if(ipal.lt.0)z=1d0-z
      endif
!
      z = int(z*ncols)
      if(z.gt.ncols-1)z=ncols-1
      icol = nint(z)+1
      z = z/(ncols-1)

      im = 1
      do i=2,ia(kpal)-1
        if(z.gt.za(i,kpal))im = i
      enddo
      ip = im+1
!
      z = (z-za(im,kpal))/(za(ip,kpal)-za(im,kpal))
      r = ra(im,kpal)+z*(ra(ip,kpal)-ra(im,kpal))
      g = ga(im,kpal)+z*(ga(ip,kpal)-ga(im,kpal))
      b = ba(im,kpal)+z*(ba(ip,kpal)-ba(im,kpal))
      return
!
999   r = 0d0 ; g = 0d0 ; b=0d0 ; icol=0
      return
      end subroutine col_steps1
!
!======================================================================
!
      subroutine col_rd_bu_11(x,xmin,xmax,nstep,r,g,b,ipal)
!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  Routine to transform 'x' into red, green and blue colours
!  using The Colorbrewer 'RdBu' spectrum.
!  Input:
!    x   -  input variable
!    xmin - minimum value (at start of first colour step)
!    xmax - maximum value (at end of last colour step)
!    nstep- number of steps (0 = smooth interpolation)
!    ipal - (optional) negative value reverses colour range.
!
!  Output:
!    r    - output red colour   (range 0d0 to 1d0)
!    g    - output green colour (range 0d0 to 1d0)
!    b    - output blue colour  (range 0d0 to 1d0)
!=======================================================================
!
      implicit none
!      integer, parameter :: DP = kind(1.0d0)
!
      real(DP),intent(in)  :: x, xmin, xmax
      integer, intent(in)  :: nstep
      integer, intent(in), optional :: ipal
      real(DP),intent(out) :: r, g, b
!
      real(DP) :: z
      integer  :: iz
!
      real(DP),parameter :: rd_bu_11(3,11)=  reshape(
     &                                              [103,    0,   31,
     &                                               178,   24,   43,
     &                                               214,   96,   77,
     &                                               244,  165,  130,
     &                                               253,  219,  199,
     &                                               247,  247,  247,
     &                                               209,  229,  240,
     &                                               146,  197,  222,
     &                                                67,  147,  195,
     &                                                33,  102,  172,
     &                                                 5,   48,   97],
     &                                                 [3,11])


      if(xmin.eq.xmax)then
        print *," Subroutine col_rd_bu_11 error"
        print *," xmin equals xmax"
        goto 999
      endif
      z = (x-xmin)/(xmax-xmin)

      if(present(ipal))then
        if(ipal.lt.0)z = 1d0-z
      endif

      if(z.lt.0d0.or.z.gt.1d0)then
        r = 0d0 ; g = 0d0; b = 0d0
        return
      endif

      if(nstep.ne.0)then
       z = int(z*nstep)
       z = z/(nstep-1)
      endif

      z  = z*10d0
      iz = min(int(z),9)
      z  = z-iz
      iz = iz+1
      r = (rd_bu_11(1,iz) + z*(rd_bu_11(1,iz+1)-rd_bu_11(1,iz)))/255d0
      g = (rd_bu_11(2,iz) + z*(rd_bu_11(2,iz+1)-rd_bu_11(2,iz)))/255d0
      b = (rd_bu_11(3,iz) + z*(rd_bu_11(3,iz+1)-rd_bu_11(3,iz)))/255d0

      return

999   print *,"Program stopping ..."
      stop

      end subroutine col_rd_bu_11
!
!======================================================================
!
      subroutine col_spec_1(x,xmin,xmax,nstep,r,g,b,ipal)
!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  Routine to transform 'x' into red, green and blue colours.
!  Input:
!    x   -  input variable
!    xmin - minimum value (at start of first colour step)
!    xmax - maximum value (at end of last colour step)
!    nstep- number of steps (0 = smooth interpolation)
!    ipal - (optional) negative value reverses colour range.
!
!  Output:
!    r    - output red colour   (range 0d0 to 1d0)
!    g    - output green colour (range 0d0 to 1d0)
!    b    - output blue colour  (range 0d0 to 1d0)
!=======================================================================
!
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer, parameter :: ncol = 7
!
      real(DP),intent(in)  :: x, xmin, xmax
      integer, intent(in)  :: nstep
      integer, intent(in), optional :: ipal
      real(DP),intent(out) :: r, g, b
!
      real(DP) :: z
      integer  :: iz
!
      real(DP),parameter :: rd_bu_11(3,ncol)=  reshape(
     &                                          [204,    0,  255,  ! lt purple
     &                                           076,  051,  255,  ! dk blue
     &                                           150,  201,  255,  ! lt blue
     &                                           000,  153,  030,  ! dk green
     &                                           255,  255,  000,  ! yellow
     &                                           255,  000,  000,  ! red
     &                                           153,  000,  000], ! dk red
     &                                             [3,ncol])


      if(xmin.eq.xmax)then
        print *," Subroutine col_rd_bu_11 error"
        print *," xmin equals xmax"
        goto 999
      endif
      z = (x-xmin)/(xmax-xmin)

      if(present(ipal))then
        if(ipal.lt.0)z = 1d0-z
      endif

      if(z.lt.0d0.or.z.gt.1d0)then
        r = 0d0 ; g = 0d0; b = 0d0
        return
      endif

      if(nstep.ne.0)then
       z = int(z*nstep)
       z = z/(nstep-1)
       if(z.gt.1d0)z=1d0
      endif

      z  = z*(ncol-1)
      iz = min(int(z),ncol-2)
      z  = z-iz
      iz = iz+1
      r = (rd_bu_11(1,iz) + z*(rd_bu_11(1,iz+1)-rd_bu_11(1,iz)))/255d0
      g = (rd_bu_11(2,iz) + z*(rd_bu_11(2,iz+1)-rd_bu_11(2,iz)))/255d0
      b = (rd_bu_11(3,iz) + z*(rd_bu_11(3,iz+1)-rd_bu_11(3,iz)))/255d0

      return

999   print *,"Program stopping ..."
      stop
      end subroutine col_spec_1
!
!======================================================================
!  Miscellaneous
!======================================================================
!
      subroutine ps_push()
      type(ps_node),pointer :: new_node
!
!      print *," Enter ps_push"

      allocate(new_node)
      if(associated(curr_node))new_node%prev_node => curr_node
      curr_node => new_node
!
      curr_node%xamin = xamin ; curr_node%xamax = xamax
      curr_node%yamin = yamin ; curr_node%yamax = yamax
      curr_node%xbmin = xbmin ; curr_node%xbmax = xbmax
      curr_node%ybmin = ybmin ; curr_node%ybmax = ybmax
      curr_node%xcmin = xcmin ; curr_node%xcmax = xcmax
      curr_node%ycmin = ycmin ; curr_node%ycmax = ycmax
      curr_node%xmmin = xmmin ; curr_node%xmmax = xmmax
      curr_node%ymmin = ymmin ; curr_node%ymmax = ymmax
!      print *," Exit  ps_push"
      return
      end subroutine
!
!======================================================================
!
      subroutine ps_pop()
      type(ps_node), pointer :: prev_node
!      print *," Enter ps_pop"
      if(.not.associated(curr_node))then
!        print *," curr_node not associated"
        return
      endif
!
      xamin = curr_node%xamin
      xamin = curr_node%xamin ; xamax = curr_node%xamax
      yamin = curr_node%yamin ; yamax = curr_node%yamax
      xbmin = curr_node%xbmin ; xbmax = curr_node%xbmax
      ybmin = curr_node%ybmin ; ybmax = curr_node%ybmax
      xcmin = curr_node%xcmin ; xcmax = curr_node%xcmax
      ycmin = curr_node%ycmin ; ycmax = curr_node%ycmax
      xmmin = curr_node%xmmin ; xmmax = curr_node%xmmax
      ymmin = curr_node%ymmin ; ymmax = curr_node%ymmax
      if(associated(curr_node%prev_node))then
!        print *," set prev_node to curr_node"
        prev_node = curr_node%prev_node
      else
!        print *," set prev_node to null"
        prev_node => null()
      endif
!      print *," deallocate curr_node"
      deallocate(curr_node)
!      print *,"  curr_node deallocated"
      if(associated(prev_node))then
        curr_node = prev_node
      else
        nullify(curr_node)
      endif
!      print *," Exit  ps_pop"
      return
      end subroutine
!
!======================================================================
!
      subroutine contc(ic,nxx,nx,ny,xe,ye,nc,rc,gc,bc)
      implicit none
!
!  subroutine to fill rectangular grid with colour using scheme to
!    reduce number of postscript commands
!
!    ic  :  array defining colours in each grid square
!    nxx :  first dimension of ic
!    nx  :  number of grids in x direction
!    ny  :  number of grids in y direction
!    xe  :  x-coordinates of grid edges (1:nx+1)
!    ye  :  y-coordinates of grid edges (1:ny+1)
!    nc  :  number of colours
!    rc  :  red array
!    gc  :  green array
!    bc  :  blue array
!
!    Warning:  The grouping algorithm works well in rectangular spaces
!    but will generate errors when used with a curved space.
!======================================================================
!
      integer, intent(in) :: ic(nxx,ny), nxx, nx, ny, nc
      real(DP), intent(in) :: xe(nx+1),ye(ny+1),rc(nc),gc(nc),bc(nc)

      integer  :: i, j, k, kl
      real(DP) :: r, g, b, rl, gl, bl
      character(len=*),parameter :: myname="contc"
!
      do j=1,ny
        r = -1d0 ; g = -1d0 ; b = -1d0 ; k = 0
        do i=1,nx
          kl = k ; rl = r ; gl = g ; bl = b
          k = ic(i,j)
          if(k.gt.nc.or.k.lt.0)then
            print *," Routine "//myname
            print *,"  Colour index error for i,j = ",i,j
            print *,"  Colour index = ",k
            print *,"  Allowed range is 0 : ",nc
            stop
          endif
          if(k.eq.0)then
            r = -1d0 ; g = -1d0 ; b = -1d0
          else
            r = rc(k); g = gc(k); b = bc(k)
          endif
          if(kl.ne.0)then
            if(i.eq.nx.or.r.ne.rl.or.g.ne.gl.or.b.ne.bl)then
              call psplot(xe(i+1),ye(j+1),1)
              call psplot(xe(i+1),ye(j),  1)
              call closepath
              call fill
            endif
          endif
          if(k.ne.0)then
            if(i.eq.1.or.r.ne.rl.or.g.ne.gl.or.b.ne.bl)then
              call setrgb(rc(k),gc(k),bc(k))
              call psplot(xe(i),ye(j),  0)
              call psplot(xe(i),ye(j+1),1)
            endif
          endif
        enddo
      enddo
      call setgrey(0d0)
      return
      end subroutine contc
!
!======================================================================
!
      subroutine pscomment(string,i1,i2)
      implicit none
!
!  subroutine to add comment line to postscript output
!
!  if i1 (or i2) is non zero, blank comment line is also added
!  before  (or after) the comment line
!
      integer i1, i2
      character*(*) string
!
!      integer iopost,ifontsize,imap0,imap1,imap2
!      logical lpsmask,lpvalid
!      common /psi/ iopost,ifontsize,imap0,imap1,imap2,lpsmask,lpvalid
!
      if(i1.ne.0)write(iopost,11)
      write(iopost,12)string
      if(i2.ne.0)write(iopost,11)
!
  11  format("%")
  12  format("% ",a)
      return
      end subroutine pscomment
!
!======================================================================
!  Private subroutines
!======================================================================
!
      subroutine pssubs_print_on_nan(x,var_name)
      implicit none
!
!  Subroutine to print values even if they are NaN values.
!  Used to debug source of NaN values.
!
!  Input:
!    x           value of variable
!    var name    string name of variable
!    var_format  format to use when printing value
!
      real(DP)    x
      character*(*) var_name
!
      if(isnan2(x))then
        print 11,var_name
      else
        print *,"    variable ",var_name," = ",x
      endif
      return
 11   format("    variable ",a," = NaN.")
!12   format("    variable ",a," = ",f10.3,".")

      end subroutine pssubs_print_on_nan
!
!======================================================================
!
      subroutine ps_plot(x,y,ii)
      implicit none
!
!   if ii == 0 then call newpath and moveto
!         else      call linetp
!   if il == 0 then the last plotted point (xl,yl) was valid
!
!  Common block /psi/ :
!
!   iopost    - fortran stream number of postscript output file
!   ifontsize - font size in points.  (See setfont)
!   imap1     - current mapping.      (See setmap)
!   lpsmask   - masking on or off.    (See psmask)
!   lpvalid   - true if last point valid (i.e. we have started to
!                       draw a line)
!
!  Limits:
!
!    xamin, ...   - limits in user units (i.e. degrees)
!    xbmin, ...   - limits in transformed units (i.e. after
!                            calling the mercator mapping function)
!    xcmin, ...   - limits of plotting region on page (pixels).
!    xmmin, ...   - limits of masking region on page (pixels).
!
      real(DP), intent(in)  :: x, y
      integer, intent(in) :: ii
!
!      real(DP)       xamin,xamax,yamin,yamax,xbmin,xbmax,ybmin,ybmax,
!     &             xcmin,xcmax,ycmin,ycmax,xmmin,xmmax,ymmin,ymmax
!      common /psr/ xamin,xamax,yamin,yamax,xbmin,xbmax,ybmin,ybmax,
!     &             xcmin,xcmax,ycmin,ycmax,xmmin,xmmax,ymmin,ymmax
!      integer :: iopost,ifontsize,imap0,imap1,imap2
!      logical lpsmask, lpvalid
!      common /psi/ iopost,ifontsize,imap0,imap1,imap2,lpsmask,lpvalid
!
      logical valid,lprint
      integer :: nerrs
      real(DP)     xl,yl,xx,yy,xxx,yyy
      save xl,yl,nerrs
      data xl,yl,nerrs /0.0,0.0,0/
!
!  call mapping if imap0 set to 1, skip if set to 0
!  to allow mix of calls from user 'a' space or page 'c' space.
!
      if(imap0.eq.0)then
        xx = x
        yy = y
      else
        call pstran(x,y,xx,yy)
      endif
!
!  check validity
!
      valid = xx.ge.xmmin.and.xx.le.xmmax
     &                       .and.yy.ge.ymmin.and.yy.le.ymmax
!
! DEBUG!
!
      lprint = .false.
      if(lprint)then
        print 51,x,y,ii,valid,lpsmask,lpvalid
        print 52,xx,yy,xmmin,xmmax,ymmin,ymmax
  51  format('ps_plot: ',2f10.3,i5,4l5)
  52  format('        ',6f10.3)
      endif
!
!  last point valid and this point valid
!
      if(.not.lpsmask.or.((lpvalid.or.ii.eq.0).and.valid))then
        if(lprint)print *,' type 1'
        if(ii.eq.0)then
          if(lpvalid)call stroke
          write(iopost,11)
          write(iopost,12) xx,yy
          if(lprint)then
            print 11
            print *, xx,yy,' -aa- ',lpsmask,ii,valid
          endif
        else
          write(iopost,21) xx,yy
        endif
        lpvalid = .true.
        xl = xx
        yl = yy
        return
!
!  last point not valid, this valid
!
      elseif(valid)then
        if(yl.lt.ymmin)then
          if(lprint)print *,' type 2'
          xxx = xl +(ymmin-yl)*(xx-xl)/(yy-yl)
          if(xxx.ge.xmmin.and.xxx.le.xmmax)then
            write(iopost,11)
            write(iopost,12) xxx,ymmin
            if(lprint)then
              print 11
              print *, xxx,ymmin,' -bb- ',lpsmask,ii,valid
            endif
            write(iopost,21) xx,yy
            lpvalid = .true.
            xl = xx
            yl = yy
            return
          endif
        elseif(yl.gt.ymmax)then
          if(lprint)print *,' type 3'
          xxx = xl +(ymmax-yl)*(xx-xl)/(yy-yl)
          if(xxx.ge.xmmin.and.xxx.le.xmmax)then
            write(iopost,11)
            write(iopost,12) xxx,ymmax
            write(iopost,21) xx,yy
            if(lprint)then
              print 11
              print *, xxx,ymmax,' -cc- ',lpsmask,ii,valid
            endif
            lpvalid = .true.
            xl = xx
            yl = yy
            return
          elseif(lprint)then
            print 121,xx,yy,xl,yl,xmmin,xmmax,xxx
 121        format('type 3:',8f10.3)
          endif
        endif
        if(xl.lt.xmmin)then
          if(lprint)print *,' type 4'
          yyy = yl +(xmmin-xl)*(yy-yl)/(xx-xl)
          if(yyy.ge.ymmin.and.yyy.le.ymmax)then
            write(iopost,11)
            write(iopost,12) xmmin,yyy
            if(lprint)then
              print 11
              print *, xmmin,yyy,' -dd- ',lpsmask,ii,valid
            endif
            write(iopost,21) xx,yy
            lpvalid = .true.
            xl = xx
            yl = yy
            return
          endif
        elseif(xl.gt.xmmax)then
          if(lprint)print *,' type 5'
          yyy = yl +(xmmax-xl)*(yy-yl)/(xx-xl)
          if(yyy.ge.ymmin.and.yyy.le.ymmax)then
            write(iopost,11)
            write(iopost,12) xmmax,yyy
            if(lprint)then
              print 11
              print 12, xmmax,yyy,' -aa- ',lpsmask,ii,valid
            endif
            write(iopost,21) xx,yy
            lpvalid = .true.
            xl = xx
            yl = yy
            return
          endif
        endif
!
!  last point valid, this not valid
!
      elseif(lpvalid)then
        if(yy.lt.ymmin)then
          if(lprint)print *,' type 6'
          xxx = xl +(ymmin-yl)*(xx-xl)/(yy-yl)
          if(xxx.ge.xmmin.and.xxx.le.xmmax)then
            write(iopost,21) xxx,ymmin
            call stroke
            xl = xx
            yl = yy
            return
          endif
        elseif(yy.gt.ymmax)then
          if(lprint)print *,' type 7'
          xxx = xl +(ymmax-yl)*(xx-xl)/(yy-yl)
          if(xxx.ge.xmmin.and.xxx.le.xmmax)then
            write(iopost,21) xxx,ymmax
            call stroke
            xl = xx
            yl = yy
            return
          endif
        endif
        if(xx.lt.xmmin)then
          if(lprint)print *,' type 8'
          yyy = yl +(xmmin-xl)*(yy-yl)/(xx-xl)
          if(yyy.ge.ymmin.and.yyy.le.ymmax)then
            write(iopost,21) xmmin,yyy
            call stroke
            xl = xx
            yl = yy
            return
          endif
        elseif(xx.gt.xmmax)then
          if(lprint)print *,' type 9'
          yyy = yl +(xmmax-xl)*(yy-yl)/(xx-xl)
          if(yyy.ge.ymmin.and.yyy.le.ymmax)then
            write(iopost,21) xmmax,yyy
            call stroke
            xl = xx
            yl = yy
            return
          endif
        endif
      endif
!
! neither point valid
!
      if(lprint)print *,' type 10'
      xl = xx
      yl = yy
      return
  11  format('n')
  12  format(f0.3,1x,f0.3,' m')
  21  format(f0.3,1x,f0.3,' l')
      end subroutine ps_plot
!
!======================================================================
!
      logical function isnan2(x)
      implicit none
!
!  function to find NaN values
!
      real(DP)  x
!
      if(x.lt.1d70.and.x.gt.-1d70)then
        isnan2 = .false.
      else
        isnan2 = .true.
      endif
      return
      end function isnan2

      end module
