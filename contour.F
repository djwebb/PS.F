!
!=======================================================================
!
!    File: contour.F
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
!
      subroutine conta(z,aa,ia,iadim,jadim,ino,jno,zx,zy,
     &                 annd,annx,anny,anna,annv,nann,mann)
!
!=======================================================================
!
!  PS.F contouring routine
!
!  Input:
!  z      - value to contour
!  aa     - data array
!  ia     - work array.  On entry set to 0 for data points
!                                        1 for points where aa undefined
!           within the routine ia(i,j) = 2 if aa(i,j).le.z
!                                      = 3 if aa(i,j).gt.z
!                                      = 4 after use for plotting
!  iadim,jadim - dimensions of aa and ia
!  ino, jno  - dimensions of region to contour (1:ino,1,jn0)
!  zx,  zy   - 1-D arrays containing coordinate of data point
!  annd      - nominal distance between contour annotations
!              distance in points (1/72 inch)
!  nann      - number of previous annotations - add new ones after this
!  mann      - maximum number of annotations to add (0=no annotations)
!
!  Output:
!  annx,anny - array page coordinate of annotation
!  anna      - angle to annotate
!  annv      - value to annotate (equals z)
!  nann      - total number of annotations
!
      use ps
      implicit none
      integer, parameter :: DP = kind(1.0d0)
!
      integer :: ia(iadim,jadim)
      real(DP)   z, aa(iadim,jadim),zx(iadim),zy(jadim)
      real(DP)   annd, annx(mann), anny(mann), anna(mann), annv(mann)
      real(DP)   bnnd, bnnx(24), bnny(24), bnna(24)
      real(DP)   xa1, ya1, sa1, xa2, ya2, sa2
      integer :: iadim, jadim, ino, jno, nann, mann, ncnn
      logical :: lann
!
      integer :: i, j, k, l, m, n, ii, jj, kk, kl, kr, ica, icu, iaa,
     &           n1, la, lb, lc, it, id, jd, i5, mst, nst, kst,
     &           i0t, i1t, ic,ml, nl, kc
      real(DP)   pi, degree, al, ar, ab, abmax, abmin, ac,
     &           bl, br,x,y,xx,yy, zc, zd
!
      integer :: land, less, more, ichk
      integer :: ida(5),jda(5),ka(5),is(4),js(4)
      Integer :: kla(12),kra(12),kca(12),ma(4),na(4),kka(4),kkr(4)
      real(DP) :: xa(5),ya(5)
!
      data land,less,more,ichk/1,2,3,4/
      data ida/1,0,-1,0,1/
      data jda/0,1,0,-1,0/
      data ka/1,2,3,4,1/
      data is/0,1,1,0/
      data js/0,0,1,1/
      data xa/0d0,1d0,1d0,0d0,0.5d0/
      data ya/0d0,0d0,1d0,1d0,0.5d0/
      data kla/12,6,8,10,7,1,9,2,11,3,5,4/
      data kra/5,7,9,11,2,12,3,6,4,8,1,10/
      data kca/0,0,0,0,3,1,4,2,1,3,2,4/
      data ma/0,1,0,-1/
      data na/-1,0,1,0/
      data kka/3,4,1,2/
      data kkr/2,3,4,1/
! DEBUG
      integer icount
      data icount /0/
!
!      print *," Enter conta",iadim,jadim,ino,jno
!
!  initialise annotations
!
      lann = mann.gt.0.and.nann.lt.mann
      bnnd = 2d0
      if(lann)then
        pi = 4d0*atan(1d0)
        degree = 180d0/pi
        ncnn = 0
        bnnd = annd/3d0
        xa1 = 0d0
        ya1 = 0d0
        sa1 = 0d0
        xa2 = 0d0
        ya2 = 0d0
        sa2 = 0d0
      endif
!
!  initialise index array ia
!
      ica = jno+ino
      icu = jno
      n1 = 1
      do 30 j=1,jno
      do 30 i=1,ino
      if(ia(i,j).ne.land)then
        if(aa(i,j).le.z)then
          ia(i,j) = less
        else
          ia(i,j) = more
        endif
      endif
  30  continue
!      print *
!      do j=jno,1,-1
!        print 31,(ia(i,j),i=1,ino)
!  31  format(i5,1x,100i1,(6x,101i1))
!      enddo
!
!  search outside of grid, then inside grid for next contour
!
      i = 1
      j = 1
      la=ia(1,1)
!
      do 110 it = 1,5
      id = ida(it)
      jd = jda(it)
      kk = ka(it)
      icu = ica-icu
      if(it.ne.5) goto 60
  40  i = 1
      j = j+1
      if(j.gt.jno) goto 110
      la = ia(i,j)
!
  60  do 100 ic = 2,icu
      i = i+id
      j = j+jd
      lb = la
      la = ia(i,j)
!  first test points indirection of search
      if(la.eq.less.and.lb.eq.more)goto 200
      if(la.eq.more.and.lb.eq.less)goto 200
      if(it.lt.5)goto 100
!  if unsuccessful and within the grid, test in negative j direction
!  this is to pick up contours which do not cross lines of constant j
      lc = ia(i-1,j-1)
      if(lb.eq.less.and.lc.eq.more)goto 150
      if(lb.eq.more.and.lc.eq.less)goto 150
      goto 100
!
!  start of contour which does not cross lines of constant j
!
 150  i5 = 1
      m = i-1
      n = j-1
      k = 4
      kl = 4
      kr = 1
      al = aa(m,j)
      ar = aa(m,n)
      ia(m,j) = ichk
      ia(m,n) = ichk
      mst = m
      nst = n
      kst = k
      goto 205
!
!  start all other contours
!
 200  i5 = 0
      k = kk
      kl = k
      kr = kkr(k)
      ii = i-id
      jj = j-jd
      al = aa(ii,jj)
      ar = aa(i,j)
      ia(ii,jj) = ichk
      ia(i,j)   = ichk
      m = ii-is(k)
      n = jj-js(k)
      mst = m
      nst = n
      kst = k
!
!  start drawing contour
!
 205  zd = ar-al
      if(zd.ne.0.0)zd=1.0/zd
      bl = zd*(z-al)
      br = zd*(ar-z)
      x = zx(m)+(zx(m+1)-zx(m))*(br*xa(kl)+bl*xa(kr))
      y = zy(n)+(zy(n+1)-zy(n))*(br*ya(kl)+bl*ya(kr))
!
      icount = icount+1
!      if(icount.lt.30) print *,'moveto ',m,n,ia(m,n),x,y
      call move2(x,y)
!
!  annotation logic at start of line
!
      if(lann)then
        call pstran(x,y,xx,yy)
        xa1=xx
        ya1=yy
        sa1=0d0
        ncnn = 0
      endif
!
!  calculate value of function at central point of next square of
!  data.  end contour if next square contains a land point where
!  function is undefined
!
 210  if(ia(m,n).eq.land)goto 410
      ac = aa(m,n)
      abmin = ac
      abmax = ac
      do 250 l=2,4
      ml = m+is(l)
      nl = n+js(l)
      if(ia(ml,nl).eq.land)goto 410
      ab = aa(ml,nl)
      ac = ac + ab
      abmax = max(ab,abmax)
      abmin = min(ab,abmin)
 250  continue
      ac = 0.5*(ac-abmin-abmax)
      kc = 5
!
!  advance baseline spanning contour
!
 260  zc = ac-z
      if(zc.gt.0.0.and.zd.gt.0.0)goto 280
      if(zc.lt.0.0.and.zd.lt.0.0)goto 280
      if(zc.eq.0.0.and.ar.eq.z)goto 280
!  turn half right
      kl = kc
      al = ac
      k = kra(k)
      goto 290
!  turn half left
 280  kr = kc
      ar = ac
      k = kla(k)
 290  zd = ar-al
      if(k.le.4)goto 310
      kc = kca(k)
      ml = m+is(kc)
      nl = n+js(kc)
      ac = aa(ml,nl)
      ia(ml,nl)=ichk
      goto 260
!  at exit of square advance contour
 310  if(zd.ne.0.0)zd=1.0/zd
      bl = zd*(z-al)
      br = zd*(ar-z)
      x = zx(m)+(zx(m+1)-zx(m))*(br*xa(kl)+bl*xa(kr))
      y = zy(n)+(zy(n+1)-zy(n))*(br*ya(kl)+bl*ya(kr))
      icount = icount+1
!      if(icount.lt.30)print *,'drawto',m,n,ia(m,n),x,y
      call draw2(x,y)
!
!  annotation logic
!
      if(lann)then
        call pstran(x,y,xx,yy)
        xa2=xa1
        ya2=ya1
        xa1=xx
        ya1=yy
        sa2=sa1
        sa1=sa2+sqrt((xa1-xa2)**2+(ya1-ya2)**2)
        if(int(sa1/bnnd).ne.int(sa2/bnnd))then
          if(ncnn.eq.24)then
            do iaa = 1,12
              bnnx(iaa)=bnnx(iaa+iaa)
              bnny(iaa)=bnny(iaa+iaa)
              bnna(iaa)=bnna(iaa+iaa)
            enddo
            ncnn = 12
          endif
          ncnn = ncnn+1
          bnnx(ncnn) = (xa1+xa2)*0.5d0
          bnny(ncnn) = (ya1+ya2)*0.5d0
          bnna(ncnn) = degree * atan2(ya1-ya2,xa1-xa2)
        endif
      endif
!
!  and find new square
      m = m+ma(k)
      n = n+na(k)
      k = kka(k)
      kl = k
      kr = kkr(k)
!  end contour if moving out of grid
      if(m.lt.1.or.m.ge.ino)goto 400
      if(n.lt.1.or.n.ge.jno)goto 400
! or back at start of contour again
      if(it.lt.5.or.n.ne.nst.or.m.ne.mst.or.k.ne.kst)goto 210
 400  call endline
!
!  annotation logic at end of line
!
      if(lann)then
        i1t = max(2,ncnn/3)
        do i0t=1,3
          if(i1t.le.ncnn.and.nann.lt.mann)then
            nann=nann+1
            annx(nann)=bnnx(i1t)
            anny(nann)=bnny(i1t)
            anna(nann)=bnna(i1t)
            annv(nann)=z
          endif
        i1t = i1t+max(4,ncnn/3)
        enddo
      endif
      goto 100
!
!  contour ends at square containing land
!
 410  call endline
!
!  annotation logic at end of line
!
      if(lann)then
        i1t = max(2,ncnn/3)
        do i0t=1,3
          if(i1t.le.ncnn.and.nann.lt.mann)then
            nann=nann+1
            annx(nann)=bnnx(i1t)
            anny(nann)=bnny(i1t)
            anna(nann)=bnna(i1t)
            annv(nann)=z
          endif
        i1t = i1t+max(4,ncnn/3)
        enddo
      endif
!  if contour started within grid crossing line of constant j in
!  positive direction, try negative direction
      if(it.lt.5.or.i5.ne.0)goto 100
      i5 = 1
      m = mst
      n = nst-1
      k = 3
      kl = k
      kr = kkr(k)
      ar = aa(m,n+1)
      al = aa(m+1,n+1)
      goto 205
!
 100  continue
      if(it.eq.5)goto 40
 110  continue
      return
      end


      subroutine mgicon(aa,ix,jy,i1,j1,minc,b1,dv,imin,inum,
     &                  itit,tita,ierr)
      implicit none
      integer, parameter :: DP = kind(1.0d0)
!
      integer ix, jy, i1, j1, minc, imin, inum, itit, ierr
      real(DP)  aa(ix,jy),b1,dv,vmask(1)
      character*(*) tita
!
        vmask(1) = b1
        call  mgicon2(aa,ix,jy,i1,j1,minc,vmask,1,dv,imin,inum,
     &                  itit,tita,ierr)
      return
      end


      subroutine mgicon2(aa,ix,jy,i1,j1,minc,vmask,nmask,dv,imin,inum,
     &                  itit,tita,ierr)
!
!  subroutine to choose interval and range for contouring
!
!  input:
!    aa    -  array of values to contour
!    ix,iy -  dimensions of a
!    i1,j1 -  region of a to scan
!    minc  -  minimum number of contours requested
!    vmask -  array of dummy values (0d0 or 1e7,2e7,...)
!    nmask -  number of dummy values
!    itit  -  if non-zero put contour information in tita
! output:
!    dv    -  contour interval
!    imin  -  minimum contour at imin*dv
!    inum  -  number of contours
!    tita  -  title
!    ierr  -  error flag
!
      implicit none
      integer, parameter :: DP = kind(1.0d0)
!
      integer  :: ix, jy, i1,j1,minc, nmask, imin, inum, itit, ierr
      logical  :: linit
      real(DP) :: aa(ix,jy),vmask(nmask)
      integer  :: i, j, k, l, iq
      real(DP) ::  p, q, r, s, dv, aamax, aamin
      character*(*) :: tita
!
      linit = .true.
!       print *,'subroutine mgicon entered'
!       print *,ix,jy,i1,j1,minc,b1,nmask,vmask
      do 20 i=1,i1
      do 20 j=1,j1
        l=0
        do k=1,nmask
          if(aa(i,j).eq.vmask(k))l=k
        end do
        if(aa(i,j).gt.1d6.and.l.eq.0)then
          print *,'subroutine mgicon error'
          print *,ix,jy,i1,j1,vmask(1)
          print *,i,j,aa(i,j)
          print *,' programme stopping ...'
          stop
        endif
        if(l.eq.0)then
          if(linit)then
            linit = .false.
            aamin = aa(i,j)
            aamax = aa(i,j)
          else
            if(aa(i,j).lt.aamin)then
              aamin = aa(i,j)
            endif
            if(aa(i,j).gt.aamax)then
               aamax = aa(i,j)
            endif
          endif
        endif
20    continue
      print 21,' mgicon - aamin aamax: ',aamin,aamax
  21  format(a,1p2e16.6)
!
      p = abs(aamax-aamin)/minc
      q = log10(p)
      iq = int(q)
      if(q.lt.0d0)iq = iq-1
      r=10d0**iq
      s = p/r
      dv = r
      if(s.gt.2d0)dv = 2d0*r
      if(s.gt.5d0)dv = 5d0*r
!
      if(abs(aamin/dv).lt.1.0d6)then
        imin = int(aamin/dv)
        if(aamin.gt.0d0)imin=imin+1
        inum = int(aamax/dv) + 1 - imin
        if(itit.ne.0)then
          write(tita,31)dv,dv*imin,dv*(imin+inum-1)
        endif
        ierr = 0
      else
        imin = 0
        inum = 0
        write(6,*)' mgicon - error.  abs(imin) is too large.'
        ierr = 1
      endif
      return
  31  format(' contour interval = ',f10.3,
     &           ' min = ',f10.3,' max = ',f10.3)
      end



      subroutine contb(ia,xe,ye,nx,ny,nxx,nyy,linter)
!
!  Subroutine to fill area for which ia(i,j) equals 1,
!  The fill colour should be set before calling contb.
!
!  Convention: 1st iaa coordinate corresponds to row, i, x or east
!              2nd iaa coordinate corresponds to column, j, y or north
!
      use ps
      implicit none
      integer, parameter :: DP = kind(1.0d0)
!
      integer,intent(inout) :: ia(nxx,nyy)     ! data array
      integer,intent(in)    :: nxx, nyy        ! size of data array
      integer,intent(in)    :: nx, ny          ! region to process
      real(DP),intent(in)   :: xe(0:nxx),ye(0:nyy) ! edge locations
      logical, intent(in) :: linter            !
!
! local variables
!
      integer :: i, i1, j1, i2, j2, i3, j3, k2, nsides
      logical :: l1,l2,lnew, log1,log2
!
!  Outer loops:  Loop over all data cells
!
      do j1 = 1,ny
      do i1 = 1,nx
!
!  If unprocessed cell found look for a block of adjacent
!  unprocessed cells.
!
      if(ia(i1,j1).eq.1)then
        l1 = .true.          !  look further in +i direction
        l2 = .true.
        nsides = 0
        j2 = j1
        i2 = i1
        i3 = i2
        j3 = j2
        do while(l1.eqv..true..and.nsides.lt.512)
!
!  search in +i direction
!
          if(l2.eqv..true.)then
            ia(i2,j2)=2                  ! Add current cell to block
!            print *,"Adding cell ",i2,j2
!  update number of sides (nsides) bounding the current block
            if(linter)then
              log1 = i2.gt.1.and.j2.gt.1
              if(log1) log1 = ia(i2-1,j2).eq.2.and.ia(i2-1,j2-1).eq.2
     &                                        .and.ia(i2  ,j2-1).ne.2
              if(log1)nsides = nsides + 2
            else
              nsides = nsides+4
              log1 = .false.
              if(i2.gt.1) log1 = ia(i2-1,j2).eq.2
              if(log1)nsides = nsides-2
              log2 = .false.
              if(j2.gt.1)log2 = ia(i2,j2-1).eq.2
              if(log2)nsides = nsides-2
            endif
!
! If necessary process end of row
!
            if(i2.eq.nx)then                ! end of row
              if(j2.eq.ny)then
                l1 = .false.
              else
                j2 = j2+1
                l2 = .false.
              endif
!
!  Otherwise see if row can be extended
!
            else
!
!  First test to ensure it will not form a hole in the block
!
              log1 = .false.
              if(i2.lt.nx.and.j2.gt.1)log1 = ia(i2+1,j2-1).eq.2
              log2 = .false.
              if(i2.lt.nx-1.and.j2.gt.1)log2 = ia(i2+2,j2-1).ne.2
!
              if(ia(i2+1,j2).eq.1.and.
     &          (j2.eq.1.or.i2.ge.nx-1.or.log1.or.log2))then
                i2 = i2+1
!
!  If row cannot be extended try a new row.
!
               else
                 if(j2.eq.ny)then
                   l1 = .false.
                 else
                   j2 = j2+1
                  l2 = .false.
                 endif
              endif
            endif
          else
!
!  search backward from (i2,j2) where (i2,j2-1) was last point
!
            i3 = 0
!           if(linter)nsides = nsides+4
            do i=i2,1,-1
              if(j2.le.1.or.i.eq.0)then
                print *,'error ',i1,j1,i2,j2,i,nsides,l1,l2
                print *,ia(i1,j1)
                print *,ia(i2,j2)
              endif
              if(ia(i,j2).eq.1.and.ia(i,j2-1).eq.2)i3=i
            end do
            if(i3.ne.0)then
              if(linter)then
                nsides = nsides + 4
                if(j2.gt.1)then
                  log1 = .false.
                  log2 = .false.
                  if(i2.gt.1 .and.j2.gt.1)log1 = ia(i2-1,j2-1).eq.0
                  if(i2.lt.nx.and.j2.gt.1)log2 = ia(i2+1,j2-1).eq.0
                  if(i2.eq.1 .or.log1)nsides=nsides-2
                  if(i2.eq.nx.or.log2)nsides=nsides-2
                endif
              endif
              i2 = i3
              if(i2.gt.1)then
                i3=i3-1
                do while(i3.ge.1.and.ia(max(i3,1),j2).eq.1
     &                          .and.nsides.lt.512)
                  ia(i3,j2) = 2
                  if(linter)then
                    if(ia(i2,j2-1).ne.2.and.ia(i2+1,j2).eq.2)
     &                nsides = nsides+2
                  else
                    nsides = nsides+2
                  endif
                  i3 = i3-1
                end do
              endif
              l2 = .true.
            else
              l1 = .false.
            endif
          endif
        end do
!
!  Outline block of cells starting on vertex at south-west corner of
!  the first cell found
!
        i2 = i1
        j2 = j1
        call move2(xe(i2-1),ye(j2-1))
!
!  Start looking eastwards (k2=1)
!
        i2 = i2+1
        k2 = 1
        lnew = .false.
!
!  Loop around outer edge of block until first point found
!
        do while(i2.ne.i1.or.j2.ne.j1)
          if(linter)then
            if(lnew)then
              call draw2(xe(i3-1),ye(j3-1))
              lnew = .false.
            endif
            i3 = i2
            j3 = j2
          else
            call draw2(xe(i2-1),ye(j2-1))
          endif
!
! Looking east from current vertex.  There must be a
! valid block north-east or north-west of current vertex
!
          if(k2.eq.1)then
            log1 = .false.
            if(i2.le.nx.and.j2.gt.1)log1 = ia(i2,j2-1).eq.2
            log2 = .false.
            if(i2.le.nx)log2 = ia(i2,j2).eq.2
            if(log1)then        !  Valid block found to south-east
              j2=j2-1          !  move south
              k2=4
              lnew=.true.
            elseif(log2)then
              i2=i2+1          ! Valid block to north-east
            else
              j2=j2+1          ! No valid block to north-east
              k2=2             ! Move north
              lnew=.true.
            endif
!
! Looking north from current vertex.  There must be a
! valid block south-west or north-west of current vertex
!
          elseif(k2.eq.2)then
            log1 = .false.
            if(j2.le.ny.and.i2.le.nx) log1 = ia(i2,j2).eq.2
            log2 = .false.
            if(i2.gt.1.and.j2.le.ny)log2 = ia(i2-1,j2).eq.2
            if(log1)then
              i2=i2+1          ! Valid block
              k2=1
              lnew=.true.
            elseif(log2)then
              j2=j2+1
            else
              i2=i2-1
              k2=3
              lnew=.true.
            endif
          elseif(k2.eq.3)then
            log1 = .false.
            if(i2.gt.1.and.j2.le.ny)log1 = ia(i2-1,j2).eq.2
            log2 = .false.
            if(i2.gt.1.and.j2.gt.1)log2 = ia(i2-1,j2-1).eq.2
            if(log1)then
              j2=j2+1
              k2=2
              lnew=.true.
            elseif(log2)then
              i2=i2-1
            else
              k2=4
              lnew=.true.
            endif
          else
            log1 = .false.
            if(i2.gt.1.and.j2.gt.1)log1 = ia(i2-1,j2-1).eq.2
            log2 = .false.
            if(j2.gt.1)log2 = ia(i2,j2-1).eq.2
            if(log1)then
              i2=i2-1
              k2=3
              lnew=.true.
            elseif(log2)then
              j2=j2-1
            else
              k2=1
              lnew=.true.
            endif
          endif
        end do
        if(linter)call draw2(xe(i3-1),ye(j3-1))
        call closepath
        call fill
!
!  set region to 3
!
        do j2 = 1,ny
        do i2 = 1,nx
          if(ia(i2,j2).eq.2)ia(i2,j2)=3
        end do
        end do
      endif
      end do
      end do
!
      return
      end

