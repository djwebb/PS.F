!
!=======================================================================
!
!    File: transeu.F
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
      subroutine transeu(rlon1, rlat1, rlon2, rlat2, deg_in, deg_out,
     &                          alpha, beta,  gamma, deg_eu)
      implicit none
      real*8  rlon2, rlat2, rlon1, rlat1, alpha, beta, gamma
      integer deg_in, deg_out, deg_eu
c
c-----------------------------------------------------------------------
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c PURPOSE To convert latitude and longitude in one spherical coordinate
c         system to latitude and longitude in a second system.  The
c         angles alpha, beta and gamma are the three Euler angles
c         defining the transform from the first system to the second.
c         The rotation is defined here is the one commonly used with
c         spherical harmonics:
c           1.  Rotate meridian by angle alpha around the pole.
c           2.  Rote pole down angle beta along the meridian.
c           3.  Rotate new meridian by angle gamma around the new pole.
c         Another commonly used definition, used in dynamics, is:
c           1.  Rotate meridian by angle alpha' around pole.
c           2.  Rote pole down angle beta', keeping the point where
c               the meridian crosses the equator stationary.  (This
c               direction defines the astronomical node).
c           3.  Rotate new meridian by angle gamma' around the new pole.
c          The angles are related by the equations:
c               alpha' = alpha + 90 degrees
c               beta'  = beta
c               gamma' = gamma - 90 degrees
c
c           The back transformation corresponding to (alpha,beta,gamma)
c           is (-gamma,-beta,-alpha).
c
c INPUT     rlat1  - model 1 latitude
c           rlon1  - model 1 longitude
c           deg_in - 1 if input values in degrees
c                    0 if input values in radians
c           deg_out- 1 if output values in degrees
c                    0 if output values in radians
c
c           alpha  - first Euler angle
c           beta   - second Euler angle (latitude)
c           gamma  - third Euler angle
c           deg_eu - 1 if Euler angles in degrees
c                    0 if Euler angles in radians
c
c OUTPUT    rlat2  - model 2 latitude
c           rlon2  - model 2 longitude
c
c  OCCAM    For the forward transform from the model 1 grid
c           to the model 2 grid, OCCAM uses:
c           alpha = 232 degrees
c           beta  =  90 degrees
c           gamma =   0 degrees
c
c  Last modified by djw, 15 June 2010
c-----------------------------------------------------------------------
c
      integer in
      real*8  a, b, c, p1, q1, p2, q2, s, t
      real*8  pi, twopi, radian, degree
      save in, pi,twopi,degree,radian
      data in/0/
c
c-----------------------------------------------------------------------
c  constants
c-----------------------------------------------------------------------
c
      if(in.eq.0)then
c        pi = 3.1415926535897932d0
        pi = 4d0*atan(1d0)
        twopi = 2d0*pi
        radian= pi/180d0
        degree= 180d0/pi
        in = 1
      endif
c
c-----------------------------------------------------------------------
c  copy across arguments and scale if necessary
c-----------------------------------------------------------------------
c
      if(deg_eu.ne.0)then
        a = alpha*radian
        b = beta*radian
        c = gamma*radian
      else
        a = alpha
        b = beta
        c = gamma
      endif
c
      if(deg_in.ne.0)then
        p1 = radian*(90d0-rlat1)
        q1 = radian*(rlon1-alpha)
      else
        p1 = 0.5d0*pi-rlat1
        q1 = rlon1-alpha
      endif
c
c  transform q1 to range -pi to +pi
c
      q1 = mod(q1,twopi)
      if(q1.gt.pi) q1=q1-twopi
      if(q1.le.-pi)q1=q1+twopi
c      print *
c      print 11,a,b,c
c      print 11,p1,q1
c  11  format(5x,5f12.6)
c
      p2 = acos(cos(b)*cos(p1)+sin(b)*sin(p1)*cos(q1))
c
      if(b.eq.0)then
        q2 = q1
      elseif(p2.eq.0d0)then
        q2 = 0d0
      else
        s = (p1+p2+b)*0.5d0
        t = sin(s-b)*sin(s-p2)/(sin(b)*sin(p2))
c  correct for small rounding errors near the limits
        if(t.lt.0d0)t = 0d0
        if(t.gt.1d0)t = 1d0
        q2 = 2d0*acos(sqrt(t))
        if(q1.lt.0d0)q2 = -q2
c
c        if(isnan2(q2))then
c          print *,"File transeu.F.  Error in routine transeu:"
c          print *," variable a = ",a
c          print *," variable b = ",b
c          print *," variable c = ",c
c          print *," variable p1 = ",q2
c          print *," variable q1 = ",q1
c          print *," variable s = ",s
c          print *," variable t = ",t
c          print *," variable p2 = ",p2
c          print *," variable q2 = ",q2
c        endif
c
c        print 11,s,t,q2
      endif
c
c-----------------------------------------------------------------------
c  rescale the output value then and return to a float
c  from the double values used in the conversion
c  program
c-----------------------------------------------------------------------
c
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
      end

