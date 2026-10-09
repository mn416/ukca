! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!     Calculates the full Jacobian matrix. Needed in diagnostic output and
!     in dense algebra. Note that O(1D) and O(3P) are the only SS species
!     treated correctly here.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Purpose
!     -------
!     Calculates the full Jacobian matrix. Needed in diagnostic output and
!     in dense algebra. Note that O(1D) and O(3P) are the only SS species
!     treated correctly here.
!
!     Interface
!     ---------
!     Called from within the Newton-Raphson solver (spimpmjp).
!
!     Arguments:
!        n_points    - No. of points calculations be done.
!
!     Method
!     ------
!
!     From the chemical reaction data, calculates the derivative w.r.t.
!     all other predicted chemicals.
!
!     Externals
!     ---------
!     (no externals)
!
!
! ######################################################################
!
MODULE asad_fuljac_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_FULJAC_MOD'

CONTAINS

SUBROUTINE asad_fuljac(s, n_points)
!$OMP DECLARE TARGET
USE asad_mod,        ONLY: ctype, frpx, jpfm, jpif, jpmsp,                     &
                           madvtr, moffam, ndepd, ndepw, nfrpx,                &
                           njcoth, nltrf, nmzjac, nmsjac, nodd,                &
                           npdfr, nsjac1, nstst, ntabpd, ntrf,                 &
                           ntrho2, ntro3, ntroh, ntrno, nzjac1,                &
                           ztabpd, jpspec, jpcspf, asad_state_type
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
TYPE(asad_state_type), INTENT(INOUT) :: s
INTEGER, INTENT(IN) :: n_points

! Local variables
REAL :: deltt
REAL :: fr
INTEGER :: irj
INTEGER :: ifamd
INTEGER :: itrd
INTEGER :: i
INTEGER :: j
INTEGER :: jc
INTEGER :: itrcr
INTEGER :: j3
INTEGER :: jn
INTEGER :: js
INTEGER :: jl
INTEGER :: i1
INTEGER :: i2
INTEGER :: iss
CHARACTER :: ityped*2
INTEGER :: ij(jpmsp)
LOGICAL :: lj(jpmsp)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_FULJAC'

!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
deltt=1.0/s%cdt
!
s%fj = 0.0
DO i=1,jpcspf
  s%fj(:,i,i)=-deltt*s%f(:,i)
END DO
!
!
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
DO jc = 1, ntrf
  itrcr = nltrf(jc)
  DO j3 = 1,nmzjac(itrcr)
    irj = nzjac1(j3,itrcr)
    !
    DO jn=1,jpmsp
      ij(jn)=njcoth(irj,jn)
      lj(jn)=ij(jn) /= 0
    END DO
    DO jn=1,2
      IF (lj(jn)) THEN
        DO jl=1,n_points
          s%fj(jl,ij(jn),itrcr) =                                              &
                        s%fj(jl,ij(jn),itrcr) - s%prk(jl,irj)
        END DO
      END IF
    END DO
    IF (nfrpx(irj) == 0) THEN
      DO jn=3,jpmsp
        IF (lj(jn)) THEN
          DO jl=1,n_points
            s%fj(jl,ij(jn),itrcr) =                                            &
                        s%fj(jl,ij(jn),itrcr) + s%prk(jl,irj)
          END DO
        END IF
      END DO
    ELSE
      DO jn=3,jpmsp
        IF (lj(jn)) THEN
          fr=frpx(nfrpx(irj)+jn-3)
          DO jl=1,n_points
            s%fj(jl,ij(jn),itrcr) =                                            &
                        s%fj(jl,ij(jn),itrcr) + fr*s%prk(jl,irj)
          END DO
        END IF
      END DO
    END IF
    !
    IF (npdfr(irj,1) /= 0) THEN
      i1 = npdfr(irj,1)
      i2 = npdfr(irj,2)
      DO jn = i1, i2
        iss = ntabpd(jn,1)
        fr = ztabpd(jn,1)
        DO jl=1,n_points
          s%fj(jl,iss,itrcr) = s%fj(jl,iss,itrcr) + fr*s%prk(jl,irj)
        END DO
      END DO
    END IF
    !
  END DO
END DO

! Go through the steady state additions to the Jacobian

DO jc = 1, nstst
  DO j3 = 1,nmsjac(jc)
    irj = nsjac1(j3,jc)
    DO jn=1,jpmsp
      ij(jn)=njcoth(irj,jn)
      lj(jn)=ij(jn) /= 0
    END DO
    DO jn=1,2
      IF (lj(jn)) THEN
        DO jl=1,n_points
          s%fj(jl,ij(jn),ntro3 ) =                                             &
             s%fj(jl,ij(jn),ntro3 ) - s%prk(jl,irj)*s%deriv(jl,jc,1)
          s%fj(jl,ij(jn),ntroh ) =                                             &
             s%fj(jl,ij(jn),ntroh ) - s%prk(jl,irj)*s%deriv(jl,jc,2)
          s%fj(jl,ij(jn),ntrho2) =                                             &
             s%fj(jl,ij(jn),ntrho2) - s%prk(jl,irj)*s%deriv(jl,jc,3)
          s%fj(jl,ij(jn),ntrno ) =                                             &
             s%fj(jl,ij(jn),ntrno ) - s%prk(jl,irj)*s%deriv(jl,jc,4)
        END DO
      END IF
    END DO
    DO jn=3,jpmsp
      IF (lj(jn)) THEN
        s%fj(1:n_points,ij(jn),ntro3) =                                        &
          s%fj(1:n_points,ij(jn),ntro3)  +                                     &
          s%prk(1:n_points,irj)*s%deriv(1:n_points,jc,1)
        s%fj(1:n_points,ij(jn),ntroh) =                                        &
          s%fj(1:n_points,ij(jn),ntroh)  +                                     &
          s%prk(1:n_points,irj)*s%deriv(1:n_points,jc,2)
        s%fj(1:n_points,ij(jn),ntrho2) =                                       &
          s%fj(1:n_points,ij(jn),ntrho2) +                                     &
          s%prk(1:n_points,irj)*s%deriv(1:n_points,jc,3)
        s%fj(1:n_points,ij(jn),ntrno ) =                                       &
          s%fj(1:n_points,ij(jn),ntrno ) +                                     &
          s%prk(1:n_points,irj)*s%deriv(1:n_points,jc,4)
      END IF
    END DO
  END DO
END DO
!
!     -----------------------------------------------------------------
!          4.  Add deposition terms to Jacobian diagonal.
!              --- ---------- ----- -- -------- ---------
!
IF ( ndepw  /=  0 .OR. ndepd  /=  0 ) THEN
  DO js = 1, jpspec
    ifamd = moffam(js)
    itrd = madvtr(js)
    ityped = ctype(js)
    !
    IF ( ifamd /= 0 ) THEN
      DO jl=1,n_points
        IF ((ityped == jpfm) .OR. (ityped == jpif .AND.                        &
        s%linfam(jl,itrd)))                                                    &
           s%fj(jl,ifamd,ifamd)=s%fj(jl,ifamd,ifamd)                           &
         - nodd(js)*(s%dpd(jl,js)+s%dpw(jl,js))                                &
          *s%y(jl,js)
      END DO
    END IF
    IF ( itrd /= 0 ) THEN
      s%fj(1:n_points,itrd,itrd) = s%fj(1:n_points,itrd,itrd)                  &
       - (s%dpd(1:n_points,js)+s%dpw(1:n_points,js))*s%y(1:n_points,js)
    END IF
  END DO
END IF
!
!     -------------------------------------------------------------
!          5.  Jacobian elements in final form
!              -------- -------- -- ----- ----
!
DO j=1,jpcspf
  s%fj(1:n_points,j,:)=s%fj(1:n_points,j,:)/s%f(1:n_points,:)
  ! filter s%f earlier!
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_fuljac
END MODULE asad_fuljac_mod
