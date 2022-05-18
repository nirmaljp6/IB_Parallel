module user

  USE petsc
#include "petsc/finclude/petsc.h"
  use parameters
  use grid
  IMPLICIT NONE


contains

!================================================================================
  FUNCTION motion_grid( it )
    
    PetscInt :: it
    PetscScalar, dimension(5) :: motion_grid
    PetscScalar ::  rx, ry
    
    ! motion_grid specifies the (possibly time-varying) u component, v component of the velocity, angular velocity, and the coordinates of the center of rotation when the body is moving with the grid.
    
    rx=0.d0
    ry=0.d0
    
    ! default is no flow
    motion_grid(1) = -1.D0   ! x-component of velocity in the body-fixed frame
    motion_grid(2) = 0.d0 ! y-component of velocity in the body-fixed frame
    motion_grid(3) = 0.d0   ! angular velocity: I think it is assumed negative in clockwise direction
    motion_grid(4) = rx    ! x-coordinate of the center of rotation
    motion_grid(5) = ry    ! y-coordinate of the center of rotation

  END FUNCTION motion_grid
  
!================================================================================

  FUNCTION bodyforcex(it, x, y, u, v ) 

    REAL(KIND(0.D0)) :: x,y,u,v, bodyforcex
    PetscInt  :: it
  
    ! specify a body force (per unit mass) on RHS of u-momentum equation,
    ! as a function of location (x,y) and possibly velocity (u,v) and possibly time

    ! the specified function should decay to zero near the edges of the (mgrid=1)
    ! inner domain
    !
    ! units are in terms of velocity (result should be consistent
    ! with U^2/L units

    bodyforcex = 0.d0 ! do not remove; overwrite values below inside if block

  END FUNCTION bodyforcex

!================================================================================

  FUNCTION bodyforcey(it, x, y, u, v ) 

    REAL(KIND(0.D0)) :: x,y,u,v, bodyforcey
    PetscInt  :: it
  
    ! specify a body force (per unit mass) on RHS of v-momentum equation,
    ! as a function of location (x,y) and possibly velocity (u,v) and possibly time

    ! the specified function should decay to zero near the edges of the (mgrid=1)
    ! inner domain
    !
    ! units are in terms of velocity (result should be consistent
    ! with U^2/L units

    bodyforcey = 0.d0 ! do not remove; overwrite values below inside if block


    !Add small body force for a small unit of time
    !xc = 0.5D0 ! center of the bodyforce
    !yc = -0.2d0

    !IF ((real(it)*dt) .ge. 115.0d0 .and. real(it)*dt .le. 116.0d0) then
    !
    !    IF ( ( (x-xc)**2.d0+(y-yc)**2.d0 ) .le. (10.D0*delta)**2.d0 ) THEN
    !        bodyforcey = 5.d-4
    !    end if
    !end IF

  END FUNCTION bodyforcey

!================================================================================

  subroutine advance_rigidbody(it, xbval, ybval)
  !This updates the position of the moving rigid body
     PetscInt :: it, i, i_bdy
     PetscScalar :: xbval(xfi:xfe), ybval(xfi:xfe)
     PetscScalar :: Ax, freq, xc, xct, alpha, rr, alpha0, beta
     
     Ax = 2.8D0
     freq = 0.25D0/2.2D0
     alpha0 = pi/2.D0
     beta = pi/4.D0
     
     xc = Ax/2.D0*cos(2.D0*pi*freq*dt*it)
     xct = Ax/2.D0*cos(2.D0*pi*freq*dt*(it-1))
     alpha = alpha0 + beta*sin(2*pi*freq*dt*it)
           
     do i=xfi,xfe
        if (codeb(i,1)==1  ) then
           rr = sqrt( (xbval(i)-xct)**2.D0 + (ybval(i))**2.D0 )
           if ( (i .ge. 1 .and. i .le. (nb-1)/2+1) .or. (i .ge. nb+1 .and. i .le. (nb-1)/2+nb+1) ) then
              !xbval(i) = xc + rr*cos(alpha)
              !ybval(i) = rr*sin(alpha)
           else
              !xbval(i) = xc - rr*cos(alpha)
              !ybval(i) = -rr*sin(alpha)
           end if
        end if
     end do
  
  end subroutine advance_rigidbody
  
!================================================================================

  subroutine body_user_prescribed(it, i_bdy)
  !This keeps track where the hinge of the torsional body is and what angle the new equilibrium position of the torsional body is inclined to w.r.t to the previous body.
  !this subroutine is local to each proc
  
     PetscInt :: i_bdy, it
     PetscScalar :: theta_prescribed, Ax, freq
     
     Ax = 1.D0
     freq = 1/(1.D0*pi)
     
     if (i_bdy==10) then 
        
     end if
     
  end subroutine body_user_prescribed

!================================================================================

  subroutine compute_rhs_prescribed(it, rhsf, xbval, ybval)
     !add -u_B*delta into the RHS. 
     
     PetscInt :: it, i
     PetscScalar :: rhsf(xfi:xfe), xbval(xfi:xfe), ybval(xfi:xfe)
     PetscScalar :: Ax, freq, alpha0, beta, xc, alpha, xcd, alphad, rr
     
     Ax = 2.8D0
     freq = 0.25D0/2.2D0
     alpha0 = pi/2.D0
     beta = pi/4.D0
     
     xc = Ax/2.D0*cos(2.D0*pi*freq*dt*it)
     alpha = alpha0 + beta*sin(2*pi*freq*dt*it)
     xcd = -1.D0*pi*Ax*freq*sin(2.D0*pi*freq*dt*it)
     alphad = 2.D0*pi*beta*freq*cos(2.D0*pi*freq*dt*it)
     
     do i=xfi,xfe
        rr = sqrt( (xbval(i)-xc)**2.D0 + (ybval(i))**2.D0 )
        
        if ( (i .ge. 1 .and. i .le. (nb-1)/2+1) .or. (i .ge. nb+1 .and. i .le. (nb-1)/2+nb+1) ) then
           if (codeb(i,3)==1  ) then
              !rhsf(i) = rhsf(i) - delta*(xcd - rr*sin(alpha)*alphad)
           else if (codeb(i,3)==2  ) then
              !rhsf(i) = rhsf(i) - delta*(rr*cos(alpha)*alphad)
           end if
        else
           if (codeb(i,3)==1  ) then
              !rhsf(i) = rhsf(i) - delta*(xcd + rr*sin(alpha)*alphad)
           else if (codeb(i,3)==2  ) then
              !rhsf(i) = rhsf(i) - delta*(-1.D0*rr*cos(alpha)*alphad)
           end if
        end if

     end do
  
  
  end subroutine compute_rhs_prescribed  

!================================================================================  

  subroutine compute_gravity(i_bdy, thet, gravity)
  
    PetscInt :: i_bdy
    PetscScalar :: gravity, rho_ratio, aspect_ratio, g_param, rr, thet
    
    rho_ratio = 20.D0
    aspect_ratio = 0.01D0
    g_param = 0.D0
    
    rr = sqrt((bdy(i_bdy)%x(bdy(i_bdy)%npts) - bdy(i_bdy)%x(1))**2 + (bdy(i_bdy)%y(bdy(i_bdy)%npts) - bdy(i_bdy)%y(1))**2)
    gravity = (rho_ratio-1)*aspect_ratio*g_param*rr/2.D0*cos(bdy(i_bdy)%gamma0 - thet)
    
  end subroutine compute_gravity  
  
!================================================================================  

  subroutine compute_dynamics(it, i_bdy, thet, dynamics)
  
    PetscInt :: i_bdy, it
    PetscScalar :: dynamics, mass_ratio, thet, xdd, Ax, freq
    
    mass_ratio = 0.2D0
    Ax = 1.D0
    freq = 1/(1.D0*pi)

    xdd = -1.D0*(2.D0*pi*freq)**2.D0*Ax/2.D0*cos(2.D0*pi*freq*dt*it)
    dynamics = -0.D0*mass_ratio/2.D0*sin(bdy(i_bdy)%gamma0 - thet)*xdd
    
  end subroutine compute_dynamics  
  
!================================================================================

end module user
