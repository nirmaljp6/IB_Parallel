MODULE parameters

  USE petsc
#include "petsc/finclude/petsc.h"

  IMPLICIT NONE
  
  PetscErrorCode :: ierr
    
  PetscInt :: nproc, rank, nproc_force
  
  ! parameters
  PetscInt :: istart                ! initial time index
  PetscInt :: istop                ! last time index
  PetscInt :: isave                 ! save data every isave steps
  PetscInt :: irestart              !start saving variables from irestart time index
  PetscInt :: m                    ! cells in x
  PetscInt :: n                   ! cells in y
  REAL(KIND(0.0D0)) :: dt     ! time step
  REAL(KIND(0.0D0)) :: Re      ! Reynolds number
  
  REAL(KIND(0.0D0)) :: fsi_tol   ! tol. for convergence of fsi algorithm
  REAL(KIND(0.0D0)) :: atol   ! absolute tol. for GMRES convergence
  REAL(KIND(0.0D0)) :: rtol   ! relative tol. for GMRES convergence
  
  
  PetscInt :: mr, mt, md, m_body          ! number of bodies (will be counted in grid.f90) 
  REAL(KIND(0.D0)) :: len   ! length scale for grid in x direction
  REAL(KIND(0.D0)) :: offsetx   ! offset for grid in x
  REAL(KIND(0.D0)) :: offsety    ! offset for grid in y
  PetscInt :: mgridlev
  REAL(KIND(0.D0)) :: pi

  REAL(KIND(0.D0)) :: rot_angle  ! rotating angle of the grid
  REAL(KIND(0.D0)) :: rox  ! x-coord. of center of rotation rotation
  REAL(KIND(0.D0)) :: roy  ! y-coord. of center of rotation rotation

  LOGICAL :: compute_pressure  ! whether to output pressure
  
  PetscInt :: dimen  !spatial dimensions of the problem
  logical :: num_stat  !whether the code is numerically stationary or not. It is true even in the case of non-inertial frame of reference where the body is stationary
  logical :: motion_prescribed  !whether motion to any rigid body is prescribed or not
  logical :: sub_domain  !whether the efficient subdomain approach is used or not
  logical :: sub_domain_full_rectangular  !whether the subdomain is fully rectangular or not. Set it to F if there is *only* one rigid body that is stationary. 
  logical :: sub_domain_precomputed  !whether the matrix has been precomputed before or not
  REAL(KIND(0.D0)) :: sdxi, sdxe, sdyi, sdye   ! 4 coordinates of the rectangular subdomain
  REAL(KIND(0.D0)) :: delta_sd    ! spacing of grid points in the subdomain. Set it to the grid spacing of the flow domain

CONTAINS
  
  SUBROUTINE input

    LOGICAL :: readinput

    NAMELIST /read_parameters/ istart,istop,isave,irestart,m,n,dt,Re, &
                                fsi_tol,atol,rtol,len,offsetx,offsety, &
                                mgridlev, compute_pressure, &
                                dimen, num_stat, motion_prescribed, sub_domain, sub_domain_full_rectangular, sub_domain_precomputed, &
                                sdxi, sdxe, sdyi, sdye, delta_sd

    pi  = 4.0d0*atan(1.0d0)
                                 
    ! read input
    INQUIRE(file='input/ib.inp',exist=readinput)

    IF (readinput) THEN
       OPEN(unit=3,file='input/ib.inp',form='formatted',status='old')
       READ(unit=3,nml=read_parameters)
       CLOSE(3)
      ELSE
       STOP 'cannot find input file'
    END IF

    if (rank==0) then
      write(*,*) 'read input file'
    end if

  END SUBROUTINE input

END MODULE parameters


