MODULE comm

  USE petsc
#include "petsc/finclude/petsc.h"
  use parameters
  IMPLICIT NONE
  
contains

!============================================================================================

  subroutine force_proc(nf, nf_sd, xfm, xsdm)
  !Identify how the force vector containing nf unknowns is distributed across processors
  
    PetscInt :: nf, proc_lim, xfm, xf_diff, nf_sd, xsdm, xsd_diff
  
    proc_lim=300  !need to determine this value by trial and error
    
    if (num_stat) nproc_force = nproc  !if numerically stationary, then divide among all processors so that Elemental library could be used   
    if (sub_domain) nproc_force = nproc  
    if (.not. num_stat .and. .not. sub_domain) nproc_force = int(nf/proc_lim+1)
    
    if (nproc_force .gt. nproc) nproc_force = nproc  !this can happen when less procs are used
  
    if (rank< nproc_force) then
      xfm = nf/nproc_force
      xf_diff = nf - xfm*nproc_force
      if (rank<xf_diff) xfm = xfm+1
      
      xsdm = nf_sd/nproc_force
      xsd_diff = nf_sd - xsdm*nproc_force
      if (rank<xsd_diff) xsdm = xsdm+1
    else
      xfm = 0
      xsdm = 0
    end if
  
  end subroutine force_proc
  
!============================================================================================
 
END MODULE comm
