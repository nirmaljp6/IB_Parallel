PROGRAM main

  USE petsc
#include "petsc/finclude/petsc.h"
  USE parameters
  USE grid
  USE variables
  USE timestepper
  USE operators_pressure
  
  implicit none
  PetscInt :: it
  integer :: nproc_mpi, rank_mpi

  !Initialize Petsc: nproc is number of processors, rank is the current processor
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc_mpi, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank_mpi, ierr)
  nproc = nproc_mpi
  rank = rank_mpi

  call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INDEX, ierr)
  call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_INDEX, ierr)
  
  call input
  call setup
  CALL setup_pressure
  call setup_variables
  
  it=istart
  call write_proc
  
  do while (it<istop)
    call advance(it)
 
    IF (it.ge.irestart) THEN 
       IF ((MOD(it-irestart,isave).eq.0).or.(it==istop)) THEN
         CALL write_variables(it)    
         
         IF (compute_pressure) THEN
             CALL calculate_pressure( it )
             call write_pressure(it)
         END IF
       END IF
    END IF
    
    call write_force(it,fvec)
    call write_force_rdst(it,f_rdstvec)
    call write_theta(it)
    if (rank==0 .and. .not. num_stat) call write_iterations(it)
  end do  

  call destroy_variables
  CALL destroy_pressure
  call destroy_grid
  
  call PetscFinalize(ierr)

END PROGRAM main
