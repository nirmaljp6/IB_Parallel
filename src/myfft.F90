module myfft

  USE petsc
#include "petsc/finclude/petsc.h"

  implicit none
  include 'fftw3-mpi.f03'  
  integer(C_INTPTR_T) :: mm, nn, alloc_local, local_n, local_n_p
  real(C_DOUBLE), pointer :: in(:,:), in_ddti(:,:)
  PetscScalar :: normalize
  type(C_PTR) :: forward, forward_ddti
  
  PetscScalar, DIMENSION(:),     ALLOCATABLE :: viscfac,vfac,con1,con2
  PetscScalar, DIMENSION(:,:),   ALLOCATABLE :: laminv, laminv_ddti
  PetscScalar, DIMENSION(:,:,:), ALLOCATABLE :: lam1, lam1i, lam1inv
  
contains

!=========================================================================

  subroutine fft_grid(ly)
  !Determine the number of rows in the flow domain to be assigned to each processor
  
    use parameters
    
    integer(C_INTPTR_T) :: local_j_offset  
    PetscScalar :: ln
    PetscInt :: ly(nproc)
    Vec :: proc_dist, proc_disty
    VecScatter :: proc_ctx
    PetscScalar, pointer :: ly_pt(:)
    
    mm=m
    nn=n
    
    call fftw_mpi_init()
    if (dimen==2) then
    alloc_local = fftw_mpi_local_size_2d(nn-1, mm-1, MPI_COMM_WORLD, local_n, local_j_offset)
    end if
    ln=local_n
    
    !Scattering local_n from all processes to all processes for creation of dmda
    call VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nproc, proc_dist, ierr)
    call VecSetValue(proc_dist, rank, ln, INSERT_VALUES, ierr)
    call VecScatterCreateToAll(proc_dist, proc_ctx, proc_disty, ierr)
    call VecScatterBegin(proc_ctx, proc_dist, proc_disty, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(proc_ctx, proc_dist, proc_disty, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayF90(proc_disty, ly_pt, ierr)
    ly=int(ly_pt)
    call VecRestoreArrayF90(proc_disty, ly_pt, ierr)
    call VecScatterDestroy(proc_ctx, ierr)
    call VecDestroy(proc_disty, ierr)
    call VecDestroy(proc_dist, ierr)    
            
  end subroutine fft_grid
  
!=========================================================================

  subroutine fft_grid_pressure(ly)
  
    use parameters
    
    integer(C_INTPTR_T) :: local_j_offset  
    PetscScalar :: ln
    PetscInt :: ly(nproc)
    Vec :: proc_dist, proc_disty
    VecScatter :: proc_ctx
    PetscScalar, pointer :: ly_pt(:)
    
    mm=m
    nn=n
    
    call fftw_mpi_init()
    if (dimen==2) then
    alloc_local = fftw_mpi_local_size_2d(nn, mm, MPI_COMM_WORLD, local_n_p, local_j_offset)
    end if
    ln=local_n_p
    
    !Scattering local_n from all processes to all processes for creation of dmda
    call VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nproc, proc_dist, ierr)
    call VecSetValue(proc_dist, rank, ln, INSERT_VALUES, ierr)
    call VecScatterCreateToAll(proc_dist, proc_ctx, proc_disty, ierr)
    call VecScatterBegin(proc_ctx, proc_dist, proc_disty, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(proc_ctx, proc_dist, proc_disty, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayF90(proc_disty, ly_pt, ierr)
    ly=int(ly_pt)
    call VecRestoreArrayF90(proc_disty, ly_pt, ierr)
    call VecScatterDestroy(proc_ctx, ierr)
    call VecDestroy(proc_disty, ierr)
    call VecDestroy(proc_dist, ierr)    
            
  end subroutine fft_grid_pressure

!=========================================================================
 
 subroutine setup_fft(xsi, xse, ysi, yse, xpi, xpe, ypi, ype, delta)
 
   use parameters
 
   type(C_PTR) :: p
   PetscScalar :: del2, del22, delta
   PetscInt :: k, i, j
   PetscInt :: xsi, ysi, xse, yse, xpi, xpe, ypi, ype
   PetscScalar, dimension(xsi:xse,ysi:yse) :: lam
   PetscScalar, dimension(xpi:xpe,ypi:ype) :: lam_ddti
   PetscScalar :: normalize_ddti
   
   p=fftw_alloc_real(alloc_local)
   call c_f_pointer(p, in, [mm-1,local_n])
   call c_f_pointer(p, in_ddti, [mm,local_n_p])
   
   ALLOCATE( viscfac(mgridlev), vfac(mgridlev), con1(mgridlev), con2(mgridlev) )
   ALLOCATE( laminv(xsi:xse,ysi:yse), &
             lam1(xsi:xse,ysi:yse,mgridlev), &
             lam1i(xsi:xse,ysi:yse,mgridlev), &
             lam1inv(xsi:xse,ysi:yse,mgridlev), &
             laminv_ddti(xpi:xpe,ypi:ype) )
             
   forward = fftw_mpi_plan_r2r_2d(nn-1,mm-1,in,in, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)
   forward_ddti = fftw_mpi_plan_r2r_2d(nn,mm,in_ddti,in_ddti, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)
              
   normalize      = 4.d0*REAL(mm*nn)
   normalize_ddti = 4.d0*REAL( (mm+1)*(nn+1) )                         
   
   ! eigenvalues for inverse of C^T C
   del2 = delta*delta
   DO k=1,mgridlev
     del22      =  del2*4.d0**(k-1)
     viscfac(k) =  dt/Re/del22
     vfac(k)    =  0.5d0*dt/Re/del22/normalize
     con1(k)    =  1.5d0*dt/del22/normalize
     con2(k)    = -0.5d0*dt/del22/normalize
   ENDDO
   
   DO j=ysi,yse
     DO i=xsi,xse
       ! turn lam into local variable, because it is only used here
       lam(i,j) = -2.d0*( COS( pi*REAL(i-1)/REAL(mm) ) + COS( pi*REAL(j-1)/REAL(nn) ) - 2.d0 )
       laminv(i,j) = 1.d0/lam(i,j)/normalize
       
       ! integrating factors for viscous terms
       DO k=1,mgridlev
         del22 = del2* 4.d0**(k-1)
         lam1(i,j,k)    =      (1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)/normalize ! original
         lam1i(i,j,k)   = 1.d0/(1.d0 + 0.5d0*dt*lam(i,j)/del22/Re) ! original
       END DO
     END DO
   END DO
   
   DO j=ypi,ype
     DO i=xpi,xpe
         lam_ddti(i,j) = 2.d0*( COS( pi*REAL(i)/REAL(mm+1) ) + COS( pi*REAL(j)/REAL(nn+1) ) - 2.d0 )
         laminv_ddti(i,j) = 1.d0/lam_ddti(i,j)/normalize_ddti        
     end do
   end do

 end subroutine setup_fft

!=========================================================================

  subroutine ainv(xsi, xse, ysi, yse, vort)
  
  PetscInt :: xsi, xse, ysi, yse
  PetscScalar :: vort(xsi:xse,ysi:yse)

  
  in = vort
  CALL fftw_mpi_execute_r2r(forward,in,in)
  in = lam1i(:,:,1) * in / normalize
  CALL fftw_mpi_execute_r2r(forward,in,in)
  vort = in

  end subroutine ainv
  
!=========================================================================

  subroutine ctci(xsi, xse, ysi, yse, vort, s)
  
    PetscInt :: xsi, xse, ysi, yse
    PetscScalar :: vort(xsi:xse,ysi:yse), s(xsi:xse,ysi:yse)
    
    in = vort
    call fftw_mpi_execute_r2r(forward,in,in)
    in = laminv * in
    call fftw_mpi_execute_r2r(forward,in,in)
    s = in
  
  end subroutine ctci
  
!=========================================================================

  subroutine ddti(xpi, xpe, ypi, ype, pressure, p)
  
    PetscInt :: xpi, xpe, ypi, ype
    PetscScalar :: pressure(xpi:xpe,ypi:ype), p(xpi:xpe,ypi:ype)
    
    in_ddti = pressure
    call fftw_mpi_execute_r2r(forward_ddti,in_ddti,in_ddti)
    in_ddti = laminv_ddti * in_ddti
    call fftw_mpi_execute_r2r(forward_ddti,in_ddti,in_ddti)
    p = in_ddti
  
  end subroutine ddti  

!=========================================================================

  function dst(xsi, xse, ysi, yse, psi)
  
    PetscInt :: xsi, xse, ysi, yse
    PetscScalar :: psi(xsi:xse,ysi:yse), dst(xsi:xse,ysi:yse)
    
    in  = psi
    CALL fftw_mpi_execute_r2r(forward,in,in)
    dst = in
  
  end function dst
!=========================================================================

  subroutine destroy_fft
  
     call fftw_destroy_plan(forward)
     call fftw_destroy_plan(forward_ddti)
     deallocate(laminv, lam1, lam1i, lam1inv)
  
  end subroutine destroy_fft

!=========================================================================

end module myfft
