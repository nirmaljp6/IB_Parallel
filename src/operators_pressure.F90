MODULE operators_pressure
 
  USE petsc
#include "petsc/finclude/petsc.h"
  USE parameters
  USE grid
  implicit none
  
  DM :: dap, dap_coarsify
  PetscInt :: xpcgi, xpcge, ypcgi, ypcge, xpcgm, ypcgm
  Vec :: pressurevec, pressurevec_local, pvec_local
  Vec, pointer :: pvec(:)
  PetscInt :: xci_p, xce_p, yci_p, yce_p, tcoarse_p
  PetscInt :: top_phi,bottom_phi,left_phi,right_phi, nbc_p
  PetscInt, allocatable :: cbcix_p1(:), cbcix_p2(:), bcix_p2(:)
  PetscInt :: tbc_p1, tbc_p2
  IS :: bcis_p2
  Vec :: pressurebcvec_local, pressurebcvec_global
  VecScatter :: ctx_bc_p2
    
contains    
 
!=======================================================================

subroutine setup_pressure()

  use myfft
  PetscInt :: ly(nproc), lpy(nproc), lpx(1)
  
  !Determine how many rows should 1 proc take
  call fft_grid_pressure(ly)
  
  if (dimen==2) then
      !Domain decomposition for pressure
      lpx(1)=m
      lpy=ly
      call DMDAcreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, m, n, ione, nproc, ione, ione, lpx, lpy, dap, ierr)
      
      call DMDAcreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, m, n, ione, nproc, ione, ithree, lpx, lpy, dap_coarsify, ierr)
  
  end if
  
  call DMSetFromOptions(dap, ierr)
  call DMSetUp(dap, ierr)
  call DMSetFromOptions(dap_coarsify, ierr)
  call DMSetUp(dap_coarsify, ierr)
  
  call DMDAGetCorners(dap, xpi, ypi, PETSC_NULL_INTEGER, xpm, ypm, PETSC_NULL_INTEGER, ierr)
  call DMDAGetGhostCorners(dap, xpgi, ypgi, PETSC_NULL_INTEGER, xpgm, ypgm, PETSC_NULL_INTEGER, ierr)
  xpi=xpi+1
  ypi=ypi+1
  xpe=xpi+xpm-1
  ype=ypi+ypm-1
  xpgi=xpgi+1
  ypgi=ypgi+1
  xpge=xpgi+xpgm-1
  ypge=ypgi+ypgm-1
  
  call DMDAGetGhostCorners(dap_coarsify, xpcgi, ypcgi, PETSC_NULL_INTEGER, xpcgm, ypcgm, PETSC_NULL_INTEGER, ierr) !c for coarsify
  xpcgi=xpcgi+1
  ypcgi=ypcgi+1
  xpcge=xpcgi+xpcgm-1
  ypcge=ypcgi+ypcgm-1

  call DMCreateGlobalVector(dap, pressurevec, ierr)
  call VecDuplicateVecsF90(pressurevec, mgridlev, pvec, ierr)
  call DMGetLocalVector(dap, pvec_local, ierr)
  
  call DMGetLocalVector(dap_coarsify, pressurevec_local, ierr)
  call setup_coarsify_pressure
  
  nbc_p = 2*(m)+2*(n)
  left_phi = 0
  right_phi = n
  bottom_phi = 2*n
  top_phi = 2*n+m
  
  call VecCreateSeq(MPI_COMM_SELF, nbc_p, pressurebcvec_local, ierr)
  call VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nbc_p, pressurebcvec_global, ierr)
  call VecZeroEntries(pressurebcvec_local, ierr)
  call VecZeroEntries(pressurebcvec_global, ierr)
  
  call setup_getbc_phi1  !for coarsest grid
  call setup_getbc_phi2  !for remaining grids
  
  
end subroutine setup_pressure

!======================================================================= 

subroutine destroy_pressure()

  call VecDestroy(pressurevec, ierr)
  call VecDestroyVecsF90(mgridlev, pvec, ierr)
  call DMRestoreLocalVector(dap, pvec_local, ierr)
  call DMDestroy(dap, ierr)

  call DMRestoreLocalVector(dap_coarsify, pressurevec_local, ierr)
  call DMDestroy(dap_coarsify, ierr)
  
  call VecScatterDestroy(ctx_bc_p2, ierr)
  call ISDestroy(bcis_p2, ierr)
  call VecDestroy(pressurebcvec_local, ierr)
  call VecDestroy(pressurebcvec_global, ierr)

end subroutine destroy_pressure

!======================================================================= 

SUBROUTINE calculate_pressure( itime )

  !***************************************************************!
  !*  Multiscale method to solve D D^T pressure = pressure_rhs   *!
  !*  D D^T pressure = D N(q) - D E^T f - d/dt(mdot) + L'mdot    *!
  !***************************************************************!
  use myfft
  use variables
  use operators_fluid
  use user
  use fsinterface
  
  PetscInt :: itime, k, i, j
  Vec :: omegaq_ubcvec_local, omegaq_vbcvec_local, omegaq_ubcvec_global, omegaq_vbcvec_global
  PetscScalar, pointer :: omega_local(:), qx_local(:), qy_local(:), omegaq_ubc(:), omegaq_vbc(:)
  Vec, pointer :: omegaqxvec(:), omegaqyvec(:), nonlinearqxvec(:), nonlinearqyvec(:), pvec_rhs_temp(:)
  PetscScalar, pointer :: omegaqx(:), omegaqy(:), nonlinearqx(:), nonlinearqy(:), fb(:), rhsfsd(:)
  Vec :: omegaqxvec_local, omegaqyvec_local
  PetscScalar :: del
  PetscScalar :: rhsfsdserialarray(nf_sd), divergence_array(xpi:xpe,ypi:ype), pressure_rhs_array(xpi:xpe,ypi:ype)
  PetscScalar, pointer :: pressure(:), p_local(:), pressurebc_array(:), pressure_rhs(:)
  PetscScalar :: omegab, uv(5)

  ! ==================================================================
  ! ================ Nonlinear term N(q)= q x omega ==================
  ! ==================================================================
  
  !temporaray vectors------------------------------
  call VecDuplicate(streambcvec_local, omegaq_ubcvec_local, ierr)  !boundary local
  call VecDuplicate(streambcvec_local, omegaq_vbcvec_local, ierr) !boundary local
  call VecDuplicate(streambcvec_global, omegaq_ubcvec_global, ierr)  !boundary global
  call VecDuplicate(streambcvec_global, omegaq_vbcvec_global, ierr) !boundary global

  call VecDuplicateVecsF90(streamvec, mgridlev, omegaqxvec, ierr)   !these will go into get_bc
  call VecDuplicateVecsF90(streamvec, mgridlev, omegaqyvec, ierr)
  call VecDuplicateVecsF90(velocityxvec, mgridlev, nonlinearqxvec, ierr)   !nonlinear terms
  call VecDuplicateVecsF90(velocityyvec, mgridlev, nonlinearqyvec, ierr)
  
  call VecDuplicate(streamvec_local, omegaqxvec_local, ierr)
  call VecDuplicate(streamvec_local, omegaqyvec_local, ierr)
  
  call VecDuplicateVecsF90(pressurevec, mgridlev, pvec_rhs_temp, ierr)
  !------------------------------------------------
  

  !k = mgridlev  ! === largest domain === !
do k=mgridlev,1,-1
  
  if (k==mgridlev) then
      call VecZeroEntries(omegaq_ubcvec_local, ierr)  ! === largest domain === !
      call VecZeroEntries(omegaq_vbcvec_local, ierr)
  else
      
      call VecGetArrayF90(omegaqxvec_local, omegaqx, ierr)
      call VecGetArrayF90(omegaqyvec_local, omegaqy, ierr)
      call get_bc(omegaqx, omegaq_ubcvec_local, omegaq_ubcvec_global, 1.d0)
      call get_bc(omegaqy, omegaq_vbcvec_local, omegaq_vbcvec_global, 1.d0)
      
      call VecAssemblyBegin(omegaq_ubcvec_global, ierr)
      call VecAssemblyBegin(omegaq_vbcvec_global, ierr)
      call VecAssemblyEnd(omegaq_ubcvec_global, ierr)
      call VecAssemblyEnd(omegaq_vbcvec_global, ierr)
      
      call VecRestoreArrayF90(omegaqxvec_local, omegaqx, ierr)
      call VecRestoreArrayF90(omegaqyvec_local, omegaqy, ierr)
      
      call VecScatterBegin(ctx_bc, omegaq_ubcvec_global, omegaq_ubcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(ctx_bc, omegaq_ubcvec_global, omegaq_ubcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterBegin(ctx_bc, omegaq_vbcvec_global, omegaq_vbcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(ctx_bc, omegaq_vbcvec_global, omegaq_vbcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
      
      
      
  end if
  
  call VecWAXPY(qhxvec, one, qxvec(k), q0xvec(k), ierr)
  call VecWAXPY(qhyvec, one, qyvec(k), q0yvec(k), ierr)
  call DMGlobalToLocalBegin(dau, qhxvec, INSERT_VALUES, qxvec_local, ierr)
  call DMGlobalToLocalBegin(dav, qhyvec, INSERT_VALUES, qyvec_local, ierr)
  call DMGlobalToLocalEnd(dau, qhxvec,    INSERT_VALUES, qxvec_local,    ierr)
  call DMGlobalToLocalEnd(dav, qhyvec,    INSERT_VALUES, qyvec_local,    ierr)  
  
  call VecGetArrayReadF90(omegavec(k), omega_local, ierr)
  call VecGetArrayReadF90(qxvec_local, qx_local, ierr)
  call VecGetArrayReadF90(qyvec_local, qy_local, ierr)
  call VecGetArrayF90(omegaqxvec(k), omegaqx, ierr)
  call VecGetArrayF90(omegaqyvec(k), omegaqy, ierr)  
  call nonlinear_corner(omega_local, qx_local, qy_local, omegaqx, omegaqy)  ! the outputs are omegaqx and omegaqy
  call VecRestoreArrayF90(omegaqxvec(k), omegaqx, ierr)  !will be used for getbc
  call VecRestoreArrayF90(omegaqyvec(k), omegaqy, ierr)
  call VecRestoreArrayReadF90(omegavec(k), omega_local, ierr)
  call VecRestoreArrayReadF90(qxvec_local, qx_local, ierr)
  call VecRestoreArrayReadF90(qyvec_local, qy_local, ierr)
  
  call DMGlobalToLocalBegin(das, omegaqxvec(k), INSERT_VALUES, omegaqxvec_local, ierr)
  call DMGlobalToLocalEnd(das, omegaqxvec(k), INSERT_VALUES, omegaqxvec_local, ierr)
  call DMGlobalToLocalBegin(das, omegaqyvec(k), INSERT_VALUES, omegaqyvec_local, ierr)
  call DMGlobalToLocalEnd(das, omegaqyvec(k), INSERT_VALUES, omegaqyvec_local, ierr)
  
  call VecGetArrayF90(nonlinearqxvec(k), nonlinearqx, ierr)
  call VecGetArrayF90(nonlinearqyvec(k), nonlinearqy, ierr)
  call VecGetArrayF90(omegaqxvec_local, omegaqx, ierr)
  call VecGetArrayF90(omegaqyvec_local, omegaqy, ierr) 
  call VecGetArrayReadF90(omegaq_ubcvec_local, omegaq_ubc, ierr)
  call VecGetArrayReadF90(omegaq_vbcvec_local, omegaq_vbc, ierr)
  
  call nonlinear_corner2edge(omegaqx, omegaqy, omegaq_ubc, omegaq_vbc, nonlinearqx, nonlinearqy)  !nonlinear
  del = delta * REAL( 2**(k-1) )
  nonlinearqx = nonlinearqx/(del**2)
  nonlinearqy = nonlinearqy/(del**2)
  
  call VecRestoreArrayF90(omegaqxvec_local, omegaqx, ierr)
  call VecRestoreArrayF90(omegaqyvec_local, omegaqy, ierr) 
  call VecRestoreArrayF90(nonlinearqxvec(k), nonlinearqx, ierr)
  call VecRestoreArrayF90(nonlinearqyvec(k), nonlinearqy, ierr)
  call VecRestoreArrayReadF90(omegaq_ubcvec_local, omegaq_ubc, ierr)
  call VecRestoreArrayReadF90(omegaq_vbcvec_local, omegaq_vbc, ierr)
  
end do

  !Doing reg to get hxvec and hyvec------------------------
  do k=1,mgridlev
      call VecZeroEntries(pvec_rhs_temp(k), ierr)
  end do
  
  call VecZeroEntries(hxvec, ierr)
  call VecZeroEntries(hyvec, ierr)
  if (.not. sub_domain) then
        call VecGetArrayReadF90(fvec, fb, ierr)
        call reg(fb/dt, hxvec, hyvec)
        call VecAssemblyBegin(hxvec, ierr)
        call VecAssemblyBegin(hyvec, ierr) 
        call VecRestoreArrayReadF90(fvec, fb, ierr)
  else 
        call VecZeroEntries(rhsfsdvec, ierr)
        call VecGetArrayReadF90(fvec, fb, ierr)
        rhsfsdserialarray(interp_indices_local(:)) = matmul(indicesserialmat, fb/dt)
        call VecSetValues(rhsfsdvec, interp_indices_local_size, interp_indices_local-1, rhsfsdserialarray(interp_indices_local(:)), ADD_VALUES, ierr)
        
        call VecAssemblyBegin(rhsfsdvec, ierr)
        call VecRestoreArrayReadF90(fvec, fb, ierr)
        call VecAssemblyEnd(rhsfsdvec, ierr)
        
        call VecGetArrayReadF90(rhsfsdvec, rhsfsd, ierr)
        call reg_sd_sparse(rhsfsd, hxvec, hyvec, size(interp_indices_global), interp_indices_global)
  
        call VecAssemblyBegin(hxvec, ierr)
        call VecAssemblyBegin(hyvec, ierr) 
        call VecRestoreArrayF90(rhsfsdvec, rhsfsd, ierr)
    
  end if
  
  call VecAssemblyEnd(hxvec, ierr)
  call VecAssemblyEnd(hyvec, ierr) 
  !--------------------------------------------
  
  !Adding rhs forcing to hxvec and hyvec
  !Adding is done inside inside rhs_forcing itself to avoid addition of pointer with array
  call VecGetArrayF90(hxvec, nonlinearqx, ierr)
  call VecGetArrayF90(hyvec, nonlinearqy, ierr)
  call VecGetArrayReadF90(qxvec_local, qx_local, ierr)
  call VecGetArrayReadF90(qyvec_local, qy_local, ierr)
  call rhs_forcing(itime, qx_local, qy_local, nonlinearqx, nonlinearqy)
  call VecRestoreArrayF90(hxvec, nonlinearqx, ierr)
  call VecRestoreArrayF90(hyvec, nonlinearqy, ierr)
  call VecRestoreArrayReadF90(qxvec_local, qx_local, ierr)
  call VecRestoreArrayReadF90(qyvec_local, qy_local, ierr)
  !------------------------------------------------
  
  !Computing divergence of all the above forcing
  call DMGlobalToLocalBegin(dau, hxvec, INSERT_VALUES, qxvec_local, ierr)
  call DMGlobalToLocalBegin(dav, hyvec, INSERT_VALUES, qyvec_local, ierr)
  call DMGlobalToLocalEnd(dau, hxvec, INSERT_VALUES, qxvec_local, ierr)
  call DMGlobalToLocalEnd(dav, hyvec, INSERT_VALUES, qyvec_local, ierr)
  
  call VecGetArrayReadF90(qxvec_local, qx_local, ierr)
  call VecGetArrayReadF90(qyvec_local, qy_local, ierr)
  call VecGetArrayF90(pvec_rhs_temp(1), pressure, ierr)
  call divergence(qx_local, qy_local, pressure)
  pressure = -pressure
  call VecRestoreArrayF90(pvec_rhs_temp(1), pressure, ierr)
  call VecRestoreArrayReadF90(qxvec_local, qx_local, ierr)
  call VecRestoreArrayReadF90(qyvec_local, qy_local, ierr)
  !----------------------------------------------------
  
  !coarsify onto bigger domain---------------
  do k=2,mgridlev
      call DMGlobalToLocalBegin(dap_coarsify, pvec_rhs_temp(k-1), INSERT_VALUES, pressurevec_local, ierr)
      call DMGlobalToLocalEnd(dap_coarsify, pvec_rhs_temp(k-1), INSERT_VALUES, pressurevec_local, ierr)
      
      call VecGetArrayReadF90(pressurevec_local, p_local, ierr)
      call coarsify_pressure(pvec_rhs_temp(k), p_local)
      call VecRestoreArrayReadF90(pressurevec_local, p_local, ierr)
  end do
  !---------------------------------------------
  
  !Adding divergence of nonlinear terms
  do k=1,mgridlev
  
      call DMGlobalToLocalBegin(dau, nonlinearqxvec(k), INSERT_VALUES, qxvec_local, ierr)
      call DMGlobalToLocalBegin(dav, nonlinearqyvec(k), INSERT_VALUES, qyvec_local, ierr)
      call DMGlobalToLocalEnd(dau, nonlinearqxvec(k), INSERT_VALUES, qxvec_local, ierr)
      call DMGlobalToLocalEnd(dav, nonlinearqyvec(k), INSERT_VALUES, qyvec_local, ierr)
  
      call VecGetArrayReadF90(qxvec_local, qx_local, ierr)
      call VecGetArrayReadF90(qyvec_local, qy_local, ierr)
      call VecGetArrayF90(pvec_rhs_temp(k), pressure, ierr)
      call divergence(qx_local, qy_local, divergence_array)
      call add_pressure_divergence(pressure, divergence_array)
      !pressure = pressure + divergence_array
      call VecRestoreArrayF90(pvec_rhs_temp(k), pressure, ierr)
      call VecRestoreArrayReadF90(qxvec_local, qx_local, ierr)
      call VecRestoreArrayReadF90(qyvec_local, qy_local, ierr)
  
  end do
  !---------------------------------------------
  
  ! =====================================================================
  ! ================ Solve for pressure =================================
  ! =====================================================================
  ! ==== DD^T (0.5*|u|^2 + P + 0.5*omegab^2*|x|^2) = RHS           ==== !
  ! == At boundary, P=0 --> BC = 0.5*|q0/del|^2 -0.5*omegab^2*|x|^2  == !
  
  !on the coarsest grid----------------------
  do k=mgridlev,1,-1  
  !k           = mgridlev
  
      !get_bc-----------
      if (k==mgridlev) then
  
          del         = 0.125D0 / ( delta * (2.d0**(k-1)) )**2
          call VecZeroEntries(pressurebcvec_local, ierr)
          call VecZeroEntries(pressurebcvec_global, ierr)
          ! user-defined angular velocity of the grid
          uv = motion_grid(itime)
          omegab = uv(3)
  
          !get_bc-----------
          call VecWAXPY(qhxvec, one, qxvec(k), q0xvec(k), ierr)
          call VecWAXPY(qhyvec, one, qyvec(k), q0yvec(k), ierr)
          call DMGlobalToLocalBegin(dau, qhxvec, INSERT_VALUES, qxvec_local, ierr)
          call DMGlobalToLocalBegin(dav, qhyvec, INSERT_VALUES, qyvec_local, ierr)
          call DMGlobalToLocalEnd(dau, qhxvec, INSERT_VALUES, qxvec_local, ierr)
          call DMGlobalToLocalEnd(dav, qhyvec, INSERT_VALUES, qyvec_local, ierr)
  
          call VecGetArrayReadF90(qxvec_local, qx_local, ierr)
          call VecGetArrayReadF90(qyvec_local, qy_local, ierr)
          call get_bc_phi1(qx_local, qy_local, pressurebcvec_local, del)
  
          call VecAssemblyBegin(pressurebcvec_local, ierr)
          call VecRestoreArrayReadF90(qxvec_local, qx_local, ierr)
          call VecRestoreArrayReadF90(qyvec_local, qy_local, ierr)
          call VecAssemblyEnd(pressurebcvec_local, ierr)
          
      else
      
          call DMGlobalToLocalBegin(dap, pvec(k+1), INSERT_VALUES, pvec_local, ierr)
          call DMGlobalToLocalEnd(dap, pvec(k+1), INSERT_VALUES, pvec_local, ierr)
          call VecGetArrayReadF90(pvec_local, pressure, ierr)
          call get_bc_phi2(pressure, pressurebcvec_local, pressurebcvec_global, 1.d0)
          call VecRestoreArrayReadF90(pvec_local, pressure, ierr)
          
          call VecAssemblyBegin(pressurebcvec_global, ierr)
          call VecAssemblyEnd(pressurebcvec_global, ierr)  
          call VecScatterBegin(ctx_bc_p2, pressurebcvec_global, pressurebcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecScatterEnd(ctx_bc_p2, pressurebcvec_global, pressurebcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
  
      end if
  
      !apply_bc---------
      call VecGetArrayF90(pvec_rhs_temp(k), pressure_rhs, ierr)
      call VecGetArrayF90(pressurebcvec_local, pressurebc_array, ierr)
      call apply_bc_phi(pressure_rhs, pressurebc_array)
      call VecRestoreArrayF90(pressurebcvec_local, pressurebc_array, ierr)
  
      !ddti------------
      call VecGetArrayF90(pvec(k), pressure, ierr)
      call ddti(xpi, xpe, ypi, ype, pressure_rhs, pressure)
      call VecRestoreArrayF90(pvec(k), pressure, ierr)
      call VecRestoreArrayF90(pvec_rhs_temp(k), pressure_rhs, ierr)
      !----------------------------------------------
      
  end do
  !solved for pressure==============================================
  
  !pressure_stag = pressure ! pressure_stag = pressure_static + pressure_dynamic [== 0.5 rho |u|^2 ]--------
  do k=1,mgridlev
      del = 0.125D0 / ( delta * (2.d0**(k-1)) )**2
      call VecWAXPY(qhxvec, one, qxvec(k), q0xvec(k), ierr)
      call VecWAXPY(qhyvec, one, qyvec(k), q0yvec(k), ierr)
      call DMGlobalToLocalBegin(dau, qhxvec, INSERT_VALUES, qxvec_local, ierr)
      call DMGlobalToLocalBegin(dav, qhyvec, INSERT_VALUES, qyvec_local, ierr)
      call DMGlobalToLocalEnd(dau, qhxvec, INSERT_VALUES, qxvec_local, ierr)
      call DMGlobalToLocalEnd(dav, qhyvec, INSERT_VALUES, qyvec_local, ierr)
      
      call VecGetArrayReadF90(qxvec_local, qx_local, ierr)
      call VecGetArrayReadF90(qyvec_local, qy_local, ierr)
      call VecGetArrayF90(pvec(k), pressure, ierr)
      call static_plus_dynamic(pressure, qx_local, qy_local, del)
      call VecRestoreArrayReadF90(qxvec_local, qx_local, ierr)
      call VecRestoreArrayReadF90(qyvec_local, qy_local, ierr)
      
      ! scale pressure to the pressure coefficient
      pressure = 2.d0*pressure
      call VecRestoreArrayF90(pvec(k), pressure, ierr)
  
  end do
  
  !destroy temporary vectors---------------------
  call VecDestroy(omegaq_ubcvec_local, ierr)
  call VecDestroy(omegaq_vbcvec_local, ierr)
  call VecDestroy(omegaq_ubcvec_global, ierr)
  call VecDestroy(omegaq_vbcvec_global, ierr)
  
  call VecDestroyVecsF90(mgridlev, omegaqxvec, ierr)
  call VecDestroyVecsF90(mgridlev, omegaqyvec, ierr)
  call VecDestroyVecsF90(mgridlev, nonlinearqxvec, ierr)
  call VecDestroyVecsF90(mgridlev, nonlinearqyvec, ierr)
  
  call VecDestroy(omegaqxvec_local, ierr)
  call VecDestroy(omegaqyvec_local, ierr)
  
  call VecDestroyVecsF90(mgridlev, pvec_rhs_temp, ierr)
  !----------------------------------------------
  
  
  
END SUBROUTINE calculate_pressure

!======================================================================= 

subroutine nonlinear_corner(omega, qx, qy, omegaqx, omegaqy)

  PetscScalar :: omega(xsi:xse,ysi:yse), qx(xugi:xuge,yugi:yuge), qy(xvgi:xvge,yvgi:yvge), omegaqx(xsi:xse,ysi:yse), omegaqy(xsi:xse,ysi:yse)
  PetscInt :: i,j
  
  do j=ysi,yse
      do i=xsi,xse
          omegaqx(i,j) =  0.5d0*( qy(i,j)+qy(i-1,j) )*omega(i,j)
          omegaqy(i,j) = -0.5d0*( qx(i,j)+qx(i,j-1) )*omega(i,j)
      end do
  end do
  
  

end subroutine nonlinear_corner

!======================================================================= 

subroutine nonlinear_corner2edge(omegaqx, omegaqy, omegaq_ubc, omegaq_vbc, nonlinearqx, nonlinearqy)

  PetscScalar :: omegaqx(xsgi:xsge,ysgi:ysge), omegaqy(xsgi:xsge,ysgi:ysge), omegaq_ubc(nbc), omegaq_vbc(nbc), nonlinearqx(xui:xue,yui:yue), nonlinearqy(xvi:xve,yvi:yve)
  PetscInt :: i,j
  
  if (xsgi==1)   then
      forall(i=xsgi:xsgi, j=ysgi:ysge);  
          omegaqx(i,j)=omegaq_ubc(left+j);   
          omegaqy(i,j)=omegaq_vbc(left+j);   
      end forall
  end if
  if (xsge==m+1) then 
      forall(i=xsge:xsge, j=ysgi:ysge)  
          omegaqx(i,j)=omegaq_ubc(right+j);   
          omegaqy(i,j)=omegaq_vbc(right+j);
      end forall 
  end if
  if (ysgi==1) then 
      forall(i=xsgi:xsge, j=ysgi:ysgi)  
          omegaqx(i,j)=omegaq_ubc(bottom+i) 
          omegaqy(i,j)=omegaq_vbc(bottom+i) 
      end forall 
  end if
  if (ysge==n+1) then 
      forall(i=xsgi:xsge, j=ysge:ysge)  
          omegaqx(i,j)=omegaq_ubc(top+i) 
          omegaqy(i,j)=omegaq_vbc(top+i) 
      end forall 
  end if
  
  do j=yui,yue
      do i=xui,xue
        nonlinearqx(i,j)=0.5d0*( omegaqx(i,j+1) + omegaqx(i,j) )
      end do
  end do
    
  do j=yvi,yve
      do i=xvi,xve
        nonlinearqy(i,j)=0.5d0*( omegaqy(i+1,j) + omegaqy(i,j) )
      end do
  end do

end subroutine nonlinear_corner2edge

!===================================================================================

  subroutine divergence(qx, qy, pressure)
  
    PetscScalar :: qx(xugi:xuge,yugi:yuge), qy(xvgi:xvge,yvgi:yvge), pressure(xpi:xpe,ypi:ype)
    PetscInt :: i, j
    
    do j=ypi,ype
      do i=xpi,xpe
        pressure(i,j) = qx(i+1,j) - qx(i,j) + qy(i,j+1) - qy(i,j)
      end do
    end do
  
  end subroutine divergence
  
!===================================================================================

  subroutine add_pressure_divergence(pressure, divergence_array)
  
    PetscScalar :: pressure(xpi:xpe,ypi:ype), divergence_array(xpi:xpe,ypi:ype)
    
    pressure = pressure + divergence_array
    
  end subroutine add_pressure_divergence
!===================================================================================

  subroutine coarsify_pressure(pvec2, p1)
  
    Vec :: pvec2
    PetscScalar :: p1(xpcgi:xpcge,ypcgi:ypcge), temp(xpcgi:xpcge,ypcgi:ypcge), yy(tcoarse_p)
    PetscInt :: i,j, ix(tcoarse_p), id, jd, iter
    
    iter = 0
    do j=yci_p,min(ypcge-1,n-2),2
        do i=xci_p,min(xpcge-1,m-2),2
            temp(i,j) =        p1(i,j)   + &
                      0.5d0*(  p1(i+1,j)   + p1(i,j+1)   + &
                               p1(i-1,j)   + p1(i,j-1) ) + &
                      0.25d0*( p1(i+1,j+1) + p1(i+1,j-1)   + &
                               p1(i-1,j-1) + p1(i-1,j+1) )
            iter = iter+1
        end do
    end do
    
    iter = 0
    do j=yci_p,yce_p,2
        do i=xci_p,xce_p,2
            iter = iter+1
            yy(iter) = (temp(i,j)   + 3.D0*temp(i+2,j) + &
                     3.D0*temp(i,j+2) + 9.D0*temp(i+2,j+2) )/16.D0
            id = m/4+i/2+1
            jd = n/4+j/2+1
            ix(iter) = m*(jd-1) + (id-1)
        end do
    end do
    call VecSetValues(pvec2, tcoarse_p, ix, yy, INSERT_VALUES, ierr)
    
    call VecAssemblyBegin(pvec2, ierr)
    call VecAssemblyEnd(pvec2, ierr)   
  
  end subroutine coarsify_pressure
  
!===================================================================================

  subroutine setup_coarsify_pressure
  
    PetscInt :: i,j
  
    !starting indices
    if (modulo(xpi,2)==0) then
      xci_p=xpi
    else
      xci_p=xpi+1
    end if
    
    if (modulo(ypi,2)==0) then
      yci_p=ypi
    else
      yci_p=ypi+1
    end if
    
    !ending indices
    if (modulo(xpe,2)==0) then
      xce_p=xpe
    else
      xce_p=xpe-1
    end if
    
    if (modulo(ype,2)==0) then
      yce_p=ype
    else
      yce_p=ype-1
    end if
    
    !some end conditions on ending indices of top and right boundaries
    if (xpe==m) xce_p = m-4
    if (ype==n) yce_p = n-4
    
    tcoarse_p = 0
    do j=yci_p,yce_p,2
        do i=xci_p,xce_p,2
            tcoarse_p = tcoarse_p+1
        end do
    end do
    
  end subroutine setup_coarsify_pressure

!===================================================================================
  
subroutine setup_getbc_phi1

  !Setting up getbc for the coarsest grid

  PetscInt :: i, j, next, tis
  
  !Finding indices of boundary points on fine grid (bcix) for scattering from global bc vec to local bc vec
    next=0
    do i=xpi,xpe
      do j=ypi,ype
      
        if (i==1) then
          next=next+1
          call grow_array(cbcix_p1, left_phi+j-1, next)
        end if
        
        if (i==m) then
          next=next+1
          call grow_array(cbcix_p1, right_phi+j-1, next)
        end if
        
        if (j==1) then
          next=next+1
          call grow_array(cbcix_p1, bottom_phi+i-1, next)
        end if
        
        if (j==n) then
          next=next+1
          call grow_array(cbcix_p1, top_phi+i-1, next)
        end if
       
      end do
    end do
    tbc_p1=next

end subroutine setup_getbc_phi1

!===================================================================================
  
subroutine setup_getbc_phi2

  !Setting up getbc for the coarsest grid

  PetscInt :: i, j, next, tis
  
  !Finding indices of boundaries on coarse grid (cbcix) that needs to be transferred to the fine grid
    next=0
    do i=xpi,xpe
      do j=ypi,ype
      
        if (j==n/4   .and. i>=m/4+1 .and. i<=3*m/4) then
          next=next+1
          call grow_array(cbcix_p2, bottom_phi+2*(i-1-m/4), next)   !-1 from i for C indexing
          next=next+1
          call grow_array(cbcix_p2, bottom_phi+2*(i-1-m/4)+1, next)   !-1 from i for C indexing
          
        end if
        
        if (j==3*n/4+1 .and. i>=m/4+1 .and. i<=3*m/4) then
          next=next+1
          call grow_array(cbcix_p2, top_phi+2*(i-1-m/4), next)   !-1 from i for C indexing
          next=next+1
          call grow_array(cbcix_p2, top_phi+2*(i-1-m/4)+1, next)   !-1 from i for C indexing
        end if
        
        if (i==m/4   .and. j>=n/4+1 .and. j<=3*n/4) then
          next=next+1
          call grow_array(cbcix_p2, left_phi+2*(j-1-n/4), next)   !-1 from j for C indexing
          next=next+1
          call grow_array(cbcix_p2, left_phi+2*(j-1-n/4)+1, next)   !-1 from j for C indexing
        end if
        
        if (i==3*m/4+1 .and. j>=n/4+1 .and. j<=3*n/4) then
          next=next+1
          call grow_array(cbcix_p2, right_phi+2*(j-1-n/4), next)   !-1 from j for C indexing
          next=next+1
          call grow_array(cbcix_p2, right_phi+2*(j-1-n/4)+1, next)   !-1 from j for C indexing
        end if
        
      end do
    end do
    
    tbc_p2=next  !total number of boundary points that needs to be transferred
  
  !Finding indices of boundary points on fine grid (bcix) for scattering from global bc vec to local bc vec
    next=0
    do i=xpi,xpe
      do j=ypi,ype
      
        if (i==1) then
          next=next+1
          call grow_array(bcix_p2, left_phi+j-1, next)
        end if
        
        if (i==m) then
          next=next+1
          call grow_array(bcix_p2, right_phi+j-1, next)
        end if
        
        if (j==1) then
          next=next+1
          call grow_array(bcix_p2, bottom_phi+i-1, next)
        end if
        
        if (j==n) then
          next=next+1
          call grow_array(bcix_p2, top_phi+i-1, next)
        end if
       
      end do
    end do
    tis=next
    call ISCreateGeneral(MPI_COMM_SELF, tis, bcix_p2, PETSC_COPY_VALUES, bcis_p2, ierr)
    call VecScatterCreate(pressurebcvec_global, bcis_p2, pressurebcvec_local, bcis_p2, ctx_bc_p2, ierr)   !Scatter context for get_bc

end subroutine setup_getbc_phi2

!===================================================================================
  
subroutine get_bc_phi1(qx, qy, pbcvec_local, del)

  PetscScalar :: qx(xugi:xuge,yugi:yuge), qy(xvgi:xvge,yvgi:yvge)
  Vec :: pbcvec_local  
  PetscInt :: i, j, next
  PetscScalar :: yy(tbc_p1), del
  
  !Finding indices of boundary points on fine grid (bcix) for scattering from global bc vec to local bc vec
    next=0
    do i=xpi,xpe
      do j=ypi,ype
      
        if (i==1) then
          next=next+1
          yy(next) = del * ( (qx(i+1,j)+qx(i,j))**2 + (qy(i,j+1)+qy(i,j))**2 )
        end if
        
        if (i==m) then
          next=next+1
          yy(next) = del * ( (qx(i+1,j)+qx(i,j))**2 + (qy(i,j+1)+qy(i,j))**2 )
        end if
        
        if (j==1) then
          next=next+1
          yy(next) = del * ( (qx(i+1,j)+qx(i,j))**2 + (qy(i,j+1)+qy(i,j))**2 )
        end if
        
        if (j==n) then
          next=next+1
          yy(next) = del * ( (qx(i+1,j)+qx(i,j))**2 + (qy(i,j+1)+qy(i,j))**2 )
        end if
       
      end do
    end do
    call VecSetValues(pbcvec_local, tbc_p1, cbcix_p1, yy, INSERT_VALUES, ierr)

end subroutine get_bc_phi1

!===================================================================================
  
subroutine get_bc_phi2(p_local, pbcvec_local, pbcvec_global, fac)

  PetscScalar :: p_local(xpgi:xpge,ypgi:ypge), fac, yy(tbc_p2)
  Vec :: pbcvec_local, pbcvec_global
  PetscInt :: i, j, next
  
  !Finding indices of boundaries on coarse grid (cbcix) that needs to be transferred to the fine grid
    next=0
    do i=xpi,xpe
      do j=ypi,ype
      
        if (j==n/4   .and. i>=m/4+1 .and. i<=3*m/4) then
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i,j+1) + p_local(i-1,j))
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i,j+1) + p_local(i+1,j))
          
        end if
        
        if (j==3*n/4+1 .and. i>=m/4+1 .and. i<=3*m/4) then
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i,j-1) + p_local(i-1,j))
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i,j-1) + p_local(i+1,j))
        end if
        
        if (i==m/4   .and. j>=n/4+1 .and. j<=3*n/4) then
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i+1,j) + p_local(i,j-1))
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i+1,j) + p_local(i,j+1))
        end if
        
        if (i==3*m/4+1 .and. j>=n/4+1 .and. j<=3*n/4) then
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i-1,j) + p_local(i,j-1))
          next=next+1
          yy(next) = 0.25*(2*p_local(i,j) + p_local(i-1,j) + p_local(i,j+1))
        end if
        
      end do
    end do
    
    yy=yy*fac
    call VecSetValues(pbcvec_global, tbc_p2, cbcix_p2, yy, INSERT_VALUES, ierr)
  
end subroutine get_bc_phi2

!===================================================================================

subroutine apply_bc_phi(pressure, pbc)

    !*****************************************************************!
    !*   given phi right outside of domain, phi_bc, (from larger,    *!
    !*   coraser mesh), add values to correct laplacian(D*D^T)of phi *!
    !*   , mass_rhs   on the (smaller, finer) domain, r.             *!
    !*****************************************************************!

  PetscScalar :: pressure(xpi:xpe,ypi:ype), pbc(1:nbc_p)
  PetscInt :: i, j
  
  !print*, rank, pressure
  if (xpi==1) then; i=1; forall(j=ypi:ype);  pressure(i,j)=pressure(i,j)-pbc(left_phi+j);   end forall; end if;
  if (xpe==m) then; i=m; forall(j=ypi:ype);  pressure(i,j)=pressure(i,j)-pbc(right_phi+j);   end forall; end if;
  if (ypi==1) then; j=1; forall(i=xpi:xpe);  pressure(i,j)=pressure(i,j)-pbc(bottom_phi+i);   end forall; end if;
  if (ype==n) then; j=n; forall(i=xpi:xpe);  pressure(i,j)=pressure(i,j)-pbc(top_phi+i);   end forall; end if; 
  
end subroutine apply_bc_phi  

!===================================================================================

subroutine static_plus_dynamic(pressure, qx, qy, del)

  PetscScalar :: pressure(xpi:xpe,ypi:ype), qx(xugi:xuge,yugi:yuge), qy(xvgi:xvge,yvgi:yvge), del
  PetscInt :: i, j
  
  do j=ypi,ype
      do i=xpi,xpe
          pressure(i,j) = pressure(i,j) - del * ( (qx(i+1,j)+qx(i,j))**2 + (qy(i,j+1)+qy(i,j))**2 )  !have not added rotation here
      end do
  end do

end subroutine static_plus_dynamic

!===================================================================================

subroutine write_pressure(it)

  PetscInt :: it, k
  character(2) :: charrank
  character(7) :: charit
  PetscScalar, pointer :: readvar(:)
  PetscScalar :: pressure_write(xpm*ypm)
  
  WRITE(charit,"(I7.7)") it
  WRITE(charrank,"(I2.2)") rank
  OPEN(unit=(rank+1)*100,file="output/pressure"//charrank//"_"//charit//".var",form="unformatted",status="unknown")
  
  WRITE((rank+1)*100) xpi, xpe, xpm, ypi, ype, ypm
  do k=1,mgridlev
      call VecGetarrayReadF90(pvec(k), readvar, ierr)
      pressure_write=readvar
      WRITE((rank+1)*100) pressure_write
      call VecRestorearrayReadF90(pvec(k), readvar, ierr)
  end do
  
  CLOSE((rank+1)*100)
  
end subroutine write_pressure
!===================================================================================
 
END MODULE operators_pressure 
