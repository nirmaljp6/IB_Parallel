module timestepper

  USE petsc
#include "petsc/finclude/petsc.h"
  USE parameters
  USE myfft
  USE grid
  USE variables
  USE operators_fluid
  USE operators_struct
  use fsinterface
  IMPLICIT NONE
  
contains  

!===================================================================================

  subroutine advance(itime)
  
    PetscInt :: itime, k
    PetscScalar, pointer :: omega_local(:), qx_local(:), qy_local(:), rhsfluid(:)
    PetscScalar, pointer :: lastbc_local(:), nlx(:), nly(:), omegabc_local(:)
    PetscScalar :: rhs(xsi:xse,ysi:yse), rhsbc(xsi:xse,ysi:yse)
    PetscScalar, pointer :: omega(:)
    PetscScalar, pointer :: velx_interface(:), vely_interface(:), rhsf(:)
    PetscScalar, pointer :: fb(:), hx(:), hy(:), vort(:), rhsfsd(:)
    PetscViewer :: viewer
    
    PetscScalar :: err_fsi, err_fsi_glob, emax, emin
    PetscInt :: fsi_count
    
    PetscInt :: ipiv(nf), info, j
    PetscScalar :: rhsfserial(xfm,1), rhsfsdserialarray(nf_sd)
    
            Mat :: a_reg, a_reg_body, dummat, dummat2, dummat3
    
    itime=itime+1
    
    IF (rank==0 .and. MOD(itime,10).eq.0 ) THEN
       WRITE(*,*) "...Advancing to itime =",itime
    END IF
    
    call motion_potential(itime-1, q0pxvec, q0pyvec)
    call motion_rotation(itime-1, q0rxvec, q0ryvec) 
    do k=1,mgridlev
      call VecWAXPY(q0xvec(k), one, q0pxvec(k), q0rxvec(k), ierr)
      call VecWAXPY(q0yvec(k), one, q0pyvec(k), q0ryvec(k), ierr)
    end do
    rot_angle=rot_angle + delta_angle(itime-1)
    
    if (itime==istart+1) call MatShellSetOperation(AB, MATOP_MULT, ab_times, ierr)

    !For 2D domain decomposition and using CG
    !Generate linear operators
    !IF (itime==istart+1) THEN
    !    call setup_linear_operators    
    !    call setup_constants        
    !    do k=1, mgridlev
    !        call MatDuplicate(Imat, MAT_SHARE_NONZERO_PATTERN, Amat(k)%matrix, ierr)
    !        call MatCopy(Imat, Amat(k)%matrix, SAME_NONZERO_PATTERN, ierr)
    !        call MatAXPY(Amat(k)%matrix, pvfac(k), CTCmat, SAME_NONZERO_PATTERN, ierr)            
    !        call MatDuplicate(Imat, MAT_SHARE_NONZERO_PATTERN, Amat2(k)%matrix, ierr)
    !        call MatCopy(Imat, Amat2(k)%matrix, SAME_NONZERO_PATTERN, ierr)
    !        call MatAXPY(Amat2(k)%matrix, mvfac(k), CTCmat, SAME_NONZERO_PATTERN, ierr)
    !        
    !        call KSPSetTolerances(vort_solver(k)%ksp, rtol, atol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER,ierr)    
    !        call KSPSetOperators(vort_solver(k)%ksp, Amat(k)%matrix, Amat(k)%matrix, ierr)
    !        call KSPSetUp(vort_solver(k)%ksp, ierr)
    !    end do       
    !    call KSPSetTolerances(ksp_sf, rtol, atol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER,ierr)    
    !    call KSPSetOperators(ksp_sf, CTCmat, CTCmat, ierr)
    !    call KSPSetUp(ksp_sf, ierr)        
    !    !call MatView(CTCmat, PETSC_VIEWER_STDOUT_WORLD, ierr)        
    !END IF
    
    IF (itime==1 .and. num_stat) THEN
      CALL preprocess !Precompute and perform Cholesky decomposition
    END IF
    IF (itime==istart+1 .and. sub_domain) THEN
      
      CALL preprocess_sd  !Precompute B matrix
      
      call VecGetArrayReadF90(xbvec, xb, ierr)
      call VecGetArrayReadF90(ybvec, yb, ierr)
      if (.not. sub_domain_precomputed) then
            call reg_prep_sd_sparse(xb, yb)
      else 
            call setup_reg_subdomain_indices(xb, yb) !creates new reg_ixs_sd_sparse and interfacex_sd_is_sparse
            call VecCreateSeq(MPI_COMM_SELF, interfacex_sd_sparse, velxvec_interface_sd_sparse, ierr)
            call VecCreateSeq(MPI_COMM_SELF, interfacey_sd_sparse, velyvec_interface_sd_sparse, ierr)
            call VecScatterCreate(velocityxvec, reg_ixs_sd_sparse, velxvec_interface_sd_sparse, interfacex_sd_is_sparse, ctx_velx_sd_sparse, ierr)
            call VecScatterCreate(velocityyvec, reg_iys_sd_sparse, velyvec_interface_sd_sparse, interfacey_sd_is_sparse, ctx_vely_sd_sparse, ierr)
      end if
      call VecRestoreArrayReadF90(xbvec, xb, ierr)
      call VecRestoreArrayReadF90(ybvec, yb, ierr)
               
    END IF
    !---------------------------------------------------------------
    
    ! Step 1: solve intermediate curl(momentum) eq. on each grid------------------------
    
    !Getting boundary conditions of previous time step-------------
    do k=1, mgridlev-1
      if (k==1) call DMGlobalToLocalBegin(das, omegavec(k+1), INSERT_VALUES, omegavec_local(k+1), ierr)
      call DMGlobalToLocalEnd(das, omegavec(k+1), INSERT_VALUES, omegavec_local(k+1), ierr)
      if (k<mgridlev-1) call DMGlobalToLocalBegin(das, omegavec(k+2), INSERT_VALUES, omegavec_local(k+2), ierr)
      
      if (k>1) call VecAssemblyEnd(streambcvecs_global(k-1), ierr)
      if (k>1) call VecScatterBegin(ctx_bc, streambcvecs_global(k-1), lastbcvec_local(k-1), INSERT_VALUES, SCATTER_FORWARD, ierr)
      
      call VecGetArrayReadF90(omegavec_local(k+1), omega_local, ierr)
      call get_bc(omega_local, lastbcvec_local(k), streambcvecs_global(k), 0.25d0)
      if (k>1) call VecScatterEnd(ctx_bc, streambcvecs_global(k-1), lastbcvec_local(k-1), INSERT_VALUES, SCATTER_FORWARD, ierr) 
      call VecAssemblyBegin(streambcvecs_global(k), ierr)
      call VecRestoreArrayReadF90(omegavec_local(k+1), omega_local, ierr)
    end do
    call VecZeroEntries(lastbcvec_local(mgridlev), ierr)
    !---------------------------------------------------------------
    
    !Computing vorticity and streamfunction predictions at all grid levels---------------------------------------------
    do k=mgridlev,1,-1
    
        if (k==mgridlev) call VecAssemblyEnd(streambcvecs_global(k-1), ierr)
        if (k==mgridlev) call VecScatterBegin(ctx_bc, streambcvecs_global(k-1), lastbcvec_local(k-1), INSERT_VALUES, SCATTER_FORWARD, ierr)
      
        !Initiate communications of boundary points for computing nonlinear terms later
        call VecWAXPY(qhxvec, one, qxvec(k), q0xvec(k), ierr)
        call VecWAXPY(qhyvec, one, qyvec(k), q0yvec(k), ierr)
        call DMGlobalToLocalBegin(dau, qhxvec, INSERT_VALUES, qxvec_local, ierr)
        call DMGlobalToLocalBegin(dav, qhyvec, INSERT_VALUES, qyvec_local, ierr)
        
        if (k==mgridlev) call VecScatterEnd(ctx_bc, streambcvecs_global(k-1), lastbcvec_local(k-1), INSERT_VALUES, SCATTER_FORWARD, ierr) 
      
        !Getting boundary conditions for the current time step from the coarser grid-------
        if (k .lt. mgridlev) then
            !we do not need to do dmglobaltolocal again here because the previous k iteration took care of it
            call DMGlobalToLocalEnd(das, omegavec(k+1), INSERT_VALUES, omegavec_local(k+1), ierr)
            call DMGlobalToLocalBegin(das, omegavec(k), INSERT_VALUES, omegavec_local(k), ierr)  !begin a different communication for the computation of the nonlinear term. Note that this routine is placed here only to optimize communication and is not needed by get_bc
            call VecGetArrayReadF90(omegavec_local(k+1), omega_local, ierr)
            call get_bc(omega_local, omegabcvec_local(k), streambcvecs_global(k), 0.25d0)
            call VecAssemblyBegin(streambcvecs_global(k), ierr)
            call VecRestoreArrayReadF90(omegavec_local(k+1), omega_local, ierr)
        else
            call DMGlobalToLocalBegin(das, omegavec(k), INSERT_VALUES, omegavec_local(k), ierr)  !begin a different communication for the computation of the nonlinear term. Note that this routine is placed here only to optimize communication and is not needed by get_bc
            call VecZeroEntries(omegabcvec_local(mgridlev), ierr)
        end if
      !---------------------------------------------------------------
      
      !Computing the nonliner term (nonlinear term mult. by dt/2)-----------------------
        call DMGlobalToLocalEnd(das, omegavec(k), INSERT_VALUES, omegavec_local(k), ierr)
        call DMGlobalToLocalEnd(dau, qhxvec,    INSERT_VALUES, qxvec_local,    ierr)
        call DMGlobalToLocalEnd(dav, qhyvec,    INSERT_VALUES, qyvec_local,    ierr)
      
        call VecGetArrayReadF90(omegavec_local(k), omega_local, ierr)
        call VecGetArrayReadF90(qxvec_local, qx_local, ierr)
        call VecGetArrayReadF90(qyvec_local, qy_local, ierr)
        call VecGetArrayReadF90(lastbcvec_local(k), lastbc_local, ierr) !this lastbc_local is used again in apply_bc
        call VecGetArrayF90(nlxvec, nlx, ierr)
        call VecGetArrayF90(nlyvec, nly, ierr)
        
        if (k.lt.mgridlev) call VecAssemblyEnd(streambcvecs_global(k), ierr)
        if (k.lt.mgridlev) call VecScatterBegin(ctx_bc, streambcvecs_global(k), omegabcvec_local(k), INSERT_VALUES, SCATTER_FORWARD, ierr)
      
        !Computing nonlinear here--
        call nonlinear(omega_local, qx_local, qy_local, lastbc_local, nlx, nly)
      !-----------------------
      !ADD user-defined RHS forcing term to momentum eq. 
      !Adding is done inside inside rhs_forcing itself to avoid addition of pointer with array
        if (k==1) then
            !call rhs_forcing(itime, qx_local, qy_local, nlx, nly)
        end if
      !-----------------------
      
        call VecRestoreArrayReadF90(omegavec_local(k), omega_local, ierr)
        call VecRestoreArrayReadF90(qxvec_local, qx_local, ierr)
        call VecRestoreArrayReadF90(qyvec_local, qy_local, ierr)
        call VecRestoreArrayF90(nlxvec, nlx, ierr)
        call VecRestoreArrayF90(nlyvec, nly, ierr)
        call DMGlobalToLocalBegin(dau, nlxvec, INSERT_VALUES, nlxvec_local, ierr)  !Initiate communications of boundary points
        call DMGlobalToLocalBegin(dav, nlyvec, INSERT_VALUES, nlyvec_local, ierr)
      !nonlinear ends----------------------------------------------------------------------------------
           
      !Applying boundary condition into the rhs---apply_bc to compute rhsbc(independent of the output of nonlinear) --------------------------------------
        rhsbc = 0.d0
        if (k.lt.mgridlev) call VecScatterEnd(ctx_bc, streambcvecs_global(k), omegabcvec_local(k), INSERT_VALUES, SCATTER_FORWARD, ierr) 
        call VecGetArrayReadF90(omegabcvec_local(k), omegabc_local, ierr)
      
        call apply_bc(rhsbc, lastbc_local,  vfac(k))
        call apply_bc(rhsbc, omegabc_local, vfac(k))
      
        call VecRestoreArrayReadF90(omegabcvec_local(k), omegabc_local, ierr)
        call VecRestoreArrayReadF90(lastbcvec_local(k), lastbc_local, ierr)
      !----------------------------------------------------------------------------------
      
        !Transpose of curl operation on nonlinear term via rot and computing rhs and rhs_old--------------------------------------------  
        !First perform on the interior terms
        call VecGetArrayReadF90(nlxvec, nlx, ierr)
        call VecGetArrayReadF90(nlyvec, nly, ierr)
        call rot_int(nlx, nly, rhs) !do rot
        call VecRestoreArrayReadF90(nlxvec, nlx, ierr)
        call VecRestoreArrayReadF90(nlyvec, nly, ierr)

        !Now perform curl on the boundary terms after ending communication
        call DMGlobalToLocalEnd(dau, nlxvec, INSERT_VALUES, nlxvec_local, ierr)
        call VecGetArrayReadF90(nlxvec_local, nlx, ierr)
        call DMGlobalToLocalEnd(dav, nlyvec, INSERT_VALUES, nlyvec_local, ierr)
        call VecGetArrayReadF90(nlyvec_local, nly, ierr)
        call rot_bc(nlx, nly, rhs) !do rot
        call VecRestoreArrayReadF90(nlxvec_local, nlx, ierr)
        call VecRestoreArrayReadF90(nlyvec_local, nly, ierr)

        ! If this is the very first time step, we need to use explicit Euler
        if (itime==1) THEN
            rhs_old(:,:,k) = rhs(:,:)
        end if
        !----------------------------------------------------------------------------------
      
        !For 2D domain decompositionaggregating rhs-----------------------------------------------------------
        !call VecGetArrayF90(rhsfluid_vec, rhsfluid, ierr)
        !call aggregate_rhs(rhsfluid, rhs, rhs_old, rhsbc, k)
        !call VecRestoreArrayF90(rhsfluid_vec, rhsfluid, ierr)
        !call MatMultAdd(Amat2(k)%matrix, omegavec(k), rhsfluid_vec, rhsfluid_vec, ierr)
      
        !!computing fluid solution----------------------------------------
        call VecGetArrayF90(omegavec(k), omega, ierr)
        !call KSPSolve(vort_solver(k)%ksp, rhsfluid_vec, omegavec(k), ierr) !For CG
        call fluid_predictor(omega, rhs, rhs_old, rhsbc, k) 
        call VecRestoreArrayF90(omegavec(k), omega, ierr)
        !----------------------------------------------------------------------------------
      
        if (k .gt. 1) call DMGlobalToLocalBegin(das, omegavec(k), INSERT_VALUES, omegavec_local(k), ierr)
      
        !update for next time step
        rhs_old(:,:,k) = rhs(:,:)

    end do  
    
    !Adjust vorticity via interpolation and compute predictions for velocity and streamfunction from vorticity
    if (mgridlev .ge. 2) then
        CALL vort2flux(qxvec, qyvec, omegavec, svec, mgridlev)
    else
        CALL vort2flux(qxvec, qyvec, omegavec, svec, mgridlev)
    end if
    call VecWAXPY(qhxvec, one, qxvec(1), q0xvec(1), ierr)
    call VecWAXPY(qhyvec, one, qyvec(1), q0yvec(1), ierr)
    !---------------------------------------------------------------

    !---step 2: compute increment in body position (dx) and forces
    IF (.not.num_stat) THEN
    
        !Initialize FSI iteration loop:
        fsi_count = 0
        err_fsi_glob = 10.d0
        ksp_iter = 0  !total number of ksp iterations in all fsi iterations
        
        !structural operators
        call get_M    !get mass matrix
        call get_JK  !get stiffness matrices
        call get_C    !get damping matrices
        call var_update  !evaluates Ftint=ktheta*theta
      
        if (motion_prescribed)  then !If motion is prescribed on any of the rigid bodies then advance the positions of those bodies
             !Scatter level 1
             call VecScatterBegin(ctx_force, fvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
             call VecGetArrayF90(xbvec, xb, ierr)
             call VecGetArrayF90(ybvec, yb, ierr)
             call advance_rigidbody(itime, xb, yb)   !advance the rigid body only if it has to physically move in the flow domain
             call VecScatterEnd(ctx_force, fvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)  
             !Scatter level 2
             call VecScatterBegin(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
             call VecRestoreArrayF90(xbvec, xb, ierr)
             call VecRestoreArrayF90(ybvec, yb, ierr) 
             call VecScatterEnd(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
             
             !Update position of any connected body to the rigid body
             call advance_bodies(itime, xbody_vec, y1body_vec, y2body_vec)   !fvec is only needed here to get the data structure      
              
              !Scatter
              call VecScatterBegin(ctx_body, y1body_vec, y1ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 start for y1ivec
              call VecScatterBegin(ctx_body2, y2body_vec, y2ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 start for y2ivec
              if (.not. sub_domain) call VecScatterDestroy(ctx_velx, ierr)
              if (.not. sub_domain) call VecScatterDestroy(ctx_vely, ierr)
              
              call VecScatterEnd(ctx_body, y1body_vec, y1ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 end for y1ivec
              call VecScatterBegin(ctx_force, y1ivec, xbvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 start for xbvec
              call VecScatterEnd(ctx_body2, y2body_vec, y2ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 end for y1ivec
              call VecScatterBegin(ctx_force2, y2ivec, ybvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 start for ybvec   
              if (.not. sub_domain) call ISDestroy(reg_ixs, ierr)
              if (.not. sub_domain) call ISDestroy(reg_iys, ierr)
              
              call VecScatterEnd(ctx_force, y1ivec, xbvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 end for xbvec
              call VecGetArrayReadF90(xbvec, xb, ierr)
              call VecScatterEnd(ctx_force2, y2ivec, ybvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 end for ybvec
              call VecGetArrayReadF90(ybvec, yb, ierr) 
              
              !Every time body position changes, E or PEo needs to be updated
              if (.not. sub_domain) call reg_prep(xb, yb) !Update E
              if (sub_domain) call reg_prep_sd_sparse(xb, yb)  !Update PEo
              
              call VecRestoreArrayReadF90(xbvec, xb, ierr)
              call VecRestoreArrayReadF90(ybvec, yb, ierr)
        end if
          
        !Need to get a consistent acceleration if at 1st time step:
        if (itime==1) then
            if (rank==0) print *, "computing consistent initial body acceleration..."
            call initial_acc
        end if
      
        !Store pos and vel at previous time step for iteration:
        if (mt_body>0) then
            theta0=theta
            thetad0=thetad
            thetadd0=thetadd
        end if
      
      !Begin FSI iterations
      do while ( (err_fsi_glob .ge. fsi_tol) .and. (fsi_count .le. 1000) )
      
          fsi_count = fsi_count+1
          
          call VecZeroEntries(rhsfvec, ierr)
          call compute_rhs(itime, ybody_vec)   !this rhs corresponds to torsional and deformable bodies
          !Scatter level 2
          call VecScatterBegin(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    
          if (motion_prescribed)  then  !Moving body will need its own boundary velocity terms
              call VecGetArrayReadF90(xbvec, xb, ierr)
              call VecGetArrayReadF90(ybvec, yb, ierr)
              call VecGetArrayF90(rhsfvec, rhsf, ierr)
              call compute_rhs_prescribed(itime, rhsf, xb, yb) !this rhs corresponds to adding the prescribed motion u_B
              call VecRestoreArrayF90(rhsfvec, rhsf, ierr)
              call VecRestoreArrayReadF90(xbvec, xb, ierr)
              call VecRestoreArrayReadF90(ybvec, yb, ierr)
          end if
          
          call VecScatterEnd(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
          !Scatter level 1
          call VecScatterBegin(ctx_force, yivec, rhsfvec, ADD_VALUES, SCATTER_REVERSE, ierr)
          call VecScatterEnd(ctx_force, yivec, rhsfvec, ADD_VALUES, SCATTER_REVERSE, ierr)  !scatter level 1 end   
          
          !--solve linear system for surface stress
          if (.not. sub_domain) then
          
              !regT----------------
              !First scatter those velocity flux qhxvec, qhyvec that are close to body points in local processors to the respective local processors
              call VecScatterBegin(ctx_velx, qhxvec, velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecScatterBegin(ctx_vely, qhyvec, velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecGetArrayF90(rhsfvec, rhsf, ierr)
              call VecScatterEnd(ctx_velx, qhxvec, velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecGetArrayReadF90(velxvec_interface, velx_interface, ierr)
              call VecScatterEnd(ctx_vely, qhyvec, velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecGetArrayReadF90(velyvec_interface, vely_interface, ierr)
          
              call regT(velx_interface, vely_interface, rhsf)!do regT; rhsfvec is added inside regT
              call VecRestoreArrayF90(rhsfvec, rhsf, ierr)
              call VecRestoreArrayReadF90(velxvec_interface, velx_interface, ierr)
              call VecRestoreArrayReadF90(velyvec_interface, vely_interface, ierr)
              !----------------------
          
              !GMRES to get fvec. Matrix free implementation via ab_times is used where a_times is used for E(CtAC)^-1Et and b_times for the remaining structural operators
              !Note that a_times and b_times are overlapped in a_times2 for improving communication
              call KSPSetOperators(ksp_fsi, AB, AB, ierr)
              call KSPSolve(ksp_fsi, rhsfvec, fvec, ierr)
              call KSPGetIterationNumber(ksp_fsi, iterations, ierr)
              ksp_iter = ksp_iter + iterations

          else
          
              !regT               
              call VecScatterBegin(ctx_velx_sd_sparse, qhxvec, velxvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecScatterBegin(ctx_vely_sd_sparse, qhyvec, velyvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecZeroEntries(rhsfsdvec, ierr)

              !Compute PBP^T+other structural matrices
              call compute_operator(lumat) 

              !regT over the subdomain------------------------------------
              call VecGetArrayF90(rhsfsdvec, rhsfsd, ierr)
              call VecScatterEnd(ctx_velx_sd_sparse, qhxvec, velxvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecGetArrayReadF90(velxvec_interface_sd_sparse, velx_interface, ierr)
              call VecScatterEnd(ctx_vely_sd_sparse, qhyvec, velyvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecGetArrayReadF90(velyvec_interface_sd_sparse, vely_interface, ierr)
    
              call regT_sd_sparse(velx_interface, vely_interface, rhsfsd, size(interp_indices_global), interp_indices_global)!do regT
              call VecRestoreArrayF90(rhsfsdvec, rhsfsd, ierr)
              call VecScatterBegin(ctx_intermediate, rhsfsdvec, intermediatelocalvec, INSERT_VALUES, SCATTER_FORWARD, ierr)
              call VecRestoreArrayReadF90(velxvec_interface_sd_sparse, velx_interface, ierr)
              call VecRestoreArrayReadF90(velyvec_interface_sd_sparse, vely_interface, ierr)
              call VecScatterEnd(ctx_intermediate, rhsfsdvec, intermediatelocalvec, INSERT_VALUES, SCATTER_FORWARD, ierr)       
              
              !using GMRES with operators instead of matrix-free---------------------------
              call VecGetArrayReadF90(intermediatelocalvec, rhsfsd, ierr)
              call VecGetArrayF90(rhsfvec, rhsf, ierr)
              do j=xfi,xfe
                  rhsf(j-xfi+1) = dot_product(interp_weights(:,j), rhsfsd(interp_indices(:,j))) + rhsf(j-xfi+1)
              end do
              call VecRestoreArrayF90(rhsfvec, rhsf, ierr)
              call VecRestoreArrayReadF90(intermediatelocalvec, rhsfsd, ierr)
              
              call KSPSetOperators(ksp_fsi_sd, lumat, lumat, ierr)
              call KSPSolve(ksp_fsi_sd, rhsfvec, fvec, ierr)
              !---------------------------------------
              
          end if
          
          !Redistributing to remove unphysical oscialltions
          call redistribute(fvec, f_rdstvec)
          !-------------------------------
          
          !Now lets compute increment in body position
          !Scatter level 1
          call VecScatterBegin(ctx_force, fvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecScatterEnd(ctx_force, fvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
          !Scatter level 2
          call VecScatterBegin(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call VecScatterEnd(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
          
          err_fsi=0.d0
          call compute_dx(itime, xbody_vec, y1body_vec, y2body_vec, err_fsi)!Update position and get maximum local deflection of the bodies in the local processor
          
          !Scatter
          call VecScatterBegin(ctx_body, y1body_vec, y1ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 start for y1ivec
          call VecScatterBegin(ctx_body2, y2body_vec, y2ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 start for y2ivec
          if (.not. sub_domain) call VecScatterDestroy(ctx_velx, ierr)
          if (.not. sub_domain) call VecScatterDestroy(ctx_vely, ierr)
          
          call VecScatterEnd(ctx_body, y1body_vec, y1ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 end for y1ivec
          call VecScatterBegin(ctx_force, y1ivec, xbvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 start for xbvec
          call VecScatterEnd(ctx_body2, y2body_vec, y2ivec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 2 end for y1ivec
          call VecScatterBegin(ctx_force2, y2ivec, ybvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 start for ybvec

          call MPI_Allreduce(err_fsi, err_fsi_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr) !Get maximum deflection among all bodies
          
          !update terms for next iteration if necessary:
          if (err_fsi_glob .ge. fsi_tol) then
              call get_JK  !get Jmat and Kmat matrices
              call var_update  !evaluates Ftint=ktheta*theta
          end if
          
          if (.not. sub_domain) call ISDestroy(reg_ixs, ierr)
          if (.not. sub_domain) call ISDestroy(reg_iys, ierr)
          call VecScatterEnd(ctx_force, y1ivec, xbvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 end for xbvec
          call VecGetArrayReadF90(xbvec, xb, ierr)
          call VecScatterEnd(ctx_force2, y2ivec, ybvec, INSERT_VALUES, SCATTER_REVERSE, ierr) !scatter level 1 end for ybvec
          call VecGetArrayReadF90(ybvec, yb, ierr) 
          
          !Every time body position changes, E or PEo needs to be updated
          if (.not. sub_domain) call reg_prep(xb, yb)
          if (sub_domain) call reg_prep_sd_sparse(xb, yb)
          call VecRestoreArrayReadF90(xbvec, xb, ierr)
          call VecRestoreArrayReadF90(ybvec, yb, ierr)
          
      end do
      
      fsi_iter = fsi_count
    
    ELSE

      !Doing regT------------------------
      call VecWAXPY(qhxvec, one, qxvec(1), q0xvec(1), ierr)
      call VecScatterBegin(ctx_velx, qhxvec, velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecWAXPY(qhyvec, one, qyvec(1), q0yvec(1), ierr) 
      call VecScatterBegin(ctx_vely, qhyvec, velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecGetArrayF90(rhsfvec, rhsf, ierr)
      rhsf = 0.d0

      call VecScatterEnd(ctx_velx, qhxvec, velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecGetArrayReadF90(velxvec_interface, velx_interface, ierr)
      call VecScatterEnd(ctx_vely, qhyvec, velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecGetArrayReadF90(velyvec_interface, vely_interface, ierr)
      
      call regT(velx_interface, vely_interface, rhsf)!do regT
      call VecRestoreArrayF90(rhsfvec, rhsf, ierr)
      call VecRestoreArrayReadF90(velxvec_interface, velx_interface, ierr)
      call VecRestoreArrayReadF90(velyvec_interface, vely_interface, ierr)
      !----------------------------------     
      
      !Add a velocity term using compute_rhs_prescribed subroutine
      
      !Force calculation using Cholesky decomposed matrix
      call KSPSolve(ksp, rhsfvec, fvec, ierr)
      
      call redistribute(fvec, f_rdstvec)
      
    END IF
   
    !step 3---complete omega on grid 1-------------------------------------------------------------

    call VecZeroEntries(hxvec, ierr)
    call VecZeroEntries(hyvec, ierr)
    
    !reg---------------------------------
    if (.not. sub_domain) then
        call VecGetArrayReadF90(fvec, fb, ierr)
        call reg(fb, hxvec, hyvec)
        call VecAssemblyBegin(hxvec, ierr)
        call VecAssemblyBegin(hyvec, ierr) 
        call VecRestoreArrayReadF90(fvec, fb, ierr)
    else 
        call VecZeroEntries(rhsfsdvec, ierr)
        call VecGetArrayReadF90(fvec, fb, ierr)
        rhsfsdserialarray(interp_indices_local(:)) = matmul(indicesserialmat, fb)
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
    !-------------------------------------------

    !rot---------------------------------
    call VecGetArrayF90(vortvec, vort, ierr)
    call VecAssemblyEnd(hxvec, ierr)
    call VecAssemblyEnd(hyvec, ierr)
    call DMGlobalToLocalBegin(dau, hxvec, INSERT_VALUES, hxvec_local, ierr)
    call DMGlobalToLocalBegin(dav, hyvec, INSERT_VALUES, hyvec_local, ierr)
    
    call VecGetArrayReadF90(hxvec, hx, ierr)
    call VecGetArrayReadF90(hyvec, hy, ierr)
    call rot_int(hx, hy, vort)
    call VecRestoreArrayReadF90(hxvec, hx, ierr)
    call VecRestoreArrayReadF90(hyvec, hy, ierr)
    
    call DMGlobalToLocalEnd(dau, hxvec, INSERT_VALUES, hxvec_local, ierr)
    call VecGetArrayReadF90(hxvec_local, hx, ierr)
    call DMGlobalToLocalEnd(dav, hyvec, INSERT_VALUES, hyvec_local, ierr)
    call VecGetArrayReadF90(hyvec_local, hy, ierr)
    call rot_bc(hx, hy, vort)
    call VecRestoreArrayReadF90(hxvec_local, hx, ierr)
    call VecRestoreArrayReadF90(hyvec_local, hy, ierr)
    !--------------------------------------------
    
    !Doing ainv----------------------------------
    !call KSPSolve(vort_solver(1)%ksp, rhsvortvec, vortvec, ierr)
    call ainv(xsi, xse, ysi, yse, vort)
    !--------------------------------------------
    call VecGetArrayF90(omegavec(1), omega, ierr)
    call fluid_corrector(omega, vort)
    call VecRestoreArrayF90(omegavec(1), omega, ierr)
    call VecRestoreArrayF90(vortvec, vort, ierr)
    !---------------------------------------------------------------------------------
    
    ! coarsify final omega and correct velocities on all grids---------------------
    CALL vort2flux(qxvec, qyvec, omegavec, svec, mgridlev)
    !----------------------------------------------------------------------------
  
  end subroutine advance
  
!===================================================================================

  subroutine preprocess
  !Computes E(CtAC)^-1Et matrix and performs Cholesky decomposition
  !The ability to save the computed matrix for restarts is not yet present
  
    PetscInt :: i, j
    Vec :: chol, z
    PetscScalar, Pointer :: cholv(:)!, cholm(:,:)
    PetscInt :: idxn(xfm), idyn(1)
    PetscViewer :: viewer

    if (rank==0) WRITE(*,*) 'precomputing body matrix for stationary geometry...please be patient'
    
    call VecDuplicate(fvec, chol, ierr)
    call VecDuplicate(fvec, z, ierr)
    !Compute the matrix column by column
    do i=1,nf      
      call VecZeroEntries(z, ierr)
      call VecSetValue(z, i-1, one, INSERT_VALUES, ierr)
      call VecAssemblyBegin(z, ierr)
      call VecAssemblyEnd(z, ierr)
      
      call a_times(z, chol)
      
      call VecGetArrayReadF90(chol, cholv, ierr)
      idxn=(/(j, j=xfi-1,xfe-1, 1)/)
      idyn(1)=i-1
      call matvecequal2(cholv, idxn, idyn)
      call VecRestoreArrayReadF90(chol, cholv, ierr)
        
    end do
    call VecDestroy(z, ierr)
    call VecDestroy(chol, ierr)   
    call MatAssemblyBegin(cholmat,MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(cholmat,MAT_FINAL_ASSEMBLY, ierr)
    
    call KSPSetOperators(ksp, cholmat, cholmat, ierr)
    call KSPSetUp(ksp, ierr)

  end subroutine preprocess

!===================================================================================

  subroutine preprocess_sd
  !Precomputes B matrix
  
    PetscInt :: i, j, k, ind
    Vec :: precompute, z
    PetscScalar, Pointer :: precomputev(:), precomputematarray(:,:)
    PetscInt :: idyn(1), nz_loc(nz)
    integer, allocatable :: idxn(:)
    integer :: reqs, status(MPI_STATUS_SIZE)
    PetscViewer :: savematrix, viewer
    PetscInt :: mmm, nnn
    integer :: sd_list_size  !keep it as integer instead of PetscInt
    PetscInt :: ranges(nproc+1), root
    Vec :: precompute0vec
    VecScatter :: ctx_pcvec
    IS :: precompute_is
    PetscScalar, allocatable :: dummy1(:)
    PetscInt, allocatable :: dummy2(:)
    PetscInt :: dummy3, dummy4(0)
    PetscScalar :: precomputeserialmat_reuse(nz)
    PetscInt :: precomputeserialindicesmat2_reuse(nz), mknz
    logical :: mk(nz)
    
         Vec :: velfullvec, velfullvec2
         Mat :: a_reg, a_ctac, a_regT, dummat, dummat2
    
    if (.not. sub_domain_precomputed) then
    if (rank==0) WRITE(*,*) 'precomputing sub domain matrix ...please be patient'
    
    call VecDuplicate(xsdvec, precompute, ierr)
    call VecDuplicate(xsdvec, z, ierr)
    call VecGetOwnershipRanges(xsdvec, ranges, ierr)
    
    call VecCreateSeq(PETSC_COMM_SELF, nf_sd, precompute0vec, ierr)
    
    !Compute column B_j one at a time
    do i=1,nf_sd      

        sd_list_size = precomputeserialnz(i)

        if (allocated(idxn)) deallocate(idxn)
        allocate(idxn(precomputeserialnz(i)))
        if (rank==0) then
            idxn = precomputeserialindicesmat(i)%value(1:precomputeserialnz(i))
        end if
        root = 0
        call MPI_Ibcast(idxn, precomputeserialnz(i), MPI_INT, root, MPI_COMM_WORLD, reqs, ierr)
        
        call VecZeroEntries(z, ierr)
        call VecSetValue(z, i-1, one, INSERT_VALUES, ierr)
        call VecAssemblyBegin(z, ierr)
        call VecAssemblyEnd(z, ierr)
        
        call MPI_Wait(reqs, status, ierr)
        
        !Get necessary weights in Eo corresponding to the sparsity of EE^T as described in get_sparsity_EEt
        if (i==1) then
            call setup_reg_subdomain(sd_list_size, idxn)
            call VecCreateSeq(MPI_COMM_SELF, interfacex_sd_sparse, velxvec_interface_sd_sparse, ierr)
            call VecCreateSeq(MPI_COMM_SELF, interfacey_sd_sparse, velyvec_interface_sd_sparse, ierr)
            call VecScatterCreate(velocityxvec, reg_ixs_sd_sparse, velxvec_interface_sd_sparse, interfacex_sd_is_sparse, ctx_velx_sd_sparse, ierr)
            call VecScatterCreate(velocityyvec, reg_iys_sd_sparse, velyvec_interface_sd_sparse, interfacey_sd_is_sparse, ctx_vely_sd_sparse, ierr)
        else
            call reg_prep_sd(sd_list_size, idxn)   !setting up regularization indices corresponding to the sparsity    
        end if
         
        if (i>1) then
            call VecScatterEnd(ctx_pcvec, precompute, precompute0vec, INSERT_VALUES, SCATTER_FORWARD, ierr) 
            call ISDestroy(precompute_is, ierr)
            call VecScatterDestroy(ctx_pcvec, ierr)
        end if

        !Compute the column B_j at idxn locations
        call a_times_sd(z, precompute, sd_list_size, idxn)

        !Drop tolerance filter and saving are performed only on the 0th processor
        if (rank==0 .and. i>1) then
            mknz = precomputeserialnz(i-1)
            mk = .FALSE.
            call VecGetArrayReadF90(precompute0vec, precomputev, ierr)
            where(abs(precomputev(precomputeserialindicesmat(i-1)%value)) > droptol*precomputev(i-1)) mk(1:mknz) = .TRUE.   !Apply drop tolerance filter
            precomputeserialnz(i-1) = count(mk(1:mknz))
            precomputeserialmat_reuse(:) = 0.d0
            precomputeserialindicesmat2_reuse(:) = 0
            precomputeserialindicesmat2_reuse(1:precomputeserialnz(i-1)) = pack(precomputeserialindicesmat(i-1)%value, mk(1:mknz))
            precomputeserialmat_reuse(1:precomputeserialnz(i-1)) = precomputev(precomputeserialindicesmat2_reuse(1:precomputeserialnz(i-1)))
            call VecRestoreArrayReadF90(precompute0vec, precomputev, ierr)
            
            !Write the column of sparsified B_j into a file
            if (i==2) OPEN(unit=100,file="input/sd_precomputed.chd",form="unformatted",status="replace")
            if (i>2) OPEN(unit=100,file="input/sd_precomputed.chd",form="unformatted",status="old",position="append")
            WRITE(100) precomputeserialnz(i-1)
            WRITE(100) precomputeserialmat_reuse(1:precomputeserialnz(i-1)), precomputeserialindicesmat2_reuse(1:precomputeserialnz(i-1))
            CLOSE(100)
        end if

        !Scatter the column B_j to processor 0 (Scattering communication is overlapped with some computation)
        if (rank==0) then
            call ISCreateGeneral(MPI_COMM_SELF, precomputeserialnz(i), precomputeserialindicesmat(i)%value(1:precomputeserialnz(i))-1, PETSC_COPY_VALUES, precompute_is, ierr)
        else
            call ISCreateGeneral(MPI_COMM_SELF, izero, dummy4, PETSC_COPY_VALUES, precompute_is, ierr) !other procs should receive nothing
        end if
        call VecScatterCreate(precompute, precompute_is, precompute0vec, precompute_is, ctx_pcvec, ierr) 
        call VecScatterBegin(ctx_pcvec, precompute, precompute0vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
        
        !Same operations for the last column
        if (i==nf_sd) then
            call VecScatterEnd(ctx_pcvec, precompute, precompute0vec, INSERT_VALUES, SCATTER_FORWARD, ierr) 
            call ISDestroy(precompute_is, ierr)
            call VecScatterDestroy(ctx_pcvec, ierr)
        end if
        if (rank==0 .and. i==nf_sd) then
            mknz = precomputeserialnz(i)
            mk = .FALSE.
            call VecGetArrayReadF90(precompute0vec, precomputev, ierr)
            where(abs(precomputev(precomputeserialindicesmat(i)%value)) > droptol*precomputev(i)) mk(1:mknz) = .TRUE.
            precomputeserialnz(i) = count(mk(1:mknz))
            precomputeserialmat_reuse(:) = 0.d0
            precomputeserialindicesmat2_reuse(:) = 0
            precomputeserialindicesmat2_reuse(1:precomputeserialnz(i)) = pack(precomputeserialindicesmat(i)%value, mk(1:mknz))
            precomputeserialmat_reuse(1:precomputeserialnz(i)) = precomputev(precomputeserialindicesmat2_reuse(1:precomputeserialnz(i)))
            call VecRestoreArrayReadF90(precompute0vec, precomputev, ierr)
            
            OPEN(unit=100,file="input/sd_precomputed.chd",form="unformatted",status="old",position="append")
            WRITE(100) precomputeserialnz(i)
            WRITE(100) precomputeserialmat_reuse(1:precomputeserialnz(i)), precomputeserialindicesmat2_reuse(1:precomputeserialnz(i))
            CLOSE(100)
        end if

        if (rank==0) print*, i, size(idxn)
        
    end do    

    call VecDestroy(z, ierr)
    call VecDestroy(precompute, ierr)
    call VecDestroy(precompute0vec, ierr)
    
    if (rank==0) then
        do i = 1,nf_sd
             deallocate(precomputeserialindicesmat(i)%value)
        end do
    end if
    if (rank==0) deallocate(precomputeserialmat, precomputeserialindicesmat, precomputeserialnz, precomputeserialindicesmat2)
    if (rank .ne. 0) deallocate(precomputeserialnz)  
        
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !Load the precomputed matrix in parallel on different processors
    !precomputemat contains values; precomputeindicesmat contains indices; precomputenz contains number of non-zeros per row
    allocate(precomputemat(xsdi:xsde), precomputeindicesmat(xsdi:xsde), precomputenz(xsdi:xsde))
    OPEN(unit=100,file="input/sd_precomputed.chd",form="unformatted",status="unknown")
    do i = 1,nf_sd
        if (i.ge.xsdi .and. i.le.xsde) then
            READ(100) precomputenz(i)
            allocate(precomputemat(i)%value(precomputenz(i)), precomputeindicesmat(i)%value(precomputenz(i)))
            READ(100) precomputemat(i)%value(:), precomputeindicesmat(i)%value(:)
        else
            READ(100) dummy3
            allocate(dummy1(dummy3),dummy2(dummy3))
            READ(100) dummy1, dummy2
            deallocate(dummy1,dummy2)
        end if
    end do
    CLOSE(100)
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
  end subroutine preprocess_sd  

!===================================================================================

  subroutine ab_times(AB, xvec, yvec)
  
    Mat :: AB
    Vec :: xvec, yvec

    call a_times2(xvec, y1vec, y2vec)  !b_times2 is included in a_times2
    call VecWAXPY(yvec, one, y1vec, y2vec, ierr)
    
  end subroutine ab_times
  
!===================================================================================

  subroutine matvecequal2(cholv, idxn, idyn)
  
    PetscInt :: idxn(xfm), idyn(1)
    PetscScalar :: cholv(xfm)
    
    call MatSetValues(cholmat, xfm, idxn, ione, idyn, cholv, ADD_VALUES, ierr)
    call MatAssemblyBegin(cholmat,MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(cholmat,MAT_FINAL_ASSEMBLY, ierr)
    
  end subroutine matvecequal2
  
!===================================================================================  
  
  subroutine write_iterations(it)
  
     PetscInt :: it
     character(2) :: charnproc
     
     WRITE(charnproc,"(I2.2)") nproc
     
     if (it==1) then
     OPEN(unit=(rank+1)*105,file="output/iter"//charnproc//".dat",form="formatted",status="replace")
     else
     OPEN(unit=(rank+1)*105,file="output/iter"//charnproc//".dat",form="formatted",status="old",position="append")
     end if
    
     WRITE((rank+1)*105,*) it, ksp_iter, fsi_iter
     CLOSE((rank+1)*105)
     
  end subroutine write_iterations

!=================================================================================== 

end module timestepper
