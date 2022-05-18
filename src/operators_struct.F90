module operators_struct

  USE petsc
#include "petsc/finclude/petsc.h"
  USE parameters
  USE grid
  use variables
  implicit none
  
contains

!===================================================================================

  subroutine b_times2(xxbody_vec, yybody_vec)
  !b_times that gets incorporated into a_times2
  !Perform all computations involving individual bodies in the inner most loop and successively scatter back
  
    Mat :: B
    Vec :: xxbody_vec, yybody_vec
    PetscInt :: i, i_bdy, iter_mt, iter_md, j, iter, k
    PetscInt :: idxm(1)
    PetscInt, DIMENSION(:), ALLOCATABLE :: idxn
    PetscScalar, DIMENSION(:), ALLOCATABLE :: value
    Vec :: xbody_vec2, ybody_vec2
    PetscScalar, pointer :: xbody_val(:), xbody_val2(:), ybody_val(:), ybody_val2(:)
    PetscInt :: xbi, xbe  !starting and ending indices of fbody_vec
    PetscScalar :: si, ci
    
    iter_mt=0
    iter_md=0
    iter=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecDuplicate(xxbody_vec, xbody_vec2, ierr)
        call VecCopy(xxbody_vec, xbody_vec2, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecCreateSeq(PETSC_COMM_SELF, 2*nb_ind(i_bdy), xbody_vec2, ierr)
        call VecGetArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecGetArrayReadF90(xxbody_vec, xbody_val, ierr)  
        xbody_val2(1:2*nb_ind(i_bdy)) = xbody_val(iter+1:iter+2*nb_ind(i_bdy))
        call VecRestoreArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecRestoreArrayReadF90(xxbody_vec, xbody_val, ierr)        
      end if
      call VecDuplicate(xbody_vec2, ybody_vec2, ierr)
      
      call VecGetOwnershipRange(xbody_vec2, xbi, xbe, ierr)
      xbi=xbi+1  !converting to fortran indices, need not do it for xbe (check VecGetOwnershipRange for details)
      
      !------------------------------------------------------------------------------
      if (bdy(i_bdy)%fam==2) then  !Torsional body equations
        iter_mt=iter_mt+1
        
        !Building sol_Jmat
        allocate(value(2*nb_ind(i_bdy)), idxn(2*nb_ind(i_bdy)))
        si = sin(bdy(i_bdy)%gamma0-theta(iter_mt))
        ci = cos(bdy(i_bdy)%gamma0-theta(iter_mt))
        do j=xbi,xbe
          idxm(1)=j-1  !-1 to convert to c indexing
          do k=1,2*nb_ind(i_bdy)
             idxn(k)=k-1   !-1 to convert to c indexing
        
             if (j .le. nb_ind(i_bdy)) then
             
                if (k .le. nb_ind(i_bdy)) then
                   value(k) = idxm(1)*idxn(k)*si**2*bdy(i_bdy)%deltas**2
                end if
                
                if (k >  nb_ind(i_bdy)) then
                   value(k) = -idxm(1)*(idxn(k)-nb_ind(i_bdy))*si*ci*bdy(i_bdy)%deltas**2
                end if
                
             end if
          
             if (j > nb_ind(i_bdy)) then
             
                if (k .le. nb_ind(i_bdy)) then
                   value(k) = -(idxm(1)-nb_ind(i_bdy))*idxn(k)*si*ci*bdy(i_bdy)%deltas**2
                end if
                
                if (k > nb_ind(i_bdy)) then
                   value(k) = (idxm(1)-nb_ind(i_bdy))*(idxn(k)-nb_ind(i_bdy))*ci**2*bdy(i_bdy)%deltas**2
                end if
                
             end if
          
          end do  !k-loop
          value = value/(4.D0/dt**2.D0*Mtmat(iter_mt) + 2.D0/dt*Ctmat(iter_mt) + Jmat(iter_mt))
          call MatSetValues(sol_Jmat(iter_mt), ione, idxm, 2*nb_ind(i_bdy), idxn, value, INSERT_VALUES, ierr)
        end do   !j-loop
        deallocate(value, idxn)
        
        call MatAssemblyBegin(sol_Jmat(iter_mt), MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(sol_Jmat(iter_mt), MAT_FINAL_ASSEMBLY, ierr)        
        call MatMult(sol_Jmat(iter_mt), xbody_vec2, ybody_vec2, ierr)
        call VecScale(ybody_vec2, (delta**2/bdy(i_bdy)%deltas/dt)*(2*bdy(i_bdy)%deltas/dt), ierr)
              
      end if
      
      !-------------------------------------------------------------------------------------
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body equations
        iter_md=iter_md+1
      end if
      
      !------------------------------------------------------------------------------------
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecCopy(ybody_vec2, yybody_vec, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecGetArrayReadF90(ybody_vec2, ybody_val2, ierr)
        call VecGetArrayF90(yybody_vec, ybody_val, ierr)  
        ybody_val(iter+1:iter+2*nb_ind(i_bdy)) = ybody_val2(1:2*nb_ind(i_bdy))
        call VecRestoreArrayReadF90(ybody_vec2, ybody_val2, ierr)
        call VecRestoreArrayF90(yybody_vec, ybody_val, ierr)        
      end if
      
      iter=iter+2*nb_ind(i_bdy)
      call VecDestroy(xbody_vec2, ierr)
      call VecDestroy(ybody_vec2, ierr)
      
    end do
      
  end subroutine b_times2

!===================================================================================

  subroutine get_M
  !Build mass matrix
  
    PetscInt :: i, i_bdy, iter_mt, iter_md
  
    iter_mt=0
    iter_md=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body mass matrix
        iter_mt=iter_mt+1
        Mtmat(iter_mt) = bdy(i_bdy)%itheta
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body matrix
        iter_md=iter_md+1
      end if
      
    end do
  
  end subroutine get_M
  
!===================================================================================

  subroutine get_JK
  !Build stiffness matrix
  
    PetscInt :: i, i_bdy, iter_mt, iter_md
  
    iter_mt=0
    iter_md=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body stiffness matrix
        iter_mt=iter_mt+1
        Jmat(iter_mt) = bdy(i_bdy)%ktheta
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body stiffness matrix
        iter_md=iter_md+1
        !Assemble Kmat matrix here
      end if
      
    end do
  
  end subroutine get_JK
  
!===================================================================================  
  
  subroutine get_C
  !Build damping matrices
  
    PetscInt :: i, i_bdy, iter_mt, iter_md
  
    iter_mt=0
    iter_md=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body mass matrix
        iter_mt=iter_mt+1
        Ctmat(iter_mt) = bdy(i_bdy)%ctheta
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body matrix
        iter_md=iter_md+1
      end if
      
    end do
  
  end subroutine get_c
  
!===================================================================================

  subroutine var_update
  
    PetscInt :: i, i_bdy, iter_mt, iter_md
  
    iter_mt=0
    iter_md=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body: Ftint=ktheta*theta
        iter_mt=iter_mt+1
        Ftint(iter_mt) = Jmat(iter_mt) * theta(iter_mt)
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body 
        iter_md=iter_md+1
        !evaluate Fdint
      end if
      
    end do
  
  
  end subroutine var_update  
  
!===================================================================================
  
  subroutine initial_acc
  
    !this evaluates initial acceleration: thetadd and udd
  
    PetscInt :: i, i_bdy, iter_mt, iter_md, it
    Vec :: fbody_vec2, BRvec
    PetscScalar, pointer :: fbody_val2(:), fbody_val(:)
    PetscInt :: iter
    PetscScalar :: BRf, gravity, dynamics
    
    !Scatter level 1
    call VecScatterBegin(ctx_force, fvec, fivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(ctx_force, fvec, fivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    
    !Scatter level 2
    call VecScatterBegin(ctx_body, fivec, fbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(ctx_body, fivec, fbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
  
    iter_mt=0
    iter_md=0
    iter=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecDuplicate(fbody_vec, fbody_vec2, ierr)
        call VecCopy(fbody_vec, fbody_vec2, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecCreateSeq(PETSC_COMM_SELF, 2*nb_ind(i_bdy), fbody_vec2, ierr)
        call VecGetArrayF90(fbody_vec2, fbody_val2, ierr)
        call VecGetArrayReadF90(fbody_vec, fbody_val, ierr)  
        fbody_val2(1:2*nb_ind(i_bdy)) = fbody_val(iter+1:iter+2*nb_ind(i_bdy))
        call VecRestoreArrayF90(fbody_vec2, fbody_val2, ierr)
        call VecRestoreArrayReadF90(fbody_vec, fbody_val, ierr)        
      end if     
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body thetadd
        iter_mt=iter_mt+1
        
        call compute_gravity(i_bdy, theta(iter_mt), gravity) !add gravity
        it = 0
        call compute_dynamics(it, i_bdy, theta(iter_mt), dynamics)  !add any additional forces such as psudo forces
        
        call VecDuplicate(fbody_vec2, BRvec, ierr)
        call BR(i_bdy, theta(iter_mt), BRvec)
        call VecDot(BRvec, fbody_vec2, BRf, ierr)
        BRf = BRf*delta/bdy(i_bdy)%deltas/dt  !multiply by scaling factor
        thetadd(iter_mt) = 1/Mtmat(iter_mt)*( -Ctmat(iter_mt)*thetad(iter_mt) -Ftint(iter_mt) + BRf*bdy(i_bdy)%deltas + gravity + dynamics)
        call VecDestroy(BRvec, ierr)
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body udd
        iter_md=iter_md+1
        !evaluate udd
      end if
      
      iter=iter+2*nb_ind(i_bdy)
      call VecDestroy(fbody_vec2, ierr)
    end do
    
  end subroutine initial_acc

!===================================================================================

  subroutine BR(i_bdy, thet, BRvec) !this should have been named QR instead
  
    !this performs Q_i R_i times xvec
    !this is also equivalent to (R_i^T Q_i^T times xvec)
  
    Vec :: BRvec
    PetscInt :: i_bdy, xbi, xbe  !starting and ending indices of fbody_vec2
    PetscScalar, pointer :: BRval(:)
    PetscScalar :: thet
    
    call VecGetOwnershipRange(BRvec, xbi, xbe, ierr)
    xbi=xbi+1  !converting to fortran indices, need not do it for xbe (check VecGetOwnershipRange for details)
    
    call VecGetArrayF90(BRvec, BRval, ierr)
    call BR_build(xbi, xbe, i_bdy, thet, BRval)
    call VecRestoreArrayF90(BRvec, BRval, ierr)
  
  end subroutine BR
  
!===================================================================================

  subroutine BR_build(xbi, xbe, i_bdy, thet, BRval)
  
    !this contructs Q_i R_i vector
    
    PetscInt :: xbi, xbe, i_bdy, j
    PetscScalar :: BRval(xbi:xbe), thet
    
    do j=xbi,xbe
        
          if (j .le. nb_ind(i_bdy)) then
            BRval(j) = (j-1)*bdy(i_bdy)%deltas * sin(bdy(i_bdy)%gamma0-thet)
          end if
          
          if (j > nb_ind(i_bdy)) then
            BRval(j) = -(j-nb_ind(i_bdy)-1)*bdy(i_bdy)%deltas * cos(bdy(i_bdy)%gamma0-thet)
          end if
          
    end do
  
  end subroutine BR_build
  
!===================================================================================

  subroutine compute_rhs(it, yybody_vec)
  !Computes the rhs terms corresponding to structural terms
  
    Vec :: yybody_vec
    Vec :: ybody_vec2
    PetscInt :: it, i, i_bdy, iter_mt, iter_md, iter
    PetscScalar, pointer :: ybody_val2(:), ybody_val(:)
    PetscScalar :: rtheta_k, Jrphi_k, gravity, dynamics
    
    iter_mt=0
    iter_md=0
    iter=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecDuplicate(yybody_vec, ybody_vec2, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecCreateSeq(PETSC_COMM_SELF, 2*nb_ind(i_bdy), ybody_vec2, ierr)  
      end if
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body thetadd
        iter_mt=iter_mt+1
        
        call compute_gravity(i_bdy, theta(iter_mt), gravity)
        call compute_dynamics(it, i_bdy, theta(iter_mt), dynamics)  
        
        rtheta_k = thetad0(iter_mt)+thetad(iter_mt) + 2.D0/dt*(theta0(iter_mt)-theta(iter_mt)) 
        Jrphi_k = 2.D0/dt/(4.D0/dt**2.D0*Mtmat(iter_mt) + 2.D0/dt*Ctmat(iter_mt) + Jmat(iter_mt))*( Mtmat(iter_mt)*( 4.D0/dt**2*theta0(iter_mt) + 4.D0/dt*thetad0(iter_mt) + thetadd0(iter_mt) ) + Ctmat(iter_mt)*(thetad0(iter_mt) + 2.D0/dt*theta0(iter_mt)) - (4.D0/dt**2*Mtmat(iter_mt) + 2.D0/dt*Ctmat(iter_mt) + Jmat(iter_mt))*theta(iter_mt) + gravity + dynamics)

        call BR(i_bdy, theta(iter_mt), ybody_vec2)
        call VecScale(ybody_vec2, delta*(rtheta_k - Jrphi_k - thetad(iter_mt)), ierr)  !here thetad is from r^c(k)
        
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body udd
        iter_md=iter_md+1
        !evaluate udd
      end if
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecCopy(ybody_vec2, yybody_vec, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecGetArrayReadF90(ybody_vec2, ybody_val2, ierr)
        call VecGetArrayF90(yybody_vec, ybody_val, ierr)  
        ybody_val(iter+1:iter+2*nb_ind(i_bdy)) = ybody_val2(1:2*nb_ind(i_bdy))
        call VecRestoreArrayReadF90(ybody_vec2, ybody_val2, ierr)
        call VecRestoreArrayF90(yybody_vec, ybody_val, ierr)
      end if
      
      iter = iter+2*nb_ind(i_bdy)
      call VecDestroy(ybody_vec2, ierr)
    end do

  end subroutine compute_rhs
  
!===================================================================================

  subroutine compute_dx(it, xxbody_vec, yy1body_vec, yy2body_vec, err_fsi)
  
    Vec :: xxbody_vec, yy1body_vec, yy2body_vec
    PetscInt :: it, i, i_bdy, iter_mt, iter_md, iter
    Vec :: xbody_vec2, BRvec
    Vec :: y1body_vec2, y2body_vec2
    PetscScalar, pointer :: xbody_val(:), xbody_val2(:), ybody_val(:), ybody_val2(:), y1val(:), y2val(:)
    PetscScalar :: Jrphi_k, BRf, err_fsi, gravity, dynamics
    PetscInt :: xbi, xbe
    PetscScalar :: d_xb(mt_body)
    
    iter_mt=0
    iter_md=0
    iter=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecDuplicate(xxbody_vec, xbody_vec2, ierr)
        call VecCopy(xxbody_vec, xbody_vec2, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecCreateSeq(PETSC_COMM_SELF, 2*nb_ind(i_bdy), xbody_vec2, ierr)
        call VecGetArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecGetArrayReadF90(xxbody_vec, xbody_val, ierr)  
        xbody_val2(1:2*nb_ind(i_bdy)) = xbody_val(iter+1:iter+2*nb_ind(i_bdy))
        call VecRestoreArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecRestoreArrayReadF90(xxbody_vec, xbody_val, ierr)        
      end if
      call VecDuplicate(xbody_vec2, y1body_vec2, ierr)
      call VecDuplicate(xbody_vec2, y2body_vec2, ierr)
      
      call VecGetOwnershipRange(xbody_vec2, xbi, xbe, ierr)
      xbi=xbi+1  !converting to fortran indices, need not do it for xbe (check VecGetOwnershipRange for details)
      
      !------------------------------------------------------------------------------
      if (bdy(i_bdy)%fam==2) then  !Torsional body equations
        iter_mt=iter_mt+1
        call compute_gravity(i_bdy, theta(iter_mt), gravity)
        call compute_dynamics(it, i_bdy, theta(iter_mt), dynamics)  
        
        Jrphi_k = 1/(4.D0/dt**2*Mtmat(iter_mt) + 2.D0/dt*Ctmat(iter_mt) + Jmat(iter_mt))*( Mtmat(iter_mt)*( 4/dt**2*theta0(iter_mt) + 4/dt*thetad0(iter_mt) + thetadd0(iter_mt) ) + Ctmat(iter_mt)*(thetad0(iter_mt) + 2.D0/dt*theta0(iter_mt)) - (4/dt**2*Mtmat(iter_mt) + 2.D0/dt*Ctmat(iter_mt) +Jmat(iter_mt))*theta(iter_mt) + gravity + dynamics )
        
        call VecDuplicate(xbody_vec2, BRvec, ierr)
        call BR(i_bdy, theta(iter_mt), BRvec)
        call VecDot(BRvec, xbody_vec2, BRf, ierr)
        BRf = bdy(i_bdy)%deltas/(4/dt**2*Mtmat(iter_mt)+Jmat(iter_mt))*BRf*(delta/bdy(i_bdy)%deltas/dt)
        dtheta(iter_mt) = Jrphi_k + BRf
        d_xb(iter_mt) = abs(dtheta(iter_mt))
        
        !Update deflection angles, velocity and acceleration
        theta(iter_mt) = theta(iter_mt) + dtheta(iter_mt)
        thetad(iter_mt) = -thetad0(iter_mt) + 2/dt*(theta(iter_mt)-theta0(iter_mt))
        thetadd(iter_mt) = -thetadd0(iter_mt) - 4/dt*thetad0(iter_mt) + 4/dt**2*(theta(iter_mt)-theta0(iter_mt))
        
        !get the new hinge location and angle of equilibrium of the rigid body to be used it for the torsional body
        call body_user_prescribed(it, i_bdy)  
        
        !---step 3 : update body positions
        call VecGetArrayF90(y1body_vec2, y1val, ierr)
        call VecGetArrayF90(y2body_vec2, y2val, ierr)
        call update_xb(xbi, xbe, i_bdy, theta(iter_mt), y1val, y2val)
        call VecRestoreArrayF90(y1body_vec2, y1val, ierr)
        call VecRestoreArrayF90(y2body_vec2, y2val, ierr)
              
        call VecDestroy(BRvec, ierr)
      end if
      
      !-------------------------------------------------------------------------------------
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body equations
        iter_md=iter_md+1
      end if
      
      !------------------------------------------------------------------------------------
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecCopy(y1body_vec2, yy1body_vec, ierr)
        call VecCopy(y2body_vec2, yy2body_vec, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecGetArrayReadF90(y1body_vec2, ybody_val2, ierr)
        call VecGetArrayF90(yy1body_vec, ybody_val, ierr)  
        ybody_val(iter+1:iter+2*nb_ind(i_bdy)) = ybody_val2(1:2*nb_ind(i_bdy))
        call VecRestoreArrayReadF90(y1body_vec2, ybody_val2, ierr)
        call VecRestoreArrayF90(yy1body_vec, ybody_val, ierr)
        
        call VecGetArrayReadF90(y2body_vec2, ybody_val2, ierr)
        call VecGetArrayF90(yy2body_vec, ybody_val, ierr)  
        ybody_val(iter+1:iter+2*nb_ind(i_bdy)) = ybody_val2(1:2*nb_ind(i_bdy))
        call VecRestoreArrayReadF90(y2body_vec2, ybody_val2, ierr)
        call VecRestoreArrayF90(yy2body_vec, ybody_val, ierr)      
      end if
      
      iter=iter+2*nb_ind(i_bdy)
      call VecDestroy(xbody_vec2, ierr)
      call VecDestroy(y1body_vec2, ierr)
      call VecDestroy(y2body_vec2, ierr)
    end do
    
    if (mrb_body>0) err_fsi = maxval(d_xb)
    if (mrb_body==0) err_fsi = 0.D0
  
  end subroutine compute_dx
  
!===================================================================================

  subroutine advance_bodies(it, xxbody_vec, yy1body_vec, yy2body_vec)
  !Update position of any torsional or deformable body that is connected to a moving rigid body
  !Only translated according to the hinge location. No rotation involved.
  
    Vec :: xxbody_vec, yy1body_vec, yy2body_vec
    PetscInt :: it, i, i_bdy, iter_mt, iter_md, iter
    Vec :: xbody_vec2, BRvec!, xivec, xbody_vec, 
    Vec :: y1body_vec2, y2body_vec2!, y1ivec, y1body_vec, y2ivec, y2body_vec
    PetscScalar, pointer :: xbody_val(:), xbody_val2(:), ybody_val(:), ybody_val2(:), y1val(:), y2val(:)
    PetscScalar :: Jrphi_k, BRf, err_fsi
    PetscInt :: xbi, xbe
    PetscScalar :: d_xb(mt_body), theta_prev

    iter_mt=0
    iter_md=0
    iter=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecDuplicate(xxbody_vec, xbody_vec2, ierr)
        call VecCopy(xxbody_vec, xbody_vec2, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecCreateSeq(PETSC_COMM_SELF, 2*nb_ind(i_bdy), xbody_vec2, ierr)
        call VecGetArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecGetArrayReadF90(xxbody_vec, xbody_val, ierr)  
        xbody_val2(1:2*nb_ind(i_bdy)) = xbody_val(iter+1:iter+2*nb_ind(i_bdy))
        call VecRestoreArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecRestoreArrayReadF90(xxbody_vec, xbody_val, ierr)        
      end if
      call VecDuplicate(xbody_vec2, y1body_vec2, ierr)
      call VecDuplicate(xbody_vec2, y2body_vec2, ierr)
      
      call VecGetOwnershipRange(xbody_vec2, xbi, xbe, ierr)
      xbi=xbi+1  !converting to fortran indices, need not do it for xbe (check VecGetOwnershipRange for details)
      
      !------------------------------------------------------------------------------
      if (bdy(i_bdy)%fam==2) then  !Torsional body equations
        iter_mt=iter_mt+1
        
        !get the new hinge location and angle of equilibrium of the rigid body to be used it for the torsional body
        call body_user_prescribed(it, i_bdy)  
        
        !---step 3 : update body positions
        call VecGetArrayF90(y1body_vec2, y1val, ierr)
        call VecGetArrayF90(y2body_vec2, y2val, ierr)
        call update_xb(xbi, xbe, i_bdy, theta(iter_mt), y1val, y2val)
        call VecRestoreArrayF90(y1body_vec2, y1val, ierr)
        call VecRestoreArrayF90(y2body_vec2, y2val, ierr)
              
      end if
      
      !-------------------------------------------------------------------------------------
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body equations
        iter_md=iter_md+1
      end if
      
      !------------------------------------------------------------------------------------
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecCopy(y1body_vec2, yy1body_vec, ierr)
        call VecCopy(y2body_vec2, yy2body_vec, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecGetArrayReadF90(y1body_vec2, ybody_val2, ierr)
        call VecGetArrayF90(yy1body_vec, ybody_val, ierr)  
        ybody_val(iter+1:iter+2*nb_ind(i_bdy)) = ybody_val2(1:2*nb_ind(i_bdy))
        call VecRestoreArrayReadF90(y1body_vec2, ybody_val2, ierr)
        call VecRestoreArrayF90(yy1body_vec, ybody_val, ierr)
        
        call VecGetArrayReadF90(y2body_vec2, ybody_val2, ierr)
        call VecGetArrayF90(yy2body_vec, ybody_val, ierr)  
        ybody_val(iter+1:iter+2*nb_ind(i_bdy)) = ybody_val2(1:2*nb_ind(i_bdy))
        call VecRestoreArrayReadF90(y2body_vec2, ybody_val2, ierr)
        call VecRestoreArrayF90(yy2body_vec, ybody_val, ierr)      
      end if
      
      iter=iter+2*nb_ind(i_bdy)
      call VecDestroy(xbody_vec2, ierr)
      call VecDestroy(y1body_vec2, ierr)
      call VecDestroy(y2body_vec2, ierr)
    end do
    
    if (mrb_body>0) err_fsi = maxval(d_xb)
  
  end subroutine advance_bodies
  
!===================================================================================

  subroutine update_xb(xbi, xbe, i_bdy, thet, y1val, y2val)
  !Update position of body
  
      PetscInt :: xbi, xbe, i_bdy, j
      PetscScalar :: thet, y1val(xbi:xbe), y2val(xbi:xbe)
      
      do j=xbi,xbe
        
          if (j .le. nb_ind(i_bdy)) then
                y1val(j) = bdy(i_bdy)%xbf0 + (j-1)*bdy(i_bdy)%deltas*cos(bdy(i_bdy)%gamma0-thet)
                y2val(j) = bdy(i_bdy)%ybf0 + (j-1)*bdy(i_bdy)%deltas*sin(bdy(i_bdy)%gamma0-thet)
          end if
          
          if (j > nb_ind(i_bdy)) then
                y1val(j) = bdy(i_bdy)%xbf0 + (j-nb_ind(i_bdy)-1)*bdy(i_bdy)%deltas*cos(bdy(i_bdy)%gamma0-thet)
                y2val(j) = bdy(i_bdy)%ybf0 + (j-nb_ind(i_bdy)-1)*bdy(i_bdy)%deltas*sin(bdy(i_bdy)%gamma0-thet)
          end if
          
      end do   
      
  end subroutine update_xb
  
!===================================================================================

  subroutine write_theta(it)
  
    PetscInt :: it
    PetscInt :: i, i_bdy, iter_mt
    character(2) :: charbody, charnproc
    
   
    iter_mt=0
    do i=1,mrb
      i_bdy=rank_body(i)
      
      !------------------------------------------------------------------------------
      if (bdy(i_bdy)%fam==2) then  !Torsional body equations
        iter_mt=iter_mt+1
        write(charbody,"(I2.2)") rank_body(iter_mt)
        write(charnproc,"(I2.2)") nproc

        if (it==1) then
                open(unit=(rank+1)*104, file="output/theta"//charnproc//"_"//charbody//".dat",form="formatted",status="replace")
        else
                open(unit=(rank+1)*104,file="output/theta"//charnproc//"_"//charbody//".dat",form="formatted",status="old",position="append")
        end if

        write((rank+1)*104,*) it, theta(iter_mt)
        close((rank+1)*104)
      end if
    end do   
  end subroutine write_theta
  
!===================================================================================

  subroutine write_force(it,f_vec)
  
    PetscInt :: it, k
    PetscScalar, pointer :: f(:)
    Vec :: f_vec
    PetscScalar :: fx(m_body), fy(m_body)
    character(2) :: charrank, charnproc
    
    call VecGetArrayReadF90(f_vec, f, ierr)
    call write_force2(f, fx, fy)
    call VecRestoreArrayReadF90(f_vec, f, ierr)
    
    WRITE(charrank,"(I2.2)") rank
    write(charnproc,"(I2.2)") nproc
    
    if (it==1) then
    OPEN(unit=(rank+1)*103,file="output/force"//charnproc//"_"//charrank//".dat",form="formatted",status="replace")
    else
    OPEN(unit=(rank+1)*103,file="output/force"//charnproc//"_"//charrank//".dat",form="formatted",status="old",position="append")
    end if
    
    WRITE((rank+1)*103,*) it, (fx(k), k=1,m_body), (fy(k), k=1,m_body)
    CLOSE((rank+1)*103)
    
  end subroutine write_force
  
!------------------------------------------------------------------------------  

  subroutine write_force_rdst(it,f_rdst_vec)
  
    PetscInt :: it, k
    PetscScalar, pointer :: f(:)
    Vec :: f_rdst_vec
    PetscScalar :: fx(m_body), fy(m_body)
    character(2) :: charrank, charnproc
    
    call VecGetArrayReadF90(f_rdst_vec, f, ierr)
    call write_force2(f, fx, fy)
    call VecRestoreArrayReadF90(f_rdst_vec, f, ierr)
    
    WRITE(charrank,"(I2.2)") rank
    write(charnproc,"(I2.2)") nproc
    
    if (it==1) then
    OPEN(unit=(rank+1)*103,file="output/forcerdst"//charnproc//"_"//charrank//".dat",form="formatted",status="replace")
    else
    OPEN(unit=(rank+1)*103,file="output/forcerdst"//charnproc//"_"//charrank//".dat",form="formatted",status="old",position="append")
    end if
    
    WRITE((rank+1)*103,*) it, (fx(k), k=1,m_body), (fy(k), k=1,m_body)
    CLOSE((rank+1)*103)
    
  end subroutine write_force_rdst
  
!------------------------------------------------------------------------------
  
  subroutine write_force2(f, fx, fy)
    
    PetscInt :: k
    PetscScalar :: f(xfi:xfe)
    PetscScalar :: fx(m_body), fy(m_body)
    
    fx=0.d0; fy=0.d0
    DO k=xfi,xfe
    
      if (codeb(k,3)==1) then
        fx(codeb(k,1))=fx(codeb(k,1))+f(k)/dt
      else if (codeb(k,3)==2) then
        fy(codeb(k,1))=fy(codeb(k,1))+f(k)/dt
      end if
      
    END DO
    
    fx=2*delta*fx
    fy=2*delta*fy
    
  end subroutine write_force2

!===================================================================================

  subroutine redistribute(f_vec, frdst_vec)
    
    use fsinterface
    
    Vec :: f_vec, frdst_vec
    PetscScalar :: one_vec(xfi:xfe)
    PetscScalar, pointer :: f(:), wghtx(:), wghty(:), frcx(:), frcy(:), finterx(:), fintery(:)
    PetscScalar, pointer :: finterx_interface(:), fintery_interface(:), frdst(:)
    
    !Define a vector of ones to get the weights:
    one_vec=1.d0
    
    !Get appropriate weights for redistribution
    call VecZeroEntries(wghtxvec, ierr)
    call VecZeroEntries(wghtyvec, ierr)
    call reg(one_vec, wghtxvec, wghtyvec)
    call VecAssemblyBegin(wghtxvec, ierr)
    call VecAssemblyBegin(wghtyvec, ierr) 
    
    
    
    !now redistribute the force: f_rdst = E*M^-1*E^T*f, where M^-1 is
    !a diagonal matrix containing the inverse of the nonzero weights wght
    call VecZeroEntries(frcxvec, ierr)
    call VecZeroEntries(frcyvec, ierr)
    call VecGetArrayReadF90(f_vec, f, ierr)
    call reg(f, frcxvec, frcyvec)
    call VecAssemblyBegin(frcxvec, ierr)
    call VecAssemblyBegin(frcyvec, ierr) 
    call VecRestoreArrayReadF90(f_vec, f, ierr)
    
    !initialize an intermediate force vector called f_inter:
    call VecSet(finterx_vec, zero, ierr)
    call VecSet(fintery_vec, zero, ierr)
    
    call VecGetArrayReadF90(finterx_vec, finterx, ierr)
    call VecGetArrayReadF90(fintery_vec, fintery, ierr)
    
    call VecAssemblyEnd(wghtxvec, ierr)
    call VecGetArrayReadF90(wghtxvec, wghtx, ierr)
    call VecAssemblyEnd(wghtyvec, ierr)
    call VecGetArrayReadF90(wghtyvec, wghty, ierr)
    call VecAssemblyEnd(frcxvec, ierr)
    call VecGetArrayReadF90(frcxvec, frcx, ierr)
    call VecAssemblyEnd(frcyvec, ierr)
    call VecGetArrayReadF90(frcyvec, frcy, ierr)
    
    call redistribute2(wghtx, wghty, frcx, frcy, finterx, fintery)
    
    call VecRestoreArrayReadF90(wghtxvec, wghtx, ierr)
    call VecRestoreArrayReadF90(wghtyvec, wghty, ierr)
    call VecRestoreArrayReadF90(frcxvec, frcx, ierr)
    call VecRestoreArrayReadF90(frcyvec, frcy, ierr)
    call VecRestoreArrayReadF90(finterx_vec, finterx, ierr)
    call VecRestoreArrayReadF90(fintery_vec, fintery, ierr)
    
    !Doing regT----------------------------------
    call VecScatterBegin(ctx_velx, finterx_vec, velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(ctx_vely, fintery_vec, velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayF90(frdst_vec, frdst, ierr)
    frdst=0.d0
    call VecScatterEnd(ctx_velx, finterx_vec, velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velxvec_interface, finterx_interface, ierr)
    call VecScatterEnd(ctx_vely, fintery_vec, velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velyvec_interface, fintery_interface, ierr)
    
    call regT(finterx_interface, fintery_interface, frdst)!do regT
    call VecRestoreArrayF90(frdst_vec, frdst, ierr)
    call VecRestoreArrayReadF90(velxvec_interface, finterx_interface, ierr)
    call VecRestoreArrayReadF90(velyvec_interface, fintery_interface, ierr)
  
  end subroutine redistribute

!----------------------------------------------------------------------------------

  subroutine redistribute2(wghtx, wghty, frcx, frcy, finterx, fintery)
  
    PetscScalar :: wghtx(xui:xue,yui:yue), frcx(xui:xue,yui:yue), finterx(xui:xue,yui:yue)
    PetscScalar :: wghty(xvi:xve,yvi:yve), frcy(xvi:xve,yvi:yve), fintery(xvi:xve,yvi:yve)
    PetscInt :: i, j
    
    do j=yui,yue
      do i=xui,xue
        if (wghtx(i,j)>1.d-10) then !Only take the reciprocal of weights if they
                                  !are above a certain tolerance
          finterx(i,j)=1.d0/wghtx(i,j)*frcx(i,j)
        end if
      end do
    end do
    
    do j=yvi,yve
      do i=xvi,xve
        if (wghty(i,j)>1.d-10) then !Only take the reciprocal of weights if they
                                  !are above a certain tolerance
          fintery(i,j)=1.d0/wghty(i,j)*frcy(i,j)
        end if
      end do
    end do
  
  end subroutine redistribute2
  
!===================================================================================

  subroutine compute_operator(lumat)
  
      !Computing P*B*P^T + the other structural matrices such as S^TJ^-1S and I^TJ^-1S
      !See subroutine setup_reg_subdomain_indices for some details about the procedure of PBP^T
      !The ouput is put in lumat
  
      Mat :: lumat, lumatserial, intermediatemat, intermediatetranposemat
      PetscScalar :: dist(nb_sd), distance(interp_points)
      PetscInt :: nearest_neighbor(interp_points), j, jj, i, k, interp_indices_global_duplicates(nf*interp_points)
      Logical, dimension(interp_points) :: mk
      PetscViewer :: viewer
      PetscScalar, Pointer :: lumatarray(:,:), lumatserialarray(:,:), intermediatearray(:), intermediatelocalarray(:)
      PetscInt :: localmm, localnn
      PetscScalar :: precomputeserialarray(nf_sd), temp(interp_points), intermediateserialarray(xsdi:xsde), intermediateserialarray2(nf_sd)
      PetscInt :: precomputeserialindicesdummy(nf_sd)
      PetscScalar, allocatable :: precomputeserialdummymat(:,:)
      PetscInt :: iter
      PetscInt, allocatable :: indices_intermediateserial(:)

      call MatZeroEntries(lumat, ierr)
      call body_mat(lumat) !the other structural matrices S^TJ^-1S and I^TJ^-1S are computed here 
      call MatAssemblyBegin(lumat, MAT_FINAL_ASSEMBLY, ierr)
      allocate(indices_intermediateserial(interp_indices_global_nproc))
      allocate(intermediateserialmat(nf,interp_indices_global_nproc))
       
      !First performing B*P^T
      iter=0
      do jj=1,size(interp_indices_global)
              
              j = interp_indices_global(jj)!row index for BP^T   
              if (j.ge.xsdi .and. j.le.xsde) then  !this if loop can be eliminated
                  iter=iter+1
                  indices_intermediateserial(iter) = j
                  precomputeserialarray(interp_indices_global(:)) = 0.d0
                  precomputeserialarray(precomputeindicesmat(j)%value(1:precomputenz(j))) = precomputemat(j)%value(1:precomputenz(j))
                  do i=1,nf
                      intermediateserialmat(i,iter) = dot_product(precomputeserialarray(interp_indices_global_sparsemat(:,i)), interp_weights_global_sparsemat(:,i))
                  end do
              end if
      end do

      !Initialize intermediate scattering required before performing P*(BP^T)
      call VecGetArrayF90(intermediatevec, intermediatearray, ierr)
      intermediatearray(indices_intermediateserial-xsdi+1) = intermediateserialmat(1,:)!this line can be optimized
      call VecRestoreArrayF90(intermediatevec, intermediatearray, ierr)
      call VecScatterBegin(ctx_intermediate, intermediatevec, intermediatelocalvec, INSERT_VALUES, SCATTER_FORWARD, ierr)          
      call MatAssemblyEnd(lumat, MAT_FINAL_ASSEMBLY, ierr) 
       
      call MatDenseGetArrayF90(lumat, lumatarray, ierr)
      do i=1,nf
          
          call VecScatterEnd(ctx_intermediate, intermediatevec, intermediatelocalvec, INSERT_VALUES, SCATTER_FORWARD, ierr)
          
          !Computing P*(BP^T) column by column with the ith column of BP^T for the previous i loop i.e. i-1
          call VecGetArrayReadF90(intermediatelocalvec, intermediatelocalarray, ierr)
          intermediateserialarray2(interp_indices_local) = intermediatelocalarray(interp_indices_local)!this is done efficiently
          call VecRestoreArrayReadF90(intermediatelocalvec, intermediatelocalarray, ierr)
          
          if (i<nf) then
          !Update intermediatevec for next scatter
          call VecGetArrayF90(intermediatevec, intermediatearray, ierr)
          intermediatearray(indices_intermediateserial-xsdi+1) = intermediateserialmat(i+1,:) !this line can be optimized
          call VecRestoreArrayF90(intermediatevec, intermediatearray, ierr)
          !Begin Scattering intermediatevec to do mat-vec mult serially once the contents of intermediatelocalvec are loaded
          call VecScatterBegin(ctx_intermediate, intermediatevec, intermediatelocalvec, INSERT_VALUES, SCATTER_FORWARD, ierr)
          end if
          
          do k=xfi,xfe
              lumatarray(k-xfi+1,i) = lumatarray(k-xfi+1,i) + dot_product(interp_weights(:,k), intermediateserialarray2(interp_indices(:,k)))
          end do
          
      end do
      call MatDenseRestoreArrayF90(lumat, lumatarray, ierr)
      deallocate(indices_intermediateserial,intermediateserialmat)
  
  end subroutine compute_operator

!===================================================================================

  subroutine body_mat(lumat)
  
    Mat :: lumat
    Vec :: xxbody_vec, yybody_vec
    PetscInt :: i, i_bdy, iter_mt, iter_md, j, iter, k
    PetscInt :: idxm(1), idxnn, idxmm
    PetscInt, DIMENSION(:), ALLOCATABLE :: idxn
    PetscScalar, DIMENSION(:), ALLOCATABLE :: value
    Vec :: xbody_vec2, ybody_vec2
    PetscScalar, pointer :: xbody_val(:), xbody_val2(:), ybody_val(:), ybody_val2(:)
    PetscInt :: xbi, xbe  !starting and ending indices of fbody_vec
    PetscScalar :: si, ci
    
    iter_mt=0
    iter_md=0
    iter=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecDuplicate(xbody_vec, xbody_vec2, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecCreateSeq(PETSC_COMM_SELF, 2*nb_ind(i_bdy), xbody_vec2, ierr)        
      end if
      call VecGetOwnershipRange(xbody_vec2, xbi, xbe, ierr)
      xbi=xbi+1  !converting to fortran indices, need not do it for xbe (check VecGetOwnershipRange for details)
      
      !------------------------------------------------------------------------------
      if (bdy(i_bdy)%fam==2) then  !Torsional body equations
        iter_mt=iter_mt+1
        
        !Building sol_Jmat
        allocate(value(2*nb_ind(i_bdy)), idxn(2*nb_ind(i_bdy)))
        si = sin(bdy(i_bdy)%gamma0-theta(iter_mt))
        ci = cos(bdy(i_bdy)%gamma0-theta(iter_mt))
        do j=xbi,xbe
          idxm(1) = nbb_ind(i_bdy) - 2*nb_ind(i_bdy) + j-1  !-1 to convert to c indexing
          idxmm = j-1
          do k=1,2*nb_ind(i_bdy)
             idxn(k) = nbb_ind(i_bdy) - 2*nb_ind(i_bdy) + k-1   !-1 to convert to c indexing
             idxnn = k-1
        
             if (j .le. nb_ind(i_bdy)) then
             
                if (k .le. nb_ind(i_bdy)) then
                   value(k) = idxmm*idxnn*si**2*bdy(i_bdy)%deltas**2
                end if
                
                if (k >  nb_ind(i_bdy)) then
                   value(k) = -idxmm*(idxnn-nb_ind(i_bdy))*si*ci*bdy(i_bdy)%deltas**2
                end if
                
             end if
          
             if (j > nb_ind(i_bdy)) then
             
                if (k .le. nb_ind(i_bdy)) then
                   value(k) = -(idxmm-nb_ind(i_bdy))*idxnn*si*ci*bdy(i_bdy)%deltas**2
                end if
                
                if (k > nb_ind(i_bdy)) then
                   value(k) = (idxmm-nb_ind(i_bdy))*(idxnn-nb_ind(i_bdy))*ci**2*bdy(i_bdy)%deltas**2
                end if
                
             end if
          
          end do  !k-loop
          value = value/(4.D0/dt**2.D0*Mtmat(iter_mt) + 2.D0/dt*Ctmat(iter_mt) + Jmat(iter_mt))*(delta**2/bdy(i_bdy)%deltas/dt)*(2*bdy(i_bdy)%deltas/dt)
          call MatSetValues(lumat, ione, idxm, 2*nb_ind(i_bdy), idxn, value, INSERT_VALUES, ierr)
        end do   !j-loop
        deallocate(value, idxn)
              
      end if
      
      !-------------------------------------------------------------------------------------
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body equations
        iter_md=iter_md+1
      end if
      
      !------------------------------------------------------------------------------------
      
      call VecDestroy(xbody_vec2, ierr)
      
    end do
   
  end subroutine body_mat

!===================================================================================

!Additional subroutines

!===================================================================================

  subroutine b_times(xvec, yvec)
  !The original b_times
  
    Mat :: B
    Vec :: xvec, yvec!, x_rdstvec
    PetscInt :: i, i_bdy, iter_mt, iter_md, j, iter, k
    PetscInt :: idxm(1)
    PetscInt, DIMENSION(:), ALLOCATABLE :: idxn
    PetscScalar, DIMENSION(:), ALLOCATABLE :: value
    Vec :: xbody_vec2, ybody_vec2!, xivec, xbody_vec, yivec, ybody_vec
    PetscScalar, pointer :: xbody_val(:), xbody_val2(:), ybody_val(:), ybody_val2(:)
    PetscInt :: xbi, xbe  !starting and ending indices of fbody_vec
    PetscScalar :: si, ci
    
    call VecZeroEntries(yvec, ierr)
    
    !Scatter level 1
    call VecScatterBegin(ctx_force, xvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(ctx_force, xvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    
    !Scatter level 2
    call VecScatterBegin(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    
    iter_mt=0
    iter_md=0
    iter=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecDuplicate(xbody_vec, xbody_vec2, ierr)
        call VecCopy(xbody_vec, xbody_vec2, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecCreateSeq(PETSC_COMM_SELF, 2*nb_ind(i_bdy), xbody_vec2, ierr)
        call VecGetArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecGetArrayReadF90(xbody_vec, xbody_val, ierr)  
        xbody_val2(1:2*nb_ind(i_bdy)) = xbody_val(iter+1:iter+2*nb_ind(i_bdy))
        call VecRestoreArrayF90(xbody_vec2, xbody_val2, ierr)
        call VecRestoreArrayReadF90(xbody_vec, xbody_val, ierr)        
      end if
      call VecDuplicate(xbody_vec2, ybody_vec2, ierr)
      
      call VecGetOwnershipRange(xbody_vec2, xbi, xbe, ierr)
      xbi=xbi+1  !converting to fortran indices, need not do it for xbe (check VecGetOwnershipRange for details)
      
      !------------------------------------------------------------------------------
      if (bdy(i_bdy)%fam==2) then  !Torsional body equations
        iter_mt=iter_mt+1
        
        !Building sol_Jmat
        allocate(value(2*nb_ind(i_bdy)), idxn(2*nb_ind(i_bdy)))
        si = sin(bdy(i_bdy)%gamma0-theta(iter_mt))
        ci = cos(bdy(i_bdy)%gamma0-theta(iter_mt))
        do j=xbi,xbe
          idxm(1)=j-1  !-1 to convert to c indexing
          do k=1,2*nb_ind(i_bdy)
             idxn(k)=k-1   !-1 to convert to c indexing
        
             if (j .le. nb_ind(i_bdy)) then
             
                if (k .le. nb_ind(i_bdy)) then
                   value(k) = idxm(1)*idxn(k)*si**2*bdy(i_bdy)%deltas**2
                end if
                
                if (k >  nb_ind(i_bdy)) then
                   value(k) = -idxm(1)*(idxn(k)-nb_ind(i_bdy))*si*ci*bdy(i_bdy)%deltas**2
                end if
                
             end if
          
             if (j > nb_ind(i_bdy)) then
             
                if (k .le. nb_ind(i_bdy)) then
                   value(k) = -(idxm(1)-nb_ind(i_bdy))*idxn(k)*si*ci*bdy(i_bdy)%deltas**2
                end if
                
                if (k > nb_ind(i_bdy)) then
                   value(k) = (idxm(1)-nb_ind(i_bdy))*(idxn(k)-nb_ind(i_bdy))*ci**2*bdy(i_bdy)%deltas**2
                end if
                
             end if
          
          end do  !k-loop
          value = value/(4.D0/dt**2.D0*Mtmat(iter_mt) + 2.D0/dt*Ctmat(iter_mt) + Jmat(iter_mt))
          call MatSetValues(sol_Jmat(iter_mt), ione, idxm, 2*nb_ind(i_bdy), idxn, value, INSERT_VALUES, ierr)
        end do   !j-loop
        deallocate(value, idxn)
        
        call MatAssemblyBegin(sol_Jmat(iter_mt), MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(sol_Jmat(iter_mt), MAT_FINAL_ASSEMBLY, ierr)
        call MatMult(sol_Jmat(iter_mt), xbody_vec2, ybody_vec2, ierr)
        call VecScale(ybody_vec2, (delta**2/bdy(i_bdy)%deltas/dt)*(2*bdy(i_bdy)%deltas/dt), ierr)
              
      end if
      
      !-------------------------------------------------------------------------------------
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body equations
        iter_md=iter_md+1
      end if
      
      !------------------------------------------------------------------------------------
      
      !Virtual scatter level 3
      if (mrb_body==1) then
        call VecCopy(ybody_vec2, ybody_vec, ierr)
      else  !i.e. if there are many bodies in the proc then we need to subdivide it even more. Note that whatever computations to be done for the case when mrb_body>1 is done only on 1 vector
        call VecGetArrayReadF90(ybody_vec2, ybody_val2, ierr)
        call VecGetArrayF90(ybody_vec, ybody_val, ierr)  
        ybody_val(iter+1:iter+2*nb_ind(i_bdy)) = ybody_val2(1:2*nb_ind(i_bdy))
        call VecRestoreArrayReadF90(ybody_vec2, ybody_val2, ierr)
        call VecRestoreArrayF90(ybody_vec, ybody_val, ierr)        
      end if
      
      iter=iter+2*nb_ind(i_bdy)
      call VecDestroy(xbody_vec2, ierr)
      call VecDestroy(ybody_vec2, ierr)
      
    end do
    
    !Scatter level 2
    call VecScatterBegin(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    call VecScatterEnd(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    
    !Scatter level 1
    call VecScatterBegin(ctx_force, yivec, yvec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    call VecScatterEnd(ctx_force, yivec, yvec, INSERT_VALUES, SCATTER_REVERSE, ierr)
   
  end subroutine b_times

!===================================================================================

end module operators_struct
