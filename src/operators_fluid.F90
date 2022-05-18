module operators_fluid

  USE petsc
#include "petsc/finclude/petsc.h"
  USE parameters
  USE grid
  use variables
  implicit none
  
contains

!===================================================================================

  subroutine a_times2(xvec, yy1vec, yy2vec)
  
  !it is atimes and btimes combined
  
    use myfft
    use fsinterface
    use operators_struct
    
    !Mat :: A
    Vec :: xvec, yy1vec, yy2vec
    PetscScalar, pointer :: x(:), hx(:), hy(:), vort(:)
    PetscScalar, pointer :: sbc_local(:), s(:), velx(:), vely(:), s_local(:)
    PetscInt :: k
    PetscScalar, pointer :: velx_interface(:), vely_interface(:), y(:)
    
    !Scatter level 1 begin
    call VecScatterBegin(ctx_force, xvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
        
    call VecZeroEntries(yy1vec, ierr)
    call VecZeroEntries(yy2vec, ierr)

    do k=1,mgridlev
      call VecZeroEntries(vortvecs(k), ierr)
    end do
    call VecZeroEntries(hxvec, ierr)
    call VecZeroEntries(hyvec, ierr)
    
    !Scatter level 1 end
    call VecScatterEnd(ctx_force, xvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !Scatter level 2 begin
    call VecScatterBegin(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    
    !reg---------------------------------
    call VecGetArrayReadF90(xvec, x, ierr)
    call reg(x, hxvec, hyvec)
    call VecAssemblyBegin(hxvec, ierr)
    call VecAssemblyBegin(hyvec, ierr)
    !Scatter level 2 end
    call VecScatterEnd(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)  
    call VecRestoreArrayReadF90(xvec, x, ierr)
    !-------------------------------------------
    
    !b_times------------------------------
    call b_times2(xbody_vec, ybody_vec)
    !Scatter level 2 begin
    call VecScatterBegin(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !-------------------------------------------

    !rot--------------------------------- 
    call VecGetArrayF90(vortvecs(1), vort, ierr)   
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
    
    !ainv----------------------------------
    call ainv(xsi, xse, ysi, yse, vort)
    call VecRestoreArrayF90(vortvecs(1), vort, ierr)
    !--------------------------------------------
    !Scatter level 2 end
    call VecScatterEnd(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !Scatter level 1 begin
    call VecScatterBegin(ctx_force, yivec, yy2vec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !--------------------------------------------
    
    !vort2flux-----------------------------
    if (mgridlev .ge. 2) then
      CALL vort2flux(velxvecs, velyvecs, vortvecs, snvecs, mgridlev)!what is the difference between if/else here?
    else
      call vort2flux(velxvecs, velyvecs, vortvecs, snvecs, mgridlev)
    end if
    !end vort2flux--------------------------------------------
    
    !regT----------------------------------    
    call VecScatterBegin(ctx_velx, velxvecs(1), velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(ctx_vely, velyvecs(1), velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayF90(yy1vec, y, ierr)
    call VecScatterEnd(ctx_velx, velxvecs(1), velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velxvec_interface, velx_interface, ierr)
    call VecScatterEnd(ctx_vely, velyvecs(1), velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velyvec_interface, vely_interface, ierr)
    
    call regT(velx_interface, vely_interface, y)!do regT
    call VecRestoreArrayF90(yy1vec, y, ierr)
    call VecRestoreArrayReadF90(velxvec_interface, velx_interface, ierr)
    call VecRestoreArrayReadF90(velyvec_interface, vely_interface, ierr)
    !--------------------------------------------
    !Scatter level 1 end
    call VecScatterEnd(ctx_force, yivec, yy2vec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !--------------------------------------------
    
  end subroutine a_times2

!===================================================================================

  subroutine a_times_sd(xvec, yvec, sd_list_size, sd_list)
  
    use myfft
    use fsinterface
    
    Vec :: xvec, yvec
    integer, intent(in) :: sd_list_size
    PetscInt, intent(in) :: sd_list(sd_list_size)
    PetscScalar, pointer :: x(:), hx(:), hy(:), vort(:)
    PetscScalar, pointer :: sbc_local(:), s(:), velx(:), vely(:), s_local(:)
    PetscInt :: k
    PetscScalar, pointer :: velx_interface(:), vely_interface(:), y(:)
        
    do k=1,mgridlev
      call VecZeroEntries(vortvecs(k), ierr)
    end do
    call VecZeroEntries(hxvec, ierr)
    call VecZeroEntries(hyvec, ierr)
    call VecZeroEntries(yvec, ierr)
    
    !Doing reg---------------------------------
    call VecGetArrayReadF90(xvec, x, ierr)
    call reg_sd_sparse(x, hxvec, hyvec, sd_list_size, sd_list) 
    call VecRestoreArrayReadF90(xvec, x, ierr)
    !-------------------------------------------

    !Doing rot--------------------------------- 
    call VecGetArrayF90(vortvecs(1), vort, ierr)   
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
    call ainv(xsi, xse, ysi, yse, vort)
    call VecRestoreArrayF90(vortvecs(1), vort, ierr)
    !--------------------------------------------
    
    !Doing vort2flux-----------------------------
    !either use vort2flux subroutine
    if (mgridlev .ge. 2) then
      CALL vort2flux(velxvecs, velyvecs, vortvecs, snvecs, mgridlev)!what is the difference between if/else here?
    else
      call vort2flux(velxvecs, velyvecs, vortvecs, snvecs, mgridlev)
    end if
    !end vort2flux--------------------------------------------
    
    !Doing regT----------------------------------
    call VecScatterBegin(ctx_velx_sd_sparse, velxvecs(1), velxvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(ctx_vely_sd_sparse, velyvecs(1), velyvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayF90(yvec, y, ierr)
    call VecScatterEnd(ctx_velx_sd_sparse, velxvecs(1), velxvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velxvec_interface_sd_sparse, velx_interface, ierr)
    call VecScatterEnd(ctx_vely_sd_sparse, velyvecs(1), velyvec_interface_sd_sparse, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velyvec_interface_sd_sparse, vely_interface, ierr)
    
    call regT_sd_sparse(velx_interface, vely_interface, y, sd_list_size, sd_list)!do regT
    call VecRestoreArrayF90(yvec, y, ierr)
    call VecRestoreArrayReadF90(velxvec_interface_sd_sparse, velx_interface, ierr)
    call VecRestoreArrayReadF90(velyvec_interface_sd_sparse, vely_interface, ierr)
    !--------------------------------------------
    
  end subroutine a_times_sd
    
!===================================================================================

  subroutine rot_int(hx, hy, vort)

   !***************************************************************!
   !*   Transpose of curl                                         *!
   !***************************************************************!
  
    PetscScalar :: hx(xui:xue,yui:yue), hy(xvi:xve,yvi:yve), vort(xsi:xse,ysi:yse)
    PetscInt :: i, j
    
    forall(i=xsi+1:xse-1,j=ysi+1:yse-1); vort(i,j) = hy(i,j) - hy(i-1,j) - hx(i,j) + hx(i,j-1); end forall
  
  end subroutine rot_int
  
!===================================================================================

  subroutine rot_bc(hx, hy, vort)
  
    PetscScalar :: hx(xugi:xuge,yugi:yuge), hy(xvgi:xvge,yvgi:yvge), vort(xsi:xse,ysi:yse)
    PetscInt :: i, j
    
    do i=xsi,xse
      j=ysi; vort(i,j) = hy(i,j) - hy(i-1,j) - hx(i,j) + hx(i,j-1)
      j=yse; vort(i,j) = hy(i,j) - hy(i-1,j) - hx(i,j) + hx(i,j-1)
    end do
    
    do j=ysi,yse
      i=xsi; vort(i,j) = hy(i,j) - hy(i-1,j) - hx(i,j) + hx(i,j-1)
      i=xse; vort(i,j) = hy(i,j) - hy(i-1,j) - hx(i,j) + hx(i,j-1)
    end do
  
  end subroutine rot_bc
      
!===================================================================================

  subroutine vort2flux(velxvec, velyvec, vortvec, snvec, mgridlev)
  
    !***************************************************************!
    !*    Multiscale method to solve C^T C s = omega               *!
    !*    and return the velocity, C s.                            *!
    !*    Results are returned in vel on each of the first nlev    *!
    !*    grids.                                                   *!
    !*    Warning: the vorticity field on all but the finest mesh  *!
    !*     is modified by the routine in the following way:        *!
    !*     the value in the center of the domain is interpolated   *!
    !*     from the next finer mesh (the value near the edge is    *!
    !*     not changed.                                            *!
    !***************************************************************!
    
    use myfft
    
    PetscInt :: mgridlev, k
    Vec, Pointer :: velxvec(:), velyvec(:), vortvec(:), snvec(:)
    PetscScalar, Pointer :: vort1(:), vort(:), s(:), s_local(:), velx(:), vely(:)
    PetscScalar, Pointer :: sbc_local(:)
    PetscScalar :: vort2(xsm*ysm)
    
    type pointertoarray2
     PetscScalar, allocatable :: array2(:)
    end type pointertoarray2
    
    type(pointertoarray2), dimension(mgridlev) :: vort2pointer
    
    !Performing coarsify---------------------------------
    !i.e. we update vorticity at coarser mesh via interpolation of vorticity at a finer mesh
    if (mgridlev .ge. 2) then
      DO k=2,mgridlev

        call DMGlobalToLocalBegin(das, vortvec(k-1), INSERT_VALUES, vortvec_local1, ierr)
        call VecGetArrayF90(velxvec(k), velxpointer(k)%array, ierr)
        call VecGetArrayF90(velyvec(k), velypointer(k)%array, ierr)        
        
        call DMGlobalToLocalEnd(das, vortvec(k-1), INSERT_VALUES, vortvec_local1, ierr)
        call VecGetArrayReadF90(vortvec_local1, vort1, ierr)
        call coarsify(vortvec(k), vort1)
        call VecAssemblyBegin(vortvec(k), ierr)
        
        call VecGetArrayF90(snvec(k), spointer(k)%array, ierr)
        spointer(k)%array = 0.d0
        call VecRestoreArrayReadF90(vortvec_local1, vort1, ierr)
        allocate(vort2pointer(k-1)%array2(xsm*ysm))
        call VecGetArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
        vort2pointer(k-1)%array2 = vortpointer(k-1)%array
        call VecRestoreArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
        call VecAssemblyEnd(vortvec(k), ierr)  
      END DO
      
      k = 1
      call VecGetArrayF90(velxvec(k), velxpointer(k)%array, ierr)
      call VecGetArrayF90(velyvec(k), velypointer(k)%array, ierr)
      call VecGetArrayF90(snvec(k), spointer(k)%array, ierr)
      spointer(k)%array = 0.d0
      
      k = mgridlev+1
      allocate(vort2pointer(k-1)%array2(xsm*ysm))
      call VecGetArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
      vort2pointer(k-1)%array2 = vortpointer(k-1)%array
      call VecRestoreArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
      
    end if
    !---------------------------------------------------
    
    !Performing ctci to obtain streamfunction from vorticity at the coarsest grid-----------------------------------------
    call ctci(xsi, xse, ysi, yse, vort2pointer(mgridlev)%array2, spointer(mgridlev)%array)
    call VecRestoreArrayF90(snvec(mgridlev), spointer(mgridlev)%array, ierr)
    !-----------------------------------------------------
    
    !Initiate boundary vectors----------------------------
    call DMGlobalToLocalBegin(das, snvec(mgridlev), INSERT_VALUES, svec_local, ierr)
    call VecGetArrayReadF90(streambcvec_local, sbc_local, ierr)
    sbc_local = 0.d0
    !-----------------------------------------------------
    
    !Performing curl to obtain velocity from streamfunction------------------------------------------    
    !First performing curl on the interior points
    call VecGetArrayReadF90(snvec(mgridlev), s, ierr)
    call curl_int(s, sbc_local, velxpointer(mgridlev)%array, velypointer(mgridlev)%array)
    call VecRestoreArrayReadF90(snvec(mgridlev), s, ierr)

    !Next curl on boundary points
    call DMGlobalToLocalEnd(das, snvec(mgridlev), INSERT_VALUES, svec_local, ierr)
    call VecGetArrayF90(svec_local, s_local, ierr)    
    call curl_bc(s_local, sbc_local, velxpointer(mgridlev)%array, velypointer(mgridlev)%array)
    
    !Get boundary conditions for the next finer grid level
    call get_bc(s_local, streambcvec_local, streambcvec_global, 1.d0)
    call VecAssemblyBegin(streambcvec_global, ierr)
    call VecRestoreArrayF90(svec_local, s_local, ierr)
    call VecRestoreArrayReadF90(streambcvec_local, sbc_local, ierr)
    deallocate(vort2pointer(mgridlev)%array2)
    !-----------------------------------------------------
    
    do k=mgridlev-1,1,-1
      
      !Performing apply_bc to include streamfunction boundary condition---------------------------
      call VecAssemblyEnd(streambcvec_global, ierr)
      call VecScatterBegin(ctx_bc, streambcvec_global, streambcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecRestoreArrayF90(velxvec(k+1), velxpointer(k+1)%array, ierr)
      call VecRestoreArrayF90(velyvec(k+1), velypointer(k+1)%array, ierr) 
      call VecScatterEnd(ctx_bc, streambcvec_global, streambcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr) 
      call VecGetArrayReadF90(streambcvec_local, sbc_local, ierr)
      call apply_bc(vort2pointer(k)%array2, sbc_local, 1.d0)!do applybc
      !------------------------------------------
      
      !ctci------------------------------
      call ctci(xsi, xse, ysi, yse, vort2pointer(k)%array2, spointer(k)%array)!do ctci
      call VecRestoreArrayF90(snvec(k), spointer(k)%array, ierr)
      call DMGlobalToLocalBegin(das, snvec(k), INSERT_VALUES, svec_local, ierr)
      !------------------------------------------
      
      !curl--------------------------------      
      call VecGetArrayReadF90(snvec(k), s, ierr)
      call curl_int(s, sbc_local, velxpointer(k)%array, velypointer(k)%array)
      call VecRestoreArrayReadF90(snvec(k), s, ierr)
      
      call DMGlobalToLocalEnd(das, snvec(k), INSERT_VALUES, svec_local, ierr)
      call VecGetArrayF90(svec_local, s_local, ierr)    !This s_local has ghost points as well
      call curl_bc(s_local, sbc_local, velxpointer(k)%array, velypointer(k)%array)!do curl
      
      if (k.gt.1) then
      call get_bc(s_local, streambcvec_local, streambcvec_global, 1.d0)
      call VecAssemblyBegin(streambcvec_global, ierr)
      end if
      
      call VecRestoreArrayF90(svec_local, s_local, ierr)
      call VecRestoreArrayReadF90(streambcvec_local, sbc_local, ierr)
      deallocate(vort2pointer(k)%array2)
      !------------------------------------------
      
    end do
    
    k = 1
    call VecRestoreArrayF90(velxvec(k), velxpointer(k)%array, ierr)
    call VecRestoreArrayF90(velyvec(k), velypointer(k)%array, ierr)  
  
  end subroutine vort2flux
  
!===================================================================================

  subroutine coarsify(vortvec2, vort1)
  
   !***************************************************************!
   !*   given vorticity on a smaller, fine mesh, (vort1) interp.    *!
   !*   values to the center region of a larger, coarser mesh     *!
   !*   (vort2).  The values outside the center region are         *!
   !*   not unmodified. Result is placed in arhs                  *!
   !***************************************************************!
  
    Vec :: vortvec2
    PetscScalar :: vort1(xsgi-1:xsge-1,ysgi-1:ysge-1)  !Note the change in indices only for convenience.
    PetscInt :: id, jd, i, j, iter
    PetscInt :: ix(tcoarse)
    PetscScalar :: yy(tcoarse)
    AO :: ao
    
    iter=0
    do j=yci,yse-1,2   !here we are looping over the finer grid, calculating the values and then
      do i=xci,xse-1,2  !setting the values at appropriate locations on the coarse grid
      
      id=m/4+i/2
      jd=n/4+j/2
      iter=iter+1
      ix(iter)=(m-1)*(jd-1)+(id-1)
      yy(iter)= vort1(i  ,j)   + &
                     0.5d0*( vort1(i+1,j)   + vort1(i  ,j+1)   + &
                             vort1(i-1,j)   + vort1(i  ,j-1) ) + &
                    0.25d0*( vort1(i+1,j+1) + vort1(i+1,j-1)   + &
                             vort1(i-1,j-1) + vort1(i-1,j+1) )
      
      end do
    end do  
    
    call DMDAGetAO(das, ao, ierr)
    call AOApplicationToPetsc(ao, tcoarse, ix, ierr)
    
    call VecSetValues(vortvec2, tcoarse, ix, yy, INSERT_VALUES, ierr)

  end subroutine coarsify

!===================================================================================

  subroutine curl_int(s_local, sbc_local, velx, vely)
  
    !***************************************************************!
    !*   returns curl(x) given x and the values of s'fun on bdy    *!
    !***************************************************************!
  
    PetscScalar :: s_local(xsi:xse,ysi:yse), sbc_local(nbc), velx(xui:xue,yui:yue), vely(xvi:xve,yvi:yve)
    PetscInt :: i,j
    
    forall(i=xui+1:xue-1,j=yui+1:yue-1); velx(i,j) = s_local(i,j+1) - s_local(i,j); end forall
    forall(i=xvi+1:xve-1,j=yvi+1:yve-1); vely(i,j) = s_local(i,j) - s_local(i+1,j); end forall
    
  end subroutine curl_int

!===================================================================================

  subroutine curl_bc(s_local, sbc_local, velx, vely)
  
    !***************************************************************!
    !*   returns curl(x) given x and the values of s'fun on bdy    *!
    !***************************************************************!
  
    PetscScalar :: s_local(xsgi:xsge,ysgi:ysge), sbc_local(nbc), velx(xui:xue,yui:yue), vely(xvi:xve,yvi:yve)
    PetscInt :: i,j
    
    if (xsgi==1)   then; forall(i=xsgi:xsgi, j=ysgi:ysge);  s_local(i,j)=sbc_local(left+j);   end forall; end if
    if (xsge==m+1) then; forall(i=xsge:xsge, j=ysgi:ysge);  s_local(i,j)=sbc_local(right+j);  end forall; end if
    if (ysgi==1)   then; forall(i=xsgi:xsge, j=ysgi:ysgi);  s_local(i,j)=sbc_local(bottom+i); end forall; end if
    if (ysge==n+1) then; forall(i=xsgi:xsge, j=ysge:ysge);  s_local(i,j)=sbc_local(top+i);    end forall; end if 
    
    do i=xui,xue
      j=yui; velx(i,j) = s_local(i,j+1) - s_local(i,j)
      j=yue; velx(i,j) = s_local(i,j+1) - s_local(i,j)
    end do
    do j=yui,yue
      i=xui; velx(i,j) = s_local(i,j+1) - s_local(i,j)
      i=xue; velx(i,j) = s_local(i,j+1) - s_local(i,j)
    end do
    
    do i=xvi,xve
      j=yvi; vely(i,j) = s_local(i,j) - s_local(i+1,j)
      j=yve; vely(i,j) = s_local(i,j) - s_local(i+1,j)
    end do
    do j=yvi,yve
      i=xvi; vely(i,j) = s_local(i,j) - s_local(i+1,j)
      i=xve; vely(i,j) = s_local(i,j) - s_local(i+1,j)
    end do
    
  end subroutine curl_bc

!===================================================================================

  subroutine get_bc(s_local, sbcvec_local, sbcvec_global, fac)
  
    !***************************************************************!
    !*   given vorticity on a larger, coarser mesh, interpolate    *!
    !*   it's values to the edge of a smaller, finer mesh          *!
    !***************************************************************!

    !The vorticity on s_local at the boundary of the finer mesh is saved to sbcvec_global, which is a global vector of size of the number of boundary points. sbcvec_global is then scattered to local boundary vectors sbcvec_local of the same size as sbcvec_global. This scattering is performed outside this subroutine.
    
    PetscScalar :: s_local(xsgi:xsge,ysgi:ysge)
    Vec :: sbcvec_local, sbcvec_global
    PetscScalar :: fac
    PetscScalar :: yy(tbc)
    PetscInt :: i, j, next
    
    next=0
    do i=xsi,xse
      do j=ysi,yse
      
        if (j==n/4+1   .and. i>=m/4+1 .and. i<=3*m/4+1) then
          next=next+1
          yy(next) = s_local(i,j)
          if (i.ne.3*m/4+1) then
            next=next+1
            yy(next) = 0.5*(s_local(i,j)+s_local(i+1,j))
          end if
        end if
        
        if (j==3*n/4+1 .and. i>=m/4+1 .and. i<=3*m/4+1) then
          next=next+1
          yy(next) = s_local(i,j)
          if (i.ne.3*m/4+1) then
            next=next+1
            yy(next) = 0.5*(s_local(i,j)+s_local(i+1,j))
          end if
        end if
        
        if (i==m/4+1   .and. j>=n/4+1 .and. j<=3*n/4+1) then
          next=next+1
          yy(next) = s_local(i,j)
          if (j.ne.3*n/4+1) then
            next=next+1
            yy(next) = 0.5*(s_local(i,j)+s_local(i,j+1))
          end if
        end if
        
        if (i==3*m/4+1 .and. j>=n/4+1 .and. j<=3*n/4+1) then
          next=next+1
          yy(next) = s_local(i,j)
          if (j.ne.3*n/4+1) then
            next=next+1
            yy(next) = 0.5*(s_local(i,j)+s_local(i,j+1))
          end if
        end if
      
      end do
    end do
    yy=yy*fac
    call VecSetValues(sbcvec_global, tbc, cbcix, yy, INSERT_VALUES, ierr)
     
  end subroutine get_bc

!===================================================================================

  subroutine apply_bc(vort, sbc_local, fac)
  
    PetscScalar :: vort(xsi:xse,ysi:yse), sbc_local(1:nbc), fac
    PetscInt :: i, j
    
    if (xsi==2) then; i=2; forall(j=ysi:yse);  vort(i,j)=vort(i,j)+fac*sbc_local(left+j);   end forall; end if;
    if (xse==m) then; i=m; forall(j=ysi:yse);  vort(i,j)=vort(i,j)+fac*sbc_local(right+j);   end forall; end if;
    if (ysi==2) then; j=2; forall(i=xsi:xse);  vort(i,j)=vort(i,j)+fac*sbc_local(bottom+i);   end forall; end if;
    if (yse==n) then; j=n; forall(i=xsi:xse);  vort(i,j)=vort(i,j)+fac*sbc_local(top+i);   end forall; end if;    
    
  end subroutine apply_bc

!===================================================================================

  subroutine nonlinear(omega, qx, qy, bc, nlx, nly)
  !Computes the nonlinear terms
  
    PetscScalar :: omega(xsgi:xsge,ysgi:ysge)
    PetscScalar :: qx(xugi:xuge,yugi:yuge)
    PetscScalar :: qy(xvgi:xvge,yvgi:yvge)
    PetscScalar :: bc(nbc), nlx(xui:xue,yui:yue), nly(xvi:xve,yvi:yve)
    PetscScalar :: uavg(xvi:xue,yui:yve), vavg(xvi:xue,yui:yve)
    PetscInt :: i, j
    
    if (xsgi==1)   then; forall(i=xsgi:xsgi, j=ysgi:ysge);  omega(i,j)=bc(left+j);   end forall; end if
    if (xsge==m+1) then; forall(i=xsge:xsge, j=ysgi:ysge);  omega(i,j)=bc(right+j);  end forall; end if
    if (ysgi==1)   then; forall(i=xsgi:xsge, j=ysgi:ysgi);  omega(i,j)=bc(bottom+i); end forall; end if
    if (ysge==n+1) then; forall(i=xsgi:xsge, j=ysge:ysge);  omega(i,j)=bc(top+i);    end forall; end if
    
    do j=yui,yve
      do i=xvi,xue
        uavg(i,j)=0.5d0*( qx(i,j)+qx(i,j-1) )
        vavg(i,j)=0.5d0*( qy(i,j)+qy(i-1,j) )
      end do
    end do
    
    do j=yui,yue
      do i=xui,xue
        nlx(i,j)=0.5d0*( vavg(i,j+1)*omega(i,j+1) + vavg(i,j)*omega(i,j) )
      end do
    end do
    
    do j=yvi,yve
      do i=xvi,xve
        nly(i,j)=-0.5d0*( uavg(i+1,j)*omega(i+1,j) + uavg(i,j)*omega(i,j) )
      end do
    end do
    
  end subroutine nonlinear

!===================================================================================

  subroutine rhs_forcing(itime, qx, qy, nlx, nly)
  
    use user
    
    PetscInt :: itime
    PetscScalar :: qx(xugi:xuge,yugi:yuge), qy(xvgi:xvge,yvgi:yvge)
    PetscScalar :: nlx(xui:xue,yui:yue), nly(xvi:xve,yvi:yve)
    PetscScalar :: dqx(xsgi:xse,ysgi:yse), dqy(xsgi:xse,ysgi:yse) !these exist on cell centers
    PetscInt :: i,j
    PetscScalar :: xx, yy, uu, vv
    PetscInt :: xuin, xuen, yuin, yuen, xvin, xven, yvin, yven
    
    dqx=0.d0
    dqy=0.d0
    
    do j=ysgi,yse
      do i=xsgi,xse
        xx = delta*(REAL(i)-0.5d0) - offsetx
        yy = delta*(REAL(j)-0.5d0) - offsety
        uu = 0.5d0 * ( qx(i,j) + qx(i+1,j) ) / delta
        vv = 0.5d0 * ( qy(i,j) + qy(i,j+1) ) / delta
        dqx(i,j)=delta* bodyforcex(itime,xx,yy,uu,vv)
        dqy(i,j)=delta* bodyforcey(itime,xx,yy,uu,vv)
      end do
    end do
    
    !Need forcing at boundary to be zero
    if (xsgi==1) then; dqx(xsgi,:)=0; dqy(xsgi,:)=0; end if
    if (xse==m)  then; dqx(xse,:)=0;  dqy(xse,:)=0;  end if
    if (ysgi==1) then; dqx(:,ysgi)=0; dqy(:,ysgi)=0; end if
    if (yse==n)  then; dqx(:,yse)=0;  dqy(:,yse)=0;  end if
    
    if (xui==1)   then; xuin=xui+1; else; xuin=xui; end if
    if (xue==m+1) then; xuen=xue-1; else; xuen=xue; end if
    if (yui==1)   then; yuin=yui+1; else; yuin=yui; end if
    if (yue==n)   then; yuen=yue-1; else; yuen=yue; end if
    
    if (xvi==1)   then; xvin=xvi+1; else; xvin=xvi; end if
    if (xve==m)   then; xven=xve-1; else; xven=xve; end if
    if (yvi==1)   then; yvin=yvi+1; else; yvin=yvi; end if
    if (yve==n+1) then; yven=yve-1; else; yven=yve; end if
    
    do j=yuin,yuen
      do i=xuin,xuen
        nlx(i,j)=nlx(i,j)+0.5d0*( dqx(i-1,j)+dqx(i,j) )
      end do
    end do
    
    do j=yvin,yven
      do i=xvin,xven
        nly(i,j)=nly(i,j)+0.5d0*( dqy(i,j-1)+dqy(i,j) )
      end do
    end do
  
  end subroutine

!===================================================================================

  subroutine fluid_predictor(omega, rhs, rhs_old, rhsbc, k)
  !Predictor step for vorticity
    
    use myfft
    
    PetscScalar :: omega(xsi:xse,ysi:yse), rhs(xsi:xse,ysi:yse), rhs_old(xsi:xse,ysi:yse,mgridlev), rhsbc(xsi:xse,ysi:yse)
    PetscInt :: k
    
    omega(:,:) = dst( xsi, xse, ysi, yse, lam1i(:,:,k) * &
     ( dst( xsi, xse, ysi, yse, con1(k)*rhs(:,:)       + &
                                con2(k)*rhs_old(:,:,k) + &
                                        rhsbc(:,:)   ) + & 
                lam1(:,:,k)*dst(xsi, xse, ysi, yse, omega(:,:)) ) )
    
  end subroutine fluid_predictor

!===================================================================================

  subroutine fluid_corrector(omega, vort)
    
    PetscScalar :: omega(xsi:xse,ysi:yse), vort(xsi:xse,ysi:yse)
    
    omega(:,:) = omega(:,:) - vort(:,:)
    
    
  end subroutine fluid_corrector

!===================================================================================

  subroutine aggregate_rhs(rhsfluid, rhs, rhs_old, rhsbc, k)
  
    PetscScalar :: rhsfluid(xsi:xse,ysi:yse), rhs(xsi:xse,ysi:yse), rhs_old(xsi:xse,ysi:yse,mgridlev), rhsbc(xsi:xse,ysi:yse)
    PetscInt :: k
    
    rhsfluid(:,:) = const1(k)*rhs(:,:) + const2(k)*rhs_old(:,:,k) + rhsbc(:,:)
  
  end subroutine aggregate_rhs

!===================================================================================

  FUNCTION delta_angle( it ) RESULT(ang)

    USE user
 
    PetscInt :: it
    PetscScalar :: ang, uv(5), k1, k2
    
    uv = motion_grid(it)
    k1 = dt*uv(3)
    uv = motion_grid(it+1)
    k2 = dt*uv(3)
    ang = 0.5d0*(k1+k2)

  END FUNCTION delta_angle
  
!===================================================================================

  subroutine setup_linear_operators
  
      PetscInt :: row, col(5), i, j, nrow, ncol, ione
      PetscReal :: value(5), one

      nrow = 1
      ione = 1
      one = 1

      do j=ysi, yse
        row = (j-ysgi)*xsgm + (xsi-xsgi) - 1
        do i=xsi, xse
          row = row+1
          !interior points
          if (i.ne.2 .and. i.ne.m .and. j.ne.2 .and. j.ne.n) then
              ncol = 5
              col(1) = row - xsgm
              col(2) = row - 1
              col(3) = row
              col(4) = row + 1
              col(5) = row + xsgm

              value(1) = -1
              value(2) = -1
              value(3) = 4
              value(4) = -1
              value(5) = -1

              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if

          if (i==2 .and. j.ne.2 .and. j.ne.n) then
              ncol = 4
              col(1) = row - xsgm
              col(2) = row
              col(3) = row + 1
              col(4) = row + xsgm

              value(1) = -1
              value(2) = 4
              value(3) = -1
              value(4) = -1

              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if

          if (i==m .and. j.ne.2 .and. j.ne.n) then
              ncol = 4
              col(1) = row - xsgm
              col(2) = row - 1
              col(3) = row
              col(4) = row + xsgm
      
              value(1) = -1
              value(2) = -1
              value(3) = 4
              value(4) = -1
      
              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if
      
          if (j==2 .and. i.ne.2 .and. i.ne.m) then
              ncol = 4
              col(1) = row - 1
              col(2) = row
              col(3) = row + 1
              col(4) = row + xsgm

              value(1) = -1
              value(2) = 4
              value(3) = -1
              value(4) = -1

              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if

          if (j==n .and. i.ne.2 .and. i.ne.m) then
              ncol = 4
              col(1) = row - xsgm
              col(2) = row - 1
              col(3) = row
              col(4) = row + 1

              value(1) = -1
              value(2) = -1
              value(3) = 4
              value(4) = -1

              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if
          
          !4 corners
          if (i==2 .and. j==2) then
              ncol = 3
              col(1) = row
              col(2) = row + 1
              col(3) = row + xsgm

              value(1) = 4
              value(2) = -1
              value(3) = -1

              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if

          if (i==2 .and. j==n) then
              ncol = 3
              col(1) = row - xsgm
              col(2) = row
              col(3) = row + 1

              value(1) = -1
              value(2) = 4
              value(3) = -1

              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if

          if (i==m .and. j==2) then
              ncol = 3
              col(1) = row - 1
              col(2) = row
              col(3) = row + xsgm

              value(1) = -1
              value(2) = 4
              value(3) = -1

              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if
      
          if (i==m .and. j==n) then
              ncol = 3
              col(1) = row - xsgm
              col(2) = row - 1
              col(3) = row

              value(1) = -1
              value(2) = -1
              value(3) = 4
      
              call MatSetValuesLocal(CTCmat, nrow, row, ncol, col, value, INSERT_VALUES, ierr)
          end if

          call MatSetValuesLocal(Imat, nrow, row, ione, row, one, INSERT_VALUES, ierr)
        end do
      end do
      
      call MatAssemblyBegin(Imat, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyBegin(CTCmat, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(Imat, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(CTCmat, MAT_FINAL_ASSEMBLY, ierr)
      
      call MatSetOption(Imat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE,ierr)
      call MatSetOption(CTCmat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE,ierr)
  
  end subroutine setup_linear_operators
  
!===================================================================================

  subroutine setup_constants
      
      PetscInt :: k
      PetscScalar :: del2, del22
  
      ALLOCATE( pvfac(mgridlev), mvfac(mgridlev), const1(mgridlev), const2(mgridlev) )
  
      del2 = delta*delta
      do k=1, mgridlev
          del22      =  del2*4.d0**(k-1)
          pvfac(k)    =  0.5d0*dt/Re/del22
          mvfac(k)    =  -0.5d0*dt/Re/del22
          const1(k)    =  1.5d0*dt/del22
          const2(k)    = -0.5d0*dt/del22
      end do
      
  end subroutine setup_constants  
  
!===================================================================================

!Additional subroutined required if using 2D domain decomposition

!===================================================================================

  subroutine a_times3(xvec, yy1vec, yy2vec)
  
  !This is the version of a_times2 but without vort2flux and fftw, to be used for CG etc. for 2D domain decomposition
  
    use myfft
    use fsinterface
    use operators_struct
    
    !Mat :: A
    Vec :: xvec, yy1vec, yy2vec
    PetscScalar, pointer :: x(:), hx(:), hy(:), vort(:)
    PetscScalar, pointer :: sbc_local(:), s(:), velx(:), vely(:), s_local(:)
    PetscInt :: k
    PetscScalar, pointer :: velx_interface(:), vely_interface(:), y(:)
    
    !Scatter level 1 begin
    call VecScatterBegin(ctx_force, xvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
        
    call VecZeroEntries(yy1vec, ierr)
    call VecZeroEntries(yy2vec, ierr)

    do k=1,mgridlev
      call VecZeroEntries(vortvecs(k), ierr)
    end do
    call VecZeroEntries(rhsvortvec, ierr)
    call VecZeroEntries(hxvec, ierr)
    call VecZeroEntries(hyvec, ierr)
    
    !Scatter level 1 end
    call VecScatterEnd(ctx_force, xvec, xivec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !Scatter level 2 begin
    call VecScatterBegin(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    
    !Doing reg---------------------------------
    call VecGetArrayReadF90(xvec, x, ierr)
    call reg(x, hxvec, hyvec)
    call VecAssemblyBegin(hxvec, ierr)
    call VecAssemblyBegin(hyvec, ierr)
    !Scatter level 2 end
    call VecScatterEnd(ctx_body, xivec, xbody_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)  
    call VecRestoreArrayReadF90(xvec, x, ierr)
    !-------------------------------------------
    
    !Doing b_times------------------------------
    call b_times2(xbody_vec, ybody_vec)
    !Scatter level 2 begin
    call VecScatterBegin(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !-------------------------------------------

    !Doing rot--------------------------------- 
    call VecGetArrayF90(rhsvortvec, vort, ierr)   
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
    call VecRestoreArrayF90(rhsvortvec, vort, ierr)
    !--------------------------------------------
    
    !Doing ainv----------------------------------
    call KSPSolve(vort_solver(1)%ksp, rhsvortvec, vortvecs(1), ierr)
    !call ainv(xsi, xse, ysi, yse, vort)
    !--------------------------------------------
    !Scatter level 2 end
    call VecScatterEnd(ctx_body, ybody_vec, yivec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !Scatter level 1 begin
    call VecScatterBegin(ctx_force, yivec, yy2vec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !--------------------------------------------
    
    !Doing vort2flux only on the finest grid-----------------------------
    call KSPSolve(ksp_sf, vortvecs(1), snvecs(1), ierr)
    
    !Initiate boundary vectors
    call DMGlobalToLocalBegin(das, snvecs(1), INSERT_VALUES, svec_local, ierr)
    call VecGetArrayReadF90(streambcvec_local, sbc_local, ierr)
    sbc_local = 0.d0
    
    !doing curl
    call VecGetArrayReadF90(snvecs(1), s, ierr)
    call VecGetArrayF90(velxvecs(1), velx, ierr)
    call VecGetArrayF90(velyvecs(1), vely, ierr)
    call curl_int(s, sbc_local, velx, vely)
    call VecRestoreArrayReadF90(snvecs(1), s, ierr)
    
    call DMGlobalToLocalEnd(das, snvecs(1), INSERT_VALUES, svec_local, ierr)
    call VecGetArrayF90(svec_local, s_local, ierr)    !This s_local has ghost points as well
    call curl_bc(s_local, sbc_local, velx, vely)
    
    call VecRestoreArrayF90(svec_local, s_local, ierr)
    call VecGetArrayF90(velxvecs(1), velx, ierr)
    call VecGetArrayF90(velyvecs(1), vely, ierr)
    call VecRestoreArrayReadF90(streambcvec_local, sbc_local, ierr)
    !end vort2flux on the finest grid--------------------------------------------
    
    !Doing regT----------------------------------
    
    call VecScatterBegin(ctx_velx, velxvecs(1), velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(ctx_vely, velyvecs(1), velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayF90(yy1vec, y, ierr)
    call VecScatterEnd(ctx_velx, velxvecs(1), velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velxvec_interface, velx_interface, ierr)
    call VecScatterEnd(ctx_vely, velyvecs(1), velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velyvec_interface, vely_interface, ierr)
    
    call regT(velx_interface, vely_interface, y)!do regT
    call VecRestoreArrayF90(yy1vec, y, ierr)
    call VecRestoreArrayReadF90(velxvec_interface, velx_interface, ierr)
    call VecRestoreArrayReadF90(velyvec_interface, vely_interface, ierr)
    !--------------------------------------------
    !Scatter level 1 end
    call VecScatterEnd(ctx_force, yivec, yy2vec, INSERT_VALUES, SCATTER_REVERSE, ierr)
    !--------------------------------------------
    
  end subroutine a_times3

!===================================================================================

  subroutine a_times(xvec, yvec)
  !this is the original a_times subroutine that does not include b_times
  
    use myfft
    use fsinterface
    
    !Mat :: A
    Vec :: xvec, yvec
    !Vec :: hxvec, hyvec, hxvec_local, hyvec_local
    PetscScalar, pointer :: x(:), hx(:), hy(:), vort(:)
    PetscScalar, pointer :: sbc_local(:), s(:), velx(:), vely(:), s_local(:)
    !Vec, pointer :: velxvecs(:), velyvecs(:)!, vortvecs(:), snvecs(:)
    PetscInt :: k
    !IS :: reg_ixs, reg_iys
    !Vec :: velxvec_interface, velyvec_interface
    !VecScatter :: ctx_velx, ctx_vely
    PetscScalar, pointer :: velx_interface(:), vely_interface(:), y(:)
    !PetscErrorCode :: ierr
        
    do k=1,mgridlev
      call VecZeroEntries(vortvecs(k), ierr)
    end do
    call VecZeroEntries(hxvec, ierr)
    call VecZeroEntries(hyvec, ierr)
    call VecZeroEntries(yvec, ierr)
    
    !Doing reg---------------------------------
    call VecGetArrayReadF90(xvec, x, ierr)
    call reg(x, hxvec, hyvec)
    call VecAssemblyBegin(hxvec, ierr)
    call VecAssemblyBegin(hyvec, ierr)  
    call VecRestoreArrayReadF90(xvec, x, ierr)
    !-------------------------------------------

    !Doing rot--------------------------------- 
    call VecGetArrayF90(vortvecs(1), vort, ierr)   
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
    call ainv(xsi, xse, ysi, yse, vort)
    call VecRestoreArrayF90(vortvecs(1), vort, ierr)
    !--------------------------------------------
    
    !Doing vort2flux-----------------------------
    !either use vort2flux subroutine
    if (mgridlev .ge. 2) then
      CALL vort2flux(velxvecs, velyvecs, vortvecs, snvecs, mgridlev)!what is the difference between if/else here?
    else
      call vort2flux(velxvecs, velyvecs, vortvecs, snvecs, mgridlev)
    end if
    !end vort2flux--------------------------------------------
    
    !Doing regT----------------------------------    
    call VecScatterBegin(ctx_velx, velxvecs(1), velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterBegin(ctx_vely, velyvecs(1), velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayF90(yvec, y, ierr)
    call VecScatterEnd(ctx_velx, velxvecs(1), velxvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velxvec_interface, velx_interface, ierr)
    call VecScatterEnd(ctx_vely, velyvecs(1), velyvec_interface, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGetArrayReadF90(velyvec_interface, vely_interface, ierr)
    
    call regT(velx_interface, vely_interface, y)!do regT
    call VecRestoreArrayF90(yvec, y, ierr)
    call VecRestoreArrayReadF90(velxvec_interface, velx_interface, ierr)
    call VecRestoreArrayReadF90(velyvec_interface, vely_interface, ierr)
    !--------------------------------------------

  end subroutine a_times

!===================================================================================

  subroutine vort2flux2(velxvec, velyvec, vortvec, snvec, mgridlev)
  !this is vort2flux for 2D domain decomposition using CG
  
    !***************************************************************!
    !*    Multiscale method to solve C^T C s = omega               *!
    !*    and return the velocity, C s.                            *!
    !*    Results are returned in vel on each of the first nlev    *!
    !*    grids.                                                   *!
    !*    Warning: the vorticity field on all but the finest mesh  *!
    !*     is modified by the routine in the following way:        *!
    !*     the value in the center of the domain is interpolated   *!
    !*     from the next finer mesh (the value near the edge is    *!
    !*     not changed.                                            *!
    !***************************************************************!
    
    use myfft
    
    PetscInt :: mgridlev, k
    Vec, Pointer :: velxvec(:), velyvec(:), vortvec(:), snvec(:)
    !Vec :: vortvec_local1, svec_local
    PetscScalar, Pointer :: vort1(:), vort(:), s(:), s_local(:), velx(:), vely(:), rhsstream(:)
    !Vec :: sbcvec_global, sbcvec_local
    PetscScalar, Pointer :: sbc_local(:)
    PetscScalar :: vort2(xsm*ysm)
    
    type pointertoarray2
     PetscScalar, allocatable :: array2(:)
    end type pointertoarray2
    
    type(pointertoarray2), dimension(mgridlev) :: vort2pointer
    
    !Doing coarsify---------------------------------
    !call DMGetLocalVector(das, vortvec_local1, ierr)
    if (mgridlev .ge. 2) then
      DO k=2,mgridlev

        call VecCopy(vortvec(k-1), streamvec_coarsify, ierr)
        call DMGlobalToLocalBegin(das_coarsify, streamvec_coarsify, INSERT_VALUES, vortvec_local1, ierr)
        call VecGetArrayF90(velxvec(k), velxpointer(k)%array, ierr)
        call VecGetArrayF90(velyvec(k), velypointer(k)%array, ierr)        
        
        call DMGlobalToLocalEnd(das_coarsify, streamvec_coarsify, INSERT_VALUES, vortvec_local1, ierr)
        call VecGetArrayReadF90(vortvec_local1, vort1, ierr)
        
        call coarsify(vortvec(k), vort1)
        call VecAssemblyBegin(vortvec(k), ierr)
        
        call VecRestoreArrayReadF90(vortvec_local1, vort1, ierr)
        allocate(vort2pointer(k-1)%array2(xsm*ysm))
        call VecGetArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
        vort2pointer(k-1)%array2 = vortpointer(k-1)%array
        call VecRestoreArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
        call VecAssemblyEnd(vortvec(k), ierr)  
        
      END DO
      
      k = 1
      call VecGetArrayF90(velxvec(k), velxpointer(k)%array, ierr)
      call VecGetArrayF90(velyvec(k), velypointer(k)%array, ierr)
      
      k = mgridlev+1
      allocate(vort2pointer(k-1)%array2(xsm*ysm))
      call VecGetArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
      vort2pointer(k-1)%array2 = vortpointer(k-1)%array
      call VecRestoreArrayReadF90(vortvec(k-1), vortpointer(k-1)%array, ierr)
      
    end if
    !---------------------------------------------------
    
    !Solving CTCs=vort+bc-----------------------------------------
    call VecGetArrayF90(rhsstream_vec, rhsstream, ierr)
    rhsstream = vort2pointer(mgridlev)%array2
    call VecRestoreArrayF90(rhsstream_vec, rhsstream, ierr)
    call KSPSolve(ksp_sf, rhsstream_vec, snvec(mgridlev), ierr)
    !-----------------------------------------------------
    
    !Initiate boundary vectors----------------------------
    call DMGlobalToLocalBegin(das, snvec(mgridlev), INSERT_VALUES, svec_local, ierr)
    call VecGetArrayReadF90(streambcvec_local, sbc_local, ierr)
    sbc_local = 0.d0
    !-----------------------------------------------------
    
    !Doing curl------------------------------------------    
    call VecGetArrayReadF90(snvec(mgridlev), s, ierr)
    call curl_int(s, sbc_local, velxpointer(mgridlev)%array, velypointer(mgridlev)%array)
    call VecRestoreArrayReadF90(snvec(mgridlev), s, ierr)
    
    call DMGlobalToLocalEnd(das, snvec(mgridlev), INSERT_VALUES, svec_local, ierr)
    call VecGetArrayF90(svec_local, s_local, ierr)    !This s_local has ghost points as well
    call curl_bc(s_local, sbc_local, velxpointer(mgridlev)%array, velypointer(mgridlev)%array)
    
    call get_bc(s_local, streambcvec_local, streambcvec_global, 1.d0)
    call VecAssemblyBegin(streambcvec_global, ierr)
    call VecRestoreArrayF90(svec_local, s_local, ierr)
    call VecRestoreArrayReadF90(streambcvec_local, sbc_local, ierr)
    deallocate(vort2pointer(mgridlev)%array2)
    !-----------------------------------------------------
    
    do k=mgridlev-1,1,-1
      
      !Doing apply_bc---------------------------
      call VecAssemblyEnd(streambcvec_global, ierr)
      call VecScatterBegin(ctx_bc, streambcvec_global, streambcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecRestoreArrayF90(velxvec(k+1), velxpointer(k+1)%array, ierr)
      call VecRestoreArrayF90(velyvec(k+1), velypointer(k+1)%array, ierr) 
      call VecScatterEnd(ctx_bc, streambcvec_global, streambcvec_local, INSERT_VALUES, SCATTER_FORWARD, ierr) 
      call VecGetArrayReadF90(streambcvec_local, sbc_local, ierr)
      call apply_bc(vort2pointer(k)%array2, sbc_local, 1.d0)!do applybc
      !------------------------------------------
      
      !Doing ctci------------------------------
      call VecGetArrayF90(rhsstream_vec, rhsstream, ierr)
      rhsstream = vort2pointer(k)%array2
      call VecRestoreArrayF90(rhsstream_vec, rhsstream, ierr)
      call KSPSolve(ksp_sf, rhsstream_vec, snvec(k), ierr)
      call DMGlobalToLocalBegin(das, snvec(k), INSERT_VALUES, svec_local, ierr)
      !------------------------------------------
      
      !Doing curl--------------------------------
      
      call VecGetArrayReadF90(snvec(k), s, ierr)
      call curl_int(s, sbc_local, velxpointer(k)%array, velypointer(k)%array)
      call VecRestoreArrayReadF90(snvec(k), s, ierr)
      
      call DMGlobalToLocalEnd(das, snvec(k), INSERT_VALUES, svec_local, ierr)
      call VecGetArrayF90(svec_local, s_local, ierr)    !This s_local has ghost points as well
      call curl_bc(s_local, sbc_local, velxpointer(k)%array, velypointer(k)%array)!do curl
      
      if (k.gt.1) then
      call get_bc(s_local, streambcvec_local, streambcvec_global, 1.d0)
      call VecAssemblyBegin(streambcvec_global, ierr)
      end if
      
      call VecRestoreArrayF90(svec_local, s_local, ierr)
      call VecRestoreArrayReadF90(streambcvec_local, sbc_local, ierr)
      deallocate(vort2pointer(k)%array2)
      !------------------------------------------
      
    end do
    
    k = 1
    call VecRestoreArrayF90(velxvec(k), velxpointer(k)%array, ierr)
    call VecRestoreArrayF90(velyvec(k), velypointer(k)%array, ierr) 
  
  end subroutine vort2flux2

!===================================================================================

end module operators_fluid
