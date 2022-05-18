MODULE variables

  USE petsc
#include "petsc/finclude/petsc.h"
  USE parameters
  USE grid
  use user
  IMPLICIT NONE

  Vec :: streamvec, streamvec_coarsify, velocityxvec, velocityyvec, streambcvec_local, streambcvec_global, velocityxvec_interface, velocityyvec_interface, rhsfluid_vec, rhsstream_vec !dummy vectors
  Vec :: streamvec_local, vortvec_local1, svec_local !dummy vectors with allocation for ghost values
  Vec, pointer :: svec(:), omegavec(:), omegavec_local(:), vortvecs(:), snvecs(:)
  Vec, pointer :: qxvec(:), q0xvec(:), q0pxvec(:), q0rxvec(:), qyvec(:), q0yvec(:), q0pyvec(:), q0ryvec(:), velxvecs(:), velyvecs(:)
  Vec, pointer :: streambcvecs_global(:), lastbcvec_local(:), omegabcvec_local(:)
  Vec :: qhxvec, qhyvec, qxvec_local, qyvec_local, nlxvec, nlyvec, nlxvec_local, nlyvec_local
  Vec :: velxvec_interface, velyvec_interface, hxvec, hyvec, vortvec, rhsvortvec, hxvec_local, hyvec_local
  PetscScalar, dimension(:,:,:), allocatable :: rhs_old
  Vec :: f_rdstvec, rhsfvec
  Mat :: cholmat
  Mat :: A, B, AB
  VecScatter :: ctx_bc, ctx_velx, ctx_vely, ctx_force, ctx_body, ctx_force2, ctx_body2
  KSP :: ksp, ksp_fsi, ksp_sf
  PC :: pc, pc_fsi, pc_sf
  Vec :: fivec !Sequential force vector for btimes
  Vec :: fbody_vec  !MPI vector on comm_body where internal matrix vector multiplication for structures happens
  Vec :: xivec, xbody_vec, yivec, ybody_vec, y1body_vec, y1ivec, y2body_vec, y2ivec, y1vec, y2vec !temporary vectors
  Vec :: frcxvec, frcyvec, wghtxvec, wghtyvec, finterx_vec, fintery_vec
  Mat, pointer :: sol_Jmat(:), sol_Kmat(:)  !This is the RBJBR matrix. It is created as a matrix pointer because a given proc may have many bodies
  Mat, pointer :: Kmat(:), Mdmat(:) !matrices for deformable modies
  PetscScalar, allocatable :: Mtmat(:), Jmat(:), Ftint(:), Ctmat(:)  !Mtmat is a single variable containing i_theta
  PetscScalar, allocatable :: theta(:), thetad(:), thetadd(:), dtheta(:)
  PetscScalar, allocatable :: theta0(:), thetad0(:), thetadd0(:)
  Vec :: dxbvec !this vector stores increment in location of body points
  PetscInt :: ksp_iter, iterations, fsi_iter, fluid_iters(2) !keeping track of iterations
  Mat :: CTCmat, Imat
  
  
  type pointertoarray
     PetscScalar, Pointer :: array(:)
  end type pointertoarray
  
  type(pointertoarray), allocatable, dimension(:) :: velxpointer, velypointer, spointer, vortpointer
  
  type matpointer
     Mat :: matrix
  end type matpointer
  
  type ksppointer
     KSP :: ksp
     PC :: pc
  end type ksppointer
    
  type(matpointer), allocatable, dimension(:) :: Amat, Amat2
  type(ksppointer), allocatable, dimension(:) :: vort_solver
  
  PetscScalar, allocatable, dimension(:) :: pvfac, mvfac, const1, const2
  
  Mat :: indicestransposemat, intermediatemat1, intermediatemat2, lumat, precomputemat2
  Vec :: velxvec_interface_sd, velyvec_interface_sd, rhsfsdvec
  VecScatter :: ctx_velx_sd, ctx_vely_sd
  KSP :: ksp_fsi_sd
  PC :: pc_fsi_sd
  PetscScalar, dimension(:,:), allocatable :: intermediateserialmat, luserialmat
  PetscInt, allocatable :: precomputeserialnz(:), precomputenz(:)

  type precomputeindicesmatpointer
     PetscInt, allocatable :: value(:)
  end type precomputeindicesmatpointer
  
  type precomputematpointer
     PetscScalar, allocatable :: value(:)
  end type precomputematpointer
  
  type(precomputeindicesmatpointer), allocatable :: precomputeserialindicesmat(:), precomputeindicesmat(:), precomputeserialindicesmat2(:)
  type(precomputematpointer), allocatable :: precomputeserialmat(:), precomputemat(:)
    
  VecScatter :: ctx_velx_sd_sparse, ctx_vely_sd_sparse
  Vec :: velxvec_interface_sd_sparse, velyvec_interface_sd_sparse
  
contains

!================================================================================

  SUBROUTINE setup_variables
  
    USE myfft
    
    PetscInt :: i
    Mat :: dummat
    PetscInt, allocatable :: lumatlocal_idx(:)
    PetscViewer :: viewer
  
    call setup_fft(xsi, xse, ysi, yse, xpi, xpe, ypi, ype, delta)
    
    !Streamfunction and vorticity vectors
    call DMCreateGlobalVector(das, streamvec, ierr)
    call DMCreateGlobalVector(das_coarsify, streamvec_coarsify, ierr)
    call VecDuplicateVecsF90(streamvec, mgridlev, svec, ierr)
    call VecDuplicateVecsF90(streamvec, mgridlev, omegavec, ierr)
    call DMGetLocalVector(das, streamvec_local, ierr) !vorticity dummy vectors which includes ghost values
    
    !Boundary condition vectors (these vectors are required only for creating scatter context in get_bc routine)
    call VecCreateSeq(MPI_COMM_SELF, nbc, streambcvec_local, ierr)
    call VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nbc, streambcvec_global, ierr)
    call VecZeroEntries(streambcvec_local, ierr)
    call VecZeroEntries(streambcvec_global, ierr)
    
    !Velocity-x vectors
    call DMCreateGlobalVector(dau, velocityxvec, ierr)
    call VecDuplicateVecsF90(velocityxvec, mgridlev, qxvec, ierr)
    call VecDuplicateVecsF90(velocityxvec, mgridlev, q0xvec, ierr)
    call VecDuplicateVecsF90(velocityxvec, mgridlev, q0pxvec, ierr)
    call VecDuplicateVecsF90(velocityxvec, mgridlev, q0rxvec, ierr)
    
    !Velocity-y vectors
    call DMCreateGlobalVector(dav, velocityyvec, ierr)
    call VecDuplicateVecsF90(velocityyvec, mgridlev, qyvec, ierr)
    call VecDuplicateVecsF90(velocityyvec, mgridlev, q0yvec, ierr)
    call VecDuplicateVecsF90(velocityyvec, mgridlev, q0pyvec, ierr)
    call VecDuplicateVecsF90(velocityyvec, mgridlev, q0ryvec, ierr)
    
    !Velocity vectors for reg and regT purposes
    call VecCreateSeq(MPI_COMM_SELF, interfacex, velocityxvec_interface, ierr)
    call VecCreateSeq(MPI_COMM_SELF, interfacey, velocityyvec_interface, ierr)
    if (sub_domain) then
        call VecCreateSeq(MPI_COMM_SELF, interfacex_sd, velxvec_interface_sd, ierr)
        call VecCreateSeq(MPI_COMM_SELF, interfacey_sd, velyvec_interface_sd, ierr)
    end if
    
    allocate(rhs_old(xsi:xse,ysi:yse,mgridlev))
    
    !variables used only as temporary placeholders
    call VecDuplicate(streamvec, rhsfluid_vec, ierr) !rhs term for fluid equation
    call VecDuplicate(streamvec, rhsstream_vec, ierr) !rhs term for getting streamfunciton from vorticity
    call VecDuplicateVecsF90(streambcvec_global, mgridlev, streambcvecs_global, ierr)  !boundary
    call VecDuplicateVecsF90(streambcvec_local, mgridlev, lastbcvec_local, ierr)  !boundary
    call VecDuplicateVecsF90(streambcvec_local, mgridlev, omegabcvec_local, ierr) !boundary
    call VecDuplicateVecsF90(streamvec_local, mgridlev, omegavec_local, ierr) !vorticity
    call VecDuplicate(velocityxvec, qhxvec, ierr)
    call VecDuplicate(velocityyvec, qhyvec, ierr)   
    call DMGetLocalVector(dau, qxvec_local, ierr) !x-velocity
    call DMGetLocalVector(dav, qyvec_local, ierr) !x-velocity     
    call VecDuplicate(velocityxvec, nlxvec, ierr) !x-nonlinear term
    call VecDuplicate(velocityyvec, nlyvec, ierr) !y-nonlinear term
    call DMGetLocalVector(dau, nlxvec_local, ierr) !x-nonlinear_local term with ghost
    call DMGetLocalVector(dav, nlyvec_local, ierr) !y-nonlinear_local term with ghost
    call VecDuplicate(fvec, rhsfvec, ierr)
    if (sub_domain) call VecDuplicate(xsdvec, rhsfsdvec, ierr)
    call VecDuplicate(velocityxvec_interface, velxvec_interface, ierr) !for regT
    call VecDuplicate(velocityyvec_interface, velyvec_interface, ierr) !for regT
    call VecDuplicate(velocityxvec, hxvec, ierr) !for reg
    call VecDuplicate(velocityyvec, hyvec, ierr) !for reg
    call VecDuplicate(streamvec, vortvec, ierr) !for rot after reg
    call VecDuplicate(streamvec, rhsvortvec, ierr) !rhs of last equation
    call DMGetLocalVector(dau, hxvec_local, ierr) !for rot after reg
    call DMGetLocalVector(dav, hyvec_local, ierr) !for rot after reg
    call VecDuplicateVecsF90(streamvec, mgridlev, vortvecs, ierr) !for a_times
    call VecDuplicateVecsF90(streamvec, mgridlev, snvecs, ierr) !for a_times
    call VecDuplicateVecsF90(velocityxvec, mgridlev, velxvecs, ierr) !for a_times
    call VecDuplicateVecsF90(velocityyvec, mgridlev, velyvecs, ierr) !for a_times
    call DMGetLocalVector(das, vortvec_local1, ierr) !for coarsify in vort2flux
    call DMGetLocalVector(das, svec_local, ierr) !for coarsify in vort2flux
    allocate(velxpointer(mgridlev))
    allocate(velypointer(mgridlev))
    allocate(spointer(mgridlev))
    allocate(vortpointer(mgridlev))
    

    !Cholesky matrices and vectors
    if (num_stat) then
       call MatCreate(MPI_COMM_WORLD, cholmat, ierr)
       call MatSetType(cholmat, MATELEMENTAL, ierr)
       call MatSetSizes(cholmat, xfm, xfm, nf, nf, ierr)
       call MatSetFromOptions(cholmat, ierr)
       call MatSetUp(cholmat, ierr)
    end if
    
    if (sub_domain) then
      
       !The entire if loop can be deleted......
       if (rank==0) then
           half_width = nint(support*2*delta/delta_sd-1)
           nz = half_width*2+1   !number of immediate neighboring elements: support*2*delta/delta_sd-1 is the number of elements on one side; it is multiplied by two to get elements on both sides and finally one is added to count the main element itself
           nz = nz*nz*2     !There are as many such clusters, therefore nz^2. Also, to incorporate corelation between and x and y off diagonal block matrices, 2 is multiplied.
           nz_indicesmat = interp_points
       else
           half_width = 0
           nz = 0
           nz_indicesmat = 0
       end if
       !Delete till here.........
       
       !lumat is the small-dimensional matrix in step 2 multiplied by fvec
       call MatCreate(MPI_COMM_WORLD, lumat, ierr)
       call MatSetType(lumat, MATDENSE, ierr)
       call MatSetSizes(lumat, xfm, xfm, nf, nf, ierr)
       call MatSetFromOptions(lumat, ierr)
       call MatSetUp(lumat, ierr)
       
       !Get the sparsity of EE^T matrix
       if (.not. sub_domain_precomputed) call get_sparsity_EEt

    end if

    call VecDuplicate(fvec, f_rdstvec, ierr)
    !Shell matrices for matrix-free implementation of GMRES
    call MatCreateShell(MPI_COMM_WORLD, xfm, xfm, nf, nf, PETSC_NULL_INTEGER, A, ierr)
    call MatCreateShell(MPI_COMM_WORLD, xfm, xfm, nf, nf, PETSC_NULL_INTEGER, B, ierr)
    call MatCreateShell(MPI_COMM_WORLD, xfm, xfm, nf, nf, PETSC_NULL_INTEGER, AB, ierr)
    
    !Body configuration variables for torsional body
    if (mt_body>0) then 
      allocate(theta(mt_body), thetad(mt_body), thetadd(mt_body), dtheta(mt_body))
      allocate(theta0(mt_body), thetad0(mt_body), thetadd0(mt_body))
    end if
    
    !initial conditions
    IF (istart==0) THEN
       CALL initial_condition
       call VecGetArrayF90(xbvec, xb, ierr)
       call VecGetArrayF90(ybvec, yb, ierr)
       call ISDestroy(reg_ixs, ierr)
       call ISDestroy(reg_iys, ierr)
       call reg_prep(xb, yb)
       call VecRestoreArrayF90(xbvec, xb, ierr)
       call VecRestoreArrayF90(ybvec, yb, ierr)
       if (mt_body>0) then
          theta=0.d0
          thetad=0.d0
       end if
       
       if (sub_domain) then
           call VecGetArrayReadF90(xsdvec, xsd, ierr)
           call VecGetArrayReadF90(ysdvec, ysd, ierr)
           call VecRestoreArrayReadF90(xsdvec, xsd, ierr)
           call VecRestoreArrayReadF90(ysdvec, ysd, ierr)
       end if
    ELSE
       call read_variables(istart)
       call VecGetArrayF90(xbvec, xb, ierr)
       call VecGetArrayF90(ybvec, yb, ierr)
       call ISDestroy(reg_ixs, ierr)
       call ISDestroy(reg_iys, ierr)
       call reg_prep(xb, yb)
       call VecRestoreArrayF90(xbvec, xb, ierr)
       call VecRestoreArrayF90(ybvec, yb, ierr)
    END IF
    
    call VecScatterCreate(streambcvec_global, bcis, streambcvec_local, bcis, ctx_bc, ierr)   !Scatter context for get_bc
    
    !call create_linear_operators !for 2D domain decomposition
    
    call setup_ksp
    call setup_body_dist
    call create_struct_matrix    
        
  END SUBROUTINE setup_variables


!================================================================================

  subroutine initial_condition

    PetscInt :: k
    LOGICAL :: readic
    PetscScalar :: uv(5)

    INQUIRE(file="input/initial.var",exist=readic)
    
    IF (readic) THEN
      !read files: left for later
    else
      
      rhs_old = 0.d0
      do k=1,mgridlev
        call VecZeroEntries(omegavec(k), ierr)
        call VecZeroEntries(qxvec(k), ierr)
        call VecZeroEntries(qyvec(k), ierr)
      end do
      
      rot_angle = 0.d0
      uv = motion_grid( 0 )
      rox = uv(4)
      roy = uv(5)
       
      call motion_potential(izero, q0pxvec, q0pyvec)
      call motion_rotation(izero, q0rxvec, q0ryvec) 
      do k=1,mgridlev
        call VecWAXPY(q0xvec(k), one, q0pxvec(k), q0rxvec(k), ierr)
        call VecWAXPY(q0yvec(k), one, q0pyvec(k), q0ryvec(k), ierr)
      end do        

    end if

  end subroutine initial_condition
  
!================================================================================

  subroutine motion_potential(it, qxref, qyref)
    
    PetscInt :: it, k
    Vec, pointer :: qxref(:), qyref(:)
    PetscScalar :: uv(5)
    PetscScalar :: fac
    
    uv=motion_grid(it)    
    do k=1,mgridlev
      fac = delta*2.d0**(k-1)
    
      call VecSet(qxref(k), -fac*uv(1), ierr)  !in new version of petsc qxref and fac*uv might be interchanged, BEWARE!!!!!!!!!!!
      call VecSet(qyref(k), -fac*uv(2), ierr)
    end do
  
  end subroutine motion_potential

!================================================================================

  subroutine motion_rotation(it, qxref, qyref)
  
    PetscInt :: it, k, iter, i, j
    Vec, pointer :: qxref(:), qyref(:)
    PetscScalar :: uv(5)
    PetscScalar :: omegab, xx, yy, fac
    PetscScalar, pointer :: qx(:), qy(:)
    
    uv=motion_grid(it) 
    omegab = uv(3)
    
    do k=1,mgridlev
      call VecGetArrayF90(qxref(k), qx, ierr)
      call VecGetArrayF90(qyref(k), qy, ierr)      
      qx=0.d0
      qy=0.d0     
      fac = delta*2.d0**(k-1)  ! cell face length on each grid
      iter = 0
      DO j=yui,yue
        DO i=xui,xue
          iter=iter+1
          xx = (REAL(i)-1-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
          yy = (REAL(j)-0.5D0-n/2)*delta*(2**(k-1)) + n/2*delta -offsety
          qx(iter) = qx(iter) + omegab*( yy - roy )
        END DO
      END DO
      iter=0
      DO j=yvi,yve
        DO i=xvi,xve
          iter=iter+1
          xx = (REAL(i)-0.5d0-m/2)*delta*(2**(k-1)) + m/2*delta - offsetx
          yy = (REAL(j)-1-n/2)*delta*(2**(k-1)) + n/2*delta -offsety
          qy(iter) = qy(iter) - omegab*( xx - rox )
        END DO
      END DO
      
      qx=fac*qx
      qy=fac*qy
      call VecRestoreArrayF90(qxref(k), qx, ierr)
      call VecRestoreArrayF90(qyref(k), qy, ierr)
    end do
      
  end subroutine motion_rotation

!================================================================================

  subroutine reg_prep(xb, yb)
  !Creates the scatter context for E and E^T
  
    PetscScalar :: xb(xfi:xfe), yb(xfi:xfe)
    
    call setup_reg(xb, yb)
    call VecScatterCreate(velocityxvec, reg_ixs, velocityxvec_interface, interfacex_is, ctx_velx, ierr)
    call VecScatterCreate(velocityyvec, reg_iys, velocityyvec_interface, interfacey_is, ctx_vely, ierr)
  
  end subroutine reg_prep  
  
!================================================================================

  subroutine reg_prep_sd_sparse(xb, yb)
  
    PetscScalar :: xb(xfi:xfe), yb(xfi:xfe)
    
    call ISDestroy(reg_ixs_sd_sparse, ierr)
    call ISDestroy(reg_iys_sd_sparse, ierr)
    call ISDestroy(interfacex_sd_is_sparse, ierr)
    call ISDestroy(interfacey_sd_is_sparse, ierr)
    call VecDestroy(velxvec_interface_sd_sparse, ierr)
    call VecDestroy(velyvec_interface_sd_sparse, ierr)
    call VecScatterDestroy(ctx_velx_sd_sparse, ierr)
    call VecScatterDestroy(ctx_vely_sd_sparse, ierr)
    
    call setup_reg_subdomain_indices(xb, yb) !creates new reg_ixs_sd_sparse and interfacex_sd_is_sparse
    call VecCreateSeq(MPI_COMM_SELF, interfacex_sd_sparse, velxvec_interface_sd_sparse, ierr)
    call VecCreateSeq(MPI_COMM_SELF, interfacey_sd_sparse, velyvec_interface_sd_sparse, ierr)
    call VecScatterCreate(velocityxvec, reg_ixs_sd_sparse, velxvec_interface_sd_sparse, interfacex_sd_is_sparse, ctx_velx_sd_sparse, ierr)
    call VecScatterCreate(velocityyvec, reg_iys_sd_sparse, velyvec_interface_sd_sparse, interfacey_sd_is_sparse, ctx_vely_sd_sparse, ierr)
  
  end subroutine reg_prep_sd_sparse
  
!================================================================================

  subroutine reg_prep_sd(sd_list_size, sd_list)
  
    integer, intent(in) :: sd_list_size
    PetscInt, intent(in) :: sd_list(sd_list_size)
    
    call ISDestroy(reg_ixs_sd_sparse, ierr)
    call ISDestroy(reg_iys_sd_sparse, ierr)
    call ISDestroy(interfacex_sd_is_sparse, ierr)
    call ISDestroy(interfacey_sd_is_sparse, ierr)
    call VecDestroy(velxvec_interface_sd_sparse, ierr)
    call VecDestroy(velyvec_interface_sd_sparse, ierr)
    call VecScatterDestroy(ctx_velx_sd_sparse, ierr)
    call VecScatterDestroy(ctx_vely_sd_sparse, ierr)
    
    
    call setup_reg_subdomain(sd_list_size, sd_list) !creates new reg_ixs_sd_sparse and interfacex_sd_is_sparse
    call VecCreateSeq(MPI_COMM_SELF, interfacex_sd_sparse, velxvec_interface_sd_sparse, ierr)
    call VecCreateSeq(MPI_COMM_SELF, interfacey_sd_sparse, velyvec_interface_sd_sparse, ierr)
    call VecScatterCreate(velocityxvec, reg_ixs_sd_sparse, velxvec_interface_sd_sparse, interfacex_sd_is_sparse, ctx_velx_sd_sparse, ierr)
    call VecScatterCreate(velocityyvec, reg_iys_sd_sparse, velyvec_interface_sd_sparse, interfacey_sd_is_sparse, ctx_vely_sd_sparse, ierr)
  
  end subroutine reg_prep_sd

!================================================================================
  
  subroutine get_sparsity_EEt
  !The sparsity of EE^T is determined here.
  !This subroutine could have been avoided entirely, but it is used to make the precomputation process faster.
  !The idea is as follows: To apply the droptol technique, we need to first compute the full column j of B. Communication of this column is required before applying the droptol technqiue. However, for very large problems, the column j can be of the size of a million and therefore, multiple evaluations can be prohibitive. Therefore, the column B_j is only evaluated at certain locations since we know in the first place that B_j is going to be sparse. These certain locations are obtained by looking at the sparsity of EE^T where these E and E^T contain delta functions with larger support. For a droptol of 7e-3, the support_nz for EE^T is set to 10 which is approx 3.5 times the original support of 3.
  
    PetscInt :: i, support_nz, nvel, nnnn, mmmm
    Mat :: regTmat, EEtmat, regTmatT
    PetscInt, pointer :: ia(:),ja(:)
    PetscInt :: nrows
    PetscBool :: done
    VecScatter :: ctx_xsd, ctx_ysd
    Vec :: xsd0vec, ysd0vec
    integer, allocatable :: pcnz_int4(:)
  
       !Let's create EE^T matrix only on rank 0-------------
       support_nz = 10!3.5*support       
       nvel = (m+1)*n + m*(n+1)
       if (rank==0) then
           call MatCreate(MPI_COMM_SELF, regTmat, ierr)
           call MatSetType(regTmat, MATSEQAIJ, ierr)
           call MatSetSizes(regTmat, nf_sd, nvel, nf_sd, nvel, ierr)
           call MatSeqAIJSetPreallocation(regTmat, (2*support_nz+1)**2, PETSC_NULL_INTEGER, ierr)
           call MatSetUp(regTmat, ierr)
       end if
       
       !Collect the entire xsdvec and ysdvec on proc 0-------
       call VecScatterCreateToZero(xsdvec, ctx_xsd, xsd0vec, ierr)
       call VecScatterCreateToZero(ysdvec, ctx_ysd, ysd0vec, ierr)
       call VecScatterBegin(ctx_xsd, xsdvec, xsd0vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterBegin(ctx_ysd, ysdvec, ysd0vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterEnd(ctx_xsd, xsdvec, xsd0vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterEnd(ctx_ysd, ysdvec, ysd0vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call VecScatterDestroy(ctx_xsd, ierr)
       call VecScatterDestroy(ctx_ysd, ierr)
       !-----------------------------------------------------
       
       if (rank==0) then
           call VecGetArrayReadF90(xsd0vec, xsd, ierr)
           call VecGetArrayReadF90(ysd0vec, ysd, ierr)      
           call setup_subdomain_sparsity(xsd, ysd, support_nz, regTmat)  
           call VecRestoreArrayReadF90(xsd0vec, xsd, ierr)
           call VecRestoreArrayReadF90(ysd0vec, ysd, ierr) 
       
           call MatMatTransposeMult(regTmat, regTmat, MAT_INITIAL_MATRIX, two, EEtmat, ierr)
           call MatDestroy(regTmat, ierr)
           call MatGetRowIJF90(EEtmat, ione, PETSC_FALSE, PETSC_FALSE, nrows, ia, ja, done, ierr) !get sparsity      
       
           allocate(precomputeserialnz(nf_sd))
           do i=1,nf_sd
               precomputeserialnz(i) = (ia(i+1) - ia(i))  !number of non-zeros
           end do
           nz = 2*maxval(precomputeserialnz)
           allocate(precomputeserialmat(nf_sd),precomputeserialindicesmat(nf_sd),precomputeserialindicesmat2(nf_sd))  !first allocate pointer
       
           do i=1,nf_sd
               allocate(precomputeserialindicesmat(i)%value(2*precomputeserialnz(i)))  
               precomputeserialindicesmat(i)%value(:) = 0

               !non-zero locations
               if (i<=nb_sd) then
                   precomputeserialindicesmat(i)%value(1:precomputeserialnz(i)) = ja(ia(i):ia(i+1)-1)
                   precomputeserialindicesmat(i)%value(precomputeserialnz(i)+1:2*precomputeserialnz(i)) = ja(ia(i):ia(i+1)-1) + nb_sd
               else if (i>nb_sd) then
                   precomputeserialindicesmat(i)%value(1:precomputeserialnz(i)) = ja(ia(i):ia(i+1)-1) - nb_sd
                   precomputeserialindicesmat(i)%value(precomputeserialnz(i)+1:2*precomputeserialnz(i)) = ja(ia(i):ia(i+1)-1)
               end if
           end do
           precomputeserialnz = 2*precomputeserialnz
       
           call MatRestoreRowIJF90(EEtmat, ione, PETSC_FALSE, PETSC_FALSE, nrows, ia, ja, done, ierr)
           call MatDestroy(EEtmat, ierr)
       end if
       
       call VecDestroy(xsd0vec, ierr)
       call VecDestroy(ysd0vec, ierr)
       !-----------------------------------------
       
       !Send the value of nz to all processors
       nz_int4 = nz
       call MPI_Bcast(nz_int4,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
       nz = nz_int4
       
       if (rank .ne. 0) allocate(precomputeserialnz(nf_sd))
       allocate(pcnz_int4(nf_sd))
       pcnz_int4 = precomputeserialnz
       call MPI_Bcast(pcnz_int4, nf_sd, MPI_INT, 0, MPI_COMM_WORLD, ierr)
       precomputeserialnz = pcnz_int4
       deallocate(pcnz_int4)
  
  end subroutine get_sparsity_EEt
            
!================================================================================

  subroutine create_linear_operators
  
      call DMCreateMatrix(das, CTCmat, ierr)
      call MatDuplicate(CTCmat, MAT_SHARE_NONZERO_PATTERN, Imat, ierr)
      
      allocate(Amat(mgridlev), Amat2(mgridlev))
  
  end subroutine create_linear_operators

!================================================================================  

  subroutine setup_ksp
  
    PetscInt :: k
  
    allocate(vort_solver(mgridlev))
    !ksp for 2D domain decomposition where CG is used
    do k=1, mgridlev
        call KSPCreate(MPI_COMM_WORLD, vort_solver(k)%ksp, ierr)
        call KSPSetType(vort_solver(k)%ksp, KSPCG, ierr)
        call KSPGeTPC(vort_solver(k)%ksp, vort_solver(k)%pc, ierr)
        call PCSetType(vort_solver(k)%pc, PCNONE, ierr)
        call KSPSetFromOptions(vort_solver(k)%ksp, ierr)
        call KSPSetInitialGuessNonzero(vort_solver(k)%ksp, PETSC_TRUE, ierr)
        call KSPSetNormType(vort_solver(k)%ksp, KSP_NORM_UNPRECONDITIONED, ierr)
    end do
    
    !ksp for 2D domain decomposition where CG is used
    call KSPCreate(MPI_COMM_WORLD, ksp_sf, ierr)
    call KSPSetType(ksp_sf, KSPCG, ierr)
    call KSPGeTPC(ksp_sf, pc_sf, ierr)
    call PCSetType(pc_sf, PCNONE, ierr)
    call KSPSetFromOptions(ksp_sf, ierr)
    call KSPSetInitialGuessNonzero(ksp_sf, PETSC_TRUE, ierr)
    call KSPSetNormType(ksp_sf, KSP_NORM_UNPRECONDITIONED, ierr)
  
    !Ksp for Cholesky decomposition
    if (num_stat) then
      call KSPCreate(MPI_COMM_WORLD, ksp, ierr)
      call KSPSetType(ksp, KSPPREONLY, ierr)
      call KSPGeTPC(ksp, pc, ierr)
      call PCSetType(pc, PCCHOLESKY, ierr)
      call PCFactorSetMatSolverTYPE(pc, MATSOLVERELEMENTAL, ierr)
      call KSPSetFromOptions(ksp, ierr)
    end if
    
    !Ksp for FSI problems
    if (.not. num_stat) then
    if (.not. sub_domain) then
      call KSPCreate(MPI_COMM_WORLD, ksp_fsi, ierr)
      call KSPSetType(ksp_fsi, KSPGMRES, ierr)
      call KSPGeTPC(ksp_fsi, pc_fsi, ierr)
      call PCSetType(pc_fsi, PCNONE, ierr)
      call KSPSetFromOptions(ksp_fsi, ierr)
      call KSPSetInitialGuessNonzero(ksp_fsi, PETSC_TRUE, ierr)
      call KSPSetNormType(ksp_fsi, KSP_NORM_UNPRECONDITIONED, ierr)
      call KSPSetTolerances(ksp_fsi, rtol, atol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER,ierr)
    else if (sub_domain) then
      
      call KSPCreate(MPI_COMM_WORLD, ksp_fsi_sd, ierr)
      call KSPSetType(ksp_fsi_sd, KSPGMRES, ierr)
      call KSPGeTPC(ksp_fsi_sd, pc_fsi_sd, ierr)
      call PCSetType(pc_fsi_sd, PCNONE, ierr)
      call KSPSetFromOptions(ksp_fsi_sd, ierr)
      call KSPSetInitialGuessNonzero(ksp_fsi_sd, PETSC_TRUE, ierr)
      call KSPSetNormType(ksp_fsi_sd, KSP_NORM_UNPRECONDITIONED, ierr)
      call KSPSetTolerances(ksp_fsi_sd, rtol, atol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER,ierr)
    end if  
    end if
    
    
    
  end subroutine setup_ksp  
  
!================================================================================  

  subroutine destroy_ksp
  
    PetscInt :: k  
  
    do k=1, mgridlev
        call KSPDestroy(vort_solver(k)%ksp, ierr)
    end do
    call KSPDestroy(ksp_sf, ierr)
    deallocate(vort_solver)
    
     if (num_stat) then
        call KSPDestroy(ksp, ierr)
     end if
    
     if (.not. num_stat) then
     if (.not. sub_domain) then
        call KSPDestroy(ksp_fsi, ierr)
     else if (sub_domain) then   
        call KSPDestroy(ksp_fsi_sd, ierr)
     end if
     end if
  
  end subroutine destroy_ksp

!================================================================================

  subroutine setup_body_dist  
    
    call VecCreateSeq(PETSC_COMM_SELF, xbm, fivec, ierr)
    call VecCreateMPI(comm_body, PETSC_DECIDE, xbm_body, fbody_vec, ierr)
    call VecScatterCreate(fvec, xbm_glob_is, fivec, xbm_is, ctx_force, ierr) !Scatter context for scatter level 1
    call VecScatterCreate(fivec, xbm_is, fbody_vec, xbm_is, ctx_body, ierr)  !Scatter context for scatter level 2
    call VecScatterCopy(ctx_force, ctx_force2, ierr)
    call VecScatterCopy(ctx_body, ctx_body2, ierr)
    
    call ISDestroy(xbm_is, ierr)
    call ISDestroy(xbm_glob_is, ierr)
    
    !temporary force vectors
    call VecDuplicate(fivec, xivec, ierr) !btimes: scatter level 1
    call VecDuplicate(fbody_vec, xbody_vec, ierr)  !btimes: Scatter level 2
    call VecDuplicate(fbody_vec, ybody_vec, ierr) !btimes: Result of multiplication on scatter level 2
    call VecDuplicate(fivec, yivec, ierr) !btimes: Result on scatter level 1
    call VecDuplicate(fbody_vec, y1body_vec, ierr) !compute_dx: Result of multiplication on scatter level 2
    call VecDuplicate(fivec, y1ivec, ierr) !compute_dx: Result on scatter level 1
    call VecDuplicate(fbody_vec, y2body_vec, ierr) !compute_dx: Result of multiplication on scatter level 2
    call VecDuplicate(fivec, y2ivec, ierr) !compute_dx: Result on scatter level 1
    call VecDuplicate(fvec, y1vec, ierr) !in ab_times
    call VecDuplicate(fvec, y2vec, ierr) !in ab_times
    
    !temporary force vectors using in redistribute
    call VecDuplicate(velocityxvec, wghtxvec, ierr) 
    call VecDuplicate(velocityyvec, wghtyvec, ierr) 
    call VecDuplicate(velocityxvec, frcxvec, ierr) 
    call VecDuplicate(velocityyvec, frcyvec, ierr) 
    call VecDuplicate(velocityxvec, finterx_vec, ierr) 
    call VecDuplicate(velocityyvec, fintery_vec, ierr) 
    
  end subroutine setup_body_dist  
  
!================================================================================  

  subroutine create_struct_matrix
  !Setup structural matrices
  
    PetscInt :: i, i_bdy, iter_mt, iter_md
    
    if (mt_body>0) allocate( sol_Jmat(mt_body), Mtmat(mt_body), Jmat(mt_body), Ftint(mt_body), Ctmat(mt_body) )
    if (md_body>0) allocate( sol_Kmat(md_body), Mdmat(md_body), Kmat(md_body) )
    
    iter_mt=0
    iter_md=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body matrix
        iter_mt=iter_mt+1
        call MatCreate(comm_body, sol_Jmat(iter_mt), ierr)
        call MatSetType(sol_Jmat(iter_mt), MATDENSE, ierr)
        call MatSetSizes(sol_Jmat(iter_mt), PETSC_DECIDE, PETSC_DECIDE, 2*nb_ind(i_bdy), 2*nb_ind(i_bdy), ierr)
        call MatSetFromOptions(sol_Jmat(iter_mt), ierr)
        call MatSetUp(sol_Jmat(iter_mt), ierr)
        
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body matrix
        iter_md=iter_md+1
        !Assign sizes for Kmat, Mdmat and sol_Kmat
      end if
      
    end do

  end subroutine create_struct_matrix
  
!================================================================================  

  subroutine destroy_struct_matrix
  
    PetscInt :: i, i_bdy, iter_mt, iter_md
    
    if (mt_body>0)  then
       call MatDestroyMatrices(mt_body, sol_Jmat, ierr)
       deallocate(Mtmat, Jmat, Ftint, Ctmat)
    end if
    
    if (md_body>0)  then
       call MatDestroyMatrices(md_body, sol_Kmat, ierr)
       call MatDestroyMatrices(md_body, Mdmat, ierr)
       call MatDestroyMatrices(md_body, Kmat, ierr)
    end if
  
  end subroutine destroy_struct_matrix

!================================================================================

  subroutine write_variables(it)
  
    PetscInt :: it, k
    character(2) :: charrank
    character(7) :: charit
    PetscScalar, pointer :: readvar(:)
    PetscScalar :: omega_write(xsm*ysm), s_write(xsm*ysm), f_write(xfm), xb_write(xfm), yb_write(xfm), velx_write(xum*yum), vely_write(xvm*yvm)
    
    if (rank==0) write(*,*) 'writing variables at it=',it
    WRITE(charit,"(I7.7)") it
    WRITE(charrank,"(I2.2)") rank
    OPEN(unit=(rank+1)*100,file="output/data"//charrank//"_"//charit//".var",form="unformatted",status="unknown")
    
    WRITE((rank+1)*100) m,n,mgridlev,nb
    WRITE((rank+1)*100) re,dt,len,offsetx,offsety
    WRITE((rank+1)*100) xsi, xse, xsm, ysi, yse, ysm
    WRITE((rank+1)*100) xui, xue, xum, yui, yue, yum
    WRITE((rank+1)*100) xvi, xve, xvm, yvi, yve, yvm
    
    do k=1,mgridlev
      call VecGetarrayReadF90(omegavec(k), readvar, ierr)
      omega_write=readvar
      WRITE((rank+1)*100) omega_write
      call VecRestorearrayReadF90(omegavec(k), readvar, ierr)
      
      call VecGetarrayReadF90(svec(k), readvar, ierr)
      s_write=readvar
      WRITE((rank+1)*100) s_write
      call VecRestorearrayReadF90(svec(k), readvar, ierr)
      
      call VecGetarrayReadF90(qxvec(k), readvar, ierr)
      velx_write=readvar
      WRITE((rank+1)*100) velx_write
      call VecRestorearrayReadF90(qxvec(k), readvar, ierr)
      
      call VecGetarrayReadF90(q0pxvec(k), readvar, ierr)
      velx_write=readvar
      WRITE((rank+1)*100) velx_write
      call VecRestorearrayReadF90(q0pxvec(k), readvar, ierr)
      
      call VecGetarrayReadF90(q0rxvec(k), readvar, ierr)
      velx_write=readvar
      WRITE((rank+1)*100) velx_write
      call VecRestorearrayReadF90(q0rxvec(k), readvar, ierr)
      
      call VecGetarrayReadF90(qyvec(k), readvar, ierr)
      vely_write=readvar
      WRITE((rank+1)*100) vely_write
      call VecRestorearrayReadF90(qyvec(k), readvar, ierr)
      
      call VecGetarrayReadF90(q0pyvec(k), readvar, ierr)
      vely_write=readvar
      WRITE((rank+1)*100) vely_write
      call VecRestorearrayReadF90(q0pyvec(k), readvar, ierr)
      
      call VecGetarrayReadF90(q0ryvec(k), readvar, ierr)
      vely_write=readvar
      WRITE((rank+1)*100) vely_write
      call VecRestorearrayReadF90(q0ryvec(k), readvar, ierr)
    end do
    
    WRITE((rank+1)*100) rhs_old
    
    WRITE((rank+1)*100) rot_angle, rox, roy
    
    WRITE((rank+1)*100) xfm, xfi, xfe
    call VecGetarrayReadF90(fvec, readvar, ierr)
    f_write = readvar
    WRITE((rank+1)*100) f_write
    call VecRestorearrayReadF90(fvec, readvar, ierr)
    
    call VecGetarrayReadF90(xbvec, readvar, ierr)
    xb_write = readvar
    WRITE((rank+1)*100) xb_write
    call VecRestorearrayReadF90(xbvec, readvar, ierr)
    
    call VecGetarrayReadF90(ybvec, readvar, ierr)
    yb_write = readvar
    WRITE((rank+1)*100) yb_write
    call VecRestorearrayReadF90(ybvec, readvar, ierr)
    
    WRITE((rank+1)*100) mt_body
    if (mt_body>0) then 
      WRITE((rank+1)*100) theta, thetad, thetadd
    end if
    
    call VecGetarrayReadF90(f_rdstvec, readvar, ierr)
    f_write = readvar
    WRITE((rank+1)*100) f_write
    call VecRestorearrayReadF90(f_rdstvec, readvar, ierr)
    
    
    CLOSE((rank+1)*100)
   
  end subroutine 
  
!================================================================================

  subroutine read_variables(it)
  
    PetscInt :: it, k
    character(2) :: charrank
    character(7) :: charit
    PetscScalar, pointer :: readvar(:)
    PetscScalar :: omega_write(xsm*ysm), s_write(xsm*ysm), f_write(xfm), xb_write(xfm), yb_write(xfm), velx_write(xum*yum), vely_write(xvm*yvm)
    PetscScalar :: dt_old
    
    if (rank==0) write(*,*) 'reading variables at it=',it
    WRITE(charit,"(I7.7)") it
    WRITE(charrank,"(I2.2)") rank
    OPEN(unit=(rank+1)*100,file="output/data"//charrank//"_"//charit//".var",form="unformatted",status="unknown")
    
    READ((rank+1)*100) m,n,mgridlev,nb
    READ((rank+1)*100) re,dt_old,len,offsetx,offsety
    READ((rank+1)*100) xsi, xse, xsm, ysi, yse, ysm
    READ((rank+1)*100) xui, xue, xum, yui, yue, yum
    READ((rank+1)*100) xvi, xve, xvm, yvi, yve, yvm
    
    do k=1,mgridlev
      call VecGetArrayF90(omegavec(k), readvar, ierr)
      omega_write=readvar
      READ((rank+1)*100) omega_write
      readvar = omega_write
      call VecRestoreArrayF90(omegavec(k), readvar, ierr)
      
      call VecGetArrayF90(svec(k), readvar, ierr)
      s_write=readvar
      READ((rank+1)*100) s_write
      readvar = s_write
      call VecRestoreArrayF90(svec(k), readvar, ierr)
      
      call VecGetArrayF90(qxvec(k), readvar, ierr)
      velx_write=readvar
      READ((rank+1)*100) velx_write
      readvar = velx_write
      call VecRestoreArrayF90(qxvec(k), readvar, ierr)
      
      call VecGetArrayF90(q0pxvec(k), readvar, ierr)
      velx_write=readvar
      READ((rank+1)*100) velx_write
      readvar = velx_write
      call VecRestoreArrayF90(q0pxvec(k), readvar, ierr)
      
      call VecGetArrayF90(q0rxvec(k), readvar, ierr)
      velx_write=readvar
      READ((rank+1)*100) velx_write
      readvar = velx_write
      call VecRestoreArrayF90(q0rxvec(k), readvar, ierr)
      
      call VecGetArrayF90(qyvec(k), readvar, ierr)
      vely_write=readvar
      READ((rank+1)*100) vely_write
      readvar = vely_write
      call VecRestoreArrayF90(qyvec(k), readvar, ierr)
      
      call VecGetArrayF90(q0pyvec(k), readvar, ierr)
      vely_write=readvar
      READ((rank+1)*100) vely_write
      readvar = vely_write
      call VecRestoreArrayF90(q0pyvec(k), readvar, ierr)
      
      call VecGetArrayF90(q0ryvec(k), readvar, ierr)
      vely_write=readvar
      READ((rank+1)*100) vely_write
      readvar = vely_write
      call VecRestoreArrayF90(q0ryvec(k), readvar, ierr)
    end do
    
    READ((rank+1)*100) rhs_old
    
    READ((rank+1)*100) rot_angle, rox, roy
    
    READ((rank+1)*100) xfm, xfi, xfe
    call VecGetArrayF90(fvec, readvar, ierr)
    f_write = readvar
    READ((rank+1)*100) f_write
    readvar = f_write
    call VecRestoreArrayF90(fvec, readvar, ierr)
    
    call VecGetArrayF90(xbvec, readvar, ierr)
    xb_write = readvar
    READ((rank+1)*100) xb_write
    readvar = xb_write
    call VecRestoreArrayF90(xbvec, readvar, ierr)
    
    call VecGetArrayF90(ybvec, readvar, ierr)
    yb_write = readvar
    READ((rank+1)*100) yb_write
    readvar = yb_write
    call VecRestoreArrayF90(ybvec, readvar, ierr)
    
    READ((rank+1)*100) mt_body
    if (mt_body>0) then 
      READ((rank+1)*100) theta, thetad, thetadd
    end if
    
    
    CLOSE((rank+1)*100)
    
    do k=1,mgridlev
      call VecWAXPY(q0xvec(k), one, q0pxvec(k), q0rxvec(k), ierr)
      call VecWAXPY(q0yvec(k), one, q0pyvec(k), q0ryvec(k), ierr)
    end do
   
  end subroutine 
  
!================================================================================

  subroutine write_proc
  
    if (rank==0) then
      OPEN(unit=(rank+1)*101,file="output/nproc.dat")
      WRITE((rank+1)*101,*) nproc
      CLOSE((rank+1)*101)
    end if
    
  end subroutine write_proc

!================================================================================

  subroutine destroy_variables
  
     USE myfft
     
     PetscInt :: k, i
     
     call destroy_fft
     
     call VecDestroy(streamvec, ierr)
     call VecDestroy(streamvec_coarsify, ierr)
     call VecDestroyVecsF90(mgridlev, omegavec, ierr)
     call VecDestroyVecsF90(mgridlev, svec, ierr)
     call DMRestoreLocalVector(das, streamvec_local, ierr)
     
     call VecDestroy(streambcvec_local, ierr)
     call VecDestroy(streambcvec_global, ierr)
     
     call VecDestroy(velocityxvec, ierr)
     call VecDestroyVecsF90(mgridlev, qxvec, ierr)
     call VecDestroyVecsF90(mgridlev, q0xvec, ierr)
     call VecDestroyVecsF90(mgridlev, q0rxvec, ierr)
     call VecDestroyVecsF90(mgridlev, q0pxvec, ierr)

     call VecDestroy(velocityyvec, ierr)
     call VecDestroyVecsF90(mgridlev, qyvec, ierr)
     call VecDestroyVecsF90(mgridlev, q0yvec, ierr)
     call VecDestroyVecsF90(mgridlev, q0ryvec, ierr)
     call VecDestroyVecsF90(mgridlev, q0pyvec, ierr)
         
     call VecDestroy(velocityxvec_interface, ierr)
     call VecDestroy(velocityyvec_interface, ierr)
     
     call VecDestroy(f_rdstvec, ierr)
     call VecDestroy(fivec, ierr)
     call VecDestroy(fbody_vec, ierr)
     call VecDestroy(xivec, ierr) !btimes
     call VecDestroy(xbody_vec, ierr) !btimes
     call VecDestroy(yivec, ierr) !btimes
     call VecDestroy(ybody_vec, ierr) !btimes
     call VecDestroy(y1ivec, ierr) !compute_dx
     call VecDestroy(y1body_vec, ierr) !compute_dx
     call VecDestroy(y2ivec, ierr) !compute_dx
     call VecDestroy(y2body_vec, ierr) !compute_dx
     call VecDestroy(wghtxvec, ierr) !redistribute
     call VecDestroy(wghtyvec, ierr) !redistribute
     call VecDestroy(frcxvec, ierr) !redistribute
     call VecDestroy(frcyvec, ierr) !redistribute
     call VecDestroy(finterx_vec, ierr) !redistribute
     call VecDestroy(fintery_vec, ierr) !redistribute
     call VecDestroy(y1vec, ierr)
     call VecDestroy(y2vec, ierr)
     
     deallocate(rhs_old)
     
     call VecDestroy(rhsfluid_vec, ierr) !rhs term of fluid equation
     call VecDestroy(rhsstream_vec, ierr) !rhs term of streamfunction equation
     call VecDestroyVecsF90(mgridlev, streambcvecs_global, ierr) !boundary
     call VecDestroyVecsF90(mgridlev, lastbcvec_local, ierr) !boundary
     call VecDestroyVecsF90(mgridlev, omegabcvec_local, ierr) !boundary
     call VecDestroyVecsF90(mgridlev, omegavec_local, ierr) !vorticity
     call VecDestroy(qhxvec, ierr) 
     call VecDestroy(qhyvec, ierr) 
     call DMRestoreLocalVector(dau, qxvec_local, ierr) !x-velocity     
     call DMRestoreLocalVector(dav, qyvec_local, ierr) !x-velocity     
     call VecDestroy(nlxvec, ierr) !x-nonlinear term
     call VecDestroy(nlyvec, ierr) !y-nonlinear term
     call DMRestoreLocalVector(dau, nlxvec_local, ierr) !x-nonlinear_local term with ghost
     call DMRestoreLocalVector(dav, nlyvec_local, ierr) !y-nonlinear_local term with ghost
     call VecDestroy(rhsfvec, ierr)
     if (sub_domain) call VecDestroy(rhsfsdvec, ierr)
     call VecDestroy(velxvec_interface, ierr) !for regT
     call VecDestroy(velyvec_interface, ierr) !for regT
     if (sub_domain) then
         call VecDestroy(velxvec_interface_sd, ierr) !for regT
         call VecDestroy(velyvec_interface_sd, ierr) !for regT
         call VecDestroy(velxvec_interface_sd_sparse, ierr)
         call VecDestroy(velyvec_interface_sd_sparse, ierr)
     end if
     call VecDestroy(hxvec, ierr) !for reg
     call VecDestroy(hyvec, ierr) !for reg
     call VecDestroy(vortvec, ierr) !for rot after reg
     call VecDestroy(rhsvortvec, ierr) !rhs of last equation
     call DMRestoreLocalVector(dau, hxvec_local, ierr) !for rot after reg
     call DMRestoreLocalVector(dav, hyvec_local, ierr) !for rot after reg
     call VecDestroyVecsF90(mgridlev, vortvecs, ierr) !for a_times
     call VecDestroyVecsF90(mgridlev, snvecs, ierr) !for a_times
     call VecDestroyVecsF90(mgridlev, velxvecs, ierr) !for a_times
    call VecDestroyVecsF90(mgridlev, velyvecs, ierr) !for a_times
    call DMRestoreLocalVector(das, vortvec_local1, ierr) !for coarsify in vort2flux
    call DMRestoreLocalVector(das, svec_local, ierr) !for coarsify in vort2flux
     
      
     if (num_stat) call MatDestroy(cholmat, ierr)
     if (sub_domain) then
       
         do i = 1,nf_sd
             if (i.ge.xsdi .and. i.le.xsde) then
                 deallocate(precomputemat(i)%value,precomputeindicesmat(i)%value)
             end if
         end do
         deallocate(precomputemat, precomputeindicesmat, precomputenz)
         deallocate(indicesserialmat)
         call MatDestroy(lumat, ierr)
     end if
     call MatDestroy(A, ierr)
     call MatDestroy(B, ierr)
     call MatDestroy(AB, ierr)
     
     call VecScatterDestroy(ctx_velx, ierr)
     call VecScatterDestroy(ctx_vely, ierr)
     if (sub_domain) then
         call VecScatterDestroy(ctx_velx_sd_sparse, ierr)
         call VecScatterDestroy(ctx_vely_sd_sparse, ierr)
     end if
     call VecScatterDestroy(ctx_bc, ierr)
     call VecScatterDestroy(ctx_force, ierr)
     call VecScatterDestroy(ctx_body, ierr)
     call VecScatterDestroy(ctx_force2, ierr)
     call VecScatterDestroy(ctx_body2, ierr)
     
     call destroy_ksp
  
     call destroy_struct_matrix
     if (mt_body>0) deallocate(theta, thetad, thetadd, dtheta, theta0, thetad0, thetadd0)
  
  end subroutine destroy_variables

!================================================================================

END MODULE variables
