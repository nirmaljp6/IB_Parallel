MODULE grid

  USE petsc
  use, intrinsic :: iso_c_binding  
#include "petsc/finclude/petsc.h"
  USE parameters
  use comm
  IMPLICIT NONE
    
  
  PetscScalar :: delta ! near field grid spacing
  PetscScalar :: support = 3.D0  ! support for smearing delta functions
  PetscInt :: interp_points = 4 !Number of interpolation points to choose for the sub-domain approach; denoted by n_p in the paper
  PetscScalar :: droptol = 7.0e-3 !drop tolerance for the sub-domain approach
  
  ! Numbers of cells, edges, etc.
  PetscInt :: ns, nb, nf, nr, nt, nd, support2
  DM :: das, das_coarsify, dau, dav
  PetscInt, dimension(:), allocatable :: nr_ind, nt_ind, nd_ind, nb_ind, nbb_ind
  
  PetscInt :: xsi, ysi, xsm, ysm, xsgi, ysgi, xsgm, ysgm, xse, yse, xsge, ysge
  PetscInt :: xui, yui, xum, yum, xugi, yugi, xugm, yugm, xue, yue, xuge, yuge
  PetscInt :: xvi, yvi, xvm, yvm, xvgi, yvgi, xvgm, yvgm, xve, yve, xvge, yvge
  PetscInt :: xpi, xpe, ypi, ype, xpm, ypm, xpgi, xpge, ypgi, ypge, xpgm, ypgm
  PetscInt :: xfi, xfe, xfm
  PetscInt :: xci, yci, tcoarse
  PetscInt :: tbc  !size of some arrays in get_bc 
  PetscInt, allocatable :: bcix(:), cbcix(:)
  IS :: bcis, interfacex_is, interfacey_is
  PetscInt :: nbc, left, right, bottom, top
  
  PetscScalar, pointer :: xb(:), yb(:)
  PetscInt, dimension(:,:), allocatable :: codeb
  
  ! index of each ib points relative to grid
    PetscInt, DIMENSION(:), ALLOCATABLE :: indexx, indexy, reg_ix, reg_iy
    IS :: reg_ixs, reg_iys
    ! arrays for smearing coefficients
    PetscScalar, DIMENSION(:,:), ALLOCATABLE :: weight
    PetscInt :: interfacex, interfacey
  
  ! a special type for bodies
  TYPE body
     PetscInt :: npts  !number of points on the body
     PetscInt :: fam   !family of the body: 1-rigid, 2-torsional, 3-deformable
     PetscInt :: bn    !assigned body number starting from 1
     PetscInt :: ipos  ! position in overall array
     REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: x,y
     PetscScalar :: gamma0  !Equilibrum angle w.r.t x-axis
     PetscScalar :: alpha0  !Equilibrum angle w.r.t previous body: always a constant no matter what
     PetscScalar :: itheta  !inertia of torsional body
     PetscScalar :: ktheta  !stiffness of torsional body
     PetscScalar :: ctheta  !damping coefficient of torsional body
     PetscScalar :: deltas  !body discretization
     PetscScalar :: xbf0, ybf0  !location at which the torsional body is hinged
  END TYPE body
  
  PetscInt, PARAMETER :: maxbodies = 999 ! a large number
  TYPE(body), DIMENSION(maxbodies) :: bdy
  
  Vec :: xbvec, ybvec, fvec
  integer :: body_color
  PetscInt, allocatable :: body_rank(:,:), rank_nb(:)  !body_rank bookkeeps which body is in which rank; rank_nb bookkeeps how many body points are there in a particular rank
  PetscInt :: xbm !It is like xfm, but the number of body unknowns to be scattered from fvec to sequential vectors
  PetscInt :: xbm_body  ! It is the value of xbm, but scattered to other processors ob comm_body for creating an MPI vector on comm_body
  MPI_Comm :: comm_body
  IS :: xbm_is, xbm_glob_is
  PetscInt :: mrb !this is the number of bodies in a given rank for btimes
  PetscInt :: mrb_body  !this is the number of bodies in a given rank for btimes which is same as mrb. However, mrb_body is a scattered mrb to other procs for enabling matrix vector multiplication in btimes.
  PetscInt, allocatable :: rank_body(:)   !This is local to each rank, and it contains information of which bodies the running processor has to handle in btimes
  integer :: nproc_comm_body_int4, rank_comm_body_int4
  PetscInt :: nproc_comm_body, rank_comm_body
  PetscInt :: mt_body, md_body !number of bodies of torsional and deformable families to be handled by each proc 
  !PetscInt :: xbi, xbe !starting and ending indices of fbody_vec
  
  PetscInt :: nb_sd, nf_sd, nr_sd !number of grid points in subdomina. nr_sd is the number of rigid body points included in the subdomain
  PetscScalar, dimension(:), allocatable :: sd_x, sd_y ! x and y coordinates of the sub domain
  PetscScalar, pointer :: xsd(:), ysd(:)
  PetscInt, dimension(:,:), allocatable :: codesd
  Vec :: xsdvec, ysdvec
  PetscInt :: xsdi, xsde, xsdm, interfacex_sd, interfacey_sd
  PetscInt, DIMENSION(:), ALLOCATABLE :: indexx_sd, indexy_sd, reg_ix_sd, reg_iy_sd
  PetscScalar, DIMENSION(:,:), ALLOCATABLE :: weight_sd
  IS :: interfacex_sd_is, interfacey_sd_is, reg_ixs_sd, reg_iys_sd
  
  PetscScalar, dimension(:,:), allocatable :: interp_weights, indicesserialmat
  Mat :: indicesmat
  PetscInt, dimension(:,:), allocatable :: interp_indices
  PetscInt, dimension(:), allocatable :: interp_indices_global, interp_indices_local
  PetscInt :: interfacex_sd_sparse, interfacey_sd_sparse
  PetscInt, DIMENSION(:), ALLOCATABLE :: reg_ix_sd_sparse, reg_iy_sd_sparse
  IS :: reg_ixs_sd_sparse, reg_iys_sd_sparse, interfacex_sd_is_sparse, interfacey_sd_is_sparse
  
  PetscInt :: nz, half_width, nz_indicesmat
  integer :: nz_int4
  
  PetscScalar, DIMENSION(:,:), ALLOCATABLE :: weight_sd_sparsity
  PetscInt :: interface_sd_sparsity
  PetscInt, DIMENSION(:), ALLOCATABLE :: reg_i_sd_sparsity
  
  Vec :: interp_indices_vec, interp_indices_global_vec, interp_weights_vec, interp_weights_global_vec
  VecScatter :: ctx_indices, ctx_weights
  PetscInt, allocatable :: interp_indices_global_sparsemat(:,:)
  PetscScalar, allocatable :: interp_weights_global_sparsemat(:,:)
  Vec :: intermediatevec, intermediatelocalvec
  VecScatter :: ctx_intermediate
  PetscInt :: interp_indices_local_size, interp_indices_global_nproc
  IS :: interp_indices_local_is
  
  PetscInt :: izero, ione, itwo, ithree
  PetscScalar :: one, zero, two
  
CONTAINS

!===================================================================================
  
  SUBROUTINE setup 

    use myfft
    PetscInt :: ly(nproc), lsy(nproc), lsx(1), luy(nproc), lux(1), lvy(nproc), lvx(1)
    
    izero=0
    ione=1
    itwo=2
    ithree=3
    one=1.d0
    zero=0.d0
    two=2.d0

    ! firt order of business is to setup the grid
    delta = len/REAL(m) 
 
    ! a few things to setup for multigrid laplace solution
    ! check that grid is divisible by 4
    IF ((MOD(m,4).ne.0).or.(MOD(n,4).ne.0)) THEN
       IF (rank==0) then
         STOP 'grid must be divisible by 4'
       END IF
    END IF
    
    ! indexing for streamfunction
    ns = (m-1)*(n-1)
    
    !Determine how many rows 1 proc should  take. In case of FFTW, we need to use FFTW routines to determine the number ly
    call fft_grid(ly)
    
    !Perform flow domain partitioning
    if (dimen==2) then
      !Domain decomposition for vorticity and streamfunction
      lsx(1)=m-1
      lsy=ly
      call DMDAcreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, m-1, n-1, ione, nproc, ione, ione, lsx, lsy, das, ierr)
      
      call DMDAcreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, m-1, n-1, ione, nproc, ione, ione, lsx, lsy, das_coarsify, ierr)
      
      !Domain decomposition for velocity-flux in x-direction
      lux(1)=m+1
      luy=lsy
      luy(nproc)=luy(nproc)+1
      call DMDAcreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, m+1, n, ione, nproc, ione, ione, lux, luy, dau, ierr)
      
      !Domain decomposition for velocity-flux in y-direction
      lvx(1)=m
      lvy=lsy
      lvy(1)=lvy(1)+1
      lvy(nproc)=lvy(nproc)+1
      call DMDAcreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, m, n+1, ione, nproc, ione, ione, lvx, lvy, dav, ierr)

      !for a 2D domain decomposition, BOX stencil needs to be used above so that nonlinear and coarsify subroutines undergo appropriate communications

    end if
    
    call DMSetFromOptions(das, ierr)
    call DMSetUp(das, ierr)
    
    call DMSetFromOptions(das_coarsify, ierr)
    call DMSetUp(das_coarsify, ierr)

    call DMSetFromOptions(dau, ierr)
    call DMSetUp(dau, ierr)

    call DMSetFromOptions(dav, ierr)
    call DMSetUp(dav, ierr)
    
    !Get indices of the corners of the various grid
    call DMDAGetCorners(das, xsi, ysi, PETSC_NULL_INTEGER, xsm, ysm, PETSC_NULL_INTEGER, ierr)
    call DMDAGetGhostCorners(das, xsgi, ysgi, PETSC_NULL_INTEGER, xsgm, ysgm, PETSC_NULL_INTEGER, ierr)
    
    call DMDAGetCorners(dau, xui, yui, PETSC_NULL_INTEGER, xum, yum, PETSC_NULL_INTEGER, ierr)
    call DMDAGetGhostCorners(dau, xugi, yugi, PETSC_NULL_INTEGER, xugm, yugm, PETSC_NULL_INTEGER, ierr)
    
    call DMDAGetCorners(dav, xvi, yvi, PETSC_NULL_INTEGER, xvm, yvm, PETSC_NULL_INTEGER, ierr)
    call DMDAGetGhostCorners(dav, xvgi, yvgi, PETSC_NULL_INTEGER, xvgm, yvgm, PETSC_NULL_INTEGER, ierr)
    
    !Converting the indices from C to Fortran
    xsi=xsi+2
    ysi=ysi+2
    xse=xsi+xsm-1
    yse=ysi+ysm-1
    xsgi=xsgi+2
    ysgi=ysgi+2
    xsge=xsgi+xsgm-1
    ysge=ysgi+ysgm-1
    
    xui=xui+1
    yui=yui+1
    xue=xui+xum-1
    yue=yui+yum-1
    xugi=xugi+1
    yugi=yugi+1
    xuge=xugi+xugm-1
    yuge=yugi+yugm-1
    
    xvi=xvi+1
    yvi=yvi+1
    xve=xvi+xvm-1
    yve=yvi+yvm-1
    xvgi=xvgi+1
    yvgi=yvgi+1
    xvge=xvgi+xvgm-1
    yvge=yvgi+yvgm-1
    
    !For multigrid boundary conditions
    nbc=2*(m+1)+2*(n+1)
    left = 0
    right = n+1
    bottom = 2*(n+1)
    top = 2*(n+1) + m+1
    
    call setup_geometry
    call setup_coarsify
    call setup_getbc
    
  END SUBROUTINE setup
  
  
!===================================================================================  

  SUBROUTINE setup_geometry

    LOGICAL :: readinput  
    CHARACTER(3) :: file_num
    PetscInt :: i, nextr, nextt, nextd, nextb
    PetscInt :: k, j

    !Look for bodies in input directory
    readinput = .TRUE.
    m_body=0
    mr = 0
    mt = 0
    md = 0
    DO WHILE (readinput)
       WRITE(file_num,"(I3.3)") m_body+1
       INQUIRE(file="input/body."//file_num//".inp",exist=readinput)
       IF (readinput) THEN
          m_body=m_body+1
          OPEN(unit=8,file="input/body."//file_num//".inp",form='formatted',status='old')
          READ(8,*) bdy(m_body)%npts
          READ(8,*) bdy(m_body)%fam
          if (bdy(m_body)%fam==2) then
            READ(8,*) bdy(m_body)%gamma0; bdy(m_body)%gamma0=bdy(m_body)%gamma0*pi/180  !converting to radians
            READ(8,*) bdy(m_body)%alpha0; bdy(m_body)%alpha0=bdy(m_body)%alpha0*pi/180  !converting to radians
            READ(8,*) bdy(m_body)%itheta
            READ(8,*) bdy(m_body)%ktheta
            READ(8,*) bdy(m_body)%ctheta
          end if
          ALLOCATE( bdy(m_body)%x(bdy(m_body)%npts), bdy(m_body)%y(bdy(m_body)%npts) )
          DO i=1,bdy(m_body)%npts
             READ(8,*) bdy(m_body)%x(i), bdy(m_body)%y(i)
          END DO
          bdy(m_body)%deltas=sqrt( (bdy(m_body)%x(2) - bdy(m_body)%x(1))**2 + (bdy(m_body)%y(2) - bdy(m_body)%y(1))**2 )   !Assuming uniform grid spacing on the body
          
          if (bdy(m_body)%fam==1) then
            mr = mr+1
            bdy(m_body)%bn=m_body
          else if (bdy(m_body)%fam==2) then
            mt = mt+1
            bdy(m_body)%bn=m_body
            bdy(m_body)%xbf0=bdy(m_body)%x(1)
            bdy(m_body)%ybf0=bdy(m_body)%y(1)
          else if (bdy(m_body)%fam==3) then
            md = md+1
            bdy(m_body)%bn=m_body
          end if
       
          CLOSE(8)
       END IF 
    END DO
    
    if (rank==0) then
      write(*,*) 'read all bodies and actuators'
    end if
    !Accumulate all bodies into global vector xb
    !Below notations are: "n"-total number of body points, "b"-body, "r"-rigid, "t"-torsional, "d"-deformable
    !nb_ind is total points for a given body, nbb_ind is cumulatively added nb_ind 
    nb = 0
    nr = 0
    nt = 0
    nd = 0
    nextb=0
    allocate(nr_ind(mr), nt_ind(mt), nd_ind(md), nb_ind(m_body), nbb_ind(0:m_body))
    nbb_ind(0) = 0
    DO i=1,m_body
       if (rank==0) then
         write(*,*) 'body no.',i,'has',bdy(i)%npts,'points.  Family?',bdy(i)%fam
       end if
       nb = nb + bdy(i)%npts
       nextb=nextb+1
       nb_ind(nextb)=bdy(i)%npts
       nbb_ind(nextb)=nbb_ind(nextb-1)+2*bdy(i)%npts
       if (bdy(i)%fam==1) then
         nr = nr+bdy(i)%npts
       else if (bdy(i)%fam==2) then
         nt = nt+bdy(i)%npts
       else if (bdy(i)%fam==3) then
         nd = nd+bdy(i)%npts
       end if
    END DO
    
    !Read the subdomain grid points
    if (sub_domain) then
        nb_sd = 0
        if (.not. sub_domain_full_rectangular)  nb_sd = nr !assumes that all rigid bodies are motionless
        nr_sd = nb_sd
        nb_sd = nb_sd + nint(((sdxe-sdxi)/delta_sd+1)*((sdye-sdyi)/delta_sd+1))  !adding number of points from the subdomain
        ALLOCATE( sd_x(nb_sd), sd_y(nb_sd)) !sd_x and sd_y contains the coordinates of all sub-domain points
        
        nextb = 0
        if (.not. sub_domain_full_rectangular) then !read the body points of the rigid stationary body
            do i=1,m_body
                if (bdy(i)%fam==1) then
                   sd_x(nextb+1:nextb+bdy(i)%npts) = bdy(i)%x(:)
                   sd_y(nextb+1:nextb+bdy(i)%npts) = bdy(i)%y(:)
                   nextb = nextb + bdy(i)%npts
                end if
            end do
        end if
        
        !Assign coordinates to the rectangular sub-domain
        do j=1, nint((sdye-sdyi)/delta_sd+1)
            do i=1, nint((sdxe-sdxi)/delta_sd+1)   !loop along x-direction first
                nextb = nextb+1
                sd_x(nextb) = sdxi + (i-1)*delta_sd
                sd_y(nextb) = sdyi + (j-1)*delta_sd
            end do
        end do
        nf_sd = 2*nb_sd  !number of forces on sub-domain points
    
    end if
    
    !Indexing for forces, positions, and velocities
    nf  = 2*nb                     ! number of forces
   
    call force_proc(nf, nf_sd, xfm, xsdm) !force vector partitioning
    call body_proc !Force vector scattering book-keeping
    call list_bodies !Force vector scattering continued

    !Initialize force vector fvec
      call VecCreateMPI(MPI_COMM_WORLD, xfm, nf, fvec, ierr)
      call VecGetOwnershipRange(fvec, xfi, xfe, ierr)
      xfi=xfi+1
      call VecDuplicate(fvec, xbvec, ierr)
      call VecDuplicate(fvec, ybvec, ierr)
      allocate(codeb(xfi:xfe,3))
      
      !Initialize xb to the positions read from the files.  If this is a restart, this will be overwritten later
      call VecGetArrayF90(xbvec, xb, ierr)
      call VecGetArrayF90(ybvec, yb, ierr)      
      call collect_bodies(xb, yb, codeb)
            
      if (rank==0) write(*,*) 'setup global positions, velocities, and forces'
      call setup_pre_reg   !only done once, just to get some max index
      CALL setup_reg(xb, yb)
      if (rank==0) write(*,*) 'setup regularization of initial geometry'
      call VecRestoreArrayF90(xbvec, xb, ierr)
      call VecRestoreArrayF90(ybvec, yb, ierr)
      
      if (sub_domain) then
          !initialize points on the subdomain
          call VecCreateMPI(MPI_COMM_WORLD, xsdm, nf_sd, xsdvec, ierr)
          call VecDuplicate(xsdvec, ysdvec, ierr)
          call VecGetOwnershipRange(xsdvec, xsdi, xsde, ierr)
          xsdi=xsdi+1
          allocate(codesd(xsdi:xsde,1))
      
         call VecGetArrayF90(xsdvec, xsd, ierr)
         call VecGetArrayF90(ysdvec, ysd, ierr)      
         call collect_subdomain(xsd, ysd, codesd)
         
         if (rank==0) write(*,*) 'setup subdomain'
         call setup_pre_reg_subdomain(xsd, ysd)   !Computes Eo
         call setup_interp_indices_scatter
         if (rank==0) write(*,*) 'setup regularization of subdomain'
         call VecRestoreArrayF90(xsdvec, xsd, ierr)
         call VecRestoreArrayF90(ysdvec, ysd, ierr) 
      end if      
    
  END SUBROUTINE setup_geometry
  
  
!===================================================================================    

  SUBROUTINE  collect_bodies(xb, yb, codeb)
  !assign body positions into xbvec and ybvec. xbvec and ybvec have the same structure as that of fvec. See the paper Nair and Goza 2021 for details about how fvec is structured. However, in xbvec and ybvec, the body coordinates are repeated, since ideally xbvec and ybvec have sizes half of that of fvec, but is constructed to have the same size of fvec. 

  !codeb is book-keeping matrix having nf rows for each body unknown in fvec. First column denotes the body number assigned to that body point. Second column denotes the family of that body point. Third column denotes whether it is the x or y coordinate.
  
    PetscScalar :: xb(xfi:xfe), yb(xfi:xfe)
    PetscInt :: codeb(xfi:xfe,3), i, j, i_bdy
    
      do i=xfi,xfe
        do j=m_body,1,-1
          if (i.le.nbb_ind(j)) then
            i_bdy=j
          end if
        end do
        
        codeb(i,2)=bdy(i_bdy)%fam
        codeb(i,1)=bdy(i_bdy)%bn
        
        if (i-nbb_ind(i_bdy-1).le.nb_ind(i_bdy)) then
          codeb(i,3)=1
          xb(i)=bdy(i_bdy)%x(i-nbb_ind(i_bdy-1))
          yb(i)=bdy(i_bdy)%y(i-nbb_ind(i_bdy-1))
        else 
          codeb(i,3)=2
          xb(i)=bdy(i_bdy)%x(i-nbb_ind(i_bdy-1)-nb_ind(i_bdy))
          yb(i)=bdy(i_bdy)%y(i-nbb_ind(i_bdy-1)-nb_ind(i_bdy))
        end if

      end do

  END SUBROUTINE  collect_bodies
  
!===================================================================================    

  SUBROUTINE  collect_subdomain(xsd, ysd, codesd)
  !same as collect_bodies but for the sub-domain
  
      PetscScalar :: xsd(xsdi:xsde), ysd(xsdi:xsde)
      PetscInt :: codesd(xsdi:xsde,1), i
    
      do i=xsdi,xsde
        if (i <= nb_sd) then
          codesd(i,1)=1
          xsd(i) = sd_x(i)
          ysd(i) = sd_y(i)
        else 
          codesd(i,1)=2
          xsd(i) = sd_x(i-nb_sd)
          ysd(i) = sd_y(i-nb_sd)
        end if
      end do

  END SUBROUTINE  collect_subdomain
    
!=================================================================================== 

  SUBROUTINE setup_reg(xb, yb)
  !Here we compute delta function weights stored in E and Et operators. The weights are given by matrix "weights" while the location index is given by reg_ix and reg_iy
    PetscScalar :: xb(xfi:xfe), yb(xfi:xfe)
    PetscScalar :: d2, x, y
    PetscInt :: k, l, p, next
    PetscInt :: iterx, itery
    AO :: aou, aov
  
  
    support2= ceiling(support)
    d2 = 0.5d0*delta
    
    IF(ALLOCATED(indexx)) THEN
        DEALLOCATE(indexx)
    END IF
    
    IF(ALLOCATED(indexy)) THEN
        DEALLOCATE(indexy)
    END IF

    IF(ALLOCATED(weight)) THEN
        DEALLOCATE(weight)
    END IF
    
    ALLOCATE(weight((2*support2+1)*(2*support2+1),xfi:xfe))
    ALLOCATE(indexx(xfi:xfe), indexy(xfi:xfe))
    
    iterx=0
    itery=0
    DO k=xfi,xfe
        
        !get index of body position relative to grid
        indexx(k)=INT((xb(k)+offsetx)/delta) + 1
        indexy(k)=INT((yb(k)+offsety)/delta) + 1
        next=0
        DO l=-support2,support2
            DO p=-support2,support2
            
                if (codeb(k,3)==1) then
                  x=delta*(indexx(k)-1+p)-offsetx ! grid location x
                  y=delta*(indexy(k)-1+l)-offsety+d2 ! grid location y
                  iterx=iterx+1
                  reg_ix(iterx) = (m+1)*(indexy(k)+l-1) + (indexx(k)+p-1)    !(i+p-1) because of 0 indexing for vecsetvalues
                else if (codeb(k,3)==2) then
                  x=delta*(indexx(k)-1+p)-offsetx+d2 ! grid location x
                  y=delta*(indexy(k)-1+l)-offsety ! grid location y
                  itery=itery+1
                  reg_iy(itery)= (m)*(indexy(k)+l-1) + (indexx(k)+p-1)
                end if
                next=next+1
                !get delta function weights
                weight(next,k)= delta * delta * &    
                deltafnc(x,xb(k),delta) * &
                deltafnc(y,yb(k),delta)
                
            end DO
        end DO
    END DO
    
    call DMDAGetAO(dau, aou, ierr)
    call AOApplicationToPetsc(aou, interfacex, reg_ix, ierr)
    
    call DMDAGetAO(dav, aov, ierr)
    call AOApplicationToPetsc(aov, interfacey, reg_iy, ierr)
    
    call ISCreateGeneral(MPI_COMM_SELF, interfacex, reg_ix, PETSC_COPY_VALUES, reg_ixs, ierr)
    call ISCreateGeneral(MPI_COMM_SELF, interfacey, reg_iy, PETSC_COPY_VALUES, reg_iys, ierr)
    
    
  END SUBROUTINE setup_reg

!=================================================================================== 

  SUBROUTINE setup_pre_reg
  !This is a precursor to setup_reg

    PetscScalar :: d2
    PetscInt :: k, l, p, i
    PetscInt, allocatable :: interfacex_ix(:), interfacey_iy(:)
  
  
    support2= ceiling(support)
    d2 = 0.5d0*delta    
    
    ! get regularized weight near ib points (u-vel points)
    interfacex=0 !This is for allocating appropriate array sizes in reg function
    interfacey=0 !Physically these represent the total number of values to send from the structure to the fluid
    DO k=xfi,xfe
        DO l=-support2,support2
            DO p=-support2,support2
            
                if (codeb(k,3)==1) then
                  interfacex=interfacex+1
                else if (codeb(k,3)==2) then
                  interfacey=interfacey+1
                end if
                
            end DO
        end DO
    END DO

    allocate(reg_ix(interfacex), reg_iy(interfacey), interfacex_ix(interfacex), interfacey_iy(interfacey))
    
    interfacex_ix=(/(i, i=0,interfacex-1, 1)/)
    interfacey_iy=(/(i, i=0,interfacey-1, 1)/)  

    call ISCreateGeneral(MPI_COMM_SELF, interfacex, interfacex_ix, PETSC_COPY_VALUES, interfacex_is, ierr)
    call ISCreateGeneral(MPI_COMM_SELF, interfacey, interfacey_iy, PETSC_COPY_VALUES, interfacey_is, ierr)    
    
  END SUBROUTINE setup_pre_reg
  
!=================================================================================== 

  SUBROUTINE setup_reg_subdomain(sd_list_size, sd_list)!some array to be equated with interp_indices_global)
  !Identifying the sparse reg and regT indices and weights using global indices that will be used for efficiently computing PEoq
  !Outputs are reg_ixs_sd_sparse and interfacex_sd_is_sparse
  
    integer, intent(in) :: sd_list_size
    PetscInt, intent(in) :: sd_list(sd_list_size)
    PetscInt :: k, l, p, next, kk, i
    PetscInt :: iterx, itery
    AO :: aou, aov
    PetscInt, allocatable :: interfacex_ix(:), interfacey_iy(:)

    interfacex_sd_sparse = count(sd_list.le.nb_sd .and. sd_list.ge.xsdi .and. sd_list.le.xsde)*(2*support2+1)**2
    interfacey_sd_sparse = count(sd_list.gt.nb_sd .and. sd_list.ge.xsdi .and. sd_list.le.xsde)*(2*support2+1)**2
    if (allocated(reg_ix_sd_sparse)) deallocate(reg_ix_sd_sparse)
    if (allocated(reg_iy_sd_sparse)) deallocate(reg_iy_sd_sparse)
    allocate(interfacex_ix(interfacex_sd_sparse), interfacey_iy(interfacey_sd_sparse), reg_ix_sd_sparse(interfacex_sd_sparse), reg_iy_sd_sparse(interfacey_sd_sparse))
    
    iterx=0
    itery=0    
    interp_indices_global_nproc = 0
    do kk=1,sd_list_size
        k = sd_list(kk)
        if (k.ge.xsdi .and. k.le.xsde) then
        next=0
        interp_indices_global_nproc = interp_indices_global_nproc+1  !this is only important when sd_list is interp_indices_global
        DO l=-support2,support2
            DO p=-support2,support2
                next=next+1
            
                if (codesd(k,1)==1) then
                  iterx=iterx+1
                  reg_ix_sd_sparse(iterx) = reg_ix_sd((k-xsdi)*(2*support2+1)**2+next)
                else if (codesd(k,1)==2) then
                  itery=itery+1
                  reg_iy_sd_sparse(itery) = reg_iy_sd((k-xsdi)*(2*support2+1)**2-interfacex_sd+next)
                end if                
            end DO
        end DO
        end if
    END DO
    
    call ISCreateGeneral(MPI_COMM_SELF, interfacex_sd_sparse, reg_ix_sd_sparse, PETSC_COPY_VALUES, reg_ixs_sd_sparse, ierr)
    call ISCreateGeneral(MPI_COMM_SELF, interfacey_sd_sparse, reg_iy_sd_sparse, PETSC_COPY_VALUES, reg_iys_sd_sparse, ierr)
      
    interfacex_ix=(/(i, i=0,interfacex_sd_sparse-1, 1)/)
    interfacey_iy=(/(i, i=0,interfacey_sd_sparse-1, 1)/)  

    call ISCreateGeneral(MPI_COMM_SELF, interfacex_sd_sparse, interfacex_ix, PETSC_COPY_VALUES, interfacex_sd_is_sparse, ierr)
    call ISCreateGeneral(MPI_COMM_SELF, interfacey_sd_sparse, interfacey_iy, PETSC_COPY_VALUES, interfacey_sd_is_sparse, ierr)

    deallocate(interfacex_ix, interfacey_iy) 
    
  END SUBROUTINE setup_reg_subdomain
  
!=================================================================================== 

  SUBROUTINE setup_pre_reg_subdomain(xsd, ysd)
  !This precomputes Eo

    PetscScalar :: xsd(xsdi:xsde), ysd(xsdi:xsde)
    PetscScalar :: d2, x, y
    PetscInt :: k, l, p, i, next
    PetscInt :: iterx, itery
    AO :: aou, aov
    
    support2= ceiling(support)
    d2 = 0.5d0*delta
    
    IF(ALLOCATED(indexx_sd)) THEN
        DEALLOCATE(indexx_sd)
    END IF
    
    IF(ALLOCATED(indexy_sd)) THEN
        DEALLOCATE(indexy_sd)
    END IF

    IF(ALLOCATED(weight_sd)) THEN
        DEALLOCATE(weight_sd)
    END IF
    
    ALLOCATE(weight_sd((2*support2+1)*(2*support2+1),xsdi:xsde))
    ALLOCATE(indexx_sd(xsdi:xsde), indexy_sd(xsdi:xsde))

    interfacex_sd = count( (/(i, i=xsdi,xsde, 1)/) .le. nb_sd)*(2*support2+1)**2
    interfacey_sd = count( (/(i, i=xsdi,xsde, 1)/) .gt. nb_sd)*(2*support2+1)**2
    
    allocate(reg_ix_sd(interfacex_sd), reg_iy_sd(interfacey_sd))  
    
    iterx=0
    itery=0
    ! get regularized weight near subdomain grid points (u-vel points)
    DO k=xsdi,xsde
        
        ! get index of body position relative to grid
        indexx_sd(k)=INT((xsd(k)+offsetx)/delta) + 1
        indexy_sd(k)=INT((ysd(k)+offsety)/delta) + 1
        next=0
        DO l=-support2,support2
            DO p=-support2,support2
            
                if (codesd(k,1)==1) then
                  x = delta*(indexx_sd(k)-1+p)-offsetx ! grid location x
                  y = delta*(indexy_sd(k)-1+l)-offsety+d2 ! grid location y
                  iterx=iterx+1
                  reg_ix_sd(iterx) = (m+1)*(indexy_sd(k)+l-1) + (indexx_sd(k)+p-1)    !(i+p-1) because of 0 indexing for vecsetvalues
                else if (codesd(k,1)==2) then
                  x = delta*(indexx_sd(k)-1+p)-offsetx+d2 ! grid location x
                  y = delta*(indexy_sd(k)-1+l)-offsety ! grid location y
                  itery=itery+1
                  reg_iy_sd(itery) = (m)*(indexy_sd(k)+l-1) + (indexx_sd(k)+p-1)
                end if
                next=next+1
                weight_sd(next,k) = delta * delta * &    
                deltafnc(x,xsd(k),delta) * &
                deltafnc(y,ysd(k),delta)
                
            end DO
        end DO
    END DO
    
  END SUBROUTINE setup_pre_reg_subdomain
  
!=================================================================================== 

  SUBROUTINE setup_subdomain_sparsity(xsd, ysd, support_nz, regTmat)
  !creating E matrix only for some sparsity estimates as described in get_sparsity_EEt

    PetscScalar :: xsd(xsdi:xsde), ysd(xsdi:xsde)
    PetscInt :: support_nz
    Mat :: regTmat
    PetscScalar :: d2, x, y
    PetscInt :: k, l, p, i, next, kk(1)
    PetscInt :: iterx, itery, temp_indexx_sd, temp_indexy_sd
    AO :: aou, aov
    
    d2 = 0.5d0*delta

    IF(ALLOCATED(weight_sd_sparsity)) THEN
        DEALLOCATE(weight_sd_sparsity)
    END IF
    
    interface_sd_sparsity = (2*support_nz+1)**2
    ALLOCATE(weight_sd_sparsity(1,interface_sd_sparsity))
    allocate(reg_i_sd_sparsity(interface_sd_sparsity))
    
    DO k=1,nf_sd
        
        ! get index of body position relative to grid
        temp_indexx_sd=INT((xsd(k)+offsetx)/delta) + 1
        temp_indexy_sd=INT((ysd(k)+offsety)/delta) + 1
        next=0
        DO l=-support_nz,support_nz
            DO p=-support_nz,support_nz
                next=next+1
                if (k .le. nb_sd) then
                  x = delta*(temp_indexx_sd-1+p)-offsetx ! grid location x
                  y = delta*(temp_indexy_sd-1+l)-offsety+d2 ! grid location y
                  reg_i_sd_sparsity(next) = (m+1)*(temp_indexy_sd+l-1) + (temp_indexx_sd+p-1)    !(i+p-1) because of 0 indexing for vecsetvalues
                else if (k .gt. nb_sd) then
                  x = delta*(temp_indexx_sd-1+p)-offsetx+d2 ! grid location x
                  y = delta*(temp_indexy_sd-1+l)-offsety ! grid location y
                  reg_i_sd_sparsity(next) = (m+1)*n-1 + (m)*(temp_indexy_sd+l-1) + (temp_indexx_sd+p-1)
                end if
                
                weight_sd_sparsity(1,next) = delta * delta * &    
                deltafnc(x,xsd(k),delta) * &
                deltafnc(y,ysd(k),delta)
                
            end DO
        end DO
        kk(1) = k-1
        call MatSetValues(regTmat, ione, kk, interface_sd_sparsity, reg_i_sd_sparsity, weight_sd_sparsity, INSERT_VALUES, ierr)
    END DO
    call MatAssemblyBegin(regTmat, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(regTmat, MAT_FINAL_ASSEMBLY, ierr)
    
    deallocate(reg_i_sd_sparsity, weight_sd_sparsity)
   
    
  END SUBROUTINE setup_subdomain_sparsity  
    
!===================================================================================

  subroutine setup_interp_indices_scatter
  !See setup_reg_subdomain_indices for details for various terms defined here. (For understanding the code, you can ignore this subroutine for now).
  
      call VecCreateMPI(PETSC_COMM_WORLD, xfm*interp_points, nf*interp_points, interp_indices_vec, ierr)
      call VecScatterCreateToAll(interp_indices_vec, ctx_indices, interp_indices_global_vec, ierr)
      
      call VecCreateMPI(PETSC_COMM_WORLD, xfm*interp_points, nf*interp_points, interp_weights_vec, ierr)
      call VecScatterCreateToAll(interp_weights_vec, ctx_weights, interp_weights_global_vec, ierr)
      
      allocate(interp_indices_global_sparsemat(interp_points,nf), interp_weights_global_sparsemat(interp_points,nf))
      
      call VecDuplicate(xsdvec, intermediatevec, ierr)
      call VecCreateSeq(MPI_COMM_SELF, nf_sd, intermediatelocalvec, ierr)
      interp_indices_local_size = 0
      allocate(interp_indices_local(interp_indices_local_size))
      call ISCreateGeneral(MPI_COMM_SELF, interp_indices_local_size, interp_indices_local, PETSC_COPY_VALUES, interp_indices_local_is, ierr)
      call VecScatterCreate(intermediatevec, interp_indices_local_is, intermediatelocalvec, interp_indices_local_is, ctx_intermediate, ierr)
  
  end subroutine setup_interp_indices_scatter
  
!===================================================================================  

  subroutine setup_reg_subdomain_indices(xb, yb)
        
  !This subroutine consists of generating P matrix to be used in the sub-domain approach. P is characterised by interp_weights and interp_indices that contain the weights and indices for interpolation. They are local to each processor. Also, Eo indices and weights compatible with P to enable efficient computation of PEoq are identified
  
  !For performing PBP^T efficiently, the non-partitioned P^T needs to be present redundantly on all processors. To facilitate this, first, the local interp_indices are collected intto MPI vecs interp_indices_vec, then scattered (using ctx_indices) fully to Seq vecs interp_indices_global_vec, which are then extracted as arrays into interp_indices_global_sparsemat. The same procedure is done on interp_weights. Finally, P^T is represented by interp_indices_global_sparsemat and interp_weights_global_sparsemat.

  !For efficiently performing P^T*f in step 3-corrector, P^T is represented by indicesserialmat.
  !Some other important variables are 
  !interp_indices_global --- contains *ALL* unique indices of P as an array, redundantly present in all processors (basically unique subroutine applied to interp_indices_global_sparsemat)
  !interp_indices_global_nproc --- this is the partitioned interp_indices_global
  !interp_indices_local --- contains unique indices of P as an array present in the local processor (basically unique subroutine applied to interp_indices)
  
  !When performing P*(BP^T), the columns of (BP^T) needs to be commnicated or scattered appropriately. These columns are first collected into intermediatevec, then scattered to intermediatelocalvec via ctx_intermediate

      
      PetscScalar :: xb(xfi:xfe), yb(xfi:xfe)
      PetscScalar :: dist(nb_sd), distance(interp_points)
      PetscInt :: i, j, i1, i2
      integer :: nearest_neighbor(interp_points), interp_indices_global_duplicates(nf*interp_points), interp_indices_local_duplicates(xfm*interp_points)
      PetscScalar :: interp_weights_local(xfm*interp_points)
      Logical, dimension(interp_points) :: mk  
      PetscInt :: k, l, p, next, iterx, itery, kk, ll(nint(sqrt(real(interp_points)))), pp(nint(sqrt(real(interp_points))))
      PetscInt, allocatable :: interfacex_ix(:), interfacey_iy(:)
      PetscScalar :: indxx_sd, indxy_sd
      Logical, dimension(nint(sqrt(real(interp_points)))) :: mkx, mky
      PetscScalar, dimension(nint(sqrt(real(interp_points)))) :: distancex, distancey
      double precision infor(MAT_INFO_SIZE), nz_a, nz_unneed, nz_use
      PetscScalar, pointer :: an_array(:)
      
      IF(.not. ALLOCATED(interp_weights)) THEN
          ALLOCATE(interp_weights(interp_points, xfi:xfe), interp_indices(interp_points, xfi:xfe))
      END IF
      
      IF(ALLOCATED(interp_indices_global)) THEN
          DEALLOCATE(interp_indices_global)
      END IF
      
      !i1 and i2 are just i1=0 and i2=1 for bilinear interpolation using np=4 points
      if (mod(interp_points,2)==0) then  !if even number
          i1 = sqrt(real(interp_points))/2 + 1 - sqrt(real(interp_points))
          i2 = i1 + sqrt(real(interp_points)) - 1
      else
          i1 = (sqrt(real(interp_points))-1)/2 + 1 - sqrt(real(interp_points))
          i2 = i1 + sqrt(real(interp_points)) - 1
      end if
      
      do j=xfi,xfe
          !Find nearest index indxx_sd and indxy_sd on the sub-domain for the body point and identify indices ll and pp for the nearest neighboring sub-domain points that form the tensor grid
          indxx_sd = (xb(j) - sdxi)/delta_sd + 1
          indxy_sd = (yb(j) - sdyi)/delta_sd + 1
          
          if (mod(interp_points,2)==0) then  !if even number
              ll = (/(int(indxx_sd)+i,i=i1,i2)/)   !ll and pp are the indices of the interpolation points in the subdomain relative to the subdomain and not the flow domain
              pp = (/(int(indxy_sd)+i,i=i1,i2)/)
          else
              ll = (/(nint(indxx_sd)+i,i=i1,i2)/)
              pp = (/(nint(indxy_sd)+i,i=i1,i2)/)
          end if
          
          if (.not. sub_domain_full_rectangular .and. codeb(j,2)==1) then  !if exact body points are included in sub-domain points---!assumes ONLY one rigid body present in the problem
              if (codeb(j,3)==1) then  !x-direction
                  nearest_neighbor = j+1  !dummy nearest neighbor if exact body points are used in sub-domain points
                  nearest_neighbor(1) = j  !choose the body point itself
              else  !y-direction
                  nearest_neighbor(:) = j+1-nb_ind(codeb(j,1))
                  nearest_neighbor(1) = j-nb_ind(codeb(j,1))
              end if
              distance(:) = 1
              distance(1) = 0  !Set interpolation weights indirectly to zero via the variable distance
          
          else
          i=0
          do l = 1,nint(sqrt(real(interp_points)))
              do p = 1,nint(sqrt(real(interp_points)))
                  i = i+1
                  IF (ll(l)<1) then
                      STOP 'Body has exited the subdomain from left'
                  else if (ll(l)>nint((sdxe-sdxi)/delta_sd+1)) then
                      STOP 'Body has exited the subdomain from right'
                  else if (pp(p)<1) then
                      STOP 'Body has exited the subdomain from bottom'
                  else if (pp(p)>nint((sdye-sdyi)/delta_sd+1) ) then
                      STOP 'Body has exited the subdomain from top'
                  END IF
                  
                  nearest_neighbor(i) = nr_sd + (pp(p)-1)*nint((sdxe-sdxi)/delta_sd+1) + ll(l)
                  distancex(l) = xb(j)-sd_x(nearest_neighbor(i))
                  distancey(p) = yb(j)-sd_y(nearest_neighbor(i))
                  distance(i) = sqrt( (distancex(l))**2.D0 + (distancey(p))**2.D0 )
              end do
          end do
          end if

          mk = .TRUE.
          !do i=1,interp_points
          i=0
          do l = 1,nint(sqrt(real(interp_points)))
          do p = 1,nint(sqrt(real(interp_points)))
              i = i+1
              !Interpolation weights and indices using linear hat delta function
              interp_weights(i,j)= delta_sd * delta_sd * deltafnc2(sd_x(nearest_neighbor(i)),xb(j),delta_sd) * deltafnc2(sd_y(nearest_neighbor(i)),yb(j),delta_sd)
              
              if (distance(i)==0) then 
                   interp_weights(i,j)=1.d0  !Set interpolation weight to 0 if the body exact exists on the sub-domain
                   mk(i) = .FALSE.
              end if
              if (codeb(j,3)==1) then
                  interp_indices(i,j) = nearest_neighbor(i)
              else
                  interp_indices(i,j) = nearest_neighbor(i) + nb_sd
              end if
          end do  !loop l
          end do  !loop p
          
          if (any(.not. mk)) where(mk .eqv. .TRUE.) interp_weights(:,j)=0.d0
      end do
      
      !need to scatter all to all interp_indices
      !in the meantime need to compute interp_indices_local
      !when scatter all to all is done, need to compute interp_indices_global which will be the same in all the procs
      !that global vec will be fed to setup_reg_subdomain
      
      !Collecting all the indices locally
      interp_indices_local_duplicates = reshape(interp_indices, (/xfm*interp_points/))
      call VecGetArrayF90(interp_indices_vec, an_array, ierr)
      an_array(1:xfm*interp_points) = interp_indices_local_duplicates(1:xfm*interp_points)
      call VecRestoreArrayF90(interp_indices_vec, an_array, ierr)
      !Scattering the local indices to a sequential global accessible to all procs
      call VecScatterBegin(ctx_indices, interp_indices_vec, interp_indices_global_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
      
      !Collecting all the weights locally
      interp_weights_local = reshape(interp_weights, (/xfm*interp_points/))
      call VecGetArrayF90(interp_weights_vec, an_array, ierr)
      an_array(1:xfm*interp_points) = interp_weights_local(1:xfm*interp_points)
      call VecRestoreArrayF90(interp_weights_vec, an_array, ierr)
      !Scattering the local weights to a sequential global accessible to all procs
      call VecScatterBegin(ctx_weights, interp_weights_vec, interp_weights_global_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)

      if (allocated(interp_indices_local)) deallocate(interp_indices_local)
      interp_indices_local = unique(interp_indices_local_duplicates)
      interp_indices_local_size = size(interp_indices_local)
      !Creating an indicesmat locally
      if (allocated(indicesserialmat)) deallocate(indicesserialmat)
      allocate(indicesserialmat(interp_indices_local_size,xfi:xfe))
      indicesserialmat(:,:) = 0.D0
      do j=xfi,xfe
          do i=1,interp_points
              where (interp_indices_local .eq. interp_indices(i,j)) indicesserialmat(:, j) = interp_weights(i,j)
          end do
      end do
      
      !Creating scatter context for intermediate matmatmult in compute_operator
      call ISDestroy(interp_indices_local_is, ierr)
      call VecScatterDestroy(ctx_intermediate, ierr)
      call ISCreateGeneral(MPI_COMM_SELF, interp_indices_local_size, interp_indices_local-1, PETSC_COPY_VALUES, interp_indices_local_is, ierr)
      call VecScatterCreate(intermediatevec, interp_indices_local_is, intermediatelocalvec, interp_indices_local_is, ctx_intermediate, ierr)  
      
      !Collecting all the global indices
      call VecScatterEnd(ctx_indices, interp_indices_vec, interp_indices_global_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecGetArrayF90(interp_indices_global_vec, an_array, ierr)
      interp_indices_global_duplicates = an_array
      call VecRestoreArrayF90(interp_indices_global_vec, an_array, ierr)
      if (allocated(interp_indices_global)) deallocate(interp_indices_global)
      interp_indices_global_sparsemat = reshape(interp_indices_global_duplicates, (/interp_points,nf/)) !this line needs to be before unique is applied, since unque modifies input
      interp_indices_global = unique(interp_indices_global_duplicates)
      
      !Identifying the sparse reg and regT indices and weights using global indices.
      CALL setup_reg_subdomain(size(interp_indices_global),interp_indices_global)
      
      call VecScatterEnd(ctx_weights, interp_weights_vec, interp_weights_global_vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecGetArrayF90(interp_weights_global_vec, an_array, ierr)
      interp_weights_global_sparsemat = reshape(an_array, (/interp_points,nf/))
      call VecRestoreArrayF90(interp_indices_global_vec, an_array, ierr)      
  
  end subroutine setup_reg_subdomain_indices

!===================================================================================  

FUNCTION deltafnc( x1, x2, dr )

real(kind(0.d0)) :: x1, x2
REAL(KIND(0.D0)) :: r,dr,deltafnc
real(kind(0.d0)) :: r1, r2, r3, r4, a1, a5, a6, a7, a8, a9

r = abs( x1 - x2 )

r1 = r/dr
r2 = r1*r1
r3 = r2*r1
r4 = r3*r1

if (r1 .le. 1.d0) then
a5 = asin((1.d0/2.d0)*sqrt(3.d0)*(2.d0*r1-1.d0))
a8 = sqrt(1.d0-12.d0*r2+12.d0*r1)

deltafnc = 0.4166666667d-1*r4+(-.1388888889+0.3472222222d-1*a8)*r3+ &
(-0.7121664902d-1-0.5208333333d-1*a8+0.2405626122*a5)*r2+&
(-.2405626122*a5-.3792313933+.1012731481*a8)*r1+0.8018753741d-1*a5 &
-0.4195601852d-1*a8+.6485698427

elseif (r1 .le. 2.d0) then
a6 = asin((1.d0/2.d0)*sqrt(3.d0)*(-3.d0+2.d0*r1))
a9 = sqrt(-23.d0+36.d0*r1-12.d0*r2)

deltafnc = -0.6250000000d-1*r4+(.4861111111-0.1736111111d-1*a9)*r3 + &
(-1.143175026+0.7812500000d-1*a9-.1202813061*a6)*r2 + &
(.8751991178+.3608439183*a6-.1548032407*a9)*r1-.2806563809*a6 + &
0.822848104d-2+.1150173611*a9

elseif (r1 .le. 3.d0 ) then
a1 = asin((1.d0/2.d0*(2.d0*r1-5.d0))*sqrt(3.d0))
a7 = sqrt(-71.d0-12.d0*r2+60.d0*r1)

deltafnc = 0.2083333333d-1*r4+(0.3472222222d-2*a7-.2638888889)*r3+ &
(1.214391675-0.2604166667d-1*a7+0.2405626122d-1*a1)*r2+ &
(-.1202813061*a1-2.449273192+0.7262731481d-1*a7)*r1 +.1523563211*a1 &
+1.843201677-0.7306134259d-1*a7
!print *, deltafnc

else
deltafnc = 0.d0

end if

deltafnc = deltafnc / dr

END FUNCTION deltafnc  

!===================================================================================  

FUNCTION deltafnc2( x1, x2, dr ) result (deltafnc)
!2-point hat function

real(kind(0.d0)) :: x1, x2
REAL(KIND(0.D0)) :: r,dr,deltafnc
real(kind(0.d0)) :: r1, r2, r3, r4, a1, a5, a6, a7, a8, a9

r = abs( x1 - x2 )

r1 = r/dr
r2 = r1*r1
r3 = r2*r1
r4 = r3*r1

if (r1 .le. 1.d0) then
deltafnc = 1-r1
else
deltafnc = 0.d0
end if

deltafnc = deltafnc / dr

END FUNCTION deltafnc2

!===================================================================================  

FUNCTION deltafnc4( x1, x2, dr ) result(deltafnc)
!piecewise cubic function

real(kind(0.d0)) :: x1, x2
REAL(KIND(0.D0)) :: r,dr,deltafnc
real(kind(0.d0)) :: r1, r2, r3, r4, a1, a5, a6, a7, a8, a9

r = abs( x1 - x2 )

r1 = r/dr
r2 = r1*r1
r3 = r2*r1
r4 = r3*r1

if (r1 .le. 1.d0) then
deltafnc = 1 - 0.5d0*r1 - r2 + 0.5d0*r3
elseif (r1 .le. 2.d0) then
deltafnc = 1 - 1.83333333d0*r1 + r2 - 0.16666666667d0*r3
else
deltafnc = 0.d0
end if

deltafnc = deltafnc / dr

END FUNCTION deltafnc4

!===================================================================================

  subroutine setup_coarsify
  !Some indices for required for communication in coarsify
  
    PetscInt :: xcoarse, ycoarse
    
    PetscInt :: iter, j, i, jd, id  !delete
    PetscInt, allocatable, dimension(:) :: ix  !delete
    AO :: ao
    
  
    if (modulo(xsi,2)==0) then
      xci=xsi
    else
      xci=xsi-1
    end if
    
    if (modulo(ysi,2)==0) then
      yci=ysi
    else
      yci=ysi-1
    end if
    
    xcoarse=(xse-xci+1)/2
    ycoarse=(yse-yci+1)/2
    tcoarse=xcoarse*ycoarse
    
  end subroutine setup_coarsify

!===================================================================================

  subroutine setup_getbc
  !Some indices for required for communication in get_bc
  
    PetscInt :: i,j, next, tis
        
    !Finding indices of boundaries on coarse grid (cbcix) that needs to be transferred to the fine grid
    next=0
    do i=xsi,xse
      do j=ysi,yse
      
        if (j==n/4+1   .and. i>=m/4+1 .and. i<=3*m/4+1) then
          next=next+1
          call grow_array(cbcix, bottom+2*(i-1-m/4), next)   !-1 from i for C indexing
          if (i.ne.3*m/4+1) then
            next=next+1
            call grow_array(cbcix, bottom+2*(i-1-m/4) + 1, next)
          end if
        end if
        
        if (j==3*n/4+1 .and. i>=m/4+1 .and. i<=3*m/4+1) then
          next=next+1
          call grow_array(cbcix, top+2*(i-1-m/4), next)   !-1 from i for C indexing
          if (i.ne.3*m/4+1) then
            next=next+1
            call grow_array(cbcix, top+2*(i-1-m/4) + 1, next)
          end if
        end if
        
        if (i==m/4+1   .and. j>=n/4+1 .and. j<=3*n/4+1) then
          next=next+1
          call grow_array(cbcix, left+2*(j-1-n/4), next)   !-1 from j for C indexing
          if (j.ne.3*n/4+1) then
            next=next+1
            call grow_array(cbcix, left+2*(j-1-n/4) + 1, next)
          end if
        end if
        
        if (i==3*m/4+1 .and. j>=n/4+1 .and. j<=3*n/4+1) then
          next=next+1
          call grow_array(cbcix, right+2*(j-1-n/4), next)   !-1 from j for C indexing
          if (j.ne.3*n/4+1) then
            next=next+1
            call grow_array(cbcix, right+2*(j-1-n/4) + 1, next)
          end if
        end if
        
      end do
    end do
    
    tbc=next  !total number of boundary points that needs to be transferred
    
    !Finding indices of boundary points on fine grid (bcix) for scattering from global bc vec to local bc vec
    next=0
    do i=xsgi,xsge
      do j=ysgi,ysge
      
        if (i==1) then
          next=next+1
          call grow_array(bcix, left+j-1, next)
        end if
        
        if (i==m+1) then
          next=next+1
          call grow_array(bcix, right+j-1, next)
        end if
        
        if (j==1) then
          next=next+1
          call grow_array(bcix, bottom+i-1, next)
        end if
        
        if (j==n+1) then
          next=next+1
          call grow_array(bcix, top+i-1, next)
        end if
       
      end do
    end do
    tis=next
    
    call ISCreateGeneral(MPI_COMM_SELF, tis, bcix, PETSC_COPY_VALUES, bcis, ierr)
    

  end subroutine setup_getbc


!===================================================================================  

  subroutine grow_array(yyy, xxx, nxt)
  !yyy is the input array needed to grow in size, xxx is the value to be appended to yyy, nxt is the length of the grown yyy
  
    PetscInt :: nxt, xxx
    PetscInt, allocatable :: temp(:), yyy(:)
  
    allocate(temp(1:nxt))
    if ( allocated (yyy) ) temp(1:nxt-1) = yyy
    temp(nxt)=xxx
    call move_alloc(temp, yyy)
  
  end subroutine grow_array
  
!================================================================================

  subroutine body_proc
  !Determine how the indiidual body operators should be divided
  !Basically three levels of scatter are used for dealing with forces on the interface.

  !Scatter level 1: The body force vector fvec and similar vectors that exist on a subset of processors are first scattered to sequential vectors local to each processor. That is, component of fvec corresponding to the first body will be scattered to proc 0 and will locally exist on proc 0. For second body, fvec will be scattered to proc 1 and so on. If there are more bodies than processors, then the round of assignment will be repeated starting from proc 0, e.g. (nproc+1)st body will be assigned to proc 0. All the scattering operations in this level is collective that occurs on MPI_COMM_WORLD.

  !Scatter level 2: After havinng scattered a parallel vector to sequential vectors, the sequential vectors may be further divided among any "extra" available procs that were not involved in scatter 1. This scattering is performed on a comm_body communicator which consists of the the main proc and the extra procs to which the force vector is scattered.

  !(Virtual) Scatter level 3: This is only performed if the main proc after scatter 1 holds more than 1 body. Truly speaking, a scatter is not practically performed here. Essentially only looping through the bodies of that main proc is done.

  !Important outputs from this subroutine are:
  !rank_body: It is an array local to each processor that contains the body index assigned to that processor after performing scatter level 1. 
  !mrb_body: Size of rank_body array 
  
    PetscInt :: rk, rki, i, rank_excess, cl
    PetscInt :: nbmax, nbmax_loc, nb_lim, nproc_body
    PetscInt, allocatable :: new_proc(:), rank_nb_dummy(:)
    PetscInt, allocatable :: xbm_ix(:), xbm_glob_ix(:)
    PetscInt :: it, body_no, ii
    Vec :: xbmvec_body, xbmvec
    VecScatter :: ctx_body
    PetscScalar :: xbm_max
    
    nb_lim=1000  !This is the limit of number of unknowns a proc should handle, required for scatter level 2 (only used when there is excess number of procs)
    
    body_color=MPI_UNDEFINED
    allocate(rank_nb(nproc), rank_nb_dummy(nproc), body_rank(mt+md,2))
    rank_nb=0
    
    IF (mt+md>nproc) then
      
      rk=0
      do i=mr+1, m_body
        
        if (rk<nproc) then
          rki=rk
        else
          rki=minloc(rank_nb,DIM=1)-1
        end if
        
        rk=rk+1
        body_rank(rk,1)=i
        body_rank(rk,2)=rki     !body number i is in rank rki
        rank_nb(rki+1)=2*nb_ind(i)+rank_nb(rki+1)
        
      end do
    
    ELSE  !(mt+md<=nproc)
      
      rk=0
      do i=mr+1, m_body
        rki=rk
        rk=rk+1
        body_rank(rk,1)=i
        body_rank(rk,2)=rki
        rank_nb(rki+1)=2*nb_ind(i)
      end do
      
      rank_nb_dummy=rank_nb
      rank_excess=nproc-(rki+1)
      cl=nproc   !make sure this number is always greater than the largest rank
      do while (rank_excess>0)
      
        nbmax=maxval(rank_nb_dummy,DIM=1)
        nbmax_loc=maxloc(rank_nb_dummy,DIM=1)
        
        if (nbmax>nb_lim) then
          nproc_body=ceiling(real(nbmax)/real(nb_lim))
          if (rank_excess<nproc_body-1) nproc_body=rank_excess+1
          allocate(new_proc(nproc_body-1))
          new_proc=(/(rki+i, i=1,nproc_body-1, 1)/)
          if (rank==body_rank(nbmax_loc,2) .or. any(new_proc==rank)) body_color=cl
          deallocate(new_proc)
          cl=cl+1
          rki=rki+nproc_body-1   !these many ranks have been considered till now
          rank_excess=nproc-(rki+1) !remaining procs
          rank_nb_dummy(nbmax_loc)=0  !zero out the element that had max number of unknonw, because we generated the appropriate colors
          
        else
          rank_excess=0  !if the maximum number of unknonws in any proc itself is smaller than the limit, then exit the while loop asap
        end if

      end do
    
    END IF
    
    if (body_color == MPI_UNDEFINED) body_color=rank !this line provides color to every processor to make comm_body
    call MPI_Comm_split(MPI_COMM_WORLD, body_color, 0, comm_body, ierr) !Define the comm_body for scatter 2
    
    do i=1,nproc
      if (rank==i-1) xbm=rank_nb(i)
    end do
    
    !xbm_ix is the indices for the sequential vector, xbm_glob_ix contains the indices of fvec which has to be scattered to the sequential vector
    allocate(xbm_ix(xbm), xbm_glob_ix(xbm))
    xbm_ix=(/(i, i=0,xbm-1, 1)/)
    call ISCreateGeneral(MPI_COMM_SELF, xbm, xbm_ix, PETSC_COPY_VALUES, xbm_is, ierr)
    
    it=0
    xbm_glob_ix=0
    do i=1,mt+md
      if (rank==body_rank(i,2)) then
        body_no=body_rank(i,1)  !body_number
        xbm_glob_ix(it+1:it+2*nb_ind(body_no))=(/(ii, ii=(nbb_ind(body_no)-2*nb_ind(body_no)+1)-1, nbb_ind(body_no)-1, 1)/)  !-1 in ii for C indexing
        it=it+2*nb_ind(body_no)
      end if
    end do
    call ISCreateGeneral(MPI_COMM_SELF, xbm, xbm_glob_ix, PETSC_COPY_VALUES, xbm_glob_is, ierr)
    
    mrb=0
    do i=1,mt+md
      if (rank==body_rank(i,2)) then
        mrb=mrb+1
        call grow_array(rank_body, body_rank(i,1), mrb)
        
      end if
    end do
    
    call MPI_COMM_SIZE(comm_body, nproc_comm_body_int4, ierr)  !Here we identify the number of procs in the respective subcommuniator
    call MPI_COMM_RANK(comm_body, rank_comm_body_int4, ierr)
    nproc_comm_body = nproc_comm_body_int4
    rank_comm_body = rank_comm_body_int4

    !Scatter the values of xbm to everyone else so that every proc knows the length of the mpi vector to be created on comm_body  (this is specifically required when one body is divided into many procs)
    call VecCreateMPI(comm_body, ione, nproc_comm_body, xbmvec, ierr)
    call VecScatterCreateToAll(xbmvec, ctx_body, xbmvec_body, ierr)
    
    call VecSetValue(xbmvec, rank_comm_body, real(xbm,8), INSERT_VALUES, ierr)
    call VecScatterBegin(ctx_body, xbmvec, xbmvec_body, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(ctx_body, xbmvec, xbmvec_body, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecMax(xbmvec_body, PETSC_NULL_INTEGER, xbm_max, ierr)
    xbm_body=int(xbm_max)
    
    call VecSetValue(xbmvec, rank_comm_body, real(mrb,8), INSERT_VALUES, ierr)
    call VecScatterBegin(ctx_body, xbmvec, xbmvec_body, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(ctx_body, xbmvec, xbmvec_body, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecMax(xbmvec_body, PETSC_NULL_INTEGER, xbm_max, ierr)
    mrb_body=int(xbm_max)
    
    if (mrb_body>0 .and. .not. allocated(rank_body)) then 
      allocate(rank_body(1))   !here we are allocating rank_body on the procs that will receive force vectors on the second level of scatter on comm_body. This is only done on those procs which has been identified to receive a part of 1 body
      rank_body=0
    end if
    
    if (mrb_body==1) then  !Now here we are assuming that the proc that that will scatter on the second level of scatter operation handles only 1 body and similarly the proc that is accepting the body is going to handle only 1 as well.
      call VecSetValue(xbmvec, rank_comm_body, real(rank_body(1),8), INSERT_VALUES, ierr)
      call VecScatterBegin(ctx_body, xbmvec, xbmvec_body, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(ctx_body, xbmvec, xbmvec_body, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecMax(xbmvec_body, PETSC_NULL_INTEGER, xbm_max, ierr)
      rank_body(1)=int(xbm_max)
    
    end if
    
    call VecDestroy(xbmvec, ierr)
    call VecScatterDestroy(ctx_body, ierr)
    call VecDestroy(xbmvec_body, ierr)
  
  end subroutine body_proc

!================================================================================

  subroutine list_bodies   
  !Identifies the number of bodies of torsional and deformable families to be handled by each proc 
  
    PetscInt :: i, i_bdy
    
    mt_body=0
    md_body=0
    do i=1,mrb_body
      i_bdy=rank_body(i)
      
      if (bdy(i_bdy)%fam==2) then  !Torsional body 
        mt_body=mt_body+1
      end if
      
      if (bdy(i_bdy)%fam==3) then  !Deformable body 
        md_body=md_body+1
      end if
      
    end do
  
  end subroutine list_bodies
    
!================================================================================  
  
    SUBROUTINE FindInVector(n,TF,npos,pos)
    
    !obtained from: https://stackoverflow.com/a/33331489/11800173
    ! Inlet variables
    INTEGER,INTENT(IN):: n      ! Dimension of logical vector
    LOGICAL,INTENT(IN):: TF(n)  ! Logical vector (True or False)
    ! Outlet variables
    INTEGER npos                ! number of "true" conditions
    INTEGER pos(n)              ! position of "true" conditions
    ! Internal variables
    INTEGER i                   ! counter
    INTEGER v(n)                ! vector of all positions

    pos = 0                     ! Initialize pos
    FORALL(i=1:n)   v(i) = i    ! Enumerate all positions
    npos  = COUNT(TF)           ! Count the elements of TF that are .True.
    pos(1:npos)= pack(v, TF)    ! With Pack function, verify position of true conditions

    END SUBROUTINE FindInVector
!================================================================================

  function find_min(dist, interp_points) result (nearest_neighbor)
  
      PetscInt :: interp_points, nearest_neighbor(interp_points)
      PetscScalar :: dist(nb_sd)
      PetscInt :: a1(1), i
      Logical, dimension(nb_sd) :: mk
      
      mk = .TRUE.
      do i=1,interp_points
          a1 = minloc(dist,mk)
          nearest_neighbor(i) = a1(1)
          mk(a1(1)) = .FALSE.      
      end do
  
  end function find_min
  
!===================================================================================  

!Unique sorting algorithm obtained from: https://stackoverflow.com/a/44242476/11800173

  Recursive Subroutine MergeSort(temp, Begin, Finish, list)
    ! 1st 3 arguments are input, 4th is output sorted list
    implicit none
    integer(kind=4),intent(inout) :: Begin,list(:),temp(:)
    integer(kind=4),intent(in) :: Finish
    integer(kind=4) :: Middle
    if (Finish-Begin<2) then    !if run size =1
       return                   !it is sorted
    else
       ! split longer runs into halves
       Middle = (Finish+Begin)/2
       ! recursively sort both halves from list into temp
       call MergeSort(list, Begin, Middle, temp)
       call MergeSort(list, Middle, Finish, temp)
       ! merge sorted runs from temp into list
       call Merge(temp, Begin, Middle, Finish, list)
     endif
  End Subroutine MergeSort

  Subroutine Merge(list, Begin, Middle, Finish, temp)
    implicit none
    integer(kind=4),intent(inout) :: list(:),temp(:)
    integer(kind=4),intent(in) ::Begin,Middle,Finish
    integer(kind=4)    :: kx,ky,kz
    ky=Begin
    kz=Middle
    !! While there are elements in the left or right runs...
    do kx=Begin,Finish-1
       !! If left run head exists and is <= existing right run head.
       if (ky.lt.Middle.and.(kz.ge.Finish.or.list(ky).le.list(kz))) then
          temp(kx)=list(ky)
          ky=ky+1
       else
          temp(kx)=list(kz)
          kz = kz + 1
       end if
    end do

  End Subroutine Merge

  Function Unique(list)
    !! usage sortedlist=Unique(list)
    implicit none
    integer(kind=4) :: strt,fin,N
    integer(kind=4), intent(inout) :: list(:)
    integer(kind=4), allocatable  :: unique(:),work(:)
    logical,allocatable :: mask(:)
    ! sort
    work=list;strt=1;N=size(list);fin=N+1
    call MergeSort(work,strt,fin,list) 
    ! cull duplicate indices
    allocate(mask(N));
    mask=.false.
    mask(1:N-1)=list(1:N-1)==list(2:N)
    unique=pack(list,.not.mask)

  End Function Unique
  
!================================================================================

  subroutine destroy_grid
  
     call DMDestroy(das, ierr)
     call DMDestroy(das_coarsify, ierr)
     call DMDestroy(dau, ierr)
     call DMDestroy(dav, ierr)
     
     call MPI_Comm_free(comm_body, ierr)
     
     call VecDestroy(fvec, ierr)
     call VecDestroy(xbvec, ierr)
     call VecDestroy(ybvec, ierr)
     
     if (sub_domain) then
         call VecDestroy(xsdvec, ierr)
         call VecDestroy(ysdvec, ierr)
     end if
     
     call ISDestroy(interfacex_is, ierr)
     call ISDestroy(interfacey_is, ierr)
     call ISDestroy(reg_ixs, ierr)
     call ISDestroy(reg_iys, ierr)
     call ISDestroy(bcis, ierr)
     
     if (sub_domain) then
         call ISDestroy(interfacex_sd_is_sparse, ierr)
         call ISDestroy(interfacey_sd_is_sparse, ierr)
         call ISDestroy(reg_ixs_sd_sparse, ierr)
         call ISDestroy(reg_iys_sd_sparse, ierr)
         
         call VecDestroy(interp_indices_global_vec, ierr)
         call VecDestroy(interp_indices_vec, ierr)
         call VecScatterDestroy(ctx_indices, ierr)
         
         call VecDestroy(interp_weights_global_vec, ierr)
         call VecDestroy(interp_weights_vec, ierr)
         call VecScatterDestroy(ctx_weights, ierr)
         
         deallocate(interp_indices_global_sparsemat, interp_weights_global_sparsemat)
         
         call VecDestroy(intermediatevec, ierr)
         call VecDestroy(intermediatelocalvec, ierr)         
         call ISDestroy(interp_indices_local_is, ierr)
         call VecScatterDestroy(ctx_intermediate, ierr)
     end if
     
  end subroutine destroy_grid  

!================================================================================

END MODULE grid
