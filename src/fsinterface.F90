module fsinterface

  USE petsc
#include "petsc/finclude/petsc.h"
  USE parameters
  USE grid
  implicit none


contains

!===================================================================================

  subroutine reg(h0, hxvec, hyvec)
  !Regularize h0 to nearby flow grid points resulting in hxvec, hyvec. Communication of hxvec and hyvec are performed outside this subroutine
  
  PetscScalar :: h0(xfi:xfe)
  Vec :: hxvec, hyvec
  PetscInt :: k, l, p, next, iterx, itery
  PetscScalar :: x(interfacex), y(interfacey)
  
  iterx=0
  itery=0
  do k=xfi,xfe

    next=0
    do l=-support2, support2
      do p=-support2, support2
        next=next+1
        
        if (codeb(k,3)==1) then
          iterx=iterx+1
          x(iterx)=weight(next, k)*h0(k)
        else if (codeb(k,3)==2) then
          itery=itery+1
          y(itery)=weight(next, k)*h0(k)
        end if
        
      end do
    end do
  
  end do
  
  call VecSetValues(hxvec, interfacex, reg_ix, x, ADD_VALUES, ierr)
  call VecSetValues(hyvec, interfacey, reg_iy, y, ADD_VALUES, ierr)  
  
  end subroutine reg

!===================================================================================

  subroutine regT(hx, hy, h0)
  !Performs interpolation i.e. action of E. Interpolates velocity field hx, hy onto the body h0
  !Before regT is performed, appropriate velocity field is scatter to obtain hx, hy since the partitioning of flow domain is different from the partitioning of the structure
  
    PetscScalar :: hx(interfacex), hy(interfacey), h0(xfi:xfe)
    PetscInt :: k, l, p, next, iterx, itery
    
    iterx=0
    itery=0
    do k=xfi,xfe

      next=0
      do l=-support2, support2
        do p=-support2, support2
          next=next+1
        
          if (codeb(k,3)==1) then
            iterx=iterx+1
            h0(k)=h0(k)+weight(next, k)*hx(iterx)
          else if (codeb(k,3)==2) then
            itery=itery+1
            h0(k)=h0(k)+weight(next, k)*hy(itery)
          end if
        
        end do
      end do
  
    end do
  
  end subroutine regT
  
!===================================================================================

  subroutine reg_sd_sparse(h0, hxsdvec, hysdvec, sd_list_size, sd_list)
  !same as reg, but on the sub-domain performed sparsely
  
  PetscScalar :: h0(xsdi:xsde)
  Vec :: hxsdvec, hysdvec
  integer, intent(in) :: sd_list_size
  PetscInt, intent(in) :: sd_list(sd_list_size)
  PetscInt :: k, l, p, next, iterx, itery, kk
  PetscScalar :: x(interfacex_sd_sparse), y(interfacey_sd_sparse)
  
  iterx=0
  itery=0
  do kk=1,sd_list_size !not looped over xsdi:xsde for efficiency
    k = sd_list(kk)
    if (k.ge.xsdi .and. k.le.xsde) then
    next=0
    do l=-support2, support2
      do p=-support2, support2
        next=next+1
        
        if (codesd(k,1)==1) then
          iterx=iterx+1
          x(iterx)=weight_sd(next, k)*h0(k)
        else if (codesd(k,1)==2) then
          itery=itery+1
          y(itery)=weight_sd(next, k)*h0(k)
        end if
        
      end do
    end do
    end if
  end do

  call VecSetValues(hxsdvec, interfacex_sd_sparse, reg_ix_sd_sparse, x, ADD_VALUES, ierr)
  call VecSetValues(hysdvec, interfacey_sd_sparse, reg_iy_sd_sparse, y, ADD_VALUES, ierr) 
  
  call VecAssemblyBegin(hxsdvec, ierr)
  call VecAssemblyBegin(hysdvec, ierr)  
  call VecAssemblyEnd(hxsdvec, ierr)
  call VecAssemblyEnd(hysdvec, ierr) 

  
  end subroutine reg_sd_sparse  
  
!===================================================================================

  subroutine regT_sd_sparse(hx, hy, h0, sd_list_size, sd_list)
  !same as regT, but on the sub-domain performed sparsely
  
    PetscScalar :: hx(interfacex_sd_sparse), hy(interfacey_sd_sparse), h0(xsdi:xsde)
    integer, intent(in) :: sd_list_size
    PetscInt, intent(in) :: sd_list(sd_list_size)
    PetscInt :: k, l, p, next, iterx, itery, kk
    
    iterx=0
    itery=0
    do kk=1,sd_list_size !not looped over xsdi:xsde for efficiency
      k = sd_list(kk)
      if (k.ge.xsdi .and. k.le.xsde) then
      next=0
      do l=-support2, support2
        do p=-support2, support2
          next=next+1
        
          if (codesd(k,1)==1) then
            iterx=iterx+1
            h0(k)=h0(k)+weight_sd(next, k)*hx(iterx)
          else if (codesd(k,1)==2) then
            itery=itery+1
            h0(k)=h0(k)+weight_sd(next, k)*hy(itery)
          end if
        
        end do
      end do
      end if
    end do
  
  end subroutine regT_sd_sparse
    
    
!===================================================================================

end module fsinterface
