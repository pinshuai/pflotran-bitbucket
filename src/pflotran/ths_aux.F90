module THS_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscBool, public :: ths_analytical_derivatives = PETSC_FALSE

  ! Toggle individual equations
  PetscBool, public :: ths_isothermal = PETSC_FALSE
  PetscBool, public :: ths_no_solvent = PETSC_FALSE
  PetscBool, public :: ths_no_solute = PETSC_FALSE

  type, public :: ths_component_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: dof
    type(ths_component_type), pointer :: next
  end type ths_component_type

  type, public :: ths_auxvar_type
    PetscInt :: istate_store(2)
    PetscInt :: phase ! (liquid = 1, gas = 2, solid = 3)
    PetscReal :: pres(2)   ! (phase pressure, capillary pressure)
    PetscReal :: sat    ! (phase)
    PetscReal :: den    ! (phase) kmol/m^3 phase
    PetscReal :: den_kg ! (phase) kg/m^3 phase
    PetscReal :: temp
    PetscReal, pointer :: xmol(:) ! (icomp)
    PetscReal :: H ! MJ/kmol
    PetscReal :: U ! MJ/kmol
    PetscReal :: mobility ! relative perm / kinematic viscosity
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: pert
    type(ths_derivative_auxvar_type), pointer :: d
  end type ths_auxvar_type

  type, public :: ths_energy_derivative_type
    PetscReal :: Ha_pa
    PetscReal :: Ha_T
    PetscReal :: Ua_T
    PetscReal :: Ua_pa
  end type ths_energy_derivative_type

  type, public :: ths_mass_derivative_type
    type(ths_component_type) :: component
    PetscReal, pointer :: xmol_P(:,:)
    PetscReal, pointer :: xmol_T(:,:)
  end type ths_mass_derivative_type

  type, public :: ths_derivative_auxvar_type
    type(ths_energy_derivative_type) :: energy_derivative
    type(ths_mass_derivative_type) :: mass_derivative
  end type ths_derivative_auxvar_type

  type, public :: ths_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:) ! (icomp)
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
  end type ths_parameter_type

  type, public :: ths_type
    PetscInt :: n_inactive_rows
    PetscInt, pointer :: inactive_rows_local(:), inactive_rows_local_ghosted(:)
    PetscInt, pointer :: row_zeroing_array(:)

    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(ths_component_type), pointer :: component
    type(ths_parameter_type), pointer :: ths_parameter
    type(ths_auxvar_type), pointer :: auxvars(:,:)
    type(ths_auxvar_type), pointer :: auxvars_bc(:)
    type(ths_auxvar_type), pointer :: auxvars_ss(:)
  end type ths_type

  public :: THSAuxCreate, &
            THSAuxDestroy, &
            THSAuxVarCompute, &
            THSAuxVarInit, &
            THSAuxVarCopy, &
            THSAuxVarDestroy, &
            THSAuxVarStrip, &
            THSAuxVarPerturb, &
            THSPrintAuxVars, &
            THSOutputAuxVars

contains

! ************************************************************************** !

function THSAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(ths_type), pointer :: THSAuxCreate

  type(ths_type), pointer :: aux
  type(ths_component_type), pointer :: component

  allocate(aux)
  component => aux%component
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0

!   if (option%ncomp > 0) then
!     do i = 1, option%ncomp
!       allocate(component)
!       component => component%next
!     enddo
!   enddo

  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_inactive_rows = 0
  nullify(aux%inactive_rows_local)
  nullify(aux%inactive_rows_local_ghosted)
  nullify(aux%row_zeroing_array)

  allocate(aux%ths_parameter)
!   allocate(aux%ths_parameter%diffusion_coefficient(option%ncomp))

  aux%ths_parameter%newton_inf_scaled_res_tol = 1.d-50
  aux%ths_parameter%check_post_converged = PETSC_FALSE

  THSAuxCreate => aux

end function THSAuxCreate

! ************************************************************************** !

subroutine THSAuxVarInit(auxvar,allocate_derivative,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  use Option_module

  implicit none

  type(ths_auxvar_type) :: auxvar
  PetscBool :: allocate_derivative
  type(option_type) :: option

  auxvar%phase = 0
  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%temp = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%pert = 0.d0


end subroutine THSAuxVarInit

! ************************************************************************** !

subroutine THSAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Michael Nole
  ! Date: 02/07/19
  ! 

  use Option_module

  implicit none
  
  type(ths_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate_store = auxvar%istate_store
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%xmol = auxvar%xmol
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%mobility = auxvar%mobility
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%pert = auxvar%pert

end subroutine THSAuxVarCopy

! ************************************************************************** !

subroutine THSAuxVarCompute()

end subroutine THSAuxVarCompute

! ************************************************************************** !

subroutine THSAuxVarPerturb()

end subroutine THSAuxVarPerturb

! ************************************************************************** !

subroutine THSPrintAuxVars()

end subroutine THSPrintAuxVars

! ************************************************************************** !

subroutine THSOutputAuxVars()

end subroutine THSOutputAuxVars

! ************************************************************************** !

subroutine THSAuxvarStrip()

end subroutine THSAuxVarStrip

! ************************************************************************** !

subroutine THSAuxDestroy()

end subroutine THSAuxDestroy

! ************************************************************************** !

subroutine THSAuxVarDestroy()

end subroutine THSAuxVarDestroy
! ************************************************************************** !

end module THS_Aux_module
