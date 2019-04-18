module THS_Common_module

#include "petsc/finclude/petscsys.h"

  use petscsys
  use THS_Aux_module
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  public :: THSAccumulation, &
            THSFlux, &
            THSBCFlux, &
            THSSrcSink
            
contains

! ************************************************************************** !

subroutine THSAccumulation(ths_auxvar,global_auxvar,option, &
                           Res, Jac, analytical_derivatives, debug_cell)
                           
  !
  ! Computes the non-fixed portion of the accumulation term for the residual.
  ! Author: Michael Nole
  ! Date: 03/12/19
  !
  
  use Option_module
  
  implicit none
  
  type(ths_auxvar_type) :: ths_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof, option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_cell
  
  PetscInt :: icomp, iphase, energy_id
  
  Res = 0.d0
  do iphase = 1,option%nphase
    do icomp = 1, option%nflowspec
!       Res(icomp) = Res(icomp) + ths_auxvar%sat(iphase) * &
!                                 ths_auxvar%den(iphase) * &
!                                 ths_auxvar%xmol(icomp,iphase)
    enddo
!     Res(energy_id) = Res(energy_id) + ths_auxvar%sat(iphase) * &
!                                       ths_auxvar%den(iphase) * &
!                                       ths_auxvar%U(iphase)
  enddo
  
!   if (analytical_derivatives) then
!   endif
  
end subroutine THSAccumulation

! ************************************************************************** !

subroutine THSFlux(ths_auxvar_up, global_auxvar_up, material_auxvar_up, &
                   ths_auxvar_dn, global_auxvar_dn, material_auxvar_dn, &
                   area, dist, upwind_direction, &
                   option, Res, Jup, Jdn, analytical_derivatives, &
                   debug_connection)
  !
  ! Computes internal flux terms for the residual
  ! Author: Michael Nole
  ! Date: 03/12/19
  !
  
  use Option_module
  use Material_Aux_class
  use Upwind_Direction_module
  
  implicit none

  type(ths_auxvar_type) :: ths_auxvar_up, ths_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area, dist(-1:3)
  PetscInt :: upwind_direction(option%nphase)
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_connection

    

end subroutine THSFlux

! ************************************************************************** !

subroutine THSBCFlux()

end subroutine THSBCFlux

! ************************************************************************** !

subroutine THSSrcSink()

end subroutine THSSrcSink

! ************************************************************************** !

end module THS_Common_module
