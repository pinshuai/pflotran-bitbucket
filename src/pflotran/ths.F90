module THS_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use THS_Aux_module
  use THS_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

!   private

  public :: THSSetup, &
            THSUpdateSolution, &
            THSTimeCut,&
            THSUpdateAuxVars, &
            THSUpdateFixedAccum, &
            THSComputeMassBalance, &
            THSResidual, &
            THSJacobian, &
            THSGetTecplotHeader, &
            THSSetPlotVariables, &
            THSMapBCAuxVarsToGlobal, &
            THSDestroy
    
contains
! ************************************************************************** !
subroutine THSSetup(realization)
  !
  !Author: Michael Nole
  !Date: 02/21/19
  !
  
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Fluid_module
  use Material_Aux_class
  use Output_Aux_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
end subroutine THSSetup

! ************************************************************************** !
subroutine THSUpdateSolution(realization)
  !
  !Author: Michael Nole
  !Date: 02/21/19
  !
  
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
end subroutine THSUpdateSolution

! ************************************************************************** !
subroutine THSTimeCut(realization)
  !
  !Author: Michael Nole
  !Date: 02/21/19
  !
  
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

end subroutine THSTimeCut
! ************************************************************************** !

subroutine THSUpdateAuxVars(realization)
  !
  !Author: Michael Nole
  !Date: 02/21/19
  !
  
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_class
  use EOS_Water_module
  use Saturation_Function_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  

end subroutine THSUpdateAuxVars

! ************************************************************************** !

subroutine THSUpdateFixedAccum(realization)
  !
  !Author: Michael Nole
  !Date: 02/21/19
  !
  
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none
  
  type(realization_subsurface_type) :: realization
  
end subroutine THSUpdateFixedAccum

! ************************************************************************** !

subroutine THSComputeMassBalance(realization, mass_balance)
  !
  !Author: Michael Nole
  !Date: 02/21/19
  !
  
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)

end subroutine THSComputeMassBalance

! ************************************************************************** !

subroutine THSResidual(snes,xx,r,realization,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !
  
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module

  use Connection_module
  use Grid_module
  use Coupler_module  
  use Debug_module
  use Material_Aux_class
  use Upwind_Direction_module
  
  implicit none
  
  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
end subroutine THSResidual
! ************************************************************************** !

subroutine THSJacobian(snes,xx,A,B,realization,ierr)
  !
  ! This routine performs a Jacobian calculation
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class
  use Upwind_Direction_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A,B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

end subroutine THSJacobian

! ************************************************************************** !

subroutine THSMapBCAuxVarsToGlobal(realization)
  !
  !Author: Michael Nole
  !Date: 02/21/19
  !
  
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  type(realization_subsurface_type) :: realization
  
end subroutine THSMapBCAuxVarsToGlobal

! ************************************************************************** !

function THSGetTecplotHeader(realization,icolumn)


  use Realization_Subsurface_class
  use Option_module
  use Field_module
    
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: THSGetTecplotHeader
  type(realization_subsurface_type) :: realization
  PetscInt :: icolumn
  
end function THSGetTecplotHeader

! ************************************************************************** !

subroutine THSSetPlotVariables(realization,list)

  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

end subroutine THSSetPlotVariables

! ************************************************************************** !

subroutine THSDestroy(realization)

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization

end subroutine THSDestroy

! ************************************************************************** !

end module THS_module
