module PM_THS_class
        
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use THS_module
  use THS_Aux_module
  use THS_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_subsurface_flow_type) :: pm_ths_type

  contains
    procedure, public :: Read => PMTHSRead
    procedure, public :: InitializeRun => PMTHSInitializeRun
    procedure, public :: InitializeTimestep => PMTHSInitializeTimestep
    procedure, public :: Residual => PMTHSResidual
    procedure, public :: Jacobian => PMTHSJacobian
    procedure, public :: UpdateTimestep => PMTHSUpdateTimestep
    procedure, public :: PreSolve => PMTHSPreSolve
    procedure, public :: PostSolve => PMTHSPostSolve
    procedure, public :: CheckUpdatePre => PMTHSCheckUpdatePre
    procedure, public :: CheckUpdatepost => PMTHSCheckUpdatePost
    procedure, public :: CheckConvergence => PMTHSCheckConvergence
    procedure, public :: TimeCut => PMTHSTimeCut
    procedure, public :: UpdateSolution => PMTHSUpdateSolution
    procedure, public :: UpdateAuxVars => PMTHSUpdateAuxVars
    procedure, public :: ComputeMassBalance => PMTHSComputeMassBalance
    procedure, public :: InputRecord => PMTHSInputRecord
    procedure, public :: CheckpointBinary => PMTHSCheckpointBinary
    procedure, public :: RestartBinary => PMTHSRestartBinary
    procedure, public :: Destroy => PMTHSDestroy
  end type pm_ths_type

  public :: PMTHSCreate

contains

! ************************************************************************** !

function PMTHSCreate()
 !
 ! Creates THS Process model shell
 !
 ! Author: Michael Nole
 ! Date: 02/06/19
 !

 use Upwind_Direction_module

 implicit none

 class(pm_ths_type), pointer :: PMTHSCreate
 
 class(pm_ths_type), pointer :: ths_pm
 
 allocate(ths_pm)
 
 call PMSubsurfaceFlowCreate(ths_pm)
 ths_pm%name = 'Thermo-Hydro-Solute'
 ths_pm%header = 'THERMO_HYDRO_SOLUTE'

 fix_upwind_direction = PETSC_FALSE

 PMTHSCreate => ths_pm

end function PMTHSCreate

! ************************************************************************** !

subroutine PMTHSRead(this, input)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !
  use THS_module
  use THS_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_ths_type) :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscReal :: tempreal

  option => this%option

  error_string = "THS Options"

  input%ierr = 0

  do
    call InputReadPflotranString(input, option)

    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case (keyword)
      case('ISOTHERMAL')
        ! Turns off energy conservation equation
        ths_isothermal = PETSC_TRUE
      case('NO_SOLUTE')
        ! Turns off solute mass conservation equation
        ths_no_solute = PETSC_TRUE
      case('NO_SOLVENT','NO_WATER')
        ! Turns off water mass conservation equation
        ths_no_solvent = PETSC_TRUE
    end select

  enddo

end subroutine PMTHSRead

! ************************************************************************** !

recursive subroutine PMTHSInitializeRun(this)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !
  
  use Realization_Base_class

  implicit none

  class(pm_ths_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr

  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMTHSInitializeRun

! ************************************************************************** !

subroutine PMTHSInitializeTimestep(this)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  use Global_module
 
  implicit none
  
  class(pm_ths_type) :: this
  
  PetscInt :: ths_newton_iteration_number
  PetscInt :: ths_ni_count
  PetscBool :: update_upwind_direction

  call PMSubsurfaceFlowInitializeTimestepA(this)

  ths_newton_iteration_number = 0
  update_upwind_direction = PETSC_TRUE

  call PMTHSUpdateFixedAccum(this%realization)

  ths_ni_count = 0

  call PMSubsurfaceFlowInitializeTimestepB(this)
  

end subroutine PMTHSInitializeTimestep

! ************************************************************************** !

subroutine PMTHSPreSolve(this)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  implicit none

  class(pm_ths_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMTHSPreSolve

! ************************************************************************** !

subroutine PMTHSPostSolve(this)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  implicit none

  class(pm_ths_type) :: this

end subroutine PMTHSPostSolve

! ************************************************************************** !

subroutine PMTHSUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                               num_newton_iterations, tfac, &
                               time_step_max_growth_factor)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !
  
  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION

  implicit none

  class(pm_ths_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ut, ux, umin
  PetscReal :: dtt
  type(field_type), pointer :: field

  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    up = this%pressure_change_governor/(this%max_pressure_change+0.1d0)
    ut = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
    ux = this%xmol_change_governor/(this%max_xmol_change+1.d-5)
    umin = min(up,ut,ux)
  endif
  ifac = max(min(num_newton_iterations,size(tfac)),1)
  dtt = fac * dt * (1.d0 + umin)
  dtt = min(time_step_max_growth_factor*dt,dtt)
  dt = min(dtt,tfac(ifac)*dt,dt_max)
  dt = max(dt,dt_min)

  !MAN: CFL limiter might need beefing up
  if (Initialized(this%cfl_governor)) then
    field => this%realization%field
    call RealizationGetVariable(this%realization,field%work, &
                                LIQUID_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               LIQUID_SATURATION,TIME_NULL)
    call RealizationGetVariable(this%realization,field%work, &
                                GAS_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               GAS_SATURATION,TIME_NULL)
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  endif

end subroutine PMTHSUpdateTimestep

! ************************************************************************** !

subroutine PMTHSResidual(this,snes,xx,r,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  implicit none
  class(pm_ths_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call THSResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTHSResidual

! ************************************************************************** !

subroutine PMTHSJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  ! 

  implicit none

  class(pm_ths_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call THSJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMTHSJacobian

! ************************************************************************** !

subroutine PMTHSCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  ! 

  use Grid_module
  use Option_module
  use Patch_module
  use THS_Aux_module
  use Global_Aux_module

  implicit none

  class(pm_THS_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

end subroutine PMTHSCheckUpdatePre

! ************************************************************************** !

subroutine PMTHSCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  ! 
  use THS_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module

  implicit none

  class(pm_ths_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

end subroutine PMTHSCheckUpdatePost

! ************************************************************************** !

subroutine PMTHSCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                     reason,ierr)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !
  
  use Convergence_module
  use THS_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use String_module

  implicit none

  class(pm_ths_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

end subroutine PMTHSCheckConvergence

! ************************************************************************** !
subroutine PMTHSTimeCut(this)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !
  
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module
  
  implicit none
  
  class(pm_ths_type) :: this
  
end subroutine PMTHSTimeCut

! ************************************************************************** !

subroutine PMTHSUpdateSolution(this)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  !
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module

  implicit none
  
  class(pm_ths_type) :: this
  
end subroutine PMTHSUpdateSolution


! ************************************************************************** !

subroutine PMTHSUpdateAuxVars(this)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19

  implicit none

  class(pm_ths_type) :: this

!   call THSUpdateAuxVars(this%realization,PETSC_FALSE)

end subroutine PMTHSUpdateAuxVars

! ************************************************************************** !

subroutine PMTHSComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  ! 

  implicit none

  class(pm_ths_type) :: this
  PetscReal :: mass_balance_array(:)

  call THSComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMTHSComputeMassBalance

! ************************************************************************** !

subroutine PMTHSInputRecord(this)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  ! 

  implicit none

  class(pm_ths_type) :: this

end subroutine PMTHSInputRecord

subroutine PMTHSCheckpointBinary(this,viewer)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_ths_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMTHSCheckpointBinary

! ************************************************************************** !

subroutine PMTHSRestartBinary(this,viewer)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  use Checkpoint_module
  use Global_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_ths_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMTHSRestartBinary

! ************************************************************************** !

subroutine PMTHSDestroy(this)
  ! 
  ! Author: Michael Nole
  ! Date: 02/06/19
  ! 

  implicit none

  class(pm_ths_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call THSDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMTHSDestroy

! ************************************************************************** !

subroutine PMTHSUpdateFixedAccum(realization)
  !
  ! Author: Michael Nole
  ! Date: 02/06/19
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none

  type(realization_subsurface_type) :: realization

end subroutine PMTHSUpdateFixedAccum


end module PM_THS_class
