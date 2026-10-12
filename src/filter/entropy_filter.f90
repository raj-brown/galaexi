#include "flexi.h"

MODULE MOD_EntropyFilter
  USE CUDAFOR
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: InitEntropyFilter
  PUBLIC :: ApplyEntropyFilter
  PUBLIC :: FinalizeEntropyFilter

  CONTAINS

  SUBROUTINE InitEntropyFilter()
  IMPLICIT NONE

  END SUBROUTINE InitEntropyFilter

  SUBROUTINE ApplyEntropyFilter(U_in, tStage, iStage)
  USE MOD_PreProc
  USE MOD_Mesh_Vars, ONLY: nElems

  IMPLICIT NONE

  REAL, DEVICE, INTENT(INOUT) :: U_in(PP_nVar, 0:PP_N, 0:PP_N,0:PP_NZ, nElems)
  REAL, INTENT(IN) :: tStage
  INTEGER, INTENT(IN) :: iStage

  END SUBROUTINE ApplyEntropyFilter


  SUBROUTINE FinalizeEntropyFilter()
    IMPLICIT NONE
END SUBROUTINE FinalizeEntropyFilter

  END MODULE MOD_EntropyFilter
