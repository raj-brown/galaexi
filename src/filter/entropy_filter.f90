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
  USE MOD_GPU, ONLY: GPU_DefaultStream
  USE MOD_Mesh_Vars, ONLY: nElems
  USE MOD_Interpolation_Vars, ONLY: d_Vdm_Leg, d_sVdm_Leg
  USE MOD_ChangeBasisByDim, ONLY: ChangeBasisVolume_GPU
  IMPLICIT NONE

  REAL, DEVICE, INTENT(INOUT) :: U_in(PP_nVar, 0:PP_N, 0:PP_N,0:PP_NZ, nElems)
  REAL, INTENT(IN) :: tStage
  INTEGER, INTENT(IN) :: iStage
  REAL, DEVICE, ALLOCATABLE :: d_U_modal(:,:,:,:,:)
  REAL, DEVICE, ALLOCATABLE :: d_U_rec(:,:,:,:,:)

  REAL, ALLOCATABLE :: U_original(:,:,:,:,:)
  REAL, ALLOCATABLE :: U_reconstructed(:,:,:,:,:)
  REAL :: error_rel, norm_ref, error_abs
  INTEGER :: nThreads, shmem_size
  INTEGER :: bytes_per_real
  INTEGER :: ierr
  IF (iStage == 1) THEN
   RETURN

  ALLOCATE(d_U_modal(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_NZ, nElems))
  ALLOCATE(d_U_rec(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_NZ, nElems))

  ALLOCATE(U_original(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_NZ, nElems))
  ALLOCATE(U_reconstructed(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_NZ, nElems))
  U_original = U_in
  nThreads = (PP_N+1)**2 * (PP_NZ+1)
  bytes_per_real = STORAGE_SIZE(0.0)/8
  shmem_size = 2*nThreads * bytes_per_real * PP_nVar

  CALL ChangeBasisVolume_GPU<<<nElems, nThreads, shmem_size>>>(&
  PP_nVar, PP_N, PP_N, nElems, d_sVdm_Leg, U_in, d_U_modal)

  CALL ChangeBasisVolume_GPU<<<nElems, nThreads, shmem_size>>>(&
  PP_nVar, PP_N, PP_N, nElems, d_Vdm_Leg, d_U_modal, d_U_rec)

  U_reconstructed = d_U_rec

  error_abs = MAXVAL(ABS(U_original - U_reconstructed))
  norm_ref = SQRT(SUM(U_original**2))

  error_rel = SQRT(SUM((U_original - U_reconstructed)**2)) / MAX(TINY, norm_ref)

  WRITE(*,'(A,I4,A,ES14.6)') &
      'EntropyFilter reconstruction: stage=',iStage, &
      ' relative_L2=',error_rel

    WRITE(*,'(A,ES14.6)') &
      'EntropyFilter reconstruction: max_abs=',error_abs

    !-------------------------------------------------------------
    ! Release temporary memory
    !-------------------------------------------------------------
    DEALLOCATE(d_U_modal,d_U_rec)
    DEALLOCATE(U_original,U_reconstructed)


  END SUBROUTINE ApplyEntropyFilter


  SUBROUTINE FinalizeEntropyFilter()
    IMPLICIT NONE
END SUBROUTINE FinalizeEntropyFilter

  END MODULE MOD_EntropyFilter
