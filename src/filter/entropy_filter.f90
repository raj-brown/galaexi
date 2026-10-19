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
  USE MOD_EOS_Vars, ONLY: Kappa
  USE MOD_Mesh_Vars, ONLY: nElems
  USE MOD_Interpolation_Vars, ONLY: d_Vdm_Leg, d_sVdm_Leg
  USE MOD_ChangeBasisByDim, ONLY: ChangeBasisVolume_GPU
  IMPLICIT NONE

  REAL, DEVICE, INTENT(INOUT) :: U_in(PP_nVar, 0:PP_N, 0:PP_N,0:PP_NZ, nElems)
  REAL, INTENT(IN) :: tStage
  INTEGER, INTENT(IN) :: iStage
  REAL, DEVICE, ALLOCATABLE :: d_U_modal(:,:,:,:,:)
  REAL, DEVICE, ALLOCATABLE :: d_U_rec(:,:,:,:,:)
  REAL, DEVICE, ALLOCATABLE :: d_rhoMin(:), d_pMin(:), d_sMIN(:)
  REAL, ALLOCATABLE :: h_rhoMin(:), h_pMin(:), h_sMIN(:)

  REAL, ALLOCATABLE :: U_original(:,:,:,:,:)
  REAL, ALLOCATABLE :: U_reconstructed(:,:,:,:,:)
  REAL :: error_rel, norm_ref, error_abs
  INTEGER :: nThreads, shmem_size
  INTEGER :: bytes_per_real
  INTEGER :: ierr
  IF (iStage /= 1) RETURN
   ALLOCATE(d_rhoMin(nElems), d_pMin(nElems), d_sMIN(nElems))
   ALLOCATE(h_rhoMin(nElems), h_pMin(nElems), h_sMIN(nElems))
   CALL EntropyMin_Kernel<<<(nElems + 255)/256, 256>>>(nElems, PP_N, U_in, Kappa, d_rhoMin, d_pMin, d_sMIN)
   h_rhoMin = d_rhoMin
   h_pMin = d_pMin
   h_sMIN = d_sMIN
   WRITE(*,'(A,ES14.6)') &
      'EntropyFilter: min_rho=',MINVAL(h_rhoMin)
    WRITE(*,'(A,ES14.6)') &
      'EntropyFilter: min_p=',MINVAL(h_pMin)
    WRITE(*,'(A,ES14.6)') &
      'EntropyFilter: min_s=',MINVAL(h_sMIN)
  DEALLOCATE(d_rhoMin, d_pMin, d_sMIN)
  DEALLOCATE(h_rhoMin, h_pMin, h_sMIN)
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
  CALL ModalDamping_Kernel<<<(nElems*nThreads + 255)/256, 256>>>(&
  nElems,PP_N, d_U_modal, 0.05)

  CALL ChangeBasisVolume_GPU<<<nElems, nThreads, shmem_size>>>(&
  PP_nVar, PP_N, PP_N, nElems, d_Vdm_Leg, d_U_modal, d_U_rec)

  U_reconstructed = d_U_rec

  error_abs = MAXVAL(ABS(U_original - U_reconstructed))
  norm_ref = SQRT(SUM(U_original**2))

  error_rel = SQRT(SUM((U_original - U_reconstructed)**2)) / MAX(1E-30, norm_ref)

  WRITE(*,'(A,I4,A,ES14.6)') &
      'EntropyFilter modal damping: stage=',iStage, &
      ' relative_L2=',error_rel

    WRITE(*,'(A,ES14.6)') &
      'EntropyFilter modal damping: max_abs=',error_abs

    !-------------------------------------------------------------
    ! Release temporary memory
    !-------------------------------------------------------------
    DEALLOCATE(d_U_modal,d_U_rec)
    DEALLOCATE(U_original,U_reconstructed)


  END SUBROUTINE ApplyEntropyFilter


  ATTRIBUTES(GLOBAL) SUBROUTINE ModalDamping_Kernel(nElems, Nloc, U_modal, zeta)
  USE MOD_PreProc
  IMPLICIT NONE
  INTEGER, VALUE, INTENT(IN)  :: nElems, Nloc
  REAL, DEVICE, INTENT(INOUT) :: U_modal(5, 0:Nloc, 0:Nloc, 0:Nloc, nElems)
  REAL, VALUE, INTENT(IN)     :: zeta
  INTEGER :: tid, nodeID, iElem
  INTEGER :: a,b,c,rest
  INTEGER :: nNodes
  REAL :: deg2, sigma

  nNodes = (Nloc+1)**2 * (Nloc+1)
  tid = (blockIdx%x - 1)* blockDim%x + threadIdx%x
  IF (tid > nElems*nNodes) RETURN
  iElem = (tid -1) / nNodes + 1
  nodeID = MOD(tid-1, nNodes)

  a = MOD(nodeID, Nloc+1)
  rest = nodeID / (Nloc+1)
  b = MOD(rest, Nloc+1)
  rest = rest / (Nloc+1)
  c = MOD(rest, Nloc+1)
  deg2 = REAL(MAX(a,b,c)**2)
  sigma = EXP(-zeta*deg2)
  U_modal(:,a,b,c,iElem) = sigma * U_modal(:,a,b,c,iElem)

  END SUBROUTINE ModalDamping_Kernel


ATTRIBUTES(GLOBAL) SUBROUTINE EntropyMin_Kernel(nElems, Nloc, U, gammaGas, rhoMin, pMin, sMIN)

IMPLICIT NONE

INTEGER, VALUE, INTENT(IN)  :: nElems, Nloc
REAL, VALUE, INTENT(IN)     :: gammaGas
REAL, DEVICE, INTENT(IN)     :: U(5, 0:Nloc, 0:Nloc, 0:Nloc, nElems)
REAL, DEVICE, INTENT(OUT)     :: rhoMin(nElems), pMin(nElems), sMIN(nElems)
INTEGER :: iElem, i,j,k
REAL :: rho, mom2, pressure, entropy
REAL :: rhoLo, pLo, sLo
iElem = (blockIdx%x - 1)* blockDim%x + threadIdx%x
IF (iElem > nElems) RETURN
rhoLo = HUGE(1.0)
pLo = HUGE(1.0)
sLo = HUGE(1.0)
DO k = 0, Nloc
  DO j = 0, Nloc
    DO i = 0, Nloc
      rho = U(1,i,j,k,iElem)
      rhoLo = MIN(rhoLo, rho)
      IF (rho <= 0.0) THEN
        pLo = MIN(pLo, -1.0)
        sLo = MIN(sLo, -1.0E30)
        CYCLE
      ENDIF

      mom2 = U(2,i,j,k,iElem)**2 + U(3,i,j,k,iElem)**2 + U(4,i,j,k,iElem)**2
      pressure = (gammaGas-1.0)*(U(5,i,j,k,iElem) - 0.5*mom2/rho)
      pLo = MIN(pLo, pressure)
      if (pressure <= 0.0) THEN
        sLo = MIN(sLo, -1.0E30)
        CYCLE
      ENDIF
      entropy = LOG(pressure) - gammaGas*LOG(rho)
      sLo = MIN(sLo, entropy)
      END DO
  END DO
END DO
rhoMin(iElem) = rhoLo
pMin(iElem) = pLo
sMIN(iElem) = sLo
END SUBROUTINE EntropyMin_Kernel


  SUBROUTINE FinalizeEntropyFilter()
    IMPLICIT NONE
END SUBROUTINE FinalizeEntropyFilter

  END MODULE MOD_EntropyFilter
