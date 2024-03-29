MODULE GWFRCHMODULE
  IMPLICIT NONE
  INTEGER,SAVE,POINTER :: NRCHOP
  INTEGER,SAVE,POINTER :: IRCHCB
  INTEGER,SAVE,POINTER :: MXNDRCH
  INTEGER,SAVE,POINTER :: NPRCH
  INTEGER,SAVE,POINTER :: IRCHPF
  INTEGER,SAVE,POINTER :: INIRCH
  INTEGER,SAVE,POINTER :: NIRCH
  REAL,SAVE,DIMENSION(:),POINTER :: RECH
  INTEGER,SAVE,DIMENSION(:),POINTER :: IRCH
  TYPE :: GWFRCHTYPE
    INTEGER,POINTER :: NRCHOP
    INTEGER,POINTER :: IRCHCB
    INTEGER,POINTER :: MXNDRCH
    INTEGER,POINTER :: NPRCH
    INTEGER,POINTER :: IRCHPF
    INTEGER,POINTER :: INIRCH
    INTEGER,POINTER :: NIRCH
    REAL,DIMENSION(:),POINTER :: RECH
    INTEGER,DIMENSION(:),POINTER :: IRCH
  CONTAINS
    PROCEDURE :: DESTROY
    PROCEDURE :: PNTSAV
    PROCEDURE :: PNTSET
  END TYPE GWFRCHTYPE
  CONTAINS
  SUBROUTINE DESTROY(THIS)
    IMPLICIT NONE
    CLASS(GWFRCHTYPE),INTENT(INOUT) :: THIS
    DEALLOCATE(THIS%NRCHOP)
    DEALLOCATE(THIS%IRCHCB)
    DEALLOCATE(THIS%MXNDRCH)
    DEALLOCATE(THIS%NPRCH)
    DEALLOCATE(THIS%IRCHPF)
    DEALLOCATE(THIS%INIRCH)
    DEALLOCATE(THIS%NIRCH)
    DEALLOCATE(THIS%RECH)
    DEALLOCATE(THIS%IRCH)
  END SUBROUTINE DESTROY
  SUBROUTINE PNTSAV(THIS)
    IMPLICIT NONE
    CLASS(GWFRCHTYPE),INTENT(INOUT) :: THIS
    THIS%NRCHOP=>NRCHOP
    THIS%IRCHCB=>IRCHCB
    THIS%MXNDRCH=>MXNDRCH
    THIS%NPRCH=>NPRCH
    THIS%IRCHPF=>IRCHPF
    THIS%INIRCH=>INIRCH
    THIS%NIRCH=>NIRCH
    THIS%RECH=>RECH
    THIS%IRCH=>IRCH
  END SUBROUTINE PNTSAV
  SUBROUTINE PNTSET(THIS)
    IMPLICIT NONE
    CLASS(GWFRCHTYPE),INTENT(INOUT) :: THIS
    NRCHOP=>THIS%NRCHOP
    IRCHCB=>THIS%IRCHCB
    MXNDRCH=>THIS%MXNDRCH
    NPRCH=>THIS%NPRCH
    IRCHPF=>THIS%IRCHPF
    INIRCH=>THIS%INIRCH
    NIRCH=>THIS%NIRCH
    RECH=>THIS%RECH
    IRCH=>THIS%IRCH
  END SUBROUTINE PNTSET
END MODULE GWFRCHMODULE
