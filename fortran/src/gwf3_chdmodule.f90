MODULE GWFCHDMODULE
  IMPLICIT NONE
  INTEGER,SAVE,POINTER :: NCHDS
  INTEGER,SAVE,POINTER :: MXCHD
  INTEGER,SAVE,POINTER :: NCHDVL
  INTEGER,SAVE,POINTER :: IPRCHD
  INTEGER,SAVE,POINTER :: NPCHD
  INTEGER,SAVE,POINTER :: ICHDPB
  INTEGER,SAVE,POINTER :: NNPCHD
  CHARACTER(LEN=16),SAVE,DIMENSION(:),POINTER :: CHDAUX
  REAL,SAVE,DIMENSION(:,:),POINTER :: CHDS
  TYPE :: GWFCHDTYPE
    INTEGER,POINTER :: NCHDS
    INTEGER,POINTER :: MXCHD
    INTEGER,POINTER :: NCHDVL
    INTEGER,POINTER :: IPRCHD
    INTEGER,POINTER :: NPCHD
    INTEGER,POINTER :: ICHDPB
    INTEGER,POINTER :: NNPCHD
    CHARACTER(LEN=16),DIMENSION(:),POINTER :: CHDAUX
    REAL,DIMENSION(:,:),POINTER :: CHDS
  CONTAINS
    PROCEDURE :: DESTROY
    PROCEDURE :: PNTSAV
    PROCEDURE :: PNTSET
  END TYPE GWFCHDTYPE
  CONTAINS
  SUBROUTINE DESTROY(THIS)
    IMPLICIT NONE
    CLASS(GWFCHDTYPE),INTENT(INOUT) :: THIS
    DEALLOCATE(THIS%NCHDS)
    DEALLOCATE(THIS%MXCHD)
    DEALLOCATE(THIS%NCHDVL)
    DEALLOCATE(THIS%IPRCHD)
    DEALLOCATE(THIS%NPCHD)
    DEALLOCATE(THIS%ICHDPB)
    DEALLOCATE(THIS%NNPCHD)
    DEALLOCATE(THIS%CHDAUX)
    DEALLOCATE(THIS%CHDS)
  END SUBROUTINE DESTROY
  SUBROUTINE PNTSAV(THIS)
    IMPLICIT NONE
    CLASS(GWFCHDTYPE),INTENT(INOUT) :: THIS
    THIS%NCHDS=>NCHDS
    THIS%MXCHD=>MXCHD
    THIS%NCHDVL=>NCHDVL
    THIS%IPRCHD=>IPRCHD
    THIS%NPCHD=>NPCHD
    THIS%ICHDPB=>ICHDPB
    THIS%NNPCHD=>NNPCHD
    THIS%CHDAUX=>CHDAUX
    THIS%CHDS=>CHDS
  END SUBROUTINE PNTSAV
  SUBROUTINE PNTSET(THIS)
    IMPLICIT NONE
    CLASS(GWFCHDTYPE),INTENT(INOUT) :: THIS
    NCHDS=>THIS%NCHDS
    MXCHD=>THIS%MXCHD
    NCHDVL=>THIS%NCHDVL
    IPRCHD=>THIS%IPRCHD
    NPCHD=>THIS%NPCHD
    ICHDPB=>THIS%ICHDPB
    NNPCHD=>THIS%NNPCHD
    CHDAUX=>THIS%CHDAUX
    CHDS=>THIS%CHDS
  END SUBROUTINE PNTSET
END MODULE GWFCHDMODULE