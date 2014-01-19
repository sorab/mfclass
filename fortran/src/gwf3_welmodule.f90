MODULE GWFWELMODULE
  IMPLICIT NONE
  INTEGER,SAVE,POINTER :: NWELLS
  INTEGER,SAVE,POINTER :: MXWELL
  INTEGER,SAVE,POINTER :: NWELVL
  INTEGER,SAVE,POINTER :: IWELCB
  INTEGER,SAVE,POINTER :: IPRWEL
  INTEGER,SAVE,POINTER :: IAFR
  INTEGER,SAVE,POINTER :: NPWEL
  INTEGER,SAVE,POINTER :: IWELPB
  INTEGER,SAVE,POINTER :: NNPWEL
  INTEGER,SAVE,POINTER :: IWELQV
  INTEGER,SAVE,POINTER :: NNPWCLN
  CHARACTER(LEN=16),SAVE,DIMENSION(:),POINTER :: WELAUX
  REAL,SAVE,DIMENSION(:,:),POINTER :: WELL
  TYPE :: GWFWELTYPE
    INTEGER,POINTER :: NWELLS
    INTEGER,POINTER :: MXWELL
    INTEGER,POINTER :: NWELVL
    INTEGER,POINTER :: IWELCB
    INTEGER,POINTER :: IPRWEL
    INTEGER,POINTER :: IAFR
    INTEGER,POINTER :: NPWEL
    INTEGER,POINTER :: IWELPB
    INTEGER,POINTER :: NNPWEL
    INTEGER,POINTER :: IWELQV
    INTEGER,POINTER :: NNPWCLN
    CHARACTER(LEN=16),DIMENSION(:),POINTER :: WELAUX
    REAL,DIMENSION(:,:),POINTER :: WELL
  CONTAINS
    PROCEDURE :: DESTROY
    PROCEDURE :: PNTSAV
    PROCEDURE :: PNTSET
  END TYPE GWFWELTYPE
  CONTAINS
  SUBROUTINE DESTROY(THIS)
    IMPLICIT NONE
    CLASS(GWFWELTYPE),INTENT(INOUT) :: THIS
    DEALLOCATE(THIS%NWELLS)
    DEALLOCATE(THIS%MXWELL)
    DEALLOCATE(THIS%NWELVL)
    DEALLOCATE(THIS%IWELCB)
    DEALLOCATE(THIS%IPRWEL)
    DEALLOCATE(THIS%IAFR)
    DEALLOCATE(THIS%NPWEL)
    DEALLOCATE(THIS%IWELPB)
    DEALLOCATE(THIS%NNPWEL)
    DEALLOCATE(THIS%IWELQV)
    DEALLOCATE(THIS%NNPWCLN)
    DEALLOCATE(THIS%WELAUX)
    DEALLOCATE(THIS%WELL)
  END SUBROUTINE DESTROY
  SUBROUTINE PNTSAV(THIS)
    IMPLICIT NONE
    CLASS(GWFWELTYPE),INTENT(INOUT) :: THIS
    THIS%NWELLS=>NWELLS
    THIS%MXWELL=>MXWELL
    THIS%NWELVL=>NWELVL
    THIS%IWELCB=>IWELCB
    THIS%IPRWEL=>IPRWEL
    THIS%IAFR=>IAFR
    THIS%NPWEL=>NPWEL
    THIS%IWELPB=>IWELPB
    THIS%NNPWEL=>NNPWEL
    THIS%IWELQV=>IWELQV
    THIS%NNPWCLN=>NNPWCLN
    THIS%WELAUX=>WELAUX
    THIS%WELL=>WELL
  END SUBROUTINE PNTSAV
  SUBROUTINE PNTSET(THIS)
    IMPLICIT NONE
    CLASS(GWFWELTYPE),INTENT(INOUT) :: THIS
    NWELLS=>THIS%NWELLS
    MXWELL=>THIS%MXWELL
    NWELVL=>THIS%NWELVL
    IWELCB=>THIS%IWELCB
    IPRWEL=>THIS%IPRWEL
    IAFR=>THIS%IAFR
    NPWEL=>THIS%NPWEL
    IWELPB=>THIS%IWELPB
    NNPWEL=>THIS%NNPWEL
    IWELQV=>THIS%IWELQV
    NNPWCLN=>THIS%NNPWCLN
    WELAUX=>THIS%WELAUX
    WELL=>THIS%WELL
  END SUBROUTINE PNTSET
END MODULE GWFWELMODULE
