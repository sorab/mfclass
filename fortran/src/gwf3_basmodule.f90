MODULE GWFBASMODULE
  IMPLICIT NONE
  INTEGER,SAVE,POINTER :: MSUM
  INTEGER,SAVE,POINTER :: IHEDFM
  INTEGER,SAVE,POINTER :: IHEDUN
  INTEGER,SAVE,POINTER :: IDDNFM
  INTEGER,SAVE,POINTER :: IDDNUN
  INTEGER,SAVE,POINTER :: IBOUUN
  INTEGER,SAVE,POINTER :: ISPCFM
  INTEGER,SAVE,POINTER :: ISPCUN
  INTEGER,SAVE,POINTER :: LBHDSV
  INTEGER,SAVE,POINTER :: LBDDSV
  INTEGER,SAVE,POINTER :: LBBOSV
  INTEGER,SAVE,POINTER :: IBUDFL
  INTEGER,SAVE,POINTER :: ICBCFL
  INTEGER,SAVE,POINTER :: IHDDFL
  INTEGER,SAVE,POINTER :: ISPCFL
  INTEGER,SAVE,POINTER :: IAUXSV
  INTEGER,SAVE,POINTER :: IBDOPT
  INTEGER,SAVE,POINTER :: IPRTIM
  INTEGER,SAVE,POINTER :: IPEROC
  INTEGER,SAVE,POINTER :: ITSOC
  INTEGER,SAVE,POINTER :: ICHFLG
  INTEGER,SAVE,POINTER :: IFRCNVG
  DOUBLE PRECISION,SAVE,POINTER :: DELT
  DOUBLE PRECISION,SAVE,POINTER :: PERTIM
  DOUBLE PRECISION,SAVE,POINTER :: TOTIM
  REAL,SAVE,POINTER :: HNOFLO
  CHARACTER(LEN=20),SAVE,POINTER :: CHEDFM
  CHARACTER(LEN=20),SAVE,POINTER :: CDDNFM
  CHARACTER(LEN=20),SAVE,POINTER :: CBOUFM
  CHARACTER(LEN=20),SAVE,POINTER :: CSPCFM
  INTEGER,SAVE,DIMENSION(:,:),POINTER :: IOFLG
  DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER :: VBVL
  CHARACTER(LEN=16),SAVE,DIMENSION(:),POINTER :: VBNM
  TYPE :: GWFBASTYPE
    INTEGER,POINTER :: MSUM
    INTEGER,POINTER :: IHEDFM
    INTEGER,POINTER :: IHEDUN
    INTEGER,POINTER :: IDDNFM
    INTEGER,POINTER :: IDDNUN
    INTEGER,POINTER :: IBOUUN
    INTEGER,POINTER :: ISPCFM
    INTEGER,POINTER :: ISPCUN
    INTEGER,POINTER :: LBHDSV
    INTEGER,POINTER :: LBDDSV
    INTEGER,POINTER :: LBBOSV
    INTEGER,POINTER :: IBUDFL
    INTEGER,POINTER :: ICBCFL
    INTEGER,POINTER :: IHDDFL
    INTEGER,POINTER :: ISPCFL
    INTEGER,POINTER :: IAUXSV
    INTEGER,POINTER :: IBDOPT
    INTEGER,POINTER :: IPRTIM
    INTEGER,POINTER :: IPEROC
    INTEGER,POINTER :: ITSOC
    INTEGER,POINTER :: ICHFLG
    INTEGER,POINTER :: IFRCNVG
    DOUBLE PRECISION,POINTER :: DELT
    DOUBLE PRECISION,POINTER :: PERTIM
    DOUBLE PRECISION,POINTER :: TOTIM
    REAL,POINTER :: HNOFLO
    CHARACTER(LEN=20),POINTER :: CHEDFM
    CHARACTER(LEN=20),POINTER :: CDDNFM
    CHARACTER(LEN=20),POINTER :: CBOUFM
    CHARACTER(LEN=20),POINTER :: CSPCFM
    INTEGER,DIMENSION(:,:),POINTER :: IOFLG
    DOUBLE PRECISION,DIMENSION(:,:),POINTER :: VBVL
    CHARACTER(LEN=16),DIMENSION(:),POINTER :: VBNM
  CONTAINS
    PROCEDURE :: DESTROY
    PROCEDURE :: PNTSAV
    PROCEDURE :: PNTSET
  END TYPE GWFBASTYPE
  CONTAINS
  SUBROUTINE DESTROY(THIS)
    IMPLICIT NONE
    CLASS(GWFBASTYPE),INTENT(INOUT) :: THIS
    DEALLOCATE(THIS%MSUM)
    DEALLOCATE(THIS%IHEDFM)
    DEALLOCATE(THIS%IHEDUN)
    DEALLOCATE(THIS%IDDNFM)
    DEALLOCATE(THIS%IDDNUN)
    DEALLOCATE(THIS%IBOUUN)
    DEALLOCATE(THIS%ISPCFM)
    DEALLOCATE(THIS%ISPCUN)
    DEALLOCATE(THIS%LBHDSV)
    DEALLOCATE(THIS%LBDDSV)
    DEALLOCATE(THIS%LBBOSV)
    DEALLOCATE(THIS%IBUDFL)
    DEALLOCATE(THIS%ICBCFL)
    DEALLOCATE(THIS%IHDDFL)
    DEALLOCATE(THIS%ISPCFL)
    DEALLOCATE(THIS%IAUXSV)
    DEALLOCATE(THIS%IBDOPT)
    DEALLOCATE(THIS%IPRTIM)
    DEALLOCATE(THIS%IPEROC)
    DEALLOCATE(THIS%ITSOC)
    DEALLOCATE(THIS%ICHFLG)
    DEALLOCATE(THIS%IFRCNVG)
    DEALLOCATE(THIS%DELT)
    DEALLOCATE(THIS%PERTIM)
    DEALLOCATE(THIS%TOTIM)
    DEALLOCATE(THIS%HNOFLO)
    DEALLOCATE(THIS%CHEDFM)
    DEALLOCATE(THIS%CDDNFM)
    DEALLOCATE(THIS%CBOUFM)
    DEALLOCATE(THIS%CSPCFM)
    DEALLOCATE(THIS%IOFLG)
    DEALLOCATE(THIS%VBVL)
    DEALLOCATE(THIS%VBNM)
  END SUBROUTINE DESTROY
  SUBROUTINE PNTSAV(THIS)
    IMPLICIT NONE
    CLASS(GWFBASTYPE),INTENT(INOUT) :: THIS
    THIS%MSUM=>MSUM
    THIS%IHEDFM=>IHEDFM
    THIS%IHEDUN=>IHEDUN
    THIS%IDDNFM=>IDDNFM
    THIS%IDDNUN=>IDDNUN
    THIS%IBOUUN=>IBOUUN
    THIS%ISPCFM=>ISPCFM
    THIS%ISPCUN=>ISPCUN
    THIS%LBHDSV=>LBHDSV
    THIS%LBDDSV=>LBDDSV
    THIS%LBBOSV=>LBBOSV
    THIS%IBUDFL=>IBUDFL
    THIS%ICBCFL=>ICBCFL
    THIS%IHDDFL=>IHDDFL
    THIS%ISPCFL=>ISPCFL
    THIS%IAUXSV=>IAUXSV
    THIS%IBDOPT=>IBDOPT
    THIS%IPRTIM=>IPRTIM
    THIS%IPEROC=>IPEROC
    THIS%ITSOC=>ITSOC
    THIS%ICHFLG=>ICHFLG
    THIS%IFRCNVG=>IFRCNVG
    THIS%DELT=>DELT
    THIS%PERTIM=>PERTIM
    THIS%TOTIM=>TOTIM
    THIS%HNOFLO=>HNOFLO
    THIS%CHEDFM=>CHEDFM
    THIS%CDDNFM=>CDDNFM
    THIS%CBOUFM=>CBOUFM
    THIS%CSPCFM=>CSPCFM
    THIS%IOFLG=>IOFLG
    THIS%VBVL=>VBVL
    THIS%VBNM=>VBNM
  END SUBROUTINE PNTSAV
  SUBROUTINE PNTSET(THIS)
    IMPLICIT NONE
    CLASS(GWFBASTYPE),INTENT(INOUT) :: THIS
    MSUM=>THIS%MSUM
    IHEDFM=>THIS%IHEDFM
    IHEDUN=>THIS%IHEDUN
    IDDNFM=>THIS%IDDNFM
    IDDNUN=>THIS%IDDNUN
    IBOUUN=>THIS%IBOUUN
    ISPCFM=>THIS%ISPCFM
    ISPCUN=>THIS%ISPCUN
    LBHDSV=>THIS%LBHDSV
    LBDDSV=>THIS%LBDDSV
    LBBOSV=>THIS%LBBOSV
    IBUDFL=>THIS%IBUDFL
    ICBCFL=>THIS%ICBCFL
    IHDDFL=>THIS%IHDDFL
    ISPCFL=>THIS%ISPCFL
    IAUXSV=>THIS%IAUXSV
    IBDOPT=>THIS%IBDOPT
    IPRTIM=>THIS%IPRTIM
    IPEROC=>THIS%IPEROC
    ITSOC=>THIS%ITSOC
    ICHFLG=>THIS%ICHFLG
    IFRCNVG=>THIS%IFRCNVG
    DELT=>THIS%DELT
    PERTIM=>THIS%PERTIM
    TOTIM=>THIS%TOTIM
    HNOFLO=>THIS%HNOFLO
    CHEDFM=>THIS%CHEDFM
    CDDNFM=>THIS%CDDNFM
    CBOUFM=>THIS%CBOUFM
    CSPCFM=>THIS%CSPCFM
    IOFLG=>THIS%IOFLG
    VBVL=>THIS%VBVL
    VBNM=>THIS%VBNM
  END SUBROUTINE PNTSET
END MODULE GWFBASMODULE