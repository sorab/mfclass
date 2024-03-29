module gwfmodule
  use ModelModule
  use PackageModule
  use global,only:gwfglotype,niunit
  use PARAMMODULE,only:gwfparamtype
  use GWFBASMODULE,only:gwfbastype
  use GWFBCFMODULE,only:gwfbcftype
  use GWFGHBMODULE,only:gwfghbtype
  use GWFWELMODULE,only:gwfweltype
  use GWFCHDMODULE,only:gwfchdtype
  use GWFDRNMODULE,only:gwfdrntype
  use GWFEVTMODULE,only:gwfevttype
  use GWFFHBMODULE,only:gwffhbtype
  use GWFHFBMODULE,only:gwfhfbtype
  use GWFRCHMODULE,only:gwfrchtype
  use GWFRIVMODULE,only:gwfrivtype
  use GWFSTRMODULE,only:gwfstrtype
  private
  public :: gwf3ar
  type, extends(modeltype) :: gwfmodeltype
    type(gwfglotype) :: gwfglodat
    type(gwfbastype) :: gwfbasdat
    type(gwfparamtype) :: gwfpardat
    type(gwfbcftype) :: gwfbcfdat
    type(gwfghbtype) :: gwfghbdat
    type(gwfweltype) :: gwfweldat
    type(gwfchdtype) :: gwfchddat
    type(gwfdrntype) :: gwfdrndat
    type(gwfevttype) :: gwfevtdat
    type(gwffhbtype) :: gwffhbdat
    type(gwfhfbtype) :: gwfhfbdat
    type(gwfrchtype) :: gwfrchdat
    type(gwfrivtype) :: gwfrivdat
    type(gwfstrtype) :: gwfstrdat
    contains
    procedure :: modelst=>gwf3st
    procedure :: modelrp=>gwf3rp
    procedure :: modelad=>gwf3ad
    procedure :: modelfmcalc=>gwf3fmcalc
    procedure :: modelfmfill=>gwf3fmfill
    procedure :: modelbd=>gwf3bd
    procedure :: modelot=>gwf3ot
    procedure :: modelbdentry=>gwf3bdentry
    procedure :: pntset=>gwf3pntset
  end type gwfmodeltype
  CHARACTER*40 VERSION
  CHARACTER*10 MFVNAM
  integer :: maxunit
  CHARACTER*80 HEADNG(2)
  PARAMETER (VERSION='GWF 3 01/18/2014')
  PARAMETER (MFVNAM='-2015') 
  INTEGER,DIMENSION(:),POINTER :: IUNIT
  CHARACTER*4 CUNIT(NIUNIT)
  DATA CUNIT/   'BCF6', 'WEL ', 'DRN ', 'RIV ', 'EVT ', '    ', 'GHB ',  & !  7
                'RCH ', '    ', '    ', '    ', 'OC  ', 'SMS ', 'PCB ',  & ! 14
                'BCT ', 'FHB ', 'RES ', 'STR ', 'IBS ', 'CHD ', 'HFB6',  & ! 21
                'LAK ', 'LPF ', 'DIS ', 'DISU', 'PVAL', '    ', 'HOB ',  & ! 28
                'CLN ', '    ', 'ZONE', 'MULT', 'DROB', 'RVOB', 'GBOB',  & ! 35
                'GNC ', '    ', 'CHOB', 'ETS ', 'DRT ', '    ', 'GMG ',  & ! 42
                'hyd ', 'SFR ', '    ', 'GAGE', 'LVDA', '    ', 'lmt6',  & ! 49
                'MNW1', '    ', '    ', 'KDEP', 'SUB ', 'UZF ', 'gwm ',  & ! 56
                'SWT ', '    ', '    ', '    ', '    ', '    ', '    ',  & ! 63
                37*'    '/
  contains

  subroutine gwf3ar(filename,id,ioutsim)
! ******************************************************************************
! gwf3ar -- GroundWater Flow Model Allocate and Read
! Subroutine: (1) creates the gwfmodel object, 
!             (2) adds the model object to modellist,
!             (3) allocates and reads packages part of this model,
!             (4) allocates memory for arrays part of this model object
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  use SimModule,only:iout
  use global,only: ia,ja
  implicit none
  character(len=*),intent(in) :: filename
  integer,intent(in) :: id,ioutsim
  integer :: inunit
  type(gwfmodeltype), pointer :: gwfmodel
! ------------------------------------------------------------------------------
!
! -- Create a new model and add it to the modellist container
  allocate(gwfmodel)
  call modellist%setmodel(gwfmodel,id)
  write(gwfmodel%name,'(a4,i1)') 'GWF_',id
!
! -- Open the gwf name file
  call freeunitnumber(inunit)
  write(iout,'(/a,a)') ' Creating model: ', gwfmodel%name
  WRITE(iout,'(a,a)') ' Using name file: ',trim(filename)
  WRITE(iout,'(a,i6)') ' On unit number: ',inunit
  open(unit=inunit,file=filename,status='old')
!
! -- Allocate and read bas and dis and then save to pointers
  CALL GLO2BAS8AR(INUNIT,CUNIT,VERSION,24,31,32,MAXUNIT,12,                    &
                    HEADNG,26,MFVNAM,29,27,30,36)
  call gwfmodel%gwfglodat%pntsav
  call gwfmodel%gwfbasdat%pntsav
  call gwfmodel%gwfpardat%pntsav
!
! -- Make a pointer to iunit for easier syntax within these routines
  iunit=>gwfmodel%gwfglodat%iunit
!
! -- Allocate and read bcf/lpf and then save to pointers
  IF(IUNIT(1).GT.0.OR.IUNIT(23).GT.0)THEN                                     
    CALL GWF2BCFU1AR(IUNIT(1),IUNIT(22),IUNIT(23))
    call gwfmodel%gwfbcfdat%pntsav
  ENDIF
!
! -- Allocate and read wel and then save to pointers
  IF(IUNIT(2).GT.0) THEN
    CALL GWF2WEL7U1AR(IUNIT(2))
    call gwfmodel%gwfweldat%pntsav
  ENDIF
!
! -- Allocate and read ghb and then save to pointers
  IF(IUNIT(7).GT.0) THEN
    CALL GWF2GHB7U1AR(IUNIT(7))
    call gwfmodel%gwfghbdat%pntsav
  ENDIF
!
! -- Allocate and read chd and then save to pointers
  IF(IUNIT(20).GT.0) THEN
    CALL GWF2CHD7U1AR(IUNIT(20))
    call gwfmodel%gwfchddat%pntsav
  ENDIF  
!
! -- Allocate and read drn and then save to pointers
  IF(IUNIT(3).GT.0) THEN
    CALL GWF2DRN7U1AR(IUNIT(3))
    call gwfmodel%gwfdrndat%pntsav
  ENDIF
!
! -- Allocate and read evt and then save to pointers
  IF(IUNIT(5).GT.0) THEN
    CALL GWF2EVT8U1AR(IUNIT(5),IUNIT(15))
    call gwfmodel%gwfevtdat%pntsav
  ENDIF
!
! -- Allocate and read fhb and then save to pointers
  IF(IUNIT(16).GT.0) THEN
    CALL GWF2FHB7U1AR(IUNIT(16))
    call gwfmodel%gwffhbdat%pntsav
  ENDIF
!
! -- Allocate and read hfb and then save to pointers
  IF(IUNIT(21).GT.0) THEN
    CALL GWF2HFB7U1AR(IUNIT(21))
    call gwfmodel%gwfhfbdat%pntsav
  ENDIF
!
! -- Allocate and read rch and then save to pointers
  IF(IUNIT(8).GT.0) THEN
    CALL GWF2RCH8U1AR(IUNIT(8))
    call gwfmodel%gwfrchdat%pntsav
  ENDIF
!
! -- Allocate and read riv and then save to pointers
  IF(IUNIT(4).GT.0) THEN
    CALL GWF2RIV7U1AR(IUNIT(4))
    call gwfmodel%gwfrivdat%pntsav
  ENDIF
!
! -- Allocate and read str and then save to pointers
  IF(IUNIT(18).GT.0) THEN
    CALL GWF2STR7U1AR(IUNIT(18))
    call gwfmodel%gwfstrdat%pntsav
  ENDIF

!
! -- Set and allocate gwfmodel variables and arrays
! -- langevin mfs015 todo: improve memory management. there is a copy
! -- of ia and ja in the model and in global
  gwfmodel%neq=gwfmodel%gwfglodat%neqs
  gwfmodel%nja=gwfmodel%gwfglodat%nja    
  allocate(gwfmodel%ia(gwfmodel%neq+1))
  allocate(gwfmodel%ja(gwfmodel%nja))
  allocate(gwfmodel%cond(gwfmodel%nja))
  allocate(gwfmodel%rhs(gwfmodel%neq))
  allocate(gwfmodel%idxglo(gwfmodel%nja))
  gwfmodel%ia(:) = ia
  gwfmodel%ja(:) = ja
!
! -- return
  return
  end subroutine gwf3ar

  subroutine gwf3st(this)
! ******************************************************************************
! gwf3st -- GroundWater Flow Model Stress Timing
! Subroutine: (1) copies head to model object x
!             (2) calls package st routines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use tdismodule,only:kper
    implicit none
    class(gwfmodeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip,i
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('gwf3st')
    call this%pntset
!
! -- Copy the starting heads into this%x
  if(kper==1) then
    do i=1,this%neq
      this%x(i)=this%gwfglodat%hnew(i)
    enddo
  endif
!
! -- Call package stress-timing subroutines
  CALL GWF2BAS8ST(kper)
  do ip=1,this%packages%npackages
    call this%packages%getpackage(p,ip)
    call p%packagest()
  enddo
!
! -- Return
    return
  end subroutine gwf3st

  subroutine gwf3rp(this)
! ******************************************************************************
! gwf3rp -- GroundWater Flow Model Read and Prepare
! Subroutine: (1) calls package read and prepare routines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(gwfmodeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine
    call this%modelsubinit('gwf3rp')
    call this%pntset
!
! -- Call package read and prepare subroutines
    IF(IUNIT(2).GT.0) CALL GWF2WEL7U1RP(IUNIT(2))
    IF(IUNIT(3).GT.0) CALL GWF2DRN7U1RP(IUNIT(3))
    IF(IUNIT(4).GT.0) CALL GWF2RIV7U1RP(IUNIT(4))
    IF(IUNIT(5).GT.0) CALL GWF2EVT8U1RP(IUNIT(5))
    IF(IUNIT(7).GT.0) CALL GWF2GHB7U1RP(IUNIT(7))
    IF(IUNIT(8).GT.0) CALL GWF2RCH8U1RP(IUNIT(8))
    IF(IUNIT(18).GT.0) CALL GWF2STR7U1RP(IUNIT(18))
    IF(IUNIT(20).GT.0) CALL GWF2CHD7U1RP(IUNIT(20))
!
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagerp()
    enddo
!
! -- Return
    return
  end subroutine gwf3rp

  subroutine gwf3ad(this)
! ******************************************************************************
! gwf3ad -- GroundWater Flow Model Time Step Advance
! Subroutine: (1) calls package advance subroutines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use tdismodule,only:kstp,kper
    implicit none
    class(gwfmodeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('gwf3ad')
    call this%pntset
!
! -- Call package read and prepare subroutines
    CALL GWF2BAS7AD(kper,kstp)
    IF(IUNIT(1).GT.0) CALL GWF2BCFU1AD(kper)
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagead()
    enddo
!
! -- return
    return
  end subroutine gwf3ad
  
  subroutine gwf3fmcalc(this)
! ******************************************************************************
! gwf3fmcalc -- GroundWater Flow Model Calculate Terms for Formulation
! Subroutine: (1) calls package fm routines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use tdismodule,only:kstp,kper
    implicit none
    class(gwfmodeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip,ipos,kkiter
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('gwf3fmcalc')
    call this%pntset
!
! -- Because global::hnew and this%x do not share memory, need to copy latest
! -- solution into hnew so that amat can be calculated using latest heads.
! -- langevin mf2015 todo: memory management of x and hnew
    this%gwfglodat%hnew(:)=this%x(:)
!
! -- Call package fm routines
    kkiter=1
    CALL GWF2BAS7U1FM
    IF(IUNIT(1).GT.0) CALL GWF2BCFU1FM(KKITER,kstp,kper)
    IF(IUNIT(21).GT.0) CALL GWF2HFB7U1FM
    IF(IUNIT(2).GT.0) CALL GWF2WEL7U1FM
    IF(IUNIT(3).GT.0) CALL GWF2DRN7U1FM
    IF(IUNIT(4).GT.0) CALL GWF2RIV7U1FM
    IF(IUNIT(5).GT.0) CALL GWF2EVT8U1FM
    IF(IUNIT(7).GT.0) CALL GWF2GHB7U1FM
    IF(IUNIT(8).GT.0) CALL GWF2RCH8U1FM(kper)
    IF(IUNIT(16).GT.0) CALL GWF2FHB7U1FM
    IF(IUNIT(18).GT.0) CALL GWF2STR7U1FM
!
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%fmcalc()
    enddo
!
! -- return
    return
  end subroutine gwf3fmcalc

  subroutine gwf3fmfill(this,amatsln,njasln)
! ******************************************************************************
! gwf3fmfill -- GroundWater Flow Model Formulation Fill
! Subroutine: (1) Fill the solution amat and rhs terms
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(gwfmodeltype) :: this
    double precision,dimension(njasln),intent(inout) :: amatsln
    integer,intent(in) :: njasln
    class(packagetype), pointer :: p
    integer :: ip,n,i,ipos
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('gwf3fmfill')
    call this%pntset
!
! -- Copy gwf amat into this model conductance
! -- Then fill this solution amat with this model conductance
! -- langevin mf2015 todo: better memory managment of amatsln,amat,cond
    do ipos=1,this%nja
        this%cond(ipos)=this%gwfglodat%amat(ipos)
        amatsln(this%idxglo(ipos))=this%cond(ipos)
    enddo
!
! -- Copy the rhs
    do ipos=1,this%neq
      this%rhs(ipos)=this%rhs(ipos)+this%gwfglodat%rhs(ipos)
    enddo
!
! -- Copy package rhs and hcof into solution amat and rhs
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        do i=1,p%nbound
            n=p%nodelist(i)
            this%rhs(n)=this%rhs(n)+p%rhs(i)
            ipos=this%ia(n)
            amatsln(this%idxglo(ipos))=amatsln(this%idxglo(ipos))+p%hcof(i)
        enddo
    enddo
!
! -- return
    return
  end subroutine gwf3fmfill

  subroutine gwf3bd(this)
! ******************************************************************************
! gwf3bd -- GroundWater Flow Model Budget
! Subroutine: (1) Tabulate Groundwater Flow Model Budget
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use tdismodule,only:kstp,kper
    use GLOBAL,only:flowja
    implicit none
    class(gwfmodeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip,icnvg
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('gwf3bd')
    call this%pntset
!      
! -- langevin mf2015 todo: get icnvg in here, also allocating flowja for budget
! -- calculations.
    icnvg=1 
    ALLOCATE(FLOWJA(this%gwfglodat%NJA))
!
! -- Push the solution results into hnew
! -- langevin mf2015 todo: memory management of x and hnew
    this%gwfglodat%hnew(:)=this%x(:)
!
! -- Output control
    CALL GWF2BAS7OC(kstp,kper,ICNVG,IUNIT(12))
    this%gwfbasdat%MSUM = 1
!
! -- Budget routines
    IF (IUNIT(1).GT.0) CALL GWF2BCFU1BDS(kstp,kper)
    IF (IUNIT(1).GT.0) THEN
      CALL GWF2BCFU1BDADJ(kstp,kper)
    ENDIF
    IF (IUNIT(1).GT.0) THEN
      CALL GWF2BCFU1BDCHWR(kstp,kper)  
      CALL GWF2BCFU1BDADJWR(kstp,kper)
    ENDIF
    DEALLOCATE(FLOWJA)
!
! -- boundary packages
    IF(IUNIT(2).GT.0) CALL GWF2WEL7U1BD(kstp,kper)
    IF(IUNIT(3).GT.0) CALL GWF2DRN7U1BD(KSTP,KPER)
    IF(IUNIT(4).GT.0) CALL GWF2RIV7U1BD(KSTP,KPER)
    IF(IUNIT(5).GT.0) CALL GWF2EVT8U1BD(KSTP,KPER,IUNIT(15))
    IF(IUNIT(7).GT.0) CALL GWF2GHB7U1BD(kstp,kper)
    IF(IUNIT(8).GT.0) CALL GWF2RCH8U1BD(KSTP,KPER)
    IF(IUNIT(16).GT.0) CALL GWF2FHB7U1BD(KSTP,KPER)
    IF(IUNIT(18).GT.0) CALL GWF2STR7U1BD(KSTP,KPER)    
!
! -- modflow2015 packages
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagebd(this%x)
      call this%modelbdentry(p%name,p%rin,p%rout)
    enddo
!
! -- return
    return
  end subroutine gwf3bd
  
  subroutine gwf3ot(this)
! ******************************************************************************
! gwf3ot -- GroundWater Flow Model Output
! Subroutine: (1) Output budget items
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use tdismodule,only:kstp,kper
    implicit none
    class(gwfmodeltype) :: this
    integer :: icnvg
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('gwf3ot')
    call this%pntset
!      
! -- langevin mf2015 todo: get icnvg in here, also allocating flowja for budget
! -- calculations.
    icnvg=1 
!
! -- Output control output
    CALL GWF2BAS7OT(kstp,kper,ICNVG,1)
!
! -- return
    return
  end subroutine gwf3ot
  
  subroutine gwf3pntset(this)
! ******************************************************************************
! gwf3pntset -- GroundWater Flow Model Pointer Set
! This subroutine should be called each time this model object is used.
! Subroutine: (1) set the pointers in global, gwfbasmodule, gwfbcfmodule, etc.
!                 using the type members of gwfglodat, gwfbasdat, etc.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(gwfmodeltype) :: this
! ------------------------------------------------------------------------------
!
    call this%gwfglodat%pntset
    iunit=>this%gwfglodat%iunit
    call this%gwfbasdat%pntset
    call this%gwfpardat%pntset
    if(iunit(1).gt.0) call this%gwfbcfdat%pntset
    if(iunit(2).gt.0) call this%gwfweldat%pntset
    if(iunit(7).gt.0) call this%gwfghbdat%pntset
!
! -- return
    return
  end subroutine gwf3pntset

  subroutine package_create(fname,filtyp,ipakid,packages,inunit,iout)
! ******************************************************************************
! NOT USED YET for GWF.  If we decide to use the package list approach, then
! we will need to get information back from basopen in order to create
! individual packages and add them to the package list.  We will also need to
! store the rhs and hcof package contributions.
! package_create -- GroundWater Flow Model Create Packages
! Subroutine: (1) add packages to this models package list
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use PackageModule
    use welmodule
    use ghbmodule
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: filtyp
    integer,intent(in) :: ipakid
    type(packagelist),intent(inout) :: packages
    integer,intent(in) :: inunit
    integer,intent(in) :: iout
    class(packagetype), pointer :: packobj
! ------------------------------------------------------------------------------
!
! -- Create and add packages to this models list of packages
    if(filtyp=='WEL') then
      call wel_create(packobj,fname,ipakid,inunit,iout)
    elseif(filtyp=='GHB') then
      call ghb_create(packobj,fname,ipakid,inunit,iout)
    endif
    call packages%setpackage(packobj,ipakid)
!
! -- return
    return
  end subroutine package_create  
 
  subroutine gwf3bdentry(this,text,rin,rout)
! ******************************************************************************
! gwf3bdentry -- GroundWater Flow Model Budget Entry
! This subroutine adds a budget entry to the flow budget.  It was added as
! a method for the gwf3 model object so that the cross object could add its
! contributions.
! Subroutine: (1) sets the pointer
!             (2) adds the entry to the vbvl array in gwfbasmodule
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    USE GWFBASMODULE,ONLY:MSUM,DELT,VBVL,VBNM
    implicit none
    class(gwfmodeltype) :: this
    character(len=*),intent(in) :: text
    real,intent(in) :: rin
    real,intent(in) :: rout
! ------------------------------------------------------------------------------
!
    call this%pntset
    VBVL(3,MSUM)=RIN
    VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
    VBVL(4,MSUM)=ROUT
    VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
    VBNM(MSUM)=TEXT
    MSUM=MSUM+1
!
! -- return
    return
  end subroutine gwf3bdentry
  
  end module gwfmodule
  