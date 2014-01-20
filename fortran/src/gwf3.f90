module gwfmodule
  use ModelModule
  use PackageModule
  use global,only:gwfglotype,niunit
  use PARAMMODULE,only:gwfparamtype
  use GWFBASMODULE,only:gwfbastype
  use GWFBCFMODULE,only:gwfbcftype
  use GWFGHBMODULE,only:gwfghbtype
  use GWFWELMODULE,only:gwfweltype
  private
  public :: gwf3ar
  type, extends(modeltype) :: gwfmodeltype
    type(gwfglotype) :: gwfglodat
    type(gwfbastype) :: gwfbasdat
    type(gwfparamtype) :: gwfpardat
    type(gwfbcftype) :: gwfbcfdat
    type(gwfghbtype) :: gwfghbdat
    type(gwfweltype) :: gwfweldat
    contains
    procedure :: modelst=>gwf3st
    procedure :: modelrp=>gwf3rp
    procedure :: modelfmcalc=>gwf3fmcalc
    procedure :: modelfmfill=>gwf3fmfill
    procedure :: modelbd=>gwf3bd
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
  print *,'Creating gwfmodel: ', id
!
! -- Open the gwf name file
  call freeunitnumber(inunit)
  print *, 'opening model namefile on unit: ', inunit
  open(unit=inunit,file=filename,status='old')
!
! -- Allocate and read bas and dis and then save to pointers
  call gwfmodel%gwfglodat%pntset
  call gwfmodel%gwfbasdat%pntset
  call gwfmodel%gwfpardat%pntset
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
    call gwfmodel%gwfbcfdat%pntset
    CALL GWF2BCFU1AR(IUNIT(1),IUNIT(22),IUNIT(23))
    call gwfmodel%gwfbcfdat%pntsav
  ENDIF
!
! -- Allocate and read wel and then save to pointers
  IF(IUNIT(2).GT.0) THEN
    call gwfmodel%gwfweldat%pntset
    CALL GWF2WEL7U1AR(IUNIT(2))
    call gwfmodel%gwfweldat%pntsav
  ENDIF
!
! -- Allocate and read ghb and then save to pointers
  IF(IUNIT(7).GT.0) THEN
    call gwfmodel%gwfghbdat%pntset
    CALL GWF2GHB7U1AR(IUNIT(7))
    call gwfmodel%gwfghbdat%pntsav
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
! -- Print routine name and set model object pointers
  print *,'gwf3st'
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
!
! -- Print routine name and set model object pointers
    print *,'gwf3rp'
    call this%pntset
!
! -- Call package read and prepare subroutines
    IF(IUNIT(2).GT.0) CALL GWF2WEL7U1RP(IUNIT(2))
    IF(IUNIT(7).GT.0) CALL GWF2GHB7U1RP(IUNIT(7))
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
! -- Print routine name and set model object pointers
    print *,'gwf3ad'
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
! -- Print routine name and set model object pointers
    print *,'gwf3fmcalc'
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
    IF(IUNIT(2).GT.0) CALL GWF2WEL7U1FM
    IF(IUNIT(7).GT.0) CALL GWF2GHB7U1FM
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
! -- Print routine name and set model object pointers
    print *,'gwf3fmfill'
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
! Subroutine: (1) Fill the solution amat and rhs terms
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
! -- Print routine name and set model object pointers
    print *,'gwf3bd'
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
    IF(IUNIT(2).GT.0) CALL GWF2WEL7U1BD(kstp,kper)
    IF(IUNIT(7).GT.0) CALL GWF2GHB7U1BD(kstp,kper)
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagebd()
    enddo
!
! -- Output control output
    CALL GWF2BAS7OT(kstp,kper,ICNVG,1)
!
! -- return
    return
  end subroutine gwf3bd
  
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

  subroutine package_create(fname,filtyp,ipakid,packages)
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
    class(packagetype), pointer :: packobj
! ------------------------------------------------------------------------------
!
! -- Create and add packages to this models list of packages
    if(filtyp=='WEL') then
      call wel_create(packobj,fname,ipakid)
    elseif(filtyp=='GHB') then
      call ghb_create(packobj,fname,ipakid)
    endif
    call packages%setpackage(packobj,ipakid)
!
! -- return
    return
  end subroutine package_create  
  
end module gwfmodule