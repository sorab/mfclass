module ModelModule
  use PackageModule
  implicit none
  private
  public :: model_list_init
  public :: modeltype
  public :: modellisttype
  public :: getmodel
  public :: model_create
  public :: modellist
!
  type :: modelcontainer
    class(modeltype), pointer :: obj
  end type modelcontainer
!
  type modellisttype
    integer :: nmodels = 0
    type(modelcontainer), allocatable, dimension(:) :: models
    contains
    procedure :: setmodel
    procedure :: getmodel
  end type modellisttype
  type(modellisttype) :: modellist
!
  type modeltype
      character(len=20) :: name
      integer :: id
      integer :: neq
      integer :: nja
      integer :: offset
      integer :: solutionid=0
      integer, allocatable, dimension(:) :: ia
      integer, allocatable, dimension(:) :: ja
      double precision, pointer, dimension(:) :: x
      double precision, pointer, dimension(:) :: rhs
      double precision, allocatable, dimension(:) :: cond
      integer, allocatable, dimension(:) :: idxglo
      type(packagelist) :: packages
      integer, dimension(10) :: models
      integer, dimension(10) :: crosses
      contains
      procedure :: disread
      procedure :: modelst
      procedure :: modelrp
      procedure :: modelad
      procedure :: modelfmcalc
      procedure :: modelfmfill
      procedure :: modelbd
      procedure :: modelot
      procedure :: modelbdentry
      procedure :: pntset=>modelpntset
      procedure :: modelsubinit
  end type modeltype

  contains

  subroutine model_list_init(nmodels)
! ******************************************************************************
! model_list_init -- Set size of modellist
! modellist is a list (actually an allocatable array) that contains all of
! the models that are part of the simulation.
! Subroutine: (1) allocate size of modellist%models
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    integer,intent(in) :: nmodels
! ------------------------------------------------------------------------------
!
    allocate(modellist%models(nmodels))
!
! -- return
    return
  end subroutine model_list_init

  subroutine setmodel(this, newmodel, ipos)
! ******************************************************************************
! setmodel -- Set Model
! For this modellist, set the pointer in position ipos to newsmodel
! Subroutine: (1) point the solution to newmodel
!             (2) increase this%nmodels to reflect the number of stored 
!                 models
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modellisttype) :: this
    class(modeltype), target, intent(in) :: newmodel
    integer, intent(in) :: ipos
! ------------------------------------------------------------------------------
!
    this%models(ipos)%obj => newmodel
    this%nmodels=max(ipos,this%nmodels)
!
! -- return
    return
  end subroutine setmodel

  subroutine getmodel(this, themodel, ipos)
! ******************************************************************************
! getmodel -- Get Model
! For this modellist, point the model to the one in position ipos
! Subroutine: (1) point the model to themodel 
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modellisttype) :: this
    class(modeltype), pointer :: themodel
! ------------------------------------------------------------------------------
!
    integer, intent(in) :: ipos
    themodel => this%models(ipos)%obj
!
! -- return
    return
  end subroutine getmodel

  subroutine model_create(filename,id)
! ******************************************************************************
! model_create -- Create a New model 
! Using the data in filename,  assign this new model an id number and store 
! the model in the modellist.
! Subroutine: (1) allocate model and assign id and name
!             (2) open the filename and read in the data 
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use SimModule,only:iout
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id
    character(len=300) :: line,fname
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,inunit,npackages,ipak
    real :: r
    type(modeltype), pointer :: model
! ------------------------------------------------------------------------------
!
! -- Create a new model and add it to the modellist container
    allocate(model)
    call modellist%setmodel(model,id)
    write(model%name,'(a4,i1)') 'MDL_',id
!
! -- Open the model name file
    call freeunitnumber(inunit)
    write(iout,'(/a,a)') ' Creating model: ', model%name
    WRITE(iout,'(a,a)') ' Using name file: ',trim(filename)
    WRITE(iout,'(a,i6)') ' On unit number: ',inunit
    open(unit=inunit,file=filename)
!
! -- Read the name file and create the packages
    ipak=1
    do
!
! -- Read line and skip if necessary
      read(inunit,'(a)',end=100) line
      if(line.eq.' ') cycle
      if(line(1:1).eq.'#') cycle
!
! -- Decode
      lloc=1
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
      filtyp=line(ityp1:ityp2)
!
! -- Simple DIS package created for this generic model
      if(filtyp.eq.'DIS') then
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
        fname=line(ityp1:ityp2)
        call model%disread(fname,id)
!
! -- Create packages
      elseif(filtyp.eq.'NPACKAGES') then
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,npackages,R,IOUT,INUNIT)
!
! -- Note that the size of model%packages is currently hardwired to 20
! -- todo: allocate(model%packages(npackages))
        ipak=1
!
! -- Create all other new-style packages (WEL, GHB, etc.)
      else
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
        fname=line(ityp1:ityp2)
        call package_create(fname,filtyp,ipak,model%packages)
        ipak=ipak+1
      endif
!
! -- Cycle loop to avoid exit
      cycle
100   exit        
    enddo
!
! -- close this file
    close(inunit)
!
! -- return
    return
  end subroutine model_create

  subroutine package_create(fname,filtyp,ipakid,packages)
! ******************************************************************************
! package_create -- Model Create Packages
! Subroutine: (1) add new-style packages to this models package list
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
! -- Now supporting new-style WEL and GHB packages
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


  subroutine disread(this,filename,id)
! ******************************************************************************
! disread -- Read Simplified DIS File
! This is not a GWF DIS input file, this is a file that contains the following:
!   neqs nodes
!   ia
!   ja
!   amat
!   rhs
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    character(len=*), intent(in) :: filename
    integer,intent(in) :: id
    character(len=300) :: line
    character(len=20) :: name
    integer :: lloc,ityp1,ityp2,n,iout,inunit
    real :: r
! ------------------------------------------------------------------------------
!
    this%id=id
    call freeunitnumber(inunit)
    print *, 'opening model dis file on unit: ', inunit
    open(unit=inunit,file=filename)
    !
    !problem size: neq, nja
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo    
    lloc=1
    CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,this%neq,R,IOUT,INUNIT)
    CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,this%nja,R,IOUT,INUNIT)
    !
    !allocate arrays
    allocate(this%ia(this%neq+1))
    allocate(this%ja(this%nja))
    allocate(this%cond(this%nja))
    allocate(this%rhs(this%neq))
    allocate(this%idxglo(this%nja))
    !
    !ia
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%ia
    !
    !ja
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%ja
    !
    !fully constructed amat matrix
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%cond
    !
    !rhs
    do
        read(inunit,'(a)') line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        exit
    enddo
    backspace(inunit)
    read(inunit,*) this%rhs
    close(inunit)
!
! -- return
    return
  end subroutine disread

  subroutine modelst(this)
! ******************************************************************************
! modelst -- Model Stress Timing
! Subroutine: (1) calls package st routines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('modelst')
    call this%pntset
!
! -- stress timing
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagest()
    enddo
!
! -- return
    return
  end subroutine modelst

  subroutine modelrp(this)
! ******************************************************************************
! modelrp -- Model Read and Prepare
! Subroutine: (1) calls package read and prepare routines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine
    call this%modelsubinit('modelrp')
    call this%pntset
!
! -- read and prepare
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagerp()
    enddo
!
! -- return
    return
  end subroutine modelrp

  subroutine modelad(this)
! ******************************************************************************
! modelad -- Model Time Step Advance
! Subroutine: (1) calls package advance subroutines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('modelad')
    call this%pntset
!
! -- advance
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagead()
    enddo
!
! -- return
    return
  end subroutine modelad

  subroutine modelfmcalc(this)
! ******************************************************************************
! modelfmcalc -- Model Calculate Terms for Formulation
! Subroutine: (1) calls package fm routines
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('modelfmcalc')
    call this%pntset
!
! -- This is where conductance should be calculated if it is dependent on x
!
! -- Calculate the rhs and hcof terms for each package
    do ip=1,this%packages%npackages
        call this%packages%getpackage(p,ip)
        call p%fmcalc()
    enddo
!
! -- return
    return
  end subroutine modelfmcalc

  subroutine modelfmfill(this,amatsln,njasln)
! ******************************************************************************
! modelfmfill -- Model Formulation Fill
! Subroutine: (1) Fill the solution amat and rhs terms
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    double precision,dimension(njasln),intent(inout) :: amatsln
    integer,intent(in) :: njasln
    class(packagetype), pointer :: p
    integer :: ip,n,i,ipos
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('modelfmfill')
    call this%pntset
!
! -- Copy the model conductance into the solution amat
    do ipos=1,this%nja
      amatsln(this%idxglo(ipos))=this%cond(ipos)
    enddo
    !
    !copy the package rhs and hcof into the rhs and amat
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
  end subroutine modelfmfill

  subroutine modelbd(this)
! ******************************************************************************
! modelbd -- Model Budget
! Subroutine: (1) Tabulate model budget
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    class(packagetype), pointer :: p
    integer :: ip
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('modelbd')
    call this%pntset
!
! -- budget
    do ip=1,this%packages%npackages
      call this%packages%getpackage(p,ip)
      call p%packagebd(this%x)
    enddo
!
! -- return
    return
  end subroutine modelbd

  subroutine modelot(this)
! ******************************************************************************
! modelot -- Model Output
! Subroutine: (1) Output budget items
! This subroutine should be overridden by a child model if necessary.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
! ------------------------------------------------------------------------------
!
! -- Call the simulation subroutine initialize routine and set pointers
    call this%modelsubinit('modelot')
    call this%pntset
!
! -- return
    return
  end subroutine modelot

  subroutine modelbdentry(this,text,rin,rout)
! ******************************************************************************
! modelbdentry -- Model Budget Entry
! Todo: need to write generic model budget capabilities
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
    character(len=*),intent(in) :: text
    real,intent(in) :: rin
    real,intent(in) :: rout
! ------------------------------------------------------------------------------
!
! -- return
    return
  end subroutine modelbdentry

  subroutine modelpntset(this)
! ******************************************************************************
! pntset -- Model Pointer Set
! This subroutine should be overridden by a child model if necessary.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(modeltype) :: this
! ------------------------------------------------------------------------------
!
! -- return
    return
  end subroutine modelpntset

  subroutine modelsubinit(this,subname)
! ******************************************************************************
! modelsubinit -- Model Subroutine Initialization
! This subroutine allows call stacks and additional information to be written
! depending on the verbosity level in simmodule.  It is called at the beginning
! of every model subroutine so it may serve other purposes as needed.
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use SimModule,only:iout,iverbose
    implicit none
    class(modeltype) :: this
    character(len=*) :: subname
    character(len=200) :: message
! ------------------------------------------------------------------------------
!
! -- construct message and send it to sim_message
    write(message,'(a,a,a,a)') '------Calling subroutine ',subname,            &
      ' on model ', this%name
    call sim_message(1,message)
!
! -- return
    return
  end subroutine modelsubinit


end module ModelModule
