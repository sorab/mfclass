
module SOLUTIONMODULE
  use ModelModule
  use CrossModule
  use sparsemodule, only:sparsematrix
  use SMSMODULE
  !use MSMSCLASS
  use XMDMODULE
  use MXMDCLASS
  use PCGUMODULE
  use MPCGUCLASS
  implicit none
  private
  public :: solutionlist
  public :: solutiontype
  public :: solution_create
  public :: solution_list_init
  public :: solutiongrouplist
  public :: solutiongrouptype
  public :: nsolgps
  public :: neq,nja,ia,ja,amat,rhs,x,active

  type :: solutioncontainer
    class(solutiontype), pointer :: obj
  end type solutioncontainer

  type solutionlisttype
    integer :: nsolutions = 0
    type(solutioncontainer), allocatable, dimension(:) :: solutions
    contains
    procedure :: setsolution
    procedure :: getsolution
  end type solutionlisttype
  type(solutionlisttype) :: solutionlist

  integer, save, pointer :: neq
  integer, save, pointer :: nja
  integer, save, pointer, dimension(:), contiguous :: ia
  integer, save, pointer, dimension(:), contiguous :: ja
  double precision, save, pointer, dimension(:), contiguous :: amat
  double precision, save, pointer, dimension(:), contiguous :: rhs
  double precision, save, pointer, dimension(:), contiguous :: x
  integer, save, pointer, dimension(:), contiguous :: active
  type solutiontype
    integer :: id
    character(len=300) :: fname
    integer :: iu
    character(len=20) :: name
    integer, pointer :: neq
    integer, pointer :: nja
    integer, pointer, dimension(:), contiguous :: ia
    integer, pointer, dimension(:), contiguous :: ja
    double precision, pointer, dimension(:), contiguous :: amat
    double precision, pointer, dimension(:), contiguous :: rhs
    double precision, pointer, dimension(:), contiguous :: x
    integer, pointer, dimension(:), contiguous :: active
    type(SMS_DATA) :: sms
    type(XMD_DATA) :: xmd
    type(PCGU_DATA) :: pcgu
    type(modellisttype) :: modellist
    type(crosslisttype) :: crosslist
    type(sparsematrix) :: sparse
    contains
    procedure :: initialize
    procedure :: reset
    procedure :: addmodel
    procedure :: addcross
    procedure :: slnassigncrosses
    procedure :: connect
    procedure :: smsinit
    procedure :: solve
    procedure :: save
    procedure :: PNTSAV => solution_pntsav
    procedure :: PNTSET => solution_pntset
  end type solutiontype

  type :: solutiongrouptype
    integer :: id
    integer :: mxiter
    integer :: nsolutions
    integer,dimension(:),allocatable :: solutionidlist
  end type solutiongrouptype

  !number of solution groups (nsolgps) in the solutiongrouplist
  integer :: nsolgps
  type(solutiongrouptype),dimension(:),pointer :: solutiongrouplist

contains

  subroutine solution_list_init(nsolutions)
! ******************************************************************************
! solution_list_init -- Set size of solutionlist
! solutionlist is a list (actually an allocatable array) that contains all of
! the solutions that are part of the simulation.
! Subroutine: (1) allocate size of solutionlist%solutions
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    integer,intent(in) :: nsolutions
    allocate(solutionlist%solutions(nsolutions))
!
! -- return
    return
  end subroutine solution_list_init

  subroutine setsolution(this, newsolution, ipos)
! ******************************************************************************
! setsolution -- Set Solution
! For this solutionlist, set the pointer in position ipos to newsolution
! Subroutine: (1) point the solution to newsolution 
!             (2) increase this%nsolutions to reflect the number of stored 
!                 solutions
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(solutionlisttype) :: this
    class(solutiontype), target, intent(in) :: newsolution;
    integer, intent(in) :: ipos
! ------------------------------------------------------------------------------
!
    this%solutions(ipos)%obj => newsolution
    this%nsolutions=max(ipos,this%nsolutions)
!
! -- return
    return
  end subroutine setsolution

  subroutine getsolution(this, thesolution, ipos)
! ******************************************************************************
! getsolution -- Get Solution
! For this solutionlist, point the solution to the one in position ipos
! Subroutine: (1) point the solution to thesolution 
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(solutionlisttype) :: this
    class(solutiontype), pointer :: thesolution
    integer, intent(in) :: ipos
! ------------------------------------------------------------------------------
!
    thesolution => this%solutions(ipos)%obj
!
! -- return
    return
  end subroutine getsolution

  subroutine solution_create(filename,id)
! ******************************************************************************
! solution_create -- Create a New Solution 
! Using the data in filename,  assign this new solution an id number and store 
! the solution in the solutionlist.
! Subroutine: (1) allocate solution and assign id and name
!             (2) open the filename for later reading 
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
    type(solutiontype), pointer :: solution
! ------------------------------------------------------------------------------
!
! -- Create a new solution and add it to the solutionlist container
    allocate(solution)
    call solutionlist%setsolution(solution,id)
    solution%id=id
    write(solution%name,'(a4,i1)') 'SLN_',id
!
! -- Open solution input file for reading later after problem size is known
    call freeunitnumber(solution%iu)
    write(iout,'(/a,a)') ' Creating solution: ', solution%name
    WRITE(iout,'(a,a)') ' Using solution input file: ',trim(filename)
    WRITE(iout,'(a,i6)') ' On unit number: ',solution%iu
    open(unit=solution%iu,file=filename)
!
! -- return
    return
  end subroutine solution_create

  subroutine initialize(this)
! ******************************************************************************
! initialize -- Initialize This Solution
! Initialize the solution.  Must be called after the models and crosses have 
! been added to this solution.
! Subroutine: (1) Allocate neq and nja
!             (2) Assign model offsets and solution ids
!             (3) Allocate and initialize the solution arrays
!             (4) Point each model's x and rhs arrays
!             (5) Initialize the sparsematrix instance
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use ModelModule
    implicit none
    class(solutiontype) :: this
    class(modeltype),pointer :: mp
    integer :: i
    integer :: modelid
    integer, allocatable, dimension(:) :: rowmaxnnz
! ------------------------------------------------------------------------------
!
! -- Allocate and initialize neq and nja
    allocate(this%neq)
    allocate(this%nja)
    this%neq = 0
    this%nja = 0
    !calculate offsets
    do i=1,this%modellist%nmodels
        call this%modellist%getmodel(mp, i)
        mp%solutionid=this%id
        mp%offset = this%neq
        this%neq=this%neq+mp%neq
        mp%solutionid=this%id
    enddo
!
! -- Allocate and initialize solution arrays
    allocate(this%ia(this%neq+1))
    allocate(this%x(this%neq))
    allocate(this%rhs(this%neq))
    do i=1,this%neq
      this%x(i) = 0.0d0
    enddo
!
! -- Go through each model and point x and rhs to solution
    do i=1,this%modellist%nmodels
      call this%modellist%getmodel(mp, i)
      mp%x=>this%x(mp%offset+1:mp%offset+mp%neq)
      mp%rhs=>this%rhs(mp%offset+1:mp%offset+mp%neq)
    enddo
!
! -- Create the sparsematrix instance
    allocate(rowmaxnnz(this%neq))
    do i=1,this%neq
        rowmaxnnz(i)=4
    enddo
    call this%sparse%init(this%neq,this%neq,rowmaxnnz)
    deallocate(rowmaxnnz)
!
! -- return
    return
  end subroutine initialize

  subroutine addmodel(this, model, ipos)
! ******************************************************************************
! addmodel -- Add Model
! Subroutine: (1) add a model to this%modellist 
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(solutiontype) :: this
    class(modeltype),intent(in) :: model
    integer,intent(in) :: ipos
! ------------------------------------------------------------------------------
!
    call this%modellist%setmodel(model,ipos)
!
! -- return
    return
  end subroutine addmodel

  subroutine addcross(this, cross, ipos)
! ******************************************************************************
! addcross -- Add Cross
! Subroutine: (1) add a cross to this%crosslist
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(solutiontype) :: this
    type(crosstype),intent(in) :: cross
    integer,intent(in) :: ipos
! ------------------------------------------------------------------------------
!
    call this%crosslist%setcross(cross,ipos)
!
! -- return
    return
  end subroutine addcross

  subroutine slnassigncrosses(this)
! ******************************************************************************
! slnassigncross -- Assign crosses to this solution
! Subroutine: (1) count the number of crosses for this solution, 
!             (2) allocate this%crosslist,
!             (3) assign the appropriate crosses to this%crosslist
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(solutiontype) :: this
    class(crosstype), pointer :: c
    integer :: ncross,ic,ipos
! ------------------------------------------------------------------------------
!
! -- Count the number of crosses for this solution.   If both models for a cross
! -- belong to this solution, then this counts as only one cross.
    ncross=0
    do ic=1,crosslist%ncrosses
      call crosslist%getcross(c,ic)
      if(c%m1%solutionid==this%id) then
        ncross=ncross+1
        cycle
      elseif(c%m2%solutionid==this%id) then
        ncross=ncross+1
        cycle
      endif
    enddo
!
! -- Allocate the cross list for this solution
    allocate(this%crosslist%crosses(ncross))
!
! -- Now add each cross to this crosslist
    ipos=1
    do ic=1,crosslist%ncrosses
      call crosslist%getcross(c,ic)
      if(c%m1%solutionid==this%id) then
        call this%addcross(c,ipos)
        ipos=ipos+1
        cycle
      elseif(c%m2%solutionid==this%id) then
        call this%addcross(c,ipos)
        ipos=ipos+1
        cycle
      endif
    enddo
!
! -- return
    return
  end subroutine slnassigncrosses

  subroutine connect(this)
! ******************************************************************************
! connect -- Assign Connections
! Main workhorse method for solution.  This goes through all the models and all
! the connections and builds up the sparse matrix.
! Subroutine: (1) Add internal model connections, 
!             (2) Add cross terms,
!             (3) Allocate solution arrays
!             (4) Create mapping arrays
!             (5) Fill cross term values if necessary
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(solutiontype) :: this
    class(modeltype), pointer :: mp,m1,m2
    class(crosstype), pointer :: cp
    integer :: im,ic,i,j,jj,iglo,jglo,n,ierror,ipos
    integer :: istart,istop
! ------------------------------------------------------------------------------
!
! -- Add internal model connections to sparse
    do im=1,this%modellist%nmodels
      call this%modellist%getmodel(mp, im)
      do i=1,mp%neq
        do jj=mp%ia(i),mp%ia(i+1)-1
          j=mp%ja(jj)
          iglo=i+mp%offset
          jglo=j+mp%offset
          call this%sparse%addconnection(iglo,jglo,1)
        enddo
      enddo
    enddo
!
! -- Add the cross terms to sparse
    do ic=1,this%crosslist%ncrosses
      call this%crosslist%getcross(cp, ic)
      if(.not. cp%implicit) cycle
      do n=1,cp%ncross
        iglo=cp%nodem1(n)+cp%m1%offset
        jglo=cp%nodem2(n)+cp%m2%offset
        call this%sparse%addconnection(iglo,jglo,1)
        call this%sparse%addconnection(jglo,iglo,1)
      enddo
    enddo
!
! -- The number of non-zero array values are now known so
! -- ia and ja can be created from sparse. then destroy sparse
    this%nja=this%sparse%nnz
    allocate(this%ja(this%nja))
    allocate(this%amat(this%nja))
    allocate(this%active(this%neq))
    call this%sparse%sort()
    call this%sparse%filliaja(this%ia,this%ja,ierror)
    call this%sparse%destroy()
!
! -- fill active
    do n = 1, this%neq
      this%active(n) = 1
    end do
!
! -- Create mapping arrays for each model
    do im=1,this%modellist%nmodels
      call this%modellist%getmodel(mp, im)
      ipos=1
      do n=1,mp%neq
        istart=this%ia(n+mp%offset)
        istop=this%ia(n+mp%offset+1)-1
        do jj=istart,istop
          j=abs(this%ja(jj))
          if(j>mp%offset .and. j<=mp%offset+mp%neq) then
            mp%idxglo(ipos)=jj
            ipos=ipos+1
          endif
        enddo
      enddo
    enddo
!
! -- Create arrays for mapping cross connections to global solution
    do ic=1,this%crosslist%ncrosses
      call this%crosslist%getcross(cp, ic)
!
! -- If models are fully coupled in a single matrix
      if(cp%implicit) then
        do n=1,cp%ncross
          iglo=cp%nodem1(n)+cp%m1%offset
          jglo=cp%nodem2(n)+cp%m2%offset
          !find jglobal value in row iglo and store in idxglo
          do ipos=this%ia(iglo),this%ia(iglo+1)-1
            if(jglo==this%ja(ipos)) then
              cp%idxglo(n)=ipos
              exit
            endif
          enddo
          do ipos=this%ia(jglo),this%ia(jglo+1)-1
            if(iglo==this%ja(ipos)) then
              cp%idxsymglo(n)=ipos
              exit
            endif
          enddo
        enddo
      else
!
! -- Not fully coupled, so put the initial x values into the cross 
! -- package stage arrays
        call cp%crossfmfill()
      endif
    enddo    
!
! -- return
    return
  end subroutine connect

  subroutine smsinit(this)
! ******************************************************************************
! smsinit -- Initialize SMS
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  implicit none
  class(solutiontype) :: this
! ------------------------------------------------------------------------------
!
  call this%PNTSET()
  call SMS7U1AR(this%iu)
  call this%sms%PNTSAV()
  if (this%sms%linmeth.eq.1) then
    call this%xmd%PNTSAV()
  else if (this%sms%linmeth.eq.2) then
    call this%pcgu%PNTSAV()
  end if
  close(this%iu)
!
! -- return
    return
  end subroutine smsinit

  subroutine reset(this)
! ******************************************************************************
! reset -- Reset This Solution
! Reset this solution by setting amat and rhs to zero
! ******************************************************************************
! 
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  implicit none
  class(solutiontype) :: this
  integer :: i
! ------------------------------------------------------------------------------
!
  do i=1,this%nja
      this%amat(i)=0.
  enddo
  do i=1,this%neq
      this%rhs(i)=0.
  enddo
!
! -- return
    return
  end subroutine reset

subroutine solve(this)
    implicit none
    class(solutiontype) :: this
    class(modeltype), pointer :: mp
    class(crosstype), pointer :: cp
    integer :: im,ic,i,nodem1sln,nodem2sln,idiagsln
    integer :: kiter
    integer :: kstp,kper
    
    kstp=1
    kper=1
    
    do kiter=1,this%sms%mxiter
      !
      !set amat and rhs to zero
      call this%reset()
      !
      !modelfmcalc each model
      do im=1,this%modellist%nmodels
          call this%modellist%getmodel(mp, im)
          call mp%modelfmcalc()
      enddo
      !
      !fill each model
      do im=1,this%modellist%nmodels
          call this%modellist%getmodel(mp, im)
          call mp%modelfmfill(this%amat,this%nja)
      enddo
      !
      !add cross conductance to solution amat if implicit or
      !or call crossfmfill to update the packages with appropriate
      !model x values
      do ic=1,this%crosslist%ncrosses
          call this%crosslist%getcross(cp, ic)
          if(cp%implicit) then
              do i=1,cp%ncross
                  this%amat(cp%idxglo(i))=cp%cond(i)
                  this%amat(cp%idxsymglo(i))=cp%cond(i)
                  nodem1sln=cp%nodem1(i)+cp%m1%offset
                  nodem2sln=cp%nodem2(i)+cp%m2%offset
                  idiagsln=this%ia(nodem1sln)
                  this%amat(idiagsln)=this%amat(idiagsln)-cp%cond(i)
                  idiagsln=this%ia(nodem2sln)
                  this%amat(idiagsln)=this%amat(idiagsln)-cp%cond(i)
              enddo
          else
              call cp%crossfmfill()
          endif
      enddo
      !
      !call sms
      call this%PNTSET()
      call this%sms%PNTSET()
      if (this%sms%linmeth.eq.1) then
        call this%xmd%PNTSET()
      else if (this%sms%linmeth.eq.2) then
        call this%pcgu%PNTSET()
      end if
      CALL GLO2SMS1AP(kiter,kstp,kper)
    !check for convergence
    if (this%sms%icnvg.eq.1) exit

    end do
    
end subroutine solve

subroutine save(this,filename)
    implicit none
    class(solutiontype) :: this
    character(len=*), intent(in) :: filename
    integer :: inunit
    call freeunitnumber(inunit)
    open(unit=inunit,file=filename,status='unknown')
    write(inunit,*) 'ia'
    write(inunit,*) this%ia
    write(inunit,*) 'ja'
    write(inunit,*) this%ja
    write(inunit,*) 'amat'
    write(inunit,*) this%amat
    write(inunit,*) 'rhs'
    write(inunit,*) this%rhs
    write(inunit,*) 'x'
    write(inunit,*) this%x
    close(inunit)
end subroutine save

subroutine solution_pntsav(this)
    !use SOLUTION
    implicit none
    class(solutiontype) :: this
    this%neq=>neq
    this%nja=>nja
    this%ia=>ia
    this%ja=>ja
    this%amat=>amat
    this%rhs=>rhs
    this%x=>x
    this%active=>active
end subroutine solution_pntsav

subroutine solution_pntset(this)
    !use SOLUTION
    implicit none
    class(solutiontype) :: this
    neq=>this%neq
    nja=>this%nja
    ia=>this%ia
    ja=>this%ja
    amat=>this%amat
    rhs=>this%rhs
    x=>this%x
    active=>this%active
end subroutine solution_pntset

end module SOLUTIONMODULE



