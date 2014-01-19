
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
  !  integer :: mxiter
  !  double precision :: xtol
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
    !integer :: nmodels=0
    !integer, dimension(10) :: imodels
      type(modellisttype) :: modellist
      type(crosslisttype) :: crosslist
    !integer :: ncrosses=0
    !integer, dimension(10) :: icrosses
    type(sparsematrix) :: sparse
    contains
    procedure :: initialize
    procedure :: reset
    procedure :: addmodel
    procedure :: addcross
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
!    double precision :: xtol
    integer :: nsolutions
    type(solutiontype),allocatable,dimension(:) :: solutions
end type solutiongrouptype
type(solutiongrouptype),dimension(1) :: solutiongrouplist

contains

subroutine setsolution(this, newsolution, ipos)
    implicit none
    class(solutionlisttype) :: this
    class(solutiontype), target, intent(in) :: newsolution;
    integer, intent(in) :: ipos
    this%solutions(ipos)%obj => newsolution
    this%nsolutions=max(ipos,this%nsolutions)
end subroutine setsolution

subroutine getsolution(this, thesolution, ipos)
    implicit none
    class(solutionlisttype) :: this
    class(solutiontype), pointer :: thesolution
    integer, intent(in) :: ipos
    thesolution => this%solutions(ipos)%obj
end subroutine getsolution

subroutine solution_list_init(nsolutions)
    implicit none
    integer,intent(in) :: nsolutions
    allocate(solutionlist%solutions(nsolutions))
end subroutine solution_list_init

subroutine solution_create(filename,id)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id
    character(len=300) :: line,fname
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,iout,inunit,npackages,ipak
    real :: r
    type(solutiontype), pointer :: solution
    !
    print *,'Creating solution: ', id
    !
    !create a new solution and add it to the solutionlist container
    allocate(solution)
    call solutionlist%setsolution(solution,id)
    solution%id=id
    !
    !open and read the solution input file
    call freeunitnumber(solution%iu)
    print *, 'opening solution namefile on unit: ', solution%iu
    open(unit=solution%iu,file=filename)
!    do
!        !read line and skip if necessary
!        read(solution%iu,'(a)',end=100) line
!        if(line.eq.' ') cycle
!        if(line(1:1).eq.'#') cycle
!        !
!        !decode
!        lloc=1
!        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
!        filtyp=line(ityp1:ityp2)
!        !
!        if(filtyp.eq.'MXITER') then
!            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,solution%mxiter,R,IOUT,INUNIT)
!        elseif(filtyp.eq.'XTOL') then
!            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,3,N,R,IOUT,INUNIT)
!            solution%xtol=r
!        else
!            print *, 'unknown file type in solution file: ', filtyp
!        endif
!    cycle
!100 exit        
!    enddo
!    close(solution%iu)
end subroutine solution_create

subroutine initialize(this)
    !this subroutine must be called after models and crosses
    !have been added
    use ModelModule
    implicit none
    class(solutiontype) :: this
    class(modeltype),pointer :: mp
    integer :: i
    integer :: modelid
    integer, allocatable, dimension(:) :: rowmaxnnz
    !allocate and initialize neq and nja
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
    !allocate and initialize solution arrays
    allocate(this%ia(this%neq+1))
    allocate(this%x(this%neq))
    allocate(this%rhs(this%neq))
    do i=1,this%neq
        this%x(i) = 0.0d0
    enddo
    !
    !go through each model and point x and rhs to solution
    do i=1,this%modellist%nmodels
        call this%modellist%getmodel(mp, i)
        mp%x=>this%x(mp%offset+1:mp%offset+mp%neq)
        mp%rhs=>this%rhs(mp%offset+1:mp%offset+mp%neq)
    enddo
    !
    !create the sparsematrix instance
    allocate(rowmaxnnz(this%neq))
    do i=1,this%neq
        rowmaxnnz(i)=4
    enddo
    call this%sparse%init(this%neq,this%neq,rowmaxnnz)
    deallocate(rowmaxnnz)
end subroutine initialize

subroutine addmodel(this, model, ipos)
    implicit none
    class(solutiontype) :: this
    class(modeltype),intent(in) :: model
    integer,intent(in) :: ipos
    call this%modellist%setmodel(model,ipos)
end subroutine addmodel

subroutine addcross(this, cross, ipos)
    implicit none
    class(solutiontype) :: this
    type(crosstype),intent(in) :: cross
    integer,intent(in) :: ipos
    call this%crosslist%setcross(cross,ipos)
end subroutine addcross

subroutine connect(this)
    implicit none
    class(solutiontype) :: this
    class(modeltype), pointer :: mp,m1,m2
    class(crosstype), pointer :: cp
    integer :: im,ic,i,j,jj,iglo,jglo,n,ierror,ipos
    integer :: istart,istop
    !
    !add internal model connections
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
    !add the cross terms
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
    !the number of non-zero array values are now known so
    !ia and ja can be created from sparse. then destroy sparse
    this%nja=this%sparse%nnz
    allocate(this%ja(this%nja))
    allocate(this%amat(this%nja))
    allocate(this%active(this%neq))
    call this%sparse%sort()
    call this%sparse%filliaja(this%ia,this%ja,ierror)
    call this%sparse%destroy()
    !fill active
    do n = 1, this%neq
      this%active(n) = 1
    end do
    !
    !create mapping arrays for each model
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
    !create arrays for mapping cross connections to global solution
    do ic=1,this%crosslist%ncrosses
        call this%crosslist%getcross(cp, ic)
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
            !put the initial x values into the cross package stage arrays
            call cp%fill()
        endif
    enddo    
end subroutine connect

subroutine smsinit(this)
    implicit none
    class(solutiontype) :: this
    call this%PNTSET()
    call SMS7U1AR(this%iu)
    call this%sms%PNTSAV()
    if (this%sms%linmeth.eq.1) then
      call this%xmd%PNTSAV()
    else if (this%sms%linmeth.eq.2) then
      call this%pcgu%PNTSAV()
    end if
    close(this%iu)
end subroutine smsinit


subroutine reset(this)
    implicit none
    class(solutiontype) :: this
    integer :: i
    do i=1,this%nja
        this%amat(i)=0.
    enddo
    do i=1,this%neq
        this%rhs(i)=0.
    enddo
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
      !add cross conductance to solution amat
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
              call cp%fill()
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



