module SolutionModule
use ModelModule
use CrossModule
use sparsemodule, only:sparsematrix
implicit none
private
public :: solutionlist
public :: solutiontype
public :: solution_create
public :: solution_list_init
public :: getsolution

type :: solutioncontainer
    class(solutiontype), pointer :: obj
end type solutioncontainer

type solutionlisttype
    integer :: numsolutions = 0
    type(solutioncontainer), allocatable, dimension(:) :: solutions
    contains
    procedure :: setsolution
    procedure :: getsolution
end type solutionlisttype
type(solutionlisttype) :: solutionlist

type solutiontype
  integer :: id
  character(len=20) :: name
  integer :: mxiter
  double precision :: xtol
  integer :: neq=0
  integer :: nja=0
  integer, allocatable, dimension(:) :: ia
  integer, allocatable, dimension(:) :: ja
  double precision, allocatable, dimension(:) :: amat
  double precision, pointer, dimension(:) :: rhs
  double precision, pointer, dimension(:) :: x
  !integer :: nmodels=0
  !integer, dimension(10) :: imodels
    type(modellisttype) :: modellist
    type(crosslisttype) :: crosslist
  integer :: ncrosses=0
  integer, dimension(10) :: icrosses
  type(sparsematrix) :: sparse
  contains
  procedure :: initialize
  procedure :: addmodel
  procedure :: addcross
  procedure :: connect
end type solutiontype

!type solutionlist
!  integer :: nval=0
!  type(solutiontype), dimension(10) :: sl
!end type solutionlist

contains

subroutine setsolution(this, newsolution, ipos)
    implicit none
    class(solutionlisttype) :: this
    class(solutiontype), target, intent(in) :: newsolution;
    integer, intent(in) :: ipos
    this%solutions(ipos)%obj => newsolution
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

subroutine initialize(self)
    !this subroutine must be called after models and crosses
    !have been added
    use ModelModule
    implicit none
    class(solutiontype) :: self
    class(modeltype),pointer :: mp
    integer :: i
    integer :: modelid
    integer, allocatable, dimension(:) :: rowmaxnnz
    do i=1,self%modellist%nmodels
        call modellist%getmodel(mp, i)
        mp%offset = self%neq
        self%neq=self%neq+mp%neq
        mp%solutionid=self%id
    enddo
    !
    !allocate solution arrays
    allocate(self%ia(self%neq+1))
    allocate(self%x(self%neq))
    allocate(self%rhs(self%neq))
    !
    !go through each model and point x and rhs to solution
    do i=1,self%modellist%nmodels
        call modellist%getmodel(mp, i)
        mp%x=>self%x(mp%offset+1:mp%offset+mp%neq)
        mp%rhs=>self%rhs(mp%offset+1:mp%offset+mp%neq)
    enddo
    !
    !create the sparsematrix instance
    allocate(rowmaxnnz(self%neq))
    do i=1,self%neq
        rowmaxnnz(i)=4
    enddo
    call self%sparse%init(self%neq,self%neq,rowmaxnnz)
    deallocate(rowmaxnnz)
end subroutine initialize

subroutine addmodel(this, model, ipos)
    implicit none
    class(solutiontype) :: this
    type(modeltype),intent(in) :: model
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

subroutine solution_create(filename,id)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in) :: id
    character(len=300) :: line,fname
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,iout,inunit,npackages,ipak
    real :: r
    type(solutiontype), pointer :: solution

    !create a new solution and add it to the solutionlist container
    allocate(solution)
    call solutionlist%setsolution(solution,id)

    print *,'Creating solution: ', id
    call freeunitnumber(inunit)
    print *, 'opening solution namefile on unit: ', inunit
    
    open(unit=inunit,file=filename)
    do
        !read line and skip if necessary
        read(inunit,'(a)',end=100) line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        !
        !decode
        lloc=1
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
        filtyp=line(ityp1:ityp2)
        !
        if(filtyp.eq.'MXITER') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,solution%mxiter,R,IOUT,INUNIT)
        elseif(filtyp.eq.'XTOL') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,3,N,R,IOUT,INUNIT)
            solution%xtol=r
        else
            print *, 'unknown file type in solution file: ', filtyp
        endif
    cycle
100 exit        
    enddo
    close(inunit)
end subroutine solution_create

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
        call modellist%getmodel(mp, im)
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
        call crosslist%getcross(cp, ic)
        if(.not. cp%implicit) cycle
        do n=1,cp%ncross
            iglo=cp%nodem1(n)+cp%m1%offset
            jglo=cp%nodem2(n)+cp%m2%offset
            call this%sparse%addconnection(iglo,jglo,1)
            call this%sparse%addconnection(jglo,iglo,1)
        enddo
    enddo
    !
    !the number of non-zero array values are now known
    this%nja=this%sparse%nnz
    allocate(this%ja(this%nja))
    allocate(this%amat(this%nja))
    call this%sparse%sort()
    call this%sparse%filliaja(this%ia,this%ja,ierror)
    !
    !create mapping arrays for each model
    do im=1,this%modellist%nmodels
        call modellist%getmodel(mp, im)
        ipos=0
        do n=1,mp%neq
            istart=this%ia(n+mp.offset)
            istop=this%ia(n+mp.offset+1)-1
            do jj=istart,istop
                j=abs(this%ja(jj))
                if(mp%offset<=j<mp%offset+mp%neq) then
                    mp%idxglo(ipos)=jj
                    ipos=ipos+1
                endif
            enddo
        enddo
    enddo

STOPPED HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

end subroutine connect

end module SolutionModule



