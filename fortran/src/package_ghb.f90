module ghbmodule
use PackageModule
private
public :: ghb_create
public :: ghbtype
type, extends(packagetype) :: ghbtype
    double precision, pointer, dimension(:) :: cond
    double precision, allocatable, dimension(:) :: stage
    contains
    procedure :: ghb_ar
    procedure :: ghb_allocate
    procedure :: packagerp=>ghbrp
    procedure :: fmcalc => ghbfmcalc
end type ghbtype

contains

subroutine ghb_create(packobj,fname,id,inunit,iout)
    !create ghb object from fname
    implicit none
    class(packagetype), pointer :: packobj
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    integer,intent(in) :: inunit
    integer,intent(in) :: iout
    type(ghbtype), pointer :: ghbobj
    allocate(ghbobj)
    packobj => ghbobj
    ghbobj%inunit=inunit
    ghbobj%iout=iout
    call ghbobj%ghb_ar(fname,id)
end subroutine ghb_create

subroutine ghb_ar(this,fname,id)
    implicit none
    class(ghbtype) :: this
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    integer :: maxbound
    write(this%iout,*) 'Creating ghb package: ', id
    this%filtyp='GHB'
    write(this%name,'(a,i1)') 'GHB_',id
    write(this%iout,*)  'opening ghb input file on unit: ', this%inunit
    open(unit=this%inunit,file=fname)
    read(this%inunit,*) maxbound
    !
    !allocate arrays in package superclass and in ghbtype
    call this%ghb_allocate(maxbound)
end subroutine ghb_ar

subroutine ghb_allocate(this,maxbound)
    !allocate ghb package of size maxbound
    implicit none
    class(ghbtype) :: this
    integer,intent(in) :: maxbound
    call this%pack_allocate(maxbound)
    allocate(this%cond(maxbound))
    allocate(this%stage(maxbound))
end subroutine ghb_allocate

subroutine ghbrp(this)
    implicit none
    class(ghbtype) :: this
    integer :: i
    if(this%inunit<=0) return
    read(this%inunit,*) this%nbound
    do i=1,this%nbound
        read(this%inunit,*) this%nodelist(i),this%stage(i),this%cond(i)
    enddo
end subroutine ghbrp

subroutine ghbfmcalc(this)
    implicit none
    class(ghbtype) :: this
    integer :: i
    do i=1,this%nbound
        this%hcof(i) = -this%cond(i)
        this%rhs(i) = -this%cond(i) * this%stage(i)
    enddo
end subroutine ghbfmcalc

end module ghbmodule
