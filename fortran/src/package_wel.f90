module welmodule
use PackageModule
private
public :: wel_create
type, extends(packagetype) :: weltype
    double precision, allocatable, dimension(:) :: q
    contains
    procedure :: wel_ar
    procedure :: read_next
    procedure :: fmcalc => welfmcalc
end type weltype

contains

subroutine wel_create(packobj,fname,id)
    implicit none
    class(packagetype), pointer :: packobj
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    type(weltype), pointer :: welobj
    integer :: maxbound
    allocate(welobj)
    packobj => welobj
    call welobj%wel_ar(fname,id)
end subroutine wel_create

subroutine wel_ar(this,fname,id)
    implicit none
    class(weltype) :: this
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    integer :: maxbound

    print *,'Creating well package: ', id
    call freeunitnumber(this%inunit)
    print *, 'opening model namefile on unit: ', this%inunit
    open(unit=this%inunit,file=fname,status='old')
    read(this%inunit,*) maxbound
    !
    !allocate package members
    call this%pack_allocate(maxbound)
    !
    !allocate well package members
    allocate(this%q(maxbound))
end subroutine wel_ar

subroutine read_next(this)
    implicit none
    class(weltype) :: this
    integer :: i
    read(this%inunit,*) this%nbound
    do i=1,this%nbound
        read(this%inunit,*) this%nodelist(i),this%q(i)
    enddo
end subroutine read_next

subroutine welfmcalc(this)
    implicit none
    class(weltype) :: this
    integer :: i
    do i=1,this%nbound
        this%hcof(i) = 0.
        this%rhs(i) = -this%q(i)
    enddo
end subroutine welfmcalc

end module welmodule

