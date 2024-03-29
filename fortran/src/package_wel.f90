module welmodule
use PackageModule
private
public :: wel_create
type, extends(packagetype) :: weltype
    double precision, allocatable, dimension(:) :: q
    contains
    procedure :: wel_ar
    procedure :: wel_allocate
    procedure :: packagerp=>welrp
    procedure :: fmcalc => welfmcalc
end type weltype

contains

subroutine wel_create(packobj,fname,id,inunit,iout)
    implicit none
    class(packagetype), pointer :: packobj
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    integer,intent(in) :: inunit
    integer,intent(in) :: iout
    type(weltype), pointer :: welobj
    integer :: maxbound
    allocate(welobj)
    packobj => welobj
    welobj%inunit=inunit
    welobj%iout=iout
    call welobj%wel_ar(fname,id)
end subroutine wel_create

subroutine wel_ar(this,fname,id)
    implicit none
    class(weltype) :: this
    character(len=*),intent(in) :: fname
    integer,intent(in) :: id
    integer :: maxbound
    write(this%iout,*) 'Creating well package: ', id
    this%filtyp='WEL'
    write(this%name,'(a,i1)') 'WEL_',id
    write(this%iout,*) 'opening well input file on unit: ', this%inunit
    open(unit=this%inunit,file=fname,status='old')
    read(this%inunit,*) maxbound
    !
    !allocate arrays in package superclass and in weltype
    call this%wel_allocate(maxbound)
end subroutine wel_ar

subroutine wel_allocate(this,maxbound)
    !allocate wel package of size maxbound
    implicit none
    class(weltype) :: this
    integer,intent(in) :: maxbound
    call this%pack_allocate(maxbound)
    allocate(this%q(maxbound))
end subroutine wel_allocate

subroutine welrp(this)
    implicit none
    class(weltype) :: this
    integer :: i
    read(this%inunit,*) this%nbound
    do i=1,this%nbound
        read(this%inunit,*) this%nodelist(i),this%q(i)
    enddo
end subroutine welrp

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

