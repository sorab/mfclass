subroutine initialize(simfile)
    use SimModule
    use SolutionModule
    use ModelModule
    use CrossModule
    implicit none
    character(len=*),intent(in) :: simfile
    character(len=300) :: line,fname,tag
    character(len=20) :: filtyp
    integer :: lloc,ityp1,ityp2,n,inunit,id,sid,mid,sgid
    integer :: nsolutions, nmodels, ncrosses, m1, m2,nsol,nmod,im
    real :: r
    class(solutiontype), pointer :: sp
    class(modeltype), pointer :: mp
    inunit=99
    iout=-1
    open(unit=inunit, file=simfile)
    do
        !read line and skip if necessary
        read(inunit,'(a)',end=100) line
        if(line.eq.' ') cycle
        if(line(1:1).eq.'#') cycle
        !
        !find number of models
        lloc=1
        CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
        filtyp=line(ityp1:ityp2)
        !
        !open list file - temporary solution???
        if (filtyp.eq.'LIST') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
            fname=line(ityp1:ityp2)
            call freeunitnumber(iout)
            open(unit=iout,file=fname,status='unknown')
        !
        !set nmodels
        elseif(filtyp.eq.'NMODELS') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,nmodels,R,IOUT,INUNIT)
            call model_list_init(nmodels)
        !
        !create models
        elseif(filtyp=='MODEL') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
            fname=line(ityp1:ityp2)
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,id,R,IOUT,INUNIT)
            call model_create(fname,id)
        !    
        !set nsolutions
        elseif(filtyp.eq.'NSOLUTIONS') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,nsolutions,R,IOUT,INUNIT)
            call solution_list_init(nsolutions)
        !
        !create solutions
        elseif(filtyp=='SOLUTION') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
            fname=line(ityp1:ityp2)
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,id,R,IOUT,INUNIT)
            call solution_create(fname,id)
        !    
        !set ncross
        elseif(filtyp.eq.'NCROSSES') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,ncrosses,R,IOUT,INUNIT)
            call cross_list_init(ncrosses)
        !
        !create solutions
        id=1
        elseif(filtyp=='XRS') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,0,N,R,IOUT,INUNIT)
            fname=line(ityp1:ityp2)
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,m1,R,IOUT,INUNIT)
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,m2,R,IOUT,INUNIT)
            call cross_create(fname,id,m1,m2)
            id=id+1
        !
        !solution groups
        elseif(filtyp=='SOLUTION_GROUP') then
            CALL URWORD(LINE,LLOC,ITYP1,ITYP2,2,sgid,R,IOUT,INUNIT)
            solutiongrouplist(1)%id=sgid
            read(inunit,*) fname,solutiongrouplist(sgid)%mxiter
            read(inunit,*) fname,nsol
            solutiongrouplist(sgid)%nsolutions=nsol
            allocate(solutiongrouplist(sgid)%solutions(nsol))
            do n=1,nsol
                read(inunit,*) fname,sid
                call solutionlist%getsolution(sp,sid)
                read(inunit,*) fname,nmod
                allocate(sp%modellist%models(nmod))
                do im=1,nmod
                    read(inunit,*) fname,mid
                    call modellist%getmodel(mp,mid)
                    call sp%addmodel(mp,im)
                enddo
            enddo
        endif
    cycle
100 exit        
    enddo
    close(inunit)
end subroutine initialize

subroutine freeunitnumber(iu)
    implicit none
    integer :: lastunitnumber
    parameter(lastunitnumber=10000)
    integer,intent(inout) :: iu
    integer, save :: nextunitnumber=1000
    integer i
    logical :: opened
    do i=nextunitnumber,lastunitnumber
        inquire(unit=i,opened=opened)
        if(.not. opened) exit
    enddo
    nextunitnumber=i
    iu=nextunitnumber
end subroutine freeunitnumber

program mf2015fort
    use SolutionModule
    use ModelModule
    use CrossModule
    implicit none
    class(modeltype), pointer :: m
    class(solutiontype), pointer :: s
    class(crosstype), pointer :: c
    integer :: is,im,nsg,ic,ipos
    print *, 'starting...'

    !read the name file and initialize instances
    call initialize('simulation.nam')
    nsg=1
    !
    !add a cross
    !allocate(s%crosslist%crosses(1))
    !call crosslist%getcross(c,1)
    !call c%initialize()
    !call s%addcross(c,1)
    !
    !initialize each solution
    do is=1,solutionlist%nsolutions
        call solutionlist%getsolution(s,is)
        call s%initialize()
    enddo
    !
    !initialize each cross
    do ic=1,crosslist%ncrosses
        call crosslist%getcross(c,ic)
        call c%initialize()
    enddo
    !
    !assign crosses to solutions
    do is=1,solutionlist%nsolutions
        call solutionlist%getsolution(s,is)
        allocate(s%crosslist%crosses(1))
        ipos=1
        do ic=1,crosslist%ncrosses
            call crosslist%getcross(c,ic)
            if(c%m1%solutionid==is) then
                call s%addcross(c,ipos)
                ipos=ipos+1
                cycle
            elseif(c%m2%solutionid==is) then
                call s%addcross(c,ipos)
                ipos=ipos+1
                cycle
            endif
        enddo
    enddo
    !
    !connect each solution
    do is=1,solutionlist%nsolutions
        call solutionlist%getsolution(s,is)
        call s%connect()
        call s%smsinit()
    enddo
    !
    !read stress period data
    do is=1,solutionlist%nsolutions
        call solutionlist%getsolution(s,is)
        do im=1,s%modellist%nmodels
            call s%modellist%getmodel(m,im)
            call m%read_prepare()
        enddo
    enddo
    !
    !solve each solution
    do is=1,solutionlist%nsolutions
        call solutionlist%getsolution(s,is)
        call s%solve()
    enddo
    !
    !save each solution
    do is=1,solutionlist%nsolutions
        call solutionlist%getsolution(s,is)
        call s%save('output.dat')
    enddo
    
    print *, 'ending...'

end program mf2015fort
