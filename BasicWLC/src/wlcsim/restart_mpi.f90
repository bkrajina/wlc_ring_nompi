subroutine pt_restart(mc,md)
! Takes wlcsim_params and wlcsim_data and restarts the MPI workers for running
! parallel-tempered MC simulations.
!
! This function takes the place of PT_override in the case of restart
! This will read from a output directory and restart multiple replicas
! Override initialization with parallel setup parameters
!  In particualar it changes: mc%AB, mc%rep, mc%mu, mc%repSufix
    use mpi
    use simMod
    Implicit none
    type(MCvar), intent(inout) :: mc
    type(MCData), intent(inout) :: md
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) id, nThreads,ierror
    integer (kind=4) error  ! error id for MIP functions
    character*64 iostrg    ! for file naming
    character*16 vNum    ! for file naming
    character*64 dir
    integer ( kind = 4 ) status(MPI_status_SIZE) ! MPI stuff
    integer, parameter :: nTerms=8  ! number of energy terms
    real(dp) mag ! magnitude for renormalizing U
    real(dp) cof(nTerms)
    integer I ! bead index

    ! file parsing
    integer ios ! read status (detect end of file)
    real(dp) temp(28)  ! values

    ! Which replica am I?
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    source=0;
    call MPI_Recv ( mc%rep, 1, MPI_integer, source, 0, &
                   MPI_COMM_WORLD, status, error )
    if (mc%rep.ne.id) then
        print*, "That's not what I expected! see restart"
    endif
    if (nThreads.lt.3) then
        print*, "don't use pt_restart for fewer than 3 treads"
        stop 1
    endif
    write(vNum,'(I4)') mc%rep
    vNum=adJustL(vNum)
    vNum="v"//trim(vNum)

    ! Where to read from
    dir="data/"

    ! read Some operation variables
    ! Many of these aren't necessary but a few are
    iostrg=trim(dir)//"out1"
    iostrg=trim(iostrg)//trim(vNum)
    print*, "reading", iostrg
    open(unit=1, file=iostrg, status ='OLD')
    read(1,*)
    do WHILE (.TRUE.)
        read(1,*,IOSTAT=ios), temp(1), temp(2), temp(3), temp(4), temp(5), &
                              temp(6), temp(7), temp(8), temp(9), temp(10), &
                              temp(11), temp(12), temp(13), temp(14), &
                              temp(15), temp(16)
        if (ios.eq.0) then
            mc%ind=nint(temp(1))
            mc%EElas(1)=temp(3)
            mc%EElas(2)=temp(4)
            mc%EElas(3)=temp(5)
            mc%ECouple=temp(6)
            mc%EKap=temp(7)
            mc%EChi=temp(8)
            mc%EField=temp(9)
            mc%ebind=temp(10)
            mc%M=temp(11)
            mc%HP1_Bind=temp(12)
            mc%chi=temp(13)
            mc%mu=temp(14)
            mc%Kap=temp(15)
            mc%h_A=temp(16)
        else
            Exit
        endif
    enddo
    close(1)
    print*, "first set from file", iostrg
    print*, temp
    ! not sure if the following if statments are necessary
    if (mc%Chi.ne.0.0) then
        mc%x_Chi=mc%EChi/mc%Chi
    endif
    if (mc%Chi.ne.0.0) then
        mc%x_Couple=mc%ECouple/mc%HP1_Bind
    endif
    if (mc%Kap.ne.0) then
        mc%x_Kap=mc%EKap/mc%Kap
    endif
    if (mc%x_Field.ne.0.0) then
        mc%x_Field=mc%EField/mc%h_A
    endif
    if (mc%Mu.ne.0.0) then
        mc%x_Mu=mc%ebind/mc%Mu
    endif

    ! read back in addaptation stuff, May make slight difference
    iostrg=trim(dir)//"out3"
    iostrg=trim(iostrg)//trim(vNum)
    print*, iostrg
    open(unit=1, file=iostrg, status ='OLD')
    read(1,*)
    do WHILE (.TRUE.)
        read(1,*,IOSTAT=ios), temp(1), temp(2), temp(3), temp(4), temp(5), &
                              temp(6), temp(7), temp(8), temp(9), temp(10), &
                              temp(11), temp(12), temp(13), temp(14), &
                              temp(15), temp(16), temp(17), temp(18), &
                              temp(19), temp(20), temp(21), temp(22), &
                              temp(23), temp(24), temp(25), temp(26), &
                              temp(27), temp(28)
        if (ios.eq.0) then
            mc%WindoW(1)=temp(3); mc%MCAMP(1)=temp(4); mc%PHIT(1)=temp(5);
            mc%WindoW(2)=temp(6); mc%MCAMP(2)=temp(7); mc%PHIT(2)=temp(8);
            mc%WindoW(3)=temp(9); mc%MCAMP(3)=temp(10); mc%PHIT(3)=temp(11);
            mc%MOVEON(4)=nint(temp(12)); mc%MCAMP(4)=temp(13); mc%PHIT(4)=temp(14);
            mc%MOVEON(5)=nint(temp(15)); mc%MCAMP(5)=temp(16); mc%PHIT(5)=temp(17);
            mc%MOVEON(6)=nint(temp(18)); mc%MCAMP(6)=temp(19); mc%PHIT(6)=temp(20);
            mc%MOVEON(7)=nint(temp(21)); mc%PHIT(7)=temp(22);
            mc%MOVEON(8)=nint(temp(23)); mc%PHIT(8)=temp(24);
            mc%MOVEON(9)=nint(temp(25)); mc%PHIT(9)=temp(26);
            mc%MOVEON(10)=nint(temp(27)); mc%PHIT(10)=temp(28)
        else
            Exit
        endif
    enddo
    close(1)
    print*, "second set from file", iostrg
    print*, temp



    ! read R and AB from file
    write(iostrg,"(I8)"), mc%ind
    iostrg=adjustL(iostrg)
    iostrg="r"//trim(iostrg)
    iostrg=trim(dir)//trim(iostrg)
    iostrg=trim(iostrg)//trim(vNum)
    print*, "reading", iostrg
    open (unit = 5, file = iostrg, status = 'OLD')
    print*, "NT=",mc%NT
    ios=0;
    do I=1,mc%NT
       if (ios.ne.0) then
           print*, "Problem while reading R, Possible incomplete file"
           stop 1
       endif
       read(5,*) md%R(I,1),md%R(I,2),md%R(I,3),md%AB(I)
    enddo
    close(5)

    ! read U
    write(iostrg,"(I8)"), mc%ind
    iostrg=adjustL(iostrg)
    iostrg="u"//trim(iostrg)
    iostrg=trim(dir)//trim(iostrg)
    iostrg=trim(iostrg)//trim(vNum)
    ! read U from file
    open (unit = 5, file = iostrg, status = 'OLD')
    do I=1,mc%NT
       read(5,*) md%U(I,1),md%U(I,2),md%U(I,3)
       mag=sqrt(md%U(I,1)**2+md%U(I,2)**2+md%U(I,3)**2)
       md%U(I,1)=md%U(I,1)/mag
       md%U(I,2)=md%U(I,2)/mag
       md%U(I,3)=md%U(I,3)/mag
    enddo
    close(5)

    ! Let head node know what cof values you read
    cof(1)=mc%chi
    cof(2)=mc%mu
    cof(3)=mc%h_A
    cof(4)=mc%HP1_Bind
    cof(5)=mc%KAP
    cof(6)=mc%para(1)
    cof(7)=mc%para(2)
    cof(8)=mc%para(3)
    dest=0
    call MPI_Send (cof,nTerms, MPI_doUBLE_PRECISION, dest,   0, &
                    MPI_COMM_WORLD,error )

    ! Make repsufix
    write(iostrg,"(I4)"), mc%rep
    iostrg=adjustL(iostrg)
    iostrg=trim(iostrg)
    iostrg="v"//trim(iostrg)
    iostrg=trim(iostrg)
    mc%repSufix=trim(iostrg)

    ! keep track of which thread you are
    mc%id=int(id)
end subroutine
