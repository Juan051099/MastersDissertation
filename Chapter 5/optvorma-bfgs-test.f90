program optvormabfgs

    use iso_c_binding, only: c_loc
    use optvormod2, only: pdata_type,evalfg,drawsol,inA,drawsolJ1,drawsolJ2,drawsolJ3, &
    evalsizeedJi2,evalangJi3,evalvolJi1,evalJg

    implicit none

    ! LOCAL SCALARS
    integer :: allocerr,i,j,inform,maxnv,n,nsites, &
    iprint,m,flag,counter,activate
    real(kind=8) :: f,seed,factr,start,finish, &
    pgtol,epsx,scale,area,marg,G1,j2,j3,j1,j4
    logical :: scaled = .true., printx = .false.

    type(pdata_type), target :: pdata

    ! LOCAL ARRAYS
    logical :: lsave(4) 
    integer :: isave(44)
    integer, allocatable :: nbd(:),iwa(:)
    real(kind=8), allocatable :: g(:),x0(:),x(:),xold(:),sites(:,:),l(:),u(:),wa(:),ji(:),ji2(:),avg(:),&
    sizes(:,:),avgang(:),ji3(:),angles(:,:),vol(:),ji1(:)
    character(len=80) :: filename,task,csave
    real(kind=8) :: dsave(29),p(2)

    ! FUNCTIONS
    real(kind=8), external :: drand

    ! Set sites

    write(*,*) 'Enter nsites>0: '
    read(*,*) nsites

    ! Set disjoint convex polygons Aj whose union define A

    pdata%npols = 1

    allocate(pdata%nvpols(pdata%npols),stat=allocerr)
    if ( allocerr .ne. 0) then
        write(*,*) 'Allocation error.'
        stop
    end if

    pdata%nvpols(1:pdata%npols) = 5

    maxnv = maxval( pdata%nvpols(1:pdata%npols))

    allocate(pdata%poledges(pdata%npols,maxnv),pdata%xpol(pdata%npols,maxnv), &
            pdata%ypol(pdata%npols,maxnv),stat=allocerr)
    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error.'
        stop
    end if

    pdata%xpol(1,1:pdata%nvpols(1)) = (/ 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0 /)
    pdata%ypol(1,1:pdata%nvpols(1)) = (/ 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0 /)

    pdata%poledges(1,1:pdata%nvpols(1)-1) = (/ 0, 0, 0, 0 /) ! 0 means frontier of A, 1 means internal

    scale = 1.0d0
    pdata%volA = 1.0d0
    pdata%drawscale = 6.0d0
    if ( scaled ) then 

        scale = sqrt( real(nsites) )
        area = scale ** 2
        write(*,*) 'area = ',area
        write(*,*) 'scale = ',scale
        pdata%volA = area
        pdata%drawscale = 1.0d0/scale

    end if

    pdata%xpol(1,1:pdata%nvpols(1)) = scale * pdata%xpol(1,1:pdata%nvpols(1))
    pdata%ypol(1,1:pdata%nvpols(1)) = scale * pdata%ypol(1,1:pdata%nvpols(1))

    pdata%xmin = minval( pdata%xpol(1,1:pdata%nvpols(1)) )
    pdata%xmax = maxval( pdata%xpol(1,1:pdata%nvpols(1)) )
    pdata%ymin = minval( pdata%ypol(1,1:pdata%nvpols(1)) )
    pdata%ymax = maxval( pdata%ypol(1,1:pdata%nvpols(1)) )

    marg = 0.1d0 * max( 1.0d0, abs( pdata%xmin ), abs( pdata%xmax ), abs( pdata%ymin ), abs( pdata%ymax ) )

    pdata%xmin = pdata%xmin - marg
    pdata%xmax = pdata%xmax + marg
    pdata%ymin = pdata%ymin - marg
    pdata%ymax = pdata%ymax + marg

    ! Goal function data

    pdata%Jweight(0) = 1.0d0 
    pdata%Jweight(1:4) = 0.0d0
    ! pdata%Jweight(1) = 1.0d0
    ! pdata%Jweight(2) = 1.0d0
    ! pdata%Jweight(3) = 1.0d0

    ! write(*,*) 'Enter ctol for J2 (between 0 and 1): '
    ! read(*,*) pdata%J2ctol
    ! write(*,*) 'Enter ctol for J3 (between 0 and 1): '
    ! read(*,*) pdata%J3ctol

    n = 2 * nsites

    allocate(x(n),g(n),x0(n),xold(n),sites(2,nsites),ji(n/2),ji2(n/2),avg(n/2),sizes(n/2,10),ji3(n/2),&
    avgang(n/2),angles(n/2,10),ji1(n/2),vol(n/2),stat=allocerr)
    if ( allocerr .ne. 0) then
        write(*,*) 'Allocation error.'
        stop
    end if

    
    seed = 123456.0d0
    if ( .not. scaled ) then
        do i=1,nsites
            sites(1,i) = 0.10d0 + 0.80d0 * drand(seed)
            sites(2,i) = 0.10d0 + 0.80d0 * drand(seed)

            x(2*i-1) = sites(1,i)
            x(2*i) = sites(2,i)
        end do
    else
        do i = 1,nsites
            p(1) = pdata%xmin + (pdata%xmax - pdata%xmin) * drand(seed)
            p(2) = pdata%ymin + (pdata%ymax - pdata%ymin) * drand(seed)

            do while ( .not. InA(p,c_loc(pdata)) )
                p(1) = pdata%xmin + ( pdata%xmax - pdata%xmin ) * drand(seed)
                p(2) = pdata%ymin + ( pdata%ymax - pdata%ymin ) * drand(seed)
            end do

            x(2*i-1) = p(1)
            x(2*i) = p(2)
        end do
    end if
    x0(1:n) = x(1:n)
    write(*,*) 'x0=',x0(1:n)
    xold(1:n) = x(1:n)
    ! write(*,*) 'x=', x(1:n)

    ! Optimize
    m = 5
    iprint = 0
    factr = 0 !1.0d+07
    pgtol = 1.0d-08
    epsx = 0.0d0 !1.0d-06 

    allocate(nbd(n),l(n),u(n),stat=allocerr)
    if ( allocerr .ne. 0) then
        write(*,*) 'Allocation error.'
        stop
    end if
    allocate(iwa(3*n),stat=allocerr)
    if ( allocerr .ne. 0) then
        write(*,*) 'Allocation error.'
        stop
    end if
    allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m), stat=allocerr)
    if ( allocerr .ne. 0) then
        write(*,*) 'Allocation error.'
        stop
    end if
    do i=1, n
        nbd(i) = 2
        l(i) = pdata%xpol(1,1)
        u(i) = pdata%xpol(1,2)
    end do

!   We start the iteration by initializing task.
 
    task = 'START'

!   The beginning of the loop
    
    call cpu_time(start)

    do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
            task.eq.'START') 
            
!   This is the call to the L-BFGS-B code.
            
        call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
                 wa, iwa, task, iprint,&
                csave, lsave, isave, dsave )
            
        if (task(1:2) .eq. 'FG') then
!   Compute f and the gradient g

            call evalfg(n,x,f,g,inform,c_loc(pdata))

        else if( task(1:5) .eq. 'NEW_X' ) then

            if ( maxval( abs( xold(1:n)-x(1:n) ) ) .le. epsx ) then 
                write(*,*) 'x and xold are very close'

                goto 100
            end if

            xold(1:n) = x(1:n)
        end if

    end do

    100 continue
    call cpu_time(finish)

    write(*,*) 'n=',nsites,'f=',f,'gnorm=',dsave(13),'it=',isave(30),'fcnt=',isave(34),'time=',finish-start

    call evalJg(n,x,G1,j2,j3,j1,j4,ji,g,flag,c_loc(pdata))

    write(*,*) 'G = ', G1, ' J2 = ', j2

    if ( .not. scaled ) then
        write(filename,"(A13,I0,A3)") 'fig-bfgs-ini-',nsites,'.mp'
        call drawsol(n,x0,filename,c_loc(pdata))

        write(filename,"(A13,I0,A3)") 'fig-bfgs-end-',nsites,'.mp'
        call drawsol(n,x,filename,c_loc(pdata))

        write(filename,"(A12,I0,A4)")'tablinebfgs-',nsites,'.txt'
        open(unit=10, file=filename)
        write(10,1000) nsites,f,dsave(13),isave(30),isave(34),finish-start
        close(10)
    else 
        write(filename,"(A15,I0,A3)") 'fig-bfgs-ini-s-',nsites,'.mp'
        call drawsol(n,x0,filename,c_loc(pdata))

        write(filename,"(A15,I0,A3)") 'fig-bfgs-end-s-',nsites,'.mp'
        call drawsol(n,x,filename,c_loc(pdata))

        write(filename,"(A19,I0,A4)")'tablinebfgs-scaled-',nsites,'.txt'
        open(unit=10,file=filename)
        write(10,1001) scale,pdata%volA,nsites,f,dsave(13),isave(30),isave(34),finish-start
        close(10)
    end if

    activate = 0
    if ( scaled ) then
        if ( pdata%J2ctol .gt. 0.0d0 .or. activate .eq. 1 ) then
            if ( n .eq. 20 ) then
                call evalsizeedJi2(n,x,avg,sizes,ji2,flag,c_loc(pdata))

                write(*,*) 'avg = ', avg(1:n/2)
                write(*,*) 'Ji(a) = ', ji2(1:n/2)
                write(*,*)
                write(*,*) 'sizes = '
                do i=1,nsites
                    write(*,*) sizes(i,:)
                end do

                write(*,*)

                write(*,*) '|E|/E_i=', sizes(3,1)/avg(3)

                write(filename,"(A12,I0,A3,F4.2,A4)") 'info-bfgsj2-',nsites,'-c-',pdata%J2ctol,'.txt'
                open(unit=20,file=filename)
                write(20,"(A)") 'avg= '
                write(20, "(1P10E12.5)") (avg(i), i=1,nsites)
                write(20,"(A)") 'Ji(a)= '
                write(20, "(1P10E12.5)") (ji2(i), i=1,nsites)
                write(20,"(A)") 'sizes= '
                do i=1,nsites
                    write(20,"(1P10E12.5)") (sizes(i,j),j=1,10)
                end do
                do i=1,nsites
                    counter = count( sizes(i,:) .gt. 0.0d0 )
                    write(20,"(1P10E12.5)") minval(sizes(i,1:counter))/avg(i)
                end do
                write(20,"(A)") '|E|/E_i= '
                write(20,"(1PE12.5)") sizes(3,1)/avg(3)

                close(20)
            end if
            
            write(filename,"(A15,I0,A3)") 'fig-bfgsJ2-end-',nsites,'.mp'
            call drawsolJ2(n,x,filename,c_loc(pdata))
            write(filename,"(A15,I0,A4)") 'tabline-bfgsJ2-',nsites,'.txt'
            open(unit=10,file=filename)
            write(10,1002) scale,pdata%volA,pdata%npols,nsites,f,dsave(13),G1,j2,isave(30),isave(34),finish-start

        else if ( pdata%J3ctol .gt. 0.0d0 .or. activate .eq. 2 ) then
            if ( n .eq. 30) then
                call evalangJi3(n,x,avgang,angles,ji3,flag,c_loc(pdata))

                write(*,*) 'avg = ', avgang(1:n/2)
                write(*,*) 'Ji(a) = ', ji3(1:n/2)
                write(*,*)
                write(*,*) 'angles = '
                do i=1,nsites
                    write(*,*) angles(i,:)
                end do
                ! write(*,*)
                ! do i=1,nsites
                !     do j=1,10
                !         write(*,*) angles(i,j)/avgang(i)
                !     end do
                !     write(*,*)
                ! end do

                write(*,*)

                write(*,*) '|theta|/theta_i=', angles(3,1)/avgang(3)

                write(filename,"(A12,I0,A3,F4.2,A4)") 'info-bfgsj3-',nsites,'-c-',pdata%J3ctol,'.txt'
                open(unit=20,file=filename)
                write(20,"(A)") 'avg= '
                write(20, "(1P10E12.5)") (avgang(i), i=1,nsites)
                write(20,"(A)") 'Ji(a)= '
                write(20, "(1P10E12.5)") (ji3(i), i=1,nsites)
                write(20,"(A)") 'angles= '
                do i=1,nsites
                    write(20,"(1P10E12.5)") (angles(i,j),j=1,10)
                end do
                write(20,"(A)")
                do i=1,nsites
                    counter = count( angles(i,:) .gt. 0.0d0 )
                    write(20,"(1P10E12.5)") minval(angles(i,1:counter))/avgang(i)
                end do
                write(20,"(A)")
                write(20,"(A)") '|theta|/theta_i= '
                write(20,"(1PE12.5)") angles(3,1)/avgang(3)

                close(20)
            end if

            write(filename,"(A15,I0,A3)") 'fig-bfgsJ3-end-',nsites,'.mp'
            call drawsolJ3(n,x,filename,c_loc(pdata))
            write(filename,"(A15,I0,A4)") 'tabline-bfgsJ3-',nsites,'.txt'
            open(unit=10,file=filename)
            write(10,1002) scale,pdata%volA,pdata%npols,nsites,f,dsave(13),G1,j3,isave(30),isave(34),finish-start
        else
            if ( n .eq. 20) then 
                call evalvolJi1(n,x,vol,ji1,flag,c_loc(pdata))

                write(*,*) 'vol = ', vol(1:n/2)
                write(*,*) 'Ji(a) = ', ji1(1:n/2)
                write(*,*)

                write(*,*)

                write(filename,"(A12,I0,A4)") 'info-bfgsj1-',nsites,'.txt'
                open(unit=20,file=filename)
                write(20,"(A)") 'vol= '
                write(20, "(1P10E12.5)") (vol(i), i=1,nsites)
                write(20,"(A)") 'Ji(a)= '
                write(20, "(1P10E12.5)") (ji1(i), i=1,nsites)

                close(20)
            end if

            write(filename,"(A15,I0,A3)") 'fig-bfgsJ1-end-',nsites,'.mp'
            call drawsolJ1(n,x,filename,c_loc(pdata))
            write(filename,"(A15,I0,A4)") 'tabline-bfgsJ1-',nsites,'.txt'
            open(unit=10,file=filename)
            write(10,1002) scale,pdata%volA,pdata%npols,nsites,f,dsave(13),G1,j1,isave(30),isave(34),finish-start
        end if
    end if

        if ( printx ) then
            write(filename,"(A14,I0,A4)") 'xoptimal-bfgs-',nsites,'.txt'
            open(unit=20,file=filename)
            do i=1, n
                write(20,"(D24.12)") x(i)
            end do
            close(20)
        end if


    deallocate(x,x0,xold,nbd,l,u,wa,iwa,g,sites,ji,ji2,avg,sizes,ji3,avgang,angles,ji1,vol)
    if ( allocerr .ne. 0 ) then
        write(*,*) 'Deallocation error.'
        stop
    end if
    stop
1000 format(I6,1X,1P,E12.5,1X,1P,E7.1,2(1X,I6),1X,0P,F14.3)
1001 format(1P,E12.5,1X,1P,E12.5,1X,I6,1X,1P,E12.5,1X,1P,E7.1,2(1X,I6),1X,0P,F14.3)
1002 format(1P,E12.5,1X,1P,E12.5,1X,I3,1X,I6,1X,1P,E12.5,1X,1P,E7.1,2(1X,1P,E12.5),2(1X,I6),1X,0P,F14.3)

end program optvormabfgs