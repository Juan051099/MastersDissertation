program optvorma
    !! USED CODE TO DO THE EXPERIMENTS CONSIDERING F_1, F_2, F_3
    use iso_c_binding, only: c_loc
    use optvormod2, only: pdata_type,proj,evalfg,approxg,drawsol,eval,drawsolJ2,drawsolJ3,evalsizeedJi2, &
    evalJg,inA,evalangJi3,evalvolJi1,drawsolJ1

    implicit none

    ! LOCAL SCALARS
    integer :: allocerr,i,j,inform,maxnv,n,nsites,maxfc,maxit, &
    spginfo,iter,iprint,fcnt,ntrials,itrial,flag                   
    real(kind=8) :: f,marg,seed,epsopt,fbest,finish,ftarget, &
    gpsupn,start,epsx,G1,j2,j3,j1,j4,scale,area

    type(pdata_type), target :: pdata

    ! LOCAL ARRAYS 
    real(kind=8) :: p(2)
    real(kind=8), allocatable :: g(:),x0(:),x(:),gtmp(:),ji(:),ji2(:),avg(:),sizes(:,:),avgang(:),ji3(:),&
                                 angles(:,:),vol(:),ji1(:)
    character(len=80) :: filename

    ! PARAMETERS

    logical, parameter :: printx = .false.                             ! use printx when need to generate a file with the final
                                                                       ! x obtained by minimizing the function   

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

    scale = sqrt( real(nsites) )
    area = scale ** 2
    write(*,*) 'area = ',area
    write(*,*) 'scale = ',scale
    pdata%volA = area

    do j = 1,pdata%npols
        pdata%xpol(j,1:pdata%nvpols(j)) = scale * pdata%xpol(j,1:pdata%nvpols(j))
        pdata%ypol(j,1:pdata%nvpols(j)) = scale * pdata%ypol(j,1:pdata%nvpols(j))
    end do

    pdata%xmin = minval( pdata%xpol(1,1:pdata%nvpols(1)) )
    pdata%xmax = maxval( pdata%xpol(1,1:pdata%nvpols(1)) )
    pdata%ymin = minval( pdata%ypol(1,1:pdata%nvpols(1)) )
    pdata%ymax = maxval( pdata%ypol(1,1:pdata%nvpols(1)) )

    marg = 0.1d0 * max( 1.0d0, abs( pdata%xmin ), abs( pdata%xmax ), abs( pdata%ymin ), abs( pdata%ymax ) )

    pdata%xmin = pdata%xmin - marg
    pdata%xmax = pdata%xmax + marg
    pdata%ymin = pdata%ymin - marg
    pdata%ymax = pdata%ymax + marg

    pdata%c(1:2) = 0.5d0 * (/ pdata%xmax + pdata%xmin, pdata%ymax + pdata%ymin /)
    pdata%rmax = norm2( pdata%c(1:2) - (/ pdata%xmin, pdata%ymin /) )
    pdata%drawscale = 1.0d0/scale

    ! Goal function data

    pdata%Jweight(0) = 1.0d0       !! Weight for the function G(a).
    pdata%Jweight(1:4) = 0.0d0
    ! pdata%Jweight(1) = 1.0d0
    ! pdata%Jweight(2) = 1.0d0
    ! pdata%Jweight(3) = 1.0d0
    ! write(*,*) 'Enter ctol for J2 (between 0 and 1): '
    ! read(*,*) pdata%J2ctol
    ! write(*,*) 'Enter ctol for J3 (between 0 and 1): '
    ! read(*,*) pdata%J3ctol
    ! pdata%J2ctol = 0.3d0

    if ( pdata%Jweight(2) .gt. 0.0d0 ) then
        write(*,*) 'Enter ctol for J2 (between 0 and 1): '
        read(*,*) pdata%J2ctol
    end if

    if ( pdata%Jweight(3) .gt. 0.0d0 ) then
        write(*,*) 'Enter ctol for J3 (between 0 and 1): '
        read(*,*) pdata%J3ctol
    end if

    ! Optimize

    write(*,*) 'Enter ntrials>0: '
    read(*,*) ntrials

    n = 2 * nsites

    allocate(x(n),g(n),gtmp(n),x0(n),ji(n/2),ji2(n/2),avg(n/2),sizes(n/2,10),ji3(n/2),&
            avgang(n/2),angles(n/2,10),ji1(n/2),vol(n/2),stat=allocerr)
    if ( allocerr .ne. 0) then
        write(*,*) 'Allocation error.'
        stop
    end if

    fbest = huge( 1.0d0 )

    do itrial=1, ntrials

        seed = 123456.0d0 * itrial

        write(*,*) 'Processing trial itrail = ', itrial

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
    !! ----------------------------------------------------------------------------------------
        !! This part is to calculate the value of the gradient with finite difference, the gradient constructed by us and compare them.

        ! call evalfg(n,x,f,g,inform,c_loc(pdata))
  
        ! write(*,*) 'f = ',f
        ! write(*,*) 'g = ',g(1:n)
   
        ! call approxg(n,x,gtmp,1.0d-06,inform,c_loc(pdata))
        ! write(*,*) 'approxg = ',gtmp(1:n)
   
        ! do i = 1,n
        !    write(*,*) i,g(i),gtmp(i),abs(g(i)-gtmp(i))
        ! end do
        
        ! write(*,*) 'Maximum absolute error = ',maxval( abs( g(1:n) - gtmp(1:n) ) )
        ! stop
    !! ----------------------------------------------------------------------------------------    

        x0(1:n) = x(1:n)
        write(*,*) 'x0= ', x0(1:n)

        epsopt = 1.0d-8 
        ftarget = 1.0d-8
        epsx = 0.0d0  
        maxit = 1000000
        maxfc = 1000000
        iprint = 0

        call cpu_time(start)

        call spg(n,x,epsopt,ftarget,maxit,maxfc,epsx,iprint,f,gpsupn,iter,fcnt,spginfo,inform, &
                evalfg,proj,c_loc(pdata))
        
        call cpu_time(finish)

        !write(*,*) 'x =',x(1:n)

        write(*,*) 'Statistics: f = ',f,' gpsupn = ',gpsupn,' it = ',iter,' fcnt = ',fcnt, &
        ' spginfo = ',spginfo,' time = ',finish-start

        call evalJg(n,x,G1,j2,j3,j1,j4,ji,g,flag,c_loc(pdata))

        ! write(*,*) 'G = ', G1, ' J1 = ', j1
        ! write(*,*) 'G = ', G1, ' J2 = ', j2
        write(*,*) 'G = ', G1, ' J3 = ', j3
        ! write(*,*) 'J_i = ', ji(1:n/2)

    !! ----------------------------------------------------------------------------------------
        !! Uncomment this part if you are minimizing the function f=G+J^2

        ! call evalsizeedJi2(n,x,avg,sizes,ji2,flag,c_loc(pdata))

        ! write(*,*) 'avg = ', avg(1:n/2)
        ! write(*,*) 'Ji(a) = ', ji2(1:n/2)
        ! write(*,*)
        ! write(*,*) 'sizes = '
        ! do i=1,nsites
        !     write(*,*) sizes(i,:)
        ! end do

        ! write(*,*)

        ! write(*,*) '|E|/E_i=', sizes(3,1)/avg(3)

        ! write(filename,"(A14,I0,A3,F4.2,A4)") 'info-squarej2-',nsites,'-c-',pdata%J2ctol,'.txt'
        ! open(unit=20,file=filename)
        ! write(20,"(A)") 'avg= '
        ! write(20, "(1P10E12.5)") (avg(i), i=1,nsites)
        ! write(20,"(A)") 'Ji(a)= '
        ! write(20, "(1P10E12.5)") (ji2(i), i=1,nsites)
        ! write(20,"(A)") 'sizes= '
        ! do i=1,nsites
        !     write(20,"(1P10E12.5)") (sizes(i,j),j=1,10)
        ! end do
        ! do i=1,nsites
        !     counter = count( sizes(i,:) .gt. 0.0d0 )
        !     write(20,"(1P10E12.5)") minval(sizes(i,1:counter))/avg(i)
        ! end do
        ! write(20,"(A)") '|E|/E_i= '
        ! write(20,"(1PE12.5)") sizes(3,1)/avg(3)

        ! close(20)
    !! ----------------------------------------------------------------------------------------

    !! ----------------------------------------------------------------------------------------
    !! Uncomment this part if you are minimizing the function f=G+J^3

        ! call evalangJi3(n,x,avgang,angles,ji3,flag,c_loc(pdata))

        ! write(*,*) 'avg = ', avgang(1:n/2)
        ! write(*,*) 'Ji(a) = ', ji3(1:n/2)
        ! write(*,*)
        ! write(*,*) 'angles = '
        ! do i=1,nsites
        !     write(*,*) angles(i,:)
        ! end do
        ! ! write(*,*)
        ! ! do i=1,nsites
        ! !     do j=1,10
        ! !         write(*,*) angles(i,j)/avgang(i)
        ! !     end do
        ! !     write(*,*)
        ! ! end do

        ! write(*,*)

        ! write(*,*) '|theta|/theta_i=', angles(3,1)/avgang(3)

        ! write(filename,"(A14,I0,A3,F4.2,A4)") 'info-squarej3-',nsites,'-c-',pdata%J3ctol,'.txt'
        ! open(unit=20,file=filename)
        ! write(20,"(A)") 'avg= '
        ! write(20, "(1P10E12.5)") (avgang(i), i=1,nsites)
        ! write(20,"(A)") 'Ji(a)= '
        ! write(20, "(1P10E12.5)") (ji3(i), i=1,nsites)
        ! write(20,"(A)") 'angles= '
        ! do i=1,nsites
        !     write(20,"(1P10E12.5)") (angles(i,j),j=1,10)
        ! end do
        ! write(20,"(A)")
        ! do i=1,nsites
        !     counter = count( angles(i,:) .gt. 0.0d0 )
        !     write(20,"(1P10E12.5)") minval(angles(i,1:counter))/avgang(i)
        ! end do
        ! write(20,"(A)")
        ! write(20,"(A)") '|theta|/theta_i= '
        ! write(20,"(1PE12.5)") angles(3,1)/avgang(3)

        ! close(20)
!! ----------------------------------------------------------------------------------------

!! ----------------------------------------------------------------------------------------
    !! Uncomment this part if you are minimizing the function f=G+J^1

        ! call evalvolJi1(n,x,vol,ji1,flag,c_loc(pdata))

        ! write(*,*) 'vol = ', vol(1:n/2)
        ! write(*,*) 'Ji(a) = ', ji1(1:n/2)
        ! write(*,*)

        ! write(*,*)

        ! write(filename,"(A14,I0,A4)") 'info-squarej1-',nsites,'.txt'
        ! open(unit=20,file=filename)
        ! write(20,"(A)") 'vol= '
        ! write(20, "(1P10E12.5)") (vol(i), i=1,nsites)
        ! write(20,"(A)") 'Ji(a)= '
        ! write(20, "(1P10E12.5)") (ji1(i), i=1,nsites)

        ! close(20)
!! ----------------------------------------------------------------------------------------


        if ( f .lt. fbest ) then
            fbest = f 

            write(*,*) 'Solution was improved. New best f= ', fbest
            write(*,*)

            if ( pdata%Jweight(2) .gt. 0.0d0 .and. pdata%Jweight(3) .eq. 0.0d0  ) then
                ! write(filename,"(A10,I0,A3,F10.2,A3)") 'testJ2ini-',nsites,'-c-',pdata%J2ctol,'.mp'
                ! call drawsol4(n,x0,filename,c_loc(pdata))
        
                write(filename,"(A17,I0,A3,F4.2,A3)") 'fig-squareJ2-end-',nsites,'-c-',pdata%J2ctol,'.mp'
                call drawsolJ2(n,x,filename,c_loc(pdata))
        
                write(filename,"(A17,I0,A3,F4.2,A4)") 'tabline-squareJ2-',nsites,'-c-',pdata%J2ctol,'.txt'
                open(unit=10,file=filename)   
                write(10,1000) scale,pdata%volA,pdata%npols,nsites,f,gpsupn,G1,j2,iter,fcnt,spginfo,finish-start
                close(10)


                if ( printx ) then
                    write(filename,"(A11,I0,A4)") 'xoptimalJ2-',nsites,'.txt'
                    open(unit=20,file=filename)
                    do i=1, n
                        write(20,"(D24.12)") x(i)
                    end do
                    close(20)
                end if
                close(10)

                if ( fbest .le. ftarget ) then
                    write(*,*) 'A global minimizer has been found!'
                    exit
                end if
            
            else if ( pdata%Jweight(3) .gt. 0.0d0 .and. pdata%Jweight(2) .eq. 0.0d0) then
                ! write(filename,"(A10,I0,A3,F10.2,A3)") 'testJ3ini-',nsites,'-c-',pdata%J3ctol,'.mp'
                ! call drawsol7(n,x0,filename,c_loc(pdata))
        
                write(filename,"(A17,I0,A3,F4.2,A3)") 'fig-squareJ3-end-',nsites,'-c-',pdata%J3ctol,'.mp'
                call drawsolJ3(n,x,filename,c_loc(pdata))
        
                write(filename,"(A17,I0,A3,F4.2,A4)") 'tabline-squareJ3-',nsites,'-c-',pdata%J3ctol,'.txt'
                open(unit=10,file=filename)   
                write(10,1000) scale,pdata%volA,pdata%npols,nsites,f,gpsupn,G1,j3,iter,fcnt,spginfo,finish-start
                close(10)

                if ( fbest .le. ftarget ) then
                    write(*,*) 'A global minimizer has been found!'
                    exit
                end if
            else if ( pdata%Jweight(1) .gt. 0.0d0 ) then
                ! write(filename,"(A15,I0,A3)") 'fig-square-ini-',nsites,'.mp'
                ! call drawsol3(n,x0,filename,c_loc(pdata))

                write(filename,"(A17,I0,A3)") 'fig-squareJ1-end-',nsites,'.mp'
                call drawsolJ1(n,x,filename,c_loc(pdata))

                write(filename,"(A17,I0,A4)") 'tabline-squareJ1-',nsites,'.txt'
                open(unit=10,file=filename)   
                write(10,1000) scale,pdata%volA,pdata%npols,nsites,f,gpsupn,G1,j1,iter,fcnt,spginfo,finish-start

            else 
                ! write(filename,"(A15,I0,A3)") 'fig-square-ini-',nsites,'.mp'
                ! call drawsolJ3(n,x0,filename,c_loc(pdata))

                write(filename,"(A15,I0,A3)") 'fig-square-end-',nsites,'.mp'
                call drawsolJ2(n,x,filename,c_loc(pdata))

                ! write(filename,"(A17,I0,A3)") 'fig-squareJ3-end-',nsites,'.mp'
                ! call drawsolJ3(n,x,filename,c_loc(pdata))

                ! write(filename,"(A17,I0,A3)") 'fig-squareJ1-end-',nsites,'.mp'
                ! call drawsolJ1(n,x,filename,c_loc(pdata))

                write(filename,"(A15,I0,A4)") 'tabline-square-',nsites,'.txt'
                open(unit=10,file=filename)   
                ! write(10,1000) scale,pdata%volA,pdata%npols,nsites,f,gpsupn,G1,j2,iter,fcnt,spginfo,finish-start
                ! write(10,1000) scale,pdata%volA,pdata%npols,nsites,f,gpsupn,G1,j3,iter,fcnt,spginfo,finish-start
                 write(10,1000) scale,pdata%volA,pdata%npols,nsites,f,gpsupn,G1,j1,iter,fcnt,spginfo,finish-start

                if ( printx ) then
                    write(filename,"(A16,I0,A4)") 'xoptimal-square-',nsites,'.txt'
                    open(unit=20,file=filename)
                    do i=1, n
                        write(20,"(D24.12)") x(i)
                    end do
                    close(20)
                end if
                close(10)

                if ( fbest .le. ftarget ) then
                    write(*,*) 'A global minimizer has been found!'
                    exit
                end if
                
            end if
        end if
    end do

    deallocate(x,x0,g,ji,ji2,avg,sizes,ji3,avgang,angles,ji1,vol,stat=allocerr)
    if ( allocerr .ne. 0 ) then
        write(*,*) 'Deallocation error.'
        stop
    end if
    stop

    1000 format(1P,E12.5,1X,1P,E12.5,1X,I3,1X,I6,1X,1P,E12.5,1X,1P,E7.1,2(1X,1P,E12.5),2(1X,I6),1X,I2,0P,F14.3)

end program optvorma