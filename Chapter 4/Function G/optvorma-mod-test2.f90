program optvorma

    use iso_c_binding, only: c_loc
    use optvormod, only: pdata_type,proj,evalf,evalg,approxg,drawsol

    implicit none

    ! LOCAL SCALARS
    integer :: allocerr,i,inform,maxnv,n,nsites,j,aux,maxfc,maxit, &
    spginfo,iter,iprint,fcnt
    real(kind=8) :: f,marg,seed,mindist,dist,epsopt,fbest,finish,ftarget, &
    gpsupn,start

    type(pdata_type), target :: pdata

    ! LOCAL ARRAYS 
    real(kind=8), allocatable :: g(:),x0(:),x(:),sites(:,:)
    character(len=80) :: filename

    ! PARAMETERS
    real(kind=8), parameter :: maxdisbs = 1.0d-01

    ! FUNCTIONS
    real(kind=8), external :: drand

    ! Set sites

    write(*,*) 'Enter nsites>0: '
    read(*,*) nsites

    !nsites = 3


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

    pdata%poledges(1,1:pdata%nvpols(1)-1) = (/ 0, 0, 0, 0 /) ! 0 means frontier of A, 1 means internarl

    pdata%volA = 1.0d0

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
    pdata%drawscale = 1.0d0

    n = 2 * nsites

    allocate(x(n),g(n),x0(n),sites(2,nsites),stat=allocerr)
    if ( allocerr .ne. 0) then
        write(*,*) 'Allocation error.'
        stop
    end if

    seed = 123456.0d0

    i=1
    do while( i .le. nsites )
        sites(1,i) = 0.10d0 + 0.80d0 * drand(seed)
        sites(2,i) = 0.10d0 + 0.80d0 * drand(seed)
        mindist = huge(1.0d0)
        do j=1,i-1
          dist = sqrt((sites(1,i) - sites(1,j))**2 + (sites(2,i) - sites(2,j))**2)
          mindist = min(dist,mindist)
        end do
        if( mindist .gt. maxdisbs ) i = i + 1
    end do

    aux = 0
    do i=1,nsites
        x(aux+1:aux+2) = (/ sites(1,i), sites(2,i) /)
        aux = aux + 2
    end do
    
    ! x(1:2) =  (/ 0.23d0, 0.65d0 /)
    ! x(3:4) =  (/ 0.87d0, 0.78d0 /)
    ! x(5:6) =  (/ 0.72d0, 0.32d0 /)

    x0(1:n) = x(1:n)
    !write(*,*) 'x0= ', x0(1:n)

    ! Optimize
    fbest = huge( 1.0d0 ) 
    epsopt = 1.0d-08
    ftarget = 1.0d-15
    maxit = 100000
    maxfc = 100000
    iprint = 0

    call cpu_time(start)

    call spg(n,x,epsopt,ftarget,maxit,maxfc,iprint,f,gpsupn,iter,fcnt,spginfo,inform, &
            evalf,evalg,proj,c_loc(pdata))
    
    call cpu_time(finish)

    !write(*,*) 'x =',x(1:n)

    write(*,*) 'Statistics: f = ',f,' gpsupn = ',gpsupn,' it = ',iter,' fcnt = ',fcnt, &
    ' spginfo = ',spginfo,' time = ',finish-start

    if ( f .lt. fbest ) then
        fbest = f 

        write(*,*) 'Solution was improved. New best f= ', fbest

        write(filename,"(A14,I0,A3)") 'fig-test-ini2-',nsites,'.mp'
        call drawsol(n,x0,filename,c_loc(pdata))

        write(filename,"(A14,I0,A3)") 'fig-test-end2-',nsites,'.mp'
        call drawsol(n,x,filename,c_loc(pdata))

        write(filename,"(A9,I0,A4)") 'tabline2-',nsites,'.txt'
        open(unit=10,file=filename)   
        write(10,1000) pdata%volA,pdata%npols,nsites,f,gpsupn,iter,fcnt,spginfo,finish-start
        close(10)

        if ( fbest .le. ftarget ) then
            write(*,*) 'A global minimizer has been found!'
        end if
    end if

    deallocate(x,x0,g,stat=allocerr)
    if ( allocerr .ne. 0 ) then
        write(*,*) 'Deallocation error.'
        stop
    end if
    stop

    1000 format(1P,E12.5,1X,I3,1X,I6,1X,1P,E12.5,1X,1P,E7.1,2(1X,I6),1X,I2,0P,F14.3)

end program optvorma