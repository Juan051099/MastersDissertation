program optvorma

    use iso_c_binding, only: c_ptr,c_loc,c_f_pointer
    use optvormod2, only: pdata_type,evalfg,proj,drawsol,InA,compscale,evalJg

    implicit none

    ! PARAMETERS
    integer, parameter :: nvert = 1000
    real(kind=8), parameter :: pi = acos( -1.0d0 )

    ! LOCAL SCALARS
    integer :: allocerr,fcnt,i,inform,iprint,iter,maxfc,maxit,maxnv,n,&
               nsites,spginfo
    real(kind=8) :: ang,area,epsopt,epsx,f,fbest,finish,ftarget,G1,gpsupn,j2,&
               j3,j1,j4,marg,scale,start
    type(pdata_type), target :: pdata
    
    ! LOCAL ARRAYS
    real(kind=8), allocatable :: g(:),gtmp(:),x(:),x0(:),ji(:)
    character(len=80) :: filename

    ! Set sites

    write(*,*) 'Enter nsites>0: '
    read(*,*) nsites

        ! Set disjoint convex polygons Aj whose union define A 
    
    pdata%npols = 1
    
    allocate(pdata%nvpols(pdata%npols),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if
  
    pdata%nvpols(1:pdata%npols) = (/ nvert + 1 /)
    
    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
    
    allocate(pdata%poledges(pdata%npols,maxnv),pdata%xpol(pdata%npols,maxnv),pdata%ypol(pdata%npols,maxnv),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if

    do i = 1,nvert + 1
        ang = 2.0d0 * pi * ( i - 1 ) / nvert
        pdata%xpol(1,i) = cos( ang )
        pdata%ypol(1,i) = sin( ang )
    end do

    pdata%poledges(1,1:nvert) = 0

    ! Scale A so cells have ideal area near unity
    
    call compscale(nsites,scale,area,c_loc(pdata))

    pdata%volA = area 
    pdata%xpol(1,1:pdata%nvpols(1)) = scale * pdata%xpol(1,1:pdata%nvpols(1))
    pdata%ypol(1,1:pdata%nvpols(1)) = scale * pdata%ypol(1,1:pdata%nvpols(1))
    
    pdata%xmin = minval( pdata%xpol(1,1:pdata%nvpols(1)) )
    pdata%xmax = maxval( pdata%xpol(1,1:pdata%nvpols(1)) )
    pdata%ymin = minval( pdata%ypol(1,1:pdata%nvpols(1)) )
    pdata%ymax = maxval( pdata%ypol(1,1:pdata%nvpols(1)) )
    
    marg = 0.1d0 * max( 1.0d0, abs( pdata%xmin ), abs( pdata%xmax ), abs( pdata%ymin), abs( pdata%ymax ) )

    pdata%xmin = pdata%xmin - marg
    pdata%xmax = pdata%xmax + marg
    pdata%ymin = pdata%ymin - marg
    pdata%ymax = pdata%ymax + marg
    

    ! pdata%c(1:2) = 0.5d0 * (/ pdata%xmax + pdata%xmin, pdata%ymax + pdata%ymin /)
    ! pdata%rmax = norm2(pdata%c - (/ pdata%xmin, pdata%ymin /))

    pdata%c(1:2) = 0.0d0
    pdata%rmax = scale
    write(*,*) 'radio= ', pdata%rmax
    pdata%drawscale = 2.5d0 / scale

    if ( nsites .ge. 10000 ) then
    pdata%drawscale = 5.0d0 / scale
    end if

    ! Goal function data
    pdata%Jweight(0) = 1.0d0
    pdata%Jweight(1:4) = 0.0d0
    pdata%Jweight(4) = 1.0d0

    pdata%J4c(1:2) = pdata%c(1:2)
    pdata%J4r      = scale
    write(*,*) 'scale=', scale
    n = 2 * nsites

    allocate(x(n),g(n),x0(n),gtmp(n),ji(n/2), stat=allocerr)
    if ( allocerr .ne. 0 ) then
        write(*,*)' Allocation error. '
        stop
    end if

    open(unit=10,file="/Users/juansebastian/Desktop/Mestrado USP/Dissertation/optlloyd/xoptimal-regpol-1000v1000.txt",status="old")

    do i=1,n
        read(10,"(D24.12)") x(i) 
    end do
    close(10)

    x0(1:n) = x(1:n)

    ! Optimize
    fbest = huge( 1.0d0 ) 
    epsopt = 1.0d-8 
    ftarget = 1.0d-12
    epsx = 0.0d0  
    maxit = 100000
    maxfc = 100000
    iprint = 0

    call cpu_time(start)

    call spg(n,x,epsopt,ftarget,maxit,maxfc,epsx,iprint,f,gpsupn,iter,fcnt,spginfo,inform, &
            evalfg,proj,c_loc(pdata))
    
    call cpu_time(finish)

    write(*,*) 'Statistics: f = ',f,' gpsupn = ',gpsupn,' it = ',iter,' fcnt = ',fcnt, &
    ' spginfo = ',spginfo,' time = ',finish-start

    call evalJg(n,x,G1,j2,j3,j1,j4,ji,g,inform,c_loc(pdata))

    write(*,*) 'G = ', G1, 'J4 = ',j4

    if ( f .lt. fbest ) then
        fbest = f

        write(*,*) 'Solution was improved. New best f = ',fbest

        ! write(filename,"(A15,I0,A3)") 'fig-regpol-ini-',nsites,'.mp'
        ! call drawsol(n,x0,filename,c_loc(pdata))
  
        write(filename,"(A17,I0,A3)") 'fig-regpolJ4-end-',nsites,'.mp'
        call drawsol(n,x,filename,c_loc(pdata))

        write(filename, "(A17,I0,A4)") 'tabline-regpolJ4-',nsites,'.txt'
        open(unit=10,file=filename)   
        write(10,1000) scale,pdata%volA,pdata%npols,nsites,f,gpsupn,G1,j4,iter,fcnt,spginfo,finish-start
        close(10)
        
        if ( fbest .le. ftarget ) then
           write(*,*) 'A global minimizer has been found!'
        end if
     end if

  
  deallocate(x,x0,g,gtmp,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if
  
  stop

1000 format(1P,E12.5,1X,1P,E12.5,1X,I3,1X,I6,1X,1P,E12.5,1X,1P,E7.1,2(1X,1P,E12.5),2(1X,I6),1X,I2,0P,F14.3)

end program optvorma