! This is a version of spg.f90 that computes f and g
! simultaneously. This version is recommended for the case in which
! the total number of functional evaluations of SPG is close to the
! number of iterations, i.e. the case in which SPG performs a single
! functional evaluation per iterationin in most iterations.

! ******************************************************************
! ******************************************************************

subroutine spg(n,x,epsopt,ftarget,maxit,maxfc,epsx,iprint,f,gpsupn,iter,fcnt, &
     spginfo,inform,evalfg,proj,pdataptr)

  use iso_c_binding, only: c_ptr
	
  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: iprint,maxfc,maxit,n
  integer, intent(out) :: fcnt,inform,iter,spginfo
  real(kind=8), intent(in) :: epsopt,epsx,ftarget
  real(kind=8), intent(out) :: gpsupn,f
  type(c_ptr), optional, intent(in) :: pdataptr

  ! ARRAY ARGUMENTS
  real(kind=8), intent(inout) :: x(n)

  ! SUBROUTINE ARGUMENTS
  external :: evalfg,proj
      
  ! PARAMETERS
  integer, parameter :: m = 10
  real(kind=8), parameter :: lmin = 1.0d-30,lmax = 1.0d+30

  ! LOCAL SCALARS
  integer :: ibest,lsinfo,nprogr
  real(kind=8) :: fbest,fnew,lambda,sts,sty

  ! LOCAL ARRAYS
  real(kind=8) :: g(n),gnew(n),gp(n),s(n),y(n),d(n),xbest(n),xnew(n),lastfv(0:m-1)

  if ( iprint .gt. 0 ) then
     write(* ,fmt=1000)
     write(10,fmt=1000)
     write(* ,fmt=1010) n
     write(10,fmt=1010) n
  end if

  inform = 0

  iter = 0
  fcnt = 0
  nprogr = 0
  
  lastfv(0:m-1) = - huge( 1.0d0 )

  call sproj(n,x,inform,pdataptr)
  if ( inform .ne. 0 ) return

  call sevalfg(n,x,f,g,inform,pdataptr)
  if ( inform .ne. 0 ) return

  fcnt = fcnt + 1

  lastfv(0) = f

  gp(1:n) = x(1:n) - g(1:n)

  call sproj(n,gp,inform,pdataptr)
  if (inform .ne. 0) return

  gp(1:n) = gp(1:n) - x(1:n)
  
  gpsupn = maxval( abs( gp(1:n) ) )

  if ( gpsupn .ne. 0.0d0) then
     lambda =  min( lmax, max( lmin, 1.0d0 / gpsupn ) )
  else
     lambda = 0.0d0
  end if

  fbest = f
  xbest(1:n) = x(1:n)
  ibest = iter
  
100 continue

  if ( iprint .gt. 0 ) then
     if ( mod(iter,10) .eq. 0 ) then
        write(* ,fmt=1020)
        write(10,fmt=1020)
     end if
     write(* ,fmt=1030) iter,f,gpsupn
     write(10,fmt=1030) iter,f,gpsupn
  end if

  if ( gpsupn .le. epsopt ) then
     spginfo = 0

     if ( iprint .gt. 0 ) then
        write(*, 1100)
        write(10,1100)
     end if
     
     go to 200
  end if

  if ( f .le. ftarget ) then
     spginfo = 1

     if ( iprint .gt. 0 ) then
        write(*, 1110)
        write(10,1110)
     end if
     
     go to 200
  end if

  if ( iter .ge. maxit ) then
     spginfo = 2

     if ( iprint .gt. 0 ) then
        write(*, 1120)
        write(10,1120)
     end if
     
     go to 200
  end if

  if ( fcnt .ge. maxfc ) then
     spginfo = 3

     if ( iprint .gt. 0 ) then
        write(*, 1130)
        write(10,1130)
     end if
     
     go to 200
  end if

!!$  if ( iter - ibest .gt. 2 * m ) then
!!$     spginfo = 4
!!$
!!$     if ( iprint .gt. 0 ) then
!!$        write(*, 1140)
!!$        write(10,1140)
!!$     end if
!!$     
!!$     go to 200
!!$  end if

  iter = iter + 1
  d(1:n) = x(1:n) - lambda * g(1:n)

  call sproj(n,d,inform,pdataptr)
  if (inform .ne. 0) return

  d(1:n) = d(1:n) - x(1:n)
  call ls(n,x,f,g,d,m,lastfv,maxfc,fcnt,xnew,fnew,gnew,lsinfo,inform,pdataptr)
  if ( inform .ne. 0 ) return
  
  if ( lsinfo .eq. 3 ) then
     spginfo = 3

     if ( iprint .gt. 0 ) then
        write(*, 1130)
        write(10,1130)
     end if

     go to 200
  end if

!   if ( maxval( abs( x(1:n) - xnew(1:n) ) ) .le. epsx ) then
!      spginfo = 6

!      if ( iprint .gt. 0 ) then
!         write(*, 1160)
!         write(10,1160)
!      end if

!      go to 200
!   end if

  if ( fnew .lt. f - sqrt( epsilon( 1.0d0 ) ) * abs( f ) .or. &
       fnew .gt. f + sqrt( epsilon( 1.0d0 ) ) * abs( f ) ) then
     nprogr = 0
  else
     nprogr = nprogr + 1
     if ( nprogr .gt. 100 ) then !era 100, solo para ver que pasa
        spginfo = 5

        if ( iprint .gt. 0 ) then
           write(*, 1150)
           write(10,1150)
        end if

        go to 200
     end if
  end if
  
  f = fnew
  lastfv(mod(iter,m)) = f

  s(1:n) = xnew(1:n) - x(1:n)
  y(1:n) = gnew(1:n) - g(1:n)
  sts = sum( s(1:n) ** 2 )
  sty = dot_product( s(1:n), y(1:n) )

  x(1:n) = xnew(1:n)
  g(1:n) = gnew(1:n)

  gp(1:n) = x(1:n) - g(1:n)
  
  call sproj(n,gp,inform,pdataptr)
  if ( inform .ne. 0 ) return

  gp(1:n) = gp(1:n) - x(1:n)
  
  gpsupn = maxval( abs( gp(1:n) ) )

  if ( sty .le. 0.0d0 ) then
     lambda = lmax
  else
     lambda = max( lmin, min( sts / sty, lmax ) )
  end if

  if ( f .lt. fbest ) then
     fbest = f
     xbest(1:n) = x(1:n)
     ibest = iter
  end if

  go to 100

200 continue

  if ( iprint .gt. 0 ) then
     write(* ,fmt=2000) iter,fcnt,f,gpsupn
     write(10,fmt=2000) iter,fcnt,f,gpsupn
  end if

  f = fbest
  x(1:n) = xbest(1:n)
  
1000 format(/,1X,78('='), &
          /,1X,'This is the SPECTRAL PROJECTED GRADIENT (SPG) for ', &
          'for convex-constrained',/,1X,'optimization. If you ', &
          'use this code, please, cite:',/, &
          /,1X,'E. G. Birgin, J. M. Martinez and M. Raydan, ', &
          'Nonmonotone spectral projected',/,1X,'gradient ', &
          'methods on convex sets, SIAM Journal on ', &
          'Optimization 10, pp.',/,1X,'1196-1211, 2000, and',/, &
          /,1X,'E. G. Birgin, J. M. Martinez and M. Raydan, ', &
          'Algorithm 813: SPG - software',/,1X,'for ', &
          'convex-constrained optimization, ACM Transactions ', &
          'on Mathematical',/,1X,'Software 27, pp. 340-349, ', &
          '2001.',/,1X,78('='))

1010 format(/,1X,'Entry to SPG.', &
          /,1X,'Number of variables: ',I7)

1020 format(/,4X,'ITER',10X,'F',8X,'GPSUPN')

1030 format(  1X,I7,1X,1P,D16.8,1X,1P,D7.1)
  
1100 format(/,1X,'Flag of SPG: Solution was found (small gradient norm).')
1110 format(/,1X,'Flag of SPG: Solution was found (target functional value reached).')
1120 format(/,1X,'Flag of SPG: Maximum of iterations reached.')
1130 format(/,1X,'Flag of SPG: Maximum of functional evaluations reached.')
!!$1140 format(/,1X,'Flag of SPG: Lack of progress.')
1150 format(/,1X,'Flag of SPG: Lack of progress (type 2).')
1160 format(/,1X,'Flag of SPG: Consecutive iterates are too close.')
2000 format(/,1X,'Number of iterations               : ',9X,I7, &
            /,1X,'Number of functional evaluations   : ',9X,I7, &
            /,1X,'Objective function value           : ',1P,D16.8, &
            /,1X,'Sup-norm of the projected gradient : ',9X,1P,D7.1)

contains
  
  ! ******************************************************************
  ! ******************************************************************

  subroutine ls(n,x,f,g,d,m,lastfv,maxfc,fcnt,xnew,fnew,gnew,lsinfo,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: maxfc,m,n
    integer, intent(out) :: lsinfo
    integer, intent(inout) :: fcnt,inform
    real(kind=8), intent(in) :: f
    real(kind=8), intent(out) :: fnew
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8),intent(in) :: d(n),g(n),lastfv(0:m-1),x(n)
    real(kind=8),intent(out) :: gnew(n),xnew(n)

    ! PARAMETERS
    real(kind=8), parameter :: gamma = 1.0d-04

    ! LOCAL SCALARS
    real(kind=8) :: alpha,atemp,fmax,gtd

    fmax = maxval( lastfv(0:m-1) )

    gtd = dot_product( g(1:n), d(1:n) )

    alpha = 1.0d0

    xnew(1:n) = x(1:n) + alpha * d(1:n)
    call sevalfg(n,xnew,fnew,gnew,inform,pdataptr)
    if ( inform .ne. 0 ) return

    fcnt = fcnt + 1

100 continue

    if ( fnew .le. fmax + gamma * alpha * gtd ) then
       lsinfo = 0
       return
    end if
    
    if ( fcnt .ge. maxfc ) then
       lsinfo = 3
       return
    end if

    if ( alpha .le. 0.1d0 ) then
       alpha = alpha / 2.0d0

    else
       atemp = ( - gtd * alpha ** 2 ) / ( 2.0d0 * ( fnew - f - alpha * gtd ) )
       
       if ( atemp .lt. 0.1d0 .or. atemp .gt. 0.9d0 * alpha ) then
          atemp = alpha / 2.0d0
       end if
       
       alpha = atemp
    end if
    
    xnew(1:n) = x(1:n) + alpha * d(1:n)
    
    call sevalfg(n,xnew,fnew,gnew,inform,pdataptr)
    if ( inform .ne. 0 ) return
    
    fcnt = fcnt + 1
    
    go to 100
    
  end subroutine ls
  
  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalfg(n,x,f,g,inform,pdataptr)
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    real(kind=8), intent(out) :: f
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    integer :: flag

    call evalfg(n,x,f,g,flag,pdataptr)

    ! This is true if f if Inf, - Inf or NaN
    if ( .not. f .gt. - 1.0d+99 .or. .not. f .lt. 1.0d+99 ) then
       f = 1.0d+99
    end if

    if ( flag .ne. 0 ) then
       inform = - 90
       call reperr(inform)
       return
    end if
    
  end subroutine sevalfg

  ! ******************************************************************
  ! ******************************************************************

  subroutine sproj(n,x,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)

    ! LOCAL SCALARS
    integer :: flag

    call proj(n,x,flag,pdataptr)
    
    if ( flag .ne. 0 ) then
       inform = - 92
       call reperr(inform)
       return
    end if
    
  end subroutine sproj
  
  ! ******************************************************************
  ! ******************************************************************

  subroutine reperr(inform)
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: inform

    if ( inform .eq. -90 ) then
       write(* ,fmt=100) 'EVALFG'
       write(10,fmt=100) 'EVALFG'

    else if ( inform .eq. -92 ) then
       write(* ,fmt=100) 'PROJ '
       write(10,fmt=100) 'PROJ '
    end if
  
100 format(/,1X,'*** There was an error in the user supplied ', &
         'subroutine ',A10,' ***',/)
  
  end subroutine reperr

end subroutine spg
