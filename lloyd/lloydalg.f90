program lloydalg

  use voro, only: voronoi,drawvor
  
  implicit none

  ! PARAMETER
  integer, parameter :: nvmax = 1000000                                    ! Max num of vertices.
  real(kind=8), parameter :: tol = 1.0d-06, maxdisbs = 1.0d-08
  logical, parameter :: perturb = .false., printx = .false.                ! By default perturb is false only use when you want change your initals random points.
                                                                           ! By default printx is false only use if you want to save the last iterate of the method. 

  ! LOCAL SCALARS
  integer :: allocerr,An,i,nsites,nv,istop,j,iter,n                        ! allocerr: To control the allocation,  An: size of the vectors which define the domain
                                                                           ! i:to generate random sites, nsites: number of sites, nv: number of vertices
  real(kind=8) :: seed,dist,mindist,start,finish                           ! seed: value to give to drand function

  ! LOCAL ARRAYS
  integer :: vflag(nvmax)                                                  ! vflag: flag for each vertex already visitted.
  integer, allocatable :: Aflag(:),sstart(:),colors(:)                     ! Aflag: flag for the points that generate the domain,   sstart?
  real(kind=8) :: vx(nvmax),vy(nvmax)                                      ! vx: coordinate x of the vertices, vy: coordinate y of the vertices
  real(kind=8), allocatable :: Ax(:),Ay(:),sites(:,:),cx(:),cy(:),Area(:),x(:) !Ax: x coordinates of the domain, Ay: y coordinates of the domain 
  
  character(len=80) :: filename                                            !sites: matrix with the coordinates of the sites
    
  ! FUNCTIONS 
  real(kind=8), external :: drand                                          ! Function that gives random points
  
  ! Set domain A 

  An = 5
  
  allocate(Ax(An+1),Ay(An+1),Aflag(An+1),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if

  Ax(1) = 0.0d0
  Ay(1) = 0.0d0
  ! Aflag(1) = 4
  Aflag(1) = 0

  Ax(2) = 1.0d0
  Ay(2) = 0.0d0
  ! Aflag(2) = 2
  Aflag(2) = 0

  Ax(3) = 1.0d0
  Ay(3) = 1.0d0
  ! Aflag(3) = 3
  Aflag(3) = 0

  Ax(4) = 0.0d0
  Ay(4) = 1.0d0
  ! Aflag(4) = 1
  Aflag(4) = 0

  Ax(5) = Ax(1)
  Ay(5) = Ay(1)
  Aflag(5) = Aflag(1)

  ! Set sites

  write(*,*) 'Enter nsites>0: '
  read(*,*) nsites

  ! nsites = 3
  
  allocate(sites(2,nsites),sstart(nsites+1),colors(nsites),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if
  
  ! open(unit=10,file="/Users/juansebastian/Desktop/Mestrado USP/Dissertation/optlloyd/xoptimal-500.txt",status="old")

  ! do j=1,nsites
  !   do i=1,2
  !     read(10,"(D24.12)") sites(i,j)
  !   end do 
  ! end do
  ! close(10)
  ! write(*,*) 'initsites='
  ! do j=1,2
  !    write(*,*) sites(j,:)
  ! end do

  ! sites(1:2,1) = (/ 0.23d0, 0.65d0 /)
  ! sites(1:2,2) = (/ 0.87d0, 0.78d0 /)
  ! sites(1:2,3) = (/ 0.72d0, 0.32d0 /)

  !!$ Sort sites within A
  
  seed = 123456.0d0

  do i = 1,nsites
     sites(1,i) = 0.10d0 + 0.80d0 * drand(seed)
     sites(2,i) = 0.10d0 + 0.80d0 * drand(seed)
  end do

  !! Uncomment this part if you want to put a constrain between the minimum distance of the sites.
  !! ---------------------------------------------------------------------------------------------
  i=1
  ! do while( i .le. nsites )
  !     sites(1,i) = 0.10d0 + 0.80d0 * drand(seed)
  !     sites(2,i) = 0.10d0 + 0.80d0 * drand(seed)
  !     mindist = huge(1.0d0)
  !     do j=1,i-1
  !       dist = sqrt((sites(1,i) - sites(1,j))**2 + (sites(2,i) - sites(2,j))**2)
  !       mindist = min(dist,mindist)
  !     end do
  !     if ( perturb ) then

  !       sites(1,i) = sites(1,i) + 2.0d0 * 1.0d-04 * drand(seed) - 1.0d-04 
  !       sites(2,i) = sites(2,i) + 2.0d0 * 1.0d-04 * drand(seed) - 1.0d-04 

  !     end if

  !     if( mindist .gt. maxdisbs ) i = i + 1
  ! end do
    !! ---------------------------------------------------------------------------------------------

  ! Sort sites within A
  
  n = 2 * nsites

  allocate(cx(nsites), cy(nsites), Area(nsites),x(2*nsites),stat=allocerr)
  if ( allocerr .ne. 0) then
    write(*,*) 'Allocation error.'
    stop
  end if

  call cpu_time(start)

  call lloyd(tol,nsites,perturb,sites,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop,iter,cx,cy,Area)

  call cpu_time(finish)

  do i=1,nsites
    x(2*i-1) = sites(1,i)
    x(2*i) = sites(2,i)
  end do

  write(*,*) 'x=',x(1:n)

  if ( perturb ) then
    write(filename,"(A9,I0,A4)") 'tablineP-',nsites,'.txt'
    open(unit=10,file=filename)
    write(10,1000) nsites,iter,finish-start
    close(10)
  else
    write(filename,"(A14,I0,A4)") 'tabline-lloyd-',nsites,'.txt'
    open(unit=10,file=filename)
    write(10,1000) nsites,iter,finish-start
    close(10)
    if ( printx ) then
      write(filename,"(A10,I0,A4)") 'xoptlloyd-',nsites,'.txt'
      open(unit=20,file=filename)
      do i=1, n
          write(20,"(D24.12)") x(i)
      end do
      close(20)
    end if

  end if

  deallocate(Ax,Ay,sites,cx,cy,Area,sstart,stat=allocerr)


  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if
  
  stop

  1000 format(1X,I6,1X,I6,1X,F14.3)

end program lloydalg


subroutine centroids(nv,nsites,vx,vy,sstart,cx,cy,Area)  ! Centroids with constant density
  
  implicit none
  
  ! SCALARS ARGUMENTS
  integer, intent(in) :: nv,nsites

  ! ARRAY ARGUMENTS
  integer, intent(in) :: sstart(nsites+1)
  real(kind=8), intent(in) :: vx(nv),vy(nv)
  real(kind=8), intent(out) :: Area(nsites),cx(nsites),cy(nsites)

  ! LOCAL ARGUMENTS
  integer :: i,j
  real(kind=8) :: a1,a2,cxaux,cyaux           ! Cx = (1/6A) * sum^{n-1}_{i=0}[(x_i + x_{i+1}) * !!$(x_i*y_{i+1}-x_{i+1}*y_i)] red part a1
                                              ! Cy = (1/6A) * sum^{n-1}_{i=0}[(y_i + y_{i+1}) * !!$(x_i*y_{i+1}-x_{i+1}*y_i)] " "
                                              ! A = 0.5 * sum^{n-1}_{i=0}[!!$(x_i*y_{i+1}-x_{i+1}*y_i)]
                                              ! a2 sum the values of a1 to get the total area of each cell                               
    do i=1, nsites
    a1 = 0.0d0
    a2 = 0.0d0
    cxaux = 0.0d0
    cyaux = 0.0d0
      do j=sstart(i), sstart(i+1)-1
        if ( j .ne. sstart(i+1)-1 ) then
          a1 = vx(j)*vy(j+1) - vx(j+1)*vy(j)
          cxaux = cxaux + (vx(j) + vx(j+1)) * a1
          cyaux = cyaux + (vy(j) + vy(j+1)) * a1
          a2 = a2 + a1
      else
          a1 = vx(j)*vy(sstart(i)) - vx(sstart(i))*vy(j)
          cxaux = cxaux + (vx(j) + vx(sstart(i))) * a1
          cyaux = cyaux + (vy(j) + vy(sstart(i))) * a1
          a2 = a2 + a1
      end if
    end do
    Area(i) = 0.5d0*a2
    cx(i) = cxaux/(6.0d0*Area(i))
    cy(i) = cyaux/(6.0d0*Area(i))
  end do

end subroutine centroids

subroutine lloyd(tol,nsites,perturb,sites,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop,iter,cx,cy,Area)
  
  use voro, only: voronoi,drawvor

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: nsites,An,nvmax
  integer, intent(out) :: istop,nv,iter
  real(kind=8), intent(in) :: tol
  logical, intent(in) :: perturb

  ! ARRAY ARGUMENTS
  integer, intent(inout) :: sstart(nsites+1),Aflag(An+1),vflag(nvmax)
  real(kind=8), intent(in) :: Ax(An+1),Ay(An+1)
  real(kind=8), intent(inout) :: sites(2,nsites),vx(nvmax),vy(nvmax),Area(nsites),cx(nsites),cy(nsites)

  ! LOCAL SCALARS 
  integer :: j
  real(kind=8) :: dist,maxdist

  ! PARAMETERS
  logical, parameter :: figeachiter = .false., debug = .false.


  ! LOCAL ARRAYS
  character(len=80) :: filename

  
  iter = 0
  if ( debug ) then
    write(*,*) 'initsites='
    do j=1,2
       write(*,*) sites(j,:)
    end do
  end if 

100 continue
  iter = iter + 1 
  call voronoi(nsites,sites,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop)
  call centroids(nv,nsites,vx,vy,sstart,cx,cy,Area)

  if ( iter .eq. 1 ) then

    if ( debug ) then
      write(*,*)
      write(*,*)'cx= ',cx
      write(*,*)'cy= ',cy
      write(*,*)'volW= ',Area
      write(*,*)
    end if

    if (.not. figeachiter) then

      if ( perturb ) then

        write(filename,"(A15,I0,A3)") 'fig-iniP-lloyd-',nsites,'.mp'
        call drawvor(nsites,sites,cx,cy,An,Ax,Ay,sstart,nv,vx,vy,filename)
      else

        write(filename,"(A14,I0,A3)") 'fig-ini-lloyd-',nsites,'.mp'
        call drawvor(nsites,sites,cx,cy,An,Ax,Ay,sstart,nv,vx,vy,filename)
      end if
    end if

  end if 

  j = 1
  ! maxdist = sqrt((cx(j)-sites(1,j))**2 + (cy(j)-sites(2,j))**2)
  maxdist = max( abs(cx(j)-sites(1,j)), abs(cy(j)-sites(2,j)) )
  do j=2,nsites
    ! dist = sqrt((cx(j)-sites(1,j))**2 + (cy(j)-sites(2,j))**2)
    dist = max( abs(cx(j)-sites(1,j)), abs(cy(j)-sites(2,j)) )
    if ( dist .gt. maxdist ) then
      maxdist = dist
    end if
  end do

  write(*,*) 'iter = ',iter,'maxdist = ',maxdist

  if ( figeachiter ) then
    write(filename,"(A4,I0,A1,I0,A3)") 'fig-',iter,'-',nsites,'.mp'
    call drawvor(nsites,sites,cx,cy,An,Ax,Ay,sstart,nv,vx,vy,filename) 
  end if
  
  sites(1, 1:nsites) = cx(1:nsites)
  sites(2, 1:nsites) = cy(1:nsites)

  if ( maxdist .gt. tol) then
    go to 100
  end if

  if (.not. figeachiter ) then

    if ( perturb ) then

      write(filename,"(A15,I0,A3)") 'fig-endP-lloyd-',nsites,'.mp'
      call drawvor(nsites,sites,cx,cy,An,Ax,Ay,sstart,nv,vx,vy,filename)

    else

      write(filename,"(A14,I0,A3)") 'fig-end-lloyd-',nsites,'.mp'
      call drawvor(nsites,sites,cx,cy,An,Ax,Ay,sstart,nv,vx,vy,filename)
    end if

  end if

  if ( debug ) then
    write(*,*) 'cx =', cx
    write(*,*) 'cy =', cy
    write(*,*) 'Area =', Area
    write(*,*)
    write(*,*) 'sstart = ',sstart(1:nsites+1)
    write(*,*) 'vx = ',vx(1:nv)
    write(*,*) 'vy = ',vy(1:nv)
    write(*,*) 'vflag = ',vflag(1:nv)
  end if
  
end subroutine lloyd


