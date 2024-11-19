module optvormod2

    use iso_c_binding, only: c_ptr,c_f_pointer
    use VorCells_Polygons, only: Polygon,VoronoiCell,voronoi_cell,voronoi_cell_old,VorCellInterConvPol, &
         VorCellInterConvPol_Collinear,voronoi_cell_destroy,polygon_destroy
  
    implicit none
    
    type :: pdata_type
       integer :: npols
       real(kind=8) :: c(2),rmax,volA,xmin,xmax,ymin,ymax,drawscale,J2ctol,J3ctol,J4c(2),J4r
       integer, allocatable :: nvpols(:),poledges(:,:)
       real(kind=8), allocatable :: xpol(:,:),ypol(:,:)
       real(kind=8) :: Jweight(0:4)
       ! 0: CVT
       ! 1: cells with identical volume given by vol(A)/m
       ! 2: each cell with equally sized edges
       ! 3: each cell with equally sized angles
       ! 4: cells with volume proportional to density function f evaluated at their sites

    end type pdata_type
    
    public :: evalfg,approxg,drawsol,proj
  
  contains
  
   ! ******************************************************************
   ! ******************************************************************

  logical function InA(p,pdataptr)
      
  implicit none

  ! SCALAR ARGUMENTS
  type(c_ptr), optional, intent(in) :: pdataptr
  
  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: p(2)

  ! LOCAL SCALARS
  integer :: allocerr,j,maxnv
  type(pdata_type), pointer :: pdata

  ! LOCAL ARRAYS
  real(kind=8), allocatable :: vertices(:,:)
  
  call c_f_pointer(pdataptr,pdata)
    
  maxnv = maxval( pdata%nvpols(1:pdata%npols) )
    
  allocate(vertices(2,maxnv),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'InA: Allocation error.'
     stop
  end if

  j = 1
  InA = .false.
  do while ( j .le. pdata%npols .and. .not. InA )
     vertices(1,1:pdata%nvpols(j)-1) = pdata%xpol(j,1:pdata%nvpols(j)-1)
     vertices(2,1:pdata%nvpols(j)-1) = pdata%ypol(j,1:pdata%nvpols(j)-1)
     call polygon_contains_point_2d_convex(pdata%nvpols(j)-1,vertices,p,inA)
     j = j + 1
  end do
      
  deallocate(vertices,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'InA: Deallocation error.'
     stop
  end if
  
end function InA

! ******************************************************************
! ******************************************************************

  subroutine comparea(scale,area,pdataptr)
      
   implicit none
 
   ! SCALAR ARGUMENTS
   real(kind=8), intent(in) :: scale
   real(kind=8), intent(out) :: area
   type(c_ptr), optional, intent(in) :: pdataptr
     
   ! LOCAL SCALARS
   type(pdata_type), pointer :: pdata
   
   ! LOCAL SCALARS
   integer :: allocerr,j,maxnv
   real(kind=8) :: areaPj
   
   ! LOCAL ARRAYS
   real(kind=8), allocatable :: poltmp(:,:)
   
   call c_f_pointer(pdataptr,pdata)

   maxnv = maxval( pdata%nvpols(1:pdata%npols) )
 
   allocate(poltmp(2,maxnv),stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'Allocation error.'
      stop
   end if
       
   area = 0.0d0
   
   do j = 1,pdata%npols
      poltmp(1,1:pdata%nvpols(j)) = scale * pdata%xpol(j,1:pdata%nvpols(j))
      poltmp(2,1:pdata%nvpols(j)) = scale * pdata%ypol(j,1:pdata%nvpols(j))
    
      call polygon_area_2d(pdata%nvpols(j)-1,poltmp(1:2,1:pdata%nvpols(j)-1),areaPj)
       
      area = area + areaPj
   end do
   
   deallocate(poltmp,stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'Deallocation error.'
      stop
   end if

  end subroutine comparea

 ! ******************************************************************
 ! ******************************************************************

   subroutine compscale(nsites,scale,area,pdataptr)
   
      implicit none
   
      ! SCALAR ARGUMENTS
      integer, intent(in) :: nsites
      real(kind=8), intent(out) :: area,scale
      type(c_ptr), optional, intent(in) :: pdataptr

      ! PARAMETERS
      real(kind=8), parameter :: eps = 0.01d0
      
      ! LOCAL SCALARS
      real(kind=8) :: a,b,phia,phib,tphi,tscale,tarea
      type(pdata_type), pointer :: pdata
      
      call c_f_pointer(pdataptr,pdata)

      tscale = 1.0d0

      call comparea(tscale,tarea,pdataptr)

      tphi = tarea / nsites - 1.0d0

      a = tscale
      phia = tphi
      
      write(*,*) 'scale factor (a) = ',tscale,' ideal area (a) = ',tarea / nsites
      
      if ( abs( tphi ) .le. eps ) then
         scale = tscale
         area = tarea
         return
      end if

      if ( tphi .lt. 0.0d0 ) then
         do while ( tphi .lt. 0.0d0 )
            tscale = 2.0d0 * tscale
            call comparea(tscale,tarea,pdataptr)
            tphi = tarea / nsites - 1.0d0
         end do
         
         b = tscale
         phib = tphi
         
      else
         do while ( tphi .gt. 0.0d0 )
            tscale = 0.5d0 * tscale
            call comparea(tscale,tarea,pdataptr)
            tphi = tarea / nsites - 1.0d0
         end do
         
         b = tscale
         phib = tphi
      end if

      write(*,*) 'scale factor (b) = ',tscale,' ideal area (b) = ',tarea / nsites
      
      if ( abs( tphi ) .le. eps ) then
         scale = tscale
         area = tarea
         return
      end if

   10  continue
      
      tscale = 0.5d0 * ( a + b )
      
      call comparea(tscale,tarea,pdataptr)

      tphi = tarea / nsites - 1.0d0

      write(*,*) 'scale factor (c) = ',tscale,' ideal area (c) = ',tarea / nsites
      
      if ( abs( tphi ) .le. eps ) then
         scale = tscale
         area = tarea
         return
      end if

      if ( sign(1.0d0,tphi) .eq. sign(1.0d0,phia) ) then
         a = tscale
         phia = tphi
      else
         b = tscale
         phib = tphi
      end if

      go to 10
      
   end subroutine compscale
     
 ! ******************************************************************
 ! ******************************************************************

   subroutine proj(n,x,flag,pdataptr)
      
      implicit none
 
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      type(c_ptr), optional, intent(in) :: pdataptr
     
      ! ARRAY ARGUMENTS
      real(kind=8), intent(inout) :: x(n)

      ! LOCAL SCALARS
      integer :: i,nsites
      real(kind=8) :: tmp
      type(pdata_type), pointer :: pdata

      call c_f_pointer(pdataptr,pdata)

      flag = 0

      nsites = n / 2

      ! write(*,*)'xprojantes=',x(1:n)
      do i = 1,nsites
         ! x(2*i-1) = max(pdata%xmin,min(x(2*i-1),pdata%xmax))
         ! x(2*i) = max(pdata%ymin,min(x(2*i),pdata%ymax))
         tmp = norm2( x(2*i-1:2*i) - pdata%c(1:2) )
         if ( tmp .gt. pdata%rmax ) then
            x(2*i-1:2*i) = pdata%c(1:2) + pdata%rmax * ( x(2*i-1:2*i) - pdata%c(1:2) ) / tmp
         end if
      end do
      ! write(*,*)'xprojdepois=',x(1:n)
   
   end subroutine proj

 ! ******************************************************************
 ! ******************************************************************

   ! subroutine proj(n,x,flag,pdataptr)
      
   !    implicit none
  
   !    ! SCALAR ARGUMENTS
   !    integer, intent(in) :: n
   !    integer, intent(out) :: flag
   !    type(c_ptr), optional, intent(in) :: pdataptr
     
   !    ! ARRAY ARGUMENTS
   !    real(kind=8), intent(inout) :: x(n)
   !    real(kind=8) :: edges(3000000)
  
   !    ! LOCAL SCALARS
   !    integer :: i,nsites
   !    real(kind=8) :: xproj,yproj
   !    type(pdata_type), pointer :: pdata
  
   !    call c_f_pointer(pdataptr,pdata)
      
   !    write(*,*) 'Enter in the projection subroutine!.'
      
   !    do i=1, pdata%nvpols(1)-1
   !       if ( pdata%xpol(1,i) .eq. pdata%xpol(1,i+1) .and. pdata%ypol(1,i) .eq. pdata%ypol(1,i+1) ) then
   !          write(*,*) 'Error in function constraint: x1=x2 and y1=y2'
   !       end if
   !       if ( pdata%ypol(1,i) .ne. pdata%ypol(1,i+1) ) then
   !          edges(3*i-2) = 1.0d0
   !          edges(3*i-1) = -( pdata%xpol(1,i+1) - pdata%xpol(1,i) ) / ( pdata%ypol(1,i+1) -  pdata%ypol(1,i) ) 
   !          edges(3*i) = -( pdata%xpol(1,i) + edges(3*i-1)*pdata%ypol(1,i) )
   !       else
   !          edges(3*i-2) = 0.0d0
   !          edges(3*i-1) = 1.0d0
   !          edges(3*i) = -pdata%ypol(1,i)
   !       end if

   !       if ( i .ne. pdata%nvpols(1) - 1 ) then

   !          if ( edges(3*i-2) * pdata%xpol(1,i+2) + edges(3*i-1) * pdata%ypol(1,i+2) + edges(3*i) .gt. 0.0d0 ) then
   !             edges(3*i-2) = - edges(3*i-2)
   !             edges(3*i-1) = - edges(3*i-1)
   !             edges(3*i) = - edges(3*i)
   !          end if
   !       else
   !          if ( edges(3*i-2) * pdata%xpol(1,2) + edges(3*i-1) * pdata%ypol(1,2) + edges(3*i) .gt. 0.0d0 ) then
   !             edges(3*i-2) = - edges(3*i-2)
   !             edges(3*i-1) = - edges(3*i-1)
   !             edges(3*i) = - edges(3*i)
   !          end if
   !       end if

   !    end do
   !    ! write(*,*) 'vx=',pdata%xpol(1,1:pdata%nvpols(1))
   !    ! write(*,*) 'vy=',pdata%ypol(1,1:pdata%nvpols(1))
   !    ! write(*,*) 'edges = ', edges(1:12)
  
   !    flag = 0
  
   !    nsites = n / 2
  
   !    do i = 1,nsites
  
   !        ! Project x(2*i-1,2*i) into the domain (Polygon)
   !        call projauxP(edges,x(2*i-1),x(2*i),xproj,yproj,flag,pdataptr)
   !        x(2*i-1) = xproj
   !        x(2*i) = yproj
          
   !    end do

   !    ! write(*,*) 'xproj = ', x(1:n)

   ! end subroutine proj

   ! ******************************************************************
   ! ******************************************************************
  
   subroutine projauxP(edges,x,y,xproj,yproj,flag,pdataptr)
  
      ! This subroutine complements proj subroutine projecting
      ! a point z_i in R^2 onto its corresponding polygon P_i.
  
      ! SCALAR ARGUMENTS
      ! integer, intent(in) :: p,base
      real(kind=8), intent(in) :: x,y
      type(c_ptr), optional, intent(in) :: pdataptr
      integer, intent(out) :: flag
      real(kind=8), intent(out) :: xproj, yproj
  
      ! PARAMETERS
      integer, parameter :: nvsmax = 3000000
      
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: edges(3000000)
  
      !LOCAL SCALARS
      real(kind=8) :: a,b,c,dist,mindis,nearvx,nearvy,px,py,v1x,v1y,&
                   v2x,v2y,vx,vy
      integer :: i 
      logical :: fact
      type(pdata_type), pointer :: pdata
  
      !LOCAL ARRAYS   
      logical :: satisf(nvsmax)

      call c_f_pointer(pdataptr,pdata)
  
      ! TEST ALL EDGES OF THE POLYGON FOR SATISFIABILITY
      
      fact = .true.
  
      do i=1, pdata%nvpols(1)-1
          a = edges(3*i-2)
          b = edges(3*i-1)
          c = edges(3*i)
          if (a*x+b*y+c.le.0.0d0) then
              satisf(i) = .true.
  
          else
              fact = .false.
              satisf(i) = .false.
          end if
      end do
  
      if ( fact ) then
          xproj = x
          yproj = y
          return
      end if
  
      ! SEARCH FOR THE CLOSEST VERTEX
  
      mindis = huge(1.0d0)
      do i=1, pdata%nvpols(1)-1
          vx = pdata%xpol(1,i)
          vy = pdata%ypol(1,i)
          dist = sqrt( (vx-x)**2 + (vy-y)**2 )
          if ( dist .le. mindis ) then
              mindis = dist
              nearvx = vx
              nearvy = vy
          end if
      end do
  
      xproj = nearvx
      yproj = nearvy
  
      ! PROJECT ONTO THE VIOLATED CONSTRAINTS
  
      do i = 1, pdata%nvpols(1)-1
          if ( .not. satisf(i) ) then
              a = edges(3*i-2)
              b = edges(3*i-1)
              c = edges(3*i)
              if ( a .ne. 0.0d0 ) then
                  py = (y-b* (c/a+x)/a)/ (b**2/a**2+1)
                  px = - (c+b*py)/a
              else 
                  if ( b .ne. 0.0d0 ) then
                  px = (x-a * (c/b+y)/b)/ (a ** 2/b ** 2+1)
                  py = - (c+a * px)/b
                  else
                      write(*,*) 'ERROR IN PROBLEM DEFINITION (a=b=0 for a constraint of the type a x + b y + c = 0)'
                      flag = 1
                  end if
              end if
  
              v1x = pdata%xpol(1,i)
              v1y = pdata%ypol(1,i)
  
              v2x = pdata%xpol(1,i+1)
              v2y = pdata%ypol(1,i+1)
  
              ! if ( i .ne. pdata%nvpols(1,1) ) then
              !     v2x = pdata%xpol(1,i+1)
              !     v2y = pdata%ypol(1,i+1)
              ! else
              !     v2x = pdata%xpol(1,1)
              !     v2y = pdata%ypol(1,1)
              ! end if
  
              if ( min(v1x,v2x) .le. px .and. px .le. max(v1x,v2x) .and. & 
                  min(v1y,v2y).le. py .and. py .le. max(v1y,v2y) ) then
                  dist = sqrt((px-x)**2+ (py-y)**2)
                  if ( dist .le. mindis ) then
                      mindis = dist
                      xproj = px
                      yproj = py
                  end if
              end if
          end if
      end do
   end subroutine projauxP

   ! ******************************************************************
   ! ******************************************************************
  
    subroutine approxg(n,x,g,h,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      real(kind=8), intent(in) :: h
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(inout) :: x(n)
      real(kind=8), intent(out) :: g(n)
  
      ! LOCAL SCALARS
      integer :: i
      real(kind=8) :: fplus,fminus,tmp
  
      ! LOCAL ARRAYS
      real(kind=8) :: gtmp(n)
  
      do i = 1,n
         tmp = x(i)
         x(i) = tmp + h
         call evalfg(n,x,fplus,gtmp,flag,pdataptr)
         x(i) = tmp - h
         call evalfg(n,x,fminus,gtmp,flag,pdataptr)
         g(i) = ( fplus - fminus ) / ( 2.0d0 * h )
         x(i) = tmp
      end do
        
    end subroutine approxg
      
   ! ******************************************************************
   ! ******************************************************************
  
    subroutine evalfg(n,x,f,g,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      real(kind=8), intent(out) :: f
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: g(n)
  
      ! PARAMETERS
      logical, parameter :: compgrad = .true., debug = .false.
      integer, parameter :: gtmpnnzmax = 100000
      
      ! LOCAL SCALARS
      integer :: allocerr,gtmpnnz,i,ierror,k,maxnv,nsites,ntriag,status
      real(kind=8) :: ftmp
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: gtmpind(gtmpnnzmax),indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: fden(n/2),gden(2,n/2),gtmpval(gtmpnnzmax),sites(2,n/2),vorxy(2,2*(n/2)),valf(n/2)
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
  
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      ! Compute density function at sites if required
      
      if ( pdata%Jweight(4) .gt. 0.0d0 ) then
         do i = 1,nsites
            call fgcompJ4(sites(1:2,i),fden(i),gden(1:2,i),compgrad,pdataptr)
         end do
      end if

      f = 0.0d0
  
      if ( compgrad ) g(1:n) = 0.0d0
      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )

      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)
  
      do i = 1,nsites
        if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
        end if
  
        call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
        if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
        end if


        if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
        else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then 
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
        end if
        
        if ( pdata%Jweight(4) .gt. 0.0d0) then
            call compJ4(nsites,sites,i,W,pdata%volA,fden,gden,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0) return

            f = f + pdata%Jweight(4) * ftmp ** 2
            if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(4) * 2.0d0 * ftmp * gtmpval(1:gtmpnnz)
        end if
            
        if ( W%deg ) then
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if
        
            f = f + 1.0d0

        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
        
            call compG(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            ! call compGalternative(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            !call compGalternative1(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            !! Comment this lines to calculate the approximation of f numerically
            !!----------------
            f = f + pdata%Jweight(0)*ftmp   
            
            if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(0)*gtmpval(1:gtmpnnz)

            if ( pdata%Jweight(1) .gt. 0.0d0) then
               call compJ1(nsites,sites,i,W,pdata%volA,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
               if ( flag .ne. 0) return

               f = f + pdata%Jweight(1) * ftmp ** 2
               if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(1) * 2.0d0 * ftmp * gtmpval(1:gtmpnnz)
            end if
            
            if ( pdata%Jweight(2) .gt. 0.0d0 ) then
               call compJ2(nsites,sites,i,W,pdata%J2ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
               if ( flag .ne. 0 ) return

               f = f + pdata%Jweight(2) * ftmp 
               if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(2) * gtmpval(1:gtmpnnz)
            end if

            if ( pdata%Jweight(3) .gt. 0.0d0 ) then
               call compJ3(nsites,sites,i,W,pdata%J3ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
               if ( flag .ne. 0 ) return

               f = f + pdata%Jweight(3) * ftmp 
               if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(3) * gtmpval(1:gtmpnnz)
            end if

            !!----------------
         end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do

      !! Use this part to calculate the approximation of f numerically
      !!------------------------------------
      ! call calcfnum(nsites,sites,valf)
      ! f = f + sum( valf(1:nsites) )
      !!------------------------------------

      f = f / nsites
      
      if ( compgrad ) g(1:n) = g(1:n) / nsites
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
  
    end subroutine evalfg

   ! ******************************************************************
   ! ******************************************************************

   !! This function works to calculate the values of G, J1, J2, J3

    subroutine evalJg(n,x,G1,j2,j3,j1,j4,ji,g,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      real(kind=8), intent(out) :: j2,G1,j3,j1,j4
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: g(n),ji(n/2)
  
      ! PARAMETERS
      logical, parameter :: compgrad = .false., debug = .false.
      integer, parameter :: gtmpnnzmax = 100000
      
      ! LOCAL SCALARS
      integer :: allocerr,gtmpnnz,i,ierror,k,maxnv,nsites,ntriag,status
      real(kind=8) :: ftmp
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: gtmpind(gtmpnnzmax),indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: fden(n/2),gden(2,n/2),gtmpval(gtmpnnzmax),sites(2,n/2),vorxy(2,2*(n/2)),valf(n/2)
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
  
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      ! Compute density function at sites if required
   
      do i = 1,nsites
         call fgcompJ4(sites(1:2,i),fden(i),gden(1:2,i),compgrad,pdataptr)
      end do

      
      G1 = 0.0d0
      j1 = 0.0d0
      j2 = 0.0d0
      j3 = 0.0d0
      j4 = 0.0d0

      if ( compgrad ) g(1:n) = 0.0d0
      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )

      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)
  
      do i = 1,nsites
        if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
        end if
  
        call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
        if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
        end if


        if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
        else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then 
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
        end if
        

        call compJ4(nsites,sites,i,W,pdata%volA,fden,gden,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
        if ( flag .ne. 0 ) return
        ! write(*,*) 'ftmpevalJg = ', ftmp
        j4 = j4 + ftmp ** 2
        if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(4) * 2.0d0 * ftmp * gtmpval(1:gtmpnnz)
            
        if ( W%deg ) then

            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if

            ! write(*,*) 'Degenerate polygon, do not anything.'
            G1 = G1 + 1.0d0
            j1 = j1 + 1.0d0
            j2 = j2 + 1.0d0
            j3 = j3 + 1.0d0

        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
        
            call compG(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            G1 =  G1 + ftmp 
         
            if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + gtmpval(1:gtmpnnz)

            call compJ1(nsites,sites,i,W,pdata%volA,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            j1 = j1 + ftmp ** 2
            ji(i) =  ftmp
            if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(1) * 2.0d0 * ftmp * gtmpval(1:gtmpnnz)
            

            call compJ2(nsites,sites,i,W,pdata%J2ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            j2 = j2 + ftmp
            if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(1) * gtmpval(1:gtmpnnz)

            call compJ3(nsites,sites,i,W,pdata%J3ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            j3 = j3 + ftmp
            if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + pdata%Jweight(2) * gtmpval(1:gtmpnnz)

         end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do

      G1 = G1 / nsites
      j2 = j2 / nsites
      j3 = j3 / nsites
      j1 = j1 / nsites
      j4 = j4 / nsites
      if ( compgrad ) g(1:n) = g(1:n) / nsites
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
  
    end subroutine evalJg

   !! This function works to calculate the volume of each cell, and the J^1_i(a) recalling that
   !! the function J^1(a)=1/nsites sum^{nsites}_{i=1}[J^1_i(a)]^2. 

    subroutine evalvolJi1(n,x,vol,ji1,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: ji1(n/2), vol(n/2)
  
      ! PARAMETERS
      logical, parameter :: compgrad = .false.,  debug = .false.
      integer, parameter :: gtmpnnzmax = 100000
      
      ! LOCAL SCALARS
      integer :: allocerr,i,ierror,k,maxnv,nsites,ntriag,status,ed,ednext,gtmpnnz
      real(kind=8) :: perim,ftmp,volW
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: gtmpind(gtmpnnzmax),indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: sites(2,n/2),vorxy(2,2*(n/2)),vsizeed(10),taued(2),gtmpval(gtmpnnzmax)
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
  
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )

      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)

      
  
      do i = 1,nsites
        if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
        end if
  
        call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
        if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
        end if


        if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
        else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then 
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
        end if
        
            
        if ( W%deg ) then
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if
            write(*,*) 'Degenerate polygon, do not anything.'
        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
            call polygon_area_2d(W%n-1,W%v(1:2,1:W%n-1),volW)
            vol(i) = volW

            call compJ1(nsites,sites,i,W,pdata%volA,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            ji1(i) = ftmp
        end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
  
   end subroutine evalvolJi1

    ! *******************************************************************
    ! *******************************************************************

   !! This function works to calculate the average of each cell, the size of the edges of each cell and the J^2_i(a) recalling that
   !! the function J^2(a)=1/nsites sum^{nsites}_{i=1}J^2_i(a). 

   subroutine evalsizeedJi2(n,x,avg,sizes,ji2,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: ji2(n/2), avg(n/2), sizes(n/2,10)
  
      ! PARAMETERS
      logical, parameter :: compgrad = .false.,  debug = .false.
      integer, parameter :: gtmpnnzmax = 100000
      
      ! LOCAL SCALARS
      integer :: allocerr,i,ierror,k,maxnv,nsites,ntriag,status,ed,ednext,gtmpnnz
      real(kind=8) :: perim,ftmp
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: gtmpind(gtmpnnzmax),indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: sites(2,n/2),vorxy(2,2*(n/2)),vsizeed(10),taued(2),gtmpval(gtmpnnzmax)
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
  
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )

      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)

      
  
      do i = 1,nsites
        if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
        end if
  
        call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
        if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
        end if


        if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
        else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then 
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
        end if
        
            
        if ( W%deg ) then
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if

            write(*,*) 'Degenerate polygon, do not anything.'

        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
            vsizeed(1:10) = 0.0d0
            perim = 0.0d0
            do ed = 1,W%n - 1
               if ( ed .eq. W%n-1) then
                  ednext = 1
               else 
                  ednext = ed + 1
               end if 
               taued(1:2) = W%v(1:2,ednext) - W%v(1:2,ed)
               vsizeed(ed) = norm2( taued(1:2) )
               perim = perim + vsizeed(ed)
            end do
            sizes(i,1:10)=vsizeed(1:10)
            avg(i) = perim / (W%n-1)

            call compJ2(nsites,sites,i,W,pdata%J2ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            ji2(i) = ftmp
        end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
  
    end subroutine evalsizeedJi2

   ! *******************************************************************
   ! *******************************************************************

   !! This function works to calculate the average of the angles without fixed angles, the angles of each cell and the J^3_i(a) recalling that
   !! the function J^3(a)=1/nsites sum^{nsites}_{i=1}J^3_i(a). 


    subroutine evalangJi3(n,x,avgang,angles,ji3,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: ji3(n/2), avgang(n/2), angles(n/2,10)
  
      ! PARAMETERS
      logical, parameter :: compgrad = .false.,  debug = .false.
      integer, parameter :: gtmpnnzmax = 100000
      
      ! LOCAL SCALARS
      integer :: allocerr,i,ierror,k,maxnv,nsites,ntriag,status,edprev,ed,ednext,gtmpnnz,aux,counter
      real(kind=8) :: sumtau,ftmp
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: gtmpind(gtmpnnzmax),indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: sites(2,n/2),vorxy(2,2*(n/2)),vang(10),taued(2),tauedprev(2),gtmpval(gtmpnnzmax)
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
  
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )

      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)

      
  
      do i = 1,nsites
        if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
        end if
  
        call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
        if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
        end if


        if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
        else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then 
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
        end if
        
            
        if ( W%deg ) then
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if

            write(*,*) 'Degenerate polygon, do not anything.'

        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
            vang(1:10) = 0.0d0
            sumtau = 0.0d0
            aux = 0
            do ed = 1,W%n - 1
               if ( ed .eq. W%n-1 ) then
                  ednext = 1
                  edprev = ed - 1
               else if ( ed .eq. 1 ) then
                  ednext = ed + 1
                  edprev = W%n-1
               else
                  ednext = ed + 1
                  edprev = ed - 1
               end if
               if ( W%e(ed) .lt. 0 .or. W%e(edprev) .lt. 0 ) then 
                  taued(1:2) = W%v(1:2,ednext) - W%v(1:2,ed)
                  tauedprev(1:2) = W%v(1:2,ed) - W%v(1:2,edprev)
                  vang(ed-aux) = acos( max( -1.0d0, min( dot_product( -tauedprev(1:2)/norm2( tauedprev(1:2) ), taued(1:2)/norm2( taued(1:2) ) ), 1.0d0) ) )
                  sumtau = sumtau + vang(ed-aux)
               else
                  aux=1
               end if

            end do
            counter = count( vang(1:10) .ne. 0.0d0 )
            angles(i,1:10)=vang(1:10)
            avgang(i) = sumtau / counter

            call compJ3(nsites,sites,i,W,pdata%J3ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            ji3(i) = ftmp
        end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
  
    end subroutine evalangJi3

   ! *******************************************************************
   ! *******************************************************************


   !! This subroutine give us the minimun value of |E| and minimum value of \theta (angle between two edges) in the whole diagram
   !! Also gives a vector with the average size of the cell's of the edges of each Voronoi cell.

    subroutine eval(n,x,e,a,avgsizeed,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      real(kind=8), intent(out) :: e,a
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: avgsizeed(n/2)

      ! PARAMETERS
      logical, parameter ::  debug = .false.
      
      ! LOCAL SCALARS
      integer :: allocerr,i,ierror,k,maxnv,nsites,ntriag,status,ed,ednext,edprev
      real(kind=8) :: perim
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: fden(n/2),gden(2,n/2),sites(2,n/2),vorxy(2,2*(n/2)),valf(n/2),taued(2),tauedprev(2), &
                     vminang(n/2),vminsizeed(n/2),vsizeed(10),vang(10)
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
  
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )

      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)
  
      do i = 1,nsites
        if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
        end if
  
        call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
        if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
        end if


        if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
        else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then 
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
        end if
        
            
        if ( W%deg ) then
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if
            write(*,*) 'Degenerate polygon, do not anything.'

        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
            perim = 0.0d0
            do ed=1,W%n-1
               if ( ed .eq. W%n-1 ) then
                  edprev = ed - 1
                  ednext = 1
               else if ( ed .eq. 1) then
                  edprev = W%n-1
                  ednext = ed + 1
               else
                  edprev = ed - 1
                  ednext = ed + 1
               end if
               taued(1:2) = W%v(1:2,ednext) - W%v(1:2,ed)
               tauedprev(1:2)= W%v(1:2,ed) - W%v(1:2,edprev)
               vsizeed(ed) = norm2( taued(1:2) ) 
               vang(ed) = acos( max( -1.0d0, min( dot_product( -tauedprev(1:2)/norm2( tauedprev(1:2) ), taued(1:2)/vsizeed(ed) ), 1.0d0) ) )
               perim = perim + vsizeed(ed)
            end do
            vminang(i) = minval( vang(1:W%n-1) )
            vminsizeed(i) = minval( vsizeed(1:W%n-1) )
            avgsizeed(i) = perim/(W%n-1)

            if ( flag .ne. 0 ) return

         end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do
      e = minval(vminsizeed(1:nsites))
      a = minval(vminang(1:nsites))
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
  
    end subroutine eval

   ! ******************************************************************
   ! ******************************************************************
   
   !! This subroutine gives the option to calculate the gradient of the CVT function, one is the gradient obtain of the integral expression
   !! and the other is the expression obtained by the explicit function. 
   !! Put changeG = .true. if you want to calculate the gradient with the explicit form or .false. otherwise

   subroutine evalfg1(n,changeG,x,f,g,flag,pdataptr)


      implicit none
    
      ! SCALAR ARGUMENTS
      logical, intent(in) :: changeG
      integer, intent(in) :: n
      integer, intent(out) :: flag
      real(kind=8), intent(out) :: f
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: g(n)
  
      ! PARAMETERS
      logical, parameter :: compgrad = .true., debug = .false.
      integer, parameter :: gtmpnnzmax = 100000
      
      ! LOCAL SCALARS
      integer :: allocerr,gtmpnnz,i,ierror,k,maxnv,nsites,ntriag,status
      real(kind=8) :: ftmp
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: gtmpind(gtmpnnzmax),indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: fden(n/2),gden(2,n/2),gtmpval(gtmpnnzmax),sites(2,n/2),vorxy(2,2*(n/2))
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
  
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      
      f = 0.0d0
  
      if ( compgrad ) g(1:n) = 0.0d0
      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )

      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)
  
      do i = 1,nsites
        if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
        end if
  
        call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
        if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
        end if


        if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
        else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then 
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
        end if
        
            
        if ( W%deg ) then
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if
        
            f = f + 1.0d0

        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
        
            if ( changeG ) then
               call compGalternative(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
               if ( flag .ne. 0 ) return
            else
               call compG(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
               if ( flag .ne. 0 ) return
               
            end if
            f = f + ftmp
            if ( compgrad ) g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + gtmpval(1:gtmpnnz)
         end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do
  
      f = f / nsites
      
      if ( compgrad ) g(1:n) = g(1:n) / nsites
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
   


   
   end subroutine evalfg1

   ! ******************************************************************
   ! ******************************************************************


   subroutine calcfnum(nsites,sites,f)

      implicit none

      ! SCALAR ARGUMENTS
      integer, intent(in) :: nsites
      real(kind=8), intent(out) :: f(nsites)
      
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: sites(2,nsites)

      ! LOCAL SCALARS
      integer :: nx,ny,i,j,k

      ! LOCAL ARRAYS
      real(kind=8) :: u(2),diff(nsites)
      integer :: minindex(1)

      ! LOCAL PARAMETERS
      real(kind=8) :: h=1.0d-4

      nx = 1/h
      ny = 1/h
      f(1:nsites) = 0.0d0

      do i=1,nx
         u(1) = ( i-0.5d0 ) * h
         do j=1,ny
            u(2) =  ( j-0.5d0 ) * h 
            do k=1,nsites
               diff(k) = norm2( u(1:2)-sites(1:2,k) )
            end do
            minindex = minloc( diff(1:nsites) )
            f(minindex(1)) = f(minindex(1)) + h * h * diff(minindex(1))**2
            !write(*,*) f(minindex(1))
         end do
      end do
   end subroutine calcfnum

   ! ******************************************************************
   ! ******************************************************************

   subroutine compG(nsites,sites,i,W,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)

      implicit none

      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad,debug
      integer, intent(in) :: gnnzmax,i,nsites
      integer, intent(out) :: gnnz
      integer, intent(inout) :: flag
      real(kind=8), intent(out) :: f
      type(polygon), intent(in) :: W

      ! ARRAY ARGUMENTS
      integer, intent(out) :: gind(gnnzmax)
      real(kind=8), intent(in) :: sites(2,nsites)
      real(kind=8), intent(out) :: gval(gnnzmax)

      ! LOCAL SCALARS
      integer :: ed,ednext,nx,ny,j,k,r
      real(kind=8) :: vdet,volW,tmp,tinteg,tintegex

      ! LOCAL ARRAYS
      real(kind=8) :: diff1(2),diff2(2),cent(2),u(2),diff(nsites),v1(2),v2(2),v3(2),v4(3)
      integer :: minindex(1)
      ! LOCAL PARAMETERS
      real(kind=8), parameter :: h = 1.0d-4

      !!--------------------------------------------------------------------------------------------------------
      
      !! Calculating the function by quadrature naive version, in this case the domain is evaluated nsites times which is expensive.
      !! See calfcnum and the modification of evalfg.
      !! Calculate the integral of the first triangle numerically
      ! nx = 1/h
      ! ny = 1/h
      ! f = 0.0d0
      ! tinteg = 0.0d0 

      ! if ( i .eq. 1) then 
      !    call normalvec(W%v(1:2,2)-W%v(1:2,1), v1)
      !    call normalvec(sites(1:2,i) - W%v(1:2,2), v2)
      !    call normalvec(W%v(1:2,1) - sites(1:2,i), v3)
      !    write(*,*) v1
      !    write(*,*) v2
      !    write(*,*) v3
      ! end if


      ! do j=1,nx
      !    u(1) = ( j-0.5d0 ) * h
      !    do k=1,ny
      !       u(2) =  ( k-0.5d0 ) * h 
      !       do r=1,nsites
      !          diff(r) = norm2( u(1:2)-sites(1:2,r) )
      !       end do
      !       if ( i .eq. 1 ) then
      !          v4(1:3) = (/ dot_product( u(1:2)-W%v(1:2,1), -v1(1:2) ), dot_product( u(1:2)-W%v(1:2,2), -v2(1:2) ),&
      !                    dot_product( u(1:2)-sites(1:2,i), -v3(1:2) ) /)
      !          if ( maxval( v4(1:3) ) .le. 0 ) tinteg = tinteg + h * h * diff(i)**2
      !       end if

      !       minindex = minloc( diff(1:nsites) )
      !       if ( minindex(1) .eq. i  ) f = f + h * h * diff(i)**2
      !    end do
      ! end do

      ! if ( i .eq. 1 ) then
      !    diff1(1:2) = W%v(1:2,1) - sites(1:2,i)
      !    diff2(1:2) = W%v(1:2,2) - sites(1:2,i)
      !    vdet = ( diff1(1) * diff2(2) ) - ( diff1(2) * diff2(1) )  

      !    tintegex = abs( vdet ) * ( sum( diff1(1:2)**2 )+ &
      !             dot_product( diff1(1:2),diff2(1:2) ) + sum( diff2(1:2)**2 ) )
                  
      !    tintegex = tintegex/12.0d0
      !    write(*,*)'tintegnum = ', tinteg
      !    write(*,*)'tintegex = ',tintegex
      !    write(*,*)'Maximum absolute error tintegnum - tintegex', abs(tinteg - tintegex)
      !    write(*,*)'tintegex/tinteg = ', tintegex/tinteg
      ! end if
      !!--------------------------------------------------------------------------------------------------------

      if ( compgrad ) then        
         
         call centroid(W,cent,volW)

         if ( debug ) then
            write(*,*)
            write(*,*) 'volWi= ',volW
            write(*,*) 'centi= ',cent(1:2)
         end if

         gnnz = 0

         if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compG: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
         end if
         
         gind(gnnz+1:gnnz+2) = (/ 2*i-1, 2*i /)
         gval(gnnz+1:gnnz+2) = 0.0
         gnnz = gnnz + 2

         gval(1:2) = gval(1:2) + 2.0d0 * volW * ( sites(1:2,i) - cent(1:2)  )

      end if
      f = 0.0d0
      do ed = 1,W%n-1

         if ( ed .eq. W%n-1 ) then
            ednext = 1
         else 
            ednext = ed + 1
         end if

         if ( debug ) then
               write(*,*)
               write(*,*) '===================='
               write(*,*) 'Processing edge ed = ',ed
               write(*,*) 'Its inner vertex, named ',ed,' as well, corresponds to the point ',W%v(1:2,ed)
               write(*,*) 'Its outer vertex, named ',ednext,'as well, corresponds to the point ',W%v(1:2,ednext)
               write(*,*) '===================='
               write(*,*)
         end if



         if ( debug ) then
               write(*,*) 'Next: ',ednext
               write(*,*) 'It corresponds to the point',W%v(1:2,ednext)
         end if

         diff1(1:2) = W%v(1:2,ed) - sites(1:2,i)
         diff2(1:2) = W%v(1:2,ednext) - sites(1:2,i)
         vdet = ( diff1(1) * diff2(2) ) - ( diff1(2) * diff2(1) )  
         ! vdet = ( W%v(1,ed) - sites(1,i) ) * ( W%v(2,ednext) - sites(2,i) ) - &
         !        ( W%v(2,ed) - sites(2,i) ) * ( W%v(1,ednext) - sites(1,i) )
         !call det(W%v(1:2,ed), sites(1:2,i), W%v(1:2,ednext),vdet)

         if ( debug ) then
               write(*,*)
               write(*,*)'v_E - a_i = ',diff1(1:2)
               write(*,*)'w_E - a_i = ',diff2(1:2)
               write(*,*)'det = ',vdet
               write(*,*)'||v_E-a_i||^2= ', sum( diff1(1:2)**2 )
               write(*,*)'(v_E-a_i)^T(w_E-a_i)', dot_product( diff1(1:2),diff2(1:2) )
               write(*,*)'||w_E-a_i||^2= ',sum( diff2(1:2)**2 ) 
               write(*,*)
         end if 

         tmp = ( abs( vdet ) ) * ( sum( diff1(1:2)**2 ) + &
         dot_product( diff1(1:2),diff2(1:2) ) + sum( diff2(1:2)**2 ) )

         f = f  + tmp

         if ( debug ) then
            write(*,*)
            write(*,*)'tmp = ', tmp
            write(*,*)'f = ', f
            write(*,*)
         end if
      
      end do
   
      f = f /12.0d0
      if ( debug ) then
         write(*,*)'f= ', f
      end if

   end subroutine compG

   ! ******************************************************************
   ! ******************************************************************

   subroutine normalvec(v,normv)

      implicit none

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: v(2)
      real(kind=8), intent(out) :: normv(2)

      normv(1:2) = (/ -v(2), v(1) /)

   end subroutine normalvec

   ! ******************************************************************
   ! ******************************************************************

   subroutine compGalternative(nsites,sites,i,W,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)
    
      implicit none

      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad,debug
      integer, intent(in) :: gnnzmax,i,nsites
      integer, intent(out) :: gnnz
      integer, intent(inout) :: flag
      real(kind=8), intent(out) :: f
      type(polygon), intent(in) :: W

      ! ARRAY ARGUMENTS
      integer, intent(out) :: gind(gnnzmax)
      real(kind=8), intent(in) :: sites(2,nsites)
      real(kind=8), intent(out) :: gval(gnnzmax)

      ! LOCAL SCALARS
      integer :: ed, edprev, ednext
      real(kind=8) :: phia,phib,phic,xi,vdet,tmp,kap,absvdet,kappe,sig

      ! LOCAL ARRAYS
      real(kind=8) :: gradphi(2),M(2,2),wrot(2),vrot(2),p(2),diff1(2),diff2(2),ai(2),q(2)
      

      xi = 1.0d0/12.0d0
      f = 0.0d0

      if ( debug ) then
         write(*,*)
         write(*,*) 'xi = ', xi
         write(*,*)
      end if

      if ( compgrad ) then
         gnnz = 0
      
         if ( gnnz + 2 .gt. gnnzmax ) then
            write(*,*) 'In compGalternative: Increase gnnzmax to at least ',gnnz+2
            flag = - 1
            return
         end if
         
         gind(gnnz+1:gnnz+2) = (/ 2*i-1, 2*i /)
         gval(gnnz+1:gnnz+2) = 0.0d0
         gnnz = gnnz + 2
      end if
         
      do ed = 1,W%n-1

         if ( debug ) then
            write(*,*)
            write(*,*) '===================='
            write(*,*) 'Processing edge ed = ',ed
            write(*,*) 'Its inner vertex, named ',ed,' as well, corresponds to the point ',W%v(1:2,ed)
            write(*,*) '===================='
            write(*,*)
         end if

         if ( ed .eq. 1 ) then
            edprev = W%n - 1
            ednext = ed + 1
         else if ( ed .eq. W%n-1 ) then
            edprev = ed - 1
            ednext = 1
         else 
            edprev = ed - 1
            ednext = ed + 1
         end if

         ! Calculating the objective function

         diff1(1:2) = W%v(1:2,ed) - sites(1:2,i)
         diff2(1:2) = W%v(1:2,ednext) - sites(1:2,i)
         vdet = ( diff1(1) * diff2(2) ) - ( diff1(2) * diff2(1) )  

         if ( debug ) then
               write(*,*)
               write(*,*)'v_E - a_i = ',diff1(1:2)
               write(*,*)'w_E - a_i = ',diff2(1:2)
               write(*,*)'det = ',vdet
               write(*,*)'||v_E-a_i||^2= ', sum( diff1(1:2)**2 )
               write(*,*)'(v_E-a_i)^T(w_E-a_i)', dot_product( diff1(1:2),diff2(1:2) )
               write(*,*)'||w_E-a_i||^2= ',sum( diff2(1:2)**2 ) 
               write(*,*)
         end if 

         kap = sum( diff1(1:2)**2 ) + dot_product( diff1(1:2),diff2(1:2) ) + sum( diff2(1:2)**2 )
         absvdet = abs(vdet) 
         tmp = absvdet * kap

         f = f + tmp
         ! f = f  + absvdet

         if ( debug ) then
            write(*,*)
            write(*,*)'tmp = ', tmp
            write(*,*)'f = ', f
            write(*,*)
         end if
         
         if ( compgrad ) then

            ! Processing edge (v_E,w_E)

            if ( debug ) then
               write(*,*) 'Previous: ', edprev, 'next: ',ednext
               write(*,*) 'They corresponds to the points ',W%v(1:2,edprev),' and ', W%v(1:2,ednext),' respectively.'
            end if

            vrot(1:2) = (/ -W%v(2,ed), W%v(1,ed) /)
            wrot(1:2) = (/ W%v(2,ednext), -W%v(1,ednext) /)
            ai(1:2) = (/ -sites(2,i), sites(1,i)/)
            kappe = sign(1.0d0,vdet) * kap
            p(1:2) = kappe*(wrot(1:2)+ai(1:2))+absvdet*(2.0d0*W%v(1:2,ed) + W%v(1:2,ednext)-3.0d0*sites(1:2,i))
            q(1:2) = kappe*(vrot(1:2)-ai(1:2))+absvdet*(W%v(1:2,ed)+2.0d0*W%v(1:2,ednext)-3.0d0*sites(1:2,i))

            sig = sign(1.0d0,vdet)


            if ( debug ) then
               write(*,*)
               write(*,*) 'vrot= ',vrot
               write(*,*) 'wrot= ',wrot
               write(*,*) 'ai*= ', ai
               write(*,*) 'k(E)= ',kappe
               write(*,*)
            end if
            
            gval(1:2) = gval(1:2) - xi*(3.0d0*absvdet * (W%v(1:2,ed) + W%v(1:2,ednext)-2.0d0*sites(1:2,i)) + &
                                    kappe*(wrot(1:2)+vrot(1:2)))
            ! if ( i .eq. 1 ) then
            !    write(*,*)'first term just depends of \delta ai.'
            !    write(*,*)'summed(1:2)= ', - xi*(absvdet * (9.0d0*W%v(1:2,ed) + 5.0d0*W%v(1:2,ednext)-14.0d0*sites(1:2,i)) + &
            !    kappe*(wrot(1:2)+vrot(1:2)))
            !    write(*,*)'gi(1:2)= ',gval(1:2)
            !    write(*,*)
            ! end if

            ! gval(1:2) = gval(1:2) - sig * xi * ( wrot(1:2) + vrot(1:2))
         
            ! Starting point v_E

            if ( debug ) then
               write(*,*)
               write(*,*) 'Processing v_E'
               write(*,*)
            end if

            if ( W%e(ed) .lt. 0 .and. W%e(edprev) .lt. 0) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(v_E)|=3. ', &
                  'Neighbour cells are: ',abs(W%e(edprev)),' and ',abs(W%e(ed))
               end if

               call matrixint(sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return

               gval(1:2) = gval(1:2) + xi * matmul( transpose( M(1:2,1:2) ), p(1:2) )
               ! gval(1:2) = gval(1:2) + sig * xi * matmul( transpose( M(1:2,1:2) ), wrot(1:2) + ai(1:2) )
               ! if ( i .eq. 1 ) then
               !    write(*,*)'Case |P(v_E)|=3.'
               !    write(*,*)'summed(1:2)= ', xi * matmul( transpose( M(1:2,1:2) ), p(1:2) ) 
               !    write(*,*)'gi(1:2)= ',gval(1:2)
               !    write(*,*)
               ! end if

               call matrixint(sites(1:2,i),sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return
            
               if ( gnnz + 2 .gt. gnnzmax) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ',gnnz + 2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), p(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), wrot(1:2) + ai(1:2) )
               gnnz = gnnz + 2

               call matrixint(sites(1:2,abs(W%e(ed))),sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ',gnnz + 2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), p(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), wrot(1:2) + ai(1:2) )
               gnnz = gnnz + 2

            else if ( W%e(ed) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(v_E)|=2 (type 1). Neighbour cell is: ',abs(W%e(ed))
               end if

               call line_exp2imp_2d(W%v(1:2,ed),W%v(1:2,edprev),phia,phib,phic)

               gradphi(1:2) = (/ phia, phib /)

               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = -gradphi(1:2)
               end if

               call matrixbd(sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               gval(1:2) = gval(1:2) + xi * matmul( transpose( M(1:2,1:2) ), p(1:2) )
               !gval(1:2) = gval(1:2) + sig * xi * matmul( transpose( M(1:2,1:2) ), wrot(1:2) + ai(1:2) )
               ! if ( i .eq. 1 ) then
               !    write(*,*)'Case |P(v_E)|=2 (type 1).'
               !    write(*,*)'summed(1:2)= ', xi * matmul( transpose( M(1:2,1:2) ), p(1:2) ) 
               !    write(*,*)'gi(1:2)= ',gval(1:2)
               !    write(*,*)
               ! end if

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), p(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), wrot(1:2) + ai(1:2) )
               gnnz = gnnz + 2

            else if ( W%e(edprev) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(v_E)|=2 (type 2). Neighbour cell is: ', abs(W%e(edprev))
               end if

               call line_exp2imp_2d(W%v(1:2,ed),W%v(1:2,ednext),phia,phib,phic)

               gradphi(1:2) = (/ phia, phib /)

               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = -gradphi(1:2)
               end if

               call matrixbd(sites(1:2,abs(W%e(edprev))),sites(1:2,i),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               gval(1:2) = gval(1:2) + xi * matmul( transpose( M(1:2,1:2) ), p(1:2) )

               !gval(1:2) = gval(1:2) + sig * xi * matmul( transpose( M(1:2,1:2) ), wrot(1:2) + ai(1:2) )

               ! if ( i .eq. 1 ) then
               !    write(*,*)'Case |P(v_E)|=2 (type 2).'
               !    write(*,*)'summed(1:2)= ', xi * matmul( transpose( M(1:2,1:2) ), p(1:2) ) 
               !    write(*,*)'gi(1:2)= ',gval(1:2)
               !    write(*,*)
               ! end if

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), p(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), wrot(1:2) + ai(1:2) )
               gnnz = gnnz + 2
            end if

            ! Ending point w_E

            if ( debug ) then
               write(*,*)
               write(*,*) 'Processing w_E'
               write(*,*)
            end if

            if (W%e(ed) .lt. 0 .and. W%e(ednext) .lt. 0) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=3. ',&
                  'Neighbour cells are: ',abs(W%e(ed)),' and ',abs(W%e(ednext))
               end if

               call matrixint(sites(1:2,abs(W%e(ed))),sites(1:2,abs(W%e(ednext))),sites(1:2,i),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0 ) return

               gval(1:2) = gval(1:2) + xi * matmul( transpose( M(1:2,1:2) ), q(1:2) )
               !gval(1:2) = gval(1:2) + sig * xi * matmul( transpose( M(1:2,1:2) ), vrot(1:2) - ai(1:2) )

               ! if ( i .eq. 1 ) then
               !    write(*,*)'Case |P(w_E)|=3.'
               !    write(*,*)'summed(1:2)= ', xi * matmul( transpose( M(1:2,1:2) ), p(1:2) ) 
               !    write(*,*)'gi(1:2)= ',gval(1:2)
               !    write(*,*)
               ! end if

               call matrixint(sites(1:2,i),sites(1:2,abs(W%e(ed))),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), q(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), vrot(1:2) - ai(1:2) )
               gnnz = gnnz + 2

               call matrixint(sites(1:2,abs(W%e(ednext))),sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), q(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), vrot(1:2) - ai(1:2) )
               gnnz = gnnz + 2

            else if ( W%e(ed) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=2 (type 1). Neighbour cell is: ',abs(W%e(ed))
               end if

               call line_exp2imp_2d(W%v(1:2,ednext),W%v(1:2,ednext+1),phia,phib,phic)

               gradphi(1:2) = (/ phia, phib /)

               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = -gradphi(1:2)
               end if
               
               call matrixbd(sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               gval(1:2) = gval(1:2) + xi * matmul( transpose( M(1:2,1:2) ), q(1:2) )
               !gval(1:2) = gval(1:2) + sig * xi * matmul( transpose( M(1:2,1:2) ), vrot(1:2) - ai(1:2) )

               ! if ( i .eq. 1 ) then
               !    write(*,*)'Case |P(w_E)|=2 (type 1).'
               !    write(*,*)'summed(1:2)= ', xi * matmul( transpose( M(1:2,1:2) ), p(1:2) ) 
               !    write(*,*)'gi(1:2)= ',gval(1:2)
               !    write(*,*)
               ! end if

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), q(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), vrot(1:2) - ai(1:2) )
               gnnz = gnnz + 2

            else if ( W%e(ednext) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=2 (type 2). Neighbour cell is: ',abs(W%e(ednext))
               end if

               call line_exp2imp_2d(W%v(1:2,ednext),W%v(1:2,ed),phia,phib,phic)

               gradphi(1:2) = (/ phia, phib /)

               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = -gradphi(1:2)
               end if
               
               call matrixbd(sites(1:2,abs(W%e(ednext))),sites(1:2,i),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               gval(1:2) = gval(1:2) + xi * matmul( transpose( M(1:2,1:2) ), q(1:2) )
               !gval(1:2) = gval(1:2) + sig * xi * matmul( transpose( M(1:2,1:2) ), vrot(1:2) - ai(1:2) )

               ! if ( i .eq. 1 ) then
               !    write(*,*)'Case |P(w_E)|=2 (type 2).'
               !    write(*,*)'summed(1:2)= ', xi * matmul( transpose( M(1:2,1:2) ), p(1:2) ) 
               !    write(*,*)'gi(1:2)= ',gval(1:2)
               !    write(*,*)
               ! end if

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compGalternative: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
               gval(gnnz+1:gnnz+2) = xi * matmul( transpose( M(1:2,1:2) ), q(1:2) )
               !gval(gnnz+1:gnnz+2) = sig * xi * matmul( transpose( M(1:2,1:2) ), vrot(1:2) - ai(1:2) )
               gnnz = gnnz + 2
            end if
         end if

      end do

      f = xi * f

      if ( debug ) then
         write(*,*) 'f = ',f
      end if
   end subroutine compGalternative

   ! ******************************************************************
   ! ******************************************************************
   
   subroutine compGalternative1(nsites,sites,i,W,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)
   
      implicit none

      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad,debug
      integer, intent(in) :: gnnzmax,i,nsites
      integer, intent(out) :: gnnz
      integer, intent(inout) :: flag
      real(kind=8), intent(out) :: f
      type(polygon), intent(in) :: W

      ! ARRAY ARGUMENTS
      integer, intent(out) :: gind(gnnzmax)
      real(kind=8), intent(in) :: sites(2,nsites)
      real(kind=8), intent(out) :: gval(gnnzmax)

      ! LOCAL SCALARS
      integer :: ed,ednext
      real(kind=8) :: vdet,sizeed,mu,volW

      ! LOCAL ARRAYS
      real(kind=8) :: p(2),diff1(2),diff2(2),zeta(2),beta(2),auxc(6),auxv(2,2),glcnst(2),cent(2)


      ! call polygon_area_2d(W%n-1,W%v(1:2,1:W%n-1),volW)

      ! if ( debug ) then
      !       write(*,*) 'vol=A',volW
      ! end if

      f = 0.0d0
      
      ! if ( compgrad ) then  

      !    glcnst(1:2) = 0.0d0
      !    do ed=1,W%n-1
      !       if ( ed .eq. W%n-1 ) then
      !          ednext = 1
      !       else 
      !          ednext = ed + 1
      !       end if
      !       call det(W%v(1:2,ed),sites(1:2,i),W%v(1:2,ednext),vdet)
      !       glcnst = glcnst + abs( vdet ) * ( 2.0d0*W%v(1:2,ed) + W%v(1:2,ednext) )
      !    end do

      !    gnnz = 0

      !    if ( gnnz + 2 .gt. gnnzmax ) then
      !       write(*,*) 'In compG: Increase gnnzmax to at least ',gnnz+2
      !       flag = - 1
      !       return
      !    end if
         
      !    gind(gnnz+1:gnnz+2) = (/ 2*i-1, 2*i /)
      !    gval(gnnz+1:gnnz+2) = 2.0d0 * ( sites(1:2,i)*volW - glcnst/6.0d0 )
      !    gnnz = gnnz + 2

      ! end if

      if ( compgrad ) then        
         
         call centroid(W,cent,volW)

         if ( debug ) then
            write(*,*)
            write(*,*) 'volWi= ',volW
            write(*,*) 'centi= ',cent(1:2)
         end if

         gnnz = 0

         if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compG: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
         end if
         
         gind(gnnz+1:gnnz+2) = (/ 2*i-1, 2*i /)
         gval(gnnz+1:gnnz+2) = 2.0d0 * volW * ( sites(1:2,i) - cent(1:2) )
         gnnz = gnnz + 2
      end if
            


      do ed = 1,W%n-1
         if ( debug ) then
               write(*,*)
               write(*,*) '===================='
               write(*,*) 'Processing edge ed = ',ed
               write(*,*) 'Its inner vertex, named ',ed,' as well, corresponds to the point ',W%v(1:2,ed)
               write(*,*) 'Its outer vertex, named ',ednext,'as well, corresponds to the point ',W%v(1:2,ednext)
               write(*,*) '===================='
               write(*,*)
         end if

         if ( ed .eq. W%n-1 ) then
               ednext = 1
         else 
               ednext = ed + 1
         end if

         if ( debug ) then
               write(*,*) 'Next: ',ednext
               write(*,*) 'It corresponds to the point',W%v(1:2,ednext)
         end if

         diff1(1:2) = W%v(1:2,ed) - sites(1:2,i)
         diff2(1:2) = W%v(1:2,ednext) - sites(1:2,i)
         call det(W%v(1:2,ed), sites(1:2,i), W%v(1:2,ednext),vdet)

         if ( debug ) then
               write(*,*)
               write(*,*)'v_E - a_i = ',diff1(1:2)
               write(*,*)'w_E - a_i = ',diff2(1:2)
               write(*,*)'det = ',vdet
               write(*,*)
         end if 

         f = f  + ( abs( vdet )/12.0d0 ) * (3.0d0 * sum( diff1(1:2)**2 ) + &
               3.0d0 * dot_product( diff1(1:2),diff2(1:2) ) + sum( diff2(1:2)**2 ) )

         if ( compgrad ) then

            if ( W%e(ed) .lt. 0 ) then
               
               sizeed = norm2( W%v(1:2,ed) - W%v(1:2,ednext) )
               mu = sizeed / norm2( sites(1:2,i)- sites(1:2,abs(W%e(ed))) )

               if ( debug ) then
                  write(*,*)
                  write(*,*)'sizeed = ',sizeed
                  write(*,*)'mu = ',mu
                  write(*,*)
               end if

               p(1:2) = 0.5d0 * ( W%v(1:2,ed) + W%v(1:2,ednext) )
               auxc(1:6) = (/ sum( W%v(1:2,ednext)**2 )/6.0d0, dot_product( W%v(1:2,ednext),W%v(1:2,ed) )/3.0d0, &
                           sum( W%v(1:2,ed)**2 )/6.0d0, -dot_product( W%v(1:2,ednext),sites(1:2,i) )/3.0d0, &
                           -dot_product( W%v(1:2,ed),sites(1:2,i) )/3.0d0, sum( sites(1:2,i)**2 ) /)
               auxv(1:2,1) = W%v(1:2,ednext) + p(1:2)
               auxv(1:2,2) = W%v(1:2,ed) + p(1:2)
               zeta(1:2) = auxc(1)*( auxv(1:2,1) - 2.0d0*sites(1:2,i) ) + auxc(2)*( p(1:2) - sites(1:2,i) ) + &
                           auxc(3)*( auxv(1:2,2) - 2.0d0*sites(1:2,i) ) + auxc(4)*( auxv(1:2,1) + p(1:2) - 3.0d0*sites(1:2,i) ) + &
                           auxc(5)*( auxv(1:2,2) + p(1:2) - 3.0d0*sites(1:2,i) ) + auxc(6)*( p(1:2) - sites(1:2,i) )

               if ( debug ) then
                  write(*,*) 'zeta = ',zeta
               end if

               gval(1:2) = gval(1:2) + mu * zeta(1:2)

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*)'In compG: Increase gnnzmax to at least ',gnnz+2
                  flag = -1
                  return
               end if

               beta(1:2) = auxc(1)*( auxv(1:2,1) - 2.0d0*sites(1:2,abs(W%e(ed))) ) + &
                           auxc(2)*( p(1:2) - sites(1:2,abs(W%e(ed))) ) + &
                           auxc(3)*( auxv(1:2,2) - 2.0d0*sites(1:2,abs(W%e(ed))) ) + & 
                           auxc(4)*( auxv(1:2,1) + p(1:2) - 3.0d0*sites(1:2,abs(W%e(ed))) ) + &
                           auxc(5)*( auxv(1:2,2) + p(1:2) - 3.0d0*sites(1:2,abs(W%e(ed))) ) + &
                           auxc(6)*( p(1:2) - sites(1:2,abs(W%e(ed))) )
               if ( debug ) then
                  write(*,*) 'beta = ',beta
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - mu * beta(1:2)
               gnnz = gnnz + 2
            else
               if ( debug ) then
                  write(*,*) 'This is a boundary edge.'
                  write(*,*) 'So, its vertices v_E and w_E are ignored.'
               end if
            end if
         end if
      end do

   end subroutine compGalternative1

   ! ******************************************************************
   ! ******************************************************************

   subroutine compallJ(n,x,f,debug,flag,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      logical, intent(in) :: debug
      integer, intent(in) :: n
      integer, intent(out) :: flag
      real(kind=8), intent(out) :: f(0:3,n/2)
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
  
      ! PARAMETERS
      logical, parameter :: compgrad = .false.
      integer, parameter :: gtmpnnzmax = 1
      
      ! LOCAL SCALARS
      integer :: allocerr,gtmpnnz,i,ierror,k,maxnv,nsites,ntriag,status
      real(kind=8) :: ftmp
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: V
      type(Polygon) :: P,W
      
      ! LOCAL ARRAYS
      integer :: gtmpind(gtmpnnzmax),indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: fden(n/2),gden(2,n/2),gtmpval(gtmpnnzmax),sites(2,n/2),vorxy(2,2*(n/2))
   
          
      flag = 0
      
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
          
      if ( nsites .lt. 3 ) then
         write(*,*) 'In evalfg: nsites must be at least 3!'
         flag = - 1
         return
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
      
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In evalfg: dtris2 returned an error different from 0 and 225: ',ierror
         flag = - 1
         return
      end if
    
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In evalfg: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         flag = - 1
         return
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In evalfg: triangle_circumcenter_2d returned NaN.'
            flag = - 1
            return
         end if
      end do

      ! Compute density function at sites if required

      do i = 1, nsites
         call fgcompJ4(sites(1:2,i),fden(i),gden(1:2,i),compgrad,pdataptr)
      end do
      
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )
        
      allocate(P%v(2,maxnv),P%e(maxnv-1),stat=allocerr)

      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Allocation error.'
         flag = - 1
         return
      end if

      P%n = pdata%nvpols(1)
      P%deg = .false.
         
      P%v(1,1:P%n) = pdata%xpol(1,1:P%n)
      P%v(2,1:P%n) = pdata%ypol(1,1:P%n)
      P%e(1:P%n-1) = pdata%poledges(1,1:P%n-1)
  
      do i = 1,nsites
         if ( debug ) then
            write(*,*)
            write(*,*) 'SITE i = ',i
            write(*,*) 'Generator: ',sites(1:2,i)
            write(*,*)
         end if
  
         call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,V,status)
         if ( status .ne. 0 ) then
            write(*,*) 'In evalfg: Error in voronoi_cell, status = ',status
            flag = - 1
            return
         end if
  
         
         if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),P,V,W)
         else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
            call VorCellInterConvPol_Collinear(i,nsites,sites,P,W)
         end if

         call compJ4(nsites,sites,i,W,pdata%volA,fden,gden,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
         if ( flag .ne. 0 ) return
         
         f(3,i) = ftmp

         if ( W%deg ) then
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i, ' IS EMPTY.'
            end if
        
            f(0,i) = 1.0d0
            f(1,i) = 1.0d0
            f(2,i) = 1.0d0
            ! f(3,i) = 1.0d0

        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
  
            call compJ1(nsites,sites,i,W,pdata%volA,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return
  
            f(0,i) = ftmp
  
            call compJ2(nsites,sites,i,W,pdata%J2ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            f(1,i) = ftmp
                        
            call compJ3(nsites,sites,i,W,pdata%J3ctol,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return
  
            f(2,i) = ftmp

         end if

         call voronoi_cell_destroy(V)
  
         call polygon_destroy(W)
      end do
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
    end subroutine compallJ

   ! ******************************************************************
   ! ******************************************************************

   subroutine compJ1(nsites,sites,i,W,volA,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad,debug
      integer, intent(in) :: gnnzmax,i,nsites
      integer, intent(out) :: gnnz
      integer, intent(inout) :: flag
      real(kind=8), intent(in) :: volA
      real(kind=8), intent(out) :: f
      type(Polygon), intent(in) :: W
        
      ! ARRAY ARGUMENTS
      integer, intent(out) :: gind(gnnzmax)
      real(kind=8), intent(in) :: sites(2,nsites)
      real(kind=8), intent(out) :: gval(gnnzmax)
  
      ! LOCAL SCALARS
      integer :: ed,ednext,j
      real(kind=8) :: mu,vol,xi
          
      ! LOCAL ARRAYS
      real(kind=8) :: p(2)
  
            
      call polygon_area_2d(W%n-1,W%v(1:2,1:W%n-1),vol)

      xi = nsites / volA
  
      if ( debug ) then
         write(*,*)
         write(*,*) 'xi = ',xi
         write(*,*)
      end if
         
      f = xi * vol - 1.0d0
  
      if ( .not. compgrad ) then
         return
      end if
      
      if ( 2 .gt. gnnzmax ) then
         write(*,*) 'In compJ1: Increase gnnzmax to at least ',2
         flag = - 1
         return
      end if
      
      gind(1:2) = (/ 2*i-1, 2*i /)
      gval(1:2) = 0.0d0
      gnnz = 2
      
      if ( .not. W%deg ) then
         if ( debug ) then
            write(*,*)
            write(*,*) 'POLYGON Wij with i = ',i,' and j = ',j
            write(*,*) 'Number of vertices: ',W%n - 1
            write(*,*)
         end if
            
         do ed = 1,W%n - 1
            if ( debug ) then
               write(*,*)
               write(*,*) '===================='
               write(*,*) 'Processing edge ed = ',ed
               write(*,*) 'Its inner vertex, named ',ed,' as well, corresponds to the point ',W%v(1:2,ed)
               write(*,*) 'Its outer vertex, named ',ednext,' as well, corresponds to the point ',W%v(1:2,ednext)
               write(*,*) '===================='
               write(*,*)
            end if
            
            if ( W%e(ed) .lt. 0 ) then
               
               if ( ed .eq. W%n - 1 ) then
                  ednext = 1
               else
                  ednext = ed + 1
               end if
               
               if ( debug ) then
                  write(*,*) 'Next: ',ednext
                  write(*,*) 'It corresponds to the point ',W%v(1:2,ednext)
               end if
               
               mu = norm2( W%v(1:2,ed) - W%v(1:2,ednext) ) / norm2( sites(1:2,i) - sites(1:2,abs(W%e(ed))) )
                  
               if ( debug ) then
                  write(*,*) 'mu = ',mu
               end if

               p(1:2) = 0.5d0 * ( W%v(1:2,ed) + W%v(1:2,ednext) )
               
               gval(1:2) = gval(1:2) + xi * mu * ( p(1:2) - sites(1:2,i) )
                  
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ1: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
         
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * mu * ( p(1:2) - sites(1:2,abs(W%e(ed))) )
               gnnz = gnnz + 2

            else
               if ( debug ) then
                  write(*,*) 'This is a fixed edge.'
                  write(*,*) 'So, its vertices v_E and w_E are ignored.'
               end if
            end if
         end do

      else
         if ( debug ) then
            write(*,*)
            write(*,*) 'POLYGON Wi with i = ',i,' IS EMPTY.'
            write(*,*)
         end if
      end if
      
   end subroutine compJ1
  
    ! ******************************************************************
    ! ******************************************************************

   subroutine compJ2(nsites,sites,i,W,ctol,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad,debug
      integer, intent(in) :: gnnzmax,i,nsites
      integer, intent(out) :: gnnz
      integer, intent(inout) :: flag
      real(kind=8), intent(in) :: ctol
      real(kind=8), intent(out) :: f
      type(Polygon), intent(in) :: W
        
      ! ARRAY ARGUMENTS
      integer, intent(out) :: gind(gnnzmax)
      real(kind=8), intent(in) :: sites(2,nsites)
      real(kind=8), intent(out) :: gval(gnnzmax)
  
      ! LOCAL SCALARS
      integer :: ed,edprev,ednext
      real(kind=8) :: glcnst,mu,perim,phia,phib,phic,sizeed,tmp,xi
          
      ! LOCAL ARRAYS
      real(kind=8) :: gradphi(2),M(2,2),tang(2)
  
      !write(*,*) 'Site i = ',i
      perim = 0.0d0
      do ed = 1,W%n - 1
         if ( ed .lt. W%n - 1 ) then
            ednext = ed + 1
         else
            ednext = 1
         end if
  
         sizeed = norm2( W%v(1:2,ednext) - W%v(1:2,ed) )
         perim = perim + sizeed
  
        !write(*,*) 'Edge ',ed,' size = ',sizeed 
      end do
     !write(*,*) 'Average size = ',perim / ( W%n - 1 )
  
      if ( debug ) then
         write(*,*)
         write(*,*) 'Perimeter = ',perim
         write(*,*)
      end if
         
      xi = 2.0d0 / perim
            
      if ( debug ) then
         write(*,*)
         write(*,*) 'xi = ',xi
         write(*,*)
      end if
  
      f = 0.0d0
  
      if ( compgrad ) then
         glcnst = 0.0d0
         do ed = 1,W%n - 1
            if ( ed .lt. W%n - 1 ) then
               ednext = ed + 1
            else
               ednext = 1
            end if
            
            sizeed = norm2( W%v(1:2,ednext) - W%v(1:2,ed) )
            glcnst = glcnst + min( 0.0d0, ( W%n - 1 ) * sizeed / perim - ctol ) * sizeed / perim
         end do
         
         if ( debug ) then
            write(*,*)
            write(*,*) 'Global constant = ',glcnst
            write(*,*)
         end if
         
         gnnz = 0
         
         if ( gnnz + 2 .gt. gnnzmax ) then
            write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
            flag = - 1
            return
         end if
         
         gind(gnnz+1:gnnz+2) = (/ 2*i-1, 2*i /)
         gval(gnnz+1:gnnz+2) = 0.0d0
         gnnz = gnnz + 2
      end if
      
      do ed = 1,W%n - 1
         if ( debug ) then
            write(*,*)
            write(*,*) '===================='
            write(*,*) 'Processing edge ed = ',ed
            write(*,*) 'Its inner vertex, named ',ed,' as well, corresponds to the point ',W%v(1:2,ed)
            write(*,*) '===================='
            write(*,*)
         end if
         
         ! Processing edge (v_E,w_E)
                     
         if ( ed .eq. 1 ) then
            edprev = W%n - 1
            ednext = ed + 1
         else if ( ed .eq. W%n - 1 ) then
            edprev = ed - 1
            ednext = 1
         else
            edprev = ed - 1
            ednext = ed + 1
         end if
                    
         if ( debug ) then
            write(*,*) 'Previous: ',edprev,' next: ',ednext
            write(*,*) 'They correspond to the points ',W%v(1:2,edprev),' and ',W%v(1:2,ednext),' respectively.'
         end if
                  
         sizeed = norm2( W%v(1:2,ednext) - W%v(1:2,ed) )
  
         if ( debug ) then
            write(*,*) 'Edge size = ',sizeed
         end if
  
         tmp = min( 0.0d0, ( W%n - 1 ) * sizeed / perim - ctol )
  
         f = f + tmp ** 2 / ( W%n - 1 )
  
         if ( .not. compgrad ) then
            cycle
         end if
         
         mu = glcnst - tmp
               
         if ( debug ) then
            write(*,*) 'mu associated wit the edge = ',mu
         end if
                  
         if ( mu .eq. 0.0d0 ) then
            cycle
         end if
         
         tang(1:2) = ( W%v(1:2,ednext) - W%v(1:2,ed) ) / sizeed
  
         if ( debug ) then
            write(*,*) 'Tangent vector: ',tang(1:2)
         end if
  
         ! Starting point v_E
  
         if ( debug ) then
            write(*,*)
            write(*,*) 'Processing v_E'
            write(*,*)
         end if
  
         if ( W%e(ed) .lt. 0 .and. W%e(edprev) .lt. 0 ) then
            if ( debug ) then
               write(*,*) 'We are in the case |P(v_E)|=3. ', &
                    'Neighbour cells are: ',abs(W%e(edprev)),' and ',abs(W%e(ed))
            end if
                  
            call matrixint(sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ed),M,flag)
            if ( flag .ne. 0 ) return
                        
            gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
                     
            call matrixint(sites(1:2,i),sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),M,flag)
            if ( flag .ne. 0 ) return
  
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
            gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
               
            call matrixint(sites(1:2,abs(W%e(ed))),sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),M,flag)
            if ( flag .ne. 0 ) return
  
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
            gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
               
         else if ( W%e(ed) .lt. 0 ) then
            if ( debug ) then
               write(*,*) 'We are in the case |P(v_E)|=2 (type 1). Neighbour cell is: ',abs(W%e(ed))
            end if
                  
            call line_exp2imp_2d(W%v(1:2,ed),W%v(1:2,edprev),phia,phib,phic)
  
            gradphi(1:2) = (/ phia, phib /)
            
            if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
               gradphi(1:2) = - gradphi(1:2)
            end if
      
            call matrixbd(sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ed),gradphi,M,flag)
            if ( flag .ne. 0 ) return
  
            gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
                     
            call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),gradphi,M,flag)
            if ( flag .ne. 0 ) return
               
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
            gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
               
         else if ( W%e(edprev) .lt. 0 ) then
            if ( debug ) then
               write(*,*) 'We are in the case |P(v_E)|=2 (type 2). Neighbour cell is: ',abs(W%e(edprev))
            end if
                  
            call line_exp2imp_2d(W%v(1:2,ed),W%v(1:2,ednext),phia,phib,phic)
  
            gradphi(1:2) = (/ phia, phib /)
                  
            if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
               gradphi(1:2) = - gradphi(1:2)
            end if
      
            call matrixbd(sites(1:2,abs(W%e(edprev))),sites(1:2,i),W%v(1:2,ed),gradphi,M,flag)
            if ( flag .ne. 0 ) return
  
            gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
                     
            call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),gradphi,M,flag)
            if ( flag .ne. 0 ) return
                        
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
            gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
         end if
                     
         ! Ending point w_E
                    
         if ( debug ) then
            write(*,*)
            write(*,*) 'Processing w_E'
            write(*,*)
         end if
            
         if ( W%e(ed) .lt. 0 .and. W%e(ednext) .lt. 0 ) then
            if ( debug ) then
               write(*,*) 'We are in the case |P(w_E)|=3. ', &
                    'Neighbour cells are: ',abs(W%e(ed)),' and ',abs(W%e(ednext))
            end if
                  
            call matrixint(sites(1:2,abs(W%e(ed))),sites(1:2,abs(W%e(ednext))),sites(1:2,i),W%v(1:2,ednext),M,flag)
            if ( flag .ne. 0 ) return
               
            gval(1:2) = gval(1:2) - xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
  
            call matrixint(sites(1:2,i),sites(1:2,abs(W%e(ed))),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),M,flag)
            if ( flag .ne. 0 ) return
                   
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
            gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
            
            call matrixint(sites(1:2,abs(W%e(ednext))),sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),M,flag)
            if ( flag .ne. 0 ) return
                   
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
            gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
            
         else if ( W%e(ed) .lt. 0 ) then 
            if ( debug ) then
               write(*,*) 'We are in the case |P(w_E)|=2 (type 1). Neighbour cell is: ',abs(W%e(ed))
            end if
                  
            call line_exp2imp_2d(W%v(1:2,ednext),W%v(1:2,ednext+1),phia,phib,phic)
  
            gradphi(1:2) = (/ phia, phib /)
            
            if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
               gradphi(1:2) = - gradphi(1:2)
            end if
      
            call matrixbd(sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ednext),gradphi,M,flag)
            if ( flag .ne. 0 ) return
  
            gval(1:2) = gval(1:2) - xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
                     
            call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),gradphi,M,flag)
            if ( flag .ne. 0 ) return
                     
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
            gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
            
         else if ( W%e(ednext) .lt. 0 ) then 
            if ( debug ) then
               write(*,*) 'We are in the case |P(w_E)|=2 (type 2). Neighbour cell is: ',abs(W%e(ednext))
            end if
            
            call line_exp2imp_2d(W%v(1:2,ednext),W%v(1:2,ed),phia,phib,phic)
               
            gradphi(1:2) = (/ phia, phib /)
                  
            if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
               gradphi(1:2) = - gradphi(1:2)
            end if
               
            call matrixbd(sites(1:2,abs(W%e(ednext))),sites(1:2,i),W%v(1:2,ednext),gradphi,M,flag)
            if ( flag .ne. 0 ) return
  
            gval(1:2) = gval(1:2) - xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            
            call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),gradphi,M,flag)
            if ( flag .ne. 0 ) return
                     
            if ( gnnz + 2 .gt. gnnzmax ) then
               write(*,*) 'In compJ3: Increase gnnzmax to at least ',gnnz+2
               flag = - 1
               return
            end if
          
            gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
            gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), tang(1:2) )
            gnnz = gnnz + 2
         end if
      end do
         
   end subroutine compJ2

   ! ******************************************************************
   ! ******************************************************************

   subroutine compJ3(nsites,sites,i,W,ctol,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad,debug
      integer, intent(in) :: gnnzmax,i,nsites
      integer, intent(out) :: gnnz
      integer, intent(inout) :: flag
      real(kind=8), intent(in) :: ctol
      real(kind=8), intent(out) :: f
      type(Polygon), intent(in) :: W
        
      ! ARRAY ARGUMENTS
      integer, intent(out) :: gind(gnnzmax)
      real(kind=8), intent(in) :: sites(2,nsites)
      real(kind=8), intent(out) :: gval(gnnzmax)
  
      ! LOCAL SCALARS
      integer :: ed,edprev,edpp,ednext,nmov
      real(kind=8) :: ang,avgang,glcnst,mu,phia,phib,phic,sizeed,sizeedprev,tmp,xi
          
      ! LOCAL ARRAYS
      real(kind=8) :: gradphi(2),M(2,2),normaled(2),normaledprev(2),taued(2),tauedprev(2),zeta(2)
  
      !write(*,*) 'Site i = ',i
      nmov = 0
      avgang = 0.0d0
      do ed = 1,W%n - 1
         if ( ed .eq. 1) then
            edprev = W%n - 1
            ednext = ed + 1
         else if ( ed .eq. W%n - 1 ) then
            edprev = ed - 1
            ednext = 1
         else
            edprev = ed - 1
            ednext = ed + 1
         end if
  
         ! If none is negative, the angle is fixed
         if ( W%e(ed) .lt. 0 .or. W%e(edprev) .lt. 0 ) then
            tauedprev(1:2) = W%v(1:2,ed) - W%v(1:2,edprev)
            taued(1:2)     = W%v(1:2,ednext) - W%v(1:2,ed)
  
            sizeedprev = norm2( tauedprev(1:2) )
            sizeed     = norm2( taued(1:2) )
         
            ang = acos( max( -1.0d0, min( dot_product( - tauedprev(1:2) / sizeedprev, taued(1:2) / sizeed ), 1.0d0 ) ) )
            !write(*,*) 'vertex = ',ed,' angle = ',ang
  
            if ( debug ) then
               write(*,*) 'vertex = ',ed,' angle = ',ang
            end if
         
            nmov = nmov + 1
            avgang = avgang + ang
         end if
      end do
      avgang = avgang / nmov
      !write(*,*) 'Average angle = ',avgang
  
      if ( debug ) then
         write(*,*)
         write(*,*) 'The polygon has nmov = ',nmov,' non-fixed angles.'
         write(*,*) 'Average angle = ',avgang
         write(*,*)
      end if
  
      xi = ( 2.0d0 / ( avgang * nmov ) )
             
      if ( debug ) then
         write(*,*)
         write(*,*) 'xi = ',xi
         write(*,*)
      end if
  
      f = 0.0d0
  
      if ( compgrad ) then
         glcnst = 0.0d0
         do ed = 1,W%n - 1
            if ( ed .eq. 1) then
               edprev = W%n - 1
               ednext = ed + 1
            else if ( ed .eq. W%n - 1 ) then
               edprev = ed - 1
               ednext = 1
            else
               edprev = ed - 1
               ednext = ed + 1
            end if
  
            ! If none is negative, the angle is fixed
            if ( W%e(ed) .lt. 0 .or. W%e(edprev) .lt. 0 ) then
               tauedprev(1:2) = W%v(1:2,ed) - W%v(1:2,edprev)
               taued(1:2)     = W%v(1:2,ednext) - W%v(1:2,ed)
               
               sizeedprev = norm2( tauedprev(1:2) )
               sizeed     = norm2( taued(1:2) )
               
               ang = acos( max( -1.0d0, min( dot_product( - tauedprev(1:2) / sizeedprev, taued(1:2) / sizeed ), 1.0d0 ) ) )
               
               glcnst = glcnst + ( ang / avgang ) * min( 0.0d0, ang / avgang - ctol ) / nmov
            end if
         end do
      
         gnnz = 0
      
         if ( gnnz + 2 .gt. gnnzmax ) then
            write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
            flag = - 1
            return
         end if
      
         gind(gnnz+1:gnnz+2) = (/ 2*i-1, 2*i /)
         gval(gnnz+1:gnnz+2) = 0.0d0
         gnnz = gnnz + 2
      end if
      
      do ed = 1,W%n - 1
         if ( debug ) then
            write(*,*)
            write(*,*) '===================='
            write(*,*) 'Processing edge ed = ',ed
            write(*,*) 'Its inner vertex, named ',ed,' as well, corresponds to the point ',W%v(1:2,ed)
            write(*,*) '===================='
            write(*,*)
         end if
         
         ! Processing edge (v_E,w_E)
                     
         if ( ed .eq. 1) then
            edprev = W%n - 1
            ednext = ed + 1
         else if ( ed .eq. W%n - 1 ) then
            edprev = ed - 1
            ednext = 1
         else
            edprev = ed - 1
            ednext = ed + 1
         end if
                    
         if ( debug ) then
            write(*,*) 'Previous: ',edprev,' next: ',ednext
            write(*,*) 'They correspond to the points ',W%v(1:2,edprev),' and ',W%v(1:2,ednext),' respectively.'
         end if
  
         ! If none is negative, the angle is fixed
         if ( W%e(ed) .lt. 0 .or. W%e(edprev) .lt. 0 ) then
            
            tauedprev(1:2) = W%v(1:2,ed) - W%v(1:2,edprev)
            taued(1:2)     = W%v(1:2,ednext) - W%v(1:2,ed)
  
            sizeedprev = norm2( tauedprev(1:2) )
            sizeed     = norm2( taued(1:2) )
         
            tauedprev(1:2) = tauedprev(1:2) / sizeedprev
            taued(1:2)     = taued(1:2) / sizeed
  
            ang = acos( max( -1.0d0, min( dot_product( - tauedprev(1:2), taued(1:2) ), 1.0d0 ) ) )
  
            tmp = min( 0.0d0, ang / avgang - ctol )
  
            f = f + tmp ** 2 / nmov
  
            if ( .not. compgrad ) then
               cycle
            end if
            
            mu = glcnst - tmp
  
            if ( debug ) then
               write(*,*) 'mu = ',mu
            end if
  
            if ( mu .eq. 0.0d0 ) then
               cycle
            end if
         
            call line_exp_normal_2d(W%v(1:2,ednext),W%v(1:2,ed),normaled)
            call line_exp_normal_2d(W%v(1:2,ed),W%v(1:2,edprev),normaledprev)
  
            if ( debug ) then
               write(*,*) 'Normal to ed     = ',normaled(1:2)
               write(*,*) 'Normal to edprev = ',normaledprev(1:2)
            end if
  
            ! Starting point v_E
  
            if ( debug ) then
               write(*,*)
               write(*,*) 'Processing v_E = w_{\hat E}'
               write(*,*)
            end if
  
          ! zeta(1:2) = ( dot_product( normaled(1:2), tauedprev(1:2) ) / sizeed ) * normaled(1:2) - &
          !             ( dot_product( normaledprev(1:2), taued(1:2) ) / sizeedprev ) * normaledprev(1:2)
            zeta(1:2) = normaled(1:2) / sizeed + normaledprev(1:2) / sizeedprev
  
            if ( debug ) then
               write(*,*) 'zeta = ',zeta(1:2)
            end if
                  
            if ( W%e(ed) .lt. 0 .and. W%e(edprev) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(v_E)|=3. ', &
                       'Neighbour cells are: ',abs(W%e(edprev)),' and ',abs(W%e(ed))
               end if
                  
               call matrixint(sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return
                        
               gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
                     
               call matrixint(sites(1:2,i),sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return
  
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
  
               call matrixint(sites(1:2,abs(W%e(ed))),sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return
  
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
  
            else if ( W%e(ed) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(v_E)|=2 (type 1). Neighbour cell is: ',abs(W%e(ed))
               end if
               
               call line_exp2imp_2d(W%v(1:2,ed),W%v(1:2,edprev),phia,phib,phic)
  
               gradphi(1:2) = (/ phia, phib /)
            
               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = - gradphi(1:2)
               end if
      
               call matrixbd(sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return
  
               gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
                     
               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return
               
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
               
            else if ( W%e(edprev) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(v_E)|=2 (type 2). Neighbour cell is: ',abs(W%e(edprev))
               end if
                  
               call line_exp2imp_2d(W%v(1:2,ed),W%v(1:2,ednext),phia,phib,phic)
  
               gradphi(1:2) = (/ phia, phib /)
                  
               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = - gradphi(1:2)
               end if
      
               call matrixbd(sites(1:2,abs(W%e(edprev))),sites(1:2,i),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return
  
               gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
                     
               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return
                        
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
            end if
                     
            ! Ending point w_E
                    
            if ( debug ) then
               write(*,*)
               write(*,*) 'Processing w_E'
               write(*,*)
            end if
                  
          ! zeta(1:2) = ( dot_product( normaled(1:2), tauedprev(1:2) ) / sizeed ) * normaled(1:2)
            zeta(1:2) = normaled(1:2) / sizeed
  
            if ( debug ) then
               write(*,*) 'zeta = ',zeta(1:2)
            end if
                  
            if ( W%e(ed) .lt. 0 .and. W%e(ednext) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=3. ', &
                       'Neighbour cells are: ',abs(W%e(ed)),' and ',abs(W%e(ednext))
               end if
                  
               call matrixint(sites(1:2,abs(W%e(ed))),sites(1:2,abs(W%e(ednext))),sites(1:2,i),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0 ) return
                     
               gval(1:2) = gval(1:2) - xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
               
               call matrixint(sites(1:2,i),sites(1:2,abs(W%e(ed))),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0 ) return
                   
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
               gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
            
               call matrixint(sites(1:2,abs(W%e(ednext))),sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0 ) return
                   
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
            
            else if ( W%e(ed) .lt. 0 ) then 
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=2 (type 1). Neighbour cell is: ',abs(W%e(ed))
               end if
                  
               call line_exp2imp_2d(W%v(1:2,ednext),W%v(1:2,ednext+1),phia,phib,phic)
  
               gradphi(1:2) = (/ phia, phib /)
            
               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = - gradphi(1:2)
               end if
      
               call matrixbd(sites(1:2,abs(W%e(ed))),sites(1:2,i),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return
  
               gval(1:2) = gval(1:2) - xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
                     
               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return
                     
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
               
            else if ( W%e(ednext) .lt. 0 ) then 
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=2 (type 2). Neighbour cell is: ',abs(W%e(ednext))
               end if
               
               call line_exp2imp_2d(W%v(1:2,ednext),W%v(1:2,ed),phia,phib,phic)
               
               gradphi(1:2) = (/ phia, phib /)
                  
               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = - gradphi(1:2)
               end if
      
               call matrixbd(sites(1:2,abs(W%e(ednext))),sites(1:2,i),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return
  
               gval(1:2) = gval(1:2) - xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
                     
               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return
                     
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
               gval(gnnz+1:gnnz+2) = - xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
            end if
  
            ! Starting point point v_{\hat E}
            
            if ( edprev .eq. 1 ) then
               edpp = W%n - 1
            else
               edpp = edprev - 1
            end if
               
            if ( debug ) then
               write(*,*)
               write(*,*) 'Processing v_{\hat E}'
               write(*,*)
            end if
                  
          ! zeta(1:2) = ( dot_product( normaledprev(1:2), taued(1:2) ) / sizeedprev ) * normaledprev(1:2)
            zeta(1:2) = - normaledprev(1:2) / sizeedprev
  
            if ( debug ) then
               write(*,*) 'zeta = ',zeta(1:2)
            end if
                  
            if ( W%e(edprev) .lt. 0 .and. W%e(edpp) .lt. 0 ) then
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=3. ', &
                       'Neighbour cells are: ',abs(W%e(edprev)),' and ',abs(W%e(edpp))
               end if
                  
               call matrixint(sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(edpp))),sites(1:2,i),W%v(1:2,edprev),M,flag)
               if ( flag .ne. 0 ) return
                     
               gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
               
               call matrixint(sites(1:2,i),sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(edpp))),W%v(1:2,edprev),M,flag)
               if ( flag .ne. 0 ) return
                   
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edpp))-1, 2*abs(W%e(edpp)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
                     
               call matrixint(sites(1:2,abs(W%e(edpp))),sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,edprev),M,flag)
               if ( flag .ne. 0 ) return
                   
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
                     
            else if ( W%e(edprev) .lt. 0 ) then 
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=2 (type 1). Neighbour cell is: ',abs(W%e(edprev))
               end if
                  
               call line_exp2imp_2d(W%v(1:2,edpp),W%v(1:2,edprev),phia,phib,phic)
  
               gradphi(1:2) = (/ phia, phib /)
            
               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = - gradphi(1:2)
               end if
      
               call matrixbd(sites(1:2,abs(W%e(edprev))),sites(1:2,i),W%v(1:2,edprev),gradphi,M,flag)
               if ( flag .ne. 0 ) return
  
               gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
                  
               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,edprev),gradphi,M,flag)
               if ( flag .ne. 0 ) return
                     
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
            
            else if ( W%e(edpp) .lt. 0 ) then 
               if ( debug ) then
                  write(*,*) 'We are in the case |P(w_E)|=2 (type 2). Neighbour cell is: ',abs(W%e(edpp))
               end if
               
               call line_exp2imp_2d(W%v(1:2,edprev),W%v(1:2,ed),phia,phib,phic)
               
               gradphi(1:2) = (/ phia, phib /)
                  
               if ( gradphi(1) * sites(1,i) + gradphi(2) * sites(2,i) + phic .gt. 0.0d0 ) then
                  gradphi(1:2) = - gradphi(1:2)
               end if
      
               call matrixbd(sites(1:2,abs(W%e(edpp))),sites(1:2,i),W%v(1:2,edprev),gradphi,M,flag)
               if ( flag .ne. 0 ) return
  
               gval(1:2) = gval(1:2) + xi * mu * matmul( transpose( M(1:2,1:2) ),zeta(1:2) )
               
               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(edpp))),W%v(1:2,edprev),gradphi,M,flag)
               if ( flag .ne. 0 ) return
                     
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ4: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
          
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edpp))-1, 2*abs(W%e(edpp)) /)
               gval(gnnz+1:gnnz+2) = xi * mu * matmul( transpose( M(1:2,1:2) ), zeta(1:2) )
               gnnz = gnnz + 2
            end if
         else
            if ( debug ) then
               write(*,*)
               write(*,*) 'The inner vertex of this edge has a fixed angle; so it will  be ignored.'
               write(*,*)
            end if
         end if
      end do
         
   end subroutine compJ3

   ! ******************************************************************
   ! ******************************************************************

   subroutine compJ4(nsites,sites,i,W,volA,fden,gden,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)

      implicit none
      
      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad,debug
      integer, intent(in) :: gnnzmax,i,nsites
      integer, intent(out) :: gnnz
      integer, intent(inout) :: flag
      real(kind=8), intent(in) :: volA
      real(kind=8), intent(out) :: f
      type(Polygon), intent(in) :: W

      ! ARRAY ARGUMENTS
      integer, intent(out) :: gind(gnnzmax)
      real(kind=8), intent(in) :: fden(nsites),gden(2,nsites),sites(2,nsites)
      real(kind=8), intent(out) :: gval(gnnzmax)

      ! LOCAL SCALARS
      integer :: ed,ednext,j
      real(kind=8) :: mu,vol,xi
          
      ! LOCAL ARRAYS
      real(kind=8) :: p(2)

      call polygon_area_2d(W%n-1,W%v(1:2,1:W%n-1),vol)

      xi = nsites / volA
  
      if ( debug ) then
         write(*,*)
         write(*,*) 'xi = ',xi
         write(*,*)
      end if

      f = xi * vol - fden(i)
  
      if ( .not. compgrad ) then
         return
      end if
      
      if ( 2 .gt. gnnzmax ) then
         write(*,*) 'In compJ4: Increase gnnzmax to at least ',2
         flag = - 1
         return
      end if

      gind(1:2) = (/ 2*i-1, 2*i /)
      gval(1:2) = -gden(1:2,i)
      gnnz = 2

      if ( .not. W%deg ) then
         if ( debug ) then
            write(*,*)
            write(*,*) 'POLYGON Wij with i = ',i,' and j = ',j
            write(*,*) 'Number of vertices: ',W%n - 1
            write(*,*)
         end if
            
         do ed = 1,W%n - 1
            if ( debug ) then
               write(*,*)
               write(*,*) '===================='
               write(*,*) 'Processing edge ed = ',ed
               write(*,*) 'Its inner vertex, named ',ed,' as well, corresponds to the point ',W%v(1:2,ed)
               write(*,*) 'Its outer vertex, named ',ednext,' as well, corresponds to the point ',W%v(1:2,ednext)
               write(*,*) '===================='
               write(*,*)
            end if
            
            if ( W%e(ed) .lt. 0 ) then
               
               if ( ed .eq. W%n - 1 ) then
                  ednext = 1
               else
                  ednext = ed + 1
               end if
               
               if ( debug ) then
                  write(*,*) 'Next: ',ednext
                  write(*,*) 'It corresponds to the point ',W%v(1:2,ednext)
               end if
               
               mu = norm2( W%v(1:2,ed) - W%v(1:2,ednext) ) / norm2( sites(1:2,i) - sites(1:2,abs(W%e(ed))) )
                  
               if ( debug ) then
                  write(*,*) 'mu = ',mu
               end if

               p(1:2) = 0.5d0 * ( W%v(1:2,ed) + W%v(1:2,ednext) )
               
               gval(1:2) = gval(1:2) + xi * mu * ( p(1:2) - sites(1:2,i) )
                  
               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJ1: Increase gnnzmax to at least ',gnnz+2
                  flag = - 1
                  return
               end if
         
               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * mu * ( p(1:2) - sites(1:2,abs(W%e(ed))) )
               gnnz = gnnz + 2

            else
               if ( debug ) then
                  write(*,*) 'This is a fixed edge.'
                  write(*,*) 'So, its vertices v_E and w_E are ignored.'
               end if
            end if
         end do

      else
         if ( debug ) then
            write(*,*)
            write(*,*) 'POLYGON Wi with i = ',i,' IS EMPTY.'
            write(*,*)
         end if
      end if

   end subroutine compJ4

   subroutine fgcompJ4(x,fden,gden,compgrad,pdataptr)

      implicit none

      ! SCALAR ARGUMENTS
      logical, intent(in) :: compgrad
      real(kind=8), intent(out) :: fden
      type(c_ptr), optional, intent(in) :: pdataptr

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(2)
      real(kind=8), intent(out) :: gden(2)

      ! LOCAL SCALARS
      type(pdata_type), pointer :: pdata
      real(kind=8) :: delta,param,a,b,scale,z1,z2,faux
      real(kind=8), parameter :: pi = 4.0d0 * atan( 1.0d0 )

      ! LOCAL ARRAY
      real(kind=8) :: z(2)

      call c_f_pointer(pdataptr,pdata)

      !! Function 1 CROSS
      ! delta = pdata%c(1)/4.0d0

      ! if ( abs( x(1) - pdata%c(1) ) .le. delta .or. abs( x(2) - pdata%c(2) ) .le. delta   ) then
      !    fden = 0.25d0
      ! else
      !    fden = 1.075d0
      ! end if

      ! if ( compgrad ) gden(1:2) = 0.0d0

      !! Function 2
      ! fden = 2.5d0 - 2.0d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 4.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      ! !! Function 2.0 Does not work for the case of the circle circunscribed in the square
      ! fden = 20.0d0 - (19.0d0 + 0.5d0) * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 39.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.1
      
      ! fden = 10.0d0 - 9.0d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 18.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.2 CINS4 CONVX2

      fden = 5.0d0 - 4.99d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      if ( compgrad ) gden(1:2) = - 9.98d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.3

      ! fden = 5.0d0 - 4.999d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 9.998d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2


      !! Function 2.4

      ! fden = 4.5d0 - 4.49d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 8.98d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.5

      ! fden = 4.0d0 - 3.99d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 7.98d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.6

      ! fden = 3.5d0 - 3.49d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 6.98d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.7

      ! fden = 3.2d0 - 3.19d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 6.38d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.8 CINS3 CONVX CONVX1

      ! fden = 2.5d0 - 2.49d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 4.98d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 2.9 Does not work squarecircle

      ! fden = 2.5d0 - 2.49d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 16 ) / pdata%J4r ** 16

      ! if ( compgrad ) gden(1:2) = - 2.49d0 * 16.0d0 * ( x(1:2) - pdata%J4c(1:2) ) ** 15 / pdata%J4r ** 16
      
      !! Function 3 CINS1

      ! fden = 0.01d0 + 20.0d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = 40.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 3.1 CINS2

      ! fden = 0.5d0 + 20.0d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = 40.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 3.2 CINS

      ! fden = 0.001d0 + 20.0d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = 40.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 3.3

      ! fden = 1.0d0 + 20.0d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = 40.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 4 Does not work

      ! fden = 4.5d0 - 2.0d0 * sum( ( x(1:2) - pdata%J4c(1:2) ) ** 2 ) / pdata%J4r ** 2

      ! if ( compgrad ) gden(1:2) = - 4.0d0 * ( x(1:2) - pdata%J4c(1:2) ) / pdata%J4r ** 2

      !! Function 5 BANANA

      ! z(1:2) = -pdata%c(1:2)/5.0d0 + 0.40d0 * x(1:2)
      ! if ( ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) ** 2 + ( 0.25d0 * z(1) - 1.0d0 ) ** 2 - 1.0d0 .le. 0.0d0 ) then
      !    fden = 0.25d0
      ! else
      !    fden = 1.075d0
      ! end if
      
      ! if ( compgrad ) gden(1:2) = 0.0d0

      !! Function 6 Works, but does not nothing to the diagram
      ! fden = exp( -10.0d0 * (x(1)**2 + x(2)**2) )

      ! if ( compgrad ) gden(1:2) = (/ -20.0d0 * x(1) * fden, -20.0d0 * x(2) * fden /)

      !! Function 7 Works, but does not nothing to the diagram

      ! fden = exp( -2.0d0 * ( x(1) + x(2) ) )

      ! if ( compgrad ) gden(1:2) = -2.0d0 * (/ fden, fden /)

      !! Function 8

      ! fden = exp( -20.0d0 * ( x(1) ** 2 + x(2) ** 2 ) )+ 0.05d0 * (sin(pi * x(1)) ** 2) * (sin(pi * x(2)) ** 2) 

      ! if ( compgrad ) then 
      !    gden(1) =  -40.0d0 * x(1) * exp( -20.0d0 * ( x(1) ** 2 + x(2) ** 2 ) ) + &
      !                0.1d0 * pi * sin( pi * x(1) ) * cos( pi * x(1) ) * (sin(pi * x(2)) ** 2) 
      !    gden(2) =  -40.0d0 * x(2) * exp( -20.0d0 * ( x(1) ** 2 + x(2) ** 2 ) ) + &
      !    0.1d0 * pi * sin( pi * x(2) ) * cos( pi * x(2) ) * (sin(pi * x(1)) ** 2) 
      ! end if

      !! Function 9

      ! if ( ( x(1) - pdata%c(1) ) ** 2 + ( x(2) - pdata%c(2) ) ** 2  .le. 9.0d0 ) then 
      !    fden = 0.01d0
      ! else if ( 9.0d0 .lt. ( x(1) - pdata%c(1) ) ** 2 + ( x(2) - pdata%c(2) ) ** 2 .and. &
      !         ( x(1) - pdata%c(1) ) ** 2 + ( x(2) - pdata%c(2) ) ** 2  .le. 18.0d0  ) then
      !    fden = 0.1d0
      ! else if ( 18.0d0 .lt. (x(1) - pdata%c(1) ) ** 2 + ( x(2) - pdata%c(2) ) ** 2 .and. & 
      !         ( x(1) - pdata%c(1) ) ** 2 + ( x(2) - pdata%c(2) ) ** 2  .le. 27.0d0 ) then
      !    fden = 1.0d0
      ! else 
      !    fden = 5.0d0
      ! end if

      ! if ( compgrad ) gden(1:2) = 0.0d0

      !! Function 10 SIN1
      ! delta = 0.5d0*sqrt( 1000.0d0 )
      ! if ( abs( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) - delta ) .lt. 0.25d0*delta  ) then
      !    fden = 0.15d0
      ! else 
      !    fden = 1.175d0
      ! end if

      ! if ( compgrad ) gden(1:2) = 0.0d0

      !! Function 11 does not working
      ! fden = exp((x(1) - pdata%c(1))**2 + (x(2) - pdata%c(2))**2 ) + 1.0d0

      ! if ( compgrad ) gden(1:2) = 2.0d0 * fden * (/ x(1) - pdata%c(1), x(2) - pdata%c(2) /)

      !! Function 12
      ! delta = 0.5d0*sqrt( 1000.d0 )
      ! param = 4.9d0/delta**2
      ! fden = 0.1d0 + param*( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) - delta ) ** 2
      ! if ( compgrad ) then
      !     gden(1) = 2.0d0 * param * ( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) ) &
      !              * ( -0.6*2.0d0*pi*cos( 2.0d0*pi*x(1)/(2.0d0*delta) )/(2.0d0*delta) )
      !     gden(2) = 2.0d0 * param * ( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) )
      ! end if
      !! Function 12.1
      ! delta = 0.5d0*sqrt( 1000.d0 )
      ! param = 3.9d0/delta**2
      ! fden = 0.1d0 + param*( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) - delta ) ** 2
      ! if ( compgrad ) then
      !     gden(1) = 2.0d0 * param * ( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) ) &
      !              * ( -0.6*2.0d0*pi*cos( 2.0d0*pi*x(1)/(2.0d0*delta) )/(2.0d0*delta) )
      !     gden(2) = 2.0d0 * param * ( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) )
      ! end if
      !! Function 12.2 SIN2
      ! delta = 0.5d0*sqrt( 1000.d0 )
      ! param = 2.9d0/delta**2
      ! fden = 0.1d0 + param*( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) - delta ) ** 2
      ! if ( compgrad ) then
      !     gden(1) = 2.0d0 * param * ( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) ) &
      !              * ( -0.6*2.0d0*pi*cos( 2.0d0*pi*x(1)/(2.0d0*delta) )/(2.0d0*delta) )
      !     gden(2) = 2.0d0 * param * ( x(2) - 0.6d0*delta*sin( 2.0d0*pi*x(1)/(2.0d0*delta) ) )
      ! end if

      ! !! Function 13 
      ! z(1:2) = -pdata%c(1:2)/5.0d0 + 0.40d0 * x(1:2)
      ! b = 0.25d0
      ! a = 0.0625d0 * (2.4375d0 - b)
      ! fden = a * ( (z(2) - ( 0.25d0 * z(1) ) ** 2 ) ** 2 + ( 0.25d0 * z(1) - 1.0d0 ) ** 2 ) + b
      ! if ( compgrad ) then
      !    gden(1) = a * ( -0.1d0 * z(1) * ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) + 0.2d0 * (0.25d0 * z(1) - 1.0d0) )
      !    gden(2) = a * ( 0.8d0 *  z(2) - ( 0.25d0 * z(1) ) ** 2 ) 
      ! end if

      !! Function 13.1 BANANA1
      ! z(1:2) = -pdata%c(1:2)/5.0d0 + 0.40d0 * x(1:2)
      ! b = 0.25d0
      ! a = 0.0625d0 * (2.0d0 - b - 0.5625d0)
      ! fden = a * ( (z(2) - ( 0.25d0 * z(1) ) ** 2 ) ** 2 + ( 0.25d0 * z(1) - 1.0d0 ) ** 2 ) + b
      ! if ( compgrad ) then
      !    gden(1) = a * ( -0.1d0 * z(1) * ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) + 0.2d0 * (0.25d0 * z(1) - 1.0d0) )
      !    gden(2) = a * ( 0.8d0 * ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) ) 
      ! end if


      !! Function 13.6 and 13.8
      ! z(1:2) = -pdata%c(1:2)/5.0d0 + 0.40d0 * x(1:2)
      ! b = 0.25d0
      ! a = 0.0625d0 * (1.9d0 - b - 0.5625d0)  !! 13.6
      ! a = 0.0625d0 * (1.75d0 - b - 0.5625d0)  !! 13.8
      ! fden = a * ( (z(2) - ( 0.25d0 * z(1) ) ** 2 ) ** 2 + ( 0.25d0 * z(1) - 1.0d0 ) ** 2 ) + b
      ! if ( compgrad ) then
      !    gden(1) = a * ( -0.1d0 * z(1) * ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) + 0.2d0 * (0.25d0 * z(1) - 1.0d0) )
      !    gden(2) = a * ( 0.8d0 * ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) ) 
      ! end if

      !! Function 14
      ! z(1:2) = -pdata%c(1:2)/5.0d0 + 0.40d0 * x(1:2)
      ! scale = sqrt( 1000.0d0 )
      ! z1 = 10.0d0 + 0.25d0 * scale
      ! z2 = 2.5d0 + 0.25d0 * scale
      ! param = ( 0.1d0 * scale + 1.0d0 ) ** 4 - z2 ** 2
      ! a = (2.0d0 - 0.25d0) / param
      ! b = 2.0d0 - a * ( ( 0.1d0 * scale + 1.0d0 ) ** 4 ) 
      ! faux = ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) ** 2 + ( 0.25d0 * z(1) - 1.0d0 ) ** 2
      ! fden = a * (  faux - x(2) ) ** 2 + b

      ! if ( compgrad ) then
      !    gden(1) = 2.0d0 * a * ( faux - x(2) ) * ( -0.1d0 * z(1) * ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) + 0.2d0 * (0.25d0 * z(1) - 1.0d0) )
      !    gden(2) = 2.0d0 * a * ( faux - x(2) ) * ( 0.8d0 * ( z(2) - ( 0.25d0 * z(1) ) ** 2 ) - 1.0d0 )
      ! end if



   end subroutine fgcompJ4




   subroutine centroid(W,cent,volW)
      
      implicit none

      ! SCALAR ARGUMENTS
      type(polygon), intent(in) :: W
      real(kind=8), intent(out) :: volW

      ! ARRAY ARGUMENTS
      real(kind=8), intent(out) :: cent(2)

      ! LOCAL SCALARS
      integer :: j
      real(kind=8) :: re

      call polygon_area_2d(W%n-1,W%v(1:2,1:W%n-1),volW)
      
      cent = 0.0d0

      do j=1,W%n-1
         if ( j .eq. W%n-1 ) then
            re = W%v(1,j)*W%v(2,1) - W%v(1,1)*W%v(2,j)
            cent = cent + re*(W%v(1:2,j) + W%v(1:2,1))
         else
            re = W%v(1,j)*W%v(2,j+1) - W%v(1,j+1)*W%v(2,j)
            cent = cent + re*(W%v(1:2,j) + W%v(1:2,j+1)) 
         end if
      end do
      cent = cent/(6.0d0*volW)

   end subroutine centroid

    ! ******************************************************************
    ! ******************************************************************

   subroutine det(x,y,z,vdet)

        implicit none

        ! SACALAR ARGUMENTS
        real(kind=8), intent(out) :: vdet

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: x(2),y(2),z(2)

        vdet = ( x(1) - y(1) ) * ( z(2) - y(2) ) - ( x(2) - y(2) ) * ( z(1) - y(1) )

   end subroutine det
  
    
    ! ******************************************************************
    ! ******************************************************************

   subroutine matrixint(x,y,z,v,M,flag)

      ! It computes the matrix used in the calculation of the gradient
      ! when v belongs to T_int.
            
      implicit none

      ! SCALAR ARGUMENTS
      integer, intent(inout) :: flag

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(2),y(2),z(2),v(2)
      real(kind=8), intent(out) :: M(2,2)

      ! LOCAL SCALARS
      real(kind=8) :: det

      ! LOCAL ARRAYS
      real(kind=8) :: tmp(2)

      det = ( y(1) - x(1) ) * ( z(2) - x(2) ) - ( y(2) - x(2) ) * ( z(1) - x(1) )

      if ( abs( det ) .lt. 1.0d+02 * epsilon( 1.0d0 ) ) then
         det = sign( 1.0d+02 * epsilon( 1.0d0 ), det )
      end if

      ! This is to avoid a warning of the compiler that says flag is not
      ! being used.
      if ( .false. ) flag = - 1

      !!$    if ( det .eq. 0.0d0 ) then
      !!$       write(*,*) 'In matrixint: Division by zero in the calculation of matrix M.'
      !!$       flag = - 1
      !!$       return
      !!$       
      !!$    else if ( abs(det) .lt. 1.0d-14 ) then
      !!$       write(*,*) 'In matrixint: Warning: Determinant is too small.'
      !!$    end if
            
      tmp(1:2) = (/ y(2) - x(2), x(1) - y(1) /)

      M(1:2,1) = tmp(1:2) * ( v(1) - z(1) )
      M(1:2,2) = tmp(1:2) * ( v(2) - z(2) )
      M(1:2,1:2) = M(1:2,1:2) / det

   end subroutine matrixint

   subroutine matrixbd(y,x,v,gradphi,M,flag)
      
      ! It computes the matrix used in the calculation of the gradient
      ! when v belongs to T_bd.
          
      implicit none
          
      ! SCALAR ARGUMENTS
      integer, intent(inout) :: flag
      
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: gradphi(2),x(2),y(2),v(2)
      real(kind=8), intent(out) :: M(2,2)
          
      ! LOCAL SCALARS
      real(kind=8) :: det
  
      ! LOCAL ARRAYS
      real(kind=8) :: tmp(2)
      
      tmp(1:2) = (/ - gradphi(2), gradphi(1) /)
      
      det = - dot_product( tmp(1:2), y(1:2) - x(1:2) )
  
      if ( abs( det ) .lt. 1.0d+02 * epsilon( 1.0d0 ) ) then
         det = sign( 1.0d+02 * epsilon( 1.0d0 ), det )
      end if
  
      ! This is to avoid a warning of the compiler that says flag is not
      ! being used.
      if ( .false. ) flag = - 1
      
  !!$    if ( det .eq. 0.0d0 ) then
  !!$       write(*,*) 'In matrixbd: Division by zero in the calculation of matrix M.'
  !!$       flag = - 1
  !!$       return
  !!$       
  !!$    else if ( abs(det) .lt. 1.0d-14 ) then
  !!$       write(*,*) 'In matrixbd: Warning: Determinant is too small.'
  !!$    end if
          
      ! Outer product of (/ - gradphi(2), gradphi(1) /) with p(1:2) - x(1:2)
      M(1:2,1) = tmp(1:2) * ( v(1) - x(1) ) 
      M(1:2,2) = tmp(1:2) * ( v(2) - x(2) )
      M(1:2,1:2) = - M(1:2,1:2) / det
      
   end subroutine matrixbd

    ! ******************************************************************
    ! ******************************************************************

  
   subroutine drawsol(n,x,filename,pdataptr)
      
      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      type(c_ptr), optional, intent(in) :: pdataptr
        
      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      character(len=80) :: filename
  
      ! LOCAL SCALARS
      logical :: debug
      integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v
      type(pdata_type), pointer :: pdata
      type(VoronoiCell) :: Vi
      type(Polygon) :: Wij,Pj,frame,Viframe
      
      ! LOCAL ARRAYS
      integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
      real(kind=8) :: f(0:3,n/2),sites(2,n/2),vorxy(2,2*(n/2))
  
      debug = .false.
      call compallJ(n,x,f,debug,flag,pdataptr)
      if ( flag .ne. 0 ) then
         write(*,*) 'In drawsol, compallJ returned a non-null flag!'
         return
      end if
      ! write(*,*) 'ftmpscompallJ= ', f(3,1:n/2)
      ! write(*,*) 'J4compall = ', sum( f(3,1:n/2) ** 2 ) / (n/2) 
      call c_f_pointer(pdataptr,pdata)
  
      nsites = n / 2
      indices(1:nsites) = (/ (i,i=1,nsites) /)
      sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
    
      if ( nsites .le. 1 ) then
         write(*,*) 'In drawsol: nsites must be at least 2!'
         stop
      end if
  
      ! Compute Delaunay triangulation
  
      call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)
  
      if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
         write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
         stop
      end if
    
      call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
      
      ! Compute vertices of the Voronoi diagram
            
      if ( ntriag .gt. 2 * nsites ) then
         write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
         write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
         write(*,*) 'cannot be larger than twice the number of sites!'
         stop
      end if
    
      do k = 1,ntriag
         call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
         if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
            write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
            write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
            stop
         end if
      end do
  
      maxnv = maxval( pdata%nvpols(1:pdata%npols) )
        
      allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In drawsol: Allocation error.'
         stop
      end if
        
      open(unit=10,file=trim(filename))                           
  
      write(10,10) 1,pdata%drawscale
  
      ! Drawing Delaunay triangulation
      
      !! For the paper it is not necessary to draw the Delaunay triangulation
      do k = 1,ntriag
         do v = 1,3
            write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
         end do
      end do
  
      ! Frame
      
      frame%n = 5
      frame%deg = .false.
      
      allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In drawsol: Allocation error.'
         stop
      end if
      
      frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
      frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
      frame%e(1:4) = (/ 0, 0, 0, 0 /)  
          
      ! Drawing Voronoi diagram
      
      do i = 1,nsites
         call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
         if ( status .ne. 0 ) then
            write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
            stop
         end if
  
         ! This is to draw the Voronoi cell restricted to the frame.
         
         if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
         else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
            call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
         end if
  
         if ( .not. Viframe%deg ) then
            do v = 1,Viframe%n - 1
               write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
            end do
         end if
  
         call polygon_destroy(Viframe)
  
         ! This is to draw the Voronoi cell restricted to A.
  
         do j = 1,pdata%npols
            Pj%n = pdata%nvpols(j)
            Pj%deg = .false.
  
            Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
            Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
            Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)
  
            if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
               call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
            else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
               call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
            end if
  
            if ( .not. Wij%deg ) then
               do v = 1,Wij%n - 1
                  if ( Wij%e(v) .lt. 1 ) then
                     write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
                  end if
               end do
            end if
  
            call polygon_destroy(Wij)
         end do
  
         ! This is to draw the sites
         write(10,40) sites(1:2,i)
         
         call voronoi_cell_destroy(Vi)
      end do
        
      deallocate(frame%v,frame%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In drawsol: Deallocation error.'
         stop
      end if
          
      ! Drawing A
      
      do j = 1,pdata%npols
         Pj%n = pdata%nvpols(j)
         Pj%deg = .false.          
  
         Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
         Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
         Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)
  
         do k = 1,Pj%n - 1
            if ( Pj%e(k) .eq. 0 ) then
               write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
            end if
         end do
      end do
  
      ! Close figure
      
      write(10,50)
      
      write(10,*) 'Value of J4i (sum/ave/min/max) = ',sum(abs(f(3,1:n/2))),sum(abs(f(3,1:n/2)))/(n/2), &
           minval(abs(f(3,1:n/2))),maxval(abs(f(3,1:n/2))),' J4 itself = ',sum(f(3,1:n/2)**2)/(n/2)
  
      close(10)
  
      deallocate(Pj%v,Pj%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In drawsol: Deallocation error.'
         stop
      end if
  
      ! NON-EXECUTABLE STATEMENTS
      
  10  format('prologues := 3;',/, &
           'outputtemplate := "%j.mps";',/, &
           'input mpcolornames;',/, &
           'beginfig(',i4,');',/, &
           'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
           'orange, plum, bluea, blueb, bluec, blued, bluee, corn, deepsaffron, ferrarired;',/, &
           'skyblue      = (135/255,206/255,235/255);',/, &
           'dollarbill   = (133/255,187/255,101/255);',/, &
           'teagreen     = (208/255,240/255,192/255);',/, &
           'desire       = (234/255, 60/255, 83/255);',/, &
           'jonquil      = (244/255,202/255, 22/255);',/, &
           'royalblue    = ( 65/255,105/255,225/255);',/, &
           'orange       = (255/255,165/255,  0/255);',/, &
           'plum         = (221/255,160/255,221/255);',/, &
           'bluea        = (  0/255,230/255,230/255);',/, &
           'blueb        = (  0/255,191/255,230/255);',/, &
           'bluec        = (  0/255,153/255,230/255);',/, &
           'blued        = (  0/255,115/255,230/255);',/, &
           'bluee        = (  0/255, 77/255,230/255);',/, &
           'corn         = (255/255,230/255,102/255);',/, &
           'deepsaffron  = (255/255,153/255, 51/255);',/, &
           'ferrarired   = (255/255, 42/255,  0/255);',/, &
           'u = ',f20.10,'cm;')
  20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
  40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
  50  format('endfig;',/,'end;')
  
end subroutine drawsol


! ******************************************************************
! ******************************************************************

! ! Draw the unbalanced edges with blue. (Here I think I have an error)
! subroutine drawsol1(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext
!    real(kind=8) :: avgsizeed, perim
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vsizeed(10)

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols
!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
!             perim = 0.0d0
!             do v = 1,Wij%n - 1
!                if ( v .eq. Wij%n-1) then
!                   vnext = 1
!                else 
!                   vnext = v + 1
!                end if 
!                tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                vsizeed(v) = norm2( tauv(1:2) )
!                perim = perim + vsizeed(v)
!             end do
!             avgsizeed = perim / (Wij%n-1)

!             do v = 1,Wij%n - 1

!                if ( Wij%e(v) .lt. 1 ) then
!                   if ( vsizeed(v)/avgsizeed .lt. pdata%J2ctol ) then
!                      write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1), 'bluep'
!                   else
!                      write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                   end if
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do

!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orange, plum, bluea, blueb, bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orange       = (255/255,165/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol1

! ! Draw the edges that form an unbalanced angle with blue.
! subroutine drawsol2(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext,vprev
!    real(kind=8) :: avgtau, sumtau
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vang(10),tauvprev(2)

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols
!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then

!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do


!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do

!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if
!       do j = 1,pdata%npols
!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
!             sumtau = 0.0d0
!             vang(1:10) = 0.0d0
!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   if ( v .eq. Wij%n-1) then
!                      vnext = 1
!                      vprev = v - 1
!                   else if ( v .eq. 1 ) then
!                      vnext = v + 1
!                      vprev = Wij%n - 1
!                   else
!                      vnext = v + 1
!                      vprev = v - 1
!                   end if 
!                   tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                   tauvprev(1:2) = Wij%v(1:2,v) - Wij%v(1:2,vprev)
!                   vang(v) = acos( max( -1.0d0, min( dot_product( -tauvprev(1:2)/norm2( tauvprev(1:2) ), tauv(1:2)/norm2( tauv(1:2) ) ), 1.0d0) ) )
!                   sumtau = sumtau + vang(v)
!                end if
!             end do
!             avgtau =  sumtau / count( vang(1:10) .ne. 0.0d0)


!             do v = 1,Wij%n - 1
!                   if ( v .eq. Wij%n-1) then
!                      vnext = 1
!                      vprev = v - 1
!                   else if ( v .eq. 1 ) then
!                      vnext = v + 1
!                      vprev = Wij%n - 1
!                   else
!                      vnext = v + 1
!                      vprev = v - 1
!                   end if 
!                   if ( vang(v)/avgtau .lt. pdata%J3ctol ) then
!                      write(10,20) Wij%v(1:2,vprev),Wij%v(1:2,v), 'royalblue'
!                      write(10,20) Wij%v(1:2,v),Wij%v(1:2,vnext), 'royalblue'
!                   end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orange, plum, bluea, blueb, bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orange       = (255/255,165/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol2

! ! Fill the cells with unbalanced edges with blue. 
! subroutine drawsol3(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext
!    real(kind=8) :: avgsizeed, perim, volW
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vsizeed(20)

!    ! LOCAL PARAMETER
!    real(kind=8), parameter :: tol = 5.0d-2

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols
!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
            
!             call polygon_area_2d(Wij%n-1,Wij%v(1:2,1:Wij%n-1),volW)
!             write(*,*) 'Vol',i,'=',volW
!             ! write(*,*) 'Vol - 1.0d0=', volW-1.0d0
!             if ( abs( volW - 1.0d0 ) .gt. tol ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'ferrarired'
!                   end if
!                end do
!             end if

!             if ( pdata%Jweight(1) .gt. 0.0d0 ) then
!                perim = 0.0d0
!                do v = 1,Wij%n - 1
!                   if ( v .eq. Wij%n-1) then
!                      vnext = 1
!                   else 
!                      vnext = v + 1
!                   end if 
!                   tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                   vsizeed(v) = norm2( tauv(1:2) )
!                   perim = perim + vsizeed(v)
!                end do
!                avgsizeed = perim / (Wij%n-1)

!                if ( any( vsizeed(1:Wij%n-1)/avgsizeed .lt. pdata%J2ctol ) ) then
!                   do v = 1,Wij%n-1
!                      if ( v .eq. 1 ) then
!                         write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                      else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                         write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                      else
!                         write(10,32) Wij%v(1:2,v), 'bluea'
!                      end if
!                   end do
!                end if
!             end if 


!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do

!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orange, plum, bluea, blueb, bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orange       = (255/255,165/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol3



! subroutine drawsol4(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext
!    real(kind=8) :: avgsizeed, perim
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vsizeed(20)

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols
!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
!             perim = 0.0d0
!             do v = 1,Wij%n - 1
!                if ( v .eq. Wij%n-1) then
!                   vnext = 1
!                else 
!                   vnext = v + 1
!                end if 
!                tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                vsizeed(v) = norm2( tauv(1:2) )
!                perim = perim + vsizeed(v)
!             end do
!             avgsizeed = perim / (Wij%n-1)

!             if ( any( vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.1 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluee'
!                   end if
!                end do
!             else if ( any( 0.1 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.2 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blued'
!                   end if
!                end do
!             else if ( any(0.2 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.3  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluec'
!                   end if
!                end do
!             else if ( any(0.3 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.4  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blueb'
!                   end if
!                end do
!             else if ( any(0.4 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.5  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluea'
!                   end if
!                end do   
!             end if


!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do

!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orange, plum, bluea, blueb, bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orange       = (255/255,165/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol4


! subroutine drawsol5(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext
!    real(kind=8) :: avgsizeed, perim
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vsizeed(20)

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols
!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
!             perim = 0.0d0
!             do v = 1,Wij%n - 1
!                if ( v .eq. Wij%n-1) then
!                   vnext = 1
!                else 
!                   vnext = v + 1
!                end if 
!                tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                vsizeed(v) = norm2( tauv(1:2) )
!                perim = perim + vsizeed(v)
!             end do
!             avgsizeed = perim / (Wij%n-1)

!             if ( any( 0.4d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed  .lt. 0.41d0 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluedd'
!                   end if
!                end do
!             else if ( any( 0.41d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.42d0 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluee'
!                   end if
!                end do
!             else if ( any( 0.42d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.43d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blued'
!                   end if
!                end do
!             else if ( any( 0.43d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.44d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blued'
!                   end if
!                end do
!             else if ( any( 0.44d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.45d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blueb'
!                   end if
!                end do
!             else if ( any( 0.45d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.46d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluea'
!                   end if
!                end do
!             else if ( any( 0.46d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.47d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'corn'
!                   end if
!                end do
!             else if ( any( 0.47d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.48d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orange'
!                   end if
!                end do
!             else if ( any( 0.48d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.49d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'plum'
!                   end if
!                end do
!             else if ( any( 0.49d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.50d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'deepsaffron'
!                   end if
!                end do
!             end if


!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do

!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orange, plum, bluea, blueb, bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orange       = (255/255,165/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol5

! subroutine drawsol6(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext
!    real(kind=8) :: avgsizeed, perim
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vsizeed(20)

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols
!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
!             perim = 0.0d0
!             do v = 1,Wij%n - 1
!                if ( v .eq. Wij%n-1) then
!                   vnext = 1
!                else 
!                   vnext = v + 1
!                end if 
!                tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                vsizeed(v) = norm2( tauv(1:2) )
!                perim = perim + vsizeed(v)
!             end do
!             avgsizeed = perim / (Wij%n-1)

!             if ( any( 0.3d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed  .lt. 0.31d0 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluedd'
!                   end if
!                end do
!             else if ( any( 0.31d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.32d0 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluee'
!                   end if
!                end do
!             else if ( any( 0.32d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.33d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blued'
!                   end if
!                end do
!             else if ( any( 0.33d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.34d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blued'
!                   end if
!                end do
!             else if ( any( 0.34d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.35d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blueb'
!                   end if
!                end do
!             else if ( any( 0.35d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.36d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluea'
!                   end if
!                end do
!             else if ( any( 0.36d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.37d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'corn'
!                   end if
!                end do
!             else if ( any( 0.37d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.38d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orange'
!                   end if
!                end do
!             else if ( any( 0.38d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.39d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'plum'
!                   end if
!                end do
!             else if ( any( 0.39d0 .le. vsizeed(1:Wij%n-1)/avgsizeed .and. vsizeed(1:Wij%n-1)/avgsizeed .lt. 0.40d0  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'deepsaffron'
!                   end if
!                end do
!             end if


!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do

!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orange, plum, bluea, blueb, bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orange       = (255/255,165/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol6

! ! Fill the cells with unbalanced angles with different colors the darker the more unbalanced the angles.
! subroutine drawsol7(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext,vprev,aux,counter
!    real(kind=8) :: sumtau,avgtau
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vang(20),tauvprev(2)

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols

!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
!             sumtau = 0.0d0
!             vang(1:20) = 0.0d0
!             aux = 0
!             do v = 1,Wij%n - 1
!                   if ( v .eq. Wij%n-1) then
!                      vnext = 1
!                      vprev = v - 1
!                   else if ( v .eq. 1 ) then
!                      vnext = v + 1
!                      vprev = Wij%n - 1
!                   else
!                      vnext = v + 1
!                      vprev = v - 1
!                   end if 

!                   if ( Wij%e(v) .lt. 0 .or. Wij%e(vprev) .lt. 0 ) then
!                      tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                      tauvprev(1:2) = Wij%v(1:2,v) - Wij%v(1:2,vprev)
!                      vang(v-aux) = acos( max( -1.0d0, min( dot_product( -tauvprev(1:2)/norm2( tauvprev(1:2) ), tauv(1:2)/norm2( tauv(1:2) ) ), 1.0d0) ) )
!                      sumtau = sumtau + vang(v-aux)
!                   else
!                      aux = 1
!                   end if

!             end do
!             counter = count( vang(1:20) .ne. 0.0d0 )
!             !write(*,*)'vang = ', vang(1:20)
!             avgtau =  sumtau / counter
!             !write(*,*)'count = ', count( vang(1:20) .ne. 0.0d0 )

!             if ( any( vang(1:counter)/avgtau .lt. 0.4 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangee'
!                   end if
!                end do
!             else if ( any( 0.4 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.5 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'oranged'
!                   end if
!                end do
!             else if ( any(0.5 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.6  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangec'
!                   end if
!                end do
!             else if ( any(0.6 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.7  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangeb'
!                   end if
!                end do
!             else if ( any(0.7 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.8  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangea'
!                   end if
!                end do   
!             end if


!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do

!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orangea, orangeb, orangec, oranged, orangee, plum, bluea, blueb, ', &
!         'bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orangea      = (255/255,165/255,  0/255);',/, &
!         'orangeb      = (223/255,139/255,  0/255);',/, &
!         'orangec      = (192/255,113/255,  0/255);',/, &
!         'oranged      = (161/255, 88/255,  0/255);',/, &
!         'orangee      = (131/255, 65/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol7

! ! Fill the cells with unbalanced angles with different colors between 0.7 and 0.8 in intervals of 0.01
! subroutine drawsol8(n,x,filename,pdataptr)
      
!    implicit none
 
!    ! SCALAR ARGUMENTS
!    integer, intent(in) :: n
!    type(c_ptr), optional, intent(in) :: pdataptr
     
!    ! ARRAY ARGUMENTS
!    real(kind=8), intent(in) :: x(n)
!    character(len=80) :: filename

!    ! LOCAL SCALARS
!    logical :: debug
!    integer :: allocerr,flag,i,ierror,j,k,maxnv,nsites,status,ntriag,v,vnext,vprev,aux,counter
!    real(kind=8) :: sumtau,avgtau
!    type(pdata_type), pointer :: pdata
!    type(VoronoiCell) :: Vi
!    type(Polygon) :: Wij,Pj,frame,Viframe
   
!    ! LOCAL ARRAYS
!    integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
!    real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2)),tauv(2),vang(20),tauvprev(2)

!    debug = .false.
!    !call compallJ(n,x,f,debug,flag,pdataptr)
!    if ( flag .ne. 0 ) then
!       write(*,*) 'In drawsol, compallJ returned a non-null flag!'
!       return
!    end if
   
!    call c_f_pointer(pdataptr,pdata)

!    nsites = n / 2
!    indices(1:nsites) = (/ (i,i=1,nsites) /)
!    sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
!    if ( nsites .le. 1 ) then
!       write(*,*) 'In drawsol: nsites must be at least 2!'
!       stop
!    end if

!    ! Compute Delaunay triangulation

!    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

!    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
!       write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
!       stop
!    end if
 
!    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
!    ! Compute vertices of the Voronoi diagram
         
!    if ( ntriag .gt. 2 * nsites ) then
!       write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
!       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
!       write(*,*) 'cannot be larger than twice the number of sites!'
!       stop
!    end if
 
!    do k = 1,ntriag
!       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
!       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
!          write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
!          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
!          stop
!       end if
!    end do

!    maxnv = maxval( pdata%nvpols(1:pdata%npols) )
     
!    allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
     
!    open(unit=10,file=trim(filename))                           

!    write(10,10) 1,pdata%drawscale

!    ! Drawing Delaunay triangulation
   
!    do k = 1,ntriag
!       do v = 1,3
!          write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
!       end do
!    end do

!    ! Frame
   
!    frame%n = 5
!    frame%deg = .false.
   
!    allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Allocation error.'
!       stop
!    end if
   
!    frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
!    frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
!    frame%e(1:4) = (/ 0, 0, 0, 0 /)  
       
!    ! Drawing Voronoi diagram
   
!    do i = 1,nsites
!       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
!       if ( status .ne. 0 ) then
!          write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
!          stop
!       end if

!       ! This is to draw the Voronoi cell restricted to the frame.
      
!       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!          call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
!       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!          call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
!       end if

!       if ( .not. Viframe%deg ) then
!          do v = 1,Viframe%n - 1
!             write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
!          end do
!       end if

!       call polygon_destroy(Viframe)

!       ! This is to draw the Voronoi cell restricted to A.

!       do j = 1,pdata%npols

!          Pj%n = pdata%nvpols(j)
!          Pj%deg = .false.

!          Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!          Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!          Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!          if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
!             call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
!          else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
!             call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
!          end if



!          if ( .not. Wij%deg ) then
!             sumtau = 0.0d0
!             vang(1:20) = 0.0d0
!             aux = 0
!             do v = 1,Wij%n - 1
!                   if ( v .eq. Wij%n-1) then
!                      vnext = 1
!                      vprev = v - 1
!                   else if ( v .eq. 1 ) then
!                      vnext = v + 1
!                      vprev = Wij%n - 1
!                   else
!                      vnext = v + 1
!                      vprev = v - 1
!                   end if 

!                   if ( Wij%e(v) .lt. 0 .or. Wij%e(vprev) .lt. 0 ) then
!                      tauv(1:2) = Wij%v(1:2,vnext) - Wij%v(1:2,v)
!                      tauvprev(1:2) = Wij%v(1:2,v) - Wij%v(1:2,vprev)
!                      vang(v-aux) = acos( max( -1.0d0, min( dot_product( -tauvprev(1:2)/norm2( tauvprev(1:2) ), tauv(1:2)/norm2( tauv(1:2) ) ), 1.0d0) ) )
!                      sumtau = sumtau + vang(v-aux)
!                   else
!                      aux = 1
!                   end if

!             end do
!             counter = count( vang(1:20) .ne. 0.0d0 )
!             !write(*,*)'vang = ', vang(1:20)
!             avgtau =  sumtau / counter
!             !write(*,*)'count = ', count( vang(1:20) .ne. 0.0d0 )

!             if ( any( 0.7 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau  .lt. 0.71 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangee'
!                   end if
!                end do
!             else if ( any( 0.71 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.72 ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'oranged'
!                   end if
!                end do
!             else if ( any( 0.72 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.73  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangec'
!                   end if
!                end do
!             else if ( any( 0.74 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.75  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangeb'
!                   end if
!                end do
!             else if ( any( 0.75 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.76  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'orangea'
!                   end if
!                end do
!             else if  ( any( 0.76 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.77  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluea'
!                   end if
!                end do
!             else if  ( any( 0.77 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.78  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blueb'
!                   end if
!                end do
!             else if  ( any( 0.78 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.79  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'bluec'
!                   end if
!                end do
!             else if  ( any( 0.79 .le. vang(1:counter)/avgtau .and. vang(1:counter)/avgtau .lt. 0.80  ) ) then
!                do v = 1,Wij%n-1
!                   if ( v .eq. 1 ) then
!                      write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else if ( v .gt. 1 .and. v .lt. Wij%n-1 ) then
!                      write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
!                   else
!                      write(10,32) Wij%v(1:2,v), 'blued'
!                   end if
!                end do
!             end if


!             do v = 1,Wij%n - 1
!                if ( Wij%e(v) .lt. 1 ) then
!                   write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
!                end if
!             end do
!          end if

!          call polygon_destroy(Wij)
!       end do

!       ! This is to draw the sites
!       write(10,40) sites(1:2,i)
      
!       call voronoi_cell_destroy(Vi)
!    end do
     
!    deallocate(frame%v,frame%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if
       
!    ! Drawing A
   
!    do j = 1,pdata%npols
!       Pj%n = pdata%nvpols(j)
!       Pj%deg = .false.          

!       Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
!       Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
!       Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

!       do k = 1,Pj%n - 1
!          if ( Pj%e(k) .eq. 0 ) then
!             write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
!          end if
!       end do
!    end do

!    ! Close figure
   
!    write(10,50)
   
!    write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
!         minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)

!    close(10)

!    deallocate(Pj%v,Pj%e,stat=allocerr)
!    if ( allocerr .ne. 0 ) then
!       write(*,*) 'In drawsol: Deallocation error.'
!       stop
!    end if

!    ! NON-EXECUTABLE STATEMENTS
   
! 10  format('prologues := 3;',/, &
!         'outputtemplate := "%j.mps";',/, &
!         'input mpcolornames;',/, &
!         'beginfig(',i4,');',/, &
!         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
!         'orangea, orangeb, orangec, oranged, orangee, plum, bluea, blueb, ', &
!         'bluec, blued, bluee, bluedd,bluep , corn, deepsaffron, ferrarired;',/, &
!         'skyblue      = (135/255,206/255,235/255);',/, &
!         'dollarbill   = (133/255,187/255,101/255);',/, &
!         'teagreen     = (208/255,240/255,192/255);',/, &
!         'desire       = (234/255, 60/255, 83/255);',/, &
!         'jonquil      = (244/255,202/255, 22/255);',/, &
!         'royalblue    = ( 65/255,105/255,225/255);',/, &
!         'orangea      = (255/255,165/255,  0/255);',/, &
!         'orangeb      = (223/255,139/255,  0/255);',/, &
!         'orangec      = (192/255,113/255,  0/255);',/, &
!         'oranged      = (161/255, 88/255,  0/255);',/, &
!         'orangee      = (131/255, 65/255,  0/255);',/, &
!         'plum         = (221/255,160/255,221/255);',/, &
!         'bluea        = (  0/255,230/255,230/255);',/, &
!         'blueb        = (  0/255,191/255,230/255);',/, &
!         'bluec        = (  0/255,153/255,230/255);',/, &
!         'blued        = (  0/255,115/255,230/255);',/, &
!         'bluee        = (  0/255, 77/255,230/255);',/, &
!         'bluedd       = (  0/255, 0/255, 100/255);',/, & 
!         'bluep        = (  0/255,100/255,255/255);',/, & 
!         'corn         = (255/255,230/255,102/255);',/, &
!         'deepsaffron  = (255/255,153/255, 51/255);',/, &
!         'ferrarired   = (255/255, 42/255,  0/255);',/, &
!         'u = ',f20.10,'cm;')
! 20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
! 30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
! 32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')
! 40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
! 50  format('endfig;',/,'end;')

! end subroutine drawsol8

!! This drawsol are the definitive for J2 I hope.

subroutine drawsolJ2(n,x,filename,pdataptr)

   implicit none
    
   ! SCALAR ARGUMENTS
   integer, intent(in) :: n
   type(c_ptr), optional, intent(in) :: pdataptr
     
   ! ARRAY ARGUMENTS
   real(kind=8), intent(in) :: x(n)
   character(len=80) :: filename

   ! LOCAL SCALARS
   integer :: allocerr,i,ierror,j,k,maxnv,nsites,status,ntriag,v
   type(pdata_type), pointer :: pdata
   type(VoronoiCell) :: Vi
   type(Polygon) :: Wij,Pj,frame,Viframe
   
   ! LOCAL ARRAYS
   integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
   real(kind=8) :: sites(2,n/2),vorxy(2,2*(n/2))

   integer :: flag,find
   logical :: debug
   real(kind=8) :: ctol, J2ctoltmp
   integer :: cellcolorJ2(n/2)
   real(kind=8) :: f(0:3,n/2)
   character(len=80) :: paintcolor

   call c_f_pointer(pdataptr,pdata) 

   J2ctoltmp = pdata%J2ctol
   
   ! The smaller the ctol for which the cell has positive f, the darker the color.
   cellcolorJ2(1:n/2) = -1       ! Means no color

   do k = 5,1,-1
      ctol = 0.1d0 * k
      pdata%J2ctol = ctol
      debug = .false.
      call compallJ(n,x,f,debug,flag,pdataptr)
      if ( flag .ne. 0 ) then
         write(*,*) 'In drawsol, compallJ returned a non-null flag!'
         return
      end if

      where ( f(1,1:n/2) .gt. 1.0d-04 )
         cellcolorJ2(1:n/2) = k
      end where
   end do

   pdata%J2ctol = J2ctoltmp

   call c_f_pointer(pdataptr,pdata)
  
   nsites = n / 2
   indices(1:nsites) = (/ (i,i=1,nsites) /)
   sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
   if ( nsites .le. 1 ) then
      write(*,*) 'In drawsol: nsites must be at least 2!'
      stop
   end if

   ! Compute Delaunay triangulation

   call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

   if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
      write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
      stop
   end if
 
   call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
   ! Compute vertices of the Voronoi diagram
         
   if ( ntriag .gt. 2 * nsites ) then
      write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
      write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
      write(*,*) 'cannot be larger than twice the number of sites!'
      stop
   end if
 
   do k = 1,ntriag
      call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
      if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
         write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
         write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
         stop
      end if
   end do

   maxnv = maxval( pdata%nvpols(1:pdata%npols) )

   allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Allocation error.'
      stop
   end if
     
   open(unit=10,file=trim(filename))                           

   write(10,10) 1,pdata%drawscale

   ! Drawing Delaunay triangulation
   
   do k = 1,ntriag
      do v = 1,3
         write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
      end do
   end do

   ! Frame
   
   frame%n = 5
   frame%deg = .false.
   
   allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Allocation error.'
      stop
   end if

   frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
   frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
   frame%e(1:4) = (/ 0, 0, 0, 0 /)  
          
      ! Drawing Voronoi diagram
      
      do i = 1,nsites
         call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
         if ( status .ne. 0 ) then
            write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
            stop
         end if
  
         ! This is to draw the Voronoi cell restricted to the frame.
         
         if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
         else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
            call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
         end if
  
         if ( .not. Viframe%deg ) then
            do v = 1,Viframe%n - 1
               write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
            end do
         end if
  
         call polygon_destroy(Viframe)

          ! This is to paint the Voronoi cell, restricted to A, if
         ! desired conditions does not hold.
  
         if ( cellcolorJ2(i) .ne. - 1 ) then
            if ( cellcolorJ2(i) .eq. 5 ) then
               paintcolor = 'bluea'
            else if ( cellcolorJ2(i) .eq. 4 ) then
               paintcolor = 'blueb'
            else if ( cellcolorJ2(i) .eq. 3 ) then
               paintcolor = 'bluec'
            else if ( cellcolorJ2(i) .eq. 2 ) then
               paintcolor = 'blued'
            else if ( cellcolorJ2(i) .eq. 1 ) then
               paintcolor = 'bluee'
            end if

            do j = 1,pdata%npols
               Pj%n = pdata%nvpols(j)
               Pj%deg = .false.         
               
               Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
               Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
               Pj%e(1:Pj%n-1)  = pdata%poledges(j,1:Pj%n-1)
  
               if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
                  call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
               else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
                  call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
               end if
  
               if ( .not. Wij%deg ) then
                  do v = 1,Wij%n - 1
                     if ( v .eq. 1 ) then
                        write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
                     else if ( v .gt. 1 .and. v .lt. Wij%n - 1 ) then 
                        write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
                     else
                        write(10,32) Wij%v(1:2,v),trim(paintcolor)
                     end if
                  end do
               end if
  
               call polygon_destroy(Wij)
            end do
         end if

         ! This is to draw the Voronoi cell restricted to A.
  
         do j = 1,pdata%npols
            Pj%n = pdata%nvpols(j)
            Pj%deg = .false.
  
            Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
            Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
            Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)
  
            if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
               call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
            else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
               call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
            end if
  
            if ( .not. Wij%deg ) then
               do v = 1,Wij%n - 1
                  if ( Wij%e(v) .lt. 1 ) then
                     write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
                  end if
               end do
            end if
  
            call polygon_destroy(Wij)
         end do



      ! This is to draw the sites
      write(10,40) sites(1:2,i)
      
      call voronoi_cell_destroy(Vi)
   end do

   deallocate(frame%v,frame%e,stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Deallocation error.'
      stop
   end if
       
   ! Drawing A
   
   do j = 1,pdata%npols
      Pj%n = pdata%nvpols(j)
      Pj%deg = .false.          

      Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
      Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
      Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

      do k = 1,Pj%n - 1
         if ( Pj%e(k) .eq. 0 ) then
            write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
         end if
      end do
   end do

   ! Close figure

   write(10,50)

   close(10)

   deallocate(Pj%v,Pj%e,stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Deallocation error.'
      stop
   end if

   ! NON-EXECUTABLE STATEMENTS
   
10  format('prologues := 3;',/, &
        'outputtemplate := "%j.mps";',/, &
        'input mpcolornames;',/, &
        'beginfig(',i4,');',/, &
        'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
        'orange, plum, bluea, blueb, bluec, blued, bluee, corn, deepsaffron, ferrarired;',/, &
        'skyblue      = (135/255,206/255,235/255);',/, &
        'dollarbill   = (133/255,187/255,101/255);',/, &
        'teagreen     = (208/255,240/255,192/255);',/, &
        'desire       = (234/255, 60/255, 83/255);',/, &
        'jonquil      = (244/255,202/255, 22/255);',/, &
        'royalblue    = ( 65/255,105/255,225/255);',/, &
        'orange       = (255/255,165/255,  0/255);',/, &
        'plum         = (221/255,160/255,221/255);',/, &
        'bluea        = (  0/255,230/255,230/255);',/, &
        'blueb        = (  0/255,191/255,230/255);',/, &
        'bluec        = (  0/255,153/255,230/255);',/, &
        'blued        = (  0/255,115/255,230/255);',/, &
        'bluee        = (  0/255, 77/255,230/255);',/, &
        'corn         = (255/255,230/255,102/255);',/, &
        'deepsaffron  = (255/255,153/255, 51/255);',/, &
        'ferrarired   = (255/255, 42/255,  0/255);',/, &
        'u = ',f20.10,'cm;')
20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')              
40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
50  format('endfig;',/,'end;')

end subroutine drawsolJ2

subroutine drawsolJ3(n,x,filename,pdataptr)

   implicit none
    
   ! SCALAR ARGUMENTS
   integer, intent(in) :: n
   type(c_ptr), optional, intent(in) :: pdataptr
     
   ! ARRAY ARGUMENTS
   real(kind=8), intent(in) :: x(n)
   character(len=80) :: filename

   ! LOCAL SCALARS
   integer :: allocerr,i,ierror,j,k,maxnv,nsites,status,ntriag,v
   type(pdata_type), pointer :: pdata
   type(VoronoiCell) :: Vi
   type(Polygon) :: Wij,Pj,frame,Viframe
   
   ! LOCAL ARRAYS
   integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
   real(kind=8) :: sites(2,n/2),vorxy(2,2*(n/2))

   integer :: flag,find
   logical :: debug
   real(kind=8) :: ctol, J3ctoltmp
   integer :: cellcolorJ3(n/2)
   real(kind=8) :: f(0:3,n/2)
   character(len=80) :: paintcolor

   call c_f_pointer(pdataptr,pdata) 

   J3ctoltmp = pdata%J3ctol

   ! The smaller the ctol for which the cell has positive f, the darker the color.
   cellcolorJ3(1:n/2) = -1       ! Means no color

   do k = 8,5,-1
      ctol = 0.1d0 * k
      pdata%J3ctol = ctol
      
      debug = .false.
      call compallJ(n,x,f,debug,flag,pdataptr)
      if ( flag .ne. 0 ) then
         write(*,*) 'In drawsol, compallJ returned a non-null flag!'
         return
      end if

      where ( f(2,1:n/2) .gt. 1.0d-04 )
         cellcolorJ3(1:n/2) = k
      end where
   end do

   pdata%J3ctol = J3ctoltmp

   call c_f_pointer(pdataptr,pdata)
  
   nsites = n / 2
   indices(1:nsites) = (/ (i,i=1,nsites) /)
   sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
   if ( nsites .le. 1 ) then
      write(*,*) 'In drawsol: nsites must be at least 2!'
      stop
   end if

   ! Compute Delaunay triangulation

   call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

   if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
      write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
      stop
   end if
 
   call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
   ! Compute vertices of the Voronoi diagram
         
   if ( ntriag .gt. 2 * nsites ) then
      write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
      write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
      write(*,*) 'cannot be larger than twice the number of sites!'
      stop
   end if
 
   do k = 1,ntriag
      call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
      if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
         write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
         write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
         stop
      end if
   end do

   maxnv = maxval( pdata%nvpols(1:pdata%npols) )

   allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Allocation error.'
      stop
   end if
     
   open(unit=10,file=trim(filename))                           

   write(10,10) 1,pdata%drawscale

   ! Drawing Delaunay triangulation
   
   do k = 1,ntriag
      do v = 1,3
         write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
      end do
   end do

   ! Frame
   
   frame%n = 5
   frame%deg = .false.
   
   allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Allocation error.'
      stop
   end if

   frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
      frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
      frame%e(1:4) = (/ 0, 0, 0, 0 /)  
          
      ! Drawing Voronoi diagram
      
      do i = 1,nsites
         call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
         if ( status .ne. 0 ) then
            write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
            stop
         end if
  
         ! This is to draw the Voronoi cell restricted to the frame.
         
         if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
         else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
            call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
         end if
  
         if ( .not. Viframe%deg ) then
            do v = 1,Viframe%n - 1
               write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
            end do
         end if
  
         call polygon_destroy(Viframe)

          ! This is to paint the Voronoi cell, restricted to A, if
         ! desired conditions does not hold.
  
         if ( cellcolorJ3(i) .ne. - 1 ) then
            if ( cellcolorJ3(i) .eq. 8 ) then
               paintcolor = 'corn'
            else if ( cellcolorJ3(i) .eq. 7 ) then
               paintcolor = 'deepsaffron'
            else if ( cellcolorJ3(i) .eq. 6 ) then
               paintcolor = 'ferrarired'
            else if ( cellcolorJ3(i) .eq. 5 ) then
               paintcolor = 'darkerred'
            end if

            do j = 1,pdata%npols
               Pj%n = pdata%nvpols(j)
               Pj%deg = .false.         
               
               Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
               Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
               Pj%e(1:Pj%n-1)  = pdata%poledges(j,1:Pj%n-1)
  
               if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
                  call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
               else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
                  call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
               end if
  
               if ( .not. Wij%deg ) then
                  do v = 1,Wij%n - 1
                     if ( v .eq. 1 ) then
                        write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
                     else if ( v .gt. 1 .and. v .lt. Wij%n - 1 ) then 
                        write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
                     else
                        write(10,32) Wij%v(1:2,v),trim(paintcolor)
                     end if
                  end do
               end if
  
               call polygon_destroy(Wij)
            end do
         end if

         ! This is to draw the Voronoi cell restricted to A.
  
         do j = 1,pdata%npols
            Pj%n = pdata%nvpols(j)
            Pj%deg = .false.
  
            Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
            Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
            Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)
  
            if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
               call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
            else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
               call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
            end if
  
            if ( .not. Wij%deg ) then
               do v = 1,Wij%n - 1
                  if ( Wij%e(v) .lt. 1 ) then
                     write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
                  end if
               end do
            end if
  
            call polygon_destroy(Wij)
         end do



      ! This is to draw the sites
      write(10,40) sites(1:2,i)
      
      call voronoi_cell_destroy(Vi)
   end do

   deallocate(frame%v,frame%e,stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Deallocation error.'
      stop
   end if
       
   ! Drawing A
   
   do j = 1,pdata%npols
      Pj%n = pdata%nvpols(j)
      Pj%deg = .false.          

      Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
      Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
      Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

      do k = 1,Pj%n - 1
         if ( Pj%e(k) .eq. 0 ) then
            write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
         end if
      end do
   end do

   ! Close figure

   write(10,50)

   close(10)

   deallocate(Pj%v,Pj%e,stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Deallocation error.'
      stop
   end if

   ! NON-EXECUTABLE STATEMENTS
   
10  format('prologues := 3;',/, &
        'outputtemplate := "%j.mps";',/, &
        'input mpcolornames;',/, &
        'beginfig(',i4,');',/, &
        'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
        'orange, plum, bluea, blueb, bluec, blued, bluee, corn, deepsaffron, ferrarired,darkerred;',/, &
        'skyblue      = (135/255,206/255,235/255);',/, &
        'dollarbill   = (133/255,187/255,101/255);',/, &
        'teagreen     = (208/255,240/255,192/255);',/, &
        'desire       = (234/255, 60/255, 83/255);',/, &
        'jonquil      = (244/255,202/255, 22/255);',/, &
        'royalblue    = ( 65/255,105/255,225/255);',/, &
        'orange       = (255/255,165/255,  0/255);',/, &
        'plum         = (221/255,160/255,221/255);',/, &
        'bluea        = (  0/255,230/255,230/255);',/, &
        'blueb        = (  0/255,191/255,230/255);',/, &
        'bluec        = (  0/255,153/255,230/255);',/, &
        'blued        = (  0/255,115/255,230/255);',/, &
        'bluee        = (  0/255, 77/255,230/255);',/, &
        'corn         = (255/255,230/255,102/255);',/, &
        'deepsaffron  = (255/255,153/255, 51/255);',/, &
        'ferrarired   = (255/255, 42/255,  0/255);',/, &
        'darkerred    = (122/255, 20/255,  0/255);',/, &
        'u = ',f20.10,'cm;')
20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')              
40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
50  format('endfig;',/,'end;')

end subroutine drawsolJ3

subroutine drawsolJ1(n,x,filename,pdataptr)

   implicit none
    
   ! SCALAR ARGUMENTS
   integer, intent(in) :: n
   type(c_ptr), optional, intent(in) :: pdataptr
     
   ! ARRAY ARGUMENTS
   real(kind=8), intent(in) :: x(n)
   character(len=80) :: filename

   ! LOCAL SCALARS
   integer :: allocerr,i,ierror,j,k,maxnv,nsites,status,ntriag,v
   type(pdata_type), pointer :: pdata
   type(VoronoiCell) :: Vi
   type(Polygon) :: Wij,Pj,frame,Viframe
   
   ! LOCAL ARRAYS
   integer :: indices(n/2),sitetv(3*2*(n/2)),start((n/2)+1),til(3,2*(n/2)),tnbr(3,2*(n/2))
   real(kind=8) :: sites(2,n/2),vorxy(2,2*(n/2))

   integer :: flag,find
   logical :: debug
   real(kind=8) :: ctol
   integer :: cellcolorJ1(n/2)
   real(kind=8) :: f(0:3,n/2)
   character(len=80) :: paintcolor

   call c_f_pointer(pdataptr,pdata) 

   ! The smaller the ctol for which the cell has positive f, the darker the color.
   cellcolorJ1(1:n/2) = -1       ! Means no color

   debug = .false.
   call compallJ(n,x,f,debug,flag,pdataptr)
   if ( flag .ne. 0 ) then
      write(*,*) 'In drawsol, compallJ returned a non-null flag!'
      return
   end if

   where ( abs( f(0,1:n/2) ) .gt. 1.0d-03 )
      cellcolorJ1(1:n/2) = 1
   end where

   call c_f_pointer(pdataptr,pdata)
  
   nsites = n / 2
   indices(1:nsites) = (/ (i,i=1,nsites) /)
   sites(1:2,1:nsites) = reshape( x(1:n), (/ 2, nsites /) )
 
   if ( nsites .le. 1 ) then
      write(*,*) 'In drawsol: nsites must be at least 2!'
      stop
   end if

   ! Compute Delaunay triangulation

   call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

   if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
      write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225: ',ierror
      stop
   end if
 
   call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
   
   ! Compute vertices of the Voronoi diagram
         
   if ( ntriag .gt. 2 * nsites ) then
      write(*,*) 'In drawsol: There is something wrong. According to the dtris2, '
      write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
      write(*,*) 'cannot be larger than twice the number of sites!'
      stop
   end if
 
   do k = 1,ntriag
      call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
      if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
         write(*,*) 'In drawsol: triangle_circumcenter_2d returned NaN.'
         write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
         stop
      end if
   end do

   maxnv = maxval( pdata%nvpols(1:pdata%npols) )

   allocate(Pj%v(2,maxnv),Pj%e(maxnv-1),stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Allocation error.'
      stop
   end if
     
   open(unit=10,file=trim(filename))                           

   write(10,10) 1,pdata%drawscale

   ! Drawing Delaunay triangulation
   
   do k = 1,ntriag
      do v = 1,3
         write(10,20) sites(1:2,til(v,k)),sites(1:2,til(mod(v,3)+1,k)),'teagreen'
      end do
   end do

   ! Frame
   
   frame%n = 5
   frame%deg = .false.
   
   allocate(frame%v(2,frame%n),frame%e(frame%n-1),stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Allocation error.'
      stop
   end if

   frame%v(1,1:5) = (/ pdata%xmin, pdata%xmax, pdata%xmax, pdata%xmin, pdata%xmin /)
   frame%v(2,1:5) = (/ pdata%ymin, pdata%ymin, pdata%ymax, pdata%ymax, pdata%ymin /)
   frame%e(1:4) = (/ 0, 0, 0, 0 /)  
          
   ! Drawing Voronoi diagram
      
   do i = 1,nsites
      call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
      if ( status .ne. 0 ) then
         write(*,*) 'In drawsol: Error in voronoi_cell, status = ',status
         stop
      end if
  
      ! This is to draw the Voronoi cell restricted to the frame.
      
      if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
         call VorCellInterConvPol(sites(1:2,i),frame,Vi,Viframe)
      else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
         call VorCellInterConvPol_Collinear(i,nsites,sites,frame,Viframe)
      end if

      if ( .not. Viframe%deg ) then
         do v = 1,Viframe%n - 1
            write(10,20) Viframe%v(1:2,v),Viframe%v(1:2,v+1),'jonquil'
         end do
      end if

      call polygon_destroy(Viframe)

         ! This is to paint the Voronoi cell, restricted to A, if
      ! desired conditions does not hold.

      if ( cellcolorJ1(i) .ne. - 1 ) then
         if ( cellcolorJ1(i) .eq. 1 ) then
            paintcolor = 'dollarbill'
         end if

         do j = 1,pdata%npols
            Pj%n = pdata%nvpols(j)
            Pj%deg = .false.         
            
            Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
            Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
            Pj%e(1:Pj%n-1)  = pdata%poledges(j,1:Pj%n-1)

            if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
               call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
            else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
               call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
            end if

            if ( .not. Wij%deg ) then
               do v = 1,Wij%n - 1
                  if ( v .eq. 1 ) then
                     write(10,30) Wij%v(1:2,v),Wij%v(1:2,v+1)
                  else if ( v .gt. 1 .and. v .lt. Wij%n - 1 ) then 
                     write(10,31) Wij%v(1:2,v),Wij%v(1:2,v+1)
                  else
                     write(10,32) Wij%v(1:2,v),trim(paintcolor)
                  end if
               end do
            end if

            call polygon_destroy(Wij)
         end do
      end if

      ! This is to draw the Voronoi cell restricted to A.

      do j = 1,pdata%npols
         Pj%n = pdata%nvpols(j)
         Pj%deg = .false.

         Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
         Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
         Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

         if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
            call VorCellInterConvPol(sites(1:2,i),Pj,Vi,Wij)
         else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
            call VorCellInterConvPol_Collinear(i,nsites,sites,Pj,Wij)
         end if

         if ( .not. Wij%deg ) then
            do v = 1,Wij%n - 1
               if ( Wij%e(v) .lt. 1 ) then
                  write(10,20) Wij%v(1:2,v),Wij%v(1:2,v+1),'desire'
               end if
            end do
         end if

         call polygon_destroy(Wij)
      end do



   ! This is to draw the sites
   write(10,40) sites(1:2,i)
   
   call voronoi_cell_destroy(Vi)
   end do

   deallocate(frame%v,frame%e,stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Deallocation error.'
      stop
   end if
       
   ! Drawing A
   
   do j = 1,pdata%npols
      Pj%n = pdata%nvpols(j)
      Pj%deg = .false.          

      Pj%v(1,1:Pj%n) = pdata%xpol(j,1:Pj%n)
      Pj%v(2,1:Pj%n) = pdata%ypol(j,1:Pj%n)
      Pj%e(1:Pj%n-1) = pdata%poledges(j,1:Pj%n-1)

      do k = 1,Pj%n - 1
         if ( Pj%e(k) .eq. 0 ) then
            write(10,20) Pj%v(1:2,k),Pj%v(1:2,k+1),'black'
         end if
      end do
   end do

   ! Close figure

   write(10,50)

   close(10)

   deallocate(Pj%v,Pj%e,stat=allocerr)
   if ( allocerr .ne. 0 ) then
      write(*,*) 'In drawsol: Deallocation error.'
      stop
   end if

   ! NON-EXECUTABLE STATEMENTS
   
10  format('prologues := 3;',/, &
        'outputtemplate := "%j.mps";',/, &
        'input mpcolornames;',/, &
        'beginfig(',i4,');',/, &
        'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
        'orange, plum, bluea, blueb, bluec, blued, bluee, corn, deepsaffron, ferrarired;',/, &
        'skyblue      = (135/255,206/255,235/255);',/, &
        'dollarbill   = (133/255,187/255,101/255);',/, &
        'teagreen     = (208/255,240/255,192/255);',/, &
        'desire       = (234/255, 60/255, 83/255);',/, &
        'jonquil      = (244/255,202/255, 22/255);',/, &
        'royalblue    = ( 65/255,105/255,225/255);',/, &
        'orange       = (255/255,165/255,  0/255);',/, &
        'plum         = (221/255,160/255,221/255);',/, &
        'bluea        = (  0/255,230/255,230/255);',/, &
        'blueb        = (  0/255,191/255,230/255);',/, &
        'bluec        = (  0/255,153/255,230/255);',/, &
        'blued        = (  0/255,115/255,230/255);',/, &
        'bluee        = (  0/255, 77/255,230/255);',/, &
        'corn         = (255/255,230/255,102/255);',/, &
        'deepsaffron  = (255/255,153/255, 51/255);',/, &
        'ferrarired   = (255/255, 42/255,  0/255);',/, &
        'u = ',f20.10,'cm;')
20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
30  format('fill ((',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
31  format('      (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u)--')
32  format('      (',f20.10,'u,',f20.10,'u)--cycle) withcolor ',a30,';')              
40  format('drawdot(',f20.10,'u,',f20.10,'u);') 
50  format('endfig;',/,'end;')

end subroutine drawsolJ1
  
end module optvormod2
