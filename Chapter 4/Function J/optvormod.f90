module optvormod

    use iso_c_binding, only: c_ptr,c_f_pointer
    use VorCells_Polygons, only: Polygon,VoronoiCell,voronoi_cell,voronoi_cell_old,VorCellInterConvPol, &
         VorCellInterConvPol_Collinear,voronoi_cell_destroy,polygon_destroy
  
    implicit none
    
    type :: pdata_type
       integer :: npols
       real(kind=8) :: c(2),rmax,volA,xmin,xmax,ymin,ymax,drawscale
       integer, allocatable :: nvpols(:),poledges(:,:)
       real(kind=8), allocatable :: xpol(:,:),ypol(:,:)
       ! 0: minimum distance between sites
       ! 1: cells with identical volume given by vol(A)/m
       ! 2: cells with volume proportional to density function f evaluated at their sites
       ! 3: each cell with equally sized edges
       ! 4: each cell with equally sized angles
       ! 5: each cell with edges being crossed at their middle point by their corresponding Delaunay edges
    end type pdata_type
    
    public :: evalfg,approxg,drawsol,proj
  
  contains
  
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

      do i = 1,nsites
         tmp = norm2( x(2*i-1:2*i) - pdata%c(1:2) )
         if ( tmp .gt. pdata%rmax ) then
            x(2*i-1:2*i) = pdata%c(1:2) + pdata%rmax * ( x(2*i-1:2*i) - pdata%c(1:2) ) / tmp
         end if
      end do
   
    end subroutine proj
  
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
        
            call compJC(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return

            
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
  
    end subroutine evalfg

    ! ******************************************************************
    ! ******************************************************************

    subroutine evalf(n,x,f,flag,pdataptr)

      implicit none

      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
      real(kind=8), intent(out) :: f
      type(c_ptr), optional, intent(in) :: pdataptr

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)

      ! PARAMETERS
      logical, parameter :: compgrad = .false., debug= .false.

      ! LOCAL ARRAYS
      real(kind=8) :: g(n)

      call evalfg(n,x,f,g,flag,pdataptr)
   
   end subroutine evalf
      
    ! ******************************************************************
    ! ******************************************************************

   subroutine evalg(n,x,g,flag,pdataptr)

      implicit none
    
      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(out) :: flag
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

  
      g(1:n) = 0.0d0
      
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
        else
            if ( debug ) then
               write(*,*)
               write(*,*) 'POLYGON Wi with i = ',i
               write(*,*) 'Number of vertices: ',W%n - 1
               write(*,*)
            end if
        
            call compJC(nsites,sites,i,W,ftmp,gtmpnnzmax,gtmpnnz,gtmpind,gtmpval,flag,compgrad,debug)
            if ( flag .ne. 0 ) return  

            g(gtmpind(1:gtmpnnz)) = g(gtmpind(1:gtmpnnz)) + gtmpval(1:gtmpnnz)
         end if

        call voronoi_cell_destroy(V)
        
        call polygon_destroy(W)

      end do
      
      g(1:n) = g(1:n) / nsites
      
      deallocate(P%v,P%e,stat=allocerr)
      if ( allocerr .ne. 0 ) then
         write(*,*) 'In evalfg: Deallocation error.'
         flag = - 1
         return
      end if
   
   end subroutine evalg

    ! ******************************************************************
    ! ******************************************************************

    subroutine compJC(nsites,sites,i,W,f,gnnzmax,gnnz,gind,gval,flag,compgrad,debug)
    
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
      real(kind=8) :: glcnst,volW,phia,phib,phic,xi,re,pd

      ! LOCAL ARRAYS
      real(kind=8) :: gradphi(2),M(2,2),wrot(2),vrot(2),cent(2),p(2),diff(2)
      
      call centroid(W,cent,volW) 

      if ( debug ) then
         write(*,*)
         write(*,*) 'volA = ', volW
         write(*,*) 'cent = ', cent(1:2)
      end if

      xi = 2.0d0/(6.0d0*volW)

      if ( debug ) then
         write(*,*)
         write(*,*) 'xi = ', xi
         write(*,*)
      end if

      diff(1:2) = sites(1:2,i) - cent(1:2)
      f = sum(diff(1:2) ** 2)

      if ( compgrad ) then  
         glcnst = 0.0d0       
         do ed = 1,W%n-1                      ! This calculates the value of the sum_{E\in E_i}(r_E*(p_E)^T(diff))
            if ( ed .lt. W%n - 1) then
               ednext = ed + 1
            else
               ednext = 1
            end if

            re = W%v(1,ed)*W%v(2,ednext) - W%v(1,ednext)*W%v(2,ed)   ! re = v_{E_1}w_{E_2} - w_{E_1}v_{E_2} 
            p(1:2) = W%v(1:2,ed) + W%v(1:2,ednext)                   ! p = v_{E} + w_{E}
            pd = dot_product(p(1:2),diff(1:2))
            glcnst = glcnst + re*pd
         end do

         gnnz = 0

         if ( gnnz + 2 .gt. gnnzmax ) then
            write(*,*) 'In compJC: Increase gnnzmax to at least ',gnnz+2
            flag = - 1
            return
         end if
         
         gind(gnnz+1:gnnz+2) = (/ 2*i-1, 2*i /)
         gval(gnnz+1:gnnz+2) = 6.0d0*volW * xi * diff(1:2)
         gnnz = gnnz + 2 
         


         do ed = 1,W%n-1
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
            else if ( ed .eq. W%n-1 ) then
               edprev = ed - 1
               ednext = 1
            else 
               edprev = ed - 1
               ednext = ed + 1
            end if

            if ( debug ) then
               write(*,*) 'Previous: ', edprev, 'next: ',ednext
               write(*,*) 'They corresponds to the points ',W%v(1:2,edprev),' and ', W%v(1:2,ednext),' respectively.'
            end if

            re = W%v(1,ed)*W%v(2,ednext) - W%v(1,ednext)*W%v(2,ed)
            p(1:2) = W%v(1:2,ed) + W%v(1:2,ednext)
            pd = dot_product(p(1:2),diff(1:2))
            vrot(1:2) = (/ -W%v(2,ed), W%v(1,ed) /)
            wrot(1:2) = (/ W%v(2,ednext), -W%v(1,ednext) /)

            if ( debug ) then
               write(*,*)
               write(*,*) 're= ',re
               write(*,*) 'p= ',p
               write(*,*) 'vrot= ',vrot
               write(*,*) 'wrot= ',wrot
               write(*,*)
            end if

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

               gval(1:2) = gval(1:2) - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), wrot(1:2)))

               call matrixint(sites(1:2,i),sites(1:2,abs(W%e(edprev))),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return
            
               if ( gnnz + 2 .gt. gnnzmax) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ',gnnz + 2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), wrot(1:2) ))
               gnnz = gnnz + 2

               call matrixint(sites(1:2,abs(W%e(ed))),sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ',gnnz + 2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), wrot(1:2) ))
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

               gval(1:2) = gval(1:2) - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
               (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), wrot(1:2)))

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), wrot(1:2) ))
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

               gval(1:2) = gval(1:2) - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
               (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), wrot(1:2)))

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(edprev))),W%v(1:2,ed),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(edprev))-1, 2*abs(W%e(edprev)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), wrot(1:2) ))
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

               gval(1:2) = gval(1:2) - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
               (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), vrot(1:2) ))

               call matrixint(sites(1:2,i),sites(1:2,abs(W%e(ed))),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), vrot(1:2) ))
               gnnz = gnnz + 2

               call matrixint(sites(1:2,abs(W%e(ednext))),sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), vrot(1:2) ))
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

               gval(1:2) = gval(1:2) - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
               (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), vrot(1:2) ))

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ed))),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ed))-1, 2*abs(W%e(ed)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), vrot(1:2) ))
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

               gval(1:2) = gval(1:2) - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
               (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), vrot(1:2) ))

               call matrixbd(sites(1:2,i),sites(1:2,abs(W%e(ednext))),W%v(1:2,ednext),gradphi,M,flag)
               if ( flag .ne. 0 ) return

               if ( gnnz + 2 .gt. gnnzmax ) then
                  write(*,*) 'In compJC: Increase gnnzmax to at least ', gnnz+2
                  flag = -1
                  return
               end if

               gind(gnnz+1:gnnz+2) = (/ 2*abs(W%e(ednext))-1, 2*abs(W%e(ednext)) /)
               gval(gnnz+1:gnnz+2) = - xi * (re * matmul( transpose( M(1:2,1:2) ), diff(1:2) ) + &
                (pd - 0.5d0*glcnst/volW) * matmul( transpose( M(1:2,1:2) ), vrot(1:2) ))
               gnnz = gnnz + 2
            end if
         end do
      end if
   
    end subroutine compJC

    ! ******************************************************************
    ! ******************************************************************
   
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
      
    ! ******************************************************************
    ! ******************************************************************
  
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
      real(kind=8) :: f(0:5,n/2),sites(2,n/2),vorxy(2,2*(n/2))
  
      debug = .false.
      !call compallJ(n,x,f,debug,flag,pdataptr)
      if ( flag .ne. 0 ) then
         write(*,*) 'In drawsol, compallJ returned a non-null flag!'
         return
      end if
      
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
      
      write(10,*) 'Value of J2i (sum/ave/min/max) = ',sum(abs(f(2,1:n/2))),sum(abs(f(2,1:n/2)))/(n/2), &
           minval(abs(f(2,1:n/2))),maxval(abs(f(2,1:n/2))),' J2 itself = ',sum(f(2,1:n/2)**2)/(n/2)
  
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
  
  end module optvormod