module voro

  use VorCells_Polygons, only: Polygon,VoronoiCell,voronoi_cell,VorCellInterConvPol, &
       VorCellInterConvPol_Collinear,voronoi_cell_destroy,polygon_destroy

  implicit none
  
  type :: pdata_type
     integer :: npols
     real(kind=8) :: c(2),rmax,volA,xmin,xmax,ymin,ymax,drawscale
     integer, allocatable :: nvpols(:),poledges(:,:)
     real(kind=8), allocatable :: xpol(:,:),ypol(:,:)
  end type pdata_type
  
  public :: voronoi,drawvor

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine voronoi(nsites,sites,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop)
	
    implicit none
  
    ! SCALAR ARGUMENTS
    integer, intent(in) :: An,nsites,nvmax
    integer, intent(out) :: istop,nv
    
    ! ARRAY ARGUMENTS
    integer, intent(in) :: Aflag(An+1)
    real(kind=8), intent(in) :: Ax(An+1),Ay(An+1),sites(2,nsites)
    integer, intent(out) :: sstart(nsites+1),vflag(nvmax)
    real(kind=8), intent(out) :: vx(nvmax),vy(nvmax)
    
    ! LOCAL SCALARS
    integer :: allocerr,i,ierror,ind,j,k,status,ntriag
    type(VoronoiCell) :: Vi
    type(Polygon) :: Wi,P
    
    ! LOCAL ARRAYS
    integer :: indices(nsites),sitetv(3*2*nsites),start(nsites+1),til(3,2*nsites),tnbr(3,2*nsites)
    real(kind=8) :: vorxy(2,2*nsites)
    
   istop = 0

    if ( nsites .lt. 3 ) then
      write(*,*) 'In voronoi: nsites must be at least 3!'
      stop
    end if

    ! Polygon representing the domain

    P%n = An + 1
    P%deg = .false.

    allocate(P%v(2,P%n),P%e(P%n-1),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'In voronoi: Allocation error.'
       stop
    end if
    
    P%v(1,1:P%n) = Ax(1:An+1)
    P%v(2,1:P%n) = Ay(1:An+1)
    P%e(1:P%n-1) = Aflag(1:An)
          
    ! Compute Delaunay triangulation

    indices(1:nsites) = (/ (i,i=1,nsites) /)

    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

    if ( ierror .ne. 0 ) then
      istop = -1
      return
    end if
    
    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
       write(*,*) 'In voronoi: dtris2 returned an error different from 0 and 225: ',ierror
       istop = -1
    end if
  
    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
    
    ! Compute Voronoi diagram
          
    if ( ntriag .gt. 2 * nsites ) then
       write(*,*) 'In voronoi: There is something wrong. According to the dtris2, '
       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
       write(*,*) 'cannot be larger than twice the number of sites!'
       stop
    end if
  
    do k = 1,ntriag
       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
          write(*,*) 'In voronoi: triangle_circumcenter_2d returned NaN.'
          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
          istop = -2
       end if
    end do

    sstart(1) = 1

    do i = 1,nsites
       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
       if ( status .ne. 0 ) then
          write(*,*) 'In voronoi: Error in voronoi_cell, status = ',status
          stop
       end if

       ! Compute the intersection of the ith Voronoi cell and the domain

       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
          call VorCellInterConvPol(sites(1:2,i),P,Vi,Wi)
       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
          call VorCellInterConvPol_Collinear(i,nsites,sites,P,Wi)
       end if
      !  write(*,*) '-------'
      !  write(*,*) Wi%deg
      !  write(*,*) Wi%v
      !  write(*,*) '-------'
       if ( Wi%deg ) then
          sstart(i+1) = sstart(i)

       else
          sstart(i+1) = sstart(i) + ( Wi%n - 1 )
          
          if ( sstart(i+1) .gt. nvmax + 1 ) then
             write(*,*) 'Increase nvmax.'
             stop
          end if
    
          do j = 1,Wi%n - 1
             ind = sstart(i) + j - 1
             vx(ind) = Wi%v(1,j)
             vy(ind) = Wi%v(2,j)
             vflag(ind) = -Wi%e(j)
          end do
       end if
       
       call polygon_destroy(Wi)

       call voronoi_cell_destroy(Vi)
    end do

    nv = sstart(nsites+1) - 1
    
    deallocate(P%v,P%e,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'In voronoi: Deallocation error.'
       stop
    end if

  end subroutine voronoi

  ! ******************************************************************
  ! ******************************************************************


  subroutine drawvor(nsites,sites,cx,cy,An,Ax,Ay,sstart,nv,vx,vy,filename)

    implicit none
 
    ! SCALAR ARGUMENTS
    integer, intent(in) :: An,nsites,nv
   
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: Ax(An+1),Ay(An+1),sites(2,nsites),vx(nv),vy(nv),cx(nsites),cy(nsites)
    integer, intent(in) :: sstart(nsites+1)
    character(len=80) :: filename
   
    ! LOCAL SCALARS
    integer :: i,j,jnext,k,ierror,ntriag,v,status,allocerr
    real(kind=8) :: marg,xmin,xmax,ymin,ymax
    type(VoronoiCell) :: Vi
    type(Polygon) :: frame,Viframe

    ! LOCAL ARRAYS 
    integer :: indices(nsites),til(3,2*nsites),tnbr(3,2*nsites),sitetv(3*2*nsites),start(nsites+1)
    real(kind=8) :: vorxy(2,2*nsites)


   ! Compute Delaunay triangulation

    indices(1:nsites) = (/ (i,i=1,nsites) /)

    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
      write(*,*) 'In drawsol: dtris2 returned an error different from 0 and 225',ierror
      stop
    end if

    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)

    ! Compute the vertices of the Voronoi diagram

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

   
    open(unit=10,file=trim(adjustl(filename)))                          

    write(10,10)

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
    marg = 0.1d0
    xmin = 0.0d0 - marg
    xmax = 1.0d0 + marg
    ymin = 0.0d0 - marg
    ymax = 1.0d0 + marg
    write(*,*) 'frame'

    frame%v(1,1:5) = (/ xmin, xmax, xmax, xmin, xmin /)
    frame%v(2,1:5) = (/ ymin, ymin, ymax, ymax, ymin /)
    frame%e(1:4) = (/ 0, 0, 0, 0/)  

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

      ! Draw the cell (intersected with the domain)

       do j = sstart(i),sstart(i+1) - 1
          if ( j .lt. sstart(i+1) - 1 ) then
             jnext = j + 1
          else
             jnext = sstart(i)
          end if
          write(10,20) vx(j),vy(j),vx(jnext),vy(jnext),'desire'
       end do
   ! Draw the ith site
      write(10,30) sites(1:2,i),'black' 
      ! write(10,31) i,sites(1:2,i)
      write(10,30) cx(i),cy(i),'orange'
    end do

    deallocate(frame%v,frame%e,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'In drawsol: Deallocation error.'
       stop
    end if

     
    ! Draw the domain
    do k = 1,An - 1
       write(10,20) Ax(k),Ay(k),Ax(k+1),Ay(k+1),'black'
    end do

    ! Close figure
   
    write(10,40)
    
    close(10)

    ! NON-EXECUTABLE STATEMENTS

    10  format('prologues := 3;',/, &
    'outputtemplate := "%j.mps";',/, &
    'input mpcolornames;',/, &
    'beginfig( 1 );',/, &
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
    'u = 6cm;')
20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
30  format('drawdot(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 2 withcolor',a12,';')
!31  format('dotlabel.bot(btex ',i2,' etex, (',f20.10,'u,',f20.10,'u));')
40  format('endfig;',/,'end;')

  end subroutine drawvor
  
end module voro
