!!$program test
!!$
!!$  ! LOCAL SCALARS
!!$  integer, parameter :: n = 20
!!$
!!$  ! LOCAL ARRAYS
!!$  integer :: ind(n),v(n)
!!$
!!$  v(1:n) = (/ 1, 3, 5, 8, 2, -10, 14, 36, 24, 3, 5, 6, 90, -17, 19, 25, 5, 3, 2, 1 /)
!!$
!!$  call mergesort(n,v,ind)
!!$
!!$  write(*,*) 'ind    = ',ind(1:n)
!!$  write(*,*) 'v(ind) = ',v(ind(1:n))
!!$  
!!$end program test

subroutine voronoi_cell_patch(nsites,ntriag,til,start,sitetv)

  implicit none
        
  ! SCALAR ARGUMENTS
  integer, intent(in) :: nsites,ntriag
  
  ! ARRAY ARGUMENTS
  integer, intent(in) :: til(3,ntriag)
  integer, intent(out) :: sitetv(3*ntriag),start(nsites+1)

  ! LOCAL SCALARS
  integer :: i,ind,itriag,iverti,j,k

  ! till(3,ntriag) contains, for each vertex of a triangle, the index
  ! of the site to which the vertex corresponds
  
  call mergesort(3*ntriag,til,sitetv)

  ! sitetv(3*ntriag). Position ind = 3 * ( i - 1 ) + j of array sitetv
  ! corresponds to vertex j of triangle i. Given ind, we can get i and
  ! j by i = ( ind - 1 ) / 3 + 1 and j = ind - 3 * ( i - 1 ).

  ! i=1 j=1,2,3 --> ind=1,2,3
  ! i=2 j=1,2,3 --> ind=4,5,6
  ! i=3 j=1,2,3 --> ind=7,8,9
  
  start(1:nsites) = - 1
  start(nsites+1) = 3 * ntriag + 1
    
  do i = 1,ntriag
     do j = 1,3
        ind = 3 * ( i - 1 ) + j

        itriag = ( sitetv(ind) - 1 ) / 3 + 1
        iverti = sitetv(ind) - 3 * ( itriag - 1 )

        if ( itriag .lt. 1 .or. itriag .gt. ntriag ) then
           write(*,*) 'There is something wrong in array itriag, computed in voronoi_cell_patch!'
           stop
        end if
           
        if ( iverti .lt. 1 .or. iverti .gt. 3  ) then
           write(*,*) 'There is something wrong in array iverti, computed in voronoi_cell_patch!'
           stop
        end if
           
        k = til(iverti,itriag)
        
        if ( k .lt. 1 .or. k .gt. nsites ) then
           write(*,*) 'There is something wrong in array til, computed by dtris2, detected in voronoi_cell_patch!'
           write(*,*) 'Value of itriag, iverti, til(ivert,itriag), nsites = ',itriag,iverti,til(iverti,itriag),nsites
           stop
        end if
          
        if ( start(k) .eq. - 1 ) then
           start(k) = ind
        end if
     end do
  end do
  
end subroutine voronoi_cell_patch

subroutine mergesort(n,v,ind)

  ! Slightly adapted from: https://github.com/Astrokiwi/simple_fortran_argsort/blob/master/sort_test.f90
  
  implicit none
        
  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  integer, intent(in) :: v(n)
  integer, intent(out) :: ind(n)

  ! LOCAL SCALARS
  integer :: i,j,left,k,ksize,stepsize

  ! LOCAL ARRAYS
  integer :: tmp(n)

  ind(1:n) = (/ (i,i=1,n) /)
            
  if ( n .eq. 1 ) return
            
  stepsize = 1
  do while ( stepsize .lt. n )
     do left = 1,n - stepsize,2 * stepsize
        i = left
        j = left + stepsize
        ksize = min(2 * stepsize,n - left + 1)
        
        k = 1
        
        do while ( i .lt. left + stepsize .and. j .lt. left + ksize )
           if ( v(ind(i)) .le. v(ind(j)) ) then
              tmp(k) = ind(i)
              i = i + 1
              k = k + 1
           else
              tmp(k) = ind(j)
              j = j + 1
              k = k + 1
           end if
        end do
        
        if ( i .lt. left + stepsize ) then
           tmp(k:ksize) = ind(i:left+stepsize-1)
        else
           tmp(k:ksize) = ind(j:left+ksize-1)
        endif

        ind(left:left+ksize-1) = tmp(1:ksize)
     end do
     
     stepsize = 2 * stepsize
  end do
  
end subroutine mergesort
 
