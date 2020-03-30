!****************************************************************** 
! Program which generates a matrix with ones on the even lines
! and increasing integers along the odd lines
!
!        1 2 3 4
!        1 1 1 1
!        5 6 7 8
!
! JLP
! Version 16/03/2001
!****************************************************************** 
program exo2
 implicit none
 integer(kind=2)			:: n, m, err, i
 integer, dimension(:,:),allocatable	:: mat

print *,"Nber of lines? :"; read(*,*) n
print *,"Nber of columns? :"; read(*,*) m

! Allocate memory
allocate(mat(n,m),stat=err)
if (err /= 0) then
  print *,"Error allocating memory"; stop 4
endif

! Fill even lines with ones:
mat(2:n:2,:) = 1

! Fill odd lines with 1,2,...
mat(1:n:2,:) = reshape((/ (i,i=1,size(mat(1:n:2,2:))) /), &
                       shape=shape(mat(1:n:2,:)), order=(/ 2,1 /))
! Order is 1,2 by default (i.e. along the columns). Here we want
! it along the lines, i.e., the lines have the priority when filling
! the array

! Now print the output:
do i=1,n
 print *,mat(i,:)
end do

end program exo2
