!-*-f90-*-
! various aux routines
 subroutine findnaninf(vector,length,nan,inf)

   implicit none
   real*8 vector(length)
   integer length
   integer i

   logical nan,inf

   nan = .false.
   inf = .false.


   do i=1,length
      if(vector(i).ne.vector(i)) then
         write(*,*) "NaN at index ",i
         nan = .true.
      endif
      
   enddo
   

 end subroutine findnaninf
