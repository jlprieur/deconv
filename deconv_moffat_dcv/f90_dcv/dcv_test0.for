c***********************************************************
c Test of minimisation routines from "Numerical Recipees"
c
c JLP
c Version 12/03/2001
c***********************************************************
      program dcv_test0
      integer i,iter,n 
      parameter(n=2)
      real p(n),ftol,fret
      external func_test,dfunc_test
      external func,dfunc

      write(6,*) "Test with Banana function"

c Starting point:
      p(1)= -1.9 
      p(2)= 2.
c Tolerance
      ftol=1.e-3
      call frprmn(p,n,ftol,iter,fret,func,dfunc)
c      call frprmn(p,n,ftol,iter,fret,func_test,dfunc_test)
      write(6,11) p(1),p(2),iter,fret
11    format(" Location of the minimum:",2(F12.3,1X),/,
     1       " Number of iterations:",I4,/,
     1       " Value of the minimum:",E12.5)
      end
!-------------------------------------------------------------------
! Function: 
! Banana function : start at X=[-1.9;2]. minimum at X=[1;1] : f(X)=0;
!-------------------------------------------------------------------
      real function func_test(x)
      real x(*)
      func = (1. - x(1))**2 + 100*(x(2) - x(1))**2
      end
!-------------------------------------------------------------------
! Gradient:
! BANANA FUNCTION
! z = 100*(x(2)-x(1))^2 + (1-x(1))^2
!  dx = [  -200*(x(2) - x(1))-2*(1-x(1)); 200*(x(2)-x(1)) ];
!-------------------------------------------------------------------
      subroutine dfunc_test(x,dx)
      real x(*),dx(*)
      dx(1) = -200*(x(2)-x(1)) - 2*(1-x(1)) 
      dx(2) = 200*(x(2)-x(1))
      return
      end
      real function func(x)
      real x(*)
      func = (1. - x(1))**2 + 100*(x(2) - x(1))**2
      end
      subroutine dfunc(x,dx)
      real x(*),dx(*)
      dx(1) = -200*(x(2)-x(1)) - 2*(1-x(1)) 
      dx(2) = 200*(x(2)-x(1))
      return
      end
