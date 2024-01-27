subroutine setdim(n,m)
	implicit none
	integer n,m,q
	n=6
	m=0
	return
end subroutine setdim

subroutine STARTP(N,X)
	integer n
	real*8 x(n)
	x(1) = 2.d0
	x(2) = 2.d0
	x(3) = 7.d0
	x(4) = 0.d0
	x(5) =-2.d0
	x(6) = 1.d0

	return
end subroutine startp

subroutine setbounds(n,bl,bu)
	implicit none
	integer :: n
	real*8 :: bl(n), bu(n)

	bl = -100.d0
	bu = +100.d0

	return
end subroutine setbounds

subroutine fconstr(n,m,x,constr)
	implicit none
	integer :: n, m, i
	real*8 :: x(n), constr(m)

	return
end subroutine fconstr

subroutine funob(n,x,f)
	implicit none
	integer n,indfun
	real*8  x(n),f
	integer i
	real*8 fapp(102),y(51),t(51)

	do i = 1,51
		t(i) = (dble(i) - 1.0d0) / 10.0d0
		y(i) = exp(t(i))/2.0d0 - exp(-2.0d0*t(i))
		y(i) = y(i) + exp(-3.0d0*t(i))/2.0d0
		y(i) = y(i) + 1.5d0*exp(-1.5d0*t(i))*sin(7.0d0*t(i))
		y(i) = y(i) +       exp(-2.5d0*t(i))*sin(5.0d0*t(i))
		fapp(i) = x(1)*exp(-x(2)*t(i))*cos(x(3)*t(i) + x(4))
		fapp(i) = fapp(i) + x(5)*exp(-x(6)*t(i)) - y(i)
	enddo

	do i = 52,102
		fapp(i) = -fapp(i-51)
	enddo
	f = maxval(fapp)
	
	return
end subroutine funob

