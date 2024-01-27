!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! prob. /18,6,3/	 & 	polak2                         &   10 &    9 &    2 &   3823 &    9.183560E+01 &    1.000901E+04 &    5.459815E+01 &    0.000000E+00 & - \\\hline
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

subroutine setdim(n,m)
	implicit none
	integer n,m,q
	n=10
	m=9
	return
end subroutine setdim

subroutine STARTP(N,X)
	integer n
	real*8 x(n)
	x(1) =100.0d0
	x(2) =  0.1d0
	x(3) =  0.1d0
	x(4) =  0.1d0
	x(5) =  0.1d0
	x(6) =  0.1d0
	x(7) =  0.1d0
	x(8) =  0.1d0
	x(9) =  0.1d0
	x(10)=  0.1d0

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

subroutine funob(n,x,f)
	implicit none
	integer n, i
	real*8 x(n),z(n),app(2),e,f

	z    = x

	z(2) = z(2) + 2.0d0
	e = (1.0d-4*z(2))**2 + z(2)**2 + z(3)**2 + (2.0d0*z(4))**2
	do i = 5,10
		e = e + z(i)**2
	enddo
	app(1) = exp(e)


	z(2) = z(2) - 4.0d0
	e = (1.0d-4*z(2))**2 + z(2)**2 + z(3)**2 + (2.0d0*z(4))**2
	do i = 5,10
		e = e + z(i)**2
	enddo
	app(2) = exp(e)

	f = maxval(app)

	return
end subroutine funob

subroutine fconstr(n,m,x,constr)
	implicit none
	integer :: n, m, i
	real*8 :: x(n), constr(m)

	do i = 1,m
		constr(i) = x(i)**2 + x(i+1)**2 + x(i)*x(i+1) - 1.0D+00
	enddo
	return
end subroutine fconstr
