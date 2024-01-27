module vincoli
	integer m
	double precision, allocatable :: eps(:),constr(:),epsiniz(:)
	double precision viol,violz
end module vincoli

program main_vincolato
	use vincoli
	implicit none

	integer :: n
	integer :: j,istop,i,hschoice
	real     :: tbegin, tend
	character*30 :: nomefun
	real*8, allocatable :: x(:),bl(:),bu(:),step(:)      


	integer ::            num_funct,num_iter
	real*8             :: f, fob, violiniz, finiz, alfa_stop
	integer            :: nf_max,iprint

!------------------------------------------------------------------------------

	call setdim(n,m)

	allocate(x(n),bl(n),bu(n),step(n))
	if(m >= 1) then
	      allocate (constr(m),epsiniz(m),eps(m))
	endif

	call setbounds(n, bl,bu)
	call startp(n, x)

	do i=1,n
		if((x(i).lt.bl(i)).or.(x(i).gt.bu(i))) then
		   write(*,*) ' starting point violates bound constraints'
		   stop
		endif
	enddo

!-----------------------------------------------------------------------
!     print starting point info
!-----------------------------------------------------------------------
  
	call funct_vinc(n,m,x,fob,constr)

	write(*,*) ' ------------------------------------------------- '
	write(1,*) ' ------------------------------------------------- '

	write(*,*) '      f(xo) = ',fob
	write(1,*) '      f(xo) = ',fob

	viol=0.d0

	do i = 1,m
		viol=viol+max(0.d0,constr(i))
	enddo

	write(*,*) '  cviol(xo) = ',viol
	write(1,*) '  cviol(xo) = ',viol


	write(*,*) ' ------------------------------------------------- '
	write(1,*) ' ------------------------------------------------- ' 		      

	do i=1,n
		write(*,*) ' xo(',i,') =',x(i)
		write(1,*) ' xo(',i,') =',x(i)
	enddo

	write(*,*) ' ------------------------------------------------- '
	write(1,*) ' ------------------------------------------------- '
       
!-----------------------------------------------------------------------
!	choice of starting penalty parameter values
!-----------------------------------------------------------------------
	do i = 1,m
		if(max(0.d0,constr(i)) < 1.d-0) then
			eps(i) = 1.d-3
		else
			eps(i) = 1.d-1
		endif
	enddo

	if(m >= 1) epsiniz     = eps
	finiz       = fob
	viol = 0.d0
	if(m >= 1) viol = max(0.d0,maxval(constr))
	violiniz    = viol

        num_funct   = 1
        num_iter    = 0
        alfa_stop=1.d-6
  	nf_max=20000
        iprint=0
        hschoice = 2
	!---------------------------------------------------
	! hschoice = 1 ---> use Halton pseudorandom seq.
	! hschoice = 2 ---> use Sobol  pseudorandom seq.
	!---------------------------------------------------

	write(*,*)
	write(*,*) ' call the optimizer ...'

	call cpu_time(tbegin)
        call sd_box(n,x,f,bl,bu,alfa_stop,nf_max,num_iter,num_funct,iprint,istop,hschoice)
	call cpu_time(tend)

	write(*,*) ' ... done'
	write(*,*)
	 
	call funct_vinc(n,m,x,fob,constr) 
        num_funct = num_funct+1
      
        write(*,*) ' ------------------------------------------------------------------------------'     
	write(1,*) ' ------------------------------------------------------------------------------'     

	write(*,*) ' total time:',tend-tbegin
	write(1,*) ' total time:',tend-tbegin
        write(*,*) ' number of function evaluations = ',num_funct 
        write(1,*) ' number of function evaluations = ',num_funct     		    
	write(*,*) ' ------------------------------------------------------------------------------'  
	write(1,*) ' ------------------------------------------------------------------------------'  

	write(*,*)
	write(1,*)

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- '

        write(*,*) '     f(xstar) = ',fob
        write(1,*) '     f(xstar) = ',fob

        if(m >= 1) then
		viol=max(0.d0,maxval(constr))

		write(*,*) ' cviol(xstar) = ',viol
		write(1,*) ' cviol(xstar) = ',viol

		write(*,*) ' ------------------------------------------------- '
		write(1,*) ' ------------------------------------------------- ' 		      
	endif
	do i=1,n
		write(*,*) ' xstar(',i,') =',x(i)
		write(1,*) ' xstar(',i,') =',x(i)
        enddo

        write(*,*) ' ------------------------------------------------- '
        write(1,*) ' ------------------------------------------------- '

	write(*,*) 
	write(1,*) 

	nomefun = 'Polak 2 + (n-1) constraints'

	write(2,135) nomefun,n,m,num_funct,finiz,violiniz,fob,viol

100  format(es15.6e3)
110  format(a15)
135  format(a30,2(' & ',i4),' & ',i6,4(' & ',es15.6),' & - \\\hline')
2002 format(2d20.10)


	deallocate(x,bl,bu,step)

	if(m.ge.1) deallocate (constr,epsiniz,eps)

end program main_vincolato

subroutine funct_vinc(n,m,x,fob,constr)
	implicit none
	integer :: n,m
	real*8 :: x(n),fob,constr(m)

	call funob(n,x,fob)
	if(m > 0) call fconstr(n,m,x,constr)

	return
end subroutine funct_vinc

subroutine i8_sobol_tmp(n8,index_sobol,d_dense)
	implicit none
	integer*8	:: n8, index_sobol
	real*8		:: d_dense(n8)

	integer	:: n,index_halton

	n = n8
	index_halton = index_sobol

	call halton(n,index_halton,d_dense)

	index_sobol = index_sobol + 2*n

	return
end subroutine i8_sobol_tmp

