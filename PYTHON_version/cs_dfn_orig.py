"""
!============================================================================================
!    CS-DFN - Derivative-Free program for Nonsmooth Nonlinear Programming
!    Copyright (C) 2013  G.Fasano, G.Liuzzi, S.Lucidi, F.Rinaldi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G.Fasano, G.Liuzzi, S.Lucidi, F.Rinaldi. A Linesearch-based Derivative-free Approach for
!    Nonsmooth Constrained Optimization, SIAM J. Optim. 24(3): 959-992, 2014
!    DOI: 10.1137/130940037
!============================================================================================
"""

import numpy as np
import ghalton
import sobol_seq

def cs_dfn(funct,m,x,bl,bu,tol,nf_max,iprint,hschoice):
	"""
	Parameters
	----------
	funct : function that computes the objective function value and
		optionally the constraint values (if m > 0).
	m : integer
		Number of inequality constraints, 0 if none.
	x : numpy array of n floats
		initial point
	bl : numpy array of floats
		lower bounds on the variables.
	bu : numpy array of floats
		upper bounds on the variables.
	tol : float
		tolerance in the stopping condition.
	nf_max : integer
		maximum number of allowed function evaluations.
	iprint : interger
		verbosity level.
	hschoice : integer
		type of dense direction used. 1 - halton, 2 - sobol.

	Returns
	-------
	x : numpy array of floats
		best point found by the optimizer.
	f : float
		best function value.

	"""

	istop = 0
	index_halton = 1000
	index_sobol  = 10000
	n = len(x)
	nf = 0
	ni = 0
	num_fal = 0
	soglia = 1.e-3
	flag_fail = [False for i in range(n)]
	fstop = np.zeros(2*n)
	xfstop= np.zeros((n,2*n))

	#####################################################
	# CHECK IF INITIAL POINT SATISFIES BOUND CONSTRAINTS
	#####################################################
	for i in range(n):
		if((x[i] < bl[i]) or (x[i] > bu[i])):
		   print('\n\n starting point violates bound constraints\n\n')
		   return x,np.Inf

	#####################################################
	# EPSILON INITIALIZATION FOR CONSTRAINTS HANDLING
	#####################################################
	if m > 0:
		fob,constr = funct(x)
	else:
		fob = funct(x)
	nf += 1
	eps = [1.0 for i in range(m)]
	for i in range(m):
		if(np.maximum(0.0,constr[i]) < 1.0):
			eps[i] = 1.e-3
		else:
			eps[i] = 1.e-1

	if m > 0:
		f = fob + np.maximum(0.0,np.max(constr/eps))
	else:
		f = fob

	sequencer = ghalton.Halton(n)
	Phalton   = sequencer.get(nf_max + 1000)

	alfa_d = np.zeros(n)
	alfa_dense = np.zeros(n)
	for i in range(n):
		alfa_d[i]   = np.maximum(np.double(1.e-3),np.minimum(np.double(1.0),np.abs(x[i])))
		for j in range(n):
			alfa_dense[j] += alfa_d[i]
		if iprint >= 1:
			print("alfainiz(%d)=%f" % (i,alfa_d[i]))

	alfa_dense /= np.double(n)
	alfa_diag   = np.copy(alfa_d)
	alfa_coord  = np.copy(alfa_d)
	alfa_max    = np.max(alfa_d)

	if n > 1:
		if hschoice == 1:
			d_dense = np.asarray(Phalton[index_halton-1], dtype = np.double)
		else:
			d_dense, index_sobol = sobol_seq.i4_sobol(n, index_sobol)

	d1 = np.zeros(n)
	direzioni = np.zeros((n,n))
	d = np.ones(n)
	for i in range(n):
		direzioni[i,i] = np.double(1.0)

	i_corr  = 0
	i_dense = 0
	j_dense = 0
	tipo_direzione = 0

	fstop   = np.zeros(2*n+1)
	xfstop  = np.zeros((n,2*n+1))
	for j in range(2*n+1):
		fstop[j] = f
	for i in range(n):
		for j in range(2*n+1):
			xfstop[i,j] = x[i]

	z = np.copy(x)

	if iprint >= 2:
		print(" ----------------------------------")
		print(" finiz =",f)
		for i in range(n):
			print(" xiniz(",i,")=",x[i])

	while True:
		if n > 1:
			alfa_max = np.max([np.max(alfa_coord),np.max(alfa_diag),np.max(alfa_dense)])
		else:
			alfa_max = np.max(alfa_d)

		istop = stop(alfa_d,alfa_max,nf,ni,tol,nf_max)

		alfa_max = np.max(alfa_d)

		if istop >= 1:
			break

		if i_corr == 0:
			dconv = np.zeros(n)
			for i in range(n):
				dconv += -direzioni[:,i]

		if iprint >= 1:
			print("----------------------------------------------")
		if iprint >= 0:
			print(" ni=%4d  nf=%5d   f=%12.5e   alfamax=%12.5e" % (ni,nf,f,alfa_max))
		if iprint >= 2:
			for i in range(n):
				print(" x(",i,")=",x[i])

		d = np.copy(direzioni[:,i_corr])

		if tipo_direzione == 0:
			alfa, fz, z1, fz1, z2, fz2, i_corr_fall, nf = linesearchbox_cont(funct,m,eps,x,f,d,alfa_d,z,i_corr,alfa_max,bl,bu,nf,iprint)

			if alfa >= 1.e-12:
				x[i_corr] = x[i_corr]+alfa*d[i_corr]
		else:

			alfa, alfatilde, fz, d, nf = linesearchbox_dense(funct,m,eps,x,f,d,alfa_d[i_corr],alfa_max,bl,bu,nf,iprint)
			alfa_d[i_corr]      = alfatilde

			if alfa >= 1.e-12:
				x = np.maximum(bl,np.minimum(bu,x+alfa*d))

		direzioni[:,i_corr] = np.copy(d)

		if alfa >= 1.e-12:
			flag_fail[i_corr] = False
			f = fz
			num_fal = 0
		else:
			flag_fail[i_corr] = True
			if i_corr_fall == 0:
				fstop[i_corr]   = fz1
				fstop[2*i_corr] = fz2
				for j in range(n):
					xfstop[j,i_corr]   = z1[j]
					xfstop[j,2*i_corr] = z2[j]
				num_fal += 1

		ni += 1
		z = np.copy(x)
		if i_corr < n-1:
			i_corr += 1
		else:
			if np.max(alfa_d) <= soglia and n > 1:
				if tipo_direzione == 0:
					fmin = fstop[0]
					fmax = fstop[0]
					imin = 0
					imax = 0
					doldalfamin = alfa_d[0]
					doldalfamax = alfa_d[0]
					iminalfa = 0
					imaxalfa = 0
					for i in range(1,n):
						if alfa_d[i] < doldalfamin:
							doldalfamin = alfa_d[i]
							iminalfa = i
						if alfa_d[i] > doldalfamax:
							doldalfamax = alfa_d[i]
							imaxalfa = i
					rapalfa = 3.0
					if doldalfamax/doldalfamin > rapalfa:
						for i in range(n):
							d1[i] = dconv[i]
						dnr = np.sqrt(np.double(n))
					else:
						for i in range(1,2*n):
							if fstop[i] < fmin:
								fmin = fstop[i]
								imin = i
							if fstop[i] > fmax:
								fmax = fstop[i]
								imax = i

						dnr = np.double(0.0)
						doldalfamedio = (doldalfamin+doldalfamax)/2.0
						for i in range(n):
							d1[i] = xfstop[i,imin]-xfstop[i,imax]
							dnr += d1[i]*d1[i]
						dnr = np.sqrt(dnr)
						if dnr <= 1.e-24:
							for i in range(n):
								d1[i] = dconv[i]
							dnr = np.sqrt(np.double(n))

					direzioni = gen_base(d1)
					direzioni = gram_schmidt(direzioni)
					tipo_direzione = 1

					alfa_coord = np.copy(alfa_d)
					if doldalfamax/doldalfamin > rapalfa:
						alfa_d = np.copy(alfa_diag)
					else:
						dnr = np.sum(alfa_d)/np.double(n)
						for i in range(n):
							alfa_d[i] = dnr

					if iprint >= 1:
						print("FINE DIR. COORDINATE")

				elif tipo_direzione == 1:
					direzioni = gen_base(d_dense)
					direzioni = gram_schmidt(direzioni)
					index_halton += 2*n
					if hschoice == 1:
						d_dense = np.asarray(Phalton[index_halton-1], dtype = np.double)
					else:
						d_dense, index_sobol = sobol_seq.i4_sobol(n, index_sobol)

					alfa_diag = np.copy(alfa_d)

					dnr = np.sum(alfa_d)/np.double(n)
					for i in range(n):
						alfa_d[i] = np.double(10.0)*dnr

					if iprint >= 1:
						print("FINE DIR. N+1")
					tipo_direzione = 2

				elif tipo_direzione == 2:

					direzioni = np.zeros((n,n))
					for i in range(n):
						direzioni[i,i] = np.double(1.0)

					tipo_direzione = 0

					alfa_dense = np.copy(alfa_d)
					alfa_d = np.copy(alfa_coord)

					if iprint >= 1:
						print("FINE DIR. DENSA")

				i_corr = 0

			i_corr = 0

		#################################################
		# EPSILON UPDATE FOR CONTRAINTS HANDLING
		#################################################
		if m > 0:
			fob,constr = funct(x)
			viol       = np.maximum(0.0,np.max(constr))
			#print(viol)
		else:
			fob = funct(x)
			viol= 0.0
		nf += 1

		if (viol > 0):
			cambio_eps = False
			maxeps = np.max(eps)

			for i in range(m):
				if(eps[i]*constr[i] > np.max(alfa_d)):
					eps[i] = 1.e-2*eps[i]
					if(iprint >= 1):
						print('**************************************')
						print('*********** aggiorno eps(',i,')=',eps[i],' *************')
						print('**************************************')
					cambio_eps = True

			if cambio_eps:
				fob,constr = funct(x)
				viol       = np.maximum(0.0,np.max(constr))
				f          = fob + np.maximum(0.0,np.max(constr/eps))
				nf += 1

				for i in range(n):
					alfa_d[i] = np.maximum(1.e-3,np.minimum(1.0,np.abs(x[i])))
					if(iprint >= 1):
						print(' alfainiz(',i,')=',alfa_d[i])

	return x, f

def gen_base(d):
	n = len(d)
	ind = np.argmax(np.abs(d))
	H = np.zeros((n,n))
	H[:,0] = d
	#print(ind)
	for i in range(1,ind+1):
		H[i-1,i] = np.double(1.0)
	for i in range(ind+1,n):
		H[i,i] = np.double(1.0)

	return H

def gram_schmidt(H):
	(n,n) = H.shape
	for i in range(1,n):
		proj = np.double(0.0)
		for j in range(i):
			proj += (np.dot(H[:,i],H[:,j])/np.dot(H[:,j],H[:,j]))*H[:,j]
		H[:,i] -= proj

	for i in range(n):
		H[:,i] = H[:,i]/np.linalg.norm(H[:,i])

	return H

def stop(alfa_d,alfa_max,nf,ni,tol,nf_max):
	istop = 0
	if alfa_max <= tol:
		istop = 1

	if nf > nf_max:
		istop = 2

	if ni > nf_max:
		istop = 3

	return istop

def linesearchbox_cont(funct,m,eps,x,f,d,alfa_d,z,j,alfa_max,bl,bu,nf,iprint):
	gamma = np.double(1.e-6)
	delta = np.double(0.5)
	delta1= np.double(0.5)
	ifront= 0
	i_corr_fall = 0
	n = len(x)

	z1 = np.zeros(n)
	z2 = np.zeros(n)
	fz1 = 0.0
	fz2 = 0.0

	if iprint >= 1:
		print("variabile continua j =",j,"    d(j) =",d[j]," alfa=",alfa_d[j])

	if np.abs(alfa_d[j]) <= 1.e-3*np.minimum(1.0,alfa_max):
		alfa = np.double(0.0)
		fz = 0.0
		if iprint >= 1:
			print("  alfa piccolo")
			print(" alfa_d(j)=",alfa_d[j],"    alfamax=",alfa_max)
		return alfa,fz, z1, fz1, z2, fz2, i_corr_fall, nf

	ifront = 0

	for ielle in [1,2]:
		if d[j] > 0.0:
			if alfa_d[j] - (bu[j]-x[j]) < -1.e-6:
				alfa = np.maximum(1.e-24,alfa_d[j])
			else:
				alfa = bu[j]-x[j]
				ifront = 1
				if iprint >= 1:
					print(" punto espan. sulla front. *")
		else:
			if alfa_d[j] - (x[j]-bl[j]) < -1.e-6:
				alfa = np.maximum(1.e-24,alfa_d[j])
			else:
				alfa = x[j]-bl[j]
				ifront = 1
				if iprint >= 1:
					print(" punto espan. sulla front. *")

		if np.abs(alfa) <= 1.e-3*np.minimum(1.0,alfa_max):
			d[j] = -d[j]
			i_corr_fall += 1
			ifront = 0
			if iprint >= 1:
				print(" direzione opposta per alfa piccolo")
				print(" j =",j,"    d(j) =",d[j])
				print(' alfa=',alfa,'    alfamax=',alfa_max)
			alfa = np.double(0.0)
		else:
			alfaex = alfa
			z[j] = x[j] + alfa*d[j]
			if m > 0:
				fob, constr = funct(z)
				fz = fob + np.maximum(0.0,np.max(constr/eps))
			else:
				fz = funct(z)
			nf    += 1

			if ielle == 1:
				z1 = np.copy(z)
				fz1 = fz
			else:
				z2 = np.copy(z)
				fz2 = fz

			if iprint >= 1:
				print(' fz =',fz,'   alfa =',alfa)
			if iprint >= 2:
				for i in range(n):
					print(' z(',i,')=',z[i])

			fpar = f - gamma * alfa**2
			if fz < fpar:
				while True:
					if ifront == 1:
						if iprint >= 1:
							print(' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa)

						alfa_d[j] = delta*alfa

						return alfa, fz, z1, fz1, z2, fz2, i_corr_fall, nf

					if d[j] > 0.0:
						if alfa/delta1 - (bu[j]-x[j]) < -1.e-6:
							alfaex = alfa/delta1
						else:
							alfaex = bu[j]-x[j]
							ifront = 1
							if iprint >= 1:
								print(' punto espan. sulla front.')
					else:
						if alfa/delta1 - (x[j]-bl[j]) < -1.e-6:
							alfaex = alfa/delta1
						else:
							alfaex = x[j]-bl[j]
							ifront = 1
							if iprint >= 1:
								print(' punto espan. sulla front.')

					z[j] = x[j] + alfaex*d[j]
					if m > 0:
						fob, constr = funct(z)
						fzdelta = fob + np.maximum(0.0,np.max(constr/eps))
					else:
						fzdelta = funct(z)
					nf    += 1

					if iprint >= 1:
						print(' fzex=',fzdelta,'  alfaex=',alfaex)
					if iprint >= 2:
						for i in range(n):
							print(' z(',i,')=',z[i])

					fpar = f - gamma * alfaex**2
					if fzdelta < fpar:
						fz = fzdelta
						alfa = alfaex
					else:
						alfa_d[j] = delta*alfa
						if iprint >= 1:
							print(' accetta punto fz =',fz,'   alfa =',alfa)
						return alfa, fz, z1, fz1, z2, fz2, i_corr_fall, nf

			else:
				d[j] = -d[j]
				ifront = 0
				if iprint >= 1:
					print(' direzione opposta')
					print(' j =',j,'    d(j) =',d[j])

	if not i_corr_fall==2:
		alfa_d[j] = delta*alfa_d[j]

	alfa = 0.0

	if iprint >= 1:
		print(' fallimento direzione')

	return alfa, fz, z1, fz1, z2, fz2, i_corr_fall, nf


def linesearchbox_dense(funct,m,eps,x,f,d,alfa_d,alfa_max,bl,bu,nf,iprint):
	n = len(x)
	gamma = np.double(1.e-6)
	delta = np.double(0.5)
	delta1= np.double(0.5)
	ifront= 0

	if iprint >= 1:
		print("direzione halton, alfa=",alfa_d)

	for ielle in [1, 2]:
		alfa   = alfa_d
		alfaex = alfa
		z      = x + alfa*d
		z      = np.maximum(bl,np.minimum(bu,z))
		if m > 0:
			fob, constr = funct(z)
			fz = fob + np.maximum(0.0,np.max(constr/eps))
		else:
			fz = funct(z)
		nf    += 1

		if iprint >= 1:
			print(" fz =",fz,"   alfa =",alfa)
		if iprint >= 2:
			for i in range(n):
				print(" z(",i,")=",z[i])

		fpar = f - gamma*alfa**2
		if fz < fpar:
			while True:
				alfaex = alfa/delta1
				z      = x + alfaex*d
				z      = np.maximum(bl,np.minimum(bu,z))
				if m > 0:
					fob, constr = funct(z)
					fzdelta = fob + np.maximum(0.0,np.max(constr/eps))
				else:
					fzdelta = funct(z)
				nf    += 1

				if iprint >= 1:
					print(" fzex=",fzdelta,"   alfaex=",alfaex)
				if iprint >= 2:
					for i in range(n):
						print(" z(",i,")=",z[i])

				fpar = f - gamma*alfaex**2
				if fzdelta < fpar:
					fz   = fzdelta
					alfa = alfaex
				else:
					alfa_d = alfa
					if iprint >= 1:
						print(" denso: accetta punto fz =",fz,"   alfa =",alfa)
					return alfa, alfa_d, fz, d, nf
		else:
			d      = -d
			ifront =  0

			if iprint >= 1:
				print("denso:  direzione opposta")

	alfa_d = delta*alfa_d
	alfa   = np.double(0.0)

	if iprint >= 1:
		print("denso: fallimento direzione")

	return alfa, alfa_d, fz, d, nf
