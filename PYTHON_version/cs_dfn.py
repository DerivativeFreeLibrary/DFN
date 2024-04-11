#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 18:05:46 2020

@author: giampo
"""

import numpy as np
import ghalton
import sobol_seq

class cs_dfn:

    x = []
    xint = []
    F = []
    eps = []
    funct = 0
    functpen = 0
    V = 0
    J = 0
    n = 0
    m = 0
    iprint = 0
    istop = 0
    nf_max = 0
    index_halton = 1000
    index_sobol  = 10000
    hschoice = 2
    num_fal = 0
    soglia = 1.e-3
    tol = 1.0
    print_format = '(c)| %5d | %5d | %5d | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %5d       |   '

    flag_fail = []
    fstop = []
    xfstop= []

    sequencer = []
    Phalton   = []

    alfa_d = []
    alfa_dense = []
    alfa_diag   = []
    alfa_coord  = []
    alfa_max    = np.inf

    d_dense = []
    d1 = []
    direzioni = []
    d = []

    i_corr  = 0
    i_dense = 0
    j_dense = 0
    tipo_direzione = 0

    f = 0
    ni = 0
    fstop   = []
    xfstop  = []
    z = []
    z1= []
    z2= []
    bl = []
    bu = []
    fz1 = 0
    fz2 = 0

    def __init__(self,functpen,f,J,V,eps,xtot,nf_max,tol,bl,bu,iprint):
        """
        Parameters
        ----------
        n : intero
            numero di var. CONTINUE
        nint : intero
            numero di var. DISCRETE
        xtot : array
            vettore completo di variabili:
                prima nint DISCRETE
                poi   n    CONTINUE
        nf_max : intero
            DESCRIPTION.
        iprint : intero
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.n = V.ncont
        n = self.n
        self.nint = V.nint
        self.m = V.m[J]
        self.x = np.copy(xtot[:n])
        self.xint = np.copy(xtot[n:])
        self.z = np.copy(self.x)
        self.bl = np.copy(bl)
        self.bu = np.copy(bu)
        print_format = '(c)| %05d | %05d | %05d | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %5d       |   '
        self.iprint = iprint
        self.tol = tol
        self.nf_max = nf_max
        self.functpen = functpen
        self.f = f
        self.V = V
        self.J = J
        self.eps = eps
        self.flag_fail = [False for i in range(n)]
        self.fstop = np.zeros(2*n)
        self.xfstop= np.zeros((n,2*n))

        self.sequencer = ghalton.Halton(n)
        self.Phalton   = self.sequencer.get(nf_max + 1000)

        self.alfa_d = np.zeros(n)
        print(self.alfa_d.shape)
        self.alfa_dense = np.zeros(n)
        for i in range(n):
            self.alfa_d[i]   = np.maximum(np.double(self.soglia*10.0),np.minimum(np.double(1.0),np.abs(self.x[i])))
            self.alfa_d[i]   = (bu[i]-bl[i])/2
            for j in range(n):
                self.alfa_dense[j] += self.alfa_d[i]
            if self.iprint >= 2:
                print("alfainiz(%d)=%f" % (i,self.alfa_d[i]))

        self.alfa_dense /= np.double(n)
        self.alfa_diag   = np.copy(self.alfa_d)
        self.alfa_coord  = np.copy(self.alfa_d)
        self.alfa_max    = np.max(self.alfa_d)

        if n > 1:
            if self.hschoice == 1:
                self.d_dense = np.asarray(self.Phalton[self.index_halton-1], dtype = np.double)
            else:
                self.d_dense, self.index_sobol = sobol_seq.i4_sobol(n, self.index_sobol)

        self.d1 = np.zeros(n)
        self.direzioni = np.zeros((n,n))
        self.d = np.ones(n)
        for i in range(n):
            self.direzioni[i,i] = np.double(1.0)

        self.fstop   = np.zeros(2*n+1)
        self.xfstop  = np.zeros((n,2*n+1))
        for j in range(2*n+1):
            self.fstop[j] = self.f
        for i in range(n):
            for j in range(2*n+1):
                self.xfstop[i,j] = self.x[i]

    def viol_constr(self,x,J,V):
        xtot = np.concatenate((x,self.xint),axis=0)
        if J == 'a':
            gloc = V.fconstr_a(xtot)
        elif J == 'b':
            gloc = V.fconstr_b(xtot)
        elif J == 'c':
            gloc = V.fconstr_c(xtot)
        elif J == 'd':
            gloc = V.fconstr_d(xtot)
        elif J == 'e':
            gloc = V.fconstr_e(xtot)
        elif J == 'f':
            gloc = V.fconstr_f(xtot)
        else:
            gloc = V.fconstr_z(xtot)

        viol = np.maximum(0.0,np.max(gloc))

        return viol, gloc

    def func_cont(self,x,eps,J,V,CACHE):
        funct_value = self.functpen(np.concatenate((x,self.xint),axis=0),eps,J,V,CACHE)
        return funct_value

    def alfa_init(self,n,x):
        self.alfa_d = np.zeros(n)
        self.alfa_dense = np.zeros(n)
        for i in range(n):
            self.alfa_d[i]   = np.maximum(np.double(self.soglia*10.0),np.minimum(np.double(1.0),np.abs(self.x[i])))
            for j in range(n):
                self.alfa_dense[j] += self.alfa_d[i]
            if self.iprint >= 2:
                print("alfainiz(%d)=%f" % (i,self.alfa_d[i]))

        self.alfa_dense /= np.double(n)
        self.alfa_diag   = np.copy(self.alfa_d)
        self.alfa_coord  = np.copy(self.alfa_d)
        self.alfa_max    = np.max(self.alfa_d)

    def step(self,alg,xint,f,eps,V,CACHE,nf):
        self.f    = f
        self.xint = np.copy(xint)
        self.eps  = np.copy(eps)
        cambio_eps = False
        while True:
            if self.n > 1:
                self.alfa_max = np.max([np.max(self.alfa_coord),np.max(self.alfa_diag),np.max(self.alfa_dense)])
            else:
                self.alfa_max = np.max(self.alfa_d)

            istop = stop(self.alfa_d,self.alfa_max,nf,self.ni,self.tol,self.nf_max)

            self.alfa_max = np.max(self.alfa_d)

            if istop >= 1:
                self.eps = np.copy(eps)
                #self.f   = f
                if self.iprint >= 2:
                    print(self.alfa_d)
                    print(self.alfa_max)
                return nf, cambio_eps
                #return nf

            if self.i_corr == 0:
                self.dconv = np.zeros(self.n)
                for i in range(self.n):
                    self.dconv += -self.direzioni[:,i]

            if self.iprint >= 2:
                print("----------------------------------------------")
            #if self.iprint >= 0:
            #	print(" ni=%4d  nf=%5d   f=%12.5e   alfamax=%12.5e" % (self.ni,self.nf,self.f,self.alfa_max))
            if self.iprint >= 2:
                #if (self.m > 0):
                #	print(self.print_format %(self.ni, self.nf, 0, self.f, self.f, 0.0, 0.0, self.alfa_max, 0))
                #else:
                viol, g = self.viol_constr(self.x,self.J,self.V)
                print(self.print_format %(self.ni, nf, CACHE.hits, self.f, self.f, np.sum(np.where(g<0,0,g)), self.alfa_max, self.i_corr))
            if self.iprint >= 2:
                for i in range(self.n):
                    print(" x(",i,")=",self.x[i])

            self.d = np.copy(self.direzioni[:,self.i_corr])

            if self.tipo_direzione == 0:
                self.alfa, self.fz, self.z1, self.fz1,self.z2, self.fz2, self.i_corr_fall, nf = linesearchbox_cont(
                                    self.func_cont,self.x,self.f,self.d,
                                    self.alfa_d,self.z,self.i_corr,self.alfa_max,
                                    self.bl,self.bu,nf,self.iprint,eps,self.J,self.V,CACHE)

                if self.alfa >= 1.e-12:
                    self.x[self.i_corr] = self.x[self.i_corr]+self.alfa*self.d[self.i_corr]
            else:

                self.alfa, self.alfatilde, self.fz, self.d, nf = linesearchbox_dense(self.func_cont,self.x,self.f,self.d,
                                                                              self.alfa_d[self.i_corr],self.alfa_max,
                                                                              self.bl,self.bu,nf,self.iprint,
                                                                              eps,self.J,self.V,CACHE)
                self.alfa_d[self.i_corr]      = self.alfatilde

                if self.alfa >= 1.e-12:
                    self.x = np.maximum(self.bl,np.minimum(self.bu,self.x+self.alfa*self.d))

            self.direzioni[:,self.i_corr] = np.copy(self.d)

            if self.alfa >= 1.e-12:
                self.flag_fail[self.i_corr] = False
                self.f = self.fz
                self.num_fal = 0
            else:
                self.flag_fail[self.i_corr] = True
                if self.i_corr_fall == 0:
                    self.fstop[self.i_corr]   = self.fz1
                    self.fstop[2*self.i_corr] = self.fz2
                    for j in range(self.n):
                        self.xfstop[j,self.i_corr]   = self.z1[j]
                        self.xfstop[j,2*self.i_corr] = self.z2[j]
                    self.num_fal += 1

            self.ni += 1
            self.z = np.copy(self.x)

            if self.iprint >= 1:
                viol, g = self.viol_constr(self.x,self.J,self.V)
                print(self.print_format %(self.ni, nf, CACHE.hits, self.f, self.f, np.sum(np.where(g<0,0,g)), self.alfa_max, self.i_corr))

            if self.i_corr < self.n-1:
                self.i_corr += 1
            else:
                if alg == 'DFN_DFL' and np.max(self.alfa_d) <= self.soglia and self.n > 1:
                    if self.tipo_direzione == 0:
                        fmin = self.fstop[0]
                        fmax = self.fstop[0]
                        imin = 0
                        imax = 0
                        doldalfamin = self.alfa_d[0]
                        doldalfamax = self.alfa_d[0]
                        iminalfa = 0
                        imaxalfa = 0
                        for i in range(1,self.n):
                            if self.alfa_d[i] < doldalfamin:
                                doldalfamin = self.alfa_d[i]
                                iminalfa = i
                            if self.alfa_d[i] > doldalfamax:
                                doldalfamax = self.alfa_d[i]
                                imaxalfa = i
                        rapalfa = 3.0
                        if doldalfamax/doldalfamin > rapalfa:
                            for i in range(self.n):
                                self.d1[i] = self.dconv[i]
                            self.dnr = np.sqrt(np.double(self.n))
                        else:
                            for i in range(1,2*self.n):
                                    if self.fstop[i] < fmin:
                                        fmin = self.fstop[i]
                                        imin = i
                                    if self.fstop[i] > fmax:
                                        fmax = self.fstop[i]
                                        imax = i

                            self.dnr = np.double(0.0)
                            doldalfamedio = (doldalfamin+doldalfamax)/2.0
                            for i in range(self.n):
                                self.d1[i] = self.xfstop[i,imin]-self.xfstop[i,imax]
                                self.dnr += self.d1[i]*self.d1[i]
                            self.dnr = np.sqrt(self.dnr)
                            if self.dnr <= 1.e-24:
                                for i in range(self.n):
                                    self.d1[i] = self.dconv[i]
                                self.dnr = np.sqrt(np.double(self.n))

                        self.direzioni = gen_base(self.d1)
                        self.direzioni = gram_schmidt(self.direzioni)
                        self.tipo_direzione = 1

                        self.alfa_coord = np.copy(self.alfa_d)
                        if doldalfamax/doldalfamin > rapalfa:
                            self.alfa_d = np.copy(self.alfa_diag)
                        else:
                            self.dnr = np.sum(self.alfa_d)/np.double(self.n)
                            for i in range(self.n):
                                self.alfa_d[i] = self.dnr

                        if self.iprint >= 2:
                            print("FINE DIR. COORDINATE")

                    elif self.tipo_direzione == 1:
                        self.direzioni = gen_base(self.d_dense)
                        self.direzioni = gram_schmidt(self.direzioni)
                        self.index_halton += 2*self.n
                        if self.hschoice == 1:
                            self.d_dense = np.asarray(self.Phalton[self.index_halton-1], dtype = np.double)
                        else:
                            self.d_dense, self.index_sobol = sobol_seq.i4_sobol(self.n, self.index_sobol)

                        self.tipo_direzione = 2
                        self.alfa_diag = np.copy(self.alfa_d)

                        self.dnr = np.sum(self.alfa_d)/np.double(self.n)
                        for i in range(self.n):
                            self.alfa_d[i] = np.double(10.0)*self.dnr

                        if self.iprint >= 2:
                            print("FINE DIR. N+1")

                    elif self.tipo_direzione == 2:

                        self.direzioni = np.zeros((self.n,self.n))
                        for i in range(self.n):
                            self.direzioni[i,i] = np.double(1.0)

                        self.tipo_direzione = 0

                        self.alfa_dense = np.copy(self.alfa_d)
                        self.alfa_d = np.copy(self.alfa_coord)

                        if self.iprint >= 2:
                            print("FINE DIR. DENSA")

                    self.i_corr = 0
                    break

                self.i_corr = 0
                break

        #####################
        #check eps
        #####################
        cambio_eps = False
        if self.m > 0:
            viol, constr = self.viol_constr(self.x,self.J,self.V)
            #print(viol,constr)
            if viol > 0.0:
                maxeps = np.max(eps)
                for i in range(self.m):
                    if(eps[i]*constr[i] > 1.0*np.max(self.alfa_d)):
                        eps[i]=1.e-2*eps[i]
                        if self.iprint >= 0:
                            print('**************************************')
                            print('*********** aggiorno eps(',i,')=',eps[i],' *************')
                            print('**************************************')
                        cambio_eps = True

                if cambio_eps:
                    self.f = self.func_cont(self.x,eps,self.J,self.V,CACHE)
                    #nf += 1

                    for i in range(self.n):
                        self.alfa_d[i]   = np.maximum(np.double(self.soglia*100.0),self.alfa_d[i])
                        self.alfa_d[i]   = (self.bu[i]-self.bl[i])/2
                    #input()

        self.eps = np.copy(eps)
        #self.f   = f

        return nf, cambio_eps
        #return nf

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

def linesearchbox_cont(funct,x,f,d,alfa_d,z,j,alfa_max,bl,bu,nf,iprint,eps,J,V,CACHE):
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

    if iprint >= 2:
        print("variabile continua j =",j,"    d(j) =",d[j]," alfa=",alfa_d[j])

    if np.abs(alfa_d[j]) <= 1.e-3*np.minimum(1.0,alfa_max):
        alfa = np.double(0.0)
        fz = 0.0
        if iprint >= 2:
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
                if iprint >= 2:
                    print(" punto espan. sulla front. *")
        else:
            if alfa_d[j] - (x[j]-bl[j]) < -1.e-6:
                alfa = np.maximum(1.e-24,alfa_d[j])
            else:
                alfa = x[j]-bl[j]
                ifront = 1
                if iprint >= 2:
                    print(" punto espan. sulla front. *")

        if np.abs(alfa) <= 1.e-3*np.minimum(1.0,alfa_max):
            d[j] = -d[j]
            i_corr_fall += 1
            ifront = 0
            if iprint >= 2:
                print(" direzione opposta per alfa piccolo")
                print(" j =",j,"    d(j) =",d[j])
                print(' alfa=',alfa,'    alfamax=',alfa_max)
            alfa = np.double(0.0)
        else:
            alfaex = alfa
            z[j] = x[j] + alfa*d[j]
            fz = funct(z,eps,J,V,CACHE)
            nf += 1

            if ielle == 1:
                z1 = np.copy(z)
                fz1 = fz
            else:
                z2 = np.copy(z)
                fz2 = fz

            if iprint >= 2:
                print(' fz =',fz,'   alfa =',alfa)
            if iprint >= 3:
                for i in range(n):
                    print(' z(',i,')=',z[i])

            fpar = f - gamma * alfa**2
            if fz < fpar:
                while True:
                    if ifront == 1:
                        if iprint >= 2:
                            print(' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa)

                        alfa_d[j] = delta*alfa

                        return alfa, fz, z1, fz1, z2, fz2, i_corr_fall, nf

                    if d[j] > 0.0:
                        if alfa/delta1 - (bu[j]-x[j]) < -1.e-6:
                            alfaex = alfa/delta1
                        else:
                            alfaex = bu[j]-x[j]
                            ifront = 1
                            if iprint >= 2:
                                print(' punto espan. sulla front.')
                    else:
                        if alfa/delta1 - (x[j]-bl[j]) < -1.e-6:
                            alfaex = alfa/delta1
                        else:
                            alfaex = x[j]-bl[j]
                            ifront = 1
                            if iprint >= 2:
                                print(' punto espan. sulla front.')

                    z[j] = x[j] + alfaex*d[j]
                    fzdelta = funct(z,eps,J,V,CACHE)
                    nf += 1

                    if iprint >= 2:
                        print(' fzex=',fzdelta,'  alfaex=',alfaex)
                    if iprint >= 3:
                        for i in range(n):
                            print(' z(',i,')=',z[i])

                    fpar = f - gamma * alfaex**2
                    if fzdelta < fpar:
                        fz = fzdelta
                        alfa = alfaex
                    else:
                        alfa_d[j] = delta*alfa
                        if iprint >= 2:
                            print(' accetta punto fz =',fz,'   alfa =',alfa)
                        return alfa, fz, z1, fz1, z2, fz2, i_corr_fall, nf

            else:
                d[j] = -d[j]
                ifront = 0
                if iprint >= 2:
                    print(' direzione opposta')
                    print(' j =',j,'    d(j) =',d[j])

    if not i_corr_fall==2:
        alfa_d[j] = delta*alfa_d[j]

    alfa = 0.0

    if iprint >= 2:
        print(' fallimento direzione')

    return alfa, fz, z1, fz1, z2, fz2, i_corr_fall, nf


def linesearchbox_dense(funct,x,f,d,alfa_d,alfa_max,bl,bu,nf,iprint,eps,J,V,CACHE):
    gamma = np.double(1.e-6)
    delta = np.double(0.5)
    delta1= np.double(0.5)
    ifront= 0

    if iprint >= 2:
        print("direzione halton, alfa=",alfa_d)

    for ielle in [1, 2]:
        alfa   = alfa_d
        alfaex = alfa
        z      = x + alfa*d
        z      = np.maximum(bl,np.minimum(bu,z))
        fz     = funct(z,eps,J,V,CACHE)
        nf    += 1

        if iprint >= 2:
            print(" fz =",fz,"   alfa =",alfa)
        if iprint >= 3:
            for i in range(n):
                print(" z(",i,")=",z[i])

        fpar = f - gamma*alfa**2
        if fz < fpar:
            while True:
                alfaex = alfa/delta1
                z      = x + alfaex*d
                z      = np.maximum(bl,np.minimum(bu,z))
                fzdelta= funct(z,eps,J,V,CACHE)
                nf    += 1

                if iprint >= 2:
                    print(" fzex=",fzdelta,"   alfaex=",alfaex)
                if iprint >= 3:
                    for i in range(n):
                        print(" z(",i,")=",z[i])

                fpar = f - gamma*alfaex**2
                if fzdelta < fpar:
                    fz   = fzdelta
                    alfa = alfaex
                else:
                    alfa_d = alfa
                    if iprint >= 2:
                        print(" denso: accetta punto fz =",fz,"   alfa =",alfa)
                    return alfa, alfa_d, fz, d, nf
        else:
            d      = -d
            ifront =  0

            if iprint >= 2:
                print("denso:  direzione opposta")

    alfa_d = delta*alfa_d
    alfa   = np.double(0.0)

    if iprint >= 2:
        print("denso: fallimento direzione")

    return alfa, alfa_d, fz, d, nf
