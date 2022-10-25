import math
import sympy as sym


##Using NIST Standard temperature and pressure: 293.15 K, 101325 Pa
stprho = 997.77 # density in kg . m^(-3)
stpmu = 1.0005E-3 # Dynamic viscosity in pa . s

#define symbols for symbolic derivative
mdot, rho, Q, diam, mu, eps, re, length, ff, k, pump = sym.symbols('mdot rho Q diam mu eps re length ff k pump', real=True)
def calcQ(mdot,rho=stprho):
    Q = mdot/rho
    return Q

def calcmdot(Q,rho=stprho):
    mdot = Q*rho
    return mdot

def calcRe(diam,mdot,mu=stpmu):
    re = (4 * abs(mdot)) / (math.pi * diam * mu)
    return re

def calcff(eps,diam,re):
    if re <= 64:
        return 1
    if re > 64 and re < 2300:
        return 64/re
    else:
        return 0.25 * math.log10(((eps/(diam*1000))/3.715) + 15/re)**-2

def calck(length,diam,ff ):
    k = (length * 8 * ff) / (9.81 * (math.pi ** 2) * diam ** 5)
    return k

def calchl(k,Q,pump=0):
    hl = k * Q**2 - pump
    return hl

def calcdhl(k,Q,pump=0):
    dhl = 2 * k * Q
    return dhl

def fullcalchl(Q,length,diam,eps,pump=0,rho=stprho,mu=stpmu):
    mdot = Q * rho
    re = (4 * abs(mdot)) / (math.pi * diam * mu)
    if (re <= 64):
        ff=1
    elif re > 64 and re <= 2300:
        ff=64/re
    elif re > 2300:
        ff=0.25 * math.log10(((eps/(diam*1000))/3.715) + 15/re)**-2
    else:
        raise ValueError('reynolds number out of bounds')

    k = (length * 8 * ff) / (9.81 * (math.pi ** 2) * diam ** 5)
    hl = k * Q**2 - pump
    return hl

def funcdhldQ(Q,length,diam,eps,pump=0):
    return sym.diff(fullcalchl(Q,length,diam,eps,pump),Q)

def dhldQ(Q_in,length_in,diam_in,eps_in,pump_in=0):
    return funcdhldQ(Q,length,diam,eps,pump).evalf(subs={Q: Q_in, length: length_in,diam: diam_in,eps: eps_in,pump: pump_in})


def numerical_dhldQ(Q,length,diam,eps,pump=0,rho=stprho,mu=stpmu,convcrit=1e-6):
    dQ = 10
    conv = 1
    d1=0
    while conv > convcrit:
        dQ = dQ * 0.1
        d1 = ( fullcalchl(Q+dQ,length,diam,eps,pump,rho,mu) - fullcalchl(Q,length,diam,eps,pump,rho,mu) ) / dQ
        d2 = (fullcalchl(Q + (dQ*0.1), length, diam, eps, pump, rho, mu) - fullcalchl(Q, length, diam, eps, pump, rho,
                                                                       mu)) / (dQ*0.1)
        conv = abs(d2-d1)
    return d1
