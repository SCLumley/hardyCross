import math

def calcQ(mdot,rho):
    Q = mdot/rho
    return Q

def calcmdot(Q,rho):
    mdot = Q*rho
    return mdot

def calcRe(mdot,diam,mu=1.0016E-3):
    re = (4 * abs(mdot)) / (math.pi * diam * mu)
    return re

def calcff(re,eps,diam):
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