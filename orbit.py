## A python code to find the orbit of a planet [2018 Nov 22]
## Email : pbera.phy@gmail.com
## d2v/dphi^2 = eta*v^2 - v + 1/(1-eps^2)
## r = b/v; b = semi-mejor axis, v = dimensionless parameter
###############################################
import numpy as np
import pylab as pl

CONST_PI = np.pi

Msun = 1.9885e30
G = 6.6740831e-11
b = 5.790905e10
c = 2.99792458e8
eps = 0.20563593

eta = 3.0*G*Msun/(b*c*c)

def func_d2v_dphi2(dv_dphi, v, phi):
	return eta*v*v - v + 1.0/(1.0-eps*eps)

def func_dv_dphi(dv_dphi, v, phi):
	return dv_dphi

Norbit = 2
orbit_steps = 1024*1024
Nsteps = Norbit*orbit_steps
dphi = 2.*CONST_PI/orbit_steps

Array_phi = [0]*(Nsteps+1)
Array_v = [0]*(Nsteps+1)
Array_dv_dphi = [0]*(Nsteps+1)

r = [0]*(Nsteps+1)
X = [0]*(Nsteps+1)
Y = [0]*(Nsteps+1)

for i in range(Nsteps+1):
	if (i == 0):
		Array_phi[0] = 0.0
		Array_v[0] = 1.0+eps
		Array_dv_dphi[0] = 0.0

	else:
		Array_phi[i] = Array_phi[i-1] + dphi

		vK1 = func_dv_dphi(Array_dv_dphi[i-1], Array_v[i-1], Array_phi[i-1])*dphi
		dv_dphiK1 = func_d2v_dphi2(Array_dv_dphi[i-1], Array_v[i-1], Array_phi[i-1])*dphi

		vK2 = func_dv_dphi(Array_dv_dphi[i-1]+0.5*dv_dphiK1, Array_v[i-1]+0.5*vK1, Array_phi[i-1]+0.5*dphi)*dphi
		dv_dphiK2 = func_d2v_dphi2(Array_dv_dphi[i-1]+0.5*dv_dphiK1, Array_v[i-1]+0.5*vK1, Array_phi[i-1]+0.5*dphi)*dphi

		vK3 = func_dv_dphi(Array_dv_dphi[i-1]+0.5*dv_dphiK2, Array_v[i-1]+0.5*vK2, Array_phi[i-1]+0.5*dphi)*dphi
		dv_dphiK3 = func_d2v_dphi2(Array_dv_dphi[i-1]+0.5*dv_dphiK2, Array_v[i-1]+0.5*vK2, Array_phi[i-1]+0.5*dphi)*dphi
	
		vK4 = func_dv_dphi(Array_dv_dphi[i-1]+dv_dphiK3, Array_v[i-1]+vK3, Array_phi[i-1]+dphi)*dphi
		dv_dphiK4 = func_d2v_dphi2(Array_dv_dphi[i-1]+dv_dphiK3, Array_v[i-1]+vK3, Array_phi[i-1]+dphi)*dphi

		Array_v[i] = Array_v[i-1] + (vK1 + 2*vK2 + 2*vK3 + vK4)/6.0
		Array_dv_dphi[i] = Array_dv_dphi[i-1] + (dv_dphiK1 + 2*dv_dphiK2 + 2*dv_dphiK3 + dv_dphiK4)/6.0
	
	r[i] = 1./Array_v[i]
	X[i] = r[i]*np.cos(Array_phi[i])
	Y[i] = r[i]*np.sin(Array_phi[i])
	if (i%(1*orbit_steps)==0):
		print Array_phi[i]/(2*CONST_PI), r[i], Array_dv_dphi[i]

print Array_phi[-1]/(2*CONST_PI), r[-1], max(Array_dv_dphi)
#pl.plot(Array_phi, Array_v)
pl.plot(Array_phi, Array_dv_dphi)
pl.plot(Array_phi, max(Array_dv_dphi)*np.sin(Array_phi))
#pl.plot(X,Y, 'x')
pl.show()

## precession angle from the phase difference of Array_dv_dphi; calculated from the sinusoildal assumption around the 'zero'
print ' Mercury\'s precession (arc-sec per Earth century): ', (2*CONST_PI/Array_phi[-1])*Array_dv_dphi[-1]/max(Array_dv_dphi)*(365.25*100/88)*(3600*180/np.pi)