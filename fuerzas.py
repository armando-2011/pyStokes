# ------------------------------------------------------------------
# Stretching Triangular Membrane   Septiembre 12 - 2010
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------
# Input: pi, pj, pk -> Nodos del elemento NO deformado i1, i2, i3
#        Pi, Pj, Pk -> Nodos del elemento deformado    I1, I2, I3
# ------------------------------------------------------------------
# Output: F -> Vector fuerza por cada nodo en la base  i1, i2, i3
# ------------------------------------------------------------------

# Librerias utilizadas 
import numpy as np


# Verificacion del algoritmo de fuerzas en el plano
pi = np.array([0,0,0])
pj = np.array([2,0,0])
pk = np.array([1,1,0])

Pi = np.array([0,0,0])
Pj = np.array([2,0,0])
Pk = np.array([5,3,0])


# Coordenadas iniciales de los tres nodos [cm] revision Sep/14/2010
xi = pi[0]
yi = pi[1]
zi = pi[2]
  
xj = pj[0]
yj = pj[1]
zj = pj[2]
  
xk = pk[0]
yk = pk[1]
zk = pk[2]

# Coordenadas del elemento deformado [cm] revision Sep/14/2010
Xi = Pi[0]
Yi = Pi[1]
Zi = Pi[2]

Xj = Pj[0]
Yj = Pj[1]
Zj = Pj[2]

Xk = Pk[0]
Yk = Pk[1]
Zk = Pk[2]

# Desplazamientos de los nodos en coordenadas x,y,z revision Sep/14/2010
ui = Xi - xi
vi = Yi - yi
wi = Zi - zi
uj = Xj - xj
vj = Yj - yj
wj = Zj - zj
uk = Xk - xk
vk = Yk - yk
wk = Zk - zk

# Coeficientes de funciones de forma revision Sep/14/2010
ai = yj-yk 
bi = xk-xj
ci = xj*yk - xk*yj
Li = ai*xi + bi*yi + ci

aj = yk - yi
bj = xi - xk
cj = xk*yi - xi*yk
Lj = aj*xj + bj*yj + cj

ak = yi - yj
bk = xj - xi
ck = xi*yj - xj*yi
Lk = ak*xk + bk*yk + ck

# Derivadas parciales para calcular el vector [G] revision Sep/14/2010
dudx = ui*ai/Li + uj*aj/Lj + uk*ak/Lk
dudy = ui*bi/Li + uj*bj/Lj + uk*bk/Lk
dvdx = vi*ai/Li + vj*aj/Lj + vk*ak/Lk
dvdy = vi*bi/Li + vj*bj/Lj + vk*bk/Lk

# Componentes del vector [G] revision Sep/14/2010
g11 = (1+dudx)**2 + (dvdx)**2
g12 = (1+dudx)*(dudy) + (1+dvdy)*(dvdx)
g21 = g12
g22 = (1+dvdy)**2 + (dudy)**2
G = np.array([[g11,g12],[g21,g22]])

# Calculo de lambda1 y lambda2 revision Sep/14/2010
l1 = np.sqrt((g11 + g22  + np.sqrt((g11-g22)**2 + 4*g12**2))/2)
l2 = np.sqrt((g11 + g22  - np.sqrt((g11-g22)**2 + 4*g12**2))/2)

# Derivadas de la funcion Strain Energy respecto a lambda 1 y lambda 2 
# Modelo de energia Skalak 1973 revision Sep/14/2010
B = 0.005	# [dyn/cm]
C = 100.0	# [dyn/cm]
I1 = (l1**2) + (l2**2) - 2
I2 = (l1**2)*(l2**2)-1
dI1dl1 = 2*l1 
dI1dl2 = 2*l2
dI2dl1 = 2*l1*(l2**2)
dI2dl2 = 2*l2*(l1**2)
dwdl1 = (B/4)*(I1*dI1dl1 + dI1dl1 - dI2dl1) + (C/4)*(I2)*(dI2dl1)
dwdl2 = (B/4)*(I1*dI1dl2 + dI1dl2 - dI2dl2) + (C/4)*(I2)*(dI2dl2)

# Calculo de diferenciales sobre l1 y l2 respecto a desplazamientos de nodos
# 1. Derivadas de [G] respecto a desplazamiento de nodos revision Sep/14/2010
dg11dui = 2*(1+ dudx)*(ai/Li)
dg11duj = 2*(1+ dudx)*(aj/Lj)
dg11duk = 2*(1+ dudx)*(ak/Lk)
dg11dvi = 2*dvdx*(ai/Li)
dg11dvj = 2*dvdx*(aj/Lj)
dg11dvk = 2*dvdx*(ak/Lk)

dg12dui = (1+dudx)*(bi/Li) + (ai/Li)*(dudy)
dg12duj = (1+dudx)*(bj/Lj) + (aj/Lj)*(dudy)
dg12duk = (1+dudx)*(bk/Lk) + (ak/Lk)*(dudy)
dg12dvi = (1+dvdy)*(ai/Li) + (bi/Li)*(dvdx)
dg12dvj = (1+dvdy)*(aj/Lj) + (bj/Lj)*(dvdx)
dg12dvk = (1+dvdy)*(ak/Lk) + (bk/Lk)*(dvdx)

dg22dui = 2*dudy*(bi/Li) 
dg22duj = 2*dudy*(bj/Lj) 
dg22duk = 2*dudy*(bk/Lk) 
dg22dvi = 2*(1+dvdy)*(bi/Li)
dg22dvj = 2*(1+dvdy)*(bj/Lj)
dg22dvk = 2*(1+dvdy)*(bk/Lk)

# 2. Calculo de las derivadas de lambda 1 y lambda 2 respecto desplazamientos
# nodales revision Sep/14/2010
if(g12!=0):
    a = (1./2.)*((g11-g22)**2 + 4*g12**2)**(-0.5)
    b = 2*(g11-g22)*(dg11dui - dg22dui) + 8*g12*dg12dui
    dl1dui = ((dg11dui) + (dg22dui) + (a*b))/(4*l1)
    dl2dui = ((dg11dui) + (dg22dui) - (a*b))/(4*l2)
    
    b = 2*(g11-g22)*(dg11duj - dg22duj) + 8*g12*dg12duj
    dl1duj = ((dg11duj) + (dg22duj) + (a*b))/(4*l1)
    dl2duj = ((dg11duj) + (dg22duj) - (a*b))/(4*l2)
    
    b = 2*(g11-g22)*(dg11duk - dg22duk) + 8*g12*dg12duk
    dl1duk = ((dg11duk) + (dg22duk) + (a*b))/(4*l1)
    dl2duk = ((dg11duk) + (dg22duk) - (a*b))/(4*l2)
    
    b = 2*(g11-g22)*(dg11dvi - dg22dvi) + 8*g12*dg12dvi
    dl1dvi = ((dg11dvi) + (dg22dvi) + (a*b))/(4*l1)
    dl2dvi = ((dg11dvi) + (dg22dvi) - (a*b))/(4*l2)
    
    b = 2*(g11-g22)*(dg11dvj - dg22dvj) + 8*g12*dg12dvj
    dl1dvj = ((dg11dvj) + (dg22dvj) + (a*b))/(4*l1)
    dl2dvj = ((dg11dvj) + (dg22dvj) - (a*b))/(4*l2)
    
    b = 2*(g11-g22)*(dg11dvk - dg22dvk) + 8*g12*dg12dvk
    dl1dvk = ((dg11dvk) + (dg22dvk) + (a*b))/(4*l1)
    dl2dvk = ((dg11dvk) + (dg22dvk) - (a*b))/(4*l2)
else:
    dl1dui = dg11dui/(2*l1)
    dl1duj = dg11duj/(2*l1)
    dl1duk = dg11duk/(2*l1)
    
    dl1dvi = dg11dvi/(2*l1)
    dl1dvj = dg11dvj/(2*l1)
    dl1dvk = dg11dvk/(2*l1)
    
    dl2dui = dg11dui/(2*l2)
    dl2duj = dg11duj/(2*l2)
    dl2duk = dg11duk/(2*l2)
    
    dl2dvi = dg11dvi/(2*l2)
    dl2dvj = dg11dvj/(2*l2)
    dl2dvk = dg11dvk/(2*l2)

# 3. Calculo de las derivadas de w respecto a los desplazamientos nodales
# revision Sep/14/2010
dwdui = dwdl1*dl1dui + dwdl2*dl2dui
dwdvi = dwdl1*dl1dvi + dwdl2*dl2dvi
dwduj = dwdl1*dl1duj + dwdl2*dl2duj
dwdvj = dwdl1*dl1dvj + dwdl2*dl2dvj
dwduk = dwdl1*dl1duk + dwdl2*dl2duk
dwdvk = dwdl1*dl1dvk + dwdl2*dl2dvk

# 4. Volumen del elemento revision Sep/14/2010
a0 = ((xj-xi)*(yk-yi) - (xk-xi)*(yj-yi))/.2
t0 = 0.5

# 5. Calculo de las componentes de fuerza revision Sep/14/2010
fxi = dwdui*a0*t0 
fyi = dwdvi*a0*t0 
fzi = 0.0
fxj = dwduj*a0*t0 
fyj = dwdvj*a0*t0 
fzj = 0.0
fxk = dwduk*a0*t0 
fyk = dwdvk*a0*t0 
fzk = 0.0

# 6. Ensamblar el vector de fuerzas a retornar por cada nodo
fi = np.array([fxi, fyi, fzi])
fj = np.array([fxj, fyj, fzj])
fk = np.array([fxk, fyk, fzk])
F = np.vstack((fi,fj))
F = np.vstack((F,fk))

# 7. Momentos en  cualquier punto 
vji = Pj-Pi
vki = Pk-Pi
M = np.cross(F[1],vji) - np.cross(F[2],vki)

print "Nodos del elemento inicial:", pi, pj, pk
print "Nodos del elemento final:", Pi, Pj, Pk
print "Tensor metrico:\n",G
print "Fuerza:\n", F
print "Equlibrio fuerzas:\n", (fxi+fxj+fxk),(fyi+fyj+fyk)
print "Momentos\n", M