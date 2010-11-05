# ------------------------------------------------------------------
# Stretching Triangular Membrane   Octubre 04 - 2010
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

def fuerzas(pi, pj, pk, Pi, Pj, Pk):
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
    B = 1.0	# [dyn/cm]
    C = 1.0	# [dyn/cm]
    I1 = (l1**2.0) + (l2**2) - 2.0
    I2 = (l1**2)*(l2**2)-1.0
    dI1dl1 = 2*l1 
    dI1dl2 = 2*l2
    dI2dl1 = 2*l1*(l2**2)
    dI2dl2 = 2*l2*(l1**2)
#    dwdl1 = (B/4)*(I1*dI1dl1 + dI1dl1 - dI2dl1) + (C/4)*(I2)*(dI2dl1)
#    dwdl2 = (B/4)*(I1*dI1dl2 + dI1dl2 - dI2dl2) + (C/4)*(I2)*(dI2dl2)
    ks = 0.022131199
    dwdl1 = (ks/12)*(2*I1*dI1dl1 + 2*dI1dl1 -2*dI2dl1) + (ks/6)*I2*dI2dl1
    dwdl2 = (ks/12)*(2*I1*dI1dl2 + 2*dI1dl2 -2*dI2dl2) + (ks/6)*I2*dI2dl2
    
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
    # Formulacion Rolling John Hopkins
    t0 = np.sqrt((g11-g22)**2 + 4.*g12*g12)
    t1 = ((g11 + g22 + t0))
    jl1 = np.sqrt(0.5*t1)
    jl2 = np.sqrt(0.5*(g11+g22-t0))
    t2 = dg11dui - dg22dui
    if(np.abs(t0) > 1.0e-3):
        t3 = 0.5/(t0*(g11-g22)*(dg11dui - dg22dui) + 8.*g12*dg12dui)
    else:
        t3=0.0
    dl1dui = 0.5*np.sqrt(0.5)/np.sqrt(t1)*(t2+t3)
    dl2dui = 0.5*np.sqrt(0.5)/np.sqrt(t1-2.*t0)*(t2-t3)
    
    t2 = dg11duj - dg22duj
    if(abs(t0) > 1.0e-03):
        t3  = + 0.5/t0*(2.*(g11-g22)*(dg11duj - dg22duj) + 8.*g12*dg12duj )
    else:
       t3 =0.0
    dl1duj = 0.5*np.sqrt(0.5)/np.sqrt(t1)*(t2+t3)
    dl2duj = 0.5*np.sqrt(0.5)/np.sqrt(t1-2.*t0)*(t2-t3)
    t2 = dg11duk + dg22duk
    if(abs(t0) > 1.0e-03):
       t3  = + 0.5/t0*(2.*(g11-g22)*(dg11duk - dg22duk) + 8.*g12*dg12duk )
    else:
       t3 = 0.0
    dl1duk = 0.5*np.sqrt(0.5)/np.sqrt(t1)*(t2+t3)
    dl2duk = 0.5*np.sqrt(0.5)/np.sqrt(t1-2.*t0)*(t2-t3)
    t2 = dg11dvi + dg22dvi
    if(abs(t0) > 1.0e-03):
       t3  = + 0.5/t0*(2.*(g11-g22)*(dg11dvi - dg22dvi) + 8.*g12*dg12dvi )
    else:
       t3 = 0.0
    dl1dvi = 0.5*np.sqrt(0.5)/np.sqrt(t1)*(t2+t3)
    dl2dvi = 0.5*np.sqrt(0.5)/np.sqrt(t1-2.*t0)*(t2-t3)
    t2 = dg11dvj + dg22dvj
    if(abs(t0) > 1.0e-03):
       t3  = + 0.5/t0*(2.*(g11-g22)*(dg11dvj - dg22dvj)+ 8.*g12*dg12dvj )
    else:
       t3 = 0.0
    dl1dvj = 0.5*np.sqrt(0.5)/np.sqrt(t1)*(t2+t3)
    dl2dvj = 0.5*np.sqrt(0.5)/np.sqrt(t1-2.*t0)*(t2-t3)
    t2 = dg11dvk + dg22dvk
    if(abs(t0) > 1.0e-03):
        t3  = + 0.5/t0*(2.*(g11-g22)*(dg11dvk - dg22dvk)+ 8.*g12*dg12dvk )
    else:
        t3 = 0.0
    dl1dvk = 0.5*np.sqrt(0.5)/np.sqrt(t1)*(t2+t3)
    dl2dvk = 0.5*np.sqrt(0.5)/np.sqrt(t1-2.*t0)*(t2-t3)
    

#    a = (1./2.)*((g11-g22)**2 + 4*g12**2)**(-1/2)
#    b = 2*(g11-g22)*(dg11dui - dg22dui) + 8*g12*dg12dui
#    dl1dui = ((dg11dui) + (dg22dui) + (a*b))/(4*l1)
#    dl2dui = ((dg11dui) + (dg22dui) - (a*b))/(4*l2)
#    
#    b = 2*(g11-g22)*(dg11duj - dg22duj) + 8*g12*dg12duj
#    dl1duj = ((dg11duj) + (dg22duj) + (a*b))/(4*l1)
#    dl2duj = ((dg11duj) + (dg22duj) - (a*b))/(4*l2)
#    
#    b = 2*(g11-g22)*(dg11duk - dg22duk) + 8*g12*dg12duk
#    dl1duk = ((dg11duk) + (dg22duk) + (a*b))/(4*l1)
#    dl2duk = ((dg11duk) + (dg22duk) - (a*b))/(4*l2)
#    
#    b = 2*(g11-g22)*(dg11dvi - dg22dvi) + 8*g12*dg12dvi
#    dl1dvi = ((dg11dvi) + (dg22dvi) + (a*b))/(4*l1)
#    dl2dvi = ((dg11dvi) + (dg22dvi) - (a*b))/(4*l2)
#    
#    b = 2*(g11-g22)*(dg11dvj - dg22dvj) + 8*g12*dg12dvj
#    dl1dvj = ((dg11dvj) + (dg22dvj) + (a*b))/(4*l1)
#    dl2dvj = ((dg11dvj) + (dg22dvj) - (a*b))/(4*l2)
#    
#    b = 2*(g11-g22)*(dg11dvk - dg22dvk) + 8*g12*dg12dvk
#    dl1dvk = ((dg11dvk) + (dg22dvk) + (a*b))/(4*l1)
#    dl2dvk = ((dg11dvk) + (dg22dvk) - (a*b))/(4*l2)
    
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
    t0 = 0.1
    
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
    
    # 6. Momentos en  cualquier punto 
    Mx = fyj*(Xj-Xi) + fyk*(Xk-Xi)
    My = fxj*(Yj-Yi) + fxk*(Yk-Yi)
    
    # 7. Ensamblar el vector de fuerzas a retornar por cada nodo
    fi = np.array([fxi, fyi, fzi])
    fj = np.array([fxj, fyj, fzj])
    fk = np.array([fxk, fyk, fzk])
    F = np.vstack((fi,fj))
    F = np.vstack((F,fk))
    return F

def main():
    
    # Comprobacipn del algoritmo de traslacion y rotacion de elementos
    print "---------------------------------------------------------\n"
    pi = np.array([0.0,0.0,0.0])
    pj = np.array([2.0,0.0,0.0])
    pk = np.array([1.0,1.0,0.0])
    Pi = np.array([0.0,0.0,0.0])
    Pj = np.array([2.0,0.0,0.0])
    Pk = np.array([1.5,2.0,0.0])
    print "Puntos del elemento inicial:", pi, pj, pk
    print "Puntos del elemento deformado:", Pi, Pj, Pk
    print "Fuerzas resultantes: \n", fuerzas(pi,pj,pk,Pi,Pj,Pk)
    print "---------------------------------------------------------\n"
    
if __name__ == '__main__':
    main()
