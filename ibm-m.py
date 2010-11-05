# ------------------------------------------------------------------
# Immersed Boundary Method         Septiembre 16 - 2010
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------

import numpy as np

#-------------------------------------------------------------------
# Function: phi_x, Funciones phi con soporte x para calcular dDirac
# ------------------------------------------------------------------
def phi_2(r):
    if (0 <= abs(r) <= 1):
        return (1.0-abs(r))
    if ( 1 <= abs(r)):
        return 0.0

def phi_3(r):
    if(0 <= abs(r) <= (1./2.)):
        return ((1./3.)*(1+np.sqrt(1-3*r**2)))
    if((1./2.) <= np.abs(r) <= (3./2.)):
        return ((1./6.)*(5-3+abs(r)-np.sqrt(-2+6*abs(r)-3*r**2)))
    if((3./2.) <= abs(r)):
        return 0.0

def phi_4(r):
    if(0 <= abs(r) <= 1):
        return ((1./8.)*(3-2*abs(r)+np.sqrt(1+4*abs(r)-4*r**2)))
    if(1 <= abs(r) <= 2):
        return ((1./8.)*(5-2*abs(r)-np.sqrt(-7+12*abs(r)-4*r**2)))
    if(2 <= abs(r)):
        return 0.0
        
#-------------------------------------------------------------------
# Function: dirac_x, Funciones dirac_x con soporte x 
# ------------------------------------------------------------------
def dirac_2(x):
    d = phi_2(x[0])*phi_2(x[1])*phi_2(x[2])
    return d

def dirac_3(x):
    d = phi_3(x[0])*phi_3(x[1])*phi_3(x[2])
    return d

def dirac_4(x):
    d = phi_4(x[0])*phi_4(x[1])*phi_4(x[2])
    return d

#-------------------------------------------------------------------
# Function: spread, propaga la fuerza de cada nodo en la membrana 
#           hacia la malla de fluido Euleriana
# ------------------------------------------------------------------
# Input: F     -> Estructura con vector fuerza en nodos de membrana
#        FV    -> Geometria de la membrana deformada
#        x,y,z -> Dimensiones de la malla euleriana 
#        dx    -> Espaciamiento uniforme de la malla 
# ------------------------------------------------------------------
# Output: f -> Estructura de fuerza Euleriana para calcular feq
#              la fuerza esta en componentes en la base x,y,z
# ------------------------------------------------------------------

def spread(F,FV,fluido,x,y,z,dx):
    
    # Estructura para la densidad de fuerza en cada nodo de fluido
    # por componentes
    f = np.zeros((x,y,z,3))
    for l in xrange(0,FV.vertices.shape[0]):
        A = int(np.rint(FV.vertices[l][0]/dx))
        B = int(np.rint(FV.vertices[l][1]/dx))
        C = int(np.rint(FV.vertices[l][2]/dx))
        i = A-3
        x = A+3
        j = B-3
        y = B+3
        k = C-3
        z = C+3
        for i in xrange(0,x):
            for j in xrange(0,y):
                for k in xrange(0,z):
                    d = dirac_4((fluido[i,j,k] - FV.vertices[l])/dx)
                    f[i,j,k,:] += F[l]*d
    return f
    

#-------------------------------------------------------------------
# Function: interpolation, propaga la velocidad del fluido hacia
#           los nodos de la membrana            
# ------------------------------------------------------------------
# Input: u_x   -> Velocidades del fluido por componentes
#        FV    -> Geometria de la membrana
#        x,y,z -> Dimensiones de la malla euleriana 
#        dx    -> Espaciamiento uniforme de la malla 
# ------------------------------------------------------------------
# Output: u -> Estructura de velocidad para cada nodo en la malla
#              lagrangiana
# ------------------------------------------------------------------
def interpolation(FV,fluido, x, y, z, dx, u_x, u_y, u_z):
    
    # Recorrido sobre los nodos de la membrana INTERPOLATION
    # Ajustado al rectangulo de interes
    u = np.zeros(FV.vertices.shape)
    for l in xrange(0,FV.vertices.shape[0]):
        A = int(np.rint(FV.vertices[l][0]/dx))
        B = int(np.rint(FV.vertices[l][1]/dx))
        C = int(np.rint(FV.vertices[l][2]/dx))
        i = A-3
        x = A+3
        j = B-3
        y = B+3
        k = C-3
        z = C+3
        for i in xrange(0,x):
            for j in xrange(0,y):
                for k in xrange(0,z):
                    d = dirac_4((fluido[i,j,k] - FV.vertices[l])/dx)
                    ux = u_x[i,j,k]*d
                    uy = u_y[i,j,k]*d
                    uz = u_z[i,j,k]*d
                    u[l,:] += (ux, uy, uz)
    return u
    
def main():
    # Crear una grafica de cada funcion de impulso con diferente soporte
    x = np.arange(-2.0,2.0,0.05)
    y_2 = np.arange(2,-2,1)
    y_3 = np.arange(2,-2,1)
    y_4 = np.arange(2,-2,1)
    for i in range(0,x.shape[0]):
        y_2 = np.append(y_2,phi_2(x[i]))
        y_3 = np.append(y_3,phi_3(x[i]))
        y_4 = np.append(y_4,phi_4(x[i]))
    print dirac_2([1,1.2,1.3])
    
if __name__ == '__main__':
    main()
    
    
    
    
    
    
    