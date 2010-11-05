# ------------------------------------------------------------------
# Funcion Delta dirac n=2,3,4      Agosto 24 - 2010
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------

from  numpy import *

def phi_2(r):
    if (0 <= abs(r) <= 1):
        return (1.0-abs(r))
    if ( 1 <= abs(r)):
        return 0.0

def phi_3(r):
    if(0 <= abs(r) <= (1./2.)):
        return ((1./3.)*(1+sqrt(1-3*r**2)))
    if((1./2.) <= abs(r) <= (3./2.)):
        return ((1./6.)*(5-3+abs(r)-sqrt(-2+6*abs(r)-3*r**2)))
    if((3./2.) <= abs(r)):
        return 0.0

def phi_4(r):
    if(0 <= abs(r) <= 1):
        return ((1./8.)*(3-2*abs(r)+sqrt(1+4*abs(r)-4*r**2)))
    if(1 <= abs(r) <= 2):
        return ((1./8.)*(5-2*abs(r)-sqrt(-7+12*abs(r)-4*r**2)))
    if(2 <= abs(r)):
        return 0.0

def main():
    
    # Crear una grafica de cada funcion de impulso con diferente soporte
    x = arange(-2.0,2.0,0.05)
    y_2 = arange(2,-2,1)
    y_3 = arange(2,-2,1)
    y_4 = arange(2,-2,1)
    for i in range(0,x.shape[0]):
        y_2 = append(y_2,phi_2(x[i]))
        y_3 = append(y_3,phi_3(x[i]))
        y_4 = append(y_4,phi_4(x[i]))
    plot(x, y_2, 'bo-')
    plot(x, y_3, 'r+-')
    plot(x, y_4, 'yx-')
    
if __name__ == '__main__':
    main()