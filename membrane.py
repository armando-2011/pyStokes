# ------------------------------------------------------------------
# Skalak Membrane Model            Agosto 22 - 2010
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------

from numpy import *
from sphere_mesh import *

def main():

    ks = 1     # Surface elastic shear modulus (Shear)
    ka = 1     # Area dilatation modulus       (Dilatation)
    
    FVe = triang(shape = 'ico', maxlevel = 3, r = 1.0 ,winding = 1)
    FVd = FVe
    WS = 0
    
    # Calcula la energia elastica para cada elemento
    for i in range(0,FV.faces.shape[0]):
        
        # Elemento en estado equilibrio
        elemento_0 = FVe.faces[i] 
        p1_0 = FVe.vertices[elemento_0[0]]
        p2_0 = FVe.vertices[elemento_0[1]]
        p3_0 = FVe.vertices[elemento_0[2]]
        l_0 = p3_0 - p1_0
        lp_0 = p2_0 -p1_0
        phi_0 = arccos((dot(l_0,lp_0)/(norm(l_0)*norm(lp_0))))
        
        # Elemento deformado 
        elemento = FVd.faces[i] 
        p1 = FVe.vertices[elemento[0]]
        p2 = FVe.vertices[elemento[1]]
        p3 = FVe.vertices[elemento[2]]
        l = p3 - p1
        lp = p2 - p1
        phi = arccos((dot(l,lp)/(norm(l)*norm(lp))))
        
        # Tomar norma de los vectores l, l_0, lp, lp_0
        l = norm(l)
        l_0 = norm(l_0)
        lp = norm(lp)
        lp_0 = norm(lp_0)
        
        # Displacement gradient tensor 
        a = l/l_0
        b = (1./sin(phi_0))*((lp/lp_0)*cos(phi)-(l/l_0)*cos(phi_0))
        c = (lp/lp_0)*(sin(phi)/sin(phi_0))
        I_1 = (a**2 + b**2 + c**2) - 2
        I_2 = (a**2)*(c**2) - 1
        
        # Areal strain energy density (Skalak Model)
        ws = (ks/12)*(I_1**2+2*I_1-2*I_2) + (ka/12)*I_2**2
        WS += ws
    print WS
    
if __name__ == '__main__':
    main()