# ------------------------------------------------------------------
# Rigid Body Movement              Octubre 04 - 2010 (Terminado)
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------
# Input: pi, pj, pk -> Nodos del elemento NO deformado x,y,z
#        Pi, Pj, Pk -> Nodos del elemento deformado    x,y,z
# ------------------------------------------------------------------
# Output: F -> Vector fuerza por cada nodo en la base  x,y,z
# ------------------------------------------------------------------

# Librerias utilizadas 
import numpy as np
import numpy.linalg as la
import triangular as t

def rotacion(pi,pj,pk,Pi,Pj,Pk):
    # Coordenadas iniciales de los tres nodos respecto al sistema 
    # global de coordenadas i1, i2, i3
    xi = pi[0]
    yi = pi[1]
    zi = pi[2]
    
    xj = pj[0]
    yj = pj[1]
    zj = pj[2]
    
    xk = pk[0]
    yk = pk[1]
    zk = pk[2]
    
    # Coordenadas del elemento deformado de los tres nodos 
    # respecto al sistema global de coordenadas i1, i2, i3
    Xi = Pi[0]
    Yi = Pi[1]
    Zi = Pi[2]
    
    Xj = Pj[0]
    Yj = Pj[1]
    Zj = Pj[2]
    
    Xk = Pk[0]
    Yk = Pk[1]
    Zk = Pk[2]
    
    # Unit vectors Local undeformed local coordinate axis 
    m1 = np.sqrt((xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)
    e1 = np.array([(xj-xi), (yj-yi), (zj-zi)])/m1
    m2 = np.sqrt((xk-xi)**2 + (yk-yi)**2 + (zk-zi)**2)
    e4 = np.array([(xk-xi), (yk-yi), (zk-zi)])/m2
    e3 = np.cross(e1,e4)
    m3 = la.norm(e3)
    e3 = np.cross(e1,e4)/la.norm(e3)
    e2 = np.cross(e3,e1)
    e = np.vstack((e1,e2))
    e = np.vstack((e,e3))
    
    # Matriz de rotacion para la configuracion NO deformada [r]
    d1 = (xj-xi)/m1
    d2 = (yj-yi)/m1
    d3 = (zj-zi)/m1
    e1 = (xk-xi)/m2
    e2 = (yk-yi)/m2
    e3 = (zk-zi)/m2
    f1 = (d2*e3-d3*e2)/m3
    f2 = (d3*e1-d1*e3)/m3
    f3 = (d1*e2-d2*e1)/m3
    g1 = f2*d3 - f3*d2
    g2 = f3*d1 - f1*d3
    g3 = f1*d2 - f2*d1
    r = np.array([[d1, d2, d3],[g1, g2, g3],[f1, f2, f3]])
    
    
    # Unit vectors Local deformed coordinate axis
    M1 = np.sqrt((Xj-Xi)**2 + (Yj-Yi)**2 + (Zj-Zi)**2)
    E1 = np.array([(Xj-Xi), (Yj-Yi), (Zj-Zi)])/M1
    M2 = np.sqrt((Xk-Xi)**2 + (Yk-Yi)**2 + (Zk-Zi)**2)
    E4 = np.array([(Xk-Xi), (Yk-Yi), (Zk-Zi)])/M2
    E3 = np.cross(E1,E4)
    M3 = la.norm(E3)
    E3 = np.cross(E1,E4)/la.norm(E3)
    E2 = np.cross(E3,E1)
    E = np.vstack((E1,E2))
    E = np.vstack((E,E3))
    
    # Matriz de rotacion para la configuracion deformada [R]
    d1 = (Xj-Xi)/M1
    d2 = (Yj-Yi)/M1
    d3 = (Zj-Zi)/M1
    e1 = (Xk-Xi)/M2
    e2 = (Yk-Yi)/M2
    e3 = (Zk-Zi)/M2
    f1 = (d2*e3-d3*e2)/M3
    f2 = (d3*e1-d1*e3)/M3
    f3 = (d1*e2-d2*e1)/M3
    g1 = f2*d3 - f3*d2
    g2 = f3*d1 - f1*d3
    g3 = f1*d2 - f2*d1
    R = np.array([[d1, d2, d3],[g1, g2, g3],[f1, f2, f3]])
    
    # Vectores posicion de los nodos en coordenadas locales estado deformado
    # tiene origen en el nodo i deformado 
    Xil = 0.0
    Yil = 0.0
    Zil = 0.0
    Xjl = Xj - Xi
    Yjl = Yj - Yi
    Zjl = Zj - Zi
    Xkl = Xk - Xi
    Ykl = Yk - Yi
    Zkl = Zk - Zi
    
    # Transformar cada coordenada local del elemento deformado
    # a la base comun i,j,k 
    Pi = np.dot(R, np.array([Xil, Yil, Zil]))
    Pj = np.dot(R, np.array([Xjl, Yjl, Zjl]))
    Pk = np.dot(R, np.array([Xkl, Ykl, Zkl]))    
    
    # Vectores posicion de los nodos en coordenadas locales estado NO deformado
    # tiene origen en el nodo i NO deformado 
    xil = 0.0
    yil = 0.0
    zil = 0.0
    xjl = xj - xi
    yjl = yj - yi 
    zjl = zj - zi
    xkl = xk - xi
    ykl = yk - yi
    zkl = zk - zi
    
    # Transformar cada coordenada local del elemento sin deformar a la base i,j,k 
    pi = np.dot(r, np.array([xil, yil, zil]))
    pj = np.dot(r, np.array([xjl, yjl, zjl]))
    pk = np.dot(r, np.array([xkl, ykl, zkl]))

    
    # Vectores de desplazamientos en el plano
    di = Pi - pi
    dj = Pj - pj
    dk = Pk - pk
    
    # Llama rutina para calcular las fuerzas en el elemento deformado 
    F = t.fuerzas(pi,pj,pk,Pi,Pj,Pk)
    
    # Transformar la fuerza de cada nodo a las coordenadas locales
    Fil = np.dot(np.transpose(R),F[0,:])
    Fjl = np.dot(np.transpose(R),F[1,:])
    Fkl = np.dot(np.transpose(R),F[2,:])
    F = np.vstack((Fil,Fjl))
    F = np.vstack((F,Fkl))
    return F

def main():
    
    # Comprobacipn del algoritmo de traslacion y rotacion de elementos
    # coordendas en DECIMALES
    
    print "---------------------------------------------------------\n"
    pi = np.array([0.0,0.0,0.0])
    pj = np.array([2.0,0.0,0.0])
    pk = np.array([1.0,1.0,0.0])
    
    Pi = np.array([0.0,0.0,0.0])
    Pj = np.array([2.0,0.0,0.0])
    Pk = np.array([1.0,0.0,2.0])
    print "Puntos del elemento inicial:", pi, pj, pk
    print "Puntos del elemento deformado:", Pi, Pj, Pk
    print "Fuerzas resultantes: \n", rotacion(pi,pj,pk,Pi,Pj,Pk)
    print "---------------------------------------------------------\n"
    
if __name__ == '__main__':
    main()