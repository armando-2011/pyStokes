import rotacion as r
import triangular as t
import numpy as np

def main():
    
    # Comprobacipn del algoritmo de traslacion y rotacion de elementos
    print "---------------------------------------------------------\n"
    pi = np.array([0.0,0.0,0.0])
    pj = np.array([2.0,0.0,0.0])
    pk = np.array([1.0,1.0,0.0])
    Pi = np.array([0.0,0.0,0.0])
    Pj = np.array([2.0,0.0,0.0])
    Pk = np.array([1.0,2.0,0.0])
    print "Puntos del elemento inicial:", pi, pj, pk
    print "Puntos del elemento deformado:", Pi, Pj, Pk
    print "Fuerzas resultantes t: \n", t.fuerzas(pi,pj,pk,Pi,Pj,Pk)
    print "Fuerzas resultantes r: \n", r.rotacion(pi,pj,pk,Pi,Pj,Pk)
    print "---------------------------------------------------------\n"
    
if __name__ == '__main__':
    main()