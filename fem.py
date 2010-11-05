# ------------------------------------------------------------------
# Forces on deformed membrane      Octubre 04 - 2010 (Terminado)
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------
# Input: FV1 -> Geometria NO deformada 
#        FV2 -> Geometria deformada o de referencia
# ------------------------------------------------------------------
# Output: fuerza -> Vector fuerza por cada nodo en la geometria
# ------------------------------------------------------------------

# Librerias utilizadas 
import numpy as np
import rotacion as rt

# Recibe estructura deformada y No deformada (Referencia)
def fem(FV1, FV2):

    # Crear estructura para almacenar las fuerzas en cada nodo
    fuerza = np.zeros((FV1.vertices.shape[0],3))
    
    # Calcular fuerzas sobre cada elemento
    for i in xrange(FV1.faces.shape[0]):
        # Encontrar los nodos del elemento
        ni = FV1.faces[i,0]
        nj = FV1.faces[i,1]
        nk = FV1.faces[i,2]
        
        # Coordenadas de cada nodo en estado inicial
        pi = FV1.vertices[ni]
        pj = FV1.vertices[nj]
        pk = FV1.vertices[nk]
    
        # Coordenadas de cada nodo en estado deformado 
        Pi = FV2.vertices[ni]
        Pj = FV2.vertices[nj]
        Pk = FV2.vertices[nk]
        
        # Calculo de la fuerza en coordenadas globales para el elemento
        F = rt.rotacion(pi,pj,pk,Pi,Pj,Pk)
        
        # Vector fuerza por componentes asociado a cada nodo
        fuerza[ni] += F[0]
        fuerza[nj] += F[1]
        fuerza[nk] += F[2]
        
    return fuerza