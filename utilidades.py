# ------------------------------------------------------------------
# Utilidades                       Septiembre 21 - 2010
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------

# Librerias utilizadas 
import lbm as lbm
import fem as fem 
import ibm as ibm
import sphere as sp
import numpy as np
import visualizacion as vs
from enthought.tvtk.api import tvtk


# ------------------------------------------------------------------
# Function:  leerVTK 
# ------------------------------------------------------------------
# Input: ruta  -> Ruta del archivo a importar 
# ------------------------------------------------------------------
# Output: estructura (Contiene coordenadas de nodos y conexiones)
# ------------------------------------------------------------------
def leerVTK(ruta):
    file = tvtk.UnstructuredGridReader()
    file.file_name = ruta
    file.update()
    s = file.output
    faces = s.get_cells
    vertices = s.points