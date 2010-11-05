# ------------------------------------------------------------------
# Visualizacion                    Septiembre 16 - 2010
# Universidad de Los Andes         arXiv:1004.2416v1
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------

# Librerias utilizadas 
from enthought.tvtk.api import tvtk
import numpy as np

# ------------------------------------------------------------------
# Function: Guardar membrana
# ------------------------------------------------------------------
# Input: FV  -> Geometria NO deformada o de referencia
#        U   -> Velocidad de la malla Lagrangiana
#        F   -> Fuerza en cada nodo de la malla Lagrangiana
#        D   -> Estructura de cambio de area por elemento
#        n   -> Tiempo de la iteracion
# ------------------------------------------------------------------
# Output: membrana.vtu
# ------------------------------------------------------------------

def guardarMembrana(FV, U, F, D, n):
    tet_type = tvtk.Triangle().cell_type
    ug = tvtk.UnstructuredGrid(points=FV.vertices)
    ug.set_cells(tet_type, FV.faces)
    
    # Asignacion de atributos a nodos
    ug.point_data.vectors = U
    ug.point_data.vectors.name = 'Velocidad'
    
    # Asignacion de atributos a caras
    ug.cell_data.scalars = D
    ug.cell_data.scalars.name = 'Cambio de area'
    
    w = tvtk.XMLUnstructuredGridWriter(input=ug, file_name='temp/membrana%d.vtu'%n)
    w.write()

# ------------------------------------------------------------------
# Function: Guardar fluido
# ------------------------------------------------------------------
# Input: FV  -> Geometria NO deformada o de referencia
#        U   -> Velocidad de la malla Lagrangiana
#        F   -> Fuerza en cada nodo de la malla Lagrangiana
#        D   -> Estructura de cambio de area por elemento
#        n   -> Tiempo de la iteracion
# ------------------------------------------------------------------
# Output: fluido.vtk
# ------------------------------------------------------------------

def guardarFluido(x,y,z,dx,rho,u_x,u_y,u_z,n):
    # Crear las coordenadas de cada punto 
    mesh = np.zeros(3)
    for i in xrange(0,x):
        for j in xrange(0,y):
            for k in xrange(0,z):
                mesh = np.vstack((mesh,(i,j,k)))
    sg = tvtk.StructuredGrid(dimensions=(x,y,z))
    sg.points = mesh[1:x*y*z+1,:]*dx
    
    # Densidad
    sg.point_data.scalars = np.ravel(rho)
    sg.point_data.scalars.name = 'Densidad'
    
    # Vectores velocidad
    u_x = np.ravel(u_x)
    u_y = np.ravel(u_y)
    u_z = np.ravel(u_z)
    velocidad = np.vstack((np.vstack((u_x,u_y)),u_z))
    sg.point_data.vectors = np.transpose(velocidad)
    sg.point_data.vectors.name = 'Velocidad'
    
    w = tvtk.StructuredGridWriter(input=sg, file_name='temp/densidad%d.vtk'%n)
    w.write()