# ------------------------------------------------------------------
# Tank-treading                    Octubre 18 - 2010
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
import time as time

# ----------------------------------------------------------------
# 0. Conversion a cantidades adimensionales y parametros de flujo
# ----------------------------------------------------------------
Re = 0.01 
ks = 0.2
u = 0.0002
tau = 1.0

# ----------------------------------------------------------------
# 1. Crear la malla Euleriana para el fluido
# ----------------------------------------------------------------
# Dimensiones
H = 100.0
x = 10
y = 10
z = 10
dx = (H/x)*0.0002
tau = 1.0
omega = 1.0/tau
R = 4*dx
current =0
other=1

# Estructuras de informacion macroscopica
f = lbm.inicializar(x,y,z)
s = np.zeros((x,y,z))

# Condiciones de frontera
s[:,:,0]=3
s[:,:,z-1]=2
df = np.zeros((x,y,z,3))
p = lbm.fluido(f,s,x,y,z)
u_x = p[0]
u_y = p[1]
u_z = p[2]
rho = p[3]
df = np.zeros((x,y,z,3))

# Coordenadas cartesianas de cada punto en el fluido
mesh = np.zeros(3)
for i in xrange(0,x):
    for j in xrange(0,y):
        for k in xrange(0,z):
            mesh = np.vstack((mesh,(i,j,k)))
fluido = mesh[1:x*y*z+1,:]*dx
fluido = fluido.reshape(x,y,z,3)

# ----------------------------------------------------------------
# 2. Crear la geometria Lagrangiana para la membrana
# ----------------------------------------------------------------
membrana = sp.triang(shape = 'ico', maxlevel = 1, r = R ,winding = 0)
referencia = sp.triang(shape = 'ico', maxlevel = 1, r = R ,winding = 0)
vertices = membrana.vertices.shape
elementos = membrana.faces.shape

# Coordenadas para trasladar el centro de la esfera
i_x = np.ones(membrana.vertices.shape[0])*((x-1)*dx)/2.0
i_y = np.ones(membrana.vertices.shape[0])*((y-1)*dx)/2.0
i_z = np.ones(membrana.vertices.shape[0])*((z-1)*dx)/2.0
incremento = np.vstack((i_x,i_y))
incremento = np.vstack((incremento,i_z))
membrana.vertices +=np.transpose(incremento)
referencia.vertices +=np.transpose(incremento)
U = np.zeros(vertices)
F = np.zeros(vertices)
D = sp.calcularArea(membrana) - sp.calcularArea(referencia)

# ----------------------------------------------------------------
# 3. Ciclo de simulacion principal 
# ----------------------------------------------------------------

vs.guardarMembrana(membrana,U,F,D,0)
vs.guardarFluido(x,y,z,dx,rho,u_x,u_y,u_z,0)

print "Iniciando...."
tic = time.clock()

for t in xrange(1,200):
    print "Paso: ", t
    # ----------------------------------------------------------------
    # 3.1 Propagar velocidad hacia la membrana INTERPOLATION
    # ----------------------------------------------------------------
    U = ibm.interpolation(membrana, fluido, x, y, z, dx, u_x, u_y, u_z)
    # ----------------------------------------------------------------
    # 3.2 Encontrar nuevas posiciones de la membrana
    # ----------------------------------------------------------------
    dt = 0.1
    membrana.vertices += U*dt
    # ----------------------------------------------------------------
    # 3.3 Calcular fuerzas en los nodos de la membrana FEM
    # ----------------------------------------------------------------
    F = fem.fem(referencia, membrana)
    # ----------------------------------------------------------------
    # 3.4 Propagar densidad fuerza hacia los nodos del fluido SPREAD
    # ----------------------------------------------------------------
    df = ibm.spread(F,membrana,fluido,x,y,z,dx)
    # ----------------------------------------------------------------
    # 3.5 Solucionar dinamica del fluido con nuevas fuerzas LBM
    # ----------------------------------------------------------------
    f = lbm.streaming(f,s,x,y,z)
    f = lbm.collide(f,s,(-1.0)*df,x,y,z,omega)
    # ----------------------------------------------------------------
    # 3.6 Calcular propiedades macroscopicas para el fluido 
    # ----------------------------------------------------------------
    p = lbm.fluido(f,s,x,y,z)
    u_x = p[0]
    u_y = p[1]
    u_z = p[2]
    rho = p[3]
    f[other]=f[current]
    # ----------------------------------------------------------------
    # 3.7 Calcular propiedades macroscopicas para la membrana
    # ----------------------------------------------------------------
    D = sp.calcularArea(membrana) - sp.calcularArea(referencia)
    # ----------------------------------------------------------------
    # 3.8 Visualizacion
    # ----------------------------------------------------------------    
    vs.guardarMembrana(membrana,U,F,D,t)
    vs.guardarFluido(x,y,z,dx,rho,u_x,u_y,u_z,t)

toc = time.clock()
print "Tiempo de %d iteraciones [s]:"%t, (toc-tic)
print "Fin - simulacion"
