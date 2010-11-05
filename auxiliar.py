

# ----------------------------------------------------------------
# 2. Crear la geometria Lagrangiana para la membrana
# ----------------------------------------------------------------
membrana = sp.triang(shape = 'ico', maxlevel = 0, r = dx ,winding = 0)
referencia = sp.triang(shape = 'ico', maxlevel = 0, r = dx ,winding = 0)
vertices = membrana.vertices.shape
elementos = membrana.faces.shape

# Coordenadas para trasladar el centro de la esfera
i_x = np.ones(membrana.vertices.shape[0])*(x*dx)/2.0
i_y = np.ones(membrana.vertices.shape[0])*(y*dx)/2.0
i_z = np.ones(membrana.vertices.shape[0])*(z*dx)/2.0
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

for t in xrange(1,100):
    print "Paso: ", t
    # ----------------------------------------------------------------
    # 3.1 Propagar velocidad hacia la membrana INTERPOLATION
    # ----------------------------------------------------------------
    U = ibm.interpolation(membrana, x, y, z, dx, u_x, u_y, u_z)
    # ----------------------------------------------------------------
    # 3.2 Encontrar nuevas posiciones de la membrana
    # ----------------------------------------------------------------
    dt = 0.05
    membrana.vertices += U*dt
    # ----------------------------------------------------------------
    # 3.3 Calcular fuerzas en los nodos de la membrana FEM
    # ----------------------------------------------------------------
    F = fem.fem(referencia, membrana)
    # ----------------------------------------------------------------
    # 3.4 Propagar densidad fuerza hacia los nodos del fluido SPREAD
    # ----------------------------------------------------------------
    df = ibm.spread(F,membrana,x,y,z,dx)
    # ----------------------------------------------------------------
    # 3.5 Solucionar dinamica del fluido con nuevas fuerzas LBM
    # ----------------------------------------------------------------
    f = lbm.streaming(f,s,x,y,z)
    f = lbm.collide(f,s,df,x,y,z,omega)
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