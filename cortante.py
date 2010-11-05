# ------------------------------------------------------------------
# Boundary Conditions              Septiembre 30 - 2010
# Universidad de Los Andes         CMES, vol.44, no.2, pp.137-155
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------

#-------------------------------------------------------------------
# Librerias utilizadas
# ------------------------------------------------------------------
import numpy as np
import visualizacion as vs
import time as time

#-------------------------------------------------------------------
# Parametros del modelo D3Q19
# ------------------------------------------------------------------

# Velocidad en componentes cartesianas 
e_x = np.array([ 0, 1, -1,  0,  0,  0,  0,  1, -1,  1, -1, 0,  0,  0,   0,  1, -1,  1, -1 ])
e_y = np.array([ 0, 0,  0,  1, -1,  0,  0,  1, -1, -1,  1, 1, -1,  1,  -1,  0,  0,  0,  0 ])
e_z = np.array([ 0, 0,  0,  0,  0,  1, -1,  0,  0,  0,  0, 1, -1, -1,   1,  1, -1, -1,  1 ])
    
# Ponderado para velocidades
w = np.array([(1./3.),
   	(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
   	(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
   	(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.)])
p = w

# Vector de indices opuestos
opuesto = np.array([0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17])
   
#-------------------------------------------------------------------
# Inicializar: Crea el dominio de fluido completo,Rev: 04/Oct/2010
# ------------------------------------------------------------------
# Input: x,y,z
# ------------------------------------------------------------------
# Output: f(x,y,z,l) Estructura con poblaciones f 
# ------------------------------------------------------------------
def inicializar(x,y,z):
    f = np.zeros((2,x,y,z,19))
    f[0:1,0:x,0:y,0:z,:]=w 
    f[1:1,0:x,0:y,0:z,:]=w 
    return f
    
#--------------------------------------------------------------------
# Funcion para calcular propiedades macroscopicas Rev: 04/Oct/2010
#--------------------------------------------------------------------
def fluido(f,x,y,z):
    rho = np.zeros((x,y,z))
    u_x = np.zeros((x,y,z))
    u_y = np.zeros((x,y,z))
    u_z = np.zeros((x,y,z))
    for i in xrange(0,x):
        for j in xrange(0,y):
            for k in xrange(0,z):
                for l in xrange(0,19):
                    rho[i,j,k]+=f[0,i,j,k,l]
                    u_x += f[1][i][j][k][l]*e_x[l]
                    u_y += f[1][i][j][k][l]*e_y[l]
                    u_z += f[1][i][j][k][l]*e_z[l]
            u_x[i,j,k] = u_x[i,j,k]/rho[i,j,k]
            u_y[i,j,k] = u_y[i,j,k]/rho[i,j,k]
            u_z[i,j,k] = u_z[i,j,k]/rho[i,j,k]
    return u_x, u_y, u_z, rho

#--------------------------------------------------------------------
# Streaming - Bounce-back
#--------------------------------------------------------------------
def streaming(f,s,x,y,z):
    for i in xrange(1,x-1):
        for j in xrange(1,y-1):
            for k in xrange(1,z-1):
                for l in xrange(0,19):
                    inv = opuesto[l]
                    if(s[i+e_x[inv],j+e_y[inv],k+e_z[inv]]==1):
                        # Bounce - back
                        f[0,i,j,k,l] = f[1,i,j,k,inv]
                    if(s[i,j,k]==3):
                        # Frontera superior
                        u = 0.0
                        v = 0.0
                        w = 0.0
                        n = nodoSuperior(f[1,i,j,k,:],u,v,w)
                        f[0,i,j,k]=n[1]
                    if(s[i,j,k]==2):
                        # Frontera inferior
                        u = 0.0
                        v = 0.0
                        w = 0.0
                        n = nodoInferior(f[1,i,j,k,:],u,v,w)
                        f[0,i,j,k]=n[1]
                    else:
                        # Streaming - normal
                        f[0,i,j,k,l] = f[1,i+e_x[inv],j+e_y[inv],k+e_z[inv],l]
    return f

#--------------------------------------------------------------------
# Consistent Boundary Conditions
#--------------------------------------------------------------------
def nodoSuperior(f,u,v,w):
    
    # Seleccion de correctores 
    f4a=f[opuesto[4]]
    f8a=f[opuesto[8]]
    f9a=f[opuesto[9]]
    f12a=f[opuesto[12]]
    f14a=f[opuesto[14]]
    
    # Var auxiliares
    rho=(f[0]+f[1]+f[2]+f[5]+f[6]+f[15]+f[16]+f[17]+f[18])/(1+v) + 2*(f[3]+f[7]+f[10]+f[11]+f[13])/(1+v)
    Qx=18*(-f[1]+f[10]-f[15]+f[16]-f[17]+f[18]+f[2]-f[7]+f8a-f9a-rho*u)
    A=v*(f[0]+f[1]+f[10]+f[11]+f12a+f[13]+f14a+f[15]+f[16]+f[17]+f[18]+f[2]+f[3]+f4a+f[5]+f[6]+f[7]+f8a+f9a)
    Qy=6*(-f[10]-f[11]+f12a-f[13]+f14a-f[3]+f4a-f[7]+f8a+f9a+A)/(1+v)
    Qz=18*(-f[11]+f12a+f[13]-f14a-f[15]+f[16]+f[17]-f[18]-f[5]+f[6]+rho* w)
    
    # Funciones desconocidas 
    f[4]=f4a-p[4]*Qy
    f[8]=f8a-p[8]*(Qx+Qy)
    f[9]=f9a+p[9]*(Qx-Qy)
    f[12]=f12a-p[12]*(Qy+Qz)
    f[14]=f14a+p[14]*(Qz-Qy)
    return rho, f
    
def nodoInferior(f,u,v,w):
    
    # Seleccion de correctores 
    f3a=f[3]
    f7a=f[7]
    f10a=f[10]
    f11a=f[11]
    f13a=f[13]
    
    # Var auxiliares
    rho=(f[0]+f[1]+f[15]+f[16]+f[17]+f[18]+f[2]+f[5]+f[6]+2*(f[9]+f[8]+f[14]+f[12]+f[4]))/(1 + v)
    Qx=18*(-f[1]+f10a-f[15]+f[16]-f[17]+f[18]+f[2]-f7a+f[8]-f[9]+rho*u)
    Qy=6*(-f10a-f11a+f[12]-f13a+f[14]-f3a+f[4]-f7a+f[8]+f[9]+rho*v)
    Qz=18*(-f11a+f[12]+f13a-f[14]-f[15]+f[16]+f[17]-f[18]-f[5]+f[6]+rho*w)
    
    # Funciones desconocidas 
    f[3]=f3a+p[3]*Qy
    f[7]=f7a+p[7]*(Qx+Qy)
    f[10]=f10a+p[10]*(Qy-Qx)
    f[11]=f11a+p[11]*(Qy+Qz)
    f[13]=f13a+p[13]*(Qy-Qz)
    return rho, f

#--------------------------------------------------------------------
# Collide - Aproximacion BGK
#--------------------------------------------------------------------
def collide(f,s,x,y,z,omega):
    for i in xrange(1,x-1):
            for j in xrange(1,y-1):
                for k in xrange(1,z-1):
                    rhop = 0.0
                    ux=0.0
                    uy=0.0
                    uz=0.0
                    if(s[i,j,k]!=1):
			        # Celda de fluido normal
		            for l in xrange(0,19):
		                fi = f[0][i][j][k][l]
		                rhop += fi
		                ux += e_x[l] * fi
		                uy += e_y[l] * fi
		                uz += e_z[l] * fi
		            
		            for l in xrange(0,19):
		                tmp = (e_x[l]*ux + e_y[l]*uy + e_z[l]*uz)
		                feq = w[l] * ( rhop - 
		                (3.0/2.0 * (ux*ux + uy*uy + uz*uz)) + 
		                (3.0 * tmp) +
		                ((9.0/2.0) * tmp*tmp ) )
		                f[0][i][j][k][l] = (1.0-omega)*f[0][i][j][k][l]+omega * feq
    return f

def main():
    # Dimensiones de la red 
    H = 100
    x = 32
    y = 32
    z = 32
    dx = H/x
    omega = 1.4
    
    #Magnitud de la velocidad en la frontera inferior y superior
    u = 0.1
    
    # Estructuras de informacion macroscopica
    s = np.zeros((x,y,z))
    
    s[0,:,:]=1
    s[x-1,:,:]=1
    s[:,0,:]=1
    s[:,y-1,:]=1
    s[:,:,0]=1
    s[:,:,z-1]=1

    s[:,:,1]=2
    
    f = inicializar(x,y,z)
    p = fluido(f,x,y,z)
    u_x = p[0]
    u_y = p[1]
    u_z = p[2]
    rho = p[3]
    
    # Ciclo principal del algoritmo LBM
    # ------------------------------------------
    # 1. Calcular Variables macroscopicas
    # 2. Streaming
    #     2.1 Bounce-back
    #     2.2 Fronteras microscopicas
    # 3. Colision 
    # 4. Actualizar 
    # ------------------------------------------
    print "Iniciando..."
    tic = time.clock()
    
    for t in xrange(0,50):
        # 1. Calcular Variables macroscopicas
        print "Paso:", t
        p = fluido(f,x,y,z)
        u_x = p[0]
        u_y = p[1]
        u_z = p[2]
        rho = p[3]
        # 2. Streaming
        f = streaming(f,s,x,y,z)
        f = collide(f,s,x,y,z,omega)
        vs.guardarFluido(x,y,z,dx,rho,u_x,u_y,u_z,t)
        f[1]=f[0]
        toc=time.clock()
    print "Fin:\n", (toc-tic)
if __name__ == '__main__':
    main()
