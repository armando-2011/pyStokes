# ------------------------------------------------------------------
# Boundary Conditions              Rev: 16/Oct/2010 (Terminado)
# Universidad de Los Andes         arXiv:0811.4593v5
# Autor: Oscar Castillo O.         ol.castillo28@uniandes.edu.co
# ------------------------------------------------------------------

#-------------------------------------------------------------------
# Librerias utilizadas
# ------------------------------------------------------------------
import numpy as np
import visualizacion as vs
import time as time

#-------------------------------------------------------------------
# Parametros del modelo D3Q19 Rev: 16/Oct/2010
# ------------------------------------------------------------------

# Velocidad en componentes cartesianas 
e_x = np.array([1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1,0,0,0,0,0], dtype=float)
e_y = np.array([0,0,1,-1,0,0,1,-1,0,0,1,-1,0,0,1,1,-1,-1,0], dtype=float)
e_z = np.array([0,0,0,0,1,-1,0,0,1,-1,0,0,1,-1,1,-1,1,-1,0], dtype=float)
    
# Ponderado para velocidades
w = np.array([(2./36.),(2./36.),(2./36.),(2./36.),(2./36.),(2./36.)
,(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
(1./36.),(1./36.),(1./36.),(12./36.)])

# Vector de indices opuestos
opuesto = np.array([1,0,3,2,5,4,11,10,13,12,7,6,9,8,17,16,15,14,18])
current = 0
other = 1

#-------------------------------------------------------------------
# Inicializar: Crea el dominio de fluido completo,Rev: 16/Oct/2010
# ------------------------------------------------------------------
# Input: x,y,z
# ------------------------------------------------------------------
# Output: f(x,y,z,l) Estructura con poblaciones f 
# ------------------------------------------------------------------
def inicializar(x,y,z):
    f = np.zeros((2,x,y,z,19))
    f[0:2,:,:,:,:]=w
    return f
    
#--------------------------------------------------------------------
# Funcion para calcular propiedades macroscopicas Rev: 16/Oct/2010
#--------------------------------------------------------------------
def fluido(f,s,x,y,z):
    rho = np.zeros((x,y,z))
    u_x = np.zeros((x,y,z))
    u_y = np.zeros((x,y,z))
    u_z = np.zeros((x,y,z))
    
    for i in xrange(0,x):
        for j in xrange(0,y):
            for k in xrange(0,z):
                rho[i,j,k]=np.sum(f[0,i,j,k,:])
                u_x[i,j,k]=np.dot(f[0,i,j,k,:],e_x)/rho[i,j,k]
                u_y[i,j,k]=np.dot(f[0,i,j,k,:],e_y)/rho[i,j,k]
                u_z[i,j,k]=np.dot(f[0,i,j,k,:],e_z)/rho[i,j,k]
    
    for i in xrange(0,x):
        for j in xrange(0,y):
            k=0
            u_x[i,j,k]=-0.1
            rho[i,j,k]=nodoInferior(f[0,i,j,k,],-0.1,0.0,0.0)[0]
                    
    for i in xrange(0,x):
        for j in xrange(0,y):
            k=z-1
            u_x[i,j,k]=0.1
            rho[i,j,k]=nodoSuperior(f[0,i,j,k,],0.1,0.0,0.0)[0]
    return u_x, u_y, u_z, rho

#--------------------------------------------------------------------
# Streaming - Bounce-back
#--------------------------------------------------------------------
def streaming(f,s,x,y,z):
                       
    for i in xrange(0,x):
        for j in xrange(0,y):
            for k in xrange(1,z-1):
                for l in xrange(0,19):
                    inv = opuesto[l]
                    a=i+e_x[inv]
                    b=j+e_y[inv]
                    c=k+e_z[inv]
                    # Periodico en x
                    if(a<0):a=x-1
                    if(a>x-1):a=0
                    
                    # Periodico en y
                    if(b<0):b=y-1
                    if(b>y-1):b=0
                    
                    if(s[a,b,c]==1):
                        # Bounce - back
                        f[current,i,j,k,l] = f[other,i,j,k,inv]                        
                    else:
                        # Streaming - normal
                        f[current,i,j,k,l] = f[other,a,b,c,l]
    for i in xrange(0,x):
        for j in xrange(0,y):
            k=0
            # Nodo inferior
            f[current,i,j,k,:] = nodoInferior(f[other,i,j,k,],-0.1,0.0,0.0)[1]
                    
    for i in xrange(0,x):
        for j in xrange(0,y):
            k=z-1
            # Nodo superior
            f[current,i,j,k,:] = nodoSuperior(f[other,i,j,k,],0.1,0.0,0.0)[1]
    return f
    
#--------------------------------------------------------------------
# Consistent Boundary Conditions
#--------------------------------------------------------------------
def nodoSuperior(f,u,v,w):
    
    A=f[0]+f[1]+f[2]+f[3]+f[6]+f[10]+f[11]+f[7]+f[18]
    B=f[4]+f[8]+f[12]+f[14]+f[16]
    rho = (A+2*B)/(w+1)
    
    Nx=(1./2.)*(f[0]+f[6]+f[7]-(f[1]+f[10]+f[11]))-(1./3.)*rho*u
    Ny=(1./2.)*(f[2]+f[6]+f[10]-(f[3]+f[7]+f[11]))-(1./3.)*rho*v
    
    f[5]=f[4]-(1./3.)*rho*w
    f[9]=f[12]+(rho/6)*(-w+u)-Nx
    f[13]=f[8]+(rho/6)*(-w-u)+Nx
    f[15]=f[16]+(rho/6)*(-w+v)-Ny
    f[17]=f[14]+(rho/6)*(-w-v)+Ny
    return rho, f
    
def nodoInferior(f,u,v,w):
    
    A=f[0]+f[1]+f[2]+f[3]+f[6]+f[7]+f[10]+f[11]+f[18]
    B=f[5]+f[9]+f[13]+f[15]+f[17]
    rho = (A+2*B)/(1-w)
    
    Nx=(1./2.)*(f[0]+f[6]+f[7]-(f[1]+f[10]+f[11]))-(1./3.)*rho*u
    Ny=(1./2.)*(f[2]+f[6]+f[10]-(f[3]+f[7]+f[11]))-(1./3.)*rho*v
    
    f[4]=f[5]+(1./3.)*rho*w
    f[8]=f[13]+(rho/6)*(w+u)-Nx
    f[12]=f[9]+(rho/6)*(w-u)+Nx
    f[14]=f[17]+(rho/6)*(w+v)-Ny
    f[16]=f[15]+(rho/6)*(w-v)+Ny
    return rho, f

#--------------------------------------------------------------------
# Collide - Aproximacion BGK
#--------------------------------------------------------------------
def collide(f,s,x,y,z,omega):
    
    rho = np.zeros((x,y,z))
    u_x = np.zeros((x,y,z))
    u_y = np.zeros((x,y,z))
    u_z = np.zeros((x,y,z))
    
    for i in xrange(0,x):
            for j in xrange(0,y):
                for k in xrange(0,z):
                    rho[i,j,k]=np.sum(f[0,i,j,k,:])
                    u_x[i,j,k]=np.dot(f[0,i,j,k,:],e_x)/rho[i,j,k]
                    u_y[i,j,k]=np.dot(f[0,i,j,k,:],e_y)/rho[i,j,k]
                    u_z[i,j,k]=np.dot(f[0,i,j,k,:],e_z)/rho[i,j,k]
                    if(s[i,j,k]!=1):
                        rhop=rho[i,j,k]
                        ux=u_x[i,j,k]
                        uy=u_y[i,j,k]
                        uz=u_z[i,j,k]
                        for l in xrange(0,19):
    		                tmp = (e_x[l]*ux + e_y[l]*uy + e_z[l]*uz)
    		                feq = w[l] * ( rhop - 
    		                (3.0/2.0 * (ux*ux + uy*uy + uz*uz)) + 
    		                (3.0 * tmp) +
    		                ((9.0/2.0) * tmp*tmp ) )
    		                f[current,i,j,k,l] = (1.0-omega)*f[1,i,j,k,l]+omega * feq
    return f

def main():
    # Dimensiones de la red 
    H = 100.0
    x = 10
    y = 10
    z = 10
    dx = H/x
    omega = 1.0

    # Estructuras de informacion macroscopica
    f = inicializar(x,y,z)
    s = np.zeros((x,y,z))
    s[:,:,0]=3
    s[:,:,z-1]=2
   
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
    
    for t in xrange(500):
        # 1. Calcular Variables macroscopicas
        print "Paso:", t
        p = fluido(f,s,x,y,z)
        u_x = p[0]
        u_y = p[1]
        u_z = p[2]
        rho = p[3]
        # 2. Streaming
        f = streaming(f,s,x,y,z)
        f = collide(f,s,x,y,z,omega)
        vs.guardarFluido(x,y,z,dx,rho,u_x,u_y,u_z,t)
        toc=time.clock()
        f[other]=f[current]
    print "Fin:\n", (toc-tic)
if __name__ == '__main__':
    main()