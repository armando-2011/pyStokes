import numpy as np
from enthought.tvtk.api import tvtk

class MeshTri:
    """Malla triangular de superficie esferica"""
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces
        

def triang(shape, maxlevel, r, winding):
    """
    triang - generate a triangle mesh approximating a sphere
    
    Usage: FV = triang(shape,Nrecurse,r,winding)
    
      shape is a string, either of the following:
      'ico'   starts with icosahedron (most even, default)
      'oct'   starts with octahedron
      'tetra' starts with tetrahedron (least even)
    
      Nrecurse is int >= 0, setting the recursions (default 0)
    
      r is the radius of the sphere (default 1)
    
      winding is 0 for clockwise, 1 for counterclockwise (default 0).  The
      matlab patch command gives outward surface normals for clockwise
      order of vertices in the faces (viewed from outside the surface).
    
      FV has fields FV.vertices and FV.faces.  The vertices 
      are listed in clockwise order in FV.faces, as viewed 
      from the outside in a RHS coordinate system.
    
    The function uses recursive subdivision.  The first
    approximation is an platonic solid, either an  icosahedron,
    octahedron or a tetrahedron.  Each level of refinement 
    subdivides each triangle face by a factor of 4 (see also 
    mesh_refine).  At each refinement, the vertices are 
    projected to the sphere surface (see sphere_project).
    
    A recursion level of 3 or 4 is a good sphere surface. """
    
    def tetrafun():
        sqrt_3 = 1./np.sqrt(3.)
        
        # +X, +Y, +Z  - PPP
        # -X, -Y, +Z  - MMP
        # -X, +Y, -Z  - MPM
        # +X, -Y, -Z  - PMM
        
        v = np.array([  [sqrt_3,  sqrt_3,  sqrt_3], \
                    [-sqrt_3, -sqrt_3,  sqrt_3],\
                    [-sqrt_3,  sqrt_3, -sqrt_3],\
                    [sqrt_3, -sqrt_3, -sqrt_3] ])

        # Structure describing a tetrahedron
        
        f = np.array([ [0, 1, 2],\
                    [0, 3, 1],\
                    [2, 1, 3],\
                    [3, 0, 2] ])

        FV = MeshTri(v,f)
        
        return FV
        
    def octfun():
        # Six equidistant points lying on the unit sphere
        v = np.array([  [1,  0,  0],\
                        [-1,  0,  0],\
                        [0,  1,  0],\
                        [0, -1,  0],\
                        [0,  0,  1],\
                        [0,  0, -1] ])

        # Join vertices to create a unit octahedron
        f = np.array([  [0, 4, 2],\
                        [2, 4, 1],\
                        [1, 4, 3],\
                        [3, 4, 0],\
                        [0, 2, 5],\
                        [2, 1, 5],\
                        [1, 3, 5],\
                        [3, 0, 5] ])
                        
        FV = MeshTri(v,f)
        return FV
        
    def icofun():
        t = (1.+np.sqrt(5.))/2.
        tau = t/np.sqrt(1.+t**2)
        one = 1./np.sqrt(1.+t**2) # Unit sphere
        
        # Twelve vertices of icosahedron on unit sphere
        
        v = np.array([[  tau,  one,    0 ],
                     [ -tau,  one,    0 ],
                     [ -tau, -one,    0 ],
                     [  tau, -one,    0 ],
                     [  one,   0 ,  tau ],
                     [  one,   0 , -tau ],
                     [ -one,   0 , -tau ],
                     [ -one,   0 ,  tau ],
                     [   0 ,  tau,  one ],
                     [   0 , -tau,  one ],
                     [   0 , -tau, -one ],
                     [   0 ,  tau, -one ]])

        # Structure for unit icosahedron
        f = np.array([  [4,  7,  8],
                        [4, 9,  7],
                        [5, 11,  6],
                        [5,  6, 10],
                        [0,  3,  4],
                        [0,  5,  3],
                        [2,  1,  7],
                        [2,  6,  1],
                        [8, 11,  0],
                        [8,  1, 11],
                        [9,  3, 10],
                        [9, 10,  2],
                        [8,  0,  4],
                        [11,  5,  0],
                        [4,  3, 9],
                        [5, 10,  3],
                        [7,  1,  8],
                        [6, 11,  1],
                        [7, 9,  2],
                        [6,  2, 10], ])
                   
        FV = MeshTri(v,f)
        return FV
        
    def errhandler():
       print "Your input has not been recognised"
        
    case = {
        "tetra" : tetrafun,
        "oct"   : octfun,
        "ico"   : icofun
    }
    
    FV = case.get(shape,errhandler)()
    
    # -----------------
    # refine the starting shapes with subdivisions
    if maxlevel:
        # Subdivide each starting triangle (maxlevel) times
        for level in xrange(maxlevel):
            # Subdivide each triangle and normalize the new points thus
            # generated to lie on the surface of a sphere radius r.
            mesh_refine_tri4(FV);
            sphere_project(FV.vertices,r);
        
    if winding:
        print('...returning counterclockwise vertex order (viewed from outside)\n')
        FV.faces = FV.faces[:,[0,2,1]]
    else:
        print('...returning clockwise vertex order (viewed from outside)\n');
    
    return FV
    
def mesh_refine_tri4(FV):
    """mesh_refine_tri4 - creates 4 triangle from each triangle of a mesh
    
    [ FV ] = mesh_refine_tri4( FV )
    
    FV.vertices   - mesh vertices (Nx3 matrix)
    FV.faces      - faces with indices into 3 rows
                    of FV.vertices (Mx3 matrix)
    
    For each face, 3 new vertices are created at the 
    triangle edge midpoints.  Each face is divided into 4
    faces and returned in FV.
    
           B
          /\
         /  \
       a/____\b       Construct new triangles
       /\    /\       [A,a,c]
      /  \  /  \      [a,B,b]
     /____\/____\     [c,b,C]
    A      c     C    [a,b,c]
    
    It is assumed that the vertices are listed in clockwise order in
    FV.faces (A,B,C above), as viewed from the outside in a RHS coordinate
    system."""
    
    print('...refining mesh (tri4)...')
    
    # NOTE
    # The centroid is located one third of the way from each vertex to 
    # the midpoint of the opposite side. Each median divides the triangle 
    # into two equal areas; all the medians together divide it into six 
    # equal parts, and the lines from the median point to the vertices 
    # divide the whole into three equivalent triangles.
    #
    # Each input triangle with vertices labelled [A,B,C] as shown
    # below will be turned into four new triangles:
    #
    # Make new midpoints
    # a = (A+B)/2
    # b = (B+C)/2
    # c = (C+A)/2
    #
    #        B
    #       /\
    #      /  \
    #    a/____\b       Construct new triangles
    #    /\    /\       [A,a,c]
    #   /  \  /  \      [a,B,b]
    #  /____\/____\     [c,b,C]
    # A        c     C    [a,b,c]
    #

    # Initialise a new vertices and faces matrix
    Nvert = FV.vertices.shape[0]
    Nface = FV.faces.shape[0]
    V2 = np.zeros((Nface*3,3),dtype=float);
    F2 = np.zeros((Nface*4,3),dtype=int);

    for f in xrange(Nface):

        # Get the triangle vertex indices
        NA, NB, NC = FV.faces[f]
        
        # Get the triangle vertex coordinates
        A = FV.vertices[NA]
        B = FV.vertices[NB]
        C = FV.vertices[NC]

        # Now find the midpoints between vertices
        a = (A + B)/2
        b = (B + C)/2
        c = (A + C)/2

        # Find the length of each median
        # A2blen = sqrt ( sum( (A - b).^2, 2 ) )
        # B2clen = sqrt ( sum( (B - c).^2, 2 ) )
        # C2alen = sqrt ( sum( (C - a).^2, 2 ) )

        # Store the midpoint vertices, while
        # checking if midpoint vertex already exists
        Na = mesh_find_vertex(FV,a)
        Nb = mesh_find_vertex(FV,b)
        Nc = mesh_find_vertex(FV,c)

        # Create new faces with orig vertices plus midpoints
        F2[(f+1)*4-4,:] = [ NA, Na, Nc ]
        F2[(f+1)*4-3,:] = [ Na, NB, Nb ]
        F2[(f+1)*4-2,:] = [ Nc, Nb, NC ]
        F2[(f+1)*4-1,:] = [ Na, Nb, Nc ]
        
    # Replace the faces matrix
    FV.faces = F2

def mesh_find_vertex(FV,vertex):
    from matplotlib.mlab import find
    
    V0 = find(FV.vertices[:,0] == vertex[0])
    V1 = find(FV.vertices[:,1] == vertex[1])
    V2 = find(FV.vertices[:,2] == vertex[2])
    V01 = np.intersect1d(V0,V1)
    V012 = np.intersect1d(V01,V2)
    if V012.shape[0] > 0:
        if V012.shape[0] == 1:
            N = V012[0]
        else:
            print 'replicated vertices'
            exit(1)
    else:
        N,M=FV.vertices.shape
        temp = np.append(FV.vertices,vertex)
        FV.vertices = temp.reshape(N+1,M)
        N = FV.vertices.shape[0]-1
    return N

def sphere_project(v,r):
    """
    sphere_project - project point X,Y,Z to the surface of sphere radius r
    
    V = sphere_project(v,r,c)
    
    Cartesian inputs:
    v is the vertex matrix, Nx3 (XYZ)
    r is the sphere radius, 1x1 (default 1)
    c is the sphere centroid, 1x3 (default 0,0,0)
    
    XYZ are converted to spherical coordinates and their radius is
    adjusted according to r, from c toward XYZ (defined with theta,phi)
     
    V is returned as Cartesian 3D coordinates
    """
    from numpy.linalg import norm
    for i in xrange(v.shape[0]):
        magn = norm(v[i])
        vecnorm = v[i]/norm(v[i])
        v[i] = vecnorm*r
        
def save_xml(ug, file_name):
    """Save the unstructured grid dataset to a VTK XML file."""
    w = tvtk.XMLUnstructuredGridWriter(input=ug, file_name=file_name)
    w.write()


def computeNormals(FV):
    """ """
    t = FV.faces.shape[0]
    normals = np.zeros((t,3))
    for i in range(0,t):
        p0 = FV.vertices[FV.faces[i][0]]
        p1 = FV.vertices[FV.faces[i][1]]
        p2 = FV.vertices[FV.faces[i][2]]
        v1 = p1 - p0
        v2 = p2 - p0
        n = np.cross(v1,v2)
        normals[i] = n/np.sqrt(np.dot(n,n))
    return normals
    
def metrica(FV):
       
    # Calcular distancia media entre nodos
    a=0
    for i in xrange(0,FV.faces.shape[0]):
        e = FV.faces[i,:]
        p1 = FV.vertices[e[0]]
        p2 = FV.vertices[e[1]]
        p3 = FV.vertices[e[2]]
        l1 = np.linalg.norm(p3-p1)
        l2 = np.linalg.norm(p3-p2)
        l3 = np.linalg.norm(p2-p1)
        a += (l1 + l2 + l3)/3
    lm = a/FV.faces.shape[0]
    dx=2*np.ceil(lm)
    print dx, lm, FV.faces.shape[0]


def calcularArea(FV):
    """Calcula el area de una esfera"""
    t = FV.faces.shape[0]
    area = np.zeros(t)
    for i in range(0,t):
        p0 = FV.vertices[FV.faces[i][0]]
        p1 = FV.vertices[FV.faces[i][1]]
        p2 = FV.vertices[FV.faces[i][2]]
        v1 = p1 - p0
        v2 = p2 - p0
        a = np.cross(v1,v2)
        a = np.sqrt(np.dot(a,a))/2
        area[i]=a
    return area

def main():
    
    FV = triang(shape = 'ico', maxlevel = 2, r = 4.0e-6 ,winding = 0)
    i_x = np.ones(FV.vertices.shape[0])*95/2
    i_y = np.ones(FV.vertices.shape[0])*95/2
    i_z = np.ones(FV.vertices.shape[0])*95/2
    incremento = np.vstack((i_x,i_y))
    incremento = np.vstack((incremento,i_z))
    print FV.vertices.shape
    print FV.faces.shape
    FV.vertices = FV.vertices + np.transpose(incremento)
    tet_type = tvtk.Triangle().cell_type
    ug = tvtk.UnstructuredGrid(points=FV.vertices)
    ug.set_cells(tet_type, FV.faces)
    velocity = np.random.random_integers(0,100,FV.faces.shape)*0.01
    
    # Asignacion de atributos a nodos
    ug.point_data.vectors = velocity
    ug.point_data.vectors.name = 'velocity'
    curvature = np.ones(FV.vertices.shape[0],dtype=float)
    ug.point_data.scalars = curvature
    ug.point_data.scalars.name = 'curvature'
    
    # Asignacion de atributos a celdas
    
    # Vector normal a cada cara
    normals = computeNormals(FV)
    ug.cell_data.normals = normals
    ug.cell_data.normals.name = 'Normal-cara'
    
    # Escalar que indica tension en cada cara
    tension = np.random.random_integers(0,100,FV.faces.shape[0])*0.01
    ug.cell_data.scalars = tension
    ug.cell_data.scalars.name = 'Tension-cara'
    
    # Atributo vectorial por cada cara
    vector = np.random.random_integers(0,100,FV.faces.shape)*0.01
    ug.cell_data.vectors = vector
    ug.cell_data.vectors.name = 'Vectores-cara'
    
    # Calcular propiedades
    metrica(FV)
    
    save_xml(ug, 'esfera.vtu')
    
if __name__ == '__main__':
    main()