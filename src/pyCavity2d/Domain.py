import numpy as _np
import matplotlib.pyplot as _plt
import ngsolve as _ng
import netgen.occ as _ngocc
from ngsolve.webgui import Draw as _Draw

c_mpers = 299792458

class DomainCreator :
    @classmethod
    def make2d_rectangular(self, cavity_freqnecy_ghz = 0.8, iris_radius = 0.025, iris_thickness = 0.01):
        b = []

        cavity_radius = 2.405*c_mpers/(cavity_freqnecy_ghz*1e9)/2/_np.pi
        cell_length = 1/(cavity_freqnecy_ghz*1e9)/2*c_mpers

        cavity_length = cell_length - 2 *iris_thickness

        b.append([-cavity_length/2-iris_thickness,0])
        b.append([-cavity_length/2-iris_thickness,iris_radius])
        b.append([-cavity_length/2,               iris_radius])
        b.append([-cavity_length/2,               cavity_radius])
        b.append([ cavity_length/2,               cavity_radius])
        b.append([ cavity_length/2,               iris_radius])
        b.append([ cavity_length/2+iris_thickness,iris_radius])
        b.append([ cavity_length/2+iris_thickness,0])
        return b

    @classmethod
    def make_boundary_repeat(self, boundary = "", n = 1):

        if type(boundary) is str :
            b = DomainLoader.load2d_boundary(boundary)
        else :
            b = boundary

        # strip the (-cavLength,0) and (cavLength,0) points
        b = b[1:-1]

        # find caviy length
        cavLength = b[-1][0] - b[0][0]

        # repeat the cavity n times
        b_repeated = []
        for i in range(0,n,1) :
            b_repeated.extend([[bi[0]+cavLength*i,bi[1]] for bi in b])

        # add back the  the (0,0) and (n*cavLength,0) points
        b_repeated.insert(0,[b_repeated[0][0],0])
        b_repeated.append([b_repeated[-1][0],0])
        b_repeated.append([b_repeated[0][0],0])

        return b_repeated

    @classmethod
    def draw(self, b):
        ba = _np.array(b)
        _plt.plot(ba[:,0],ba[:,1])

class DomainWriter :

    def __init__(self):
        pass

    @classmethod
    def write(cls, fileName, boundary):
        f = open(fileName,"w")

        for p in boundary :
            f.write(f"{p[0]} {p[1]}\n")
        f.close()

class DomainLoader :
    def __init__(self):
        pass

    @classmethod
    def load2d_boundary(self, fileName = ""):
        f = open(fileName,"r")

        boundary = []
        for l in f :
            boundary.append([float(v) for v in l.split()])

        return boundary

class Domain2D :

    def __init__(self, boundary, maxh=0.001, single_cell = True) :
    #def __init__(self, boundary = [[-0.02,0], [-0.02,0.005], [-0.015,0.005], [-0.015,0.047], [0.015,0.047],[0.015,0.005], [0.02,0.005], [0.02,0], [-0.02,0]]) :
    #def __init__(self, boundary = [[0,0], [-0.015,0.047], [0.015,0.047], [0.015,-0.047], [-0.015,-0.047], [-0.015,0]]) :
    #def __init__(self, boundary = [[-0.015,0], [-0.015,0.047], [0.015,0.047], [0.015,0],[-0.015,0]]): # SF example
    #def __init__(self, boundary = [[-0.03,0], [-0.03,0.047], [0.03,0.047], [0.03,0],[-0.03,0]]):
        
        self.boundary = boundary
        
        wp = _ngocc.WorkPlane()
        wp.MoveTo(*self.boundary[0])
        for p in self.boundary[1:]:
            wp.LineTo(*p)
        wp.Close().Reverse()
        self.domain = wp.Face()

        if single_cell :
            self.domain.edges.Min(_ngocc.X).name = "zmin"
            self.domain.edges.Max(_ngocc.X).name = "zmax"
            self.domain.edges.Min(_ngocc.X).col = (1, 0, 0) #TODO what is this col?
            self.domain.edges.Max(_ngocc.X).col = (1, 0, 0)

        self.domain.edges.Min(_ngocc.Y).name = "rmin"
        self.domain.edges.Min(_ngocc.Y).col = (1, 0, 0)
        
        geo = _ngocc.OCCGeometry(self.domain, dim=2)

        # mesh
        self.ngmesh = geo.GenerateMesh(maxh=maxh)
        self.mesh = _ng.Mesh(self.ngmesh)

        print(self.mesh.GetBoundaries())

    def draw(self):
        _Draw(self.mesh)

class Domain3D :
    def __init__(self, boundary, maxh=0.005, single_cell = True) :
        coil = _ngocc.Cylinder(_ngocc.Axes((0, 0, 1), _ngocc.Z), r=0.047, h=0.03)

        self.domain = coil

        self.domain.edges.Min(_ngocc.X).name = "zmin"
        self.domain.edges.Max(_ngocc.X).name = "zmax"
        self.domain.edges.Min(_ngocc.X).col = (1, 0, 0)  # TODO what is this col?
        self.domain.edges.Max(_ngocc.X).col = (1, 0, 0)

        geo = _ngocc.OCCGeometry(self.domain, dim=3)

        self.ngmesh = geo.GenerateMesh(maxh=maxh)
        self.mesh = _ng.Mesh(self.ngmesh)
    def draw(self):
        _Draw(self.mesh)