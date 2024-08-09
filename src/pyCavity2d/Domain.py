import ngsolve as _ng
import netgen.occ as _ngocc
from ngsolve.webgui import Draw as _Draw

class DomainCreator :
    pass

class DomainLoader :
    def __init__(self):
        pass

    @classmethod
    def load2d_boundary(self, fileName = ""):
        f = open(fileName,"r")

        boundary = []
        for l in f :
            boundary.append(map(float,l.split()))

        return boundary


class Domain2D :

    def __init__(self, boundary = [[-0.15,0], [-0.15,0.03], [-0.1,0.03], [-0.1,0.3], [0.1,0.3], [0.1,0.03], [0.15,0.03], [0.15,0], [-0.15,0]]) :
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

        #self.domain.edges.Min(_ngocc.X).name = "zmin"
        #self.domain.edges.Max(_ngocc.X).name = "zmax"
        #self.domain.edges.Min(_ngocc.X).col = (1, 0, 0) #TODO what is this col?
        #self.domain.edges.Max(_ngocc.X).col = (1, 0, 0)
        self.domain.edges.Min(_ngocc.Y).name = "rmin"
        self.domain.edges.Min(_ngocc.Y).col = (1, 0, 0)
        
        geo = _ngocc.OCCGeometry(self.domain, dim=2)

        # mesh
        self.ngmesh = geo.GenerateMesh(maxh=0.01)
        self.mesh = _ng.Mesh(self.ngmesh)

        print(self.mesh.GetBoundaries())
        
        _Draw(self.mesh)
