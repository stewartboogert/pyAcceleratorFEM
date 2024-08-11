import ngsolve as _ng
import numpy as _np
from ngsolve import VTKOutput as _VTKOutput
from ngsolve.webgui import Draw as _Draw

q0   = 1.60217663e-19
m0   = 9.1093837e-31
mu0  = 4 * _np.pi * 1e-7
eps0 = 8.85418782e-12
c0   = 299792458

class FieldCalculator :

    def __init__(self, domain):
        self.domain = domain

    def compute(self, nmodes = 10, maxit = 200, axisymmetric = True):
        # define finite element space
        fes = _ng.HCurl(self.domain.mesh, order=1, dirichlet='default')
        u, v = fes.TnT()

        if axisymmetric :
            a = _ng.BilinearForm(_ng.y * _ng.curl(u) * _ng.curl(v) * _ng.dx).Assemble()
            m = _ng.BilinearForm(_ng.y * u * v * _ng.dx).Assemble()
            apre = _ng.BilinearForm(_ng.y * _ng.curl(u) * _ng.curl(v) * _ng.dx + _ng.y * u * v * _ng.dx).Assemble()
            pre = _ng.Preconditioner(apre, "direct", inverse="sparsecholesky")
        else :
            a = _ng.BilinearForm(_ng.curl(u) * _ng.curl(v) * _ng.dx).Assemble()
            m = _ng.BilinearForm(u * v * _ng.dx).Assemble()

            apre = _ng.BilinearForm(_ng.curl(u) * _ng.curl(v) * _ng.dx + u * v * _ng.dx).Assemble()
            pre = _ng.Preconditioner(apre, "direct", inverse="sparsecholesky")

        with _ng.TaskManager():
            a.Assemble()
            m.Assemble()
            apre.Assemble()

            # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
            gradmat, fesh1 = fes.CreateGradient()
            gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
            math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
            math1[0, 0] += 1  # fix the 1-dim kernel
            invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())
            # build the Poisson projector with operator Algebra:
            proj = _ng.IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat
            projpre = proj @ pre.mat

            self.K = a.mat
            self.M = m.mat
            self.precond = pre.mat
            self.eigenvals, self.eigenvecs = _ng.solvers.PINVIT(a.mat, m.mat, pre=projpre, num=nmodes,
                                                                maxit=maxit,
                                                                printrates=False)

        # print out eigenvalues
        self.eigen_freq = []
        for i, lam in enumerate(self.eigenvals):
            self.eigen_freq.append(c0 * _np.sqrt(lam) / (2 * _np.pi) * 1e-6)
            print(i, lam, 'freq: ', c0 * _np.sqrt(lam) / (2 * _np.pi) * 1e-6, "MHz")

        # plot results
        self.gfu_E = []
        self.gfu_H = []
        for i in range(len(self.eigenvecs)):
            w = 2 * _ng.pi * self.eigen_freq[i] * 1e6
            gfu = _ng.GridFunction(fes)
            gfu.vec.data = self.eigenvecs[i]

            self.gfu_E.append(gfu)
            self.gfu_H.append(1j / (mu0 * w) * _ng.curl(gfu)) # H is 1j out of phase with curl of E


    def draw(self, imode = 1, field = "E") :
        if field == "E" :
            _Draw(self.gfu_E[imode], self.domain.mesh,
                  vectors = {"grid_size" : 25, "offset" : 5},
                  # draw_surf=False,
                  order=2)

        elif field == "B" :
            _Draw(_ng.Norm(self.gfu_H[imode]), self.domain.mesh, order=2)
        else :
            print("field = E|B")


    def vtk(self, imode = 1, field = "E"):
        if field == "E" :
            vtk = _VTKOutput(self.domain.mesh,
                             coefs=[self.gfu_E[imode]],
                             names=[f"E_{imode}"],
                             filename=f"E_{imode}", subdivision=2)
        else :
            vtk = _VTKOutput(self.domain.mesh,
                             coefs=[self.gfu_H[imode]],
                             names=[f"B_{imode}"],
                             filename=f"B_{imode}", subdivision=2)
        vtk.Do()

    def plotAxialFields(self, imode=1, field = "E", npt = 100) :

        zmin = self.domain.boundary[0][0]+0.005
        zmax = self.domain.boundary[-2][0]-0.005
        z = _np.linspace(zmin,zmax,npt)
        mip = self.domain.mesh(z,0)
        if field == "E" :
            return [z, self.gfu_E[imode](mip)]
        elif field == "B" :
            return [z, self.gfu_H[imode](mip)]
