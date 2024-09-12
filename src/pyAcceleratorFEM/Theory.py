import numpy as _np

def CavityFrequencyFromRadius(radius) :
    pass

def CavityRadriusFromFrequency(freq) :
    pass

def AxialElectricField(t = 0, n = [-1,0], b = [1,1], l = 1, m = 1 , cell_length = 0.2, cavity_length = 1.0, npts = 200, E0 = 1.0) :
    #
    # t time
    # bp [0... \infty]
    # bn [-1, -2, ... \infty]
    # psi = l/m pi (phase advance
    # cell_length
    # cavity_length

    psi = l/m * _np.pi

    z  = _np.linspace(-cavity_length/2,cavity_length/2, npts)
    Ez = 1j*_np.zeros(npts, dtype=_np.complex128)

    for ni, bi in zip(n,b) :
        kn = (psi + 2*_np.pi*ni)/cell_length
        Ez += bi * _np.exp(-1j * kn*z)

    Ez = Ez/len(b)
    return z, Ez

def ComputeErAndBpFromEz(z, Ez, r = 0.01 ,freq = 802e6) :
    Er = _np.gradient(Ez, z)/2 * r
    Bp = Ez / 2 / 3e8 ** 2 * (2 * _np.pi * freq) * r

    return Er, Bp

def CavityBodyTransverseMatrix(gammaI, gammaF, L, alpha, eta, deltaPhi) :
    '''
    Calculates RF cavity body transverse matrix.

    Example:

    >>> deltaPhi = 0
    >>> gammaI = 1
    >>> gammaF = 3
    >>> eta = CavityBodyEta([1],[0],deltaPhi)
    >>> alpha = CavityBodyAlpha(gammaI,gammaF,eta,deltaPhi)
    >>> m = CavityFringeTransverseMatrix(gammaI, gammaF, alpha, eta, deltaPhi)
        array([[ 0.92550932,  0.53559779],
               [-0.0892663 ,  0.30850311]])

    +-----------------+---------------------------------------------------------+
    | **Parameters**  | **Description**                                         |
    +-----------------+---------------------------------------------------------+
    | gammaI          | Incoming energy (gamma units do not matter)             |
    +-----------------+---------------------------------------------------------+
    | gammaF          | Outgoing energy (or gamma, units do not matter)         |
    +-----------------+---------------------------------------------------------+
    | alpha           | Entrance (True) or exit (False)                         |
    +-----------------+---------------------------------------------------------+
    | eta             | Energy gain (or gamma, units do not matter). If None    |
    +-----------------+---------------------------------------------------------+
    | deltaPhi        | Energy gain (or gamma, units do not matter). If None    |
    +-----------------+---------------------------------------------------------+
    '''

    gammaPrime = (gammaF - gammaI)/L
    return _np.array([[_np.cos(alpha), _np.sqrt(8/eta)*gammaI/gammaPrime*_np.cos(deltaPhi)*_np.sin(alpha)],
                      [-_np.sqrt(eta/8)*gammaPrime/(gammaF*_np.cos(deltaPhi))*_np.sin(alpha), gammaI/gammaF*_np.cos(alpha)]])

def CavityBodyEta(bn = [0], bmn = [1], deltaPhi = 0) :
    bn = _np.array(bn)
    bmn = _np.array(bmn)
    return (bn**2 + bmn**2 + 2*bn*bmn*_np.cos(2*deltaPhi)).sum()

def CavityBodyAlpha(gammaI, gammaF, eta, deltaPhi) :
    return _np.sqrt(eta/8.0)/_np.cos(deltaPhi)*_np.log(gammaF/gammaI)

def CavityFringeTransverseMatrix(gammaI, gammaF, L = 1, inward = True, gammaPrime = None) :
    '''
    Calculates RF cavity fringe transverse matrix.

    Example:

    >>> cf = CavityFringeTransverseMatrix(1,2 inward = True)
        array([[1.   , 0.   ],
               [0.375, 1.   ]])

    +-----------------+---------------------------------------------------------+
    | **Parameters**  | **Description**                                         |
    +-----------------+---------------------------------------------------------+
    | gammaI          | Incoming energy (gamma units do not matter)             |
    +-----------------+---------------------------------------------------------+
    | gammaF          | Outgoing energy (or gamma, units do not matter)         |
    +-----------------+---------------------------------------------------------+
    | inward          | Entrance (True) or exit (False)                         |
    +-----------------+---------------------------------------------------------+
    | gammaPrime      | Energy gain (or gamma, units do not matter). If None    |
    |                 | it is calculated from gammaF - gammaI                   |
    +-----------------+---------------------------------------------------------+
    '''

    if not gammaPrime and L != 0:
        gammaPrime = (gammaF - gammaI)/L

    if inward:
        gammaPrime = -gammaPrime
        gamma = gammaI
    else:
        gamma = gammaF

    return _np.array([[1,                   0],
                      [gammaPrime/(2*gamma),1]])


def CavityPiModeComplete(gammaI, gammaF, L, alpha) :
    gammaPrime = (gammaF-gammaI)/L
    return _np.array([[_np.cos(alpha )-_np.sqrt(2)*_np.sin(alpha),_np.sqrt(8)*gammaI/gammaPrime*_np.sin(alpha)],
                      [-3*gammaPrime/(_np.sqrt(8)*gammaF)*_np.sin(alpha),gammaI/gammaF*(_np.cos(alpha)+_np.sqrt(2)*_np.sin(alpha))]])

def CavityGammaPrime(E0, deltaPhi, q = 1, m0 = 0.511) :
    '''
    Calculates RF cavity energy gain (inits are energy MeV, momentum MeV/c, mass MeV/c**2)

    gammaPrime = q*E0*cos(deltaPhi)/m

    +-----------------+---------------------------------------------------------+
    | **Parameters**  | **Description**                                         |
    +-----------------+---------------------------------------------------------+
    | E0              | Peak field (MV)                                         |
    +-----------------+---------------------------------------------------------+
    | deltaPhi        | Phase difference compared to crest phase                |
    +-----------------+---------------------------------------------------------+
    | charge          | Particle charge in units of e                           |
    +-----------------+---------------------------------------------------------+
    | m0              | Rest mass (MeV/c^2)                                     |
    +-----------------+---------------------------------------------------------+
    '''

    return q*E0*_np.cos(deltaPhi)/m0
