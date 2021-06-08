import scipy as sp
import matplotlib
import matplotlib.pyplot as py
import math

#Some global settings for matplotlib plots
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['font.weight'] = 'bold'

# [Starter code for Problem 8]

# ------------------------------------------------------------
# Solving the bulk with an inhomogeneous iterative method
# ------------------------------------------------------------

#Define some useful objects
class Latticesite:
    def __init__(self,siteNumber):
        self.index = siteNumber
        self.coordinate = []
        self.NNs = []
        self.NNNs = []
        self.potential = 0.0
        self.density_current = 0.0
        self.density_previous = 0.0
    
    def update(self):
        '''Update the density. Remember to account for divergences under iteration'''
        if self.density_current<1.0:
            self.density_previous = self.density_current
        elif self.density_current>1.0:
            self.density_previous = 1.0
        else:
            self.density_previous = 0.0

def iterate(sites,k,mu,beta):
    '''Perform a single iteration of eq. 46 for the particle at site k, given the sitelist sites'''
    nDens = 0.0
    nnDens = 0.0
    for neighbor in sites[k].NNs:
        nDens = nDens + sites[neighbor].density_previous
    for neighbor in sites[k].NNNs:
        nnDens = nnDens + sites[neighbor].density_previous
    return (1-sites[k].density_previous)*sp.exp(beta*(mu+nDens+0.25*nnDens-sites[k].potential)) #Here we assume T,mu,V are all in units of the interaction strength 


#The lattice is a (Lx x Ly square lattice)
Lx = 5
Ly = 6

#Initialize the lattice
def initialize(Lx,Ly):

    def getIndex(Lx,Ly):
        '''Converts from coordinates (x,y) to lattice index'''
        return lambda x,y: Lx*y + x

    nSites = Lx*Ly
    sites = [Latticesite(k) for k in range(nSites)]  
    pos = getIndex(Lx,Ly)     
    for site in sites:
        x,y = [site.index%Lx,site.index//Lx]                                #Get x and y coordinates and from those the coordinates of the neighbors
        nns = set([((x+1)%Lx,y),((x-1)%Lx,y),(x,(y+1)%Ly),(x,(y-1)%Ly)])
        nnns = set([((x+1)%Lx,(y+1)%Ly),((x-1)%Lx,(y-1)%Ly),((x+1)%Lx,(y-1)%Ly),((x-1)%Lx,(y-1)%Ly)])
        #nns = set([((x+1)%Lx,y),((x-1)%Lx,y),(x,(y+1)%Ly),(x,(y-1)%Ly)])    
        #nnns = set([( (x+1)%Lx,(y+1)%Ly ),( (x-1)%Lx, (y-1)%Ly),((x-1)%Lx,(y+1)%Ly),((x+1)%Lx,(y-1)%Ly)])
        site.NNs = [pos(x[0],x[1]) for x in nns]           #Store the neighbor indices as instance variables
        site.NNNs = [pos(x[0],x[1]) for x in nnns]
        site.density_previous = 0.2                        #Initialize the system in the low density limit
    return sites

#Now we iterate the solver until the density is converged
def run(mu,T,Lx,Ly,cTol=10**-8,mixing=0.1,iterMax=1000,show=True):
    'Calculates the density profile at a given mu,T'
    sites = initialize(Lx,Ly)
    convergence = 0.1
    iteration = 0.0
    while (convergence>cTol) and (iteration<iterMax):
        for k in range(len(sites)):
            sites[k].density_current = sites[k].density_previous*(1-mixing) + mixing*iterate(sites,k,mu,1/T) #Calculate new state of the system from the old state of the system
        two_norm = sum([(site.density_current-site.density_previous)**2 for site in sites])
        convergence = math.sqrt(two_norm)
        iteration = iteration +1
        for site in sites:
            site.update()
    
    'Can then return an image of the density profile'
    z = []
    for site in sites:
        z.append(site.density_previous)
    Z = sp.array(z)
    Z = sp.reshape(Z,(Lx,Ly))
    x,y = sp.meshgrid(range(Ly),range(Lx))
    contourlevels = py.MaxNLocator(nbins=10).tick_values(Z.min(), Z.max())
    cmap = py.get_cmap('hot')
    norm = matplotlib.colors.BoundaryNorm(contourlevels, ncolors=cmap.N, clip=True)
    if show==True:
        print(sites[5].density_previous)
        fig,ax = py.subplots()
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.plot()
        raw = ax.pcolormesh(x,y,Z,cmap=cmap,norm=norm)
        fig.colorbar(raw)
        py.show()
    else:
        #print(sites[5].density_previous)
        return py.pcolormesh(x,y,Z,cmap=cmap,norm=norm)

#Run a few examples
for mu in [-1,-1.5,-4]:
    run(mu,1.5,Lx,Ly)
#figs, ax = py.subplots()
#ax = [run(k,0.5,4,4,show=False) for k in [-1,-4]]
