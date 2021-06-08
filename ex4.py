import scipy as sp
import matplotlib.pyplot as py

'Part 1: Solve the simple two-body problem with one-body fixed'

#For motion in a plane we have four coupled differential equations
#for the (x,y) co-ordinate and (vx,vx) velocities/momenta
#To solve this problem we implement the following algorithm:
#1.) Initialize with (x0,y0), (vx0,vy0)
#2.) We have set up the RK4 equations which evalutate kn = fn(x0,y0,vx0,vy0)
#    need to now map the f's across various inputs to get the kn at each step

class body:
    ''' Defines a massive body object with various properties as defined in __init__
        I would add in here additional functions under init which could be used to get kinetic
        energy etc. '''
    def __init__(self,initState,M,F):
        self.state = initState      #Current state-point
        self.M = M                  #Mass
        self.F = F                  #Force acting on body
    def getR(self):
        '''Return position co-ordinates'''
        return (self.state[1],self.state[2])
    def getV(self):
        '''Return velocity co-ordinates'''
        return (self.state[3],self.state[4])
    def propagate(self,F,dx):
        '''Perform one RK4 update for a single body in motion around another'''
        while True:
            var = self.state[1:]  
            K1 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K1)]
            K2 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K2)]
            K3 = [f(var) for f in F]
            var = [s+dx*k for (s,k) in zip(var,K3)]
            K4 = [f(var) for f in F]
            #Update the state. First list is the timestep, second list is the RK update. Does not yet update internally
            yield [self.state[0]+dx] + [x+dx*(k1+2*(k2+k3)+k4)/6 for (x,k1,k2,k3,k4) in zip(self.state[1:],K1,K2,K3,K4)],dx

    def propagateMulti(self,bodies,dx):
        '''Perform one RK4 update for a single body in motion around another'''
        def forces(bodies):
            '''Update the force on the body due to other bodies'''    
            dSquared = lambda X:-1*sum(b.M()*sum((x0-x)**2 for x0,x in zip(X,b.getR())) for b in bodies)         
            return [lambda var:var[k]*F for k in range(1,len(self.getR()+1))]       
        F = [lambda var:var[2],lambda var: var[3]].extend(forces(bodies))
        while True:
            var = self.state[1:]  
            K1 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K1)]
            K2 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K2)]
            K3 = [f(var) for f in F]
            var = [s+dx*k for (s,k) in zip(var,K3)]
            K4 = [f(var) for f in F]
            #Update the state. First list is the timestep, second list is the RK update. Does not yet update internally
            yield [self.state[0]+dx] + [x+dx*(k1+2*(k2+k3)+k4)/6 for (x,k1,k2,k3,k4) in zip(self.state[1:],K1,K2,K3,K4)],dx
'Part 1: Simple orbit around a massive object'
G = 1       #I have used units MG=1
M = 1/G
F = [lambda var:var[2],lambda var: var[3], lambda var:-G*M*var[0]/((var[0]**2+var[1]**2)**1.5),lambda var:-G*M*var[1]/((var[0]**2+var[1]**2)**1.5)]

moon = body([0,0,7,0.36,0],1/G)         #Set up initial state of the moon. These parameters give a decently stable orbit

newState = moon.propagate(F,10**-4)     #newState is a "propagate" function. Calling next(newState) will generate the next state point

trajectory = []                         #Store all position co-ordinates. Could also choose to store only a few
dev = []                                #Honestly, should probably use scipy arrays here
k=0
while len(dev)<5:
    newPoint,dx = next(newState)           #Calculate next state
    if (moon.getR()[0]<0 and moon.getR()[1]>0) and newPoint[1]>0:
        dev.append(newPoint[1:3])
    moon.state = newPoint               #Propagate one timestep and update the state
    if k%10==0:                         #Store each 10th co-ordinate point
        trajectory.append(moon.getR())
    k+=1
trajectory = sp.array(trajectory)
#dev = [(d-7)*100/7 for d in dev]
#py.plot(dev)
#py.ylabel('Percent deviation ')
py.plot(trajectory[:,0],trajectory[:,1])
py.plot(0,0,'o')
py.show()
