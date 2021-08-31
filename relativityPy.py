from sympy import *
import numpy as np

class ReimannManifold:
    
    # set hyperspace and its variables 
    
    def __init__(self, coor, metric):
    
        self.coor = coor
        self.metric = metric
        self.d = len(coor)
    
    
    def pGammaString(self,u,v,i,val):
        gamma = "Γ {}".format(u)
        gamma = gamma.translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
        gamma = gamma + "{}{}".format(v,i)
        gamma = gamma.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
        print(gamma + " = {}".format(val))

    
    def getChristofellSymbols(self, string=False):
        # N x N x N 
        christofell_symbols = [[[0] * self.d for i in range(self.d)] for j in range(self.d)]
        
        for i in range(self.d):
            for j in range(self.d):
                for k in range(self.d):
                    
                    gamma = 0
                    for d in range(self.d):
                        gamma += 0.5 * Matrix(self.metric).inv(method="LU")[i,d] * (diff(self.metric[k][d], self.coor[j])\
                                                                               + diff(self.metric[d][j],self.coor[k])\
                                                                               - diff(self.metric[j][k],self.coor[d]))

                    gamma = simplify(gamma)
        
                    christofell_symbols[i][j][k] = gamma
                    
                    if string == True and gamma != 0:
                        self.pGammaString(i,j,k,gamma)               
        
        return christofell_symbols

    
    def pRiemannTensorString(self,p,i,j,k,val):
        R = "R {}".format(p)
        R = R.translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
        R = R + "{}{}{}".format(i,j,k)
        R = R.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
        print(R + " = {}".format(val))

        
    def getRiemannTensor(self, string=False):
        
        # N x N x N x N
        tensor_components = [[[[0] * self.d for i in range(self.d)] for j in range(self.d)] for k in range(self.d)]
        
        cs = self.getChristofellSymbols()

        for p in range(self.d):
           
            for i in range(self.d):
                for j in range(self.d):
                    for k in range(self.d):
                        
                        part1 = diff(cs[p][k][i], self.coor[j]) - diff(cs[p][j][i], self.coor[k])
                        
                        # summing over dummy a
                        s1 = 0
                        for a in range(self.d):
                            s1+=cs[a][k][i] * cs[p][j][a]
                        
                        s2=0
                        for a in range(self.d):
                            s2+=cs[a][j][i] * cs[p][k][a]
                                                
                        component = simplify(part1 + s1 - s2)
                        
                        tensor_components[p][i][j][k] = component 
                    
                        if string == True and component != 0:
                            self.pRiemannTensorString(p,i,j,k,component)
                            
        return tensor_components


    def pRicciTensorString(self, i,j, val):
        R = "R{}{}".format(i,j,k)
        R = R.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
        print(R + " = {}".format(val))
        
        
    def getRicciTensor(self, string=False):
        
        # N x N 
        tensor_components = [[0] * self.d for i in range(self.d)]
        
        rt = self.getRiemannTensor(string=True)

        for i in range(self.d):
            for j in range(self.d):

                    # summing over dummy a
                    s = 0
                    for p in range(self.d):
                        s += rt[p][i][p][j]

                    s = simplify(s)
                    tensor_components[i][j] = s

                    if string == True and s != 0:
                        self.pRicciTensorString(i,j,s)
                            
        return tensor_components
    
    
    def getLowerIndecies(self, tensor, n):
        
        # N x N 
        tensor_components = [[0] * self.d for i in range(self.d)]
        
        rt = self.getRiemannTensor(string=True)

        for i in range(self.d):
            for j in range(self.d):

                    # summing over dummy a
                    s = 0
                    for p in range(self.d):
                        s += rt[p][i][p][j]

                    s = simplify(s)
                    tensor_components[i][j] = s

                    if string == True and s != 0:
                        self.pRicciTensorString(i,j,s)
                            
        return tensor_components


    
    def vis():
        
        rspace = 10.0 #space
        
        anglespace = 2*pi
        
        nx = 1000 #discretize
        
        dx = xspace/nx #increment
        dtheta = anglespace/nx
        dphi = anglespace/nx
        
        r = 2
        
        reimann_tensor = rm.getRiemannTensor(string=True)
        
        r_vals = np.linspace(0,rspace,nx) 
        theta_vals = np.linspace(0,anglespace,nx) 
        phi_vals = np.linspace(0,anglespace,nx) 
        
        
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)
        
        
        lam_x = lambdify([r, phi, theta], x, modules=['numpy'])
        lam_y = lambdify([r, phi, theta], x, modules=['numpy'])
        lam_z = lambdify([r, phi, theta], x, modules=['numpy'])
        
        print(lam_x, lam_y , lam_z)
        
        


# swchizler
#M, t, theta, r, phi= symbols(" M t theta r phi")
#metric = [[-(1-2*M/r),0,0,0], [0,1/(1-2*M/r),0,0], [0,0,r**2,0], [0,0,0,r**2*sin(theta)**2]]
#coor = [t,r,theta,phi]


# rindler
#x, t, D = symbols(" x t D")
#metric = [[(x/D)**2, 0] , [0, -1]]
#coor = [t, x]


# spher
r, theta, phi = symbols(" r theta phi")
metric = [[1, 0, 0] , [0, r**2 , -1], [0,0, r**2*sin(theta)**2]]
coor = [r, phi, theta]


# spher
#r, theta = symbols(" r theta")
#metric = [[r**2 , 0], [0, r**2*sin(theta)**2]]
#coor = [r, theta]


rm = ReimannManifold(coor, metric)
#gammas = rm.getChristofellSymbols(string=True)
#ricci_tensor = rm.getRicciTensor(string=True)
reimann_tensor = rm.getRiemannTensor(string=True)


e0=0
e1=0
e2=0
n=len(coor)
for i in range(n):
    for j in range(n):
        for k in range(n):
            e0 += reimann_tensor[0][i][j][k]

for i in range(n):
    for j in range(n):
        for k in range(n):
            e1 += reimann_tensor[1][i][j][k]

            
e0 = simplify(e0)
e1 = simplify(e1)
e2 = simplify(e2)

print(e0)
print(e1)
print(e2)