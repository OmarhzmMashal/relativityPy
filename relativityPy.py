from sympy import *
import numpy as np

class RiemannManifold:
    
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
        
        tensor_components = [[[[0] * self.d for i in range(self.d)] for j in range(self.d)] for k in range(self.d)]
        
        cs = self.getChristofellSymbols()

        for p in range(self.d):   
            for i in range(self.d):
                for j in range(self.d):
                    for k in range(self.d):
                        
                        part1 = diff(cs[p][k][i], self.coor[j]) - diff(cs[p][j][i], self.coor[k])
                        
                        # summing over dummy a
                        summation1 = 0
                        for a in range(self.d):
                            summation1 += cs[a][k][i] * cs[p][j][a]
                        
                        summation2 = 0
                        for a in range(self.d):
                            summation2 += cs[a][j][i] * cs[p][k][a]
                                                
                        component = simplify(part1 + summation1 - summation2)
                        tensor_components[p][i][j][k] = component 
                    
                        if string == True and component != 0:
                            self.pRiemannTensorString(p,i,j,k,component)
                            
        return tensor_components


    def pRicciTensorString(self, i,j, val):
        
        R = "R{}{}".format(i,j,k)
        R = R.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
        print(R + " = {}".format(val))
        
        
    def getRicciTensor(self, string=False):
        
        tensor_components = [[0] * self.d for i in range(self.d)]
        
        rt = self.getRiemannTensor(string=True)

        for i in range(self.d):
            for j in range(self.d):

                    # summing over dummy p
                    summation = 0
                    for p in range(self.d):
                        summation += rt[p][i][p][j]

                    summation = simplify(summation)
                    tensor_components[i][j] = summation

                    if string == True and summation != 0:
                        self.pRicciTensorString(i,j,summation)
                            
        return tensor_components
    

# Schwarzschild 
#M, t, theta, r, phi= symbols(" M t theta r phi")
#metric = [[-(1-2*M/r),0,0,0], [0,1/(1-2*M/r),0,0], [0,0,r**2,0], [0,0,0,r**2*sin(theta)**2]]
#coor = [t,r,theta,phi]


# Rindler
#x, t, D = symbols(" x t D")
#metric = [[(x/D)**2, 0] , [0, -1]]
#coor = [t, x]


# Sphere
r, theta, phi = symbols(" r theta phi")
metric = [[1, 0, 0] , [0, r**2 , -1], [0,0, r**2*sin(theta)**2]]
coor = [r, phi, theta]


rm = RiemannManifold(coor, metric)
gammas = rm.getChristofellSymbols(string=True)
riemann_tensor = rm.getRiemannTensor(string=True)
ricci_tensor = rm.getRicciTensor(string=True)

