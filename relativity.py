from sympy import *
import numpy as np
import math 

class RiemannManifold: 
   
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

    
    def getChristoffelSymbols(self, printstr=False):
        
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
                    
                    if printstr == True and gamma != 0:
                        self.pGammaString(i,j,k,gamma)               
        
        return christofell_symbols

    
    def pRiemannTensorString(self,p,i,j,k,val):
        
        R = "R {}".format(p)
        R = R.translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
        R = R + "{}{}{}".format(i,j,k)
        R = R.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
        print(R + " = {}".format(val))

        
    def getRiemannTensor(self, printstr=False):
        
        tensor_components = [[[[0] * self.d for i in range(self.d)] for j in range(self.d)] for k in range(self.d)]
        
        cs = self.getChristoffelSymbols()

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
                    
                        if printstr == True and component != 0:
                            self.pRiemannTensorString(p,i,j,k,component)
                            
        return tensor_components


    def pRicciTensorString(self, i,j, val):
        
        R = "R{}{}".format(i,j,k)
        R = R.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
        print(R + " = {}".format(val))
        
        
    def getRicciTensor(self, printstr=False):
        
        tensor_components = [[0] * self.d for i in range(self.d)]
        
        rt = self.getRiemannTensor()

        for i in range(self.d):
            for j in range(self.d):

                # summing over dummy p
                summation = 0
                for p in range(self.d):
                    summation += rt[p][i][p][j]

                summation = simplify(summation)
                tensor_components[i][j] = summation

                if printstr == True and summation != 0:
                    self.pRicciTensorString(i,j,summation)
                            
        return tensor_components
    
    
    def getRicciScalar(self, printstr=False):
        
        ricci_scalar = 0
        
        rt = self.getRicciTensor()

        for i in range(self.d):
            for j in range(self.d):

                ricci_scalar += Matrix(self.metric).inv(method="LU")[i,j] * rt[i][j] 


        if printstr == True:
            R = "R = {}".format(ricci_scalar)
            print(R)  
                          
        return ricci_scalar
    
    
class RiemannManifoldSolution:
    
    def __init__(self):
        pass
    
    
    def schwarzschild(self, returned_tensor="riemann", printstr=False):

        M, t, theta, r, phi= symbols(" M t theta r phi")
        metric = [[-(1-2*M/r),0,0,0], [0,1/(1-2*M/r),0,0], [0,0,r**2,0], [0,0,0,r**2*sin(theta)**2]]
        coor = [t,r,theta,phi]
        
        rm = RiemannManifold(coor, metric)

        if returned_tensor == "riemann":
            return rm.getRiemannTensor(printstr)
        
        elif returned_tensor == "gammas":
            return  rm.getChristoffelSymbols(printstr)
            
        elif returned_tensor == "ricci":
            return rm.getRicciTensor(printstr)

        elif returned_tensor=="ricciscalar":
            return rm.getRicciScalar(printstr)
        

    def rindler(self, returned_tensor="riemann", printstr=False):

        # where D = speed of light squared / proper acceleration
        x, t, D = symbols(" x t D")
        metric = [[(x/D)**2, 0] , [0, -1]]
        coor = [t, x]
        
        rm = RiemannManifold(coor, metric)

        if returned_tensor == "riemann":
            return rm.getRiemannTensor(printstr)
        
        elif returned_tensor == "gammas":
            return  rm.getChristoffelSymbols(printstr)
            
        elif returned_tensor == "ricci":
            return rm.getRicciTensor(printstr)

        elif returned_tensor=="ricciscalar":
            return rm.getRicciScalar(printstr)
        
    
    def sphere(self, returned_tensor="riemann", printstr=False):

        # where D = speed of light squared / proper acceleration
        r, theta, phi = symbols(" r theta phi")
        metric = [[1, 0, 0] , [0, r**2 , -1], [0,0, r**2*sin(theta)**2]]
        coor = [r, phi, theta]
        
        rm = RiemannManifold(coor, metric)

        if returned_tensor == "riemann":
            return rm.getRiemannTensor(printstr)
        
        elif returned_tensor == "gammas":
            return  rm.getChristoffelSymbols(printstr)
            
        elif returned_tensor == "ricci":
            return rm.getRicciTensor(printstr)

        elif returned_tensor=="ricciscalar":
            return rm.getRicciScalar(printstr)


    def euclidean(self, returned_tensor="riemann", printstr=False):

        # where D = speed of light squared / proper acceleration
        x, y, z = symbols(" x y z")
        metric = [[1, 0, 0] , [0, 1 , 0], [0,0, 1]]
        coor = [x, y, z]
        
        rm = RiemannManifold(coor, metric)

        if returned_tensor == "riemann":
            return rm.getRiemannTensor(printstr)
        
        elif returned_tensor == "gammas":
            return  rm.getChristoffelSymbols(printstr)
            
        elif returned_tensor == "ricci":
            return rm.getRicciTensor(printstr) 
        
        elif returned_tensor=="ricciscalar":
            return rm.getRicciScalar(printstr)


class CoordinateTransformation:
    
    def __init__(self, velocity):
        
        self.v1 = velocity[1]
        self.v2 = velocity[2]
        self.v3 = velocity[3]
    
    def lorentzT(self, vector):
        
        c = 299792458
        V = math.sqrt(self.v1**2 + self.v2**2 + self.v3**2)
        lf = 1/math.sqrt(1-(V**2/c**2))
        
        B1 = self.v1/c
        B2 = self.v2/c
        B3 = self.v3/c
                
        lB = [[lf,-lf*self.v1,-lf*self.v2,-lf*self.v3],
              [-lf*self.v1, 1+(lf-1)*(self.v1**2/V**2), (lf-1)*(self.v1*self.v2/V**2), (lf-1)*(self.v1*self.v3/V**2)],
              [-lf*self.v2, (lf-1)*(self.v1*self.v2/V**2), 1+(lf-1)*(self.v2**2/V**2), (lf-1)*(self.v2*self.v3/V**2)],
              [-lf*self.v3, (lf-1)*(self.v1*self.v3/V**2), (lf-1)*(self.v2*self.v3/V**2), 1+(lf-1)*(self.v3**2/V**2)]] 
           
        return Matrix(lB) * Matrix(vector) 
         
    def galileanT(self, vector, vector_field="spacetime"):
        
        gB = [[1, 0, 0, 0],
              [self.v1, 1, 0, 0],
              [self.v2, 0, 1, 0],
              [self.v3, 0, 0, 1]] 
        
        return Matrix(gB) * Matrix(vector) 
