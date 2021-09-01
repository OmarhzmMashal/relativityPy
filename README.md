# relativityPy

- Christoffel Symbols (Torsion free)
- Riemann Tensor
- Ricci Tensor

```
# define paramaters
M, t, theta, r, phi= symbols(" M t theta r phi")
metric = [[-(1-2*M/r),0,0,0], [0,1/(1-2*M/r),0,0], [0,0,r**2,0], [0,0,0,r**2*sin(theta)**2]]
coor = [t,r,theta,phi]

# initilaze manifold; coordinates order must match the metric's
rm = RiemannManifold(coor, metric)

# calculate christoffel symbols (Levi-Civita connection)
christoffel_symbols = rm.getChristoffelSymbols(printstr=True)

# calculate riemann tensor
riemann_tensor = rm.getRiemannTensor(printstr=True)

# calculate ricci tensor
ricci_tensor = rm.getRicciTensor(printstr=True)

# calculate ricci scalar
ricci_scalar = rm.getRicciScalar(printstr=True)
```

- Stored Solutions (Schwarzschild, Rindler, & Sphere)
```
# sphere solution
riemann_tensor_for_sphere = RiemannManifoldSolution().sphere(returned_tensor="ricciscalar", printstr=True) 
```

- Spacetime Coordinate Transformation (Galilean & Lorentz)
```
# init trnasformer with a velocity vector
c = 299792458

# input velocity vector [Vx, Vy, Vz]
coor_transformer = CoordinateTransformation(velocity=[0.9*c, 0, 0])

# input spacetime vector [t, x, y, z]
lorentz_transformed_vector = coor_transformer.lorentzT(vector=[1,0,0,0])
galilean_transformed_vector = coor_transformer.galileanT(vector=[1,0,0,0])
```

![alt text](https://github.com/OmarhzmMashal/relativityPy/blob/main/tensors.png)
