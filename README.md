# relativityPy

- Christoffel Symbols (Torsion free)
- Riemann Curvature Tensor (rank 4 tensor)
- Ricci Curvature Tensor (rank 2 tensor)
- Ricci Curvature Scalar (rank 0 tensor)

```
# define paramaters
M, t, theta, r, phi= symbols(" M t theta r phi")
metric = [[-(1-2*M/r),0,0,0], [0,1/(1-2*M/r),0,0], [0,0,r**2,0], [0,0,0,r**2*sin(theta)**2]]
coor = [t,r,theta,phi]

# init manifold; coordinates order must match the metric's
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
riemann_tensor_for_sphere = RiemannManifoldSolution().sphere(returned_tensor="riemann", printstr=True) 
```

![alt text](https://github.com/OmarhzmMashal/relativityPy/blob/main/tensors.png)
