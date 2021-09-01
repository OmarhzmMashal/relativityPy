# relativityPy

- Christoffel symbols (Torsion free)
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
```

- Stored Solutions (Schwarzschild, Rindler, & Sphere)
- Spacetime Coordinate Transformation (Galilean & Lorentz)

![alt text](https://github.com/OmarhzmMashal/relativityPy/blob/main/tensors.png)
