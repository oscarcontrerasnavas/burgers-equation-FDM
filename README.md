# Viscous Burgers' Equation
The Burgers' equation or Bateman-Burgers' equation was first introduced in 1915
and later studied in 1948 as a simplification of the Navier-Stock equation to
understand its main mathematical properties such as the tendency of the viscous
term to cero or the nonlinerity of the convective.

Navier-Stokes equation (1):

![](http://latex.codecogs.com/png.latex?%5Crho%5Cdfrac%7B%5Cpartial%5Cmathbf%7Bu%7D%7D%7B%5Cpartial%20t%7D%20&plus;%20%5Crho%28%5Cmathbf%7Bu%7D%5Ccdot%5Cnabla%29%5Cmathbf%7Bu%7D%20%3D%20-%5Cnabla%20p%20&plus;%20%5Cmu%5Cnabla%5E2%5Cmathbf%7Bu%7D%20&plus;%20%5Crho%5Cmathbf%7Bg%7D)

Thus, the Burgers' equation neglects the pressure and gravity terms from
(1)transforms it into a quasi-linear parabolic PDE. We
have written down the expanded equation for 2D Cartesian coordinates below assuming
that ![](http://latex.codecogs.com/png.latex?%5Cinline%20%5Cmathbf%7Bu%7D%3D%5Cleft%3C%20u%28t%2Cx%2Cy%29%2C%20v%28t%2Cx%2Cy%29%20%5Cright%3E) and letting ![](http://latex.codecogs.com/png.latex?%5Cinline%20%5Cnu%3D%5Cmu/%5Crho) :

Burgers' Equation:

![](http://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20t%7D%20&plus;%20u%5Cdfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20x%7D%20&plus;%20v%5Cdfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20y%7D%20%3D%20%5Cnu%20%5Cleft%28%5Cdfrac%7B%5Cpartial%5E2%7Bu%7D%7D%7B%5Cpartial%7Bx%5E2%7D%7D%20&plus;%20%5Cdfrac%7B%5Cpartial%5E2%7Bu%7D%7D%7B%5Cpartial%7By%5E2%7D%7D%5Cright%29)

![](http://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20t%7D%20&plus;%20u%5Cdfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20x%7D%20&plus;%20v%5Cdfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20y%7D%20%3D%20%5Cnu%20%5Cleft%28%5Cdfrac%7B%5Cpartial%5E2%7Bv%7D%7D%7B%5Cpartial%7Bx%5E2%7D%7D%20&plus;%20%5Cdfrac%7B%5Cpartial%5E2%7Bv%7D%7D%7B%5Cpartial%7By%5E2%7D%7D%5Cright%29)

And this is the equation we are solving for two particular situations: first, we
solve it for a high value of the kinetic viscosity and, second, we drop that
value making the convective source predominates until the *shocks* appears. All of
that using a **python** script and the libraries `numpy` and `matplotlib` within
a beautiful animated plot.

<hr />

For more information about the mathematical background and how it was coded,
please follow my blog [here](https://oscarcontrerasnavas.github.io/python/burgers-equation-with-finite-difference-method.html)