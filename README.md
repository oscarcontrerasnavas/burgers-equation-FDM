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

## The Finite Differences

The finite difference method or *FDM* for short is the simplest and more
accessible of the *numeric method* employed to find a solution of systems with partial
differential equations such it is the Burgers' that we are already familiar
with. The method implies to change such a system for an easier algebraic one using
for that approximations based on the **Taylor's Expansion** of a function called
finite difference, in a very similar way of how to the trapezoids method for
integrals does.

![Finite Differences](https://upload.wikimedia.org/wikipedia/commons/thumb/9/90/Finite_difference_method.svg/320px-Finite_difference_method.svg.png)

<small class="label-text">
<em>Figure 1. Finite Differences Stencil and its geometrical representation.</em>
</small>

As the image above shows, three different schemes could be formulated at the
point ![](http://latex.codecogs.com/png.latex?x):

**Forward Differences**:

![](http://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%7Bf%7D%7D%7B%5Cpartial%7Bx%7D%7D%5Capprox%20%5Cdfrac%7Bf%28x&plus;%5CDelta%7Bx%7D%29-f%28x%29%7D%7B%5CDelta%7Bx%7D%7D)

**Backward Differences**:

![](http://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%7Bf%7D%7D%7B%5Cpartial%7Bx%7D%7D%5Capprox%20%5Cdfrac%7Bf%28x%29-f%28x-%5CDelta%7Bx%7D%29%7D%7B%5CDelta%7Bx%7D%7D)

**Central Differences**:

![](http://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%7Bf%7D%7D%7B%5Cpartial%7Bx%7D%7D%5Capprox%20%5Cdfrac%7Bf%28x&plus;%5CDelta%7Bx%7D%29-f%28x-%5CDelta%7Bx%7D%29%7D%7B2%5CDelta%7Bx%7D%7D)

Also, we have a formula to approximate the second partial derivatives, you can
demonstrate its veracity replacing with the previous differences.

**Second Order Central Differences**:

![](http://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%5E2%7Bf%7D%7D%7B%5Cpartial%7Bx%5E2%7D%7D%5Capprox%20%5Cdfrac%7Bf%28x&plus;%5CDelta%7Bx%7D%29-2f%28x%29&plus;f%28x-%5CDelta%7Bx%7D%29%7D%7B%5CDelta%7Bx%7D%5E2%7D)

## Discretizing the Domain. The Numerical Scheme

Assuming we have a PDE over the region ![](http://latex.codecogs.com/png.latex?%5Cinline%20%5C%7B%20%5Cmathbf%7Bx%7D%20%5Cin%20%28a%2C%20b%29%5C%7D) and we
divided it into finite numbers of points *Nx* in a grid where every node can be
written with a index *ni* then, we have now a way of represent our region as a
consecutive number of points from *i=0* to *i=Nx*.

With our already set up grid it is necessary to define our numerical model in
the discretized domain as well as which numerical scheme we are going to use for
each of the derivatives. Following there is the scheme for the *x-coordinate* with
time index *n*, x-direction index *i* and y-direction index *j*.

**Forward Differences for partial derivative with respect to time**:

$$
\dfrac{\partial{u}}{\partial{t}} =
\dfrac{u^{n+1}_{i,j}-u^{n}_{i,j}}{\Delta{t}}
$$

**Backward Differences for partial derivative with respect to x**:

$$
\dfrac{\partial{u}}{\partial{x}} =
\dfrac{u^{n}_{i,j}-u^{n}_{i-1,j}}{\Delta{x}}
$$

**Backward Differences for partial derivative with respect to y**:

$$
\dfrac{\partial{u}}{\partial{y}} =
\dfrac{u^{n}_{i,j}-u^{n}_{i,j-1}}{\Delta{y}}
$$

**2nd Order Central Differences for 2nd order partial derivative with respect to x**:

$$
\dfrac{\partial^2{u}}{\partial{x}^2} =
\dfrac{u^{n}_{i+1,j}-2u^{n}_{i,j}+u^{n}_{i-1,j}}{\Delta{x}^2}
$$

**2nd Order Central Differences for 2nd order partial derivative with respect to y**:

$$
\dfrac{\partial^2{u}}{\partial{y}^2} =
\dfrac{u^{n}_{i,j+1}-2u^{n}_{i,j}+u^{n}_{i,j-1}}{\Delta{y}^2}
$$

<div class='alert alert-info'>
<strong>Note:</strong> Don't worry if you feel lost at this point. You probably
need to keep reading or go to the reference at the bottom of the page for detailed
information.
</div>

## Algebraic Model
{: #anchor-4}

Still here? Good. From here, we will transpose our PDE's —one for the x-direction and
other for the y-direction— and solved for the next step in time, if you read
carefully, you would have noted that next time step is our only unknown value.
This scheme is known as the explicit scheme, easy for working with but it also
have some stability issues of which you can read later at the references
section.

With the previous settings, we get:

$$
\dfrac{u^{n+1}_{i,j}-u^{n}_{i,j}}{\Delta{t}} +
u^{n}_{i,j} \dfrac{u^{n}_{i,j}-u^{n}_{i-1,j}}{\Delta{x}} +
v^{n}_{i,j} \dfrac{u^{n}_{i,j}-u^{n}_{i,j-1}}{\Delta{y}}
$$

$$=$$

$$
\nu \dfrac{u^{n}_{i+1,j}-2u^{n}_{i,j}+u^{n}_{i-1,j}}{\Delta{x}^2} +
\nu \dfrac{u^{n}_{i,j+1}-2u^{n}_{i,j}+u^{n}_{i,j-1}}{\Delta{y}^2}
$$

And,

$$
\dfrac{v^{n+1}_{i,j}-u^{n}_{i,j}}{\Delta{t}} +
u^{n}_{i,j} \dfrac{v^{n}_{i,j}-v^{n}_{i-1,j}}{\Delta{x}} +
v^{n}_{i,j} \dfrac{v^{n}_{i,j}-v^{n}_{i-1,j}}{\Delta{y}}
$$

$$=$$

$$
\nu \dfrac{v^{n}_{i+1,j}-2v^{n}_{i,j}+v^{n}_{i-1,j}}{\Delta{x}^2} +
\nu \dfrac{v^{n}_{i,j+1}-2v^{n}_{i,j}+v^{n}_{i,j-1}}{\Delta{y}^2}
$$

**Now we rearrange these equations for $$u, v$$**:

$$
u^{n+1}_{i,j} = u^{n}_{i,j} - u^{n}_{i,j}
\dfrac{\Delta{t}}{\Delta{x}}(u^{n}_{i,j}-u^{n}_{i-1,j}) \\ - v^{n}_{i,j}
\dfrac{\Delta{t}}{\Delta{y}}(u^{n}_{i,j}-u^{n}_{i,j-1})
+\nu\dfrac{\Delta{t}}{\Delta{x}^2}(u^{n}_{i+1,j}-2u^{n}_{i,j}+u^{n}_{i-1,j}) \\
+\nu\dfrac{\Delta{t}}{\Delta{y}^2}(u^{n}_{i,j+1}-2u^{n}_{i,j}+u^{n}_{i,j-1})
\tag{9}$$

And,

$$
v^{n+1}_{i,j} = v^{n}_{i,j} - u^{n}_{i,j}
\dfrac{\Delta{t}}{\Delta{x}}(v^{n}_{i,j}-v^{n}_{i-1,j}) \\ - v^{n}_{i,j}
\dfrac{\Delta{t}}{\Delta{y}}(v^{n}_{i,j}-v^{n}_{i,j-1})
+\nu\dfrac{\Delta{t}}{\Delta{x}^2}(v^{n}_{i+1,j}-2v^{n}_{i,j}+v^{n}_{i-1,j}) \\
+\nu\dfrac{\Delta{t}}{\Delta{y}^2}(v^{n}_{i,j+1}-2v^{n}_{i,j}+v^{n}_{i,j-1})
\tag{10}$$
