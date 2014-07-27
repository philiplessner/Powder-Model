## Powder Model ##

### CV/g of Anodes Made from Cylinders ###

#### CV/g of Cylinders ####

The capacitance of a cylindrical capacitor is given by:

![eq 1](./images/eq1.png)

where ![ep0](./images/epsilon0.png) is the permittivity of free space, ![epr](./images/epsilonr.png) is the relative dielectric constant of the metal oxide,  L is the length of the cylinder, b is the outer radius and a is the inner radius. The thickness of the dielectric is then b-a. The cylinder with dielectric was made by anodizing a cylinder with an orginal radius of R.

The thickness of the dielectric can be written in terms of the native oxide thickness and the formation constant as:

![eq 2](./images/eq2.png)

where \\( \alpha \\) is the formation constant \\( (\mu m/V) \\) and \\( b_0 \\) is the thickness of the native oxide film.

The orginal weight of the cylinder is \\( \pi R^2 L \rho_{M} \\). From the above equations, we can derive CV/g:

![eq 3](./images/eq3.png)

Using stoichiometry and rearranging, we can derive an expression for \\( R^2 \\):

![eq 4](./images/eq4.png)

where \\( \gamma = \frac{\rho_{oxide} {MW}_{M}}{\rho_{M} {MW}_{oxide}} \\) and \\( X \\) is the number of moles of metal per mole of oxide.

Subsituting the expression for \\( R^2 \\) into the equation for CV/g gives:

![eq 5](./images/eq5.png)

#### Corrected CV/g of Cylinders ####

In the anode made of cylinders, the cylinders touch as shown in this figure (for prismatic cylinders):

![cylinders crossing](https://github.com/philiplessner/Powder-Model/blob/master/images/cyl-crossing.png)
This reduces the surface area of the anode. If we define the parameter:

\\[ n = \frac{L}{d} \\]

where \\( L \\) is the length of the prism and \\( d \\) is its 'diameter'.

We can write the density of the anode in terms of the density of the metal and \\( n \\):

\\[ D_s = \left(\frac{3}{n^2} - \frac{2}{n^3}\right)\rho_M \\]

or if we define the fractional density \\( D_f = \frac{D_s}{\rho_M} \\), then:

\\[ D_f = \left(\frac{3}{n^2} - \frac{2}{n^3}\right) \\]

We can also write the corrected \\( \frac{CV}{g} \\) in terms of n:

\\[ \frac{CV}{g}_{corrected} = \frac{(n-1)}{n}\frac{CV}{g}_{ucorrected} \\]

We can then plot the correction factor \\( \frac{n-1}{n} \\) vs. \\( D_f \\):

![correction factor](https://github.com/philiplessner/Powder-Model/blob/master/images/correction.svg)
Over the portion of the curve that is nearly linear (say ~0.18 to ~0.5) we can fit a line and then write an expression for \\( \frac{CV}{g}_{corrected} \\) vs \\( D_f \\):

\\[ \frac{CV}{g}_{corrected} = -0.712D_f + 0.8529 \\]
