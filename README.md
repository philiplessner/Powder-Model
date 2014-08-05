## Powder Model ##

### CV/g of Anodes Made from Cylinders ###

#### CV/g of Cylinders ####

The capacitance of a cylindrical capacitor is given by:

![eq 1](./images/eq1.png)

where ![ep0](./images/epsilon0.png) is the permittivity of free space, ![epr](./images/epsilonr.png) is the relative dielectric constant of the metal oxide,  L is the length of the cylinder, b is the outer radius and a is the inner radius. The thickness of the dielectric is then b-a. The cylinder with dielectric was made by anodizing a cylinder with an orginal radius of R.

The thickness of the dielectric can be written in terms of the native oxide thickness and the formation constant as:

![eq 2](./images/eq2.png)

where ![alpha](./images/alpha.png) is the formation constant in ![umpV](./images/micronpervolt.png)and  b<sub>0</sub>  is the thickness of the native oxide film.

The orginal weight of the cylinder is ![weightcyl](./images/weightcyl.png). From the above equations, we can derive CV/g:

![eq 3](./images/eq3.png)

Using stoichiometry and rearranging, we can derive an expression for R<sup>2</sup> :

![eq 4](./images/eq4.png)

where ![gamma](./images/gamma.png) and X is the number of moles of metal per mole of oxide.

Subsituting the expression for R<sup>2</sup> into the equation for CV/g gives:

![eq 5](./images/eq5.png)

#### Corrected CV/g of Cylinders ####

In the anode made of cylinders, the cylinders touch as shown in this figure (for prismatic cylinders):

![cylinders crossing](./images/cyl-crossing.png)

This reduces the surface area of the anode. If we define the parameter:

![eq 6](./images/eq6.png)

where L is the length of the prism and d is its 'diameter'.

We can write the density of the anode in terms of the density of the metal and n:

![eq 7](./images/eq7.png)

or if we define the fractional density ![fracden](./images/fracden.png) then:

![eq 8](./images/eq8.png)

We can also write the corrected CV/g in terms of n:

![eq 9](./images/eq9.png)

We can then plot the correction factor ![corrfac](./images/corrfac.png) vs. D<sub>f</sub> 

![correction factor](./images/correction.png)

Over the portion of the curve that is nearly linear (say ~0.18 to ~0.5) we can fit a line and then write an expression for CV/g<sub>corrected</sub> vs D<sub>f</sub>:

![eq 10](./images/eq10.png)
