DDMDD
=====

Distorted-distance models for directional dispersal
(with Bram van Putten, Marco D. Visser, Helene C. Muller-Landau,
  and Patrick A. Jansen)

The published code for fitting distorted-distance models from 
our paper in [MEE](http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2012.00208.x/abstract).


Download the [zip](https://github.com/MarcoDVisser/DDMDD/zipball/master) 
or [tar ball](https://github.com/MarcoDVisser/DDMDD/tarball/master).

Summary
-------

1. Seed and pollen dispersal is often directionally biased, because of the inherent directionality of wind and many other dispersal vectors. Nevertheless, the vast majority of studies of seed and pollen dispersal fit isotropic dispersal kernels to data, implicitly assuming that dispersal is equally likely in all directions.

2. Here, we offer a flexible method for stochastic modelling of directional dispersal data. We show how anisotropic models can be constructed by combining standard dispersal functions with ‘distorted-distance functions’ that transform the circular contour lines of any isotropic dispersal kernel into non-circular shapes. Many existing anisotropic phenomenological models of seed and pollen dispersal are special cases of our framework.

3. We present functional forms for the specific case of elliptic distorted-distance functions, under which contour lines of the seed shadow become non-concentric, nested ellipses, and show how models using these functions can be constructed and parameterized. R-code is provided.

4. We applied the elliptic anisotropic models to characterize seed dispersal in the wind-dispersed Neotropical tree Luehea seemannii (Malvaceae) on Barro Colorado Island, Panama. We used inverse modelling to fit alternative models to data of seed rain into seed traps, the locations of seed traps and adult trees, and tree size.

5. Our anisotropic model performed considerably better than commonly applied isotropic models, revealing that seed dispersal of L. seemannii was strongly directional. The best-fitting model combined a 3-parameter elliptic distorted-distance function that captured the strong directional biases with a 1-parameter exponential dispersal kernel, a 1-parameter negative binomial probability distribution describing the clumping of seed rain and a 1-parameter function relating tree fecundity to tree diameter.

6. The framework presented in this paper enables more flexible and accurate modelling of directional dispersal data. It is applicable not only to studies of seed dispersal, but also to a wide range of other problems in which large numbers of particles disperse from one or more point sources.


![](https://raw.github.com/MarcoDVisser/DDMDD/master/images/DDMDD.png)
