These codes concern the computation and study of the spectrum of operators derived from quantum mechanics. More
specifically, Schrödinger operators with periodic potentials, particularly relevant in solid-state physics for describing the dynamics of
electrons in a crystal through the probability wave function. The eigenvalues
then correspond to the different energy levels that can exist within the crystal.

The spectrum of these operators is well known in the case of real potentials.
In this case, the spectrum coincides with a union of bands on the real axis. 
The properties of these bands vary significantly depending on the dimensionality
of the problem. In the one-dimensional case, these bands become increasingly
longer, and the gaps between them decrease rapidly, but the number of gaps
remains infinite. In the case of a 2D or 3D lattice, the band structure is more
complex. The bands can overlap, and the gaps between them may be reduced or
even disappear entirely, leading to a continuous spectrum rather than distinct
bands. Conjectured in the 1930s by the physicists A. Sommerfeld and H. Bethe,
it is now proved since 2007 that the number of gaps in all higher dimensions than
one is always finite. 
For the non-selfadjoint case (complex potential), there are significantly fewer results available. Developing numerical methods for the rapid
calculation of such spectra is important for formulating or disproving new conjectures.

In the 1D case, introducing a complex potential does not radically
change the nature of the spectrum, as it still consists of bands (although they are displaced from the real axis) with boundaries that are relatively
easy to determine. However, in the 2D case, introducing a complex potential profoundly changes the nature of the spectrum. It
is no longer composed of bands but rather of surfaces, and the contours of these
surfaces are very challenging to delineate.

-----------------------------------------------------------------------------------------------------------------------------

The operators we study are : 

1D:	H1D = -d^2\dx^2 + V(x)			
	
	V is 2pi periodic   

2D:	H2D = -d^2\dx^2 - d^2\dy^2 + V(x,y) 	
	
	V(.,y) is 2pi periodic
	V(x,.) is 2pi periodic

-----------------------------------------------------------------------------------------------------------------------------

DESCRIPTION OF THE CODES
To compute the spectrum of the Schrödinger operator with periodic potentials, we use the Floquet-Bloch transform, which provides a decomposition of this spectrum into a family of sub-problems formulated over one period of the potential, with boundary conditions depending on a parameter k. This family is uncountable (indexed by a continuous k). Numerically, only a finite number of k values can be taken.


pot1D -> To enter V(x)

pot2D -> To enter V(x,y)

main1D -> Compute the spectrum of H1D by using a regular discretization of the half Brillouin zone (0,0.5)

Remark: For the case of a real V, we could only use k=0 and k=0.5 to get the opposite endpoints of each band and then connect them with a continuous line on the real axis. This is not possible for a complex V, where the two endpoints are connected by a curve that extends beyond the real axis, and whose shape is unknown.

main2D -> Compute the spectrum of H2D by using a set of n_rand k values randomly sampled with the Brillouin zone (0.1)*(0.1)

Remark: It also draws in dark the eigenvalues obtained by variating k along the edges of a quarter of the Brillouin Zone (square (0,0)->(0.5,0.5))

This gives some boundaries of the spectrum... but not all of them. Why?

In some special cases like:
- V real 
- V(x,y) = V(-x,y) and V(x,y) = V(x,-y)
- V(x,y) = V(x)+V(y)

it can be showed that using the quarter of the Brillouin Zone is enough to get all the spectrum.

The eigenvalues are continuous fonctions of the parameter k=(k1,k2) so the image of this square by a given eigenvalue has to be a connected component of the spectrum. However this does not implies that some of the inside of this square can not map outside the image of its boundaries...

-----------------------------------------------------------------------------------------------------------------------------
For more infos:

These codes are part of those developed during a research internship supervised by Mr. Marletta at Cardiff University. The report of this internship is provided in this folder, as well as the slides of a presentation given in French at my engineering school ENSTA Paris.
