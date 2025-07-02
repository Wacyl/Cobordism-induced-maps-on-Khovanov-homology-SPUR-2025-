# SPUR-2025
Khovanov homology implementation of movie cobordisms based on Sagemath's library


The code in this project heavily relies on the SageMath library's implementation of both Khovanov Homology for knots as well as Chain Complexes and other mathematical objects. 
The goal is to be able to automate cobordism movie moves to compute the induced map on homology which is well-defined up to sign modulo smooth isotopies of surfaces in B^4 bounding links in S^3.

I) CONVENTIONS:
---------------
1) As per Sagemath's documentation for 10.1>=, the planar diagram notation goes counterclockwise (instead of clockwise in older versions)
2) The chain complex who'se homology is invariant under R1/R2/R3 for Links is the one presented in Bar-Natan's paper for a categorification of the Jones polynomial.[1]
3) For ease of implementation, we allow planar diagram notations to not be continuous (eg: [5,10,10,5] represents the unknot with a twist).
4) TYPO : For the _enhanced_states method in the Link object in SageMath, the second tuple represents negative labelled circles and the third tuple **positive** labelled circles.

II) ENHANCED STATES EXPLANATION: 
--------------------------------
Computation of the Khovanov homology for an arbitrary link in Sagemath is done by computing the enhanced Kauffman states of the link, as described in [2] by Oleg Viro.
The chain complexes are built degree-wise in Sagemath, and we will follow their approach for the cobordism movie moves. 
The basis elements are all enhanced states of a specific homological and quantum degree. The maps are described as matrices with +-1 going to the possible adjacent states.
The signs are determined by the usual Bar-Natan's convention [1] (i.e (-1)^{# of 1's before the digit shift}. Notice that this differs from Sage's implementation (number of 1's PAST the digit shift).

III) REIDEMEISTER MOVES:
-----------------------
TODO


IV) BIRTH/DEATH:
---------------
TODO



V) SADDLE MOVES:
---------------
TODO



VI) OTHER DETAILS:
------------------
TODO




VII) REFERENCES:
-------------
- [1] : On Khovanov’s categorification of the Jones polynomial, Dror Bar-Natan.
- [2] : REMARKS ON DEFINITION OF KHOVANOV HOMOLOGY, Oleg Viro.
