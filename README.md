# PSTM
## Kirchhoff Migration (PSTM) tests 

Testing forward and adjoint operators for PSTM. 

The subsurface image is given by $m(z,x)$.  The image an estimator of the subsurface reflectivity. The data are given by $d(t,k)$ where $k$ is trace number. Receiver and source positions are $xs(k)$ and $xg(k)$ which are deployed at $z=0$. 

The code *main.jl* shows how to use the two operators

1) Demigration $d = L m$
2) Migration $m' = L' d$ 
3) How to use CGLS to solve the LS Migration problem for PSTM. In other words
CGLS is used to minimize the cost function J. The latter requires to use implicit form 
operators L and 'L

$$J = \| L m - d\|_2^2 + \mu \| m\|_2^2$$



*main.jl* runs an example where I use the "invere problem crime" to model data and then retrieve the model via CGLS

*dot_product_test.jl* checks that $L$ and $L'$ pass the dot product test

Results should look like:

![image](figure1.png)

![image](figure2.png)

![image](figure3.png)

![image](figure4.png)