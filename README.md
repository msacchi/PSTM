# PSTM
## Kirchhoff Migration (PSTM) tests 

Testing forward and adjoint operators for Prestack Time Migration (PSTM)

The subsurface image is given by $m(z,x)$.  The image is an estimator of the subsurface reflectivity. The data are given by $d(t,k)$ where $k$ is the trace number. Receiver and source positions are $xs(k)$ and $xg(k)$ which are deployed at $z=0$. 

The code *main.jl* shows how to use the two operators

1) Demigration $d =W L m$  where $W$ is time convolution and $L$ is the de-migration operator:

  `d  = Operator_Conv(Operator_PSTM(m,false; Param_PSTM...), false; Param_Conv...)`

2) Migration $m' = L'W' d$ where $L'$ is the migration operator and $W'$ is the adjoint of time domain 
convilution equivalent to time domain cross-correlation:

`mig = Operator_PSTM(Operator_Conv(d,true; Param_Conv...), true; Param_PSTM...)`

3) In Least-squares Migration (LSM) we estimate the reflectivity by solving
the following problem

$$m_{ls} = argmin{J}$$ with cost $$J \| WL m - d\|_2^2 + \mu \| m\|_2^2$$
 The cost function is minimized using the Conjugate Gradient methodimplicit which requires the  
operators $WL$ and $L'W'$. In the code, this is given by 


`m_ls,J = ConjugateGradients(d,[Operator_Conv, Operator_PSTM],[Param_Conv, Param_PSTM]; Niter=20,mu=0.01)`


*main.jl* runs an example where I use the "invere problem crime" to model data and then retrieve the model via CGLS

*dot_product_test.jl* checks that $L$ and $L'$ pass the dot product test

The results are shown below. First, we show the data which is computed via the demigration operator $WLm$ 

![image](figure1.png)

Now, I show the migrated image which is computed via $L'W'd$

![image](figure2.png)

And finally, the least-sqaures migration after $20$ CG iterations:

![image](figure3.png)

I also show the cost $J$ vs iteration number. Clearly, the algorithn must converge because I am solving
a linear problem by minimizing a quadratic cost. 

![image](figure4.png)