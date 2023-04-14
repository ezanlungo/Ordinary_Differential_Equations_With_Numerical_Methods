# Ordinary_Differential_Equations_With_Numerical_Methods

Systems with more than one chemical reaction associated in series-parallel can be modeled using ordinary differential equations (ODEs). 
In this case, the following system of ODEs is considered:

𝑑𝑥1/dt = −𝑘1 ∗ 𝑥1 + 𝑘2 ∗ 𝑥2 ∗ 𝑥3
𝑑𝑥2/𝑑𝑡 = 𝑘1 ∗ 𝑥1 − 𝑘2 ∗ 𝑥2 ∗ 𝑥3 − 𝑘3 ∗ 𝑥2^2
𝑑𝑥3/𝑑𝑡 = 𝑘3 ∗ 𝑥2^2

Normally, solving these systems is cumbersome, but with numerical methods, it can be made simpler while maintaining a high level of accuracy. 
In this study, the system of ODEs will be solved using different numerical methods with various step sizes, 
and then compared against an analytical solution to draw conclusions.

Vectors are constructed using the lines of code shown in the attached .py file, with one vector for each step size used. 
In this study, four step sizes were used, corresponding to 10, 100, 500, and 1000 subdivisions in the working interval. 
The following equations were used for each numerical method:

-Euler
-Taylor
-Runge-Kutta 4
-Predictor-Corrector

All the obtained graphs have in common that as more subdivisions are added, the curves become smoother and visually appear to provide a better fit.

Next, the solutions of each method are contrasted against an analytical solution. 
Since an analytical solution is not available, the solution obtained from a high-precision method such as Runge-Kutta 4 with a very small step size (corresponding to 20000 subdivisions) is designated as the reference solution.

The results show that both the Euler and Taylor methods increase in accuracy as the number of subdivisions in the working interval is increased and the step size is decreased. 
On the other hand, the Predictor-Corrector method increases in error as the step size is increased, although the differences are very small, resulting in similar results compared to the previous two methods.

It is observed that the use of the Predictor-Corrector method with a large step size produces a lower global error compared to Taylor or Euler with a large number of subdivisions, as shown in the last graph. 
It should be noted that the total trajectory of the curves is not being analyzed, but rather the global error, which corresponds to the difference in the last position of the vectors.

Finally, the results of the Runge-Kutta 4 method are the best by a noticeable margin, due to two main factors. 
The first reason is that it is indeed the method that offers the best numerical solution for our system of ODEs, 
and the second reason is that it is logical for the result to fit better, as the analytical solution is calculated using the same method but with a much smaller step size.

It is worth noting that the accuracy of the Runge-Kutta 4 method slightly improves with an increase in the number of subdivisions or a decrease in the step size.

In conclusion, one should always be aware of the resources available for computation. 
If the computer allows for it, it is always better to work with high-precision methods such as RK-4 or Adam-Bashforth of order 4 or 5. 
However, in this study, it has been shown that mathematically simpler methods such as Euler or Taylor can provide acceptable results with lower resource demands and computational operations.
