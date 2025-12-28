# Simulation for Binary Response Forecasting under a Factor-Augmented Framework
In this repository, we put all of the simulation code and result.

The result can be found in the document "Updated Result of Simulation.pdf", where I have concluded all of the result in a clear way.

In the folder “Code_for_Example_1”, we include the data-generating processes corresponding to Example 1 in the paper, where the error terms follow a standard normal distribution. In this example, we estimate all coefficients using the standard normal likelihood function.

The three files — simulation_CASE1_DGP_1, simulation_CASE1_DGP_2, and simulation_CASE1_DGP_3 — correspond to the three error-term specifications described in the paper:

· simulation_CASE1_DGP_1 implements the case where the error terms are independent and identically distributed (i.i.d.) standard normal variables.

· simulation_CASE1_DGP_2 implements the case with weak autocorrelation, where the error terms follow an AR(1) process with coefficient 0.3.

· simulation_CASE1_DGP_3 implements the case with strong autocorrelation, where the error terms follow an AR(1) process with coefficient 0.7.

In the folder “Code_for_Example_2”, we include the data-generating processes corresponding to Example 2 in the paper, where the error terms follow a standard logistic distribution. In this example, we estimate all coefficients using the standard logistic likelihood function.

The three files — simulation_CASE2_DGP_1, simulation_CASE2_DGP_2, and simulation_CASE2_DGP_3 — correspond to the three error-term specifications described in the paper:

· simulation_CASE2_DGP_1 implements the case where the error terms are independent and identically distributed (i.i.d.) standard logistic variables.

· simulation_CASE2_DGP_2 implements the case with weak autocorrelation, where the error terms follow an AR(1) process with coefficient 0.3.

· simulation_CASE2_DGP_3 implements the case with strong autocorrelation, where the error terms follow an AR(1) process with coefficient 0.7.

In the folder “Code_for_Example_3”, we include the data-generating processes corresponding to Example 3 in the paper, where the error terms follow a standard logistic distribution. In this example, we estimate all coefficients using the standard normal likelihood function. This is the example we adopt the Quasi-MLE idea.

The three files — simulation_CASE3_DGP_1, simulation_CASE3_DGP_2, and simulation_CASE3_DGP_3 — correspond to the three error-term specifications described in the paper:

· simulation_CASE3_DGP_1 implements the case where the error terms are independent and identically distributed (i.i.d.) standard logistic variables.

· simulation_CASE3_DGP_2 implements the case with weak autocorrelation, where the error terms follow an AR(1) process with coefficient 0.3.

· simulation_CASE3_DGP_3 implements the case with strong autocorrelation, where the error terms follow an AR(1) process with coefficient 0.7.
