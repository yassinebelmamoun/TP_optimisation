function [fU,dfU]=f_df_quasi_newton(U,B,S)
fU = transpose(U) * S * U - transpose(B) * U;
dfU = 2 * S * U - B;

