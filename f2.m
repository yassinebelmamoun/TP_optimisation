function fU=f2(U,B,S)
       
fU = transpose(U) * S * U + transpose(U) * exp(U);
