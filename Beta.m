function pen=Beta(U)
I=[1;1;1;1;1];
pen=transpose(max(0,U-I))*max(0,U-I) + transpose(max(0,-U))*max(0,-U);




