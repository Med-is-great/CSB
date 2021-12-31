function bis= hill_function(u,gamma,psi,ec50)
a=(u*psi)^gamma;
b=a+((ec50)^gamma);
bis=10*a/b;
end