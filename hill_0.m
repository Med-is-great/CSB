function bis= hill_function(gamma,Ce,ec50)
a=(Ce)^gamma;
b=a+((ec50)^(gamma+0.5));
bis=10*a/b;
end