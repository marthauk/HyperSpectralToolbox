function Kn = inverse_cordic_growth_constant(niter)
    % By looping, it is seen that Kn will converge against 0.6073. niter =7
    % is decent.
    Kn = 1/prod(sqrt(1+2.^(-2*(0:double(niter)-1))))
   end