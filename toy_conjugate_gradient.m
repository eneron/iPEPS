function [ x ] = conjugate_gradient( matrix_A,vector_b,EPSILON )
% finds minimum using the CG method

% In real life, this corresponds to the principle that one need to try
% different approaches that are "orthogonal" to each other
% rather than just optimal at the current local point
% to deal with given problems (or situations)

r=vector_b;
p=r;
x=0; % initial guess x0=0

k=1;

while norm(r)>EPSILON
    r_old=r;
    p_old=p;
    x_old=x;
    
    alpha=transpose(r_old)*r/(transpose(p)*matrix_A*p);
    x=x_old+alpha*p;
    
    r=r_old-alpha*matrix_A*p_old;
    p=r+(transpose(r)*r/(transpose(r_old)*r_old))*p_old;
    k=k+1;
end

end
