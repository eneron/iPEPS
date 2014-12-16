function [ exp_value,exp_value_A,exp_value_B ] =...
    exp_value_right( tensor_A,tensor_B,tensor_Env,Observable )
% calculates the expectation value of the observable.

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% environment tensor E_a for the reduced tensor a
tensor_E_A=ncon({tensor_Env,tensor_b},...
    {[-1,-5,-2,-6,-3,-7,1,2,3,4,5,6],[1,-4,3,5,2,-8,4,6]});
rho_A=ncon({tensor_E_A,tensor_A,conj(tensor_A)},...
    {[1,2,3,4,5,6,7,8],[-1,1,2,3,4],[-2,5,6,7,8]});
% environment tensor E_b for the reduced tensor b
tensor_E_B=ncon({tensor_Env,tensor_a},...
    {[1,2,3,4,5,6,-1,-5,-3,-7,-4,-8],[1,3,5,-2,2,4,6,-6]});
rho_B=ncon({tensor_E_B,tensor_B,conj(tensor_B)},...
    {[1,2,3,4,5,6,7,8],[-1,1,2,3,4],[-2,5,6,7,8]});

if abs(trace(rho_A))<10^(-12) || abs(trace(rho_B))<10^(-12)
    error('tr(rho_A)=0 or tr(rho_B)=0: CHECK THEM.')
%     exp_value_A=0;
%     exp_value_B=0;
%     exp_value=0;
else
    exp_value_A=ncon({rho_A,Observable},{[1,2],[2,1]})/...
        trace(rho_A);
    exp_value_B=ncon({rho_B,Observable},{[1,2],[2,1]})/...
        trace(rho_B);
    
    exp_value=(exp_value_A+exp_value_B)/2;
end
end

