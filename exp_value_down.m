function [ exp_value,exp_value_A,exp_value_B ] =...
    exp_value_down( tensor_A,tensor_B,tensor_Env,Observable )
% calculates the expectation value of the observable.

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% environment tensor E_a for the reduced tensor a
tensor_E_A=ncon({tensor_Env,tensor_b},...
    {[1,2,3,4,5,6,-1,-2,-4,-5,-6,-8],[-3,1,2,3,-7,4,5,6]});
rho_A=ncon({tensor_E_A,tensor_A,conj(tensor_A)},...
    {[1,2,3,4,5,6,7,8],[-1,1,2,3,4],[-2,5,6,7,8]});
% environment tensor E_b for the reduced tensor b
tensor_E_B=ncon({tensor_Env,tensor_a},...
    {[-2,-3,-4,-6,-7,-8,1,2,3,4,5,6],[1,2,-1,3,4,5,-5,6]});
rho_B=ncon({tensor_E_B,tensor_B,conj(tensor_B)},...
    {[1,2,3,4,5,6,7,8],[-1,1,2,3,4],[-2,5,6,7,8]});

if abs(trace(rho_A))<10^(-12) || abs(trace(rho_B))<10^(-12)
    disp('tr(rho_A)=0 or tr(rho_B)=0: CHECK THEM.')
    exp_value_A=0;
    exp_value_B=0;
    exp_value=0;
else
    exp_value_A=ncon({rho_A,Observable},{[1,2],[2,1]})/...
        trace(rho_A);
    exp_value_B=ncon({rho_B,Observable},{[1,2],[2,1]})/...
        trace(rho_B);
    
    exp_value=(exp_value_A+exp_value_B)/2;
end
end

