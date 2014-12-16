% environment tensor E
tensor_Env=tensor_Env_right(...
        tensor_Cell{1},tensor_Cell{2},tensor_C,tensor_T);

tensor_a=ncon({tensor_Cell{1},conj(tensor_Cell{1})},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_Cell{2},conj(tensor_Cell{2})},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});    
    
% environment tensor E_a for the reduced tensor a
tensor_E_a=ncon({tensor_Env,tensor_b},...
    {[-1,-5,1,2,3,4,5,6,-3,-7,-2,-6],[1,-4,5,3,2,-8,6,4]});
tensor_rho_A=ncon({tensor_E_a,tensor_Cell{1},conj(tensor_Cell{1})},...
    {[1,2,3,4,5,6,7,8],[-1,1,2,3,4],[-2,5,6,7,8]});
tensor_rho_A=tensor_rho_A/trace(tensor_rho_A)
eig(tensor_rho_A)

% environment tensor E_b for the reduced tensor b
tensor_E_b=ncon({tensor_Env,tensor_a},...
    {[1,2,-1,-5,-4,-8,-3,-7,5,6,3,4],[1,3,5,-2,2,4,6,-6]});
tensor_rho_B=ncon({tensor_E_b,tensor_Cell{2},conj(tensor_Cell{2})},...
    {[1,2,3,4,5,6,7,8],[-1,1,2,3,4],[-2,5,6,7,8]});
tensor_rho_B=tensor_rho_B/trace(tensor_rho_B)
eig(tensor_rho_B)

if trace(tensor_rho_A)<10^(-12)||trace(tensor_rho_B)<10^(-12)
    sprintf('tr(rho_A)=%.3g, tr(rho_B)=.3g',...
        trace(tensor_E_A),trace(tensor_E_B))
end

% sprintf('%s= %.3g, %.3g','E_A',eig(tensor_E_A))
% sprintf('%s= %.3g, %.3g','E_B',eig(tensor_E_B))