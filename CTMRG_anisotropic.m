function [tensor_C,tensor_T,iteration] =...
    CTMRG_anisotropic( EPSILON_ENV,CONVERGENCE_REPEAT,ITERATION_MAX_ENV,...
    tensor_Cell,tensor_C,tensor_T,varargin)
% Updates CTM tensors within the convergence limit 'epsilon'

% Determine convergence by a series of epsilon values rather than by two
% values
convergence_repeat_count=0;
converged_iteration=0;

epsilon_env_old=10.0^5;
epsilon_env=10.0^5;
iteration=0;

while convergence_repeat_count<CONVERGENCE_REPEAT &&...
        iteration<ITERATION_MAX_ENV
    
    iteration=iteration+1;
    
    epsilon_env_old=epsilon_env;
    
    tensor_C_old=tensor_C;
    tensor_T_old=tensor_T;
    
    % horizontal CTMRG process
    [tensor_C,tensor_T]=...
        CTMRG_horizontal(tensor_Cell,tensor_C,tensor_T);
    % vertical CTMRG process
    [tensor_C,tensor_T]=...
        CTMRG_vertical(tensor_Cell,tensor_C,tensor_T);
    
    if nargin==7 && strcmp(varargin{1},'sv')
        epsilon_env=...
            sum(abs(svd(tensor_C{1})-svd(tensor_C_old{1}))+...
            abs(svd(tensor_C{2})-svd(tensor_C_old{2}))+...
            abs(svd(tensor_C{3})-svd(tensor_C_old{3}))+...
            abs(svd(tensor_C{4})-svd(tensor_C_old{4})));
    else
        epsilon_env=...
            convergence_CTM(tensor_C_old,tensor_T_old,tensor_C,tensor_T,...
            tensor_A,tensor_B);
    end
    
    % convergence test: convergence in a series of m values.
    if abs(epsilon_env-epsilon_env_old)<EPSILON_ENV
        if convergence_repeat_count==0
            convergence_repeat_count=1;
            converged_iteration=iteration;
        elseif converged_iteration==iteration-1
            convergence_repeat_count=convergence_repeat_count+1;
            converged_iteration=iteration;
        else
            convergence_repeat_count=1;
            converged_iteration=iteration;
        end
    end
    
end
if iteration>=ITERATION_MAX_ENV
    disp('Maximum iteration reached.')
end
sprintf('CTMRG anisotropic: iteration=%d',iteration)
end

