function [tensor_C,tensor_T,iteration] =...
    CTMRG_directional...
    (EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
    tensor_A,tensor_B,tensor_C,tensor_T,varargin)
% Updates CTM tensors within the convergence limit 'epsilon'

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% Determine convergence by a series of epsilon values rather than by two
% values
convergence_repeat_count=0;
converged_iteration=0;

epsilon_env_old=10.0^5;
epsilon_env=10.0^5;
iteration=0;

while convergence_repeat_count<CONVERGENCE_REPEAT_CTMRG &&...
        iteration<ITERATION_MAX_ENV
    
    iteration=iteration+1;
    
    epsilon_env_old=epsilon_env;
    
    tensor_C_old=tensor_C;
    tensor_T_old=tensor_T;
    
    % left CTMRG process
    [tensor_C{1},tensor_C{4},...
        tensor_T{7},tensor_T{8}]=...
        CTMRG_left(tensor_C,tensor_T,tensor_a,tensor_b);
    % right CTMRG process
    [tensor_C{2},tensor_C{3},...
        tensor_T{3},tensor_T{4}]=...
        CTMRG_right(tensor_C,tensor_T,tensor_a,tensor_b);
    % up CTMRG process
    [tensor_C{1},tensor_C{2},...
        tensor_T{1},tensor_T{2}]=...
        CTMRG_up(tensor_C,tensor_T,tensor_a,tensor_b);
    % down CTMRG process
    [tensor_C{3},tensor_C{4},...
        tensor_T{5},tensor_T{6}]=...
        CTMRG_down(tensor_C,tensor_T,tensor_a,tensor_b);
    
    if nargin==8 && strcmp(varargin{1},'sv')
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
    disp('directional CTMRG: Maximum iteration reached.')
end
sprintf('directional CTMRG: iteration=%d',iteration)
end

