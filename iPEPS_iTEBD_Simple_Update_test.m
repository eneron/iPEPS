%% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!

PHYS_DIM=2;
BOND_DIM=2;
CHI=20;

J_ISING=-1.0;
ISING_DIRECTION='X';
H_FIELD=-2.7;

DELTA_T=10^(-1); % time step for iTEBD and Simple update
EPSILON_ENV=10^(-3); % convergence limit for CTMRG
EPSILON_BULK=10^(-4); % convergence limit for the simple update
ITERATION_MAX_SIMPLE=10^3; % maximum number of iteration for 'simple_update'
ITERATION_MAX_ENV=500; % maximum number of iteration for 'CTMRG_directional'

CONVERGENCE_REPEAT_CTMRG=2;

% gate tensor generation
gate_tensor=gate_TFIM(J_ISING,H_FIELD,DELTA_T,ISING_DIRECTION);

% From Coroboz's talk, I learned that exploiting simple updates for a few
% times initially will speed up convergence because those steps would
% prepare a good initial state.
%% INITIAL STATE & CTM GENERATION (RANDOM; Normalized A,B, CTM)
% initial preparation of (random) tensors A, B, a, b
[Gamma_A,Gamma_B,lambda]=...
    unit_cell_generation_simple_update(PHYS_DIM,BOND_DIM);
% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% SIMPLE UPDATE
epsilon_bulk=10^5;
iteration=0;
while epsilon_bulk>=EPSILON_BULK && iteration<ITERATION_MAX_SIMPLE
    [Gamma_A,Gamma_B,lambda]=...
        simple_update(Gamma_A,Gamma_B,lambda,gate_tensor);
    iteration=iteration+1;
    multiWaitbar('Simple updating...',iteration/ITERATION_MAX_SIMPLE);
end
if iteration==ITERATION_MAX_SIMPLE
    disp('Simple Update: Maximum iteration reached.')
end
multiWaitbar('Simple updating...','close');
sprintf('Simple Update:\n iteration=%d, epsilon=%.5g',...
    iteration,epsilon_bulk)

tensor_A=ncon({Gamma_A,sqrt(lambda{1}),sqrt(lambda{2}),...
    sqrt(lambda{3}),sqrt(lambda{4})},...
    {[-1,1,2,3,4,],[-2,1],[-3,2],[3,-4],[4,-5]});
tensor_A=tensor_A/max(abs(tensor_A(:)));
% tensor_A(abs(tensor_A)<CUTOFF)=0;
tensor_B=ncon({Gamma_B,sqrt(lambda{3}),sqrt(lambda{4}),...
    sqrt(lambda{1}),sqrt(lambda{2})},...
    {[-1,1,2,3,4,],[-2,1],[-3,2],[3,-4],[4,-5]});
tensor_B=tensor_B/max(abs(tensor_B(:)));
% tensor_B(abs(tensor_B)<CUTOFF)=0;

[tensor_C,tensor_T]=...
    CTMRG_directional(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
    tensor_A,tensor_B,tensor_C,tensor_T,'sv');

% tensor_Cell={tensor_A,tensor_B,tensor_B,tensor_A};
% [tensor_C,tensor_T]=...
%         CTMRG_anisotropic(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
%         tensor_Cell,tensor_C,tensor_T,'sv');

tensor_Env=tensor_Env_right(tensor_A,tensor_B,tensor_C,tensor_T);
temp_m=exp_value_right(tensor_A,tensor_B,tensor_Env,[0.0,1.0;1.0,0.0])

% % save workspace variables in a file:
% save(strcat('simple_update_',ISING_DIRECTION,...
%     sprintf(...
%     '_J_%.2f_H_%.2f_delT_%.2e_eps_%.2e_d_%d_D_%d_CHI_%d_iter_%.1e.mat',...
%     J_ISING,H_FIELD,DELTA_T,EPSILON_BULK,PHYS_DIM,BOND_DIM,CHI,iteration)))
% % load workspace variables in a file:
% clear
% load(strcat('simple_update_',ISING_DIRECTION,...
%     sprintf('_d_%d_D_%d_CHI_%d_iter_%d.mat',...
%     PHYS_DIM,BOND_DIM,CHI,10^5)))

% % initial preparation of a cell of (random) tensors C, T
% [tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%
% [tensor_C,tensor_T]=...
%     full_update_CTM( EPSILON_ENV,ITERATION_MAX_ENV,...
%     tensor_A,tensor_B,tensor_C,tensor_T);
%
% shortcut_exp_values

load chirp
sound(y,Fs)
clear('y','Fs')