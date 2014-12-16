%% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!

PHYS_DIM=2;
BOND_DIM=3;
CHI=12;

J_ISING=-1.0;
ISING_DIRECTION='X'; % magnetic field along the 'Z'-direction
H_FIELD=0.0;

DELTA_T=10^(-1); % time step for iTEBD and Simple update
EPSILON_ENV=10^(-2); % convergence limit for CTMRG
EPSILON_CELL=10^(-2); % convergence limit for the update_cell
EPSILON_M=10^(-1);
% ITERATION_MAX_SIMPLE=10^5; % maximum number of iteration for 'simple_update'
ITERATION_MAX_ENV=100; % maximum number of iteration for 'full_update_CTM'
CONVERGENCE_REPEAT=3;
%ITERATION_MAX_CELL=100; % maximum number of iteration for 'update_tensor'
ITERATION_MAX_FULL_UPDATE=100;

% Hamiltonian components
pauli_X=[0.0,1.0;1.0,0.0];
pauli_Y=[0.0,-1.0j;1.0j,0.0];
pauli_Z=[1.0,0.0;0.0,-1.0];

% gate tensor generation
gate_tensor=gate_TFIM(J_ISING,H_FIELD,DELTA_T,ISING_DIRECTION);
%% INITIAL STATE & CTM GENERATION (RANDOM; Normalized A,B, CTM)
% initial preparation of (random) tensors A, B, a, b
[tensor_A,tensor_B,tensor_a,tensor_b]=...
    unit_cell_generation(PHYS_DIM,BOND_DIM);
% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |+>_X|+>_X & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% GATE APPLICATION & UPDATE
temp_m_old=0;
temp_m=0;

m_list=[];

convergence_repeat_count=0;
converged_iteration=0;

iteration_full_update=0;
while iteration_full_update<ITERATION_MAX_FULL_UPDATE &&...
        convergence_repeat_count<CONVERGENCE_REPEAT
    
    iteration_full_update=iteration_full_update+1;
    temp_m_old=temp_m;
    
    % A-rightward-B gate application & CTMRG & update tensor A, B
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        gate_application...
        (tensor_A,tensor_B,gate_tensor,'rightward');
    [tensor_C,tensor_T]=...
        full_update_CTM( EPSILON_ENV,CONVERGENCE_REPEAT,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T);
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        update_cell_right(tensor_A,tensor_B,tensor_a,tensor_b,tensor_C,tensor_T,...
        EPSILON_CELL);
    % A-leftward-B gate application & CTMRG & update tensor A, B
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        gate_application...
        (tensor_A,tensor_B,gate_tensor,'leftward');
    [tensor_C,tensor_T]=...
        full_update_CTM( EPSILON_ENV,CONVERGENCE_REPEAT,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T);
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        update_cell_left(tensor_A,tensor_B,tensor_a,tensor_b,tensor_C,tensor_T,...
        EPSILON_CELL);
    % A-upward-B gate application & CTMRG & update tensor A, B
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        gate_application...
        (tensor_A,tensor_B,gate_tensor,'upward');
    [tensor_C,tensor_T]=...
        full_update_CTM( EPSILON_ENV,CONVERGENCE_REPEAT,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T);
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        update_cell_up(tensor_A,tensor_B,tensor_a,tensor_b,tensor_C,tensor_T,...
        EPSILON_CELL);
    % A-downward-B gate application & CTMRG & update tensor A, B
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        gate_application...
        (tensor_A,tensor_B,gate_tensor,'downward');
    [tensor_C,tensor_T]=...
        full_update_CTM( EPSILON_ENV,CONVERGENCE_REPEAT,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T);
    [tensor_A,tensor_B,tensor_a,tensor_b]=...
        update_cell_down(tensor_A,tensor_B,tensor_a,tensor_b,tensor_C,tensor_T,...
        EPSILON_CELL);
    
    temp_m=...
        exp_value(tensor_A,tensor_B,tensor_C,tensor_T,...
        [0.0,1.0;1.0,0.0]);
    m_list=[m_list,temp_m];
    
    % convergence test: convergence in a series of m values.
    if abs(temp_m-temp_m_old)<EPSILON_M
        if convergence_repeat_count==0
            convergence_repeat_count=1;
            converged_iteration=iteration_full_update;
        elseif converged_iteration==iteration_full_update-1
            convergence_repeat_count=convergence_repeat_count+1;
            converged_iteration=iteration_full_update;
        else
            convergence_repeat_count=1;
            converged_iteration=iteration_full_update;
        end
    end
    
end

figure(1);
plot(1:iteration_full_update,m_list,'r')
xlabel('Iteration','FontSize',13,'FontWeight','bold')
ylabel('\left<\sigma_X\right>')
legend('\left<\sigma_X\right>','Location','northwest')