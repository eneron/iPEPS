% BACK-UP: adding new functions environment_tensor_left, etc.

% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!
% path('.\conjugate_gradient_Poblano')

PHYS_DIM=2;
BOND_DIM=2;
CHI=20;

J_ISING=-1.0;
ISING_DIRECTION='X'; % magnetic field along the 'Z'-direction
H_FIELD_MAX=3.4;
H_FIELD_MIN=2.9; % 3.0 For reference, QMC predicts 3.044
H_FIELD_STEP=0.01;
H_LIST=H_FIELD_MIN:H_FIELD_STEP:H_FIELD_MAX;

DELTA_T=10^(-1); % time step for iTEBD and Simple update
EPSILON_ENV=10^(-5); % convergence limit for CTMRG
EPSILON_CELL=10^(-8); % convergence limit for the update_cell
EPSILON_M=10^(-7);
% ITERATION_MAX_SIMPLE=10^5; % maximum number of iteration for 'simple_update'
ITERATION_MAX_ENV=100; % maximum number of iteration for 'full_update_CTM'
CONVERGENCE_REPEAT=3;
% ITERATION_MAX_CELL=100; % maximum number of iteration for 'update_tensor'
ITERATION_MAX_TEBD=1000;

% Hamiltonian components
pauli_X=[0.0,1.0;1.0,0.0];
pauli_Y=[0.0,-1.0j;1.0j,0.0];
pauli_Z=[1.0,0.0;0.0,-1.0];

% Skip the input-error checks by the tensor contraction module 'ncon':
% global ncon_skipCheckInputs;
% ncon_skipCheckInputs = true;

% Skip the input-error checks by the tensor contraction module 'ncon':
global ncon_skipCheckInputs;
ncon_skipCheckInputs = true;

% Turn off the warning 'Complex values of exp_values bulabula~'
warning('off','all');
%% iTEBD iteration

tic

% INITIAL STATE |+>_X|+>_X & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);

m_list=zeros(1,round((H_FIELD_MAX-H_FIELD_MIN)/H_FIELD_STEP)+1);
iteration_TEBD_list=[];

for H_FIELD=H_FIELD_MIN:H_FIELD_STEP:H_FIELD_MAX
    sprintf('H_FIELD=%.4f',H_FIELD)
    % gate tensor generation
    gate_tensor=gate_TFIM(J_ISING,H_FIELD,DELTA_T,ISING_DIRECTION);
    
    temp_m_old=0;
    temp_m=0;
    
    convergence_repeat_count=0;
    converged_iteration=0;
    
    iteration_TEBD=0;
    while iteration_TEBD<ITERATION_MAX_TEBD &&...
            convergence_repeat_count<CONVERGENCE_REPEAT
        
        iteration_TEBD=iteration_TEBD+1;
        temp_m_old=temp_m;
        
        % GATE APPLICATION & UPDATE
        % A-rightward-B gate application & CTMRG & update tensor A, B
        sprintf('iTEBD iteration=%d, Rightward',iteration_TEBD)
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
        sprintf('iTEBD iteration=%d, Leftward',iteration_TEBD)
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
        sprintf('iTEBD iteration=%d, Upward',iteration_TEBD)
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
        sprintf('iTEBD iteration=%d, Downward',iteration_TEBD)
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
        
        % convergence test: convergence in a series of m values.
        if abs(temp_m-temp_m_old)<EPSILON_M
            if convergence_repeat_count==0
                convergence_repeat_count=1;
                converged_iteration=iteration_TEBD;
            elseif converged_iteration==iteration_TEBD-1
                convergence_repeat_count=convergence_repeat_count+1;
                converged_iteration=iteration_TEBD;
            else
                convergence_repeat_count=1;
                converged_iteration=iteration_TEBD;
            end
        end
        
    end
    m_list(round((H_FIELD-H_FIELD_MIN)/H_FIELD_STEP)+1)=temp_m;
    iteration_TEBD_list=[iteration_TEBD_list,iteration_TEBD];
end
elapsed_time=toc
%% Plot the average m values
% figure('visible','off')
figure('visible','off')
plot(H_LIST,m_list,'r-x',H_LIST,iteration_TEBD_list,'b-x')
xlabel('H_Z FIELD')
ylabel('<X>_{avg}')
legend('<X>_{dCTMRG}','iter_TEBD/Max_iter','Location','northwest')
title(sprintf(...
    'TFIM 2D quantum [dCTMRG] epsm:\n%.2g, BOND DIM:%d, CHI:%d, crepeat: %d, max iter (TEBD)=%d',...
    EPSILON_M,BOND_DIM,CHI,CONVERGENCE_REPEAT,ITERATION_MAX_TEBD))
% ylim([0,1.1])
print('-depsc',...
    sprintf('H_%.2g_epsm_%.2g_CHI_%d_crepeat_%d_iter_CTMRG_%d_iter_TEBD_%d.eps',...
    H_FIELD,EPSILON_M,CHI,CONVERGENCE_REPEAT,ITERATION_MAX_ENV,ITERATION_MAX_TEBD));
save(sprintf('H_%.2g_epsm_%.2g_CHI_%d_crepeat_%d_iter_CTMRG_%d_iter_TEBD_%d.mat',...
    H_FIELD,EPSILON_M,CHI,CONVERGENCE_REPEAT,ITERATION_MAX_ENV,ITERATION_MAX_TEBD));

% Send 'job-completed' notification email
setpref('Internet','E_mail','hojoon.kim@kias.re.kr');
setpref('Internet','SMTP_Server','mail.kias.re.kr');
setpref('Internet','SMTP_Username','hojoon.kim');
setpref('Internet','SMTP_Password','kias7321!');
sendmail('hojoon.kim@kias.re.kr','[abacus4] iTEBD completed',...
    'iTEBD completed at Abacus4.');
