% 20141201: Wrong understanding on the TEBD algorithm; gate application
% goes after CTMRG and the construction of tensor_Env.

% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!
path('./conjugate_gradient_Poblano',path)

PHYS_DIM=2;
BOND_DIM=2;
CHI=20;

J_ISING=-1.0; % H= J_ISING X_i*X_j +H_FIELD/4 (Z_i+Z_j)
ISING_DIRECTION='X'; % magnetic field along the 'Z'-direction

H_FIELD=-3.3; % 3.0 For reference, QMC predicts 3.044

DELTA_T=10^(-1); % time step for iTEBD and Simple update
EPSILON_ENV=10^(-2); % convergence limit for CTMRG
EPSILON_CELL=10^(-10); % convergence limit for the update_cell
EPSILON_M=10^(-5);
% ITERATION_MAX_SIMPLE=10^5; % maximum number of iteration for 'simple_update'
ITERATION_MAX_ENV=1000; % maximum number of iteration for 'CTMRG_directional'
CONVERGENCE_REPEAT_CTMRG=2;
CONVERGENCE_REPEAT_TEBD=2;
% ITERATION_MAX_CELL=100; % maximum number of iteration for 'update_tensor'
ITERATION_MAX_TEBD=10;

% Skip the input-error checks by the tensor contraction module 'ncon':
% global ncon_skipCheckInputs;
% ncon_skipCheckInputs = true;

% Turn off the warning 'Complex values of exp_values bulabula~'
warning('off','all');
%% iTEBD iteration
tic

% RANDOM INITIAL STATE & RANDOM CTM GENERATION
[tensor_A,tensor_B,~,~]=unit_cell_generation(PHYS_DIM,BOND_DIM);

% INITIAL STATE |+>_X|+>_X & RANDOM CTM GENERATION
% tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
% tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

% tensor_A_old=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
% tensor_B_old=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

% INITIAL STATE |0>_Z|0>_Z & RANDOM CTM GENERATION
% tensor_A=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
% tensor_B=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
% tensor_A(1,:,:,:,:)=1;
% tensor_B(1,:,:,:,:)=1;

% tensor_a=ncon({tensor_A,conj(tensor_A)},...
%     {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
% tensor_b=ncon({tensor_B,conj(tensor_B)},...
%     {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);

m_list=[];

% gate tensor generation
gate_tensor=gate_TFIM(J_ISING,H_FIELD,DELTA_T,ISING_DIRECTION);

temp_m_old=0;
temp_m=0;

convergence_repeat_count=0;
converged_iteration=0;

iteration_TEBD=0;
while iteration_TEBD<ITERATION_MAX_TEBD &&...
        convergence_repeat_count<CONVERGENCE_REPEAT_TEBD
    multiWaitbar('iTEBD processing...',iteration_TEBD/ITERATION_MAX_TEBD,...
        'Color','g');
    iteration_TEBD=iteration_TEBD+1;
    temp_m_old=temp_m;
    
    % GATE APPLICATION & UPDATE
    % A-rightward-B gate application & CTMRG & update tensor A, B
    sprintf('H_FIELD=%.4f, iTEBD iteration=%d, Rightward',...
        H_FIELD,iteration_TEBD)
    [tensor_C,tensor_T]=...
        CTMRG_directional(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T,'sv');
    tensor_Env=tensor_Env_right(tensor_A,tensor_B,tensor_C,tensor_T);
    %             [tensor_A,tensor_B]=...
    %                 gate_application...
    %                 (tensor_A,tensor_B,gate_tensor,'rightward');
    [tensor_A,tensor_B]=...
        update_cell_right(tensor_A,tensor_B,gate_tensor,...
        tensor_Env,EPSILON_CELL);
    % A-leftward-B gate application & CTMRG & update tensor A, B
    sprintf('H_FIELD=%.4f, iTEBD iteration=%d, Leftward',...
        H_FIELD,iteration_TEBD)
    [tensor_C,tensor_T]=...
        CTMRG_directional(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T,'sv');
    tensor_Env=tensor_Env_left(tensor_A,tensor_B,tensor_C,tensor_T);
    %             [tensor_A,tensor_B]=...
    %                 gate_application...
    %                 (tensor_A,tensor_B,gate_tensor,'leftward');
    [tensor_A,tensor_B]=...
        update_cell_left(tensor_A,tensor_B,gate_tensor,...
        tensor_Env,EPSILON_CELL);
    % A-upward-B gate application & CTMRG & update tensor A, B
    sprintf('H_FIELD=%.4f, iTEBD iteration=%d, Upward',...
        H_FIELD,iteration_TEBD)
    [tensor_C,tensor_T]=...
        CTMRG_directional(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T,'sv');
    tensor_Env=tensor_Env_up(tensor_A,tensor_B,tensor_C,tensor_T);
    %             [tensor_A,tensor_B]=...
    %                 gate_application...
    %                 (tensor_A,tensor_B,gate_tensor,'upward');
    [tensor_A,tensor_B]=...
        update_cell_up(tensor_A,tensor_B,gate_tensor,...
        tensor_Env,EPSILON_CELL);
    % A-downward-B gate application & CTMRG & update tensor A, B
    sprintf('H_FIELD=%.4f, iTEBD iteration=%d, Downward',...
        H_FIELD,iteration_TEBD)
    [tensor_C,tensor_T]=...
        CTMRG_directional(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T,'sv');
    tensor_Env=tensor_Env_down(tensor_A,tensor_B,tensor_C,tensor_T);
    %             [tensor_A,tensor_B]=...
    %                 gate_application...
    %                 (tensor_A,tensor_B,gate_tensor,'downward');
    [tensor_A,tensor_B]=...
        update_cell_down(tensor_A,tensor_B,gate_tensor,...
        tensor_Env,EPSILON_CELL);
    
    % Calculation of the new environment and the expectation value
    [tensor_C,tensor_T]=...
        CTMRG_directional(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
        tensor_A,tensor_B,tensor_C,tensor_T,'sv');
    tensor_Env=tensor_Env_right(tensor_A,tensor_B,tensor_C,tensor_T);
    temp_m=exp_value_right(tensor_A,tensor_B,tensor_Env,...
        [0.0,1.0;1.0,0.0]);
    
    if ~isfinite(temp_m)
        disp('NaN or Inf occurred in the magnetization: CHECK IT')
        iteration_TEBD=ITERATION_MAX_TEBD;
    end
    
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
    
    m_list=[m_list,temp_m];
end
multiWaitbar('iTEBD processing...','close');
elapsed_time=toc
%% Plot the average m values
% Data folder creation
% DATA_FOLDER=data_folder_create('DATA_TEBD_dCTMRG');
% save('20141201_temp')
figure(1)
% figure('visible','off')
plot(1:iteration_TEBD,m_list,'r-x')
xlabel('Iteration TEBD')
ylabel('<X>_{avg}')
legend('<X>_{dCTMRG}','Location','northwest')
title(sprintf(...
    'TFIM 2D quantum [dCTMRG]\n H FIELD=%.2g, BOND DIM=%d, CHI=%d, epsm=%.2g,\n crepeat CTMRG=%d, crepeat TEBD=%d, max iter TEBD=%d, max iter CTMRG=%d',...
    H_FIELD,BOND_DIM,CHI,EPSILON_M,CONVERGENCE_REPEAT_CTMRG,CONVERGENCE_REPEAT_TEBD,ITERATION_MAX_TEBD,ITERATION_MAX_ENV))
% print('-depsc',strcat('./',DATA_FOLDER,...
%     sprintf(...
%     '/H_%.2g_epsm_%.2g_CHI_%d_crepeat_CTMRG_%d_crepeat_TEBD_%d_iter_CTMRG_%d_iter_TEBD_%d.eps',...
%     H_FIELD,EPSILON_M,CHI,CONVERGENCE_REPEAT_CTMRG,CONVERGENCE_REPEAT_TEBD,...
%     ITERATION_MAX_ENV,ITERATION_MAX_TEBD)));
%
% save(strcat('./',DATA_FOLDER,sprintf(...
%     '/H_%.2g_epsm_%.2g_CHI_%d_crepeat_CTMRG_%d_crepeat_TEBD_%d_iter_CTMRG_%d_iter_TEBD_%d.mat',...
%     H_FIELD,EPSILON_M,CHI,CONVERGENCE_REPEAT_CTMRG,CONVERGENCE_REPEAT_TEBD,...
%     ITERATION_MAX_ENV,ITERATION_MAX_TEBD)));

% % Send 'job-completed' notification email
% setpref('Internet','E_mail','hojoon.kim@kias.re.kr');
% setpref('Internet','SMTP_Server','mail.kias.re.kr');
% setpref('Internet','SMTP_Username','hojoon.kim');
% setpref('Internet','SMTP_Password','kias7321!');
% sendmail('hojoon.kim@kias.re.kr','[abacus4] iTEBD completed',...
%     'iTEBD completed at Abacus4.');

load chirp
sound(y,Fs)
clear('y','Fs')
