% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!
path('./conjugate_gradient_Poblano',path)

PHYS_DIM=2;
BOND_DIM=2;
CHI=20;

J_ISING=-1.0;
ISING_DIRECTION='X'; % magnetic field along the 'Z'-direction
H_FIELD_MAX=3.2;
H_FIELD_MIN=2.8; % 3.0 For reference, QMC predicts 3.044
H_FIELD_STEP=0.1;
H_LIST=H_FIELD_MIN:H_FIELD_STEP:H_FIELD_MAX;

DELTA_T=10^(-1); % time step for iTEBD and Simple update
EPSILON_ENV=10^(-3); % convergence limit for CTMRG
EPSILON_CELL=10^(-10); % convergence limit for the update_cell
EPSILON_M=10^(-5);
% ITERATION_MAX_SIMPLE=10^5; % maximum number of iteration for 'simple_update'
ITERATION_MAX_ENV=500; % maximum number of iteration for 'CTMRG_directional'
CONVERGENCE_REPEAT_CTMRG=2;
CONVERGENCE_REPEAT_TEBD=2;
% ITERATION_MAX_CELL=100; % maximum number of iteration for 'update_tensor'
ITERATION_MAX_TEBD=100;

% Skip the input-error checks by the tensor contraction module 'ncon':
global ncon_skipCheckInputs;
ncon_skipCheckInputs = true;

% Turn off the warning 'Complex values of exp_values bulabula~'
warning('off','all');

sprintf('ISING_DIRECTION=%c, J_ISING=%.2f, H_FIELD=%.4f:%.4f:%.4f',...
    ISING_DIRECTION,J_ISING,H_FIELD_MAX,H_FIELD_STEP,H_FIELD_MIN)
sprintf('DELTA_T=%.1g, EPSILON_M=%.1g, EPSILON_ENV=%.1g, EPSILON_CELL=%.1g',...
    DELTA_T,EPSILON_M,EPSILON_ENV,EPSILON_CELL)
sprintf('ITERATION_MAX_ENV=%d, CONVERGENCE_REPEAT_CTMRG=%d, ITERATION_MAX_TEBD=%d, CONVERGENCE_REPEAT_TEBD=%d',...
    ITERATION_MAX_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_TEBD,CONVERGENCE_REPEAT_TEBD)
%% iTEBD iteration
tic
repeat=1;
REPEAT_MAX=1;
while repeat<=REPEAT_MAX
    % RANDOM INITIAL STATE & RANDOM CTM GENERATION
    % [tensor_Cell{1:4}]=unit_cell_generation(PHYS_DIM,BOND_DIM);
    
    % INITIAL STATE |+>_X|+>_X & RANDOM CTM GENERATION
    % initial preparation of trivial tensors A, B, a, b
    tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
    tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
    tensor_Cell{1}=tensor_A;
    tensor_Cell{2}=tensor_B;
    tensor_Cell{3}=tensor_B;
    tensor_Cell{4}=tensor_A;
    
    % tensor_a=ncon({tensor_A,conj(tensor_A)},...
    %     {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
    % tensor_b=ncon({tensor_B,conj(tensor_B)},...
    %     {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
    
    % initial preparation of a cell of (random) tensors C, T
    [tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
    
    m_list=zeros(1,round((H_FIELD_MAX-H_FIELD_MIN)/H_FIELD_STEP)+1);
    iteration_TEBD_list=[];
    iteration_aCTMRG_avg_list=[];
    
    for H_FIELD=H_FIELD_MIN:H_FIELD_STEP:H_FIELD_MAX
        sprintf('H_FIELD=%.4f',H_FIELD)
        % gate tensor generation
        gate_tensor=gate_TFIM(J_ISING,H_FIELD,DELTA_T,ISING_DIRECTION);
        
        temp_m_old=0;
        temp_m=0;
        
        convergence_repeat_count=0;
        converged_iteration=0;
        
        iteration_TEBD=0;
        iteration_aCTMRG=0;
        iteration_aCTMRG_list=[];
        while iteration_TEBD<ITERATION_MAX_TEBD &&...
                convergence_repeat_count<CONVERGENCE_REPEAT_TEBD
            
            iteration_TEBD=iteration_TEBD+1;
            temp_m_old=temp_m;
            sprintf('repeat=%d/%d',repeat,REPEAT_MAX)
            
            % GATE APPLICATION & UPDATE
            % A-rightward-B gate application & CTMRG & update tensor A, B
            sprintf('H_FIELD=%.4f (%.1f%%), iTEBD iteration=%d, Rightward',...
                H_FIELD,100.0*(H_FIELD-H_FIELD_MIN)/(H_FIELD_MAX-H_FIELD_MIN),iteration_TEBD)
            [tensor_C,tensor_T,iteration_aCTMRG]=...
                CTMRG_anisotropic(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
                tensor_Cell,tensor_C,tensor_T,'sv');
            iteration_aCTMRG_list=[iteration_aCTMRG_list,iteration_aCTMRG];
            tensor_Env=tensor_Env_right(tensor_Cell{1},tensor_Cell{2},tensor_C,tensor_T);
%             [tensor_Cell{1},tensor_Cell{2}]=...
%                 gate_application...
%                 (tensor_Cell{1},tensor_Cell{2},gate_tensor,'rightward');
            [tensor_Cell{1},tensor_Cell{2}]=...
                update_cell_right(tensor_Cell{1},tensor_Cell{2},gate_tensor,...
                tensor_Env,EPSILON_CELL);
            tensor_Cell{3}=tensor_Cell{2};
            tensor_Cell{4}=tensor_Cell{1};
            % A-leftward-B gate application & CTMRG & update tensor A, B
            sprintf('H_FIELD=%.4f (%.1f%%), iTEBD iteration=%d, Leftward',...
                H_FIELD,100.0*(H_FIELD-H_FIELD_MIN)/(H_FIELD_MAX-H_FIELD_MIN),iteration_TEBD)
            [tensor_C,tensor_T,iteration_aCTMRG]=...
                CTMRG_anisotropic(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
                tensor_Cell,tensor_C,tensor_T,'sv');
            iteration_aCTMRG_list=[iteration_aCTMRG_list,iteration_aCTMRG];
            tensor_Env=tensor_Env_left(tensor_Cell{1},tensor_Cell{2},tensor_C,tensor_T);
%             [tensor_Cell{1},tensor_Cell{2}]=...
%                 gate_application...
%                 (tensor_Cell{1},tensor_Cell{2},gate_tensor,'leftward');
            [tensor_Cell{1},tensor_Cell{2}]=...
                update_cell_left(tensor_Cell{1},tensor_Cell{2},gate_tensor,...
                tensor_Env,EPSILON_CELL);
            tensor_Cell{3}=tensor_Cell{2};
            tensor_Cell{4}=tensor_Cell{1};
            % A-upward-B gate application & CTMRG & update tensor A, B
            sprintf('H_FIELD=%.4f (%.1f%%), iTEBD iteration=%d, Upward',...
                H_FIELD,100.0*(H_FIELD-H_FIELD_MIN)/(H_FIELD_MAX-H_FIELD_MIN),iteration_TEBD)
            [tensor_C,tensor_T,iteration_aCTMRG]=...
                CTMRG_anisotropic(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
                tensor_Cell,tensor_C,tensor_T,'sv');
            iteration_aCTMRG_list=[iteration_aCTMRG_list,iteration_aCTMRG];
            tensor_Env=tensor_Env_up(tensor_Cell{1},tensor_Cell{2},tensor_C,tensor_T);
%             [tensor_Cell{1},tensor_Cell{2}]=...
%                 gate_application...
%                 (tensor_Cell{1},tensor_Cell{2},gate_tensor,'upward');
            [tensor_Cell{1},tensor_Cell{2}]=...
                update_cell_up(tensor_Cell{1},tensor_Cell{2},gate_tensor,...
                tensor_Env,EPSILON_CELL);
            tensor_Cell{3}=tensor_Cell{2};
            tensor_Cell{4}=tensor_Cell{1};
            % A-downward-B gate application & CTMRG & update tensor A, B
            sprintf('H_FIELD=%.4f (%.1f%%), iTEBD iteration=%d, Downward',...
                H_FIELD,100.0*(H_FIELD-H_FIELD_MIN)/(H_FIELD_MAX-H_FIELD_MIN),iteration_TEBD)
            [tensor_C,tensor_T,iteration_aCTMRG]=...
                CTMRG_anisotropic(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
                tensor_Cell,tensor_C,tensor_T,'sv');
            iteration_aCTMRG_list=[iteration_aCTMRG_list,iteration_aCTMRG];
            tensor_Env=tensor_Env_down(tensor_Cell{1},tensor_Cell{2},tensor_C,tensor_T);
%             [tensor_Cell{1},tensor_Cell{2}]=...
%                 gate_application...
%                 (tensor_Cell{1},tensor_Cell{2},gate_tensor,'downward');
            [tensor_Cell{1},tensor_Cell{2}]=...
                update_cell_down(tensor_Cell{1},tensor_Cell{2},gate_tensor,...
                tensor_Env,EPSILON_CELL);
            tensor_Cell{3}=tensor_Cell{2};
            tensor_Cell{4}=tensor_Cell{1};
            
            % Calculation of the new environment and the expectation value
            [tensor_C,tensor_T]=...
                CTMRG_anisotropic(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
                tensor_Cell,tensor_C,tensor_T,'sv');
            tensor_Env=tensor_Env_right(tensor_Cell{1},tensor_Cell{2},tensor_C,tensor_T);
            temp_m=exp_value_right(tensor_Cell{1},tensor_Cell{2},tensor_Env,...
                [0.0,1.0;1.0,0.0]);
            
            if ~isfinite(temp_m)
                error('NaN or Inf occurred in the magnetization: CHECK IT')
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
            
        end
        m_list(round((H_FIELD-H_FIELD_MIN)/H_FIELD_STEP)+1)=temp_m;
        iteration_TEBD_list=[iteration_TEBD_list,iteration_TEBD];
        iteration_aCTMRG_avg_list=[iteration_aCTMRG_avg_list,...
            mean(iteration_aCTMRG_list)];
    end
    elapsed_time=toc
    
    %% Plot the average m values
    % Data folder creation
    DATA_FOLDER=data_folder_create('DATA_TEBD_aCTMRG');
    
    figure('visible','off')
    % figure(1)
    plot(H_LIST,m_list,'r-x',...
        H_LIST,iteration_TEBD_list/ITERATION_MAX_TEBD,'b-x',...
        H_LIST,iteration_aCTMRG_avg_list/ITERATION_MAX_ENV,'k-x')
    xlabel('H_Z FIELD')
    ylabel('<X>_{avg}')
    legend('<X>_{avg}','iter TEBD/iter Max','iter aCTMRG/iter Max','Location','northwest')
    title(sprintf(...
        'TFIM 2D quantum [aCTMRG]\n epsm=%.2g, BOND DIM=%d, CHI=%d,\n crepeat CTMRG=%d, crepeat TEBD=%d, max iter CTMRG=%d, max iter TEBD=%d',...
        EPSILON_M,BOND_DIM,CHI,CONVERGENCE_REPEAT_CTMRG,CONVERGENCE_REPEAT_TEBD,ITERATION_MAX_ENV,ITERATION_MAX_TEBD))
    % ylim([0,1.1])
    print('-depsc',strcat('./',DATA_FOLDER,...
        sprintf('/epsm_%.2g_CHI_%d_crepeat_CTMRG_%d_crepeat_TEBD_%d_iter_CTMRG_%d_iter_TEBD_%d.eps',...
        EPSILON_M,CHI,CONVERGENCE_REPEAT_CTMRG,CONVERGENCE_REPEAT_TEBD,ITERATION_MAX_ENV,ITERATION_MAX_TEBD)));
    save(strcat('./',DATA_FOLDER,sprintf('/epsm_%.2g_CHI_%d_crepeat_CTMRG_%d_crepeat_%d_iter_CTMRG_%d_iter_TEBD_%d.mat',...
        EPSILON_M,CHI,CONVERGENCE_REPEAT_CTMRG,CONVERGENCE_REPEAT_TEBD,ITERATION_MAX_ENV,ITERATION_MAX_TEBD)));
    
    % Send 'job-completed' notification email
%     setpref('Internet','E_mail','hojoon.kim@kias.re.kr');
%     setpref('Internet','SMTP_Server','mail.kias.re.kr');
%     setpref('Internet','SMTP_Username','hojoon.kim');
%     setpref('Internet','SMTP_Password','kias7321!');
%     sendmail('hojoon.kim@kias.re.kr','[abacus4] iTEBD_aCTMRG completed',...
%         strcat(sprintf('Data folder for repeat=%d/%d is  \b',repeat,REPEAT_MAX),DATA_FOLDER));
    
    repeat=repeat+1;
end

load chirp
sound(y,Fs)
clear('y','Fs')
