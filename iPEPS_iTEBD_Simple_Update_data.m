% 20141103 (MON) KIAS
% 20141127 (Thu) KIAS
% 20141128 (Fri) KIAS
%% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!

PHYS_DIM=2;
BOND_DIM=2;
CHI=20;

J_ISING=-1.0;
ISING_DIRECTION='X';
H_FIELD_MIN=-3.2; % H=J_ISING*X_i*X_j+H_FIELD/4*(Z_i+Z_j)
H_FIELD_STEP=-0.01; % H_FIELD_C=-3.044
H_FIELD_MAX=-2.8;

DELTA_T=10^(-2); % time step for iTEBD and Simple update
EPSILON_SIMPLE=10^(-8); % convergence limit for the simple update
ITERATION_MAX_SIMPLE=10^4; % maximum number of iteration for 'simple_update'
CONVERGENCE_REPEAT_SIMPLE_UPDATE=3;
CONVERGENCE_REPEAT_CTMRG=3;

repeat=1;
REPEAT_MAX=5;
while repeat<=REPEAT_MAX
    % In HCJiang, ITERATION_MAX_SIMPLE=10^5~10^6.
    %% Simple update iteration & Data collection
    sprintf('PHYS_DIM=%d, BOND_DIM=%d, CHI=%d',PHYS_DIM,BOND_DIM,CHI)
    sprintf('ISING_DIRECTION=%s, J_ISING=%.2f, DELTA_T=%.2e, EPSILON_SIMPLE=%.2e',...
        ISING_DIRECTION,J_ISING,DELTA_T,EPSILON_SIMPLE)
    sprintf('H_FIELD=%.2f:%.2f:%.2f',H_FIELD_MIN,H_FIELD_STEP,H_FIELD_MAX)
    sprintf('ITERATION_MAX_SIMPLE=%.2e',ITERATION_MAX_SIMPLE)
    
    % % parameter check query
    % if strcmp(input('Right parameters? Shoot? (yes/no) ','s'),'no')
    %     error('Adjust the parameters.')
    % end
    
    % elapsed time check
    tic
    % Data folder creation
    DATA_FOLDER=data_folder_create('DATA_Simple_Update');
    iteration_list=[];
    
    for H_FIELD=H_FIELD_MAX:H_FIELD_STEP:H_FIELD_MIN
        %         multiWaitbar('H_FIELD iterating...',...
        %             (H_FIELD_MAX-H_FIELD)/(H_FIELD_MAX-H_FIELD_MIN-H_FIELD_STEP),...
        %             'Color','g');
        
        sprintf('repeat=%d/%d, H_FIELD= %.4f (%.1f%%)',repeat,REPEAT_MAX,H_FIELD,...
            100*(H_FIELD_MAX-H_FIELD)/(H_FIELD_MAX-H_FIELD_MIN-H_FIELD_STEP))
        
        % gate tensor generation
        gate_tensor=gate_TFIM(J_ISING,H_FIELD,DELTA_T,ISING_DIRECTION);
        
        % initial preparation of (random) tensors A, B, a, b
        [Gamma_A,Gamma_B,lambda]=...
            unit_cell_generation_simple_update(PHYS_DIM,BOND_DIM);
        
        % SIMPLE UPDATE
        epsilon_simple=10^5;
        epsilon_simple_old=0;
        convergence_repeat_count=0;
        converged_iteration=0;
        iteration=0;
        while iteration<ITERATION_MAX_SIMPLE &&...
                convergence_repeat_count<=CONVERGENCE_REPEAT_SIMPLE_UPDATE
                
            lambda_old=lambda;
            epsilon_simple_old=epsilon_simple;
            
            [Gamma_A,Gamma_B,lambda]=...
                simple_update(Gamma_A,Gamma_B,lambda,gate_tensor);
            
            epsilon_simple=sum(abs(diag(lambda{1})-diag(lambda_old{1})))+...
                sum(abs(diag(lambda{2})-diag(lambda_old{2})))+...
                sum(abs(diag(lambda{3})-diag(lambda_old{3})))+...
                sum(abs(diag(lambda{4})-diag(lambda_old{4})));
            
            % convergence test: convergence in a series of m values.
            if abs(epsilon_simple-epsilon_simple_old)<EPSILON_SIMPLE
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
            
            iteration=iteration+1;
            %         multiWaitbar('Simple updating...',iteration/ITERATION_MAX_SIMPLE,...
            %             'Color','b');
        end
        
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
        
        iteration_list=[iteration_list,iteration];
        
        % Save workspace variables in a file:
        save(strcat(DATA_FOLDER,'/Simple_Update_',ISING_DIRECTION,...
            sprintf(...
            '_J_%.1f_H_%.4f_delT_%.1e_eps_%.1e_d_%d_D_%d_CHI_%d_iter_%d.mat',...
            J_ISING,H_FIELD,DELTA_T,EPSILON_SIMPLE,PHYS_DIM,BOND_DIM,CHI,ITERATION_MAX_SIMPLE)))
    end
    % multiWaitbar('Simple updating...','close');
    % multiWaitbar('H_FIELD iterating...','close');
    
    elapsedTime=toc
    save(strcat(DATA_FOLDER,'/elapsedTime.txt'),'elapsedTime','-ascii')
    
    % Send 'job-completed' notification email
    setpref('Internet','E_mail','hojoon.kim@kias.re.kr');
    setpref('Internet','SMTP_Server','mail.kias.re.kr');
    setpref('Internet','SMTP_Username','hojoon.kim');
    setpref('Internet','SMTP_Password','kias7321!');
    sendmail('hojoon.kim@kias.re.kr','[abacus4] Simple_Update completed',...
        strcat(sprintf('ITERATION_MAX_SIMPLE=%d, Data for repeat=%d/%d is  \b',...
        ITERATION_MAX_SIMPLE,repeat,REPEAT_MAX),DATA_FOLDER));
    
    repeat=repeat+1;
end
