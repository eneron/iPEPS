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

DELTA_T=10^(-3); % time step for iTEBD and Simple update
ITERATION_MAX_SIMPLE=10^4; % maximum number of iteration for 'simple_update'

repeat=1;
REPEAT_MAX=5;
while repeat<=REPEAT_MAX
    % In HCJiang, ITERATION_MAX_SIMPLE=10^5~10^6.
    %% Simple update iteration & Data collection
    sprintf('PHYS_DIM=%d, BOND_DIM=%d, CHI=%d',PHYS_DIM,BOND_DIM,CHI)
    sprintf('ISING_DIRECTION=%s, J_ISING=%.2f, DELTA_T=%.2e',...
        ISING_DIRECTION,J_ISING,DELTA_T)
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
        iteration=0;
        while iteration<ITERATION_MAX_SIMPLE
            [Gamma_A,Gamma_B,lambda]=...
                simple_update(Gamma_A,Gamma_B,lambda,gate_tensor);
            iteration=iteration+1;
        end
        
        tensor_A=ncon({Gamma_A,sqrt(lambda{1}),sqrt(lambda{2}),...
            sqrt(lambda{3}),sqrt(lambda{4})},...
            {[-1,1,2,3,4,],[-2,1],[-3,2],[3,-4],[4,-5]});
        tensor_A=tensor_A/max(abs(tensor_A(:)));

        tensor_B=ncon({Gamma_B,sqrt(lambda{3}),sqrt(lambda{4}),...
            sqrt(lambda{1}),sqrt(lambda{2})},...
            {[-1,1,2,3,4,],[-2,1],[-3,2],[3,-4],[4,-5]});
        tensor_B=tensor_B/max(abs(tensor_B(:)));
        
        % Save workspace variables in a file:
        save(strcat(DATA_FOLDER,'/Simple_Update_',ISING_DIRECTION,...
            sprintf(...
            '_J_%.1f_H_%.4f_delT_%.1e_d_%d_D_%d_CHI_%d_iter_%d.mat',...
            J_ISING,H_FIELD,DELTA_T,PHYS_DIM,BOND_DIM,CHI,ITERATION_MAX_SIMPLE)))
    end
    
    elapsedTime=toc
    save(strcat(DATA_FOLDER,'/elapsedTime.txt'),'elapsedTime','-ascii')
    
    % Send 'job-completed' notification email
    setpref('Internet','E_mail','hojoon.kim@kias.re.kr');
    setpref('Internet','SMTP_Server','mail.kias.re.kr');
    setpref('Internet','SMTP_Username','hojoon.kim');
    setpref('Internet','SMTP_Password','kias7321!');
    sendmail('hojoon.kim@kias.re.kr','[abacus3] Simple_Update (iteration) completed',...
        strcat(sprintf('ITERATION_MAX_SIMPLE=%d, Data for repeat=%d/%d is  \b',...
        ITERATION_MAX_SIMPLE,repeat,REPEAT_MAX),DATA_FOLDER));
    
    repeat=repeat+1;
end
