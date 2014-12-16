%% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!

PHYS_DIM=2;
EPSILON_ENV=10^(-3); % convergence limit for CTMRG
ITERATION_MAX_ENV=300; % maximum number of iteration for 'full_update_CTM'
CONVERGENCE_REPEAT=5;

% Skip the input-error checks by the tensor contraction module 'ncon':
global ncon_skipCheckInputs;
ncon_skipCheckInputs = true;

% Turn off the warning 'Complex values of exp_values bulabula~'
warning('off','all');

% h=waitbar(0,'butts up!');
% steps=40;

DATA_FOLDER=data_folder_create('DATA_dCTMRG_convergence_check');

tic;

for BOND_DIM=2:4
    for CHI=BOND_DIM:(BOND_DIM^2+5)
        %         waitbar(((BOND_DIM-2)*10+(CHI-BOND_DIM^2+1))/steps)
        %% INITIAL STATE & CTM GENERATION (RANDOM; Normalized A,B, CTM)
        % initial preparation of (random) tensors A, B, a, b
        [tensor_A,tensor_B,~,~]=...
            unit_cell_generation(PHYS_DIM,BOND_DIM);
        tensor_a=ncon({tensor_A,conj(tensor_A)},...
            {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
        tensor_b=ncon({tensor_B,conj(tensor_B)},...
            {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
        % initial preparation of a cell of (random) tensors C, T
        [tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
        %% CTMRG PROCESS CONJUGATE (METHOD 1: Follow Orus.2009)
        sprintf('BOND_DIM=%d, CHI=%d',BOND_DIM,CHI)
        epsilon_list=[];
        
        m_z_A_list=[];
        m_z_B_list=[];
        m_z_list=[];
        m_x_A_list=[];
        m_x_B_list=[];
        m_x_list=[];
        m_y_A_list=[];
        m_y_B_list=[];
        m_y_list=[];
        
        tic;
        % Determine convergence by a series of epsilon values rather than by two
        % values
        convergence_repeat_count=0;
        converged_iteration=0;
        
        epsilon_env_old=10^5;
        epsilon_env=10^5;
        iteration=0;
        
        while convergence_repeat_count<CONVERGENCE_REPEAT &&...
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
            
            % Calculate expectation values
            tensor_Env=tensor_Env_right(tensor_A,tensor_B,tensor_C,tensor_T);
            
            [temp_m,temp_m_A,temp_m_B]=...
                exp_value_right(tensor_A,tensor_B,tensor_Env,[1.,0.;0.,-1.]);
            m_z_A_list=[m_z_A_list,temp_m_A];
            m_z_B_list=[m_z_B_list,temp_m_B];
            m_z_list=[m_z_list,temp_m];
            
            [temp_m,temp_m_A,temp_m_B]=...
                exp_value_right(tensor_A,tensor_B,tensor_Env,[0.,1.;1.,0.]);
            m_x_A_list=[m_x_A_list,temp_m_A];
            m_x_B_list=[m_x_B_list,temp_m_B];
            m_x_list=[m_x_list,temp_m];
            
            [temp_m,temp_m_A,temp_m_B]=...
                exp_value_right(tensor_A,tensor_B,tensor_Env,[0.,-1.j;1.j,0.]);
            m_y_A_list=[m_y_A_list,temp_m_A];
            m_y_B_list=[m_y_B_list,temp_m_B];
            m_y_list=[m_y_list,temp_m];
            
            epsilon_env=...
                convergence_CTM(tensor_C_old,tensor_T_old,tensor_C,tensor_T,...
                tensor_A,tensor_B);
            epsilon_list=[epsilon_list,epsilon_env];
            
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
        
        elapsedTime=toc;
        
        if iteration>5
            m_max=max(abs([m_x_list(end-5:end),m_y_list(end-5:end),...
                m_z_list(end-5:end)]))+0.05;
        else
            m_max=max(abs([m_x_list,m_y_list,...
                m_z_list]));
        end
        
        figure('visible','off')
        plot(epsilon_list)
        ylim([0,0.1])
        xlabel('Iteration')
        legend('epsilon\_Env')
        title(sprintf('D=%d, CHI=%d: epsilon=[%.5g, %.5g]',...
            BOND_DIM,CHI,min(epsilon_list),max(epsilon_list)))
        print('-depsc',...
            strcat('./',DATA_FOLDER,sprintf('/D_%d_CHI_%d_eps.eps',...
            BOND_DIM,CHI)));
        
        x=1:iteration;
        
        plot(x,m_x_list,'b',x,m_y_list,'g',x,m_z_list,'r')
        ylim([-m_max,m_max])
        xlabel('Iteration')
        legend('<X>_{avg}','<Y>_{avg}','<Z>_{avg}','Location','northwest')
        title(sprintf('D=%d, CHI=%d',...
            BOND_DIM,CHI))
        print('-depsc',...
            strcat('./',DATA_FOLDER,sprintf('/D_%d_CHI_%d_m_value.eps',...
            BOND_DIM,CHI)));
        
        plot(x,epsilon_list,'k',x,m_x_list,'b',x,m_y_list,'g',x,m_z_list,'r')
        ylim([-m_max,m_max])
        xlabel('Iteration')
        legend('eps','<X>_{avg}','<Y>_{avg}','<Z>_{avg}','Location','northwest')
        title(sprintf('D=%d, CHI=%d',...
            BOND_DIM,CHI))
        print('-depsc',...
            strcat('./',DATA_FOLDER,sprintf('/D_%d_CHI_%d_m_eps.eps',...
            BOND_DIM,CHI)));
        
        save(strcat('./',DATA_FOLDER,sprintf('/D_%d_CHI_%d_m_eps',...
            BOND_DIM,CHI)))
    end
    % Send 'job-completed' notification email
    setpref('Internet','E_mail','hojoon.kim@kias.re.kr');
    setpref('Internet','SMTP_Server','mail.kias.re.kr');
    setpref('Internet','SMTP_Username','hojoon.kim');
    setpref('Internet','SMTP_Password','kias7321!');
    sendmail('hojoon.kim@kias.re.kr','[abacus3] dCTMRG check',...
        strcat(sprintf('Data folder for BOND_DIM=%d is  \b',BOND_DIM),DATA_FOLDER));
end

% close(h)
warning('on','all');
