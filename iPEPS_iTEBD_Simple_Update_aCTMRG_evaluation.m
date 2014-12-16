%% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!

% PHYS_DIM=2;
% BOND_DIM=2;
% CHI=20;
%
% J_ISING=-1.0;
% ISING_DIRECTION='X';
% H_FIELD_MIN=-3.2; % H=J_ISING*X_i*X_j+H_FIELD/4*(Z_i+Z_j)
% H_FIELD_STEP=-0.01; % H_FIELD_C=-3.044
% H_FIELD_MAX=-2.8;
%
% DELTA_T=10^(-1); % time step for iTEBD and Simple update
% EPSILON_BULK=10^(-7); % convergence limit for the simple update
% ITERATION_MAX_SIMPLE=10^3; % maximum number of iteration for 'simple_update'
% CONVERGENCE_REPEAT_SIMPLE_UPDATE=3;


% In HCJiang, ITERATION_MAX_SIMPLE=10^5~10^6.
%% Evaluate magnetic fields from save data
m_x_A_list=[];
m_x_B_list=[];
m_x_list=[];
% m_y_A_list=[];
% m_y_B_list=[];
% m_y_list=[];
% m_z_A_list=[];
% m_z_B_list=[];
% m_z_list=[];

rng('shuffle'); % to get the "shuffled" random numbers!
% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
EPSILON_ENV=10^(-3); % convergence limit for CTMRG
CONVERGENCE_REPEAT_CTMRG=3;
ITERATION_MAX_ENV=500; % maximum number of iteration for 'CTMRG_directional'

DATA_FOLDER='DATA_Simple_Update_20141205_1';

tic;

for H_FIELD=H_FIELD_MAX:H_FIELD_STEP:H_FIELD_MIN
    sprintf('H_FIELD=%.4f',H_FIELD)
    temp_file=dir(strcat(DATA_FOLDER,sprintf(...
        '/*H_%.4f*.mat',H_FIELD)));
    load(strcat(DATA_FOLDER,'/',temp_file.name),...
        'tensor_A','tensor_B')
    
    tensor_Cell{1}=tensor_A;
    tensor_Cell{2}=tensor_B;
    tensor_Cell{3}=tensor_B;
    tensor_Cell{4}=tensor_A;
    
    [tensor_C,tensor_T]=...
        CTMRG_anisotropic(EPSILON_ENV,CONVERGENCE_REPEAT_CTMRG,ITERATION_MAX_ENV,...
        tensor_Cell,tensor_C,tensor_T,'sv');
    
    tensor_Env=tensor_Env_right(tensor_A,tensor_B,tensor_C,tensor_T);
    temp_m=exp_value_right(tensor_A,tensor_B,tensor_Env,[0.0,1.0;1.0,0.0]);
    m_x_list=[m_x_list,temp_m];
    
    %         % Calculate expectation values
    %         [temp_m,temp_m_A,temp_m_B]=...
    %             exp_value_right(tensor_A,tensor_B,...
    %             tensor_Env_right(tensor_A,tensor_B,tensor_C,tensor_T),...
    %             [1.0,0.0;0.0,-1.0]);
    %         m_z_A_list=[m_z_A_list,temp_m_A];
    %         m_z_B_list=[m_z_B_list,temp_m_B];
    %         m_z_list=[m_z_list,temp_m];
end

elapsed_time=toc

x=H_FIELD_MAX:H_FIELD_STEP:H_FIELD_MIN;
figure(1)
plot(x,real(m_x_list),'r-x')
xlabel('H FIELD')
legend('<X>_{avg}','Location','northwest')

load chirp
sound(y,Fs)
clear('y','Fs')

% figure(2);
% plot(x,m_x_list,'b',x,m_x_A_list,'b--',x,m_x_B_list,'b:')
% xlabel('H_FIELD')
% legend('<X>','<X_{A}>','<X_{B}>','Location','northwest')