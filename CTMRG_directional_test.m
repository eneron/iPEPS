%% PARAMETERS
rng('shuffle'); % to get the "shuffled" random numbers!

PHYS_DIM=2;
BOND_DIM=3;
CHI=4;

EPSILON_ENV=10^(-2); % convergence limit for CTMRG
ITERATION_MAX_ENV=300; % maximum number of iteration for 'full_update_CTM'
CONVERGENCE_REPEAT=3;
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
%% INITIAL STATE |00> (SINGULAR) & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(1,1,1,1,1)=1;
tensor_B(1,1,1,1,1)=1;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |00> & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(1,:,:,:,:)=1;
tensor_B(1,:,:,:,:)=1;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |11> & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(2,:,:,:,:)=1;
tensor_B(2,:,:,:,:)=1;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

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
%% INITIAL STATE |++>_X (SINGULAR) & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(1,1,1,1,1)=1;
tensor_B(1,1,1,1,1)=1;
tensor_A(2,1,1,1,1)=1;
tensor_B(2,1,1,1,1)=1;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |->_X|->_X & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(2,:,:,:,:)=-1;
tensor_B(2,:,:,:,:)=-1;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |+>_Y|+>_Y & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(2,:,:,:,:)=1j;
tensor_B(2,:,:,:,:)=1j;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |->_Y|->_Y & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(2,:,:,:,:)=-1j;
tensor_B(2,:,:,:,:)=-1j;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |0+> & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=zeros(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(1,:,:,:,:)=1;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |+>_X|+>_Y & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_B(2,:,:,:,:)=1j;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% INITIAL STATE |->_X|+>_Y & RANDOM CTM GENERATION
% initial preparation of trivial tensors A, B, a, b
tensor_A=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
tensor_B=ones(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);

tensor_A(2,:,:,:,:)=-1;
tensor_B(2,:,:,:,:)=1j;

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

% initial preparation of a cell of (random) tensors C, T
[tensor_C,tensor_T]=init_CTM_generation(BOND_DIM,CHI);
%% CTMRG PROCESS CONJUGATE (METHOD 1: Follow Orus.2009)
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
    multiWaitbar('dCTMRG iterating...',iteration/ITERATION_MAX_ENV,...
        'Color','g');
    
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
multiWaitbar('dCTMRG iterating...','close');
if iteration>=ITERATION_MAX_ENV
    disp('Maximum iteration reached.')
end

elapsedTime=toc

x=1:iteration;

if iteration>10
    m_max=max(abs([m_x_list(end-10:end),m_y_list(end-10:end),...
        m_z_list(end-10:end)]))+0.05;
else
    m_max=max(abs([m_x_list(end-3:end),m_y_list(end-3:end),...
        m_z_list(end-3:end)]))+0.05;
end

figure(1)
plot(x,epsilon_list)
xlabel('Iteration')
ylabel('epsilon Env')
legend('epsilon Env','Location','northwest')
title(sprintf('D=%d, CHI=%d: epsilon=[%.5g, %.5g]',...
    BOND_DIM,CHI,min(epsilon_list),max(epsilon_list)))
ylim([0,max(epsilon_list(end-(iteration-5):end))+0.01])

figure(2)
plot(x,m_x_list,'b',x,m_y_list,'g',x,m_z_list,'r')
xlabel('Iteration')
legend('<X>_{avg}','<Y>_{avg}','<Z>_{avg}','Location','northwest')
ylim([-m_max,m_max])

figure(3)
plot(x,epsilon_list,'k',...
    x,m_x_list,'b',x,m_y_list,'g',x,m_z_list,'r')
xlabel('Iteration')
ylim([-m_max,m_max])
legend('epsilon Env','<X>_{avg}','<Y>_{avg}','<Z>_{avg}',...
    'Location','northwest')
% figure(3);
% plot(x,m_z_list,'r',x,m_z_A_list,'r--',x,m_z_B_list,'r:')
% xlabel('Iteration')
% legend('m_z','m_z_A','m_z_B','Location','northwest')

% shortcut_rho_eig

sprintf('%s=%.5g','epsilon_min',min(epsilon_list))
sprintf('%s=%.5g','epsilon_max',max(epsilon_list))

load chirp
sound(y,Fs)
clear('y','Fs')