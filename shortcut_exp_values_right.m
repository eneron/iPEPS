[temp_exp_I,temp_exp_I_A,temp_exp_I_B]=...
    exp_value_right(tensor_A,tensor_B,tensor_Env,[1,0;0,1]);
[temp_exp_X,temp_exp_X_A,temp_exp_X_B]=...
    exp_value_right(tensor_A,tensor_B,tensor_Env,[0,1;1,0]);
[temp_exp_Y,temp_exp_Y_A,temp_exp_Y_B]=...
    exp_value_right(tensor_A,tensor_B,tensor_Env,[0,-1j;1j,0]);
[temp_exp_Z,temp_exp_Z_A,temp_exp_Z_B]=...
    exp_value_right(tensor_A,tensor_B,tensor_Env,[1,0;0,-1]);

sprintf('%s =%.2g, %s =%.2g, %s =%.2g, %s =%.2g',...
    '<I>',temp_exp_I,...
    '<X>',temp_exp_X,...
    '<Y>',temp_exp_Y,...
    '<Z>',temp_exp_Z)
sprintf('%s =%.2g, %s =%.2g, %s =%.2g, %s =%.2g',...
    '<I_A>',temp_exp_I_A,...
    '<X_A>',temp_exp_X_A,...
    '<Y_A>',temp_exp_Y_A,...
    '<Z_A>',temp_exp_Z_A)
sprintf('%s =%.2g, %s =%.2g, %s =%.2g, %s =%.2g',...
    '<I_B>',temp_exp_I_B,...
    '<X_B>',temp_exp_X_B,...
    '<Y_B>',temp_exp_Y_B,...
    '<Z_B>',temp_exp_Z_B)