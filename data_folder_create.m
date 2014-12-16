function [ DATA_FOLDER ] = data_folder_create( folder_name )
% Creat a DATA_FOLDER (in the current folder)
% named as the input string plus current date and a natural number;
% Return the name as the output string.
temp_i=1;
while exist...
        (strcat(folder_name,'_',datestr(now,'yyyymmdd'),...
        '_',num2str(temp_i)),'dir')
    temp_i=temp_i+1;
end
mkdir(strcat(folder_name,'_',datestr(now,'yyyymmdd'),...
    '_',num2str(temp_i)));
DATA_FOLDER=strcat(folder_name,'_',datestr(now,'yyyymmdd'),...
    '_',num2str(temp_i));
end

