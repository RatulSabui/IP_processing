function [first_scan,last_scan,file_list_run,sav_nam] = file_nam(directory,run_num)
    
    dir_list = dir(directory);
    file_list = string(zeros(1,length(dir_list)));
    %file_tif = [];

    for i = 1:length(dir_list)
        name = dir_list(i).name;
        file_list(i) = name;
    end
    % selecting the names with tiff
    TF = contains(file_list,'.tif');
    file_list = file_list(TF);
    %TF_sc1 = contains(file_list,'_sc1');
    %file_list_sc1 = file_list(TF_sc1);
    TF_run = contains(file_list,run_num);
    file_list_run = file_list(TF_run);
    for i = 1:length(file_list_run)
        file_list_run(i) = strcat(directory,file_list_run(i));
    end
    last_scan = file_list_run(end);
    first_scan = file_list_run(1);
    sav_nam = erase(first_scan,["_sc1",".tif"]);
    sav_nam = strcat(sav_nam,".csv");
end
