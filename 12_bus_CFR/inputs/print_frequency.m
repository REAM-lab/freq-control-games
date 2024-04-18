
%% Get path to save results
mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end)-1);
newdir = newdir + "/outputs/";

%% Export results
load("trad2_results.mat")

variables = ["SG1_f", "SG2_f" , "SG3_f", "SG4_f", "SG5_f", "SG6_f"];
file_name = "frequency.csv";
export_results_csv(out, variables, newdir, file_name)