function  export_results_csv(struct_to_open, variables, target_folder, file_name)

nvars = size(variables',1);

for k=1:nvars
    if k==1
        tiss = struct_to_open.(variables(k)) ;
    else
        tiss = [tiss, struct_to_open.(variables(k))];
    end
end

time_vector = tiss.Time ; % Get vector of time to get rid of 'sec' string
TT = timeseries2timetable( tiss );
TT = renamevars(TT,TT.Properties.VariableNames, variables);

TT = timetable2table(TT); % table is easier to handle

TT.Time = time_vector ; % replce time column with time as float

writetable(TT, target_folder + file_name, 'Delimiter',',')

fprintf("Ok - File printed %s \n", target_folder + file_name)

end
