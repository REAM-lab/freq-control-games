% We use this plot to calculate the settling time. Uncomment the load
% command to load the results for the first time.

%load("nash_results.mat")

t = out.SG1_f.Time;
plot(t, out.SG1_f.Data, t , out.SG2_f.Data, t, out.SG3_f.Data, t, out.SG4_f.Data, t, out.SG5_f.Data, t, out.SG6_f.Data)
grid
xticks(12:0.5:32)
legend
