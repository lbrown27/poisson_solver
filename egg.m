clear all;
close all;
clc;

dat_num = load('sol_num.dat');

dat_ref = load('sol_ref.dat');

h = pcolor(dat_num);
set(h, 'EdgeColor', 'none');
title("Numerical Solution");
colorbar
caxis([-.5 .5])
figure();

h = pcolor(dat_ref);
set(h, 'EdgeColor', 'none');
colorbar
caxis([-.5 .5])
title("Analytical (Exact) Solution");

figure();

h= pcolor(dat_num - dat_ref);
set(h, 'EdgeColor', 'none');
colorbar
%caxis([-.5 .5])
title("Error");