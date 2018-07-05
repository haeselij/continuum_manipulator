syms C_ohne_l_bar k_spline k;
%%
D_damp_spline_ = 0;
m_ = 0.265;
g_ = - 9.81;
q_0_ = 0.3014;
A_ = 0.0183*0.0183*pi;
D_damp_ = 2000;
r_b_ = 0.026;


%%
A_total = zeros(3*size(q_,1), 3);
b_total = zeros(3*size(q_,1), 1);

for i=1:size(q_,1)
    position_in_matrix_tot = (i + 3*(i-1));
    [A_total(position_in_matrix_tot:position_in_matrix_tot+2,:), b_total(position_in_matrix_tot:position_in_matrix_tot+2,:)] = SolveSystemParamID(D_damp_spline_, m_, g_, q_0_, A_, D_damp_, r_b_, q_(i,1), q_(i,2), q_(i,3), 0, 0, 0, 0, 0, 0, 0, 0, 0, pressure_static(i), 0, 0, theta_draw_wire_(i), pi);
end
%%
sol = A_total\b_total;