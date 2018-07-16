%%
m_ = 0.265;
g_ = - 9.81;
q_0_ = 0.3014;
A_ = 0.0183*0.0183*pi;
r_b_ = 0.026;

%%

for i=1:size(q_, 2)

    position_in_matrix_tot = (i + 2*(i - 1));
    [A_total(position_in_matrix_tot:position_in_matrix_tot + 2, :), b_total(position_in_matrix_tot:position_in_matrix_tot + 2, 1)] = SolveSystemParamID(m_, g_, q_0_, A_, r_b_, q_(1,i), q_(2,i), q_(3,i), pressure_static(1, i), pressure_static(2, i), pressure_static(3, i), theta_draw_wire_(i), phi_muslce_(i));

end

%%
sol = A_total\b_total;
