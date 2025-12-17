d = 3;
R,a = rational_function_field(QQ, :a => (1:d,1:d));
C = axis_core_tensor(d);
s = matrix_tensor_congruence(Array(a),C);
lambda23 = (s[d,d-1,d-1]- s[d-1,d,d-1])*inv(-s[d, d-1, d]+s[d-1, d, d]);
H1 = matrix_tensor_congruence(Array(R[1 0 0; 0 1 lambda23; 0 0 1]),s);
@assert H1[2,3,2] == H1[3,2,2]
lambda13 = (H1[d,d-1,d-2]- H1[d-1,d,d-2])*inv(-H1[d, d-1, d]+H1[d-1, d, d]);
H2 = matrix_tensor_congruence(Array(R[1 0 lambda13; 0 1 0; 0 0 1]),H1);
lambda12 = (H2[d-1,d-2,d-2]- H2[d-2,d-1,d-2])*inv(-H2[d-1, d-2,d-1]+H2[d-2, d-1, d-1]);
H3 = matrix_tensor_congruence(Array(R[1 lambda12 0; 0 1 0; 0 0 1]),H2);
@assert all([H3[1,i,1] == H3[i,1,1] for i in (1:3)])
@assert all([H3[i,1,1]*H3[1,j,1] == H3[i,j,1]*H3[1,1,1] for i in (1:3) for j in (1:3)])
