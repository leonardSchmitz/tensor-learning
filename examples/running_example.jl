d = 4;
C = axis_core_tensor(d);
A = Array(QQ[0 -1 0 -1 ; -1 1 0 1 ; 1 0 0 1 ; -1 0 1 1])
G = matrix_tensor_congruence(A,C);
Q = tensor_learning(G);
@assert matrix_tensor_congruence(Q,G) == C
@assert inv(Q) == A
# algorithm in every step 
# s = 1
U1 = Array(QQ[1 1 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1]);
H1 = matrix_tensor_congruence(U1,G);
L1 = Array(QQ[1 0 0 0; -1 1 0 0 ; 1 0 1 0 ; -1 0 0 1]);
J1 = matrix_tensor_congruence(L1,H1);
D1 = Array(QQ[-1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1]);
G1 = matrix_tensor_congruence(D1,J1);
# s = 2
Q1 = Array(QQ[1 0 0 0; 0 1 0 0 ; 0 0 1 (G1[4,3,3] - G1[3,4,3])*inv(G1[3,4,4] - G1[4,3,4]) ; 0 0 0 1]);
G1p = matrix_tensor_congruence(Q1,G1);
Q2 = Array(QQ[1 0 0 0; 0 1 0 (G1p[4,3,2] - G1p[3,4,2])*inv(G1p[3,4,4] - G1p[4,3,4]) ; 0 0 1 0 ; 0 0 0 1]);
G1pp = matrix_tensor_congruence(Q2,G1p);
Q3 = Array(QQ[1 0 0 0; 0 1 (G1pp[3,2,2] - G1pp[2,3,2])*inv(G1pp[2,3,3] - G1pp[3,2,3]) 0 ; 0 0 1 0 ; 0 0 0 1]);
G2 = matrix_tensor_congruence(Q3,G1pp);
# s = 3
L3 = Array(QQ[1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 1 1]);
J3 = matrix_tensor_congruence(L3,G2);
D3 = Array(QQ[1 0 0 0; 0 1 0 0 ; 0 0 -1 0 ; 0 0 0 1]);
G3 = matrix_tensor_congruence(D3,J3);
@assert G3 == C
@assert inv(D3*L3*Q3*Q2*Q1*D1*L1*U1) == A





