d = 2;
# case A_22 = 0
R, a, g = polynomial_ring(QQ, :a => (1:d,1:d), :g => (1:d,1:d,1:d));
C = axis_core_tensor(2);
h = matrix_tensor_congruence(Array(a),C);
I = ideal(R,vec(g-h));
@assert ideal_membership(g[1,2,2]-g[2,1,2],I+ideal(R,[a[2,2]]))
f = det(matrix(R,a))*a[2,2] - QQ(1,3)*(g[1,2,2] - g[2,1,2])
@assert ideal_membership(f,I)
@assert ideal_membership(g[2,2,2]-a[2,1]^3,I+ideal(R,[a[2,2]]))
@assert ideal_membership(g[2,1,2]*g[1,2,2]-g[2,2,2]*g[1,1,2],I+ideal(R,[a[2,2]]))
@assert ideal_membership(g[2,2,1]*g[2,1,2]-g[1,2,1]*g[2,2,2],I+ideal(R,[a[2,2]]))
# show H_111 \not= 0
J = ideal(R,[g[1,2,1]-g[2,1,1],
             g[1,2,1]*g[2,1,1]-g[1,1,1]*g[2,2,1], 
             g[1,1,2]*g[1,2,1]-g[2,1,2]*g[1,1,1], 
             g[1,1,1]]);
@assert radical_membership(a[1,2]*g[2,2,2],J+I)
@assert radical_membership(a[1,1]*g[2,2,2],J+I)
#case A_22 \not= 0 
R,a = rational_function_field(QQ, :a => (1:d,1:d));
C = axis_core_tensor(d);
s = matrix_tensor_congruence(Array(a),C);
## determinant, so case study in a[2,2] remains. Assume a[2,2]!=0:
x = (s[2,1,1]-s[1,2,1])*inv(s[1,2,2]-s[2,1,2]);
@assert x == - a[1,2]*inv(a[2,2])
H = matrix_tensor_congruence(Array(R[1 x; 0 1]),s);
@assert H[1,2,1] == H[2,1,1]
@assert H[1,2,1]*H[2,1,1] == H[1,1,1]*H[2,2,1]
@assert H[1,1,2]*H[1,2,1] == H[2,1,2]*H[1,1,1]
