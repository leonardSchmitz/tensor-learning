d =4;
R, a = polynomial_ring(QQ, :a => (1:d,1:d));
C = axis_core_tensor(d);
time_total = 0.0;
time_min = 10000.0;
time_max = 0.0;
nr_rec = 10;
for t in (1:nr_rec)
  if divides(t,1)[1]
    print('*')
  end
  A = generic_transform(d);
  G = matrix_tensor_congruence(A,C);
  F = matrix_tensor_congruence(Array(a),C);
  I = ideal(R,vec(F-G));
  time_cur = @elapsed groebner_basis(I,algorithm=:f4);
  time_total = time_total +  time_cur
  time_min = min(time_min ,  time_cur)
  time_max = max(time_max ,  time_cur)
end
println(time_min)
println(time_total/nr_rec)
println(time_max)


