d =5;
C = axis_core_tensor(d);
time_total = 0.0;
time_min = 10000.0;
time_max = 0.0;
nr_rec = 10000;
for t in (1:nr_rec)
  if divides(t,100)[1]
    print('*')
  end
  A = generic_transform(d);
  G = matrix_tensor_congruence(A,C);
  time_cur = @elapsed tensor_learning(G);
  Q = tensor_learning(G);
  time_total = time_total +  time_cur
  time_min = min(time_min ,  time_cur)
  time_max = max(time_max ,  time_cur)
  @assert matrix_tensor_congruence(Q,G) == C
  @assert inv(Q) == A
end
println(time_min)
println(time_total/nr_rec)
println(time_max)
