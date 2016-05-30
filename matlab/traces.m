
ic = [-1.00	1.00	0.00	0.00	0.00	0.00	-1.00	1.00	
1.00	-1.00	-1.00	1.00	0.00	0.00	0.00	0.00	
0.00	0.00	1.00	-1.00	-1.00	1.00	0.00	0.00	
0.00	0.00	0.00	0.00	1.00	-1.00	1.00	-1.00	];

rs = [ -1.3187	 -0.3093	  6.4103	
  4.8223	  3.9712	  3.8978	
 -7.1369	 -5.0653	 -1.6420	
 -3.0547	 -2.6071	  5.1663	];

rk = [ -5.0493	 -3.4499	  3.2456	
-14.8055	 -5.1366	 -0.8397	
 -0.2254	  3.2372	 -3.7783	
  3.6110	  6.4743	 -7.2285	];

model = struct('ic', ic, 'rs', rs, 'rk', rk);
y0 = initial_state(model, -120);


dt = 0.01;
vms = -35:15:20;
t_end = 2.0;

G = zeros(length(0:dt:t_end), length(vms));
for i = 1:length(vms)
    Q = transition_matrix(model, vms(i));
    ex = expm(Q * dt);
    y = y0; G(1, i) = y(1);
    for j = 1:length(dt:dt:t_end)
        y = ex * y; G(j, i+1) = y(1);
    end
end

vm = -20;
Q = transition_matrix(model, vms(i));
ex = expm(Q * dt);
Y = zeros(length(0:dt:t_end), 4);
Y(1, :) = y0;
for i = 1:length(dt:dt:t_end)
    Y(i+1, :) = ex * Y(i, :)';
end

plot(0:dt:t_end, Y)