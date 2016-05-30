function [ y0 ] = initial_state( model, vm )

vars = [1; vm/100; tanh((vm+20)/50.0)];
y0 = exp(model.rs * vars);
y0 = y0 / sum(y0);


end

