function [ Q ] = transition_matrix( model , vm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ic = model.ic;
rs = model.rs;
rk = model.rk;

args = [1; vm/100; tanh((vm+20)/50.0)];
N = size(rs, 1);
E = size(rk, 1);

D = zeros(E, 2*E);
D(:, 1:2:end) = eye(E);
D(:, 2:2:end) = -eye(E);
D = [D;abs(D)];

r_vec = exp(D^-1 * [ic(:, 1:2:end)'*rs; rk] * args);
Q = ic * diag(r_vec(:)) * (-min(ic, 0)');

end

