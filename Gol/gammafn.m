% This function is implementing the Gamma Function used to integrate Golomb
% model

function [V] = gammafn(vv,theta,sigma)

V = 1/(1 + exp(-(vv-theta)/sigma));
