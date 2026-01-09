function pstruct = prc_sutton_k1_binary_namep(pvec)
% made by tom

pstruct = struct;

pstruct.mu      = pvec(1);
pstruct.Rhat    = pvec(2);
pstruct.vhat_1  = pvec(3);
pstruct.h_1     = pvec(4);


end