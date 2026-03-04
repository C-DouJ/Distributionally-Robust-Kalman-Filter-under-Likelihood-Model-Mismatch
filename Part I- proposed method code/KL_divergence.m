function kl = KL_divergence(u0, S0, u1, S1)
    u0 = u0(:);
    u1 = u1(:);

    S1_inv_S0 = S1\ S0;  % S1^{-1} * S0
    trace_term = trace(S1_inv_S0);
    diff_u = u1 - u0;
    quad_term = diff_u' / S1 * diff_u;
    log_det_S1 = log(det(S1));
    log_det_S0 = log(det(S0));  
    log_det_term = log_det_S1 - log_det_S0;
    kl = 0.5 * (trace_term - length(u0) + quad_term + log_det_term);
end