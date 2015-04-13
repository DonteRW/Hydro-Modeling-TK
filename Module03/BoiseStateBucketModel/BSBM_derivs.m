function [dVdt,qdro,qb,qtro] = BSBM_derivs(P,E,T,phi,V,dt)

    % Retrieve parameters from input vector phi
    beta  = phi(1);
    V1max = phi(2);
    k12   = phi(3);
    V2max = phi(4);
    k23   = phi(5);
    V3max = phi(6);
    kb    = phi(7);
    
    % Retrieve states from input vector V
    V1 = V(1);
    V2 = V(2);
    V3 = V(3);
    
    % Compute potential change in storage deficit over dt
    qdro  = P*(V1/V1max)^beta;
    qf    = P - qdro;
    SD1dt = (V1max - V1)/dt - (E + k12*V1) + qf;
    SD2dt = (V2max - V2)/dt - (T + k23*V2) + k12*V1;
    SD3dt = (V3max - V3)/dt - kb*V3 + k23*V2;

    % Compute percolation rates
    qf  = min(P-qdro,SD1dt);
    q12 = min(k12*V1,SD2dt);
    q23 = min(k23*V2,SD3dt);
    qb  = kb*V3;
    
    % Compute derivatives
    dV1dt = qf - (E + q12);
    dV2dt = q12 - (T + q23);
    dV3dt = q23 - qb;
    
    % Compute total runoff
    qtro = qb + qdro;
    
    % Store computed derivatives in dVdt
    dVdt(1,1) = dV1dt;
    dVdt(2,1) = dV2dt;
    dVdt(3,1) = dV3dt;
    
end