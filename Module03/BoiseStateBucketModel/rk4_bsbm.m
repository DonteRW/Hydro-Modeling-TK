function [Vnplus1,qdro,qb,qtro,mberr] = rk4_bsbm(Vn,phi,U,dt)

    % Retrieve forcings from input vector
    P = U(1);
    E = U(2);
    T = U(3);
    
    % Retrieve parameters needed here
    beta  = phi(1);
    V1max = phi(2);
    V2max = phi(4);
    V3max = phi(6);
    kb    = phi(7);
    
    [k1,qdro,qb,qtro] = BSBM_derivs(P,E,T,phi,Vn,0.0); %#ok<*NASGU,*ASGLU>
    [k2,qdro,qb,qtro] = BSBM_derivs(P,E,T,phi,(Vn + (dt/2)*k1),0.5*dt);
    [k3,qdro,qb,qtro] = BSBM_derivs(P,E,T,phi,(Vn + (dt/2)*k2),0.5*dt);
    [k4,qdro,qb,qtro] = BSBM_derivs(P,E,T,phi,(Vn + dt*k3),dt);

    Vnplus1 = Vn + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % Check for overshoots as mass balance errors
    mberr = 0.0;
    
    % Upper reservoir
    if(Vnplus1(1)<0.0)
        mberr = mberr + abs(Vnplus1(1));
        Vnplus(1) = 0.0;
    elseif(Vnplus1(1)>V1max)
        mberr = mberr + Vnplus1(1) - V1max;
        Vnplus(1) = V1max;
    end
    
    % Middle reservoir
    if(Vnplus1(2)<0.0)
        mberr = mberr + abs(Vnplus1(2));
        Vnplus(2) = 0.0;
    elseif(Vnplus1(2)>V2max)
        mberr = mberr + Vnplus1(2) - V2max;
        Vnplus(2) = V2max;
    end

    % Middle reservoir
    if(Vnplus1(3)<0.0)
        mberr = mberr + abs(Vnplus1(3));
        Vnplus(3) = 0.0;
    elseif(Vnplus1(3)>V3max)
        mberr = mberr + Vnplus1(3) - V3max;
        Vnplus(3) = V3max;
    end
    
    % Diagnose runoff from final states
    V1 = Vnplus1(1);
    V3 = Vnplus1(3);
    qdro = P*(V1/V1max)^beta;
    qb   = kb*V3;
    qtro = qb + qdro;
    
return