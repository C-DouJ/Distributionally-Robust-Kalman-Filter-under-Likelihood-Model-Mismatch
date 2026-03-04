function [xkk,Pkk]=DRKF(xkk,Pkk,F,H,z,Q,R)

iter = 20;
SearchIter = 5;
%%%%%Time update

xk1k=F*xkk;

Pk1k=F*Pkk*F'+Q;

%%%%%Measurement update

xkk = xk1k;
Pkk = Pk1k;
Pk1k_inv = Pk1k \ eye(size(Pk1k));

for nn = 1:iter
    Pkk_inv = Pkk \ eye(size(Pkk));
    S = H*Pkk*H'+R;
    S_inv = S \ eye(size(S));
    gradx = Pk1k_inv*(xkk-xk1k) - H'*S_inv*(z-H*xkk);

    gradP = 0.5*(-Pkk_inv + Pk1k_inv + H'*S_inv*H - H'*S_inv*(z-H*xkk)*(z-H*xkk)'*S_inv*H);
    ghat = 2*Pkk*gradP*Pkk;
    t = 1/nn;
    L_cost_now = L_cost(xkk,Pkk,Pk1k, xk1k, H, R, z);
    for i = 1:SearchIter
        xkk_new = xkk - t*Pkk*gradx;
        Pkk_new = Pkk - t*ghat + 0.5*t^2*ghat*Pkk_inv*ghat;
        Pkk_new = (Pkk_new+Pkk_new')/2;
        L_cost_new = L_cost(xkk_new,Pkk_new,Pk1k, xk1k, H, R, z);
        if L_cost_new>L_cost_now - 1e-3 * t * norm(gradx)^2
            t = 0.2*t;
        else
            break
        end
    end

    if KL_divergence(xkk_new,Pkk_new,xkk,Pkk) < 1e-3
        xkk = xkk_new ;Pkk = Pkk_new;
        break;
    end
    xkk = xkk_new ;Pkk = Pkk_new;

end

