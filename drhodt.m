function [drho]=drhodt(i0,Slambda,rho_t,A_0)
    drho=reshape(Slambda(i0,:,:),2,2)*rho_t*A_0-A_0*reshape(Slambda(i0,:,:),2,2)*rho_t ...
        + A_0*rho_t*reshape(Slambda(i0,:,:),2,2)'-rho_t*reshape(Slambda(i0,:,:),2,2)'*A_0;
end

