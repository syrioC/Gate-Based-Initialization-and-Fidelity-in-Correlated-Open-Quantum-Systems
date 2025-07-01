function [GB,dGB,TB]=timedDkT(T,bcf,C)
% 
% C=reshape(C,1,numel(C));
% nf=numel(C);
% N=numel(T)-1;
% dt=T(2)-T(1);
% bcf=bcf(:);T=T(:);
% bcend=bcf(end);
% kr=512;
% krb=kr/2;
% kc=N/kr;
% tp=T(kr+1);
% tr=dt*(0:kr)';
% bcf=reshape(bcf(1:end-1),kr,kc);
% bcf(kr+1,1:kc-1)=bcf(1,2:kc);
% bcf(kr+1,kc)=bcend;
% 
% % GB=zeros(krb,kc,nf);
% % dGB=zeros(krb,kc,nf);
% % %find SD
% % parfor ic=1:kc
% %     z=bcf(:,ic).*exp(1i*(tp*(ic-1)+tr).*C);
% %     z=(z(1:2:end-2,:)+4*z(2:2:end-1,:)+z(3:2:end,:))*dt/3;    
% %     res=0;
% %     zb=zeros(krb,nf);
% %     for j=1:krb
% %         res=res+z(j,:);
% %         zb(j,:)=res;
% %     end
% %     GB(:,ic,:)=zb;
% % end
% % 
% % for  ic=2:kc
% %     GB(:,ic,:)=GB(end,ic-1,:)+GB(:,ic,:);
% % end
% % 
% % % find derivative 
% % parfor ic=1:kc
% %     yb=zeros(krb,nf)
% %     gb=reshape(GB(:,ic,:),krb,nf);
% %     yb(1,:)=gb(1,:);
% %     for j=2:krb
% %         yb(j,:)=gb(j,:)+yb(j-1,:);
% %     end
% %     dGB(:,ic,:)=reshape(yb,krb,1,nf);
% % end
% % 
% % for  ic=2:kc
% %     dGB(:,ic,:)=dGB(end,ic-1,:)+dGB(:,ic,:);
% % end
% % 
% %         
% % 
% % N1=N/2;
% % GB=reshape(GB,krb*kc,nf);
% % dGB=reshape(dGB,krb*kc,nf);
% % 
% % 
% % dGB=1i*2*dt*((1:N1)'.*GB-dGB);
% % 
% % GB=[zeros(1,nf); GB];
% % dGB=[zeros(1,nf); dGB];
% % TB=T(1:2:end);
% % end

%%%ALTERNATIVE%%%

C=reshape(C,1,numel(C));
nf=numel(C);
N=numel(T)-1;
dt=T(2)-T(1);
bcf=bcf(:);T=T(:);
bcend=bcf(end);
kr=512;
krb=kr/2;
kc=N/kr;
tp=T(kr+1);
tr=dt*(0:kr)';
bcf=reshape(bcf(1:end-1),kr,kc);
bcf(kr+1,1:kc-1)=bcf(1,2:kc);
bcf(kr+1,kc)=bcend;

GB=zeros(krb,kc,nf);
dGB=zeros(krb,kc,nf);
parfor ic=1:kc
    z=bcf(:,ic).*exp(1i*(tp*(ic-1)+tr).*C);
    y=1i*(tp*(ic-1)+tr).*z; %derivative with omega
    z=(z(1:2:end-2,:)+4*z(2:2:end-1,:)+z(3:2:end,:))*dt/3;
    y=(y(1:2:end-2,:)+4*y(2:2:end-1,:)+y(3:2:end,:))*dt/3;
    
    res=0;
    re=0;
    zb=zeros(krb,nf);
    yb=zeros(krb,nf)
    for j=1:krb
        res=res+z(j,:);
        zb(j,:)=res;
        re=re+y(j,:);
        yb(j,:)=re;
    end
    GB(:,ic,:)=zb;
    dGB(:,ic,:)=yb;
end

for  ic=2:kc
    GB(:,ic,:)=GB(end,ic-1,:)+GB(:,ic,:);
    dGB(:,ic,:)=dGB(end,ic-1,:)+dGB(:,ic,:);
end

GB=reshape(GB,krb*kc,nf);
GB=[zeros(1,nf); GB];
dGB=reshape(dGB,krb*kc,nf);
dGB=[zeros(1,nf); dGB];
TB=T(1:2:end);
end
