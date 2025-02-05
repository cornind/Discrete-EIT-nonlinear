n=2;   %%方形网格边长
AA = repmat(0.4*(1:n)',1,n+1);
BB = repmat(0.2*(1:n+1),n,1);
CC = (sin(AA+BB)).^2;
conduct_horizonal= (CC+0.5)./1.5; %%水平电导分布


AA = repmat(2*(1:n)'/11,1,n+1);
BB = repmat(4*(1:n+1)/11,n,1);
CC = (sin(AA+BB)).^2;
conduct_vertical= (CC'+0.5)./1.5;  %%垂直电导分布
% conduct_horizonal= ones(n,n+1); %%水平电导分布
% conduct_vertical= ones(n+1,n);  %%垂直电导分布
bb = repmat(1:4*n,4*n,1);aa = bb';
phi = eye(4*n,4*n);



I= zeros(4*n,4*n);U=zeros(n^2+4*n,4*n);
[I,U]=experi_func(phi,n,conduct_horizonal,conduct_vertical);


% parfor id=1:4*n
%     [II(:,id),UU(:,id)]=solvebyCurrent(I(:,id),n,conduct_horizonal,conduct_vertical);
% end
errI = I.*randn(4*n,4*n)*0.01/100*0;
I = I+errI;
boundayrate = 1;innerrate = 1;
% [gradL_condhorz,gradL_condvert] = generateGrad(n,boundayrate);
%     save("gradL_cond.mat","gradL_condvert","gradL_condhorz");

save("diriclet_Neum3.mat","I","n","conduct_vertical","conduct_horizonal","U","phi","errI");
save("boundayrate.mat","boundayrate","innerrate")

% 
% 
% n=5;   %%方形网格边长
% AA = repmat(2*(1:n)'/n,1,n+1);
% BB = repmat(2*(1:n+1)/(n+1),n,1);
% CC = (sin(AA+BB)).^2;
% % conduct_horizonal= (CC+0.5)./1.5; %%水平电导分布
% % conduct_vertical= (CC'+0.5)./1.5;  %%垂直电导分布
% conduct_horizonal= ones(n,n+1); %%水平电导分布
% conduct_vertical= ones(n+1,n);  %%垂直电导分布
% phi = eye(4*n);
% I= zeros(4*n,4*n);U=zeros(n^2+4*n,4*n);
% parfor id=1:4*n
%     [I(:,id),U(:,id)]=experi_func(phi(:,id),n,conduct_horizonal,conduct_vertical);
% end
% % I = I+I.*randn(4*n,4*n)/10000;
% boundayrate = n;
% I_vain = I;U_vain = U;
% save("diriclet_Neum_vain.mat","I_vain","U_vain");

