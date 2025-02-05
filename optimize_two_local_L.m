function [form,gradPhi,gradCond_vert,gradCond_horz,HessPhi,HessCond,form1,form2] = optimize_two_local_L(Phi_opt,conduct_vertopt,conduct_horzopt,boundaryPhi,boundaryI,gradL_condhorz,gradL_condvert)
%OPTIMIZE_TWO 此处显示有关此函数的摘要
%   此处显示详细说明
Num = size(boundaryPhi,2);
smoothpara = 5000*0;
gap=0.4;
boundpara =100*0;
boundayrate=10;
innerrate=1;


n=size(conduct_vertopt,2);
% J_vert=(sign(conduct_vertopt-0.1)-1)/2;
% J_horz=(sign(conduct_horzopt-0.1)-1)/2;
% K_vert=(sign(conduct_vertopt-1.2)+1)/2;
% K_horz=(sign(conduct_horzopt-1.2)+1)/2;


% diff_vert = conduct_vertopt(1:end-1,:)-conduct_vertopt(2:end,:);
% diff_horz = conduct_horzopt(:,1:end-1)-conduct_horzopt(:,2:end);
% diff_vert_horz = conduct_vertopt(:,1:end-1)-conduct_vertopt(:,2:end);
% diff_horz_vert = conduct_horzopt(1:end-1,:)-conduct_horzopt(2:end,:);


% J_diffvert = (sign(diff_vert+gap)-1)/2;
% J_diffhorz = (sign(diff_horz+gap)-1)/2;
% K_diffvert = (sign(diff_vert-gap)+1)/2;
% K_diffhorz = (sign(diff_horz-gap)+1)/2;
% smoothform = norm(J_diffvert.*(diff_vert+gap),2)^2+norm(J_diffhorz.*(diff_horz+gap),2)^2+norm(K_diffvert.*(diff_vert-gap),2)^2+norm(K_diffhorz.*(diff_horz-gap),2)^2;
% 
% J_diffvert_horz=(sign(diff_vert_horz+gap)-1)/2;
% J_diffhorz_vert=(sign(diff_horz_vert+gap)-1)/2;
% K_diffvert_horz=(sign(diff_vert_horz-gap)+1)/2;
% K_diffhorz_vert=(sign(diff_horz_vert-gap)+1)/2;
% smoothform = smoothform+norm(J_diffvert_horz.*(diff_vert_horz+gap),2)^2+norm(J_diffhorz_vert.*(diff_horz_vert+gap),2)^2+norm(K_diffvert_horz.*(diff_vert_horz-gap),2)^2+norm(K_diffhorz_vert.*(diff_horz_vert-gap),2)^2;


Phi_opt1=vertcat(Phi_opt*boundaryPhi,boundaryPhi);
Phi_opt2=vertcat(Phi_opt*boundaryPhi,zeros(4*n,Num));


conduct_horizonal= conduct_horzopt; %%水平电导分布
conduct_vertical= conduct_vertopt; %%垂直电导分布


 L = generateL(conduct_horizonal,conduct_vertical,boundayrate,innerrate);



beforeform = L*Phi_opt1-innerrate*Phi_opt2.^3+vertcat(zeros(n^2,Num),boundaryI*boundayrate);
form1=norm(reshape(beforeform,[],1),2)^2;
% form2=smoothpara*smoothform/2;
form2 = 0;
form = form1 + form2;
gradPhi = 2*L'*beforeform*boundaryPhi'-6*(innerrate*Phi_opt2.^2.*beforeform)*boundaryPhi'; gradPhi=gradPhi(1:n^2,1:4*n);
% L_cut = L(1:end,1:n^2);
% if nargout >4
% A = L_cut'*L_cut;B=boundaryPhi*boundaryPhi';
% A=reshape(A,[n^2,1,n^2,1]);
% B= reshape(B,[1,4*n,1,4*n]);
% HessPhi = A.*B;
% HessPhi = reshape(HessPhi,4*n*n*n,4*n*n*n);
% end
% gradPhi = gradPhi / norm(gradPhi,2);


% afterform = 2*pagemtimes(gradL_condvert,Phi_opt1);


% beforeform1 = repmat(beforeform,1,1,n+1,n);afterform = beforeform1;
gradCond_vert = zeros(n+1,n);

% gradCond_vert1 = zeros(n+1,n);

for i = 1:n+1
    for j = 1:n
    G = reshape(gradL_condvert(:,(j-1)*(n+1)+i),n^2+4*n,n^2+4*n);
    [row,col,v]=find(G);
    
%     B=zeros(n^2+4*n,4*n);
    for nu = 1:size(row)
%          B(row(nu),:) = B(row(nu),:)+v(nu)*Phi_opt1(col(nu),:);
        gradCond_vert(i,j) = gradCond_vert(i,j)+2*v(nu)*dot(Phi_opt1(col(nu),:),beforeform(row(nu),:));
    end
%     afterform(:,:,i,j)=2*B;
    
%      gradCond_vert1(i,j)=sum(2*B.*beforeform,'all');
    end
end
% beforeform1 = repmat(beforeform,1,1,n+1,n);
% dotprod = afterform.*beforeform1;
% gradCond_vert = sum(dotprod,[1,2]);
% gradCond_vert = reshape(gradCond_vert,n+1,n);
% gradCond_vert =gradCond_vert -boundpara*(J_vert).*conduct_vertopt+boundpara*(K_vert).*(conduct_vertopt-1.2)...
%     -vertcat(smoothpara*(J_diffvert).*(diff_vert+gap),zeros(1,n))+vertcat(zeros(1,n),smoothpara*(J_diffvert).*(diff_vert+gap))...
%     +vertcat(smoothpara*(K_diffvert).*(diff_vert-gap),zeros(1,n))-vertcat(zeros(1,n),smoothpara*(K_diffvert).*(diff_vert-gap))...
%     -horzcat(smoothpara*(J_diffvert_horz).*(diff_vert_horz+gap),zeros(n+1,1))+horzcat(zeros(n+1,1),smoothpara*(J_diffvert_horz).*(diff_vert_horz+gap))...
%     +horzcat(smoothpara*(K_diffvert_horz).*(diff_vert_horz-gap),zeros(n+1,1))-horzcat(zeros(n+1,1),smoothpara*(K_diffvert_horz).*(diff_vert_horz-gap));
% Jane = reshape(afterform,(n^2+4*n)*Num,n^2+n);

% gradCond_vert = gradCond_vert/norm(gradCond_vert,2);

% afterform = 2*pagemtimes(gradL_condhorz,Phi_opt1);
 
% beforeform1 = repmat(beforeform,1,1,n,n+1);afterform = beforeform1;
gradCond_horz = zeros(n,n+1);
for i = 1:n
    for j = 1:n+1
    G = reshape(gradL_condhorz(:,(j-1)*(n)+i),n^2+4*n,n^2+4*n);
    [row,col,v]=find(G);
    

    for nu = 1:size(row)
        gradCond_horz(i,j) = gradCond_horz(i,j)+2*v(nu)*dot(Phi_opt1(col(nu),:),beforeform(row(nu),:));
    end
    
    end
end

% dotprod = afterform.*beforeform1;
% gradCond_horz = sum(dotprod,[1,2]);
% gradCond_horz = reshape(gradCond_horz,n,n+1);

% gradCond_horz =gradCond_horz-boundpara*(J_horz).*conduct_horzopt+boundpara*(K_horz).*(conduct_horzopt-1.2)...
%     -horzcat(smoothpara*(J_diffhorz).*(diff_horz+gap),zeros(n,1))+horzcat(zeros(n,1),smoothpara*(J_diffhorz).*(diff_horz+gap))...
%     +horzcat(smoothpara*(K_diffhorz).*(diff_horz-gap),zeros(n,1))-horzcat(zeros(n,1),smoothpara*(K_diffhorz).*(diff_horz-gap))...
%     -vertcat(smoothpara*(J_diffhorz_vert).*(diff_horz_vert+gap),zeros(1,n+1))+vertcat(zeros(1,n+1),smoothpara*(J_diffhorz_vert).*(diff_horz_vert+gap))...
%     +vertcat(smoothpara*(K_diffhorz_vert).*(diff_horz_vert-gap),zeros(1,n+1))-vertcat(zeros(1,n+1),smoothpara*(K_diffhorz_vert).*(diff_horz_vert-gap));

% if nargout > 5
%     Janey=reshape(afterform,(n^2+4*n)*Num,n^2+n);
%     
%     Jack = cat(2,Jane,Janey);
%     
%     J_vert = reshape(J_vert,n^2+n,1); 
%     J_horz = reshape(J_horz,n^2+n,1); 
%     K_vert = reshape(K_vert,n^2+n,1); 
%     K_horz = reshape(K_horz,n^2+n,1); 
%     
%     Jd_vert = reshape(vertcat(J_diffvert,zeros(1,n))+vertcat(zeros(1,n),J_diffvert),n^2+n,1);
%     Jd_horz = reshape(horzcat(zeros(n,1),J_diffhorz)+horzcat(J_diffhorz,zeros(n,1)),n^2+n,1);
%     Kd_vert = reshape(vertcat(K_diffvert,zeros(1,n))+vertcat(zeros(1,n),K_diffvert),n^2+n,1);
%     Kd_horz = reshape(horzcat(zeros(n,1),K_diffhorz)+horzcat(K_diffhorz,zeros(n,1)),n^2+n,1);
%     Jd_vert_horz = reshape(horzcat(zeros(n+1,1),J_diffvert_horz)+horzcat(J_diffvert_horz,zeros(n+1,1)),n^2+n,1);
%     Jd_horz_vert = reshape(vertcat(J_diffhorz_vert,zeros(1,n+1))+vertcat(zeros(1,n+1),J_diffhorz_vert),n^2+n,1);
%     Kd_vert_horz = reshape(horzcat(zeros(n+1,1),K_diffvert_horz)+horzcat(K_diffvert_horz,zeros(n+1,1)),n^2+n,1);
%     Kd_horz_vert = reshape(vertcat(K_diffhorz_vert,zeros(1,n+1))+vertcat(zeros(1,n+1),K_diffhorz_vert),n^2+n,1);
%     
%     Jind_vert = find(vertcat(J_diffvert,zeros(1,n))~=0);Jend_vert = Jind_vert +1;
%     Jind_horz = find(horzcat(J_diffhorz,zeros(n,1))~=0)+n*(n+1);Jend_horz = Jind_horz +n;
%     Kind_vert = find(vertcat(K_diffvert,zeros(1,n))~=0);Kend_vert = Kind_vert +1;
%     Kind_horz = find(horzcat(K_diffhorz,zeros(n,1))~=0)+n*(n+1);Kend_horz = Kind_horz +n;
%     Jind_vert_horz = find(horzcat(J_diffvert_horz,zeros(n+1,1))~=0);Jend_vert_horz = Jind_vert_horz +n;
%     Jind_horz_vert = find(vertcat(J_diffhorz_vert,zeros(1,n+1))~=0)+n*(n+1);Jend_horz_vert = Jind_horz_vert +1;
%     Kind_vert_horz = find(horzcat(K_diffvert_horz,zeros(n+1,1))~=0);Kend_vert_horz = Kind_vert_horz +n;
%     Kind_horz_vert = find(vertcat(K_diffhorz_vert,zeros(1,n+1))~=0)+n*(n+1);Kend_horz_vert = Kind_horz_vert +1;
%     
%     Jd_interact = zeros(2*n*(n+1));
%     Jd_interact(Jind_vert,Jend_vert)=2*eye(size(Jind_vert,1));
%     Jd_interact(Jend_vert,Jind_vert)=2*eye(size(Jind_vert,1));
%     Jd_interact(Jind_horz,Jend_horz)=2*eye(size(Jind_horz,1));
%     Jd_interact(Jend_horz,Jind_horz)=2*eye(size(Jind_horz,1));
%     Kd_interact = zeros(2*n*(n+1));
%     Kd_interact(Kind_vert,Kend_vert)=2*eye(size(Kind_vert,1));
%     Kd_interact(Kend_vert,Kind_vert)=2*eye(size(Kind_vert,1));
%     Kd_interact(Kind_horz,Kend_horz)=2*eye(size(Kind_horz,1));
%     Kd_interact(Kend_horz,Kind_horz)=2*eye(size(Kind_horz,1));
%     
%     Jdc_interact = zeros(2*n*(n+1));
%     Jdc_interact(Jind_vert_horz,Jend_vert_horz)=2*eye(size(Jind_vert_horz,1));
%     Jdc_interact(Jend_vert_horz,Jind_vert_horz)=2*eye(size(Jind_vert_horz,1));
%     Jdc_interact(Jind_horz_vert,Jend_horz_vert)=2*eye(size(Jend_horz_vert,1));
%     Jdc_interact(Jend_horz_vert,Jind_horz_vert)=2*eye(size(Jend_horz_vert,1));
%     Kdc_interact = zeros(2*n*(n+1));
%     Kdc_interact(Kind_vert_horz,Kend_vert_horz)=2*eye(size(Kind_vert_horz,1));
%     Kdc_interact(Kend_vert_horz,Kind_vert_horz)=2*eye(size(Kind_vert_horz,1));
%     Kdc_interact(Kind_horz_vert,Kend_horz_vert)=2*eye(size(Kend_horz_vert,1));
%     Kdc_interact(Kend_horz_vert,Kind_horz_vert)=2*eye(size(Kend_horz_vert,1));
%     
%     J= vertcat(J_vert,J_horz);K= vertcat(K_vert,K_horz);
%     Jd = vertcat(Jd_vert,Jd_horz);Kd = vertcat(Kd_vert,Kd_horz);
%     Jdc = vertcat(Jd_vert_horz,Jd_horz_vert);Kdc = vertcat(Kd_vert_horz,Kd_horz_vert);
%     HessCond = Jack'*Jack/2-boundpara*diag(J)+boundpara*diag(K)...
%         -smoothpara*diag(Jd)-smoothpara*Jd_interact...
%         +smoothpara*diag(Kd)-smoothpara*Kd_interact...
%         -smoothpara*diag(Jdc)-smoothpara*Jdc_interact...
%         +smoothpara*diag(Kdc)-smoothpara*Kdc_interact;
% end



end

 

