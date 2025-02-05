
format long
newData1 = load('-mat', "diriclet_Neum3");

% 在基础工作区中从这些字段创建新变量。
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
clear i vars newData1 

load boundayrate.mat

boundaryPhi = phi(:,3*n+2:end);

boundaryI = I(:,3*n+2:end);

Num = size(phi,2);
Firsttime = true;
load boundayrate.mat
    
if Firsttime
    x=ones(n^2*4*n +2*n^2+2*n,1);
   
    save("X.mat","x");
    %关于电导梯度生成
%     load diriclet_Neum_vain.mat
% 
%     Phi_vain = U_vain(1:n^2,:); 
    load boundayrate.mat

    [gradL_condhorz,gradL_condvert ] = generateGrad(n,boundayrate);
    save("gradL_cond.mat","gradL_condvert","gradL_condhorz");



end
load X.mat; 
load gradL_cond.mat



% boundayrate = 10*n;save("boundayrate.mat","boundayrate")


    Phi_optres = (reshape(TransformFunc1(x(1:n^2*4*n )),n^2,4*n ));
    conduct_vertoptres = (reshape(TransformFunc(x(n^2*4*n +1:n^2*4*n +1*n^2+n)),n+1,n));
    conduct_horzoptres =  (reshape(TransformFunc(x(n^2*4*n +n^2+n+1:n^2*4*n +2*n^2+2*n)),n,n+1));
% fun = @(x)getObjValue1(x,conduct_vertoptres,conduct_horzoptres,boundaryPhi,boundaryI,gradL_condhorz,gradL_condvert);
% [fval]= fun(x(1:n^2*4*n ))
y=reshape(conduct_vertical,[],1);
y=vertcat(y,reshape(conduct_horizonal,[],1));
y = TransformFuncInv(y);
y=vertcat(reshape(U(1:n^2,:),[],1),y);
         m = zeros(size(x));
    v = zeros(size(x));

     fun = @(x)getObjValueAll(x,conduct_vertoptres,conduct_horzoptres,boundaryPhi,boundaryI,gradL_condhorz,gradL_condvert );
[val]=fun(y)
Num=1000;
Record = zeros(0,0);
MAX=1e30;
OST=0;
for id = 1:100000
% xx = TransformFunc(x);
% conduct_vertoptres = reshape(abs(xx(1:1*n^2+n)),n+1,n);
% conduct_horzoptres = reshape(abs(xx(n^2+n+1:2*n^2+2*n)),n,n+1);
% testforadam
% d_N_eig
% boundaryPhi=phi;boundaryI=I*phi;
     fun = @(x)getObjValueAll(x,conduct_vertoptres,conduct_horzoptres,boundaryPhi,boundaryI,gradL_condhorz,gradL_condvert );

% [val,grad] = fun(x);
% for i = 1:size(x,1)
%     xx = x; xx(i) = xx(i) + 1e-8;
%     vali = fun(xx);
%     gg(i)= (vali-val)*1e8;
%     gg = gg';
% 
% end
% grad = grad(1:end);
% figure
%     plot(1:size(gg,1),gg,1:size(grad,1),grad);

%      [val,x,record,m,v] = adam_optimizer(0.002, 1-0.1,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
     [val,x,record,m,v] = adam_optimizer(0.001, 1-0.01,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
     % [val,x,record,m,v] = adam_optimizer(0.0005, 1-0.01,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
%       [val,x,record,m,v] = adam_optimizer(0.00025, 1-0.001,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
%        [val,x,record,m,v] = adam_optimizer(0.00012, 1-0.001,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
%        [val,x,record,m,v] = adam_optimizer(0.0001, 1-0.001,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
% [val,x,record,m,v] = adam_optimizer(0.00005, 1-0.001,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
%  [val,x,record,m,v] = adam_optimizer(0.000001, 1-0.001,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
%  [val,x,record,m,v] = adam_optimizer(0.0000001, 1-0.001,  1-0.001, 1e-30,  fun,  x, 1000,m,v);
%   [val,x,record,m,v] = adam_optimizer(0.00000005, 1-0.001,  1-0.001, 1e-30,  fun,  x, 1000,m,v);

conduct_vertoptres = (reshape(TransformFunc(x(n^2*4*n +1:n^2*4*n +1*n^2+n)),n+1,n));
        conduct_horzoptres =  (reshape(TransformFunc(x(n^2*4*n +n^2+n+1:n^2*4*n +2*n^2+2*n)),n,n+1));
         Record=vertcat(Record,record);
        MIN = min(record);
        
%         TEMP = MAX;
%         
%         MAX = max(record);
%         checker =sum(abs(record(2:end)-record(1:end-1)))/(MAX-MIN);
       [MIN]
%         if checker>10
%             OST=OST+1
%             m=zeros(size(m));
%             v=zeros(size(v));
%             
%         end
        

end

 


        
Phi_optres = (reshape(TransformFunc1(x(1:n^2*4*n )),n^2,4*n ));
    conduct_vertoptres = (reshape(TransformFunc(x(n^2*4*n +1:n^2*4*n +1*n^2+n)),n+1,n));
    conduct_horzoptres =  (reshape(TransformFunc(x(n^2*4*n +n^2+n+1:n^2*4*n +2*n^2+2*n)),n,n+1));
figure
surfc(conduct_horzoptres);

Phi_optres_slice = reshape(Phi_optres(:,1),n,n);
err_vertical = conduct_vertical-conduct_vertoptres;
err_horizonal = conduct_horizonal-conduct_horzoptres;
figure
surfc((err_horizonal));
figure
surfc(conduct_horizonal);  
save("X.mat","x");

temp = Record;
load Record.mat
Record = vertcat(Record,temp);
save("Record.mat","Record")
createfigure(Record)

% 目标函数 
function [objValue,grad,Hess] = getObjValueAll(parameter,conduct_vertoptres,conduct_horzoptres,boundaryPhi,boundaryI,gradL_condhorz,gradL_condvert )
    n = size(conduct_vertoptres,2);
    NUM_data = size(boundaryI,2);
    % 需要优化的参数
    parameter1 = parameter(1:n^2*4*n );
    parameter2 = parameter(n^2*4*n +1:n^2*4*n +2*n*n+2*n);
    [parameter2,gradout,Hessout] = TransformFunc(parameter2);
    [parameter1,gradout1] = TransformFunc1(parameter1);


    ph =  reshape(parameter1,n^2,4*n );
    conduct_vertopt =  (reshape(parameter2(1:n^2+n),n+1,n));
    conduct_horzopt =  (reshape(parameter2(n^2+n+1:2*n^2+2*n),n,n+1));
    
    % 读取训练数据
    
    % 准确率（非显式过程，需要调用子函数）
   if nargout >2
    [acc,grad1,grad_vert,grad_horz,Hess1,Hess] = optimize_two_local_L(ph,conduct_vertopt,conduct_horzopt,boundaryPhi,boundaryI,gradL_condhorz,gradL_condvert ); Hess1=2*Hess1;
    grad = zeros(2*n^2+2*n,1);
    grad(1:n^2+n)=reshape(grad_vert,n^2+n,1);
    grad(n^2+n+1:2*n^2+2*n)=reshape(grad_horz,n^2+n,1);
    Hess = Hess.*(gradout*gradout')+diag(grad).*Hessout;
    else
        [acc,grad1,grad_vert,grad_horz] = optimize_two_local_L(ph,conduct_vertopt,conduct_horzopt,boundaryPhi,boundaryI,gradL_condhorz,gradL_condvert );
    end

    grad = zeros(2*n^2+2*n,1);
    grad(1:n^2+n)=reshape(grad_vert,n^2+n,1);
    grad(n^2+n+1:2*n^2+2*n)=reshape(grad_horz,n^2+n,1);
    grad =  grad.*gradout;
    % 目标函数值
    grad1= reshape(grad1,n^2*4*n ,1).*gradout1;

    grad = vertcat(grad1,grad);
    objValue = acc;
   
end

function [out,gradout] = TransformFunc1(in)
    
        out = (in);
    gradout = ones(size(in));
%     out = (atan(in)+pi/2)/2;
%     gradout = 1./(in.*in+1)./2;
%     Hessout = diag(-2*in.*(1./(in.*in+1)).^2./2);
end

function [out] = TransformFuncInv(in)
    
        out =  tan(2*(in - 0.1)-pi/2);
end

function [out,gradout,Hessout] = TransformFunc(in)
%     out = abs(in);
%     gradout = sign(in);
%     Hessout = 0;
    out = (atan(in)+pi/2)/2+0.1;
    gradout = 1./(in.*in+1)./2;
    Hessout = diag(-2*in.*(1./(in.*in+1)).^2./2);
end

function [val,theta,record,m,v] = adam_optimizer(alpha, beta1, beta2, epsilon, f, theta,Times,m,v)
    % Initialize variables
%     m = zeros(size(theta));
%     v = zeros(size(theta));


 if nargin > 2
       record = zeros(Times,1); 
    end
    t = 0;
    [val]=f(theta);
    % Loop until convergence
    while t<Times% define your convergence criterion
        % Increment timestep
        t = t + 1;

        % Compute gradients
        [val,g]=f(theta);val;

        % Update biased first moment estimate
        m = beta1 * m + (1 - beta1) * g;

        % Update biased second raw moment estimate
        v  = beta2 * v  + (1 - beta2) * (g  .^ 2);

        % Compute bias-corrected first moment estimate
        m_hat  = m  / (1 - beta1^t*0);

        % Compute bias-corrected second raw moment estimate
        v_hat  = v  / (1 - beta2^t*0);

        % Update parameters
        theta = theta  - alpha * m_hat  ./ (sqrt(v_hat) + epsilon);

        % Update variables for next iteration
%         if mod(t,60)==0
%             alpha=alpha/1.5;
%         end

         if nargin > 2
            record(t) = val; 
        end
    end
            [val ]=f(theta);

end
