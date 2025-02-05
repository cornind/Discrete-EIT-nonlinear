function [I,U] = experi_func(phi,n,conduct_horizonal,conduct_vertical,tol)
%EXPERI_FUNC 此处显示有关此函数的摘要
%   此处显示详细说明

    L=zeros(n^2+4*n,n^2+4*n);
    L(n^2+1:end,n^2+1:end)=eye(4*n);

    L(1,1)=-conduct_vertical(1,1)-conduct_vertical(2,1)-conduct_horizonal(1,1)-conduct_horizonal(1,2);
    L(1,2)=conduct_horizonal(1,2); L(1,n+1)= conduct_vertical(2,1);L(1,n^2+3*n+1)=conduct_vertical(1,1);L(1,n^2+3*n)=conduct_horizonal(1,1);
    L(n,n)=-conduct_vertical(1,n)-conduct_vertical(2,n)-conduct_horizonal(1,n)-conduct_horizonal(1,n+1);
    L(n,n^2+1)= conduct_horizonal(1,n+1);L(n,2*n)=conduct_vertical(2,n);L(n,n^2+4*n)=conduct_vertical(1,n);L(n,n-1)=conduct_horizonal(1,n);

    L(n^2-n+1,n^2-n+1)=-conduct_vertical(n,1)-conduct_vertical(n+1,1)-conduct_horizonal(n,1)-conduct_horizonal(n,2);
    L(n^2-n+1,n^2-2*n+1)=conduct_vertical(n,1);L(n^2-n+1,n^2+2*n)=conduct_vertical(n+1,1);L(n^2-n+1,n^2+2*n+1)=conduct_horizonal(n,1);L(n^2-n+1,n^2-n+2)=conduct_horizonal(n,2);

    L(n^2,n^2)=-conduct_vertical(n,n)-conduct_vertical(n+1,n)-conduct_horizonal(n,n)-conduct_horizonal(n,n+1);
    L(n^2,n^2-n)=conduct_vertical(n,n);L(n^2,n^2+n+1)=conduct_vertical(n+1,n);L(n^2,n^2-1)=conduct_horizonal(n,n);L(n^2,n^2+n)=conduct_horizonal(n,n+1);

    L(2:n-1,2:n-1)=L(2:n-1,2:n-1)+diag(-conduct_vertical(1,2:n-1)-conduct_vertical(2,2:n-1)-conduct_horizonal(1,2:n-1)-conduct_horizonal(1,3:n));
    L(2:n-1,n^2+3*n+2:n^2+4*n-1)=L(2:n-1,n^2+3*n+2:n^2+4*n-1)+diag(conduct_vertical(1,2:n-1));
    L(2:n-1,n+2:2*n-1)=L(2:n-1,n+2:2*n-1)+diag(conduct_vertical(2,2:n-1));
    L(2:n-1,1:n-2)=L(2:n-1,1:n-2)+diag(conduct_horizonal(1,2:n-1));
    L(2:n-1,3:n)=L(2:n-1,3:n)+diag(conduct_horizonal(1,3:n));

    L(n^2-n+2:n^2-1,n^2-n+2:n^2-1)=L(n^2-n+2:n^2-1,n^2-n+2:n^2-1)+diag(-conduct_vertical(n,2:n-1)-conduct_vertical(n+1,2:n-1)-conduct_horizonal(n,2:n-1)-conduct_horizonal(n,3:n));
    L(n^2-n+2:n^2-1,n^2-2*n+2:n^2-n-1)=L(n^2-n+2:n^2-1,n^2-2*n+2:n^2-n-1)+diag(conduct_vertical(n,2:n-1));
    L(n^2-n+2:n^2-1,n^2+2*n-1:-1:n^2+n+2)=L(n^2-n+2:n^2-1,n^2+2*n-1:-1:n^2+n+2)+diag(conduct_vertical(n+1,2:n-1));
    L(n^2-n+2:n^2-1,n^2-n+1:n^2-2)=L(n^2-n+2:n^2-1,n^2-n+1:n^2-2)+diag(conduct_horizonal(n,2:n-1));
    L(n^2-n+2:n^2-1,n^2-n+3:n^2)=L(n^2-n+2:n^2-1,n^2-n+3:n^2)+diag(conduct_horizonal(n,3:n));

    for i = 2:n-1
        L(i*n-n+1,i*n-n+1)=-conduct_vertical(i,1)-conduct_vertical(i+1,1)-conduct_horizonal(i,1)-conduct_horizonal(i,2);
        L(i*n-n+1,i*n-2*n+1)=conduct_vertical(i,1);
        L(i*n-n+1,i*n+1)=conduct_vertical(i+1,1);
        L(i*n-n+1,n^2+3*n-i+1)=conduct_horizonal(i,1);
        L(i*n-n+1,i*n-n+2)=conduct_horizonal(i,2);

        L(i*n,i*n)=-conduct_vertical(i,n)-conduct_vertical(i+1,n)-conduct_horizonal(i,n)-conduct_horizonal(i,n+1);
        L(i*n,i*n-n)=conduct_vertical(i,n);
        L(i*n,i*n+n)=conduct_vertical(i+1,n);
        L(i*n,i*n-1)=conduct_horizonal(i,n);
        L(i*n,n^2+i)=conduct_horizonal(i,n+1);

        L(i*n-n+2:i*n-1,i*n-n+2:i*n-1)=L(i*n-n+2:i*n-1,i*n-n+2:i*n-1)+diag(-conduct_vertical(i,2:n-1)-conduct_vertical(i+1,2:n-1)-conduct_horizonal(i,2:n-1)-conduct_horizonal(i,3:n));
        L(i*n-n+2:i*n-1,i*n-2*n+2:i*n-n-1)=L(i*n-n+2:i*n-1,i*n-2*n+2:i*n-n-1)+diag(conduct_vertical(i,2:n-1));
        L(i*n-n+2:i*n-1,i*n+2:i*n+n-1)=L(i*n-n+2:i*n-1,i*n+2:i*n+n-1)+diag(conduct_vertical(i+1,2:n-1));  
        L(i*n-n+2:i*n-1,i*n-n+1:i*n-2)=L(i*n-n+2:i*n-1,i*n-n+1:i*n-2)+diag(conduct_horizonal(i,2:n-1));
        L(i*n-n+2:i*n-1,i*n-n+3:i*n)=L(i*n-n+2:i*n-1,i*n-n+3:i*n)+diag(conduct_horizonal(i,3:n));
    end



    phi = vertcat(zeros(n*n,size(phi,2)),phi);
    U = zeros(n^2+4*n,size(phi,2));
    
    
    for  i = 1:size(phi,2)
        err = 1;
    while err > tol
        
        w = (L-diag(vertcat(3*U(1:n^2,i).^2,zeros(4*n,1))))\(phi(:,i)-L*U(:,i)+vertcat(U(1:n^2,i).^3,zeros(4*n,1)));
        err = norm(w,2);
        U(:,i) = U(:,i)+w;

    end
    end
    I = zeros(4*n,size(U,2));
    I(1:n,:)=-(conduct_horizonal(1:n,n+1)).*(U((1:n)*n,:)-U(n^2+1:n^2+n,:));
    I(n+1:2*n,:)=-(conduct_vertical(n+1,n:-1:1)').*(U(n^2:-1:n^2-n+1,:)-U(n^2+1+n:n^2+2*n,:));
    I(2*n+1:3*n,:)=-(conduct_horizonal(n:-1:1,1)).*(U((n:-1:1)*n-n+1,:)-U(n^2+2*n+1:n^2+3*n,:));
    I(3*n+1:4*n,:)=-(conduct_vertical(1,1:n)').*(U(1:n,:)-U(n^2+3*n+1:n^2+4*n,:));



end

