function [gradL_condhorz,gradL_condvert] = generateGrad(n,boundayrate)
%GENERATEGRAD 此处显示有关此函数的摘要
%   此处显示详细说明
    gradL_condhorz = zeros(n^2+4*n,n^2+4*n,n,n+1);
    gradL_condvert = zeros(n^2+4*n,n^2+4*n,n+1,n);
    
    J=n+1;
parfor id = 1:n
    for j=1:J
    conduct_horizonal= zeros(n,n+1); %%水平电导分布
    conduct_vertical= zeros(n+1,n); %%垂直电导分布
    conduct_horizonal(id,j)=1;

    L=zeros(n^2+4*n,n^2+4*n);
    L(n^2+1:n^2+n,n^2+1:n^2+n)=diag(-conduct_horizonal(1:n,n+1)*boundayrate);
    L(n^2+1:n^2+n,(1:n).*n)=diag(conduct_horizonal(1:n,n+1)*boundayrate);
    L(n^2+n+1:n^2+2*n,n^2+n+1:n^2+2*n)=diag(-conduct_vertical(n+1,n:-1:1)*boundayrate);
    L(n^2+n+1:n^2+2*n,n^2:-1:n^2-n+1)=diag(conduct_vertical(n+1,n:-1:1)*boundayrate);
    L(n^2+2*n+1:n^2+3*n,n^2+2*n+1:n^2+3*n)=diag(-conduct_horizonal(n:-1:1,1)*boundayrate);
    L(n^2+2*n+1:n^2+3*n,(n:-1:1).*n-(n-1))=diag(conduct_horizonal(n:-1:1,1)*boundayrate);
    L(n^2+3*n+1:n^2+4*n,n^2+3*n+1:n^2+4*n)=diag(-conduct_vertical(1,1:n)*boundayrate);
    L(n^2+3*n+1:n^2+4*n,1:n)=diag(conduct_vertical(1,1:n)*boundayrate);

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
        gradL_condhorz(:,:,id,j) = L;
    end
end
parfor id = 1:n+1
    for j=1:n
    conduct_horizonal= zeros(n,n+1); %%水平电导分布
    conduct_vertical= zeros(n+1,n); %%垂直电导分布
    conduct_vertical(id,j)=1;

    L=zeros(n^2+4*n,n^2+4*n);
    L(n^2+1:n^2+n,n^2+1:n^2+n)=diag(-conduct_horizonal(1:n,n+1)*boundayrate);
    L(n^2+1:n^2+n,(1:n).*n)=diag(conduct_horizonal(1:n,n+1)*boundayrate);
    L(n^2+n+1:n^2+2*n,n^2+n+1:n^2+2*n)=diag(-conduct_vertical(n+1,n:-1:1)*boundayrate);
    L(n^2+n+1:n^2+2*n,n^2:-1:n^2-n+1)=diag(conduct_vertical(n+1,n:-1:1)*boundayrate);
    L(n^2+2*n+1:n^2+3*n,n^2+2*n+1:n^2+3*n)=diag(-conduct_horizonal(n:-1:1,1)*boundayrate);
    L(n^2+2*n+1:n^2+3*n,(n:-1:1).*n-(n-1))=diag(conduct_horizonal(n:-1:1,1)*boundayrate);
    L(n^2+3*n+1:n^2+4*n,n^2+3*n+1:n^2+4*n)=diag(-conduct_vertical(1,1:n)*boundayrate);
    L(n^2+3*n+1:n^2+4*n,1:n)=diag(conduct_vertical(1,1:n)*boundayrate);

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
        gradL_condvert(:,:,id,j)=L
    end
end
gradL_condhorz=sparse(reshape(gradL_condhorz,(n^2+4*n)^2,n*(n+1)));
gradL_condvert=sparse(reshape(gradL_condvert,(n^2+4*n)^2,n*(n+1)));
end

