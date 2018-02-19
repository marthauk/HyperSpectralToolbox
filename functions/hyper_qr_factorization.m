function [Q,R,Q_built_in, R_built_in] = hyper_qr_factorization(M)
%HYPER_QR_FACTORIZATION Summary of this function goes here
%   Detailed explanation goes here

%MATLABS inbuilt qr_factorization function; returning an upper triangular
%matrix R and a unitary matrix Q such that A = Q*R
[Q_built_in, R_built_in] = qr(M);

% for verification:Verify that A = Q*R using isAlways: isAlways(A == Q*R)

[m,n]=size(M);
R=M; %Start with R=A
Q=eye(m); %Set Q as the identity matrix
% m= P_BANDS in my case....
for k=1:m-1
    x=zeros(m,1);
    x(k:m,1)=R(k:m,k);
    g=norm(x);
    v=x; v(k)=x(k)+g;
    %Orthogonal transformation matrix that eliminates one element
    %below the diagonal of the matrix it is post-multiplying:
    s=norm(v);
    if s~=0 
        w=v/s; u=2*R'*w;
        R=R-w*u'; %Product HR
        Q=Q-2*Q*w*w'; %Product QR
    end
end

end

