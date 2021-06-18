function [ t1,t2 ] = ex2_7( M )
% matlab book example 2.7
%     M = 500;
    matlabpool('open','local')
    N = 100000;
    a = rand(M,N);
    b1 = zeros(N,1);
    display(strcat('datasize:',num2str(N * M / 1024 / 1024),'M doubles'));
    
    tic;
    parfor kk = 1:N
        b1(kk) = sum(sqrt(a(:,kk)));
    end
    t1 = toc;
    display(strcat('parfor:',num2str(t1)));
    matlabpool close;
    
    tic
    b2 = zeros(N,1);
    for kk = 1:N
        b2(kk) = sum(sqrt(a(:,kk)));
    end
    t2 = toc;
    display(strcat('for:',num2str(t2)));
    

end

