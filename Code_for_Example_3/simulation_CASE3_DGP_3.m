%simulation_CASE3_DGP:3
%observable model
rng(1234);
%设置模拟轮次，对于每一组(T,N)都进行500轮模拟;
R=500;
%设置模拟的时间跨度;
Time=[100,200,400];
%设置用于构建因子估计量的可观测变量个数;
Number=[100,200,300];
%保存结果，待集中输出
data=zeros(9,2);
result = zeros(9,6);
colNames={'Mean','Standard Deviation'};
colNames1={'RMSE_all','RMSE_cons','RMSE_F1','RMSE_F2','RMSE_W1','RMSE_W2'};
%保留每一轮模拟中AUC的计算结果
auc_output=zeros(9,3);
%分别报告上述九种情况中多轮模拟AUC结果的平均值,中位数和标准差
colnames_auc={'AUC-mean','AUC-median','AUC-std'};
group=0;
item=0;
Betake=zeros(R,5);
tic;
for N=Number
    for T=Time
        item=item+1;
        NOR=zeros(R,6);
        %加入留存每次AUC计算结果的一维数组
        AUC=zeros(R,1);
        for k=1:R
            %Observable variables(用于构造因子估计量的可观测变量)
            W_obv=zeros(T,N);
            %Observable regressors
            W1=zeros(T,1);
            W2=zeros(T,1);
            %Unobservable regressors F1 and F2
            F1=zeros(T,1);
            F2=zeros(T,1);
            %factor loading
            lambda=zeros(N,2);
            %random term in regression model
            E=zeros(T+1,1);
            Ys=zeros(T+1,1);
            Y_interm=zeros(T+1,1);
            F1_0=unifrnd(0,2);
            F2_0=unifrnd(0,2);
            E(1) =  random('Logistic',0,1);
            %generating factor loadings lambda_i,where i=1,2,...,N
            for c=1:N
                lambda(c,1)=unifrnd(0,6);
                lambda(c,2)=unifrnd(0,6);
            end
            for t=1:T
                W1(t)=unifrnd(0,3);
                W2(t)=unifrnd(-3,3);
                %Generating unobservable factor regressors
                if(t<2)
                    F1(t)=0.8*F1_0+sqrt((1-0.8*0.8))*normrnd(0,1);
                    F2(t)=0.64*F2_0+sqrt((1-0.64*0.64))*normrnd(0,1);
                else
                    F1(t)=0.8*F1(t-1)+sqrt((1-0.8*0.8))*normrnd(0,1);
                    F2(t)=0.64*F2(t-1)+sqrt((1-0.64*0.64))*normrnd(0,1);
                end
                % 生成 DGP2 的误差项
                rho_e = 0.7;
                % 对应 DGP2 的 rho 值
                E(t+1) = rho_e * E(t) + sqrt(1 - rho_e^2) * random('Logistic',0,1);  % 使用 AR(1) 模型生成误差项
                %Constructing observable variables
                for j=1:N
                    W_obv(t,j)=normrnd(0,1)+lambda(j,1)*F1(t)+lambda(j,2)*F2(t);
                end
                %Theoretical result
                Ys(t+1)=-2+F1(t)+F2(t)+W1(t)+W2(t)+E(t+1);
                Y_interm(t+1)=(Ys(t+1)>=0);
            end
            W = [F1,F2];
            %实际上，我们用于训练模型的结果是Y数组的第2至第(T+1)个成分
            Y=Y_interm(2:T+1);
            [F,V]=estimateF(W_obv);
            H = V\(F'*W/T)*(lambda'*lambda/N);
            R1=F(:,1);
            R2=F(:,2);
            %用于估计\hat{r}的自变量矩阵
            X=[ones(T,1),R1,R2,W1,W2];
            %优化起始点
            x0=[0,0,0,0,1];
            A=[];
            b=[];
            Aeq=[];
            beq=[];
            %优化取值范围上下限
            lb=[-10,-10,-10,-10,-10];
            ub=[10,10,10,10,10];
            beta=fmincon(@(beta)Fun3(X,Y,beta),x0,A,b,Aeq,beq,lb,ub);
            Betake(k,1)=beta(1);
            Betake(k,2)=beta(2);
            Betake(k,3)=beta(3);
            Betake(k,4)=beta(4);
            Betake(k,5)=beta(5);
            hatr=[beta(1);beta(2);beta(3);beta(4);beta(5)];
            idealr=[-2;H'\[1;1];1;1];
            NOR(k,1) = norm((hatr-idealr), 'fro');
            %得到数值型的预测结果
            Result_forecast=X*hatr;
            Result_forecast=1 ./ (1 + exp(-Result_forecast));
            %使用matlab内置的AUC计算函数计算每次模拟的AUC数值
            [~,~,~,AUC(k)]=perfcurve(Y, Result_forecast, 1);
            for i = 2:6
                NOR(k,i) = hatr(i-1)-idealr(i-1);
            end
        end
        d=3*group+item;
        data(d,1)=mean(NOR(:,1));
        data(d,2)=std(NOR(:,1));
        %留存AUC计算结果
        auc_output(d,1)=mean(AUC);
        auc_output(d,2)=median(AUC);
        auc_output(d,3)=std(AUC);
        for i = 1:6
            result(d,i) = sqrt(mean(NOR(:,i).^2));
        end
    end
    item=0;
    group=group+1;
end
%将结果整理为表格输出
rowNames = {'(N=100,T=100)', '(N=100,T=200)', '(N=100,T=400)','(N=200,T=100)','(N=200,T=200)','(N=200,T=400)','(N=300,T=100)', '(N=300,T=200)', '(N=300,T=400)'};
T = array2table(data, 'RowNames', rowNames, 'VariableNames', colNames);
TT = array2table(result, 'RowNames', rowNames, 'VariableNames', colNames1);
disp(T)
disp(TT)
%AUC输出结果
TTT=array2table(auc_output,"RowNames",rowNames,"VariableNames",colnames_auc);
disp(TTT)
writetable(TT, 'Results_CASE3_DGP_3.xlsx', 'Sheet','CASE3_DGP2','WriteRowNames', true);
writetable(TTT, 'Results_AUC_CASE3_DGP_3.xlsx', 'Sheet','CASE3_DGP2','WriteRowNames', true); 
t1 = toc;  % 计算经过的时间
disp(['经过时间：', num2str(t1), '秒']);
%Function of MLE
function output=Fun3(X,Y,beta)
leng=length(Y);
part1=Y'*log(normcdf(X*[beta(1);beta(2);beta(3);beta(4);beta(5)],0,1));
part2=(ones(leng,1)-Y)'*log(ones(leng,1)-normcdf(X*[beta(1);beta(2);beta(3);beta(4);beta(5)],0,1));
output=-(part1+part2);
end
%Function of Factor estimators
function [F,V]=estimateF(X)
[T,N]=size(X);
CovMatrix = (1/(N*T))*(X*X');
[eigenvectors, eigenvalues] = eig(CovMatrix);
eigenvalues = diag(eigenvalues);
[sortedEigenvalues, idx] = sort(eigenvalues, 'descend');
sortedEigenvectors = eigenvectors(:, idx);
V = diag(sortedEigenvalues(1:2));
F = sortedEigenvectors(:, 1:2);
F = sqrt(T)*F;
end