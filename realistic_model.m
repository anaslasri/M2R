% generating our data for ODE
f = @(t,y) [-0.00001*y(1)*y(2);
    0.00001*y(1)*y(2)+0.1*y(3)-0.00002*y(1)*y(2);
    0.00002*y(1)*y(2)-0.1*y(3)];
[t,y] = ode45(f,[0 10],[40000 10000 0]);
z=zeros(85,3);
for i=1:85
    z(i,1)=random('normal',y(i,1),500);
    z(i,2)=random('normal',y(i,2),500);
    z(i,3)=random('normal',y(i,3),500);
end

%diffz gives us derivative of all populations so we can work out constants
%for ODEs
diffz = [0 0 0; diff(z)];

% initial guess at variables
init_var = [0.00001 0.00002 0.1];
new_var = init_var;


% initialising vector
acc_var = [];

%calculating initial loglikelihood
S_ODE = -init_var(1)*((z(:,1).*z(:,2)));
Z_ODE = init_var(1)*z(:,1).*z(:,2)+init_var(3)*z(:,3) -init_var(2)*(z(:,1).*z(:,2));
R_ODE = init_var(2)*(z(:,1).*z(:,2))-init_var(3)*z(:,3);

Old_Lhood = sum(log(normpdf(diffz(:,1),S_ODE,500)))+sum(log(normpdf(diffz(:,2),Z_ODE,500)))+sum(log(normpdf(diffz(:,3),R_ODE,500)));

for i=1:10000
    new_var = [normrnd(0.00001,0.000000000001) new_var(2) new_var(3)];
    S_ODE = -new_var(1)*((z(:,1).*z(:,2)));
    Z_ODE = new_var(1)*z(:,1).*z(:,2)+new_var(3)*z(:,3) -new_var(2)*(z(:,1).*z(:,2));
    R_ODE = new_var(2)*(z(:,1).*z(:,2))-new_var(3)*z(:,3);
    New_Lhood = sum(log(normpdf(diffz(:,1),S_ODE,500)))+sum(log(normpdf(diffz(:,2),Z_ODE,500)))+sum(log(normpdf(diffz(:,3),R_ODE,500)));
    if New_Lhood - Old_Lhood > log(rand)
        acc_var = [acc_var; new_var];
        Old_Lhood = New_Lhood;
    end
    
    new_var = [new_var(1) normrnd(0.00002,0.000000003) new_var(3)];
    S_ODE = -new_var(1)*((z(:,1).*z(:,2)));
    Z_ODE = new_var(1)*z(:,1).*z(:,2)+new_var(3)*z(:,3) -new_var(2)*(z(:,1).*z(:,2));
    R_ODE = new_var(2)*(z(:,1).*z(:,2))-new_var(3)*z(:,3);
    New_Lhood = sum(log(normpdf(diffz(:,1),S_ODE,500)))+sum(log(normpdf(diffz(:,2),Z_ODE,500)))+sum(log(normpdf(diffz(:,3),R_ODE,500)));
    if New_Lhood - Old_Lhood > log(rand)
        acc_var = [acc_var; new_var];
        Old_Lhood = New_Lhood;
    end

    new_var = [new_var(1) new_var(2) normrnd(0.1,0.0000001)];
    S_ODE = -new_var(1)*((z(:,1).*z(:,2)));
    Z_ODE = new_var(1)*z(:,1).*z(:,2)+new_var(3)*z(:,3) -new_var(2)*(z(:,1).*z(:,2));
    R_ODE = new_var(2)*(z(:,1).*z(:,2))-new_var(3)*z(:,3);
    New_Lhood = sum(log(normpdf(diffz(:,1),S_ODE,500)))+sum(log(normpdf(diffz(:,2),Z_ODE,500)))+sum(log(normpdf(diffz(:,3),R_ODE,500)));
    if New_Lhood - Old_Lhood > log(rand)
        acc_var = [acc_var; new_var];
        Old_Lhood = New_Lhood;
    end
end

scatter3(acc_var(:,1), acc_var(:,2), acc_var(:,3))
