%Generating data
f = @(t,y) [-0.00001*y(1)*y(2);
    0.00001*y(1)*y(2)+0.1*y(3)-0.00002*y(1)*y(2);
    0.00002*y(1)*y(2)-0.1*y(3)];
[t,y] = ode45(f,[linspace(0,10,50)],[40000 10000 0]);
z=zeros(50,3);
for i=1:50
    z(i,1)=random('normal',y(i,1),500);
    z(i,2)=random('normal',y(i,2),500);
    z(i,3)=random('normal',y(i,3),500);
end
%z(:,1)=s(t) data for human population
%z(:,2)=z(t) data for zombie population
%z(:,3)=r(t) data for removed zombie population

%True value for alpha, beta, gamma
true_alphabetagamma = [0.00002, 0.00001, 0.1];

%Initialising likelihood
old_lhood = sum(log(normpdf(z(:,1),y(:,1), 500)))+sum(log(normpdf(z(:,3),y(:,3), 500)))+sum(log(normpdf(z(:,2),y(:,2), 500)));

%Initialising variables 
old_alphabetagamma = true_alphabetagamma;


%Initialising vector for accepted sampling of parameters
acc_alphabetagamma = [];


%MCMC Iteration using Gibbs Samplings
for i=1:10000
     new_alphabetagamma=[old_alphabetagamma(1), random('normal', old_alphabetagamma(2),0.0000003), old_alphabetagamma(3)];
     f = @(t_new,y_new) [-new_alphabetagamma(2)*y_new(1)*y_new(2);
     new_alphabetagamma(2)*y_new(1)*y_new(2)+new_alphabetagamma(3)*y_new(3)-new_alphabetagamma(1)*y_new(1)*y_new(2);
     new_alphabetagamma(1)*y_new(1)*y_new(2)-new_alphabetagamma(3)*y_new(3)];
     [t,y_new] = ode45(f,[linspace(0,10,50)],[40000 10000 0]);
     new_lhood = sum(log(normpdf(z(:,1),y_new(:,1), 500)))+sum(log(normpdf(z(:,3),y_new(:,3), 500)))+sum(log(normpdf(z(:,2),y_new(:,2), 500)));
     u=random('uniform', 0, 1);
     if new_lhood-old_lhood>log(u)
         acc_alphabetagamma=[acc_alphabetagamma; new_alphabetagamma];
         old_lhood=new_lhood;
         old_alphabetagamma=new_alphabetagamma;
     end
     
     new_alphabetagamma=[random('normal', old_alphabetagamma(1),0.0000003),old_alphabetagamma(2),old_alphabetagamma(3)];
     f = @(t_new,y_new) [-new_alphabetagamma(2)*y_new(1)*y_new(2);
     new_alphabetagamma(2)*y_new(1)*y_new(2)+new_alphabetagamma(3)*y_new(3)-new_alphabetagamma(1)*y_new(1)*y_new(2);
     new_alphabetagamma(1)*y_new(1)*y_new(2)-new_alphabetagamma(3)*y_new(3)];
     [t_new,y_new] = ode45(f,[linspace(0,10,50)],[40000 10000 0]);
     new_lhood = sum(log(normpdf(z(:,1),y_new(:,1), 500)))+sum(log(normpdf(z(:,3),y_new(:,3), 500)))+sum(log(normpdf(z(:,2),y_new(:,2), 500)));
     u=random('uniform', 0, 1);
     if new_lhood-old_lhood>log(u)
         acc_alphabetagamma=[acc_alphabetagamma; new_alphabetagamma];
         old_lhood=new_lhood;
         old_alphabetagamma=new_alphabetagamma;
     end
     
     new_alphabetagamma=[old_alphabetagamma(1),  old_alphabetagamma(2), random('normal',old_alphabetagamma(3),0.0000003)];
     f = @(t_new,y_new) [-new_alphabetagamma(2)*y_new(1)*y_new(2);
     new_alphabetagamma(2)*y_new(1)*y_new(2)+new_alphabetagamma(3)*y_new(3)-new_alphabetagamma(1)*y_new(1)*y_new(2);
     new_alphabetagamma(1)*y_new(1)*y_new(2)-new_alphabetagamma(3)*y_new(3)];
     [t_new,y_new] = ode45(f,[linspace(0,10,50)],[40000 10000 0]);
     new_lhood = sum(log(normpdf(z(:,1),y_new(:,1), 500)))+sum(log(normpdf(z(:,3),y_new(:,3), 500)))+sum(log(normpdf(z(:,2),y_new(:,2), 500)));
     u=random('uniform', 0, 1);
     if new_lhood-old_lhood>log(u)
         acc_alphabetagamma=[acc_alphabetagamma; new_alphabetagamma];
         old_lhood=new_lhood;
         old_alphabetagamma=new_alphabetagamma;
     end
end

%Plot for accepted alpha, beta, gamma
figure
subplot(2,2,1)
scatter3(acc_alphabetagamma(:,1),acc_alphabetagamma(:,2),acc_alphabetagamma(:,3))
%Plot for distribution of alpha, beta, gamma separately 
rng 'default'
subplot(2,2,2)
histfit(acc_alphabetagamma(:,1),50,'normal')
subplot(2,2,3)
histfit(acc_alphabetagamma(:,2),50,'normal')
subplot(2,2,4)
histfit(acc_alphabetagamma(:,3),50,'normal')

