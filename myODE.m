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

plot(t,y(:,3))
hold on 
scatter(t,z(:,3))
hold off