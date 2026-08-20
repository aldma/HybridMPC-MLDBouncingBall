function plotting(T_p,T_f,u_f,x_1,x_2,u_)
j=[];
j(1)=1;
for i=1:T_p*T_f
    if u_f(i)==1
        j(i+1)=j(i);
    else
        j(i+1)=j(i)+1;
    end
end
t0=0:(T_p*(T_f-1));
t=0:(T_p*(T_f));


subplot(4,1,1)
plotflows(t,j,x_1)
grid on
ylabel('$\hat{x}_1$','interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
xlim([0 t(end)])

subplot(4,1,2)
plotflows(t,j,x_2)
grid on
ylabel('$\hat{x}_2$','interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
xlim([0 t(end)])

subplot(4,1,3)
stairs(t0,u_f,'b-','LineWidth', 1)
ylabel('$\rho(u_{f1},u_{f2})$','interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
ylim([-.02,1.02])
xlim([0 t(end)])
xlabel('$k$','interpreter','latex')

subplot(4,1,4)
stairs(t0,u_,'b-','LineWidth', 1) % plot
xlim([0 t(end)])
xlabel('$\ell$','interpreter','latex')
ylabel('$u$','interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)

end