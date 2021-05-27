%alpha=1/6
p=Pr{1,1};
h1=plot(p(:,2),p(:,1),'-sb');hold on %beta=0
p=Pr{1,2};
h2=plot(p(:,2),p(:,1),'-xb');hold on%beta=1
p=Pr{1,3};
h3=plot(p(:,2),p(:,1),'-ob');hold on%beta=2
p=Pr{1,4};
h4=plot(p(:,2),p(:,1),'-db');hold on%beta=3
p=Pr{1,5};
h5=plot(p(:,2),p(:,1),'-^b');hold on%beta=4
p=Pr{1,6};
h6=plot(p(:,2),p(:,1),'-vb');hold on%beta=6
p=Pr{1,7};
h7=plot(p(:,2),p(:,1),'-*b');hold on%beta=12
p=Pr{1,8};
h8=plot(p(:,2),p(:,1),'-.b');hold on%beta=33

%alpha=2/6
p=Pr{2,1};
h11=plot(p(:,2),p(:,1),'-sr');hold on%beta=1/6
p=Pr{2,2};
h21=plot(p(:,2),p(:,1),'-xr');hold on%beta=2/6
p=Pr{2,3};
h31=plot(p(:,2),p(:,1),'-or');hold on%beta=3/6
p=Pr{2,4};
h41=plot(p(:,2),p(:,1),'-dr');hold on%beta=4/6
p=Pr{2,5};
h51=plot(p(:,2),p(:,1),'-^r');hold on%beta=5/6
p=Pr{2,6};
h61=plot(p(:,2),p(:,1),'-vr');hold on%beta=6/6
p=Pr{1,7};
h71=plot(p(:,2),p(:,1),'-*r');hold on%beta=12
p=Pr{1,8};
h81=plot(p(:,2),p(:,1),'-.r');hold on%beta=33

%alpha=3/6
p=Pr{3,1};
h12=plot(p(:,2),p(:,1),'-sg');hold on%beta=1/6
p=Pr{3,2};
h22=plot(p(:,2),p(:,1),'-xg');hold on%beta=2/6
p=Pr{3,3};
h32=plot(p(:,2),p(:,1),'-og');hold on%beta=3/6
p=Pr{3,4};
h42=plot(p(:,2),p(:,1),'-dg');hold on%beta=4/6
p=Pr{3,5};
h52=plot(p(:,2),p(:,1),'-^g');hold on%beta=5/6
p=Pr{3,6};
h62=plot(p(:,2),p(:,1),'-vg');hold on%beta=6/6
p=Pr{1,7};
h72=plot(p(:,2),p(:,1),'-*g');hold on%beta=12
p=Pr{1,8};
h82=plot(p(:,2),p(:,1),'-.g');hold on%beta=33

%alpha=4/6
p=Pr{4,1};
h13=plot(p(:,2),p(:,1),'-sm');hold on%beta=1/6
p=Pr{4,2};
h23=plot(p(:,2),p(:,1),'-xm');hold on%beta=2/6
p=Pr{4,3};
h33=plot(p(:,2),p(:,1),'-om');hold on%beta=3/6
p=Pr{4,4};
h43=plot(p(:,2),p(:,1),'-dm');hold on%beta=4/6
p=Pr{4,5};
h53=plot(p(:,2),p(:,1),'-^m');hold on%beta=5/6
p=Pr{4,6};
h63=plot(p(:,2),p(:,1),'-vm');hold on%beta=6/6
p=Pr{1,7};
h73=plot(p(:,2),p(:,1),'-*m');hold on%beta=12
p=Pr{1,8};
h83=plot(p(:,2),p(:,1),'-.m');hold on%beta=33

%alpha=5/6
p=Pr{5,1};
h14=plot(p(:,2),p(:,1),'-sc');hold on%beta=1/6
p=Pr{5,2};
h24=plot(p(:,2),p(:,1),'-xc');hold on%beta=2/6
p=Pr{5,3};
h34=plot(p(:,2),p(:,1),'-oc');hold on%beta=3/6
p=Pr{5,4};
h44=plot(p(:,2),p(:,1),'-dc');hold on%beta=4/6
p=Pr{5,5};
h54=plot(p(:,2),p(:,1),'-^c');hold on%beta=5/6
p=Pr{5,6};
h64=plot(p(:,2),p(:,1),'-vc');hold on%beta=6/6
p=Pr{1,7};
h74=plot(p(:,2),p(:,1),'-*c');hold on%beta=12
p=Pr{1,8};
h84=plot(p(:,2),p(:,1),'-.c');hold on%beta=33

%alpha=1
p=Pr{6,1};
h15=plot(p(:,2),p(:,1),'-sy');hold on%beta=1/6
p=Pr{6,2};
h25=plot(p(:,2),p(:,1),'-xy');hold on%beta=2/6
p=Pr{6,3};
h35=plot(p(:,2),p(:,1),'-oy');hold on%beta=3/6
p=Pr{6,4};
h45=plot(p(:,2),p(:,1),'-dy');hold on%beta=4/6
p=Pr{6,5};
h55=plot(p(:,2),p(:,1),'-^y');hold on%beta=5/6
p=Pr{6,6};
h65=plot(p(:,2),p(:,1),'-vy');hold on%beta=6/6
p=Pr{1,7};
h75=plot(p(:,2),p(:,1),'-*b');hold on%beta=12
p=Pr{1,8};
h85=plot(p(:,2),p(:,1),'-.b');hold on%beta=33

l1= legend([h1,h2,h3,h4,h5,h6,h7,h8,h11,h21,h31,h41,h51,h61,h71,h81,h12,h22,h32,h42,h52,h62,h72,h82,h13,h23,h33,h43,h53,h63,h73,h83,h14,h24,h34,h44,h54,h64,h74,h84,h15,h25,h35,h45,h55,h65,h75,h85],'$\alpha=1/6$, $\beta=0$','$\alpha=1/6$, $\beta=1$','$\alpha=1/6$, $\beta=2$','$\alpha=1/6$, $\beta=3$','$\alpha=1/6$, $\beta=4$','$\alpha=1/6$, $\beta=6$','$\alpha=1/6$, $\beta=12$','$\alpha=1/6$, $\beta=33$','$\alpha=2/6$, $\beta=0$','$\alpha=2/6$, $\beta=1$','$\alpha=2/6$, $\beta=2$','$\alpha=2/6$, $\beta=3$','$\alpha=2/6$, $\beta=4$','$\alpha=2/6$, $\beta=6$','$\alpha=2/6$, $\beta=12$','$\alpha=2/6$, $\beta=33$','$\alpha=3/6$, $\beta=0$','$\alpha=3/6$, $\beta=1$','$\alpha=3/6$, $\beta=2$','$\alpha=3/6$, $\beta=3$','$\alpha=3/6$, $\beta=4$','$\alpha=3/6$, $\beta=6$','$\alpha=3/6$, $\beta=12$','$\alpha=3/6$, $\beta=33$','$\alpha=4/6$, $\beta=0$','$\alpha=4/6$, $\beta=1$','$\alpha=4/6$, $\beta=2$','$\alpha=4/6$, $\beta=3$','$\alpha=4/6$, $\beta=4$','$\alpha=4/6$, $\beta=6$','$\alpha=4/6$, $\beta=12$','$\alpha=4/6$, $\beta=33$','$\alpha=5/6$, $\beta=0$','$\alpha=5/6$, $\beta=1$','$\alpha=5/6$, $\beta=2$','$\alpha=5/6$, $\beta=3$','$\alpha=5/6$, $\beta=4$','$\alpha=5/6$, $\beta=6$','$\alpha=5/6$, $\beta=12$','$\alpha=5/6$, $\beta=33$','$\alpha=1$, $\beta=0$','$\alpha=1$, $\beta=1$','$\alpha=1$, $\beta=2$','$\alpha=1$, $\beta=3$','$\alpha=1$, $\beta=4$','$\alpha=1$, $\beta=6$','$\alpha=1$, $\beta=12$','$\alpha=1$, $\beta=33$');hold on
set(l1,'interpreter','latex')
%ylim([3 10])

ylabel('Error(%)') % x-axis label
xlabel('Number of lePrned examples')
hold on

