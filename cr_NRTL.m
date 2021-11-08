function [x]=cr_NRTL(x10,x20,x30,l)
%%
%para sistema Acido láctico(1)-N-butanol(2)-(3)-Lactato(4)-Water

%x10=composición en la alimentaciónn de A.L
%x20=composición en la alimentaciónn de n-butanol
%x30=composición en la alimentaciónn de Lactato de butilo
%l=l­mide limite de integración

%Para llamar la grafica
tspan=[0 l];

fig=quatplot(1);
hold on

M=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0];

%Opciones (configurar del método de solución DAE´s-ODE 15s)
opts=odeset('Refine',5,'MassSingular','yes','Mass',M,'AbsTol',1e-6,'RelTol',1e-3);


x4=@(x1,x2,x3)1-x1-x2-x3;

c0=[x10; x20; x30];


for k=1:length(x10)
    
    [~,T0]=vle(c0(:,k)'); %para darle una buena temperatura inicial
    
    [~,x] = ode15s(@cambio,tspan,[c0(:,k); T0],opts);
    plot3(x(:,1),x(:,2),x4(x(:,1),x(:,2),x(:,3)),'.k');
    hold on
    [~,x] = ode15s(@cambio,-tspan,[c0(:,k); T0],opts); 
    plot3(x(:,1),x(:,2),x4(x(:,1),x(:,2),x(:,3)),'.k');
    
end

% quatlabel('Lactic Acid','1-Butanol','1-Butyl Lactate','Water',fig);
quatlabel('Lactic Acid','1-Butanol','1-Butyl Lactate','Water',fig);
[~,T]=vle([1 0 0 0]); C = {num2str(T,'%.1f'),'K'};%|Para el Comp1
str = strjoin(C); %|
text(1, 0,-0.1,str,'horizontalalignment','center');
[~,T]=vle([0 1 0 0]);C = {num2str(T,'%.1f'),'K'};%|Para el Comp2
str = strjoin(C);
text(0, 1, -0.1, str, 'HorizontalAlignment', 'center');
[~,T]=vle([0 0 1 0]);C = {num2str(T,'%.1f'),'K'};%|Para el Comp3
str = strjoin(C);
text(0, 0, -0.1, str, 'HorizontalAlignment', 'center');
[~,T]=vle([0 0 0 1]);C = {num2str(T,'%.1f'),'K'};%|Para el Comp4
str = strjoin(C);
text(0,0, 1.1,str,'horizontalalignment','center');%|


hold off
%%
end

function dxdt=cambio(~,x)
%%
T=x(4);
[y,~]=vle([x(1) x(2) x(3)],T); %para conseguir las composiciones en el vapor

dxdt=[x(1)-y(1);
    x(2)-y(2);
    x(3)-y(3);
    1-y(1)-y(2)-y(3)-y(4)];
%%
end

function [y,T]=vle(x,t)
%%
%para sistema Acido láctico(1)-N-butanol(2)-(3)-Lactato(4)-Water

%la x de entrada es:
%x(1)=fracción molar del componente 1
%x(2)=fracción molar del componente 2
%x(3)=fracción molar del componente 3

%la t es la temperatura

c=[x(1) x(2) x(3) 1-x(1)-x(2)-x(3)];

if nargin<2
%     T=(78+80+81+100)/4;
%     T=T+273.15;%T=(Tb(Et)+Tb(MEK)+Tb(ciclohexano)+Tb(agua))/4
 T=((448.07)+(390.90)+(458.25)+(373.17))/4;
%  T=((273.15)+(300)+(300)+(300))/4;
else
    T=t;
end

% P=5000; %presión en Pa
P=101325;
% P=94600;
% P=10000;

%Antoine extendido
C=[7.51107 -1965.70 -91.021 0 0 0 0 %Lactic ácid (k, log, kPa)
    99.3822 -9866.4 0 0 -11.655 1.08e-17 6;%n-butanol (experimental) k,ln,kPa
    7.142 -2160.134 235.46 0 0 0 0;%butyl lactate (experimental) °C,log,kPa
    73.649 -7258.2 0 0 -7.3037 4.1653e-6 2];%water K, Ln, Pa

p_sat1=@(T) (10.^(C(1,1)+(C(1,2)/(T+C(1,3)))+C(1,4)*T+C(1,5)*log(T)+C(1,6)*T^C(1,7)))*1000;%p_sat Lactic ácid
p_sat2=@(T) (exp(C(2,1)+(C(2,2)/(T+C(2,3)))+C(2,4)*T+C(2,5)*log(T)+C(2,6)*T^C(2,7)))*1000;%p_sat n-butanol
p_sat3=@(T) (10.^(C(3,1)+(C(3,2)/((T-273.15)+C(3,3)))+C(3,4)*T+C(3,5)*log(T)+C(3,6)*T^C(3,7)))*1000;%p_sat butyl lactate
p_sat4=@(T) exp(C(4,1)+(C(4,2)/(T+C(4,3)))+C(4,4)*T+C(4,5)*log(T)+C(4,6)*T^C(4,7));%p_sat water


%P en Pa

%vector de actividades:
g=g_NRTL(c,T);
   
vapor1=@(T,g,c)((c(1).*p_sat1(T)).*g(1))./P;
vapor2=@(T,g,c)((c(2).*p_sat2(T)).*g(2))./P;
vapor3=@(T,g,c)((c(3).*p_sat3(T)).*g(3))./P;
vapor4=@(T,g,c)((c(4).*p_sat4(T)).*g(4))./P;

if nargin<2 %es decir, si se utiliza la funciónn vle(x,t) por primera vez para calcular
            %una buena condición inicial
    
    dT=0.1; 
    f=@(T,g,c) vapor1(T,g,c)+vapor2(T,g,c)+vapor3(T,g,c)+vapor4(T,g,c)-1;
    df=@(T,dT,g,c) f(T+dT,g,c)-f(T,g,c);

    while (abs(f(T,g,c))>(1e-8))

        Tn=T-f(T,g,c)/(df(T,dT,g,c)/dT);
        dT=Tn-T;
 
        T=Tn;

        g=g_NRTL(c,T); %se actualizan los valores de las actividades
    
    end 
    
    y=[vapor1(T,g,c), vapor2(T,g,c), vapor3(T,g,c), vapor4(T,g,c)];

else
    
    y=[vapor1(T,g,c), vapor2(T,g,c), vapor3(T,g,c), vapor4(T,g,c)];

end
%%

end

function [g]=g_NRTL(x,T)
%%
%cÃ¡culo de las matrices para el modelo NRTL aplicado al sistema
%sistema Lactic acid(1)-N-butanol(2)-Butyl Lactate(3)-Water(4)

%--------matrices para el cÃ¡lculo de las actividades------%
%Datos para temperatura en K


%%
% NRTL Coeficientes de interacción binaria matriz a
a=[ 
    0 0 0 0;
    0 0 0.791708572 0;
    0 -0.72642125 0 0;
    0 0 0 0;
    ];

% NRTL Coeficientes de interacción binaria matriz b
b = [
     0 -654.6 327.590125 -148.07 ;
     1537.6 0 -158.78125 -274.66;
     -130.300539 134.427168 0 0; 
     -128.98 1442.4 0 0;
    ];

% NRTL Coeficientes de interacción binaria matriz alpha
alfa = [
      0 0.2 0.3 0.2;
      0.2 0 0.3 0.2;
      0.3  0.3 0 0;
      0.2 0.2  0 0
    ];
%%
%matriz tau
tau=a+(1/T).*b; 

%matriz de Gibbs
G=exp(-alfa.*tau);

%---matrices de la descomposiciÃ³n de la expresiÃ³n de NRTL para el cÃ¡culo
%del logaritmo neperiano de la actividad.Ver pÃ¡gina 29 del cuaderno para
%nomenclatura
A=[(dot(tau(:,1),(G(:,1).*x')))/dot(G(:,1),x);
    (dot(tau(:,2),(G(:,2).*x')))/dot(G(:,2),x);
    (dot(tau(:,3),(G(:,3).*x')))/dot(G(:,3),x);
    (dot(tau(:,4),(G(:,4).*x')))/dot(G(:,4),x)];

%ahora para la matriz B

c11=tau(1,1)-dot(tau(:,1),x'.*G(:,1))/dot(G(:,1),x);
c12=tau(1,2)-dot(tau(:,2),x'.*G(:,2))/dot(G(:,2),x);
c13=tau(1,3)-dot(tau(:,3),x'.*G(:,3))/dot(G(:,3),x);
c14=tau(1,4)-dot(tau(:,4),x'.*G(:,4))/dot(G(:,4),x);

c21=tau(2,1)-dot(tau(:,1),x'.*G(:,1))/dot(G(:,1),x);
c22=tau(2,2)-dot(tau(:,2),x'.*G(:,2))/dot(G(:,2),x);
c23=tau(2,3)-dot(tau(:,3),x'.*G(:,3))/dot(G(:,3),x);
c24=tau(2,4)-dot(tau(:,4),x'.*G(:,4))/dot(G(:,4),x);

c31=tau(3,1)-dot(tau(:,1),x'.*G(:,1))/dot(G(:,1),x);
c32=tau(3,2)-dot(tau(:,2),x'.*G(:,2))/dot(G(:,2),x);
c33=tau(3,3)-dot(tau(:,3),x'.*G(:,3))/dot(G(:,3),x);
c34=tau(3,4)-dot(tau(:,4),x'.*G(:,4))/dot(G(:,4),x);

c41=tau(4,1)-dot(tau(:,1),x'.*G(:,1))/dot(G(:,1),x);
c42=tau(4,2)-dot(tau(:,2),x'.*G(:,2))/dot(G(:,2),x);
c43=tau(4,3)-dot(tau(:,3),x'.*G(:,3))/dot(G(:,3),x);
c44=tau(4,4)-dot(tau(:,4),x'.*G(:,4))/dot(G(:,4),x);

C=[c11 c12 c13 c14; c21 c22 c23 c24; c31 c32 c33 c34; c41 c42 c43 c44];

D=zeros(4,4);

D(:,1)=G(:,1)./dot(G(:,1),x);
D(:,2)=G(:,2)./dot(G(:,2),x);
D(:,3)=G(:,3)./dot(G(:,3),x);
D(:,4)=G(:,4)./dot(G(:,4),x);

B(1)=dot((D(1,:).*x),C(1,:));
B(2)=dot((D(2,:).*x),C(2,:));
B(3)=dot((D(3,:).*x),C(3,:));
B(4)=dot((D(4,:).*x),C(4,:));

g=[exp(A(1)+B(1)), exp(A(2)+B(2)),exp(A(3)+B(3)),exp(A(4)+B(4))];
%%

end

%%------------funciones para representar los diagramas cuaternarios---------%%
function [fig]=quatplot(n)
fig = figure(n);
fig.Color='white';
axis off
ax=gca;
V=[0 0 0; 1 0 0; 0 1 0; 0 0 1]; %V=[V1; V2; V3; V4];
F=[1 2 3; 1 2 4; 1 3 4; 2 3 4]; %F=[Cara1; Cara2;..]=[V1-V2-V3; V1-V2-V4...]

patch('Vertices',V,'Faces',F,'FaceAlpha',0);

view(76, 28);

view(76, 28);
ax.CameraTarget=[0 0 0];
% text(0.5,-0.1,'0.5');
% text(-0.1,0.5,'0.5');
% text(-0.1,-0.1,0.5,'0.5');

end

function  []=quatlabel(A,B,C,D,fig)
%los nomnbres tienen que ir de mÃ¡s ligero a mÃ¡s pesado

N = {'Sistema',A,'-',B,'-',C,'-',D};
N = strjoin(N);
fig.Name=N;


text(1, 0,-0.05, A, 'HorizontalAlignment', 'center')
text(0, 1, -0.05, B, 'HorizontalAlignment', 'center')
text(0, 0, -0.05, C, 'HorizontalAlignment', 'center')
text(0, 0, 1.05, D, 'HorizontalAlignment', 'center')

end
