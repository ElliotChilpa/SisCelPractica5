% Elaboró Carlos E Nuñez C 
clc
clear all
close all

disp("Se considerara que la celda central esta en la coordena 0,0 y que el radio es de 2 km")
num_usuario=input("Ingrese el numero de usuarios que dentra cada celda: ");

r=2000;
apotema=sqrt(3).*r/2;


colores(1,:)="#0072BD";
colores(2,:)="#D95319";
colores(3,:)="#EDB120";
colores(4,:)="#7E2F8E";
colores(5,:)="#77AC30";
colores(6,:)="#4DBEEE";
colores(7,:)="#A2142F";


[Centros_x Centros_y]=Obtener_Centros_de_Celdas(r,apotema);

[Hexagonos_x Hexagonos_y]=Crear_Hexagonos(Centros_x, Centros_y,r);

[Usuarios_x, Usuarios_y] = Crear_Usuarios(Hexagonos_x, Hexagonos_y, num_usuario);

figure()
for i=7:-1:1
plot(Hexagonos_x(i,:),Hexagonos_y(i,:),'LineWidth',2,'Color',colores(i,:))
grid on
hold on
plot(Usuarios_x(i,:),Usuarios_y(i,:),'.','LineWidth',2,'Color',colores(i,:))
end
title('Usuarios de cada estacion base. Pérdidas por distancia únicamente')

contador=1;
for alpha= [4,6,8,10];

Potencias_de_usuarios = Calcular_Potencias(Centros_x, Centros_y, Usuarios_x, Usuarios_y,alpha);


figure()
for i=7:-1:1
plot(Hexagonos_x(i,:),Hexagonos_y(i,:),'LineWidth',2,'Color',colores(i,:))
grid on
hold on
end
idx = find( Potencias_de_usuarios(:,10)==1);
Tabla_celda_0{contador}=Potencias_de_usuarios(idx,:); 

%primera columna = coordena x
%segunda columna = coordena y
%tercera columna = Potencia recibida de la estacion base 0
%cuarta columna =  Potencia recibida de la estacion base 1
%quinta columna =  Potencia recibida de la estacion base 2
%sexta columna =   Potencia recibida de la estacion base 3
%septima columna = Potencia recibida de la estacion base 4
%octava columna =  Potencia recibida de la estacion base 5
%novema columna =  Potencia recibida de la estacion base 6
%decima columna =  Numero de la estacion base que provee mayor potencia




plot(Tabla_celda_0{contador}(:,1),Tabla_celda_0{contador}(:,2),'.','LineWidth',2,'Color',colores(1,:))

title(strcat('Usuarios de la estacion base central, modelo lognormal.  \alpha = ',' ',num2str(alpha)));

figure()
for i=7:-1:1
plot(Hexagonos_x(i,:),Hexagonos_y(i,:),'LineWidth',2,'Color',colores(i,:))
grid on
hold on
idx = find( Potencias_de_usuarios(:,10)==i);
Potencias_graficar=Potencias_de_usuarios(idx,:);
plot(Potencias_graficar(:,1),Potencias_graficar(:,2),'.','LineWidth',2,'Color',colores(i,:))
end
title(strcat('Usuarios de cada estacion base, modelo lognormal. \alpha = ',' ',num2str(alpha)))
contador=contador+1;
end

function [Centros_x, Centros_y]= Obtener_Centros_de_Celdas(r,apotema)
Centros_x(1)=0;
Centros_y(1)=0;

Centros_x(2)=3*r/2;
Centros_y(2)=apotema;

Centros_x(3)=0;
Centros_y(3)=2*apotema;

Centros_x(4)=-3*r/2;
Centros_y(4)=apotema;

Centros_x(5)=-3*r/2;
Centros_y(5)=-apotema;

Centros_x(6)=0;
Centros_y(6)=-2*apotema;

Centros_x(7)=3*r/2;
Centros_y(7)=-apotema;
end

function [Hexagonos_x, Hexagonos_y] = Crear_Hexagonos(Centros_x, Centros_y,r)
t = linspace(0,2*pi,7);
for i=1:length(Centros_x)
Hexagonos_x(i,:) = Centros_x(i)+r*cos(t);
Hexagonos_y(i,:) = Centros_y(i)+r*sin(t);
end
end


function [Usuarios_x, Usuarios_y] = Crear_Usuarios(Hexagonos_x, Hexagonos_y, num_usuario)

for i=1:length(Hexagonos_x)
%En general, puede generar números aleatorios N en el intervalo (a,b) con la fórmula r = a + (b-a).*rand(N,1).
usuario_x_aux = min(Hexagonos_x(i,:)) + ((max(Hexagonos_x(i,:))-(min(Hexagonos_x(i,:))))).*rand(num_usuario,1);
usuario_y_aux = min(Hexagonos_y(i,:)) + ((max(Hexagonos_y(i,:))-(min(Hexagonos_y(i,:))))).*rand(num_usuario,1);
p = inpolygon(usuario_x_aux,usuario_y_aux,Hexagonos_x(i,:),Hexagonos_y(i,:));
usuario_x_aux = usuario_x_aux(p);
usuario_y_aux = usuario_y_aux(p);

while(length(usuario_x_aux)< num_usuario)

usuario_x_aux_2 = min(Hexagonos_x(i,:)) + ((max(Hexagonos_x(i,:))-(min(Hexagonos_x(i,:))))).*rand(1,1);
usuario_y_aux_2 = min(Hexagonos_y(i,:)) + ((max(Hexagonos_y(i,:))-(min(Hexagonos_y(i,:))))).*rand(1,1);

    while(inpolygon(usuario_x_aux_2,usuario_y_aux_2,Hexagonos_x(i,:),Hexagonos_y(i,:))==0)
usuario_x_aux_2 = min(Hexagonos_x(i,:)) + ((max(Hexagonos_x(i,:))-(min(Hexagonos_x(i,:))))).*rand(1,1);
usuario_y_aux_2 = min(Hexagonos_y(i,:)) + ((max(Hexagonos_y(i,:))-(min(Hexagonos_y(i,:))))).*rand(1,1);
    end

usuario_x_aux(end+1)=usuario_x_aux_2;
usuario_y_aux(end+1)=usuario_y_aux_2;
end

Usuarios_x(i,:)=usuario_x_aux;
Usuarios_y(i,:)=usuario_y_aux;

end
end




function Potencias_de_usuarios = Calcular_Potencias(Centros_x, Centros_y, Usuarios_x, Usuarios_y,alpha)
P_tx=10;
P_tx_db = 10*log10(P_tx*1000);
G_tx_db = 12;
G_rx_db = 2;
desviacion_estandar= 7;
mu=0;
constante_del_modelo=P_tx_db + G_tx_db + G_rx_db;

for i=1:length(Centros_x)
    for j=1:length(Usuarios_x(i,:))
distancia=sqrt((Centros_x-Usuarios_x(i,j)).^2+(Centros_y-Usuarios_y(i,j)).^2);

Potencia_aux=constante_del_modelo-10*alpha*log10(distancia)-normrnd(mu,desviacion_estandar,[1,7]);

[Potencia_max, Indice]=max(Potencia_aux);

        if(i==1 && j==1)
Potencias_de_usuarios(1,:)=[Usuarios_x(i,j),Usuarios_y(i,j),Potencia_aux, Indice];
        else
Potencias_de_usuarios(end+1,:)=[Usuarios_x(i,j),Usuarios_y(i,j),Potencia_aux, Indice];   
        end
    end
end
end


