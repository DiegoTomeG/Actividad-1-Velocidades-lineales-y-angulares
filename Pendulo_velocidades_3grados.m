%Limpieza de pantalla
clear all
close all
clc

%Declaración de variables simbólicas (no tienen valor especifico)
syms th1(t) l1 
syms th2(t) l2
syms th3(t) l3

RP=[0 0 0]; %configuracion del robot, 0 para junta rotacional y 1 para prismatica 

%Vector de coordenadas articulares
Q = [th1 th2 th3]; 
disp('Coordenadas articulares'); 
pretty(Q);

%vector de velocidad articulares
Qp = diff(Q,t); %se utiliza diff para derivadas cuya variable no depende de otras
disp('Velocidades articulares'); 
pretty(Qp); 

%numero de grado de libertad del robot
GDL = size(RP,2); %se coloca 2 para referirse a las columnas de la matriz 
GDL_str = num2str(GDL); %convertir numero a string


%Articulacion 1 
%posicion de la junta 1 respecto a 0
P(:,:,1) = [l1*cos(th1);
            l1*sin(th1); 
            0];  %vector de posicion indexado por pagina

%Articulación 2
%posicion de la junta 2 respecto a 1
P(:,:,2) = [l2*cos(th2);
            l2*sin(th2); 
            0];  %vector de posicion indexado por pagina


%Articulación 3
%posicion de la junta 3 respecto a 2
P(:,:,3) = [l3*cos(th3);
            l3*sin(th3); 
            0];  %vector de posicion indexado por pagina

%matriz de rotacion de articulacion 1
R(:,:,1) = [cos(th1) -sin(th1) 0; %analisis de robot pendulo
            sin(th1) cos(th1)  0; 
            0        0         1];    

%matriz de rotacion de articulacion 2
R(:,:,2) = [cos(th2) -sin(th2) 0; %analisis de robot pendulo
            sin(th2) cos(th2)  0; 
            0        0         1];    

%matriz de rotacion de articulacion 3
R(:,:,3) = [cos(th3) -sin(th3) 0; %analisis de robot pendulo
            sin(th3) cos(th3)  0; 
            0        0         1];   

Vector_Zeros = zeros(1,3); 

%inicializamos matriz de transformacion homogenea local
A(:,:,GDL) = simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);

%inicializamos matriz de transformacion homogenea global ()
T(:,:,GDL) = simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);

PO(:,:,GDL)= P(:,:,GDL); %vectores de posicion vistos desde el marco de referencia inercial

RO(:,:,GDL)= R(:,:,GDL); %matrices de rotacion distas desde el marco de referencia inercial

for i= 1: GDL
    i_str = num2str(i);
    %locales
    disp(strcat('Matriz de transformacion local A',i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
    pretty (A(:,:,i));

    %globales
    try 
        T(:,:,i) = T(:,:,i-1)*A(:,:,i);
    catch 
        T(:,:,i) = A(:,:,i); %caso especifico cuando i= 1 nos marcaria error en try
    end 
    disp(strcat('Matriz de transformacion global T',i_str));
    T(:,:,i) = simplify(T(:,:,i));
    pretty(T(:,:,i));

    %obtenemos la matriz de rotacion RO y el vector de translacion PO de la
    %matriz de transformacion homogenea global T(:,:,GDL)

    RO(:,:,i) = T(1:3,1:3,i);
    PO(:,:,i) = T(1:3,4,i);
    pretty(RO(:,:,i));
    pretty(PO(:,:,i));
end 

%calculamos el jacobiano lineal de forma diferencial
disp('Jacobiano lineal obtenido de forma diferencial');
%derivadas parciales de x con respecto a th1 y th2
Jv11= functionalDerivative(PO(1,1,GDL), th1);
Jv12= functionalDerivative(PO(1,1,GDL), th2);
Jv13= functionalDerivative(PO(1,1,GDL), th3);
%derivadas de y con respecto a th1 y th2
Jv21= functionalDerivative(PO(2,1,GDL), th1);
Jv22= functionalDerivative(PO(2,1,GDL), th2);
Jv23= functionalDerivative(PO(2,1,GDL), th3);
%parciales de z con respecto a th1 y th2
Jv31= functionalDerivative(PO(3,1,GDL), th1);
Jv32= functionalDerivative(PO(3,1,GDL), th2);
Jv33= functionalDerivative(PO(3,1,GDL), th3);

%creamos la matriz del jacobiano lineal
jv_d = simplify([Jv11 Jv12; ...
                Jv21 Jv22 ; ...
                Jv31 Jv32]); 

pretty(jv_d)


%------------------------------------------------
%Calculamos el jacobiano lineal y angular de forma analitica
%inicializamod jacobianos analiticos (lineal y angular)
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);

for k=1:GDL
    if((RP(k)==0)|(RP(k)==1))

        try
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1)); %Producto cruz
            Jw_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));
            Jw_a(:,k)=[0,0,1]; %No hay matriz de rotacion anterior, se imprime la matriz de identidad
        end
    else
        %articulaciones prismaticas
       try
           Jv_a(:,k)= RO(:,3,k-1);
       catch
           Jv_a(:,k)=[0,0,1];
       end
           Jw_a(:,k)=[0,0,0];
    end
end

Jv_a= simplify (Jv_a);
Jw_a= simplify (Jw_a);
disp('Jacobiano lineal obtenido de forma analítica');
pretty (Jv_a);
disp('Jacobiano ángular obtenido de forma analítica');
pretty (Jw_a);

disp('Velocidad lineal obtenida mediante el Jacobiano Lineal');
V=simplify (Jv_a*Qp');
pretty(V);
disp('Velocidad angular obtenido mediante el Jacobiano angular')
W=simplify (Jw_a*Qp');
pretty(W);
