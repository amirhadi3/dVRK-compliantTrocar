classdef shaft
    properties
        Ixx
        Iyy
        Izz
        Area
    end
    properties
        L
        ri
        ro
        E
        ks
        G
        tpF
        csF
        Ls
        Lp
        Nn
    end
    properties (Constant)
        A =[0.0002    0.0142   -0.0001   -0.0002   -0.0000    0.0004;
            0.0063    0.0101   -0.0004   -0.0003    0.0002    0.0000;
            -0.0107   -0.0044    0.0002    0.0001   -0.0002    0.0004;
            -0.0129    0.0021   -0.0000    0.0000   -0.0003   -0.0000;
            0.0096   -0.0076   -0.0001    0.0001    0.0002    0.0004;
            0.0065   -0.0126   -0.0002    0.0002    0.0002    0.0000];
    end
    methods
        function obj = shaft(length,innerRadi,outerRadi,modElasticity,...
                shModElasticity,trocarStiffness)
            obj.L = length;
            obj.ri = innerRadi;
            obj.ro = outerRadi;
            obj.E = modElasticity;
            obj.G = shModElasticity;
            obj.ks = trocarStiffness;
            
            obj.Ixx = 1/4*pi*(obj.ro^4-obj.ri^4);
            obj.Iyy = 1/4*pi*(obj.ro^4-obj.ri^4);
            obj.Izz = 1/2*pi*(obj.ro^4-obj.ri^4);
            obj.Area = pi*(obj.ro^2-obj.ri^2);
        end
        function obj = ft_to_F(obj,ft,colLs,Lp)
            % n is the number of measrement points
            % "ft" is an nx6 matrix of forces at the tip at every measurement
            % point
            % "colLs" is an nx1 column vector
            % Lp is the location at which the sensor is mounted
            obj.Ls = colLs;
            obj.tpF = ft;
            obj.Lp = Lp;
            for counter = 1:numel(colLs)
                ls = obj.Ls(counter);
                c = 6*obj.E*obj.Ixx/obj.ks;
                g = ls^2/(c+2*ls^3);
                F11 = 1 - (3*obj.L-ls)*g;
                F15 = -3*g;
                F42 = (3*obj.L-ls)*(ls-Lp)*g-(obj.L-Lp);
                F44 = 1-3*(ls-Lp)*g;
                omega = [0 1;-1 0];
                %{
                B is a 6 by 6 matrix mapping forces at the tip to
                the forces at a cross section located at Lp
                
                F(6x1) = B(6x6) * ft(6x1)
                F(nx6) = ft(nx6) * B(6x6)'
                %}
                B = [F11*eye(2) zeros(2,1) F15*omega zeros(2,1);
                    zeros(1,2) 1 zeros(1,3);
                    F42*omega zeros(2,1) F44*eye(2) zeros(2,1);
                    zeros(1,5) 1];
                
                obj.csF(counter,:) = ft(counter,:)*B.';
            end
        end
        function obj = F_to_Nn(obj)
            % Nn is an nx6 matrix
            % F is an Nx6 matrix
            % A is a 6x6 matrix where Nn = A*F
            obj.Nn = obj.csF*obj.A.';
        end
    end
end