classdef shaft
    properties
        Ixx       % area moment of inertia about the x axis (m4)
        Iyy       % area moment of inertia about the y axis (m4)
        Izz       % area moment of inertia about the z axis (m4)
    end
    properties
        L         % length of the shaft (m)
        ri        % inner diameter of the shaft (m)
        ro        % outer diameter of the shaft (m)
        E         % modulus of elasticity of the shaft (N/m2)
        ks        % stiffness at the trocar (N/m)
        G         % shear modulus of elasticity (N/m2)
        tpF       % the wrench vector at the tip (6x1 N and N.m)
        csF       % the wrench vector at Lp cross section (6x1 N and N.m)
        Ls        % insertion depth into the trocar (m)
        Lp        % location of the cross section from the instrument base (m)
        Nn        % the normalized transducer signals
    end
    properties (Constant)
        % taken from the TRO paper on sensor calibration
        A =[0.0002    0.0142   -0.0001   -0.0002   -0.0000    0.0004;
            0.0063    0.0101   -0.0004   -0.0003    0.0002    0.0000;
            -0.0107   -0.0044    0.0002    0.0001   -0.0002    0.0004;
            -0.0129    0.0021   -0.0000    0.0000   -0.0003   -0.0000;
            0.0096   -0.0076   -0.0001    0.0001    0.0002    0.0004;
            0.0065   -0.0126   -0.0002    0.0002    0.0002    0.0000]*diag([1,1,1,1000,1000,1000]);
    end
    methods
        % constructor
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
        end
        function obj = ft_to_csF(obj,ft,colLs,Lp)
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
                cx = 6*obj.E*obj.Ixx/obj.ks;
                cy = 6*obj.E*obj.Iyy/obj.ks;
                
                gx = ls^2/(cx+2*ls^3);
                gy = ls^2/(cy+2*ls^3);
                
                F11 = 1 - (3*obj.L-ls)*gy;
                F15 = -3*gy;
                
                F22 = 1 - (3*obj.L-ls)*gx;
                F24 = 3*gx;
                
                F42 = (3*obj.L-ls)*(ls-Lp)*gx-(obj.L-Lp);
                F44 = 1-3*(ls-Lp)*gx;
                
                F51 = (obj.L-Lp)-(3*obj.L-ls)*(ls-Lp)*gy;
                F55 = 1-3*(ls-Lp)*gy;
                omega = [0 1;-1 0];
                %{
                B is a 6 by 6 matrix mapping forces at the tip to
                the forces at a cross section located at Lp
                
                F(6x1) = B(6x6) * ft(6x1)
                F(nx6) = ft(nx6) * B(6x6)'
                %}
                B = [F11 0 0 0 F15 0;
                     0 F22 0 F24 0 0;
                     0 0 1 0 0 0;
                     0 F42 0 F44 0 0;
                     F51 0 0 0 F55 0;
                     0 0 0 0 0 1];
                 
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