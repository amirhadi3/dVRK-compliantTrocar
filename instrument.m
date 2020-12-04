classdef instrument
    properties
        L         % length of the instrument (m)
        cx        %
        cy        %
        Lp        % sensor position (m)
        Los       % offset in Ls
        A         % calibration matrix
        wt
    end
    methods
        % constructor
        function obj = instrument(x,m)
            obj.A = reshape(x(1:6*m),6,m);
            obj.L = x(6*m+1);
            obj.cx = x(6*m+2);
            obj.cy = x(6*m+3);
            obj.Lp = x(6*m+4);
            obj.Los = x(6*m+5);
        end
        function obj = Nn_to_wt(obj,Nn,possig)
            csf = Nn*obj.A';                %% cross-section force = n*m*m*6 = n*6
            
            Ls = obj.Los-possig;
            gx = Ls.^2./(obj.cx+2*Ls.^3);
            gy = Ls.^2./(obj.cy+2*Ls.^3);
            
            F11 = 1-(3*obj.L-Ls).*gy;
            F15 = -3*gy;
            F22 = 1-(3*obj.L-Ls).*gx;
            F24 = 3*gx;
            F42 = (3*obj.L-Ls).*(Ls-obj.Lp).*gx-(obj.L-obj.Lp);
            F44 = 1-3*(Ls-obj.Lp).*gx;
            F51 = (obj.L-obj.Lp)-(3*obj.L-Ls).*(Ls-obj.Lp).*gy;
            F55 = 1-3*(Ls-obj.Lp).*gy;
            
            det1 = F11.*F55-F15.*F51;
            det2 = F22.*F44-F24.*F42;
            
            obj.wt(:,1) = (1./det1).*(F55.*csf(:,1)-F15.*csf(:,5));
            obj.wt(:,5) = (1./det1).*(-F51.*csf(:,1)+F11.*csf(:,5));
            obj.wt(:,3) = csf(:,3);
            obj.wt(:,2) = (1./det2).*(F44.*csf(:,2)-F24.*csf(:,4));
            obj.wt(:,4) = (1./det2).*(-F42.*csf(:,2)+F22.*csf(:,4));
            obj.wt(:,6) = csf(:,6);
        end
    end
end