classdef djc_ldv_vitruvian_man_shapes
 % djc_ldv_vitruvian_man_shapes  main class for the Leonardo da Vinci vitruvian Man shapes
 %
 % see djm_ldv_vitruvian_man_shapes for examples of use
 %
 % See also  djc_cc_vitruvian_man_shapes, 
 %
 % (c) Djamal Boukerroui, 2022
 %
 % Revisions
 
    
    properties
        gama
        w
        r
        dx = 1e-4;

    end
    properties (Hidden)
        centerC = [0, 0];
        centerS;
        corners;
        opts;
    end
    
    methods ( Access = 'public' )
        function obj = djc_ldv_vitruvian_man_shapes(w, dx, gama)
            obj.gama = 274/225;
            obj.w = w;
            if nargin  > 1 && ~isempty(dx)
                obj.dx = dx;
            end
            if nargin  == 3 
                if gama >=1 && gama <= 5/4
                obj.gama = gama;
                else
                    warning('Non valid gama value');
                end
            end
            obj.r = obj.w * obj.gama;
            obj = obj.update_square_corners();
        end
		
		
        function [L, L0]= get_square_lengths(obj)
            L(1) = 1;
            L(2) = obj.gama - sqrt(obj.gama^2 -1);
            L(3) = sqrt(obj.gama^2 -1);            
            L(4) = L(3);
            L(5) = 2 - obj.gama - sqrt(obj.gama^2 -1);
            L(6) = 1 - 2*sqrt(obj.gama -1);
            L(7) = 2*sqrt(obj.gama -1);
            L    = L*obj.w;
            L0 = (2- obj.gama) * obj.w;
            
            assert(abs(L(7) + L(6) - obj.w) < 1e-6);
            assert(abs(L(4) + L(5) - L0) < 1e-6);
            assert(abs(L(3) + L(2) - obj.r) < 1e-6);
            assert(abs(L0 + obj.r - 2*obj.w) < 1e-6);            
        end
		
		
        function alpha = get_alphas(obj, L, L0)
            if nargin == 1
                [L, L0]= get_square_lengths(obj);
            end
            tmp = (1 - obj.gama)/sqrt(2)/obj.gama;
            alpha(1) = asin(tmp) - pi/4;
            alpha(2) = - asin(L(4)/obj.r);
            alpha(3) = 0;
            alpha(4) = - alpha(2);
            alpha(5) = acos(-tmp) - pi/4;
            alpha(6) = asin(L0/obj.r);
        end
        
        % Ploting functions
        function plot(obj)
            nbPoints = ceil(2*pi*obj.r/obj.dx);
            theta = linspace(0, 2*pi, nbPoints);
            x = cos(theta) *obj.r + obj.centerC(1); x(end+1) = x(1);
            y = sin(theta)*obj.r + obj.centerC(2); y(end+1) = y(1);
            hc = plot(x, y, 'r', 'Tag', 'VCircle'); hold on;
            
            hs(1) = line([obj.corners(1,1), obj.corners(2,1)],  [obj.corners(1,2), obj.corners(2,2)]); 
            hs(2) = line([obj.corners(2,1), obj.corners(3,1)],  [obj.corners(2,2), obj.corners(3,2)]);
            hs(3) = line([obj.corners(3,1), obj.corners(4,1)],  [obj.corners(3,2), obj.corners(4,2)]);
            hs(4) = line([obj.corners(4,1), obj.corners(1,1)],  [obj.corners(4,2), obj.corners(1,2)]);
            set(hs, 'Color', 'g', 'Tag', 'VSquare');   
            axis image; hold off;
            grid on;
        end
		
        function plot_intersections(obj)
            [L, l0] = get_square_lengths(obj);
            p(1,:) = [L(1), l0 ];
            p(2,:) = [obj.w, L(4)];
            p(3,:) = [obj.w, -L(5)];
            hold on; plot(p(:,1), p(:,2),'+m'); 
            hold off
        end
		
        function plot_alpha_segments(obj)
            R = obj.r*sqrt(2)*1.05;
            alpha = get_alphas(obj);
            ca = cos(alpha)*R;
            sa = sin(alpha)*R;
            hold on;
            for i=1:length(alpha)
                h = line([0 ca(i)], [0 sa(i)]);
                set(h, 'Color', 'b', 'Tag', 'VAngle', 'LineStyle', '--');
            end
            hold off;
        end
        
        % get distances along the segment of the square. See figure xx for
        % parametrisation
        function st = get_da(obj, seg)            
            [L, fds] = get_da_information(obj); 
            Lc =[0 cumsum(L)];
            
            if nargin == 1
                seg = 1:7;
            end            
         
            st = struct();
            for i=1:length(seg)
                 k = seg(i);
                 st(i).seg = k;
                 st(i).h = 0: obj.dx: L(k);
                 st(i).fda = fds{k};
                 st(i).da = st(i).fda(st(i).h);  
                 %These to get the non fliped h and da
                 switch (k)
                     % Non fliped segment just add the origin
                    case {1,4,5}
                        st(i).fget_ch  =  @(x) Lc(k) + x;
                        st(i).fget_cda =  @(x) x;
                       %Fliped segment start from the end and flip it.
                    case {2,3,6,7}                    
                        st(i).fget_ch  =  @(x) Lc(k+1) - x(end:-1:1);
                        st(i).fget_cda =  @(x) x(end:-1:1);                                              
                 end               
            end                            
        end
        
        % Returns segments lenghts and handles for distance functions
        function [L, fds, fds_primitive, fds_inverse, monotony] = get_da_information(obj) 
            [L, L0] = get_square_lengths(obj);
            fds = {...
                @(h) sqrt(h.^2 +obj.r^2) - obj.r;...
                @(h) sqrt((h + L(3)).^2 + obj.w^2) - obj.r; ...
                @(h) obj.r - sqrt(h.^2 + obj.w^2); ...
                @(h) obj.r - sqrt(h.^2 + obj.w^2); ...
                @(h) sqrt((h + L(4)).^2 + obj.w^2) - obj.r ;...
                @(h) sqrt((h + L(7)).^2 + L0^2)  - obj.r ;...
                @(h) obj.r - sqrt(h.^2 + L0^2)};
                
            fds_primitive ={...                     
                @(h)  -obj.r*h + 1/2* h.*sqrt(h.^2 + obj.r^2 )...
                      + 1/2* obj.r^2.*log(h + sqrt(h.^2 + obj.r^2));...                     
                @(h) -obj.r*h + 1/2* (h + L(3)).*sqrt((h + L(3)).^2 + obj.w^2 )...
                     + 1/2* obj.w^2.*log(h + L(3) + sqrt((h + L(3)).^2 + obj.w^2));...
                @(h)  obj.r*h - 1/2* h.*sqrt(h.^2 + obj.w^2 )...
                      - 1/2* obj.w^2.*log(h + sqrt(h.^2 + obj.w^2));...                     
                @(h)  obj.r*h - 1/2* h.*sqrt(h.^2 + obj.w^2 )...
                      - 1/2* obj.w^2.*log(h + sqrt(h.^2 + obj.w^2));...                     
                @(h) -obj.r*h + 1/2* (h + L(4)).*sqrt((h + L(4)).^2 + obj.w^2 )...
                     + 1/2* obj.w^2.*log(h + L(4) + sqrt((h + L(4)).^2 + obj.w^2));...
                @(h) -obj.r*h + 1/2* (h + L(7)).*sqrt((h + L(7)).^2 + L0^2 )...
                     + 1/2* L0^2.*log(h + L(7) + sqrt((h + L(7)).^2 + L0^2));...
                @(h)  obj.r*h - 1/2* h.*sqrt(h.^2 + L0^2 )...
                      - 1/2* L0^2.*log(h + sqrt(h.^2 + L0^2));...
                };

            fds_inverse = {...                                  
                 @(y) sqrt((y + obj.r).^2 - obj.r^2);...
                 @(y) sqrt((y + obj.r).^2 - obj.w^2) - L(3);...
                 @(y) sqrt((y - obj.r).^2 - obj.w^2);...
                 @(y) sqrt((y - obj.r).^2 - obj.w^2);...
                 @(y) sqrt((y + obj.r).^2 - obj.w^2) - L(4);...
                 @(y) sqrt((y + obj.r).^2 - L0^2) - L(7);...
                 @(y) sqrt((y - obj.r).^2 - L0^2);...
                 };
             monotony = [+1 +1 -1 -1 +1 +1 -1]';             
        end
		
        % get distance profiles of segments for -pi/2 <= alpha  <pi/2
        function st = get_db(obj, seg)
            [L, L0] = get_square_lengths(obj);
            alphas = get_alphas(obj);
            if nargin == 1
                seg = 1:7;
            end
            st = struct();
            [a, b, fds] = get_db_information(obj);
            dalpha =  obj.dx/obj.r;
            for i=1:length(seg)
                k = seg(i);  %for clarity
                st(i).fdbb =  @(alpha)  fds{k}(alpha, a(k), b(k));
                st(i).seg = k;
                
                switch (k)
                    case 1                        
                        st(i).fdb = @(alpha)  obj.r*(1 + sin(alpha));                        
                        st(i).alpha = -pi/2:dalpha: alphas(k);                                          
                    case 2
                        st(i).fdb = @(alpha)  obj.w - obj.r*cos(alpha);
                        st(i).alpha = alphas(k-1): dalpha: alphas(k);
                    case 3
                        st(i).fdb = @(alpha)  -(obj.w - obj.r*cos(alpha));
                         st(i).alpha = alphas(k-1): dalpha: alphas(k);
                    case 4
                        st(i).fdb = @(alpha)  -(obj.w - obj.r*cos(alpha));
                         st(i).alpha = alphas(k-1): dalpha: alphas(k);
                    case 5
                        st(i).fdb = @(alpha)  obj.w - obj.r*cos(alpha);
                         st(i).alpha = alphas(k-1): dalpha: alphas(k);
                    case 6
                        st(i).fdb = @(alpha)  L0 - obj.r*sin(alpha);
                         st(i).alpha = alphas(k-1): dalpha: alphas(k);
                    case 7
                        st(i).fdb = @(alpha)  obj.r*sin(alpha) - L0;
                        st(i).alpha = alphas(6):dalpha: pi/2;                       
                end
                st(i).db = st(i).fdb(st(i).alpha);
            end
        end
		
        % This is a compact formulation of the distance functions.
        function [a, b, fds, fds_primitive, fds_inverse, monotony] = get_db_information(obj) 
            a =  [obj.gama, 1, -1, -1, 1 2-obj.gama, -(2-obj.gama)];
            b =  [1, -1, 1, 1, -1, -1, 1]*obj.gama;
            monotony = [+1 -1 +1 -1 +1 -1 +1];
            a =  a*obj.w;
            b =  b*obj.w;
                    
            fds= { ...
                @(x, a, b)  a + b*sin(x);...
                @(x, a, b)  a + b*cos(x); ...
                @(x, a, b)  a + b*cos(x); ...
                @(x, a, b)  a + b*cos(x); ...
                @(x, a, b)  a + b*cos(x); ...
                @(x, a, b)  a + b*sin(x);...
                @(x, a, b)  a + b*sin(x)};
            fds_primitive = { ...
                @(x, a, b)  a*x - b*cos(x);...
                @(x, a, b)  a*x + b*sin(x);...
                @(x, a, b)  a*x + b*sin(x);...
                @(x, a, b)  a*x + b*sin(x);...
                @(x, a, b)  a*x + b*sin(x);...
                @(x, a, b)  a*x - b*cos(x);...
                @(x, a, b)  a*x - b*cos(x)};
            fds_inverse = {...
                @(y, a, b) asin((y - a)/b);...
                @(y, a, b) -acos((y - a)/b);...
                @(y, a, b) -acos((y - a)/b);...
                @(y, a, b) acos((y - a)/b);...
                @(y, a, b) acos((y - a)/b);...
                @(y, a, b) asin((y - a)/b);...
                @(y, a, b) asin((y - a)/b)};
        end
		
        %Using the basis function to compute the Average distance
        function [theta,f_B1, f_B2] = get_db_basis(obj) 
            theta(:,1) =  [obj.gama; 1; -1;-1; 1; 2-obj.gama; -(2-obj.gama)];
            theta(:,2) =  [0; -1; 1; 1; - 1; 0; 0].*obj.gama;
            theta(:,3) =  [1; 0; 0; 0; 0; -1; 1].*obj.gama;
            theta = theta*obj.w;            
                        
            f_B1 = {@(x)  ones(size(x)) ; @cos; @sin};
            f_B2 = {@(x)  x ; @sin; @(x) -cos(x)};            
        end
       
        %CDF  Square to Circle cumulative distribution
        function cdf = cdf_da(obj, d)
            [L, fds, fds_primitive, fds_inverse, monotony] = get_da_information(obj);
            % Here the def interval is always [0 L(i)]
            prior = L/sum(L);
            cdf = zeros(size(d));
            for i = 1:length(L)
                if monotony(i) == 1
                    out = obj.cdf_increasing_uniform(d, 0, L(i), fds{i}, fds_inverse{i});
                else
                    out = obj.cdf_decreasing_uniform(d, 0, L(i), fds{i}, fds_inverse{i});
                end
                cdf = cdf + prior(i)*out;
            end
        end
        
        % Circle to square cumulative distribution
        function cdf = cdf_db(obj, d)
            [a, b, fds, fds_primitive, fds_inverse, monotony] = get_db_information(obj);
                       
            alphas = get_alphas(obj);
            alpha_left = [-pi/2 alphas];  % Interval
            alpha_right = [alphas pi/2];
            prior = (alpha_right - alpha_left)/pi;
            cdf = zeros(size(d));
            for i = 1:length(a)
                fd  = @(x) fds{i}(x, a(i), b(i));
                fdi = @(x) fds_inverse{i}(x, a(i), b(i));                
                if monotony(i) == 1
                    out = obj.cdf_increasing_uniform(d, alpha_left(i), alpha_right(i), fd, fdi);                    
                else
                    out = obj.cdf_decreasing_uniform(d, alpha_left(i), alpha_right(i), fd, fdi); 
                end
                cdf = cdf + prior(i)*out;
            end
        end
        
        % Mean Distance db
        function mu = mean_db(obj)
            [theta, f_B1, f_B2] = get_db_basis(obj);
            alphas = get_alphas(obj);
            alpha_im1 = [-pi/2 alphas];  % Interval
            alpha_i = [alphas pi/2];
            
            B2_im1 = [f_B2{1}(alpha_im1); f_B2{2}(alpha_im1);f_B2{3}(alpha_im1)];
            B2_i = [f_B2{1}(alpha_i); f_B2{2}(alpha_i);f_B2{3}(alpha_i)];
            mu  =  trace(theta*(B2_i - B2_im1))/pi;
        end
        
        % Mean Distance da
        function   mu = mean_da(obj) 
            [L, fds, fds_primitive, fds_inverse, monotony] = get_da_information(obj);
            mu = 0;
            for i=1:length(L)
                mu = mu + fds_primitive{i}(L(i)) -  fds_primitive{i}(0);
            end
            mu = mu / sum(L);
        end
        
        % Hausdorff distance
        function d = hd(obj)
            d = max([sqrt(1+ obj.gama^2)- obj.gama, 2*(obj.gama -1) ]);
            d =  d*obj.w;
        end
                
        %Medians and 95%HD
        %A
        function M = hd50_da(obj)
           M = obj.distance_percentile(.5, 1);
        end
         function M = hd95_da(obj)
           M = obj.distance_percentile( .95, 1);
         end
         %B
         function M = hd50_db(obj)
           M = obj.distance_percentile( .5, 2);
        end
         function M = hd95_db(obj)
           M = obj.distance_percentile( .95, 2);
         end
                
         % Lengths, APL and normalised APL
         function out = get_lengths(obj, epsilon)
             out.lengthA = obj.w*2*4;
             out.lengthB = 2*pi*obj.r;
             for i=1:length(epsilon)
                 out.naplA(i) = 1 - obj.cdf_da(epsilon(i));
                 out.naplB(i) = 1 - obj.cdf_db(epsilon(i));
             end
             out.aplA = out.naplA*out.lengthA;
             out.aplB = out.naplB*out.lengthB;             
         end
         
        % Volumes
        function out = get_vols(obj)
           alphas = get_alphas(obj);
           [L, L0] = get_square_lengths(obj);
           
            out.volA = obj.w^2*4;
            out.volB = pi*obj.r^2;
            out.volABc = 2*obj.w*(obj.r + L0 - L(3))...
                - obj.r^2 * (alphas(2)+ pi/2 + alphas(6) - alphas(4)) ...
                -  L0*L(7);
                        
            out.volBAc = obj.r^2*(pi/2 - alphas(6) +2*alphas(4)) ...
                -L0*L(7) - 2*obj.w*L(4);
            
            out.volAB  = out.volA - out.volABc;
            out.DSC = 2*out.volAB/(out.volA + out.volB);
        end
    end
    
    methods (Access = 'private')
        function obj = update_square_corners(obj)
            obj.centerS = [0, obj.w-obj.r ];
            obj.corners = [-1, +1 ; 1, 1; 1, -1; -1, -1];
            obj.corners = obj.corners.*obj.w;
            obj.corners = bsxfun(@plus, obj.corners, obj.centerS);                               
        end
		
        % proposition 2
        function out = cdf_increasing_uniform(obj, x, a, b, f, finverse)
            out = zeros(size(x));
            out(x >= f(b)) = 1;
            valid =  (x <= f(b)) & (x >= f(a));
            out(valid) = (finverse(x(valid))- a)/(b-a);
        end
		
        % proposition 4
        function out = cdf_decreasing_uniform(obj, x, a, b, f, finverse)
            out = zeros(size(x));
            out(x >= f(a)) = 1;
            valid =  (x <= f(a)) & (x >= f(b));
            out(valid) = 1 - (finverse(x(valid))- a)/(b-a);
        end            
        
        % Needed for Median, APL and 95% HD
        function M = distance_percentile(obj, p, ctype)
             if ctype == 1
                 fx = @(x) (cdf_da(obj, x) - p);
             else
                 fx = @(x) (cdf_db(obj, x) - p);
             end
            x2 = hd(obj);
            %M  =   fminbnd(@(x) abs(fx(x)), x1, x2, obj.opts);
            M = fzero(fx, x2/2);
         end
            
        % Average Distance same as mean_db %not used
        function mu = mean_db0(obj)            
           [a, b, fds, fds_primitive] = get_db_information(obj);
           alphas = get_alphas(obj);
           alpha_left = [-pi/2 alphas];  % Interval
           alpha_right = [alphas pi/2];
                      
           for i = 1:length(a)
               mu(i) =  fdsp{i}(alpha_right(i), a(i), b(i)) - fds_primitive{i}(alpha_left(i), a(i), b(i));
           end
           mu =  sum(mu)/pi;
        end        
    end    
end
