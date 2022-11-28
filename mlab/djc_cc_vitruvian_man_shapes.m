classdef djc_cc_vitruvian_man_shapes
 % djc_cc_vitruvian_man_shapes  main class for the vitruvian_man shapes
 % centred Cylinder vs Cube shapes
 %
 % see djm_cc_vitruvian_man_shapes for examples of use
 %
 % See also  djc_ldc_vitruvian_man_shapes, 
 %
 % (c) Djamal Boukerroui, 2022
 %
 % Revisions
    
    properties
        w
        r
        dx = 1e-4;        
    end
    properties (Hidden)
        opts;
    end
    
    methods ( Access = 'public' )
        function obj = djc_cc_vitruvian_man_shapes(w, r, dx)
            obj.w = w;
            obj.r = r;
            
            if w < r/sqrt(2) - 1e-6
                error('sqrt(2) w >= r failed');
            end
            if nargin  > 2
                obj.dx = dx;
            end
        end
        function out = get_case(obj)
            if obj.w < obj.r/sqrt(2)- 1e-6
                error('sqrt(2) w >= r  failed');
            end
            if abs( obj.w*sqrt(2) - obj.r)<1e-6
                out =  1;
            elseif obj.w >= obj.r
                out = 2;
            else
                out = 3;
            end
        end
        function h0 = get_h0(obj)
            if obj.w >= obj.r
                h0 = 0;
            else
                h0 = sqrt(obj.r*obj.r - obj.w*obj.w);
            end
        end
        
        function alpha0 = get_alpha0(obj)
            if obj.w >= obj.r
                alpha0 = 0;
            else
                alpha0 = acos(obj.w/obj.r);
            end
        end
        
        function [da, h] = get_da(obj, h)
            if nargin == 1
                h= 0: obj.dx: obj.w;
            end
            da = abs(sqrt(obj.w*obj.w + h.*h) - obj.r);
        end
        
        function [db, alpha] = get_db(obj, alpha)
            if nargin == 1
                alpha= 0: obj.dx/obj.r: pi/4;
            end
            db = abs(obj.w - obj.r*cos(alpha));
        end
        
        %CDF da
        function cdf = cdf_da(obj, d)
            valid =  d > 0;
            dv =  d(valid);
            cdf1  = obj.cdf_uniform_square( (obj.r + dv).^2 - obj.w*obj.w, obj.w);
            cdf2  = obj.cdf_uniform_square( (obj.r - dv).^2 - obj.w*obj.w, obj.w);
            
            cdf = zeros(size(d));
            cdf(valid) = cdf1 - cdf2;
        end
        
        function cdf = cdf_da_case_1(obj, d)
            assert(abs(obj.w*sqrt(2) - obj.r) < 1e-6, ' w and r values do not correspond to case 1');
            
            cdf = zeros(size(d));
            valid = (d >= 0) & (d <= (sqrt(2) - 1)*obj.w);
            cdf(valid) = 1 - 1/obj.w*sqrt( (sqrt(2)*obj.w - d).^2 - obj.w*obj.w);
            cdf (d >= (sqrt(2) -1)*obj.w) = 1;
            
        end
        function cdf = cdf_da_case_2(obj, d)
            assert(obj.w > obj.r, ' w and r values do not correspond to case 2 w >= r');
            cdf = zeros(size(d));
            valid = (d >= obj.w - obj.r) & (d <= sqrt(2)*obj.w - obj.r );
            cdf(valid) =  1/obj.w*sqrt( (obj.r + d).^2 - obj.w*obj.w);
            cdf (d >= (sqrt(2)*obj.w - obj.r)) = 1;
        end
        
        %CDF db
        function cdf = cdf_db(obj, d)
            valid =  d > 0;
            dv =  d(valid);
            cdf1  = obj.cdf_cos_uniform( (obj.w + dv)./obj.r);
            cdf2  = obj.cdf_cos_uniform( (obj.w - dv)./obj.r);
            
            cdf = zeros(size(d));
            cdf(valid) = cdf1 - cdf2;
        end
        function cdf = cdf_db_case_1(obj, d)
            valid =  d > 0;
            dv =  d(valid);
            cdf1  = obj.cdf_cos_uniform( (obj.w + dv)./obj.r);
            
            cdf = zeros(size(d));
            cdf(valid) = cdf1;
        end
        function cdf = cdf_db_case_2(obj, d)
            valid =  d > 0;
            dv =  d(valid);
            cdf1  = obj.cdf_cos_uniform( (obj.w - dv)./obj.r);
            
            cdf = zeros(size(d));
            cdf(valid) = 1 - cdf1;
        end
        
        % AverageDistance
        function mu = mean_da(obj)
            I1 = @(h) 1/2*h.*sqrt(h.*h + obj.w^2) + 1/2 * obj.w^2 *log( h + sqrt(h.*h + obj.w^2)) - obj.r*h;
            h0 = obj.get_h0();
            mu = 1/obj.w*(I1(obj.w) + I1(0) - 2*I1(h0));
        end
        function mu = mean_da_case_1(obj)
            if obj.get_case() == 1
                mu =  (sqrt(2) - log(1 + sqrt(2)))*obj.w/2;
            else
                error('obj is not in case 1');
            end
        end
        function mu = mean_da_case_2(obj)
            if obj.get_case() == 2
                mu =  (sqrt(2) + log(1 + sqrt(2)))*obj.w/2 - obj.r;
            else
                error('obj is not in case 2');
            end
        end
        function mu = mean_db(obj)
            I2 = @(alpha) obj.w*alpha - obj.r*sin(alpha);
            alpha0 = obj.get_alpha0();
            mu = 4/pi*(I2(pi/4) + I2(0) - 2*I2(alpha0));
        end
        function mu = mean_db_case_1(obj)
            if obj.get_case() == 1
                mu =  obj.w * (4/pi -1);
            else
                error('obj is not in case 1');
            end
        end
        function mu = mean_db_case_2(obj)
            if obj.get_case() == 2
                mu =  obj.w - 2*sqrt(2)/pi * obj.r;
            else
                error('obj is not in case 2');
            end
        end
        
        function d = hd(obj)
            d = max([abs(obj.w - obj.r), abs(sqrt(2)*obj.w - obj.r), abs(obj.w - obj.r/sqrt(2))]);
        end
        
        %Medians
        %A
        function M = hd50_da(obj)
            M = obj.distance_percentile( .5, 1);
        end
        function M = hd50_da_case_1(obj)
            if obj.get_case() == 1
                M  =  obj.w* (sqrt(2) - sqrt(5)/2);
            else
                error('obj is not in case 1');
            end
        end
        function M = hd50_da_case_2(obj)
            if obj.get_case() == 2
                M  =  sqrt(5)/2*obj.w - obj.r;
            else
                error('obj is not in case 2');
            end
        end
        %B
        function M = hd50_db(obj)
            M = obj.distance_percentile( .5, 2);
        end
        function M = hd50_db_case_1(obj)
            if obj.get_case() == 1
                M = obj.w*(sqrt(2) * cos(pi/8) -1);
            else
                error('obj is not in case 1');
            end
        end
        function M = hd50_db_case_2(obj)
            if obj.get_case() == 2
                M = obj.w - obj.r*cos(pi/8);
            else
                error('obj is not in case 1');
            end
        end
        
        %95Hausdorff
        function M = hd95_da(obj)
            M = obj.distance_percentile( .95, 1);
        end
        function M = hd95_da_case_1(obj)
            if obj.get_case() == 1
                M  =  obj.r - obj.w*sqrt(1 + 0.05^2);
            else
                error('obj is not in case 1');
            end
        end
        function M = hd95_da_case_2(obj)
            if obj.get_case() == 2
                M  =  obj.w*sqrt(1+ 0.95^2) - obj.r;
            else
                error('obj is not in case 2');
            end
        end
        
        %B
        function M = hd95_db(obj)
            M = obj.distance_percentile( .95, 2);
        end
        function M = hd95_db_case_1(obj)
            if obj.get_case() == 1
                M = obj.r*cos(0.05*pi/4) - obj.w;
            else
                error('obj is not in case 1');
            end
        end
        function M = hd95_db_case_2(obj)
            if obj.get_case() == 2
                M = obj.w - obj.r*cos(0.95*pi/4);
            else
                error('obj is not in case 2');
            end
        end

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
         
        function out = get_vols(obj)
            h0 = get_h0(obj);
            alpha0 = get_alpha0(obj);
            out.volA = (1/2*obj.w^2)*8;
            out.volB = (pi/8*obj.r^2)*8;
            out.volABc = 1/2*(obj.w^2 - h0*obj.w - (pi/4-alpha0)*obj.r^2)*8;
            out.volBAc = 1/2*(alpha0*obj.r^2 - h0*obj.w)*8;
            out.volAB  = out.volA - out.volABc;
            out.DSC = 2*((h0*obj.w)/obj.r^2 + (pi/4 - alpha0))/((obj.w/obj.r)^2 + pi/4);
        end
    end
    
    methods (Access = 'private')
        function out = cdf_cos_uniform(obj, x)
            out = zeros(size(x));
            out(x >= 1) = 1;
            valid =  (x < 1) & (x > 1/sqrt(2));
            out(valid) = 1 - 4/pi*acos(x(valid));
        end
        function out = cdf_uniform_square(obj, h, w)
            out = zeros(size(h));
            out(h > w*w) = 1;
            valid = (h < w*w) & (h >0);
            out(valid ) = sqrt(h(valid))/w;
        end
        function M = distance_percentile(obj, p, ctype)
             if ctype == 1
                 fx = @(x) (cdf_da(obj, x) - p);
             else
                 fx = @(x) (cdf_db(obj, x) - p);
             end
            x2 = hd(obj);
            %M  =   fminbnd(@(x) abs(fx(x)), 0, x2, obj.opts);
            M = fzero(fx, x2/2);
         end        
    end    
end
