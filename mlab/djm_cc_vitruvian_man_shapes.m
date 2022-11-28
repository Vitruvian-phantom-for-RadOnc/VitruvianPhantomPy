function djm_cc_vitruvian_man_shapes(job, opts)
%  djm_cc_vitruvian_man_shapes  main file to generate GT numbers and figures
%  Centred Cylinder vs Cube Shapes
%
%  djm_cc_vitruvian_man_shapes(job, opt)
%
%  job : integer value
%       1 : generate figures and results
%       2 : generate json file to be used by the Python shape generator
%           requires JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
% 
%
%  opts : a structure of option for job 1
%       .dx     step length for empirical integration (in mm)
%       .debug  option to compare theoretical results with empirical integration
%       .figsave a handle to a function to save figures  @f(filename)
%
%
% see also  djm_cc_vitruvian_man_shapes
%
% (c) Djamal Boukerroui, 2022
%
% Revisions

if nargin == 0
    job =1;
end

if nargin < 2
    opts = struct();
end

if isempty(which('djf_save_fig'))
    figsave1 = @(x) djf_save_fig(x, [], fullfile(getenv('MATLAB_SAVE_FIGURES'), 'ShapePaper'));
else
    figsave1 = [];
end
figsave2 = @(x) x;

%  Needed only for job 1
dx      = getoptions(opts, 'dx', 1e-4);
debug   = getoptions(opts, 'debug', false);
figsave = getoptions(opts, 'figsave', figsave1);


wr_ratios = fliplr([1/sqrt(2),  sqrt(pi)/2,  1, sqrt(2)]);
ws = [10 20 40 80];

% keep all figures if not saving figs
if isempty(figsave)
    figsave = figsave2;
    figsNb = 1:2:2*length(ws);
else
    figsNb = ones(size(ws));
end

if job == 1
    colors = {'r', 'm', 'b', 'c'};
    
    results = [];
    fields = {'w', 'r', 'hd', 'hd95da', 'hd95db', 'hd50da', 'hd50db', 'mua', 'mub'};
    fmergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
    
    for i=1:length(ws)
        ifig = figsNb(i);
        figure(ifig)  ; clf;   set(gcf, 'unit','normalized'); set(gcf, 'Position', [.1 .1 .3 .5]);
        figure(ifig+1); clf;   set(gcf, 'unit','normalized'); set(gcf, 'Position', [.1 .1 .3 .5]);
        
        for j= 1:length(wr_ratios)
            w = ws(i);
            r =  ws(i)./wr_ratios(j);
            r = str2num(num2str(r, '%.10g'));
            
            disp([datestr(now), 'Processing w = ', num2str(w, '%4.2f'), ',   r = ', num2str(r, '%4.2f')]);
            
            shape = djc_cc_vitruvian_man_shapes(w, r, dx);
            
            [da, h_values] = shape.get_da();
            [db, alpha_values] = shape.get_db();
            ca   = shape.get_case();
            hd   = shape.hd();
            
            % HD95
            hd95da_num =  shape.hd95_da();
            hd95db_num =  shape.hd95_db();
            % Check theory with numerical values
            if ca == 1
                hd95da =  shape.hd95_da_case_1();
                hd95db =  shape.hd95_db_case_1();
                if debug
                    my_assert(hd95da, hd95da_num, 1e-6, {'hd95 Da case 1 failed', 'hd95 Da case 1 ok'});
                    my_assert(hd95db, hd95db_num, 1e-6, {'hd95 Db case 1 failed', 'hd95 Db case 1 ok'});
                end
            elseif ca == 2
                hd95da =  shape.hd95_da_case_2();
                hd95db =  shape.hd95_db_case_2();
                if debug
                    my_assert(hd95da, hd95da_num, 1e-6, {'hd95 Da case 2 failed', 'hd95 Da case 2 ok'});
                    my_assert(hd95db, hd95db_num, 1e-6, {'hd95 Db case 2 failed', 'hd95 Db case 2 ok'});
                end
            else
                hd95da = hd95da_num;
                hd95db = hd95db_num;
            end
            
            % Median
            hd50da_num =  shape.hd50_da();
            hd50db_num =  shape.hd50_db();
            % Check theory with numerical values
            if ca == 1
                hd50da =  shape.hd50_da_case_1();
                hd50db =  shape.hd50_db_case_1();
                
                if debug
                    my_assert(hd50da, hd50da_num, 1e-6, {'hd50 Da case 1 failed', 'hd50 Da case 1 ok'});
                    my_assert(hd50db, hd50db_num, 1e-6, {'hd50 Db case 1 failed', 'hd50 Db case 1 ok'});
                end
            elseif ca == 2
                hd50da =  shape.hd50_da_case_2();
                hd50db =  shape.hd50_db_case_2();
                if debug
                    my_assert(hd50da, hd50da_num, 1e-6, {'hd50 Da case 2 failed', 'hd50 Da case 2 ok'});
                    my_assert(hd50db, hd50db_num, 1e-6, {'hd50 Db case 2 failed', 'hd50 Db case 2 ok'});
                end
            else
                hd50da = hd50da_num;
                hd50db = hd50db_num;
            end
            
            % Mean
            mua = shape.mean_da();
            mub = shape.mean_db();
            if debug
                % Check simplified expressions
                if ca == 1
                    my_assert(shape.mean_da_case_1(), mua, 1e-6, {'Average Da case 1 failed', 'Average Da case 1 ok'});
                    my_assert(shape.mean_db_case_1(), mub, 1e-6, {'Average Db case 1 failed', 'Average Db case 1 ok'});
                elseif ca == 2
                    my_assert(shape.mean_da_case_2(), mua, 1e-6, {'Average Da case 2 failed', 'Average Da case 2 ok'});
                    my_assert(shape.mean_db_case_2(), mub, 1e-6, {'Average Db case 2 failed', 'Average Db case 2 ok'});
                end
                % Check empirical estimates with theory
                my_assert(mean(da), mua, dx, {'Average Da Theory-Empi failed', 'Average Da Theory-Empi ok'});
                my_assert(mean(db), mub, dx, {'Average Db Theory-Empi failed', 'Average Db Theory-Empi ok'});
            end
            
            if j == 1
                d_values = -1:.01:hd + 1;
            end
            cdfb = shape.cdf_db(d_values);
            cdfa = shape.cdf_da(d_values);
            
            
            % Plots  -----------------------------------------------------
            figure(ifig)
            subplot(211); hold on; hand_da(j) = plot(h_values, da, 'Color', colors{j});
            xlabel('h (in mm)'); ylabel('d_a(h) (in mm)')
            subplot(212); hold on; hand_cdfa(j)= plot(d_values, cdfa, 'Color', colors{j});
            xlabel('d (in mm)'); ylabel('P(d_a \leq d)')
            
            figure(ifig+1)
            subplot(211);hold on; hand_db(j) = plot(alpha_values, db, 'Color', colors{j});
            xlabel('\alpha'); ylabel('d_b(\alpha) (in mm)')
            subplot(212);hold on; hand_cdfb(j)= plot(d_values, cdfb, 'Color', colors{j});
            xlabel('d (in mm)'); ylabel('P(d_b \leq d)')
            
            % Volumes and Dice --------------------------------------------
            volumes = shape.get_vols(); out = [];
            for f=1:length(fields)
                out.(fields{f}) = eval(fields{f});
            end
            
            %  APL and Lengths   ------------------------------------------
            lengths = shape.get_lengths([1.0 2.0]);
            apls.lengthA = lengths.lengthA;
            apls.nAPLA1mm  = lengths.naplA(1);
            apls.nAPLA2mm  = lengths.naplA(2);
            apls.lengthB = lengths.lengthB;
            apls.nAPLB1mm  = lengths.naplB(1);
            apls.nAPLB2mm  = lengths.naplB(2);
            
            out = fmergestructs(out, volumes);
            out = fmergestructs(out, apls);
            if isempty(results)
                results = out;
            else
                results(end+1) = out;
            end
            %
        end
        figure(ifig)
        subplot(211); set(gca, 'YMinorTick', 'on'); grid on;
        subplot(212);  set(gca, 'Ylim', [0 1.1]);
        legend({'1', '2', '3', '4'}, 'Location', 'SouthEast');
        set(gca, 'YTick', 0:.2:1, 'YMinorTick', 'on'); grid on;
        figsave(['profiles_da_w', num2str(w)]);
        
        figure(ifig+1)
        subplot(211); set(gca, 'YMinorTick', 'on'); grid on;
        subplot(212);  set(gca, 'Ylim', [0 1.1]);
        legend({'1', '2', '3', '4'}, 'Location', 'SouthEast');
        set(gca, 'YTick', 0:.2:1, 'YMinorTick', 'on'); grid on;
        figsave(['profiles_db_w', num2str(w)]);
    end
    
    writetable(struct2table(results), 'CC_VitruvianMan_shapePaper_theory_results.csv');
end

if job == 2
    [ww, wr] = meshgrid(ws, wr_ratios);
    
    w = ww(:);
    r = w./wr(:);
    
    namesA = arrayfun(@(x) ['CC_VM_ShapeA', num2str(x)], 1:length(w), 'UniformOutput', false);
    namesB = arrayfun(@(x) ['CC_VM_ShapeB', num2str(x)], 1:length(w), 'UniformOutput', false);
    
    CylinderCuboid = struct();
    for i=1:length(w)
        CylinderCuboid(i).w = w(i);
        CylinderCuboid(i).r = r(i);
        CylinderCuboid(i).nameA = namesA{i};
        CylinderCuboid(i).nameB = namesB{i};
        CylinderCuboid(i).theta = 0.0;
    end
    %savejson('CylinderCuboid',CylinderCuboid,'FileName', 'CylinderCuboid.json', 'FloatFormat','\t%3.8f');
    savejson('CC_VitruvianMan',CylinderCuboid,'CC_VitruvianMan.json');
end
end


function my_assert(a, b, tol, message)
tmp =  abs(a - b) > tol;
if tmp
    disp( message{1});
else
    disp( message{2});
end
end


function v = getoptions(opt, name, v)
if isfield(opt, name)
    v =opt.(name);
end
end




