function djm_ldv_vitruvian_man_shapes(job, opts)
%  djm_ldv_vitruvian_man_shapes main file to generate GT numbers and figures 
%  Leonardo da Vinci Vitruvian Man Shapes
%
%  djm_ldv_vitruvian_man_shapes(job, opt)
% 
%  job : integer value
%       1 : generate figures and results
%       2 : generate json file to be used by the Python shape generator
%           requires JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%       
%  opts : a structure of option for job 1
%       .dx     step length for empirical integration (in mm)
%       .debug  option to compare theoretical results with empirical integration
%       .figsave a handle to a function to save figures  @f(filename)  
%   
%
% see also  djc_ldv_vitruvian_man_shapes
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

gammas = fliplr([1., 1.1 274/225 5/4]);
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
    fields = {'w', 'r', 'gamma', 'hd', 'hd95da', 'hd95db', 'hd50da', 'hd50db', 'mua', 'mub'};
    fmergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
      
    for i=1:length(ws)
        ifig = figsNb(i);
        figure(ifig)  ; clf;   set(gcf, 'unit','normalized'); set(gcf, 'Position', [.1 .1 .3 .5]);
        figure(ifig+1); clf;   set(gcf, 'unit','normalized'); set(gcf, 'Position', [.1 .1 .3 .5]);
        for j= 1:length(gammas)
            w = ws(i);
            gamma = str2num(num2str(gammas(j), '%.10g'));
            r =  ws(i)*gamma;
            
            disp([datestr(now), 'Processing w = ', num2str(w, '%4.2f'), ',   r = ', num2str(r, '%4.2f'),  ',   gamma = ', num2str(gamma, '%4.2f')]);
            
            shape = djc_ldv_vitruvian_man_shapes(w, dx, gamma);
            
            da_st = shape.get_da();
            h_values = []; da = []; % Concatenates all segments
            for k = 1:length(da_st)
                h_values =  [h_values; da_st(k).fget_ch(da_st(k).h)'];
                da =  [da; da_st(k).fget_cda(da_st(k).da)'];
            end
            
            db_st= shape.get_db();
            alpha_values = []; db = []; % Concatenates all segments
            for k = 1:length(db_st)
                alpha_values = [alpha_values; db_st(k).alpha(:)];
                db = [db; db_st(k).db(:)];
            end
            
            hd   = shape.hd();
            %HD95
            hd95da =  shape.hd95_da();
            hd95db =  shape.hd95_db();
            
            %Median
            hd50da =  shape.hd50_da();
            hd50db =  shape.hd50_db();
            
            % Mean
            mua = shape.mean_da();
            mub = shape.mean_db();
            
            % Check theory with numerical values
            if debug
                num_a =  djf_percentile(da, [ 0.5 .95 1]);
                num_b =  djf_percentile(db, [ 0.5 .95 1]);
                
                my_assert(hd, max(num_a(3), num_b(3)), 1e-6, {'hd failed', 'hd ok'});
                
                my_assert(hd95da, num_a(2), 1e-6, {'hd95 Da failed', 'hd95 Da ok'});
                my_assert(hd95db, num_b(2), 1e-6, {'hd95 Db failed', 'hd95 Db ok'});
                
                my_assert(hd50da, num_a(1), 1e-6, {'hd50 Da failed', 'hd50 Da ok'});
                my_assert(hd50db, num_b(1), 1e-6, {'hd50 Db failed', 'hd50 Db ok'});
                
                my_assert(mean(da), mua, dx, {'Average Da Theory-Empi failed', 'Average Da Theory-Empi ok'});
                my_assert(mean(db), mub, dx, {'Average Db Theory-Empi failed', 'Average Db Theory-Empi ok'});
            end
            
            if j == 1
                d_values = -hd*.05:dx:hd*1.05;
            end
            cdfb = shape.cdf_db(d_values);
            cdfa = shape.cdf_da(d_values);
            
            
            % Plots ------------------------------------------------------
            figure(ifig);
            subplot(211); hold on; hand_da(j) = plot(h_values, da, 'Color', colors{j});
            xlabel('h (in mm)'); ylabel('d_a(h) (in mm)')
            subplot(212); hold on; hand_cdfa(j)= plot(d_values, cdfa, 'Color', colors{j});
            xlabel('d (in mm)'); ylabel('P(d_a \leq d)')
            
            figure(ifig+1);
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
        figure(ifig);
        subplot(211); set(gca, 'YMinorTick', 'on'); grid on;
        subplot(212);
        set(gca, 'Ylim', [0 1.1], 'Xlim', [d_values(1), d_values(end)]);
        legend({'1', '2', '3', '4'}, 'Location', 'SouthEast');
        set(gca, 'YTick', 0:.2:1, 'YMinorTick', 'on'); grid on;
        figsave(['vitruvian_profiles_da_w', num2str(w)]);
        
        figure(ifig+1);
        subplot(211); set(gca, 'Xlim', [-1.6 1.6], 'YMinorTick', 'on'); grid on;
        
        subplot(212);
        set(gca, 'Ylim', [0 1.1], 'Xlim', [d_values(1), d_values(end)]);
        legend({'1', '2', '3', '4'}, 'Location', 'SouthEast');
        set(gca, 'YTick', 0:.2:1, 'YMinorTick', 'on'); grid on;
        figsave(['vitruvian_profiles_db_w', num2str(w)]);
    end
    
    %writetable(struct2table(results), 'VitruvianMan_shapePaper_theory_results.txt');
    writetable(struct2table(results), 'LdV_VitruvianMan_shapePaper_theory_results.csv');
end

% Create and save a json file to be used by the python synthetic generator
if job == 2
    [ww, wr] = meshgrid(ws, gammas);
    
    w = ww(:);
    r = w.*wr(:);
    
    namesA = arrayfun(@(x) ['LdV_VM_ShapeA', num2str(x)], 1:length(w), 'UniformOutput', false);
    namesB = arrayfun(@(x) ['LdV_VM_ShapeB', num2str(x)], 1:length(w), 'UniformOutput', false);
    
    VitruvianMan = struct();
    for i=1:length(w)
        VitruvianMan(i).w = w(i);
        VitruvianMan(i).r = r(i);
        VitruvianMan(i).gamma = wr(i);
        VitruvianMan(i).nameA = namesA{i};
        VitruvianMan(i).nameB = namesB{i};
        VitruvianMan(i).theta = 0.0;
    end
    %savejson('CylinderCuboid',CylinderCuboid,'FileName', 'CylinderCuboid.json', 'FloatFormat','\t%3.8f');
    savejson('LdV_VitruvianMan',VitruvianMan,'LdV_VitruvianMan.json');
end

% Just to fuse both json file into a single file.
if job == 3
    a = loadjson('CC_VitruvianMan.json');
    b = loadjson('LdV_VitruvianMan.json');
    a.LdV_VitruvianMan = b.LdV_VitruvianMan;
    savejson('',a,'AllShapes.json');
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
    

