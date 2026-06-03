%% ============================================================
%  Speckle DIC — GRISHMA Project, IIT Kharagpur
%  Prof. Shibayan Roy | Materials Science Center
%  MATLAB R2026a | Image Processing Toolbox
%
%  Method:
%   - Coarse NCC search (normxcorr2)
%   - Parabolic sub-pixel refinement
%   - Grip-referenced displacement
%   - LSQ window strain calculation
%   - Full-field strain maps + summary plots
% ============================================================

clc; clear; close all;

%% ── USER SETTINGS ───────────────────────────────────────────

IMAGE_FOLDER   = '/MATLAB Drive';        % folder with images
IMAGE_PATTERN  = 'img*.JPG';             % filename pattern
RESULTS_FOLDER = '/MATLAB Drive/dic_results'; % output folder

% ROI — gauge section only (x1, y1, x2, y2)
ROI = [973, 853, 1596, 982];

% DIC parameters
SUBSET_SIZE  = 29;     % subset size in pixels (odd number)
STEP_SIZE    = 5;      % grid spacing in pixels
MAX_SHIFT    = 25;     % coarse search range in pixels

% Post-processing
STRAIN_WINDOW  = 5;    % LSQ strain window (grid points, odd)
STRAIN_SMOOTH  = 1.2;  % Gaussian smoothing sigma (grid points)
OUTLIER_NSIGMA = 2.5;  % outlier rejection threshold
MIN_CORR       = 0.70; % minimum correlation threshold
REF_GRIP_COLS  = 3;    % leftmost grid columns = fixed grip

% Physical scale (set if known)
GAUGE_LENGTH_MM = [];  % e.g. 5.0 — leave empty if unknown

%% ── SETUP ───────────────────────────────────────────────────

% Create results folder
if ~exist(RESULTS_FOLDER, 'dir')
    mkdir(RESULTS_FOLDER);
end

% Get image list
files = dir(fullfile(IMAGE_FOLDER, IMAGE_PATTERN));
files = sort_nat({files.name});   % natural sort: img1,img2,...img10
nFiles = length(files);

fprintf('=======================================================\n');
fprintf('  DIC Pipeline — GRISHMA Project, IIT Kharagpur\n');
fprintf('=======================================================\n');
fprintf('Found      : %d images\n', nFiles);
fprintf('Reference  : %s\n', files{1});
fprintf('ROI        : [%d %d %d %d] → %d × %d px\n', ...
        ROI(1),ROI(2),ROI(3),ROI(4), ...
        ROI(3)-ROI(1), ROI(4)-ROI(2));
fprintf('Subset     : %dpx | Step: %dpx | Max shift: %dpx\n\n', ...
        SUBSET_SIZE, STEP_SIZE, MAX_SHIFT);

%% ── LOAD REFERENCE IMAGE ────────────────────────────────────

ref_raw  = load_gray(fullfile(IMAGE_FOLDER, files{1}));
ref_img  = crop_roi(ref_raw, ROI);
[h_roi, w_roi] = size(ref_img);

% Build DIC grid
half = floor(SUBSET_SIZE/2);
xs   = (half + MAX_SHIFT + 1) : STEP_SIZE : (w_roi - half - MAX_SHIFT);
ys   = (half + MAX_SHIFT + 1) : STEP_SIZE : (h_roi - half - MAX_SHIFT);
[grid_x, grid_y] = meshgrid(xs, ys);
[ny, nx] = size(grid_x);

fprintf('Grid size  : %d × %d = %d points\n\n', nx, ny, nx*ny);

%% ── RESULT STORAGE ──────────────────────────────────────────

rec_frame    = zeros(nFiles-1, 1);
rec_corr     = zeros(nFiles-1, 1);
rec_mean_exx = zeros(nFiles-1, 1);
rec_mean_eyy = zeros(nFiles-1, 1);
rec_max_exx  = zeros(nFiles-1, 1);
rec_max_ux   = zeros(nFiles-1, 1);
rec_nu       = nan(nFiles-1, 1);

%% ── MAIN LOOP ───────────────────────────────────────────────

for idx = 2 : nFiles
    fname   = files{idx};
    frameNo = idx - 1;

    fprintf('Frame %03d/%03d  %-16s', frameNo, nFiles-1, fname);

    % Load deformed image
    def_raw = load_gray(fullfile(IMAGE_FOLDER, fname));
    def_img = crop_roi(def_raw, ROI);

    % Normalize both images
    ref_n = normalize_img(ref_img);
    def_n = normalize_img(def_img);

    % Allocate displacement and correlation arrays
    disp_x = zeros(ny, nx);
    disp_y = zeros(ny, nx);
    corr   = zeros(ny, nx);

    % ── Correlate each grid point ────────────────────────────
    for i = 1:ny
        for j = 1:nx
            cx = grid_x(i,j);
            cy = grid_y(i,j);

            % Extract reference subset
            ref_sub = ref_n(cy-half : cy+half, ...
                            cx-half : cx+half);

            % Search window in deformed image
            sy1 = cy - half - MAX_SHIFT;
            sy2 = cy + half + MAX_SHIFT;
            sx1 = cx - half - MAX_SHIFT;
            sx2 = cx + half + MAX_SHIFT;

            % Clamp to image bounds
            sy1 = max(sy1, 1); sy2 = min(sy2, h_roi);
            sx1 = max(sx1, 1); sx2 = min(sx2, w_roi);

            search = def_n(sy1:sy2, sx1:sx2);

            % NCC using normxcorr2
            C = normxcorr2(ref_sub, search);

            % Valid region of C (remove padding)
            C_valid = C(2*half+1 : end-2*half, ...
                        2*half+1 : end-2*half);

            % Find peak
            [max_val, max_idx] = max(C_valid(:));
            [r_peak, c_peak]   = ind2sub(size(C_valid), max_idx);

            % Integer displacement
            int_dx = (sx1 - 1) + c_peak - (cx - half);
            int_dy = (sy1 - 1) + r_peak - (cy - half);

            % Parabolic sub-pixel refinement
            [sub_dx, sub_dy] = parabolic_subpixel(C_valid, r_peak, c_peak);

            disp_x(i,j) = int_dx + sub_dx;
            disp_y(i,j) = int_dy + sub_dy;
            corr(i,j)   = max_val;
        end
    end

    % ── Outlier rejection ────────────────────────────────────
    [disp_x, mask_x] = reject_outliers(disp_x, OUTLIER_NSIGMA);
    [disp_y, mask_y] = reject_outliers(disp_y, OUTLIER_NSIGMA);
    bad_mask = (corr < MIN_CORR) | mask_x | mask_y;

    % ── Grip-referenced displacement ─────────────────────────
    ref_ux = mean(mean(disp_x(:, 1:REF_GRIP_COLS)));
    ref_uy = mean(mean(disp_y(:, 1:REF_GRIP_COLS)));
    disp_x = disp_x - ref_ux;
    disp_y = disp_y - ref_uy;

    % ── Smooth displacements ─────────────────────────────────
    ux = imgaussfilt(disp_x, 1.0);
    uy = imgaussfilt(disp_y, 1.0);

    % Mask bad points
    ux(bad_mask) = NaN;
    uy(bad_mask) = NaN;

    % ── Strain via LSQ window ────────────────────────────────
    [exx, eyy, exy] = compute_strain_lsq(ux, uy, STEP_SIZE, STRAIN_WINDOW);

    % Smooth strain fields
    exx = smooth_field(exx, STRAIN_SMOOTH);
    eyy = smooth_field(eyy, STRAIN_SMOOTH);
    exy = smooth_field(exy, STRAIN_SMOOTH);

    % ── Statistics ───────────────────────────────────────────
    mc      = mean(corr(~bad_mask), 'omitnan');
    m_exx   = mean(exx(:), 'omitnan');
    m_eyy   = mean(eyy(:), 'omitnan');
    p95_exx = prctile(exx(:), 95);
    max_ux  = max(ux(:), [], 'omitnan');
    if abs(m_exx) > 5e-4
        nu = -m_eyy / m_exx;
    else
        nu = NaN;
    end

    fprintf('ZNCC=%.4f  maxU=%+.2fpx  exx=%+.3f%%', ...
            mc, max_ux, m_exx*100);
    if ~isnan(nu)
        fprintf('  nu=%.3f', nu);
    end
    fprintf('\n');

    % Store results
    rec_frame(frameNo)    = frameNo;
    rec_corr(frameNo)     = mc;
    rec_mean_exx(frameNo) = m_exx;
    rec_mean_eyy(frameNo) = m_eyy;
    rec_max_exx(frameNo)  = p95_exx;
    rec_max_ux(frameNo)   = max_ux;
    rec_nu(frameNo)       = nu;

    % ── Von Mises strain ─────────────────────────────────────
    e_vm = sqrt(exx.^2 + eyy.^2 - exx.*eyy + 3*exy.^2);

    % ── Save frame plot ──────────────────────────────────────
    save_frame_plot(frameNo, ref_img, grid_x, grid_y, ...
                    ux, uy, exx, eyy, exy, e_vm, corr, bad_mask, ...
                    ref_ux, ref_uy, mc, RESULTS_FOLDER);

    % ── Save data ────────────────────────────────────────────
    save(fullfile(RESULTS_FOLDER, sprintf('data_%03d.mat', frameNo)), ...
         'grid_x','grid_y','ux','uy','exx','eyy','exy','corr','bad_mask');
end

%% ── SUMMARY PLOTS ───────────────────────────────────────────

fprintf('\nGenerating summary plots...\n');

fn      = rec_frame;
exx_pct = rec_mean_exx * 100;
eyy_pct = rec_mean_eyy * 100;
p95_pct = rec_max_exx  * 100;

% Valid Poisson frames
valid = ~isnan(rec_nu) & rec_nu > -0.1 & rec_nu < 0.9;
med_nu = median(rec_nu(valid), 'omitnan');

fig = figure('Units','pixels','Position',[0 0 1400 900], ...
             'Color','white','Visible','off');
set(fig,'InvertHardcopy','off');
sgtitle('DIC Summary — GRISHMA Project, IIT Kharagpur', ...
        'FontSize',14, 'FontWeight','bold', 'Color','black');

% ── Axial strain ──────────────────────────────────────────
ax1 = subplot(2,2,1);
set(ax1,'Color','white','XColor','black','YColor','black', ...
        'GridColor',[0.8 0.8 0.8],'MinorGridColor',[0.9 0.9 0.9]);
hold(ax1,'on');
fill(ax1,[fn; flipud(fn)],[exx_pct; zeros(size(exx_pct))], ...
     [0.2 0.4 0.8],'FaceAlpha',0.15,'EdgeColor','none');
plot(ax1,fn,p95_pct,'--','LineWidth',1,'Color',[0.5 0.6 1.0]);
plot(ax1,fn,exx_pct,'b-o','MarkerSize',4,'LineWidth',1.8, ...
     'MarkerFaceColor','b');
xlabel(ax1,'Frame'); ylabel(ax1,'Strain (%)');
title(ax1,'Axial Strain \epsilon_{xx}','Color','black');
legend(ax1,'','95th pct','Mean \epsilon_{xx}','Location','northwest');
grid(ax1,'on'); grid(ax1,'minor');

% ── Transverse strain ─────────────────────────────────────
ax2 = subplot(2,2,2);
set(ax2,'Color','white','XColor','black','YColor','black', ...
        'GridColor',[0.8 0.8 0.8],'MinorGridColor',[0.9 0.9 0.9]);
hold(ax2,'on');
fill(ax2,[fn; flipud(fn)],[eyy_pct; zeros(size(eyy_pct))], ...
     [0.8 0.1 0.1],'FaceAlpha',0.15,'EdgeColor','none');
plot(ax2,fn,eyy_pct,'r-o','MarkerSize',4,'LineWidth',1.8, ...
     'MarkerFaceColor','r');
yline(0,'k--','LineWidth',1.0,'Parent',ax2);
xlabel(ax2,'Frame'); ylabel(ax2,'Strain (%)');
title(ax2,'Transverse Strain \epsilon_{yy}  (negative = Poisson contraction)', ...
     'Color','black');
grid(ax2,'on'); grid(ax2,'minor');

% ── Gauge elongation ──────────────────────────────────────
ax3 = subplot(2,2,3);
set(ax3,'Color','white','XColor','black','YColor','black', ...
        'GridColor',[0.8 0.8 0.8],'MinorGridColor',[0.9 0.9 0.9]);
hold(ax3,'on');
fill(ax3,[fn; flipud(fn)],[rec_max_ux; zeros(size(rec_max_ux))], ...
     [0.1 0.6 0.1],'FaceAlpha',0.15,'EdgeColor','none');
plot(ax3,fn,rec_max_ux,'Color',[0 0.5 0],'Marker','o', ...
     'MarkerSize',4,'LineWidth',1.8,'MarkerFaceColor',[0 0.5 0]);
xlabel(ax3,'Frame'); ylabel(ax3,'Displacement (px)');
title(ax3,'Gauge Elongation (max U relative to left grip)','Color','black');
grid(ax3,'on'); grid(ax3,'minor');
if ~isempty(GAUGE_LENGTH_MM)
    scale = GAUGE_LENGTH_MM / (ROI(3) - ROI(1));
    yyaxis(ax3,'right');
    plot(ax3,fn,rec_max_ux*scale,'g--','LineWidth',1);
    ylabel(ax3,'Elongation (mm)');
end

% ── Poisson ratio ─────────────────────────────────────────
ax4 = subplot(2,2,4);
set(ax4,'Color','white','XColor','black','YColor','black', ...
        'GridColor',[0.8 0.8 0.8],'MinorGridColor',[0.9 0.9 0.9]);
hold(ax4,'on');
scatter(ax4,fn(valid),rec_nu(valid),40,[0.5 0 0.8],'filled');
yline(med_nu,'k--','LineWidth',1.5,'Label', ...
      sprintf('Median \\nu = %.3f',med_nu),'LabelHorizontalAlignment','left');
xlabel(ax4,'Frame'); ylabel(ax4,'Poisson''s ratio \nu');
title(ax4,'Apparent Poisson''s Ratio Evolution','Color','black');
ylim(ax4,[-0.1 0.7]);
grid(ax4,'on'); grid(ax4,'minor');

print(fig, fullfile(RESULTS_FOLDER,'summary.png'),'-dpng','-r150');
close(fig);

%% ── FINAL REPORT ────────────────────────────────────────────

fprintf('\n=======================================================\n');
fprintf('  Median Poisson ratio  : %.3f\n', med_nu);
fprintf('  Peak axial strain     : %.3f%%\n', max(exx_pct));
fprintf('  Max gauge elongation  : %.2f px\n', max(rec_max_ux));
if ~isempty(GAUGE_LENGTH_MM)
    s = GAUGE_LENGTH_MM / (ROI(3)-ROI(1));
    fprintf('  Max elongation (mm)   : %.4f mm\n', max(rec_max_ux)*s);
end
fprintf('  Results saved to      : %s\n', RESULTS_FOLDER);
fprintf('=======================================================\n\n');


%% ════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════

function img = load_gray(fpath)
    raw = imread(fpath);
    if size(raw,3) == 3
        img = double(rgb2gray(raw));
    else
        img = double(raw);
    end
end

function out = crop_roi(img, roi)
    x1=roi(1); y1=roi(2); x2=roi(3); y2=roi(4);
    out = img(y1:y2, x1:x2);
end

function out = normalize_img(img)
    m = mean(img(:));
    s = std(img(:));
    out = (img - m) / (s + 1e-10);
end

function [sub_dx, sub_dy] = parabolic_subpixel(C, r, c)
    [nr, nc] = size(C);
    sub_dx = 0; sub_dy = 0;
    % X direction
    if c > 1 && c < nc
        fl = C(r, c-1); fc = C(r, c); fr = C(r, c+1);
        d  = 2*(fl - 2*fc + fr);
        if abs(d) > 1e-10
            sub_dx = (fl - fr) / d;
            sub_dx = max(-1, min(1, sub_dx));
        end
    end
    % Y direction
    if r > 1 && r < nr
        fu = C(r-1, c); fc = C(r, c); fd = C(r+1, c);
        d  = 2*(fu - 2*fc + fd);
        if abs(d) > 1e-10
            sub_dy = (fu - fd) / d;
            sub_dy = max(-1, min(1, sub_dy));
        end
    end
end

function [cleaned, mask] = reject_outliers(field, nsigma)
    med  = median(field(:), 'omitnan');
    mad  = median(abs(field(:) - med), 'omitnan') * 1.4826 + 1e-10;
    mask = abs(field - med) > nsigma * mad;
    cleaned = field;
    % Replace outliers with local median
    med_field = medfilt2(field, [5 5], 'symmetric');
    cleaned(mask) = med_field(mask);
    med_cleaned = medfilt2(cleaned, [5 5], 'symmetric');
    cleaned(mask) = med_cleaned(mask);
end

function [exx, eyy, exy] = compute_strain_lsq(ux, uy, step_size, window)
    [ny, nx] = size(ux);
    exx = zeros(ny, nx);
    eyy = zeros(ny, nx);
    exy = zeros(ny, nx);
    hw  = floor(window/2);

    for i = 1:ny
        for j = 1:nx
            i0 = max(1, i-hw); i1 = min(ny, i+hw);
            j0 = max(1, j-hw); j1 = min(nx, j+hw);

            ux_p = ux(i0:i1, j0:j1);
            uy_p = uy(i0:i1, j0:j1);

            % Local coordinates
            [JJ, II] = meshgrid(j0-j:j1-j, i0-i:i1-i);
            II_f = II(:); JJ_f = JJ(:);
            ux_f = ux_p(:); uy_f = uy_p(:);

            % Remove NaN points
            valid = ~isnan(ux_f) & ~isnan(uy_f);
            if sum(valid) < 4, continue; end

            A = [ones(sum(valid),1), JJ_f(valid), II_f(valid)];
            % Skip if rank deficient (e.g. all points in same row/col)
            if rank(A) < 3, continue; end
            cx_u = A \ ux_f(valid);
            cx_v = A \ uy_f(valid);

            exx(i,j) = cx_u(2) / step_size;
            eyy(i,j) = cx_v(3) / step_size;
            exy(i,j) = 0.5 * (cx_u(3) + cx_v(2)) / step_size;
        end
    end
end

function out = smooth_field(f, sigma)
    % Gaussian smoothing ignoring NaNs
    nan_mask = isnan(f);
    f(nan_mask) = 0;
    out = imgaussfilt(f, sigma);
    out(nan_mask) = NaN;
end

function save_frame_plot(frameNo, ref_img, gx, gy, ...
                          ux, uy, exx, eyy, exy, e_vm, ...
                          corr, bad_mask, rx, ry, mc, outdir)

    % Mask arrays
    ux_plot   = ux;   ux_plot(bad_mask)  = NaN;
    uy_plot   = uy;   uy_plot(bad_mask)  = NaN;
    exx_plot  = exx;  exx_plot(bad_mask) = NaN;
    eyy_plot  = eyy;  eyy_plot(bad_mask) = NaN;
    exy_plot  = exy;  exy_plot(bad_mask) = NaN;
    evm_plot  = e_vm; evm_plot(bad_mask) = NaN;
    corr_plot = corr; corr_plot(bad_mask)= NaN;

    % Plot each panel separately and save as one combined image
    % Use tall figure with 4 rows x 2 cols for better aspect ratio
    fig = figure('Units','pixels','Position',[0 0 1600 1800], ...
                 'Color','white','Visible','off');

    ttl = sprintf('DIC Results — Frame %03d  |  ZNCC=%.4f  |  Left-grip ref: (%.2f, %.2f) px', ...
                  frameNo, mc, rx, ry);
    annotation('textbox',[0 0.97 1 0.03],'String',ttl, ...
               'FontSize',12,'FontWeight','bold', ...
               'HorizontalAlignment','center','EdgeColor','none');

    % Subplot positions: [left bottom width height]
    % 4 rows, 2 cols with good spacing
    W = 0.44; H = 0.20;
    L1 = 0.04; L2 = 0.52;
    B4 = 0.76; B3 = 0.52; B2 = 0.28; B1 = 0.04;

    % Row 4: Reference image + Correlation
    ax = axes('Position',[L1 B4 W H]);
    imagesc(ax, ref_img); colormap(ax,'gray'); axis(ax,'image','off');
    hold(ax,'on');
    scatter(ax, gx(~bad_mask), gy(~bad_mask), 2, 'g', 'filled');
    if any(bad_mask(:))
        scatter(ax, gx(bad_mask), gy(bad_mask), 2, 'r', 'filled');
    end
    title(ax, sprintf('ROI + Grid  |  valid=%d  masked=%d', ...
          sum(~bad_mask(:)), sum(bad_mask(:))), 'FontSize',10);

    ax = axes('Position',[L2 B4 W H]);
    pcolor_plot(gx, gy, corr_plot, 'ZNCC Correlation', ...
                'rdylgn', false, 0.7, 1.0);

    % Row 3: U and V displacement
    ax = axes('Position',[L1 B3 W H]);  %#ok
    pcolor_plot(gx, gy, ux_plot, 'U — Horizontal displacement (px)  [left grip = 0]', ...
                'rdbu_r', true, [], []);

    ax = axes('Position',[L2 B3 W H]);  %#ok
    pcolor_plot(gx, gy, uy_plot, 'V — Vertical displacement (px)  [left grip = 0]', ...
                'rdbu_r', true, [], []);

    % Row 2: exx and eyy
    ax = axes('Position',[L1 B2 W H]);  %#ok
    pcolor_plot(gx, gy, exx_plot, '\epsilon_{xx} — Axial Strain', ...
                'jet', true, [], []);

    ax = axes('Position',[L2 B2 W H]);  %#ok
    pcolor_plot(gx, gy, eyy_plot, '\epsilon_{yy} — Transverse Strain', ...
                'jet', true, [], []);

    % Row 1: Von Mises and exy
    ax = axes('Position',[L1 B1 W H]);  %#ok
    pcolor_plot(gx, gy, evm_plot, 'Von Mises Strain', ...
                'hot', false, 0, []);

    ax = axes('Position',[L2 B1 W H]);  %#ok
    pcolor_plot(gx, gy, exy_plot, '\epsilon_{xy} — Shear Strain', ...
                'jet', true, [], []);

    % Save at high resolution
    print(fig, fullfile(outdir, sprintf('frame_%03d.png', frameNo)), ...
          '-dpng', '-r150');
    close(fig);
end

function pcolor_plot(gx, gy, data, ttl, cmap_name, symmetric, vmin, vmax)
    s = pcolor(gx, gy, data);
    s.EdgeColor = 'none';
    shading interp;
    % Use MATLAB built-in colormaps
    switch lower(cmap_name)
        case 'rdbu_r'
            cm = make_rdbu(256);
        case 'rdylgn'
            cm = make_rdylgn(256);
        case 'jet'
            cm = jet(256);
        case 'hot'
            cm = hot(256);
        otherwise
            cm = jet(256);
    end
    colormap(gca, cm);
    if symmetric
        lim = max(abs(data(:)), [], 'omitnan');
        lim = max(lim, 1e-8);
        clim([-lim, lim]);
    else
        if ~isempty(vmin) && ~isempty(vmax)
            clim([vmin, vmax]);
        else
            hi = max(data(:), [], 'omitnan');
            if isnan(hi), hi = 1; end
            clim([0, hi + 1e-8]);
        end
    end
    colorbar; axis image;
    set(gca,'YDir','reverse');
    title(ttl, 'FontSize',9);
end

function cm = make_rdbu(n)
    % Red-White-Blue diverging colormap
    top    = [0.7 0.0 0.1];
    mid    = [1.0 1.0 1.0];
    bot    = [0.1 0.3 0.8];
    h1     = round(n/2);
    h2     = n - h1;
    cm1    = [linspace(bot(1),mid(1),h1)', ...
              linspace(bot(2),mid(2),h1)', ...
              linspace(bot(3),mid(3),h1)'];
    cm2    = [linspace(mid(1),top(1),h2)', ...
              linspace(mid(2),top(2),h2)', ...
              linspace(mid(3),top(3),h2)'];
    cm     = [cm1; cm2];
end

function cm = make_rdylgn(n)
    % Red-Yellow-Green colormap for correlation
    colors = [0.7 0.0 0.0;   % red
              1.0 0.5 0.0;   % orange
              1.0 1.0 0.0;   % yellow
              0.5 0.8 0.2;   % light green
              0.0 0.5 0.0];  % green
    cm = interp1(linspace(0,1,size(colors,1)), colors, ...
                 linspace(0,1,n));
end

function sorted = sort_nat(names)
    % Natural sort for filenames (img1, img2, ..., img10)
    n = length(names);
    nums = zeros(n,1);
    for i = 1:n
        tok = regexp(names{i}, '\d+', 'match');
        if ~isempty(tok)
            nums(i) = str2double(tok{end});
        end
    end
    [~, idx] = sort(nums);
    sorted = names(idx);
end