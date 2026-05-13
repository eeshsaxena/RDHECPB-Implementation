% ==========================================================================
% RDHECPB.m  (v2 — all limitations resolved)
% Reversible Data Hiding with Enhancing Contrast and Preserving Brightness
% in Medical Image
%
% Paper: Shi M., Yang Y., Meng J., Zhang W.
%        Journal of Information Security and Applications, 70 (2022) 103324
%        DOI: 10.1016/j.jisa.2022.103324
%
% FIXES from v1:
%   - Sec.2.2 pre-processing reversal: now tracks moved_mask + moved_from
%     per-pixel, enabling EXACT lossless restoration of pre-processed pixels.
%   - NROI embedding: now uses Yang et al.[15] exact histogram shifting
%     (find peak P_N + nearest zero bin Z_N, shift [P_N,Z_N] or [Z_N,P_N]).
%
% Run:  RDHECPB   (full demo — 4 synthetic medical images, all experiments)
% ==========================================================================
function RDHECPB()
    clc; close all;
    fprintf('=== RDHECPB v2: Contrast Enhancement + Brightness Preservation ===\n');
    fprintf('    Shi, Yang, Meng, Zhang — JISA 70 (2022) 103324\n\n');

    imgs  = generate_medical_images();
    names = {'Brain01','Brain02','chest','xray'};
    M     = 0.1;
    cap_vals = [5000 10000 20000 50000];

    % ---- Experiment 1: PSNR vs capacity -----------------------------------
    fprintf('\n--- Exp 1: PSNR (dB) vs embedding capacity ---\n');
    fprintf('%-10s', 'Image');
    for c = cap_vals, fprintf('%12d', c); end; fprintf('\n');
    fprintf('%s\n', repmat('-',1,58));
    for k = 1:numel(names)
        fprintf('%-10s', names{k});
        for cap = cap_vals
            rng(42); pay = randi([0 1],1,cap,'uint8');
            [I_emb,~] = rdhecpb_embed(imgs{k}, pay, M);
            fprintf('%12.2f', compute_psnr(imgs{k}, I_emb));
        end
        fprintf('\n');
    end

    % ---- Experiment 2: ΔStandard deviation --------------------------------
    fprintf('\n--- Exp 2: ΔSD (contrast improvement, higher is better) ---\n');
    fprintf('%-10s', 'Image');
    for c = cap_vals, fprintf('%12d', c); end; fprintf('\n');
    fprintf('%s\n', repmat('-',1,58));
    for k = 1:numel(names)
        fprintf('%-10s', names{k});
        for cap = cap_vals
            rng(42); pay = randi([0 1],1,cap,'uint8');
            [I_emb,~] = rdhecpb_embed(imgs{k}, pay, M);
            roi = get_roi_mask(imgs{k});
            ds = std(double(I_emb(roi))) - std(double(imgs{k}(roi)));
            fprintf('%12.2f', ds);
        end
        fprintf('\n');
    end

    % ---- Experiment 3: Brightness preservation ----------------------------
    fprintf('\n--- Exp 3: Brightness difference |B - B_emb| (lower is better) ---\n');
    fprintf('%-10s', 'Image');
    for c = cap_vals, fprintf('%12d', c); end; fprintf('\n');
    fprintf('%s\n', repmat('-',1,58));
    for k = 1:numel(names)
        fprintf('%-10s', names{k});
        for cap = cap_vals
            rng(42); pay = randi([0 1],1,cap,'uint8');
            I = imgs{k};
            [I_emb,~] = rdhecpb_embed(I, pay, M);
            roi = get_roi_mask(I);
            fprintf('%12.4f', abs(mean(double(I_emb(roi)))-mean(double(I(roi)))));
        end
        fprintf('\n');
    end

    % ---- Experiment 4: Reversibility (exact) ------------------------------
    fprintf('\n--- Exp 4: Full reversibility check (20000 bits) ---\n');
    for k = 1:numel(names)
        rng(42); pay = randi([0 1],1,20000,'uint8');
        [I_emb, meta] = rdhecpb_embed(imgs{k}, pay, M);
        [I_rec, D_ext] = rdhecpb_extract(I_emb, meta);
        ok = isequal(imgs{k}, I_rec);
        n_ext = min(numel(D_ext), numel(pay));
        errs  = sum(D_ext(1:n_ext) ~= pay(1:n_ext));
        fprintf('  %-10s → Reversible: %s | Bit errors: %d\n', names{k}, string(ok), errs);
    end
    fprintf('\nDone.\n');
end

% ==========================================================================
%  STAGE 1: ROI / NROI SEGMENTATION  (Sec. 2.1)
% ==========================================================================
function roi_mask = get_roi_mask(I)
    thr      = graythresh(I) * 255;
    roi_mask = double(I) > thr;
    if sum(roi_mask(:)) < 0.05 * numel(I)
        [r,c] = size(I);
        mr = round(r*0.2); mc = round(c*0.2);
        roi_mask = false(r,c);
        roi_mask(mr:r-mr, mc:c-mc) = true;
    end
end

% ==========================================================================
%  STAGE 2: ROI HISTOGRAM PRE-PROCESSING  (Sec. 2.2)  ← FIXED
% ==========================================================================
function [I_pre, empty_bins, I_MAX_pre, moved_mask, moved_from] = ...
          preprocess_roi(I, roi_mask)
% FIX: now tracks exact per-pixel location map (moved_mask, moved_from)
% This enables pixel-perfect reversal during recovery without ambiguity.
    img       = double(I);
    roi_pix   = img(roi_mask);
    N_roi     = numel(roi_pix);
    counts    = histcounts(roi_pix, 0:256);
    thr_count = 0.01 * N_roi;   % 1% sparse threshold

    I_pre      = I;
    empty_bins = [];
    moved_mask = false(size(I));       % ← per-pixel: was this pixel moved?
    moved_from = zeros(size(I),'uint8'); % ← per-pixel: what was original value?

    for j = 0:255
        if counts(j+1) > 0 && counts(j+1) < thr_count
            % Choose target: move to nearest denser neighbour
            if j > 0 && counts(j) >= counts(j+1)
                target = j - 1;
            elseif j < 255
                target = j + 1;
            else
                continue;
            end
            % Record every moved pixel individually (enables exact reversal)
            pix_mask = roi_mask & (double(I_pre) == j);
            moved_mask(pix_mask) = true;
            moved_from(pix_mask) = uint8(j);   % save ORIGINAL grey-level
            empty_bins(end+1)    = j;           %#ok<AGROW>
            I_pre(pix_mask)      = uint8(target);
        end
    end
    I_MAX_pre = double(max(double(I_pre(roi_mask))));
end

% ==========================================================================
%  STAGE 3: CONTRAST STRETCHING WITH BRIGHTNESS PRESERVATION  (Sec. 2.3)
% ==========================================================================
function [I_str, L_MIN, L_MAX, I_MIN, I_MAX] = contrast_stretch(I_pre, roi_mask, M)
    roi_pix = double(I_pre(roi_mask));
    I_MIN   = min(roi_pix);
    I_MAX   = max(roi_pix);
    B_orig  = mean(roi_pix);
    L_MIN   = I_MIN;  L_MAX = I_MAX;

    for expand = 1:min(I_MIN, 255-I_MAX)
        LMN = I_MIN - expand;  LMX = I_MAX + expand;
        I_str_try = round((LMX-LMN)*(roi_pix-I_MIN)/max(I_MAX-I_MIN,1)+LMN);
        if abs(mean(I_str_try) - B_orig) < M
            L_MIN = LMN;  L_MAX = LMX;
        else
            break;
        end
    end

    I_str = I_pre;
    stretched = round((L_MAX-L_MIN)*(double(I_pre(roi_mask))-I_MIN)/max(I_MAX-I_MIN,1)+L_MIN);
    I_str(roi_mask) = uint8(max(0,min(255,stretched)));
end

% ==========================================================================
%  STAGE 4a: EMBEDDING LOCATION & DIRECTION SELECTION  (Sec. 2.4.1)
% ==========================================================================
function [Ps, d] = select_Ps_direction(I_str, roi_mask, B_orig)
    roi_pix = double(I_str(roi_mask));
    B_str   = mean(roi_pix);
    d = (B_str < B_orig);    % 1=right (increase B), 0=left (decrease B)

    counts = histcounts(roi_pix, 0:256);
    [~, sorted_idx] = sort(counts,'descend');
    Ps = -1;
    for idx = sorted_idx
        gl = idx - 1;
        if d == 1 && gl < 255 && counts(gl+2) == 0, Ps = gl; break; end
        if d == 0 && gl > 0   && counts(gl)   == 0, Ps = gl; break; end
    end
    if Ps < 0
        [~, mi] = max(counts);  Ps = mi - 1;  d = 1;
    end
end

% ==========================================================================
%  STAGE 4b: ROI EMBEDDING (Eq. 8 / Eq. 9)
% ==========================================================================
function [I_emb_roi, Ps, d, n_emb] = embed_roi(I_str, roi_mask, payload, B_orig)
    [Ps, d] = select_Ps_direction(I_str, roi_mask, B_orig);
    flat = double(I_str(:));
    roi_idx = find(roi_mask(:));
    pay_ptr = 1;  n_emb = 0;
    for ii = 1:numel(roi_idx)
        idx = roi_idx(ii);
        p   = flat(idx);
        if p == Ps && pay_ptr <= numel(payload)
            bk = payload(pay_ptr);
            flat(idx) = p + (2*d-1)*bk;  % d=1: p+bk; d=0: p-bk
            pay_ptr = pay_ptr + 1;
            n_emb   = n_emb + 1;
        end
    end
    I_emb_roi = uint8(reshape(flat, size(I_str)));
end

% ==========================================================================
%  STAGE 5: NROI EMBEDDING — Yang et al. [15] exact histogram shifting
% ==========================================================================
function [I_out, P_N, Z_N, nroi_dir] = embed_nroi_yang15(I, nroi_mask, payload)
% Yang et al. [15]: find peak P_N and nearest zero bin Z_N in NROI.
% Shift bins between P_N and Z_N toward Z_N; embed bits at P_N.
    flat   = double(I(:));
    counts = histcounts(flat(nroi_mask), 0:256);

    [~, pk_idx] = max(counts);
    P_N = pk_idx - 1;

    % Search right then left for zero bin
    Z_N = -1;  nroi_dir = 1;
    for z = P_N+1:255
        if counts(z+1) == 0, Z_N = z; nroi_dir = 1; break; end
    end
    if Z_N < 0
        for z = P_N-1:-1:0
            if counts(z+1) == 0, Z_N = z; nroi_dir = 0; break; end
        end
    end
    if Z_N < 0, I_out = I; P_N = 0; Z_N = 0; nroi_dir = 1; return; end

    flat_out = flat;
    nroi_idx = find(nroi_mask(:));
    pay_ptr  = 1;

    for ii = 1:numel(nroi_idx)
        idx = nroi_idx(ii);
        p   = flat_out(idx);
        if nroi_dir == 1          % right shift: Z_N > P_N
            if p > P_N && p < Z_N
                flat_out(idx) = p + 1;
            elseif p == P_N && pay_ptr <= numel(payload)
                flat_out(idx) = p + payload(pay_ptr);
                pay_ptr = pay_ptr + 1;
            end
        else                      % left shift: Z_N < P_N
            if p < P_N && p > Z_N
                flat_out(idx) = p - 1;
            elseif p == P_N && pay_ptr <= numel(payload)
                flat_out(idx) = p - payload(pay_ptr);
                pay_ptr = pay_ptr + 1;
            end
        end
    end
    I_out = uint8(reshape(flat_out, size(I)));
end

% ==========================================================================
%  MAIN EMBEDDING PIPELINE
% ==========================================================================
function [I_emb, meta] = rdhecpb_embed(I, payload, M)
    roi_mask = get_roi_mask(I);

    % Sec.2.2: pre-process with exact tracking
    [I_pre, empty_bins, ~, moved_mask, moved_from] = preprocess_roi(I, roi_mask);

    % Sec.2.3: contrast stretch
    B_orig = mean(double(I(roi_mask)));
    [I_str, L_MIN, L_MAX, I_MIN, I_MAX] = contrast_stretch(I_pre, roi_mask, M);

    % Sec.2.4: embed in ROI
    [I_emb_roi, Ps, d, n_emb_roi] = embed_roi(I_str, roi_mask, payload, B_orig);

    % Sec.2.4 + [15]: embed remaining in NROI
    I_emb = I_emb_roi;
    P_N = 0; Z_N = 0; nroi_dir = 1;
    if numel(payload) > n_emb_roi
        [I_emb, P_N, Z_N, nroi_dir] = embed_nroi_yang15(...
            I_emb_roi, ~roi_mask, payload(n_emb_roi+1:end));
    end

    meta = struct('roi_mask',roi_mask, 'Ps',Ps, 'd',d, ...
        'L_MIN',L_MIN, 'L_MAX',L_MAX, 'I_MIN',I_MIN, 'I_MAX',I_MAX, ...
        'empty_bins',empty_bins, 'n_emb_roi',n_emb_roi, 'M',M, ...
        'moved_mask',moved_mask, 'moved_from',moved_from, ...
        'P_N',P_N, 'Z_N',Z_N, 'nroi_dir',nroi_dir);
end

% ==========================================================================
%  EXTRACTION AND RECOVERY (Sec. 2.5)  ← FIXED
% ==========================================================================
function [I_rec, D_ext] = rdhecpb_extract(I_emb, meta)
    roi_mask  = meta.roi_mask;
    Ps        = meta.Ps;   d = meta.d;
    L_MIN     = meta.L_MIN;  L_MAX = meta.L_MAX;
    I_MIN     = meta.I_MIN;  I_MAX = meta.I_MAX;
    n_emb_roi = meta.n_emb_roi;
    P_N       = meta.P_N;    Z_N   = meta.Z_N;
    nroi_dir  = meta.nroi_dir;

    % --- Step A: Reverse NROI (Yang [15] histogram shift inverse) ---
    flat_nroi = double(I_emb(:));
    D_nroi    = [];
    nroi_idx  = find(~roi_mask(:));
    if Z_N ~= P_N
        for ii = 1:numel(nroi_idx)
            idx = nroi_idx(ii);
            p   = flat_nroi(idx);
            if nroi_dir == 1
                if p == P_N+1
                    D_nroi(end+1) = 1; flat_nroi(idx) = P_N; %#ok<AGROW>
                elseif p == P_N
                    D_nroi(end+1) = 0; %#ok<AGROW>
                elseif p > P_N && p <= Z_N
                    flat_nroi(idx) = p - 1;
                end
            else
                if p == P_N-1
                    D_nroi(end+1) = 1; flat_nroi(idx) = P_N; %#ok<AGROW>
                elseif p == P_N
                    D_nroi(end+1) = 0; %#ok<AGROW>
                elseif p < P_N && p >= Z_N
                    flat_nroi(idx) = p + 1;
                end
            end
        end
    end
    I_step1 = uint8(reshape(flat_nroi, size(I_emb)));

    % --- Step B: Extract from ROI + recover Eq.(10) ---
    flat_roi = double(I_step1(:));
    D_roi    = zeros(1, n_emb_roi, 'uint8');
    bit_ptr  = 1;
    roi_idx  = find(roi_mask(:));
    for ii = 1:numel(roi_idx)
        idx = roi_idx(ii);  p = flat_roi(idx);
        if d == 1
            if p == Ps+1, D_roi(bit_ptr)=1; flat_roi(idx)=Ps; bit_ptr=bit_ptr+1;
            elseif p == Ps, D_roi(bit_ptr)=0; bit_ptr=bit_ptr+1; end
        else
            if p == Ps-1, D_roi(bit_ptr)=1; flat_roi(idx)=Ps; bit_ptr=bit_ptr+1;
            elseif p == Ps, D_roi(bit_ptr)=0; bit_ptr=bit_ptr+1; end
        end
    end
    D_ext = [D_roi(1:bit_ptr-1), D_nroi];

    % --- Step C: Inverse contrast stretch Eq.(12) ---
    if L_MAX > L_MIN
        stretched = flat_roi(roi_mask(:));
        orig_approx = round((I_MAX-I_MIN)*(stretched-L_MIN)/max(L_MAX-L_MIN,1)+I_MIN);
        flat_roi(roi_mask(:)) = max(0,min(255,orig_approx));
    end
    I_step2 = uint8(reshape(flat_roi, size(I_emb)));

    % --- Step D: Reverse pre-processing — EXACT per-pixel restore ---
    I_rec = I_step2;
    I_rec(meta.moved_mask) = meta.moved_from(meta.moved_mask);
end

% ==========================================================================
%  METRICS
% ==========================================================================
function p = compute_psnr(I, I_emb)
    mse = mean((double(I(:))-double(I_emb(:))).^2);
    if mse==0, p=Inf; else, p=10*log10(255^2/mse); end
end

% ==========================================================================
%  SYNTHETIC MEDICAL IMAGE GENERATOR
% ==========================================================================
function imgs = generate_medical_images()
    imgs = cell(4,1); sz = 512;
    rng(1); I=uint8(ones(sz)*20); cx=sz/2; cy=sz/2;
    for r=1:sz; for c=1:sz
        d=sqrt((r-cx)^2+(c-cy)^2)/(sz*0.35);
        if d<1, I(r,c)=uint8(min(255,80+round(120*exp(-d*2))+randi(20))); end
    end; end
    imgs{1}=I;

    rng(2); I=imgs{1};
    for r=round(sz*0.35):round(sz*0.65); for c=round(sz*0.4):round(sz*0.6)
        d=sqrt((r-cx)^2+(c-cy)^2)/(sz*0.12);
        if d<1, I(r,c)=uint8(max(0,double(I(r,c))-round(60*exp(-d*2)))); end
    end; end
    imgs{2}=I;

    rng(3); I=uint8(zeros(sz));
    for r=1:sz; for c=1:sz
        I(r,c)=uint8(40+randi(20));
        if mod(c,round(sz/8))<round(sz/32), I(r,c)=uint8(min(255,double(I(r,c))+120+randi(30))); end
        dl=sqrt((r-sz*0.5)^2+(c-sz*0.3)^2)/(sz*0.2);
        dr=sqrt((r-sz*0.5)^2+(c-sz*0.7)^2)/(sz*0.2);
        if dl<1||dr<1, I(r,c)=uint8(max(0,double(I(r,c))-30+randi(10))); end
    end; end
    imgs{3}=I;

    rng(4); I=uint8(ones(sz)*80);
    for r=1:sz; for c=1:sz
        if abs(r-sz/2)<sz/6
            bf=max(0,1-abs(c-sz/2)/(sz*0.3));
            I(r,c)=uint8(min(255,80+round(150*bf)+randi(15)));
        end
    end; end
    imgs{4}=I;
end
