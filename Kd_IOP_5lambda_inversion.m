function [chl, aphi, adg, a, bbp, bb, MSE] = ...
    Kd_IOP_5lambda_inversion(Kd, chl_initial, plotlogical, S_CDOM, mu)

% ========================================================================
% Invert spectral Kd for:
%   chl
%   a_phi(lambda)
%   a_dg(lambda)
%   bbp(lambda)
%
% Wavelengths:
%   412, 443, 490, 510, 555 nm
%
% Iterative scheme:
%
%   1. Evaluate Bricaud a_phi*(lambda) at current chl
%   2. Evaluate MM01 bbp(lambda) at current chl
%   3. Solve linear least squares for [chl, adg_443]
%   4. Recompute FULL model using candidate chl
%   5. Compute actual model MSE
%   6. Accept candidate only if actual MSE decreases
%   7. Otherwise retain previous accepted solution and stop
%
% Requirements:
%
%   Kd at lambda = 412, 443, 490, 510, 555 nm
%   chl_initial initial chl concentration guess required (~0.5 mg/m3)
%   plotlogical = 1 plots the results of the inversion (0) 
%   S_CDOM custom or 0.012-0.016 typical in oceanic waters
%   mu 0.76 (~0.7-0.9)
% 
%
% version 2026/10/16 
%
% Terence O'Brien
% ========================================================================


%% ------------------------------------------------------------------------
% Constants
% -------------------------------------------------------------------------

lambda = [412, 443, 490, 510, 555];

% Force Kd to same orientation as lambda
Kd = Kd(:).';

Temp = 15;
Salt = 33;

[aw, bbw] = IOP_of_sw(lambda, Temp, Salt);

aw  = aw(:).';
bbw = bbw(:).';


%% ------------------------------------------------------------------------
% Model parameters
% -------------------------------------------------------------------------

S = S_CDOM;

lambda_0 = 443;

maxIter = 20;


%% ------------------------------------------------------------------------
% Observed RHS
%
%   mu Kd = aw + aphi + adg + bbw + bbp
%
% therefore:
%
%   h = mu Kd - aw - bbw
%     = aphi + adg + bbp
% -------------------------------------------------------------------------

h = mu .* Kd - bbw - aw;


%% ------------------------------------------------------------------------
% CDOM/detrital spectral shape
%
% adg(lambda) = adg_443 * a_dg_star(lambda)
% -------------------------------------------------------------------------

a_dg_star = exp(-S .* (lambda - lambda_0));


%% ------------------------------------------------------------------------
% Initialize iteration
% -------------------------------------------------------------------------

chl_old = chl_initial;

MSE_old = Inf;

best = struct();

best.chl      = NaN;
best.adg_443  = NaN;
best.aphi     = nan(size(lambda));
best.adg      = nan(size(lambda));
best.bbp      = nan(size(lambda));
best.bb       = nan(size(lambda));
best.MSE      = Inf;


%% ========================================================================
% ITERATIVE SOLUTION
% ========================================================================

for it = 1:maxIter


    %% --------------------------------------------------------------------
    % Bricaud a_phi* evaluated at CURRENT chl
    % ---------------------------------------------------------------------



    [~, ~, ~, a_phi_star] = ...
       get_bricaud_aph(chl_old, lambda);

    a_phi_star = a_phi_star(:);


    %% --------------------------------------------------------------------
    % MM01 bbp evaluated at CURRENT chl
    % ---------------------------------------------------------------------

    [bbp_current, ~] = ...
        MM01bbpchl(lambda, chl_old);

    bbp_current = bbp_current(:);


    %% --------------------------------------------------------------------
    % Remove bbp contribution before solving absorption components
    %
    % h2 = aphi + adg
    % ---------------------------------------------------------------------

    h2 = h(:) - bbp_current;


    %% --------------------------------------------------------------------
    % Linear system evaluated around CURRENT chl
    %
    % aphi(lambda) = chl * a_phi_star(lambda)
    %
    % adg(lambda)  = adg_443 * a_dg_star(lambda)
    %
    % h2 = H * [chl ; adg_443]
    % ---------------------------------------------------------------------

    H = [a_phi_star(:), a_dg_star(:)];

    x = H \ h2;

    chl_new     = x(1);
    adg_443_new = x(2);


    %% --------------------------------------------------------------------
    % Candidate QC before evaluating parameterizations
    % ---------------------------------------------------------------------

    if ~isfinite(chl_new) || chl_new <= 0 || ...
       ~isfinite(adg_443_new) || adg_443_new < 0

        break

    end


    %% ====================================================================
    % RECOMPUTE COMPLETE MODEL USING CANDIDATE CHL
    %
    % This is the important correction.
    %
    % The candidate solution is tested using:
    %
    %   aphi(chl_new)
    %   bbp(chl_new)
    %   adg(adg_443_new)
    %
    % rather than using spectral quantities evaluated at chl_old.
    % ====================================================================


    %% --------------------------------------------------------------------
    % Bricaud evaluated at CANDIDATE chl
    % ---------------------------------------------------------------------

    [~, ~, ~, a_phi_star_new] = ...
         get_bricaud_aph(chl_new, lambda);


    a_phi_star_new = a_phi_star_new(:);


    %% --------------------------------------------------------------------
    % Candidate aphi
    % ---------------------------------------------------------------------

    aphi_new = ...
        chl_new .* a_phi_star_new;


    %% --------------------------------------------------------------------
    % Candidate adg
    % ---------------------------------------------------------------------

    adg_new = ...
        adg_443_new .* a_dg_star(:);


    %% --------------------------------------------------------------------
    % MM01 evaluated at CANDIDATE chl
    % ---------------------------------------------------------------------

    [bbp_new, bb_new] = ...
        MM01bbpchl(lambda, chl_new);

    bbp_new = bbp_new(:);
    bb_new  = bb_new(:);


    %% --------------------------------------------------------------------
    % Complete candidate forward model
    %
    % h_model = aphi + adg + bbp
    % ---------------------------------------------------------------------

    h_modeled_new = ...
        aphi_new ...
        + adg_new ...
        + bbp_new;


    residuals_new = ...
        h(:) - h_modeled_new;


    MSE_new = ...
        mean(residuals_new.^2);


    %% --------------------------------------------------------------------
    % Accept ONLY if FULL model MSE improves
    % ---------------------------------------------------------------------
tolRel = 1e-6;

if isinf(MSE_old) || ...
   MSE_new < MSE_old * (1 - tolRel)

   % if isfinite(MSE_new) && MSE_new < MSE_old

        % Store complete accepted state
        best.chl      = chl_new;
        best.adg_443  = adg_443_new;
        best.aphi     = aphi_new;
        best.adg      = adg_new;
        best.bbp      = bbp_new;
        best.bb       = bb_new;
        best.MSE      = MSE_new;

        MSE_old = MSE_new;

        % Next iteration is constructed around accepted chl
        chl_old = chl_new;

    else

        % Candidate made the actual forward-model fit worse.
        %
        % Reject it and retain previous accepted state.
        break

    end

end


%% ========================================================================
% FINAL ACCEPTED SOLUTION
% ========================================================================

chl     = best.chl;
adg_443 = best.adg_443;

aphi = best.aphi;
adg  = best.adg;
bbp  = best.bbp;
bb   = best.bb;

MSE = best.MSE;


%% ------------------------------------------------------------------------
% Total absorption
% -------------------------------------------------------------------------

a = aw(:) + aphi(:) + adg(:);


%% ------------------------------------------------------------------------
% Final modeled h
% -------------------------------------------------------------------------

h_modeled = ...
    aphi(:) ...
    + adg(:) ...
    + bbp(:);


%% ------------------------------------------------------------------------
% QC
% -------------------------------------------------------------------------

if ~isfinite(chl) || chl <= 0 || ...
   any(~isfinite(aphi),'all') || ...
   any(~isfinite(adg),'all')  || ...
   any(~isfinite(bbp),'all')  || ...
   any(aphi < 0,'all') || ...
   any(adg < 0,'all')  || ...
   any(bbp < 0,'all')  || ...
   ~isfinite(MSE)

    chl  = NaN;
    aphi = nan(size(lambda));
    adg  = nan(size(lambda));
    a    = nan(size(lambda));
    bbp  = nan(size(lambda));
    bb   = nan(size(lambda));
    MSE  = NaN;

end


%% ========================================================================
% OPTIONAL PLOTTING
% ========================================================================

if plotlogical == 1 && isfinite(chl)

    figure;
    set(gcf,'Units','inches','Position',[4 4 12 5]);


    % ---------------------------------------------------------------------
    % Observed versus modeled h
    % ---------------------------------------------------------------------

    subplot(1,2,1)

    plot(lambda,h,'o-', ...
        'LineWidth',1.5, ...
        'DisplayName','Observed h(\lambda)');

    hold on

    plot(lambda,h_modeled,'s--', ...
        'LineWidth',1.5, ...
        'DisplayName','Modeled h(\lambda)');

    xlabel('Wavelength (nm)')
    ylabel('h(\lambda) (m^{-1})')

    title('Observed vs Modeled h(\lambda)')

    legend
    grid on

    ax = gca;
    set(ax,'LineWidth',1,'FontSize',16)


    % ---------------------------------------------------------------------
    % Retrieved IOP components
    % ---------------------------------------------------------------------

    subplot(1,2,2)

    plot(lambda,aphi,'-o', ...
        'LineWidth',1.5, ...
        'DisplayName','a_{\phi}(\lambda)');

    hold on

    plot(lambda,adg,'-s', ...
        'LineWidth',1.5, ...
        'DisplayName','a_{dg}(\lambda)');

    plot(lambda,aw,'-^', ...
        'LineWidth',1.5, ...
        'DisplayName','a_w(\lambda)');

    plot(lambda,bbp,'-d', ...
        'LineWidth',1.5, ...
        'DisplayName','b_{bp}(\lambda)');

    xlabel('Wavelength (nm)')
    ylabel('(m^{-1})')

    title('IOP Components')

    legend
    grid on

    ax = gca;
    set(ax,'LineWidth',1,'FontSize',16)

end

end