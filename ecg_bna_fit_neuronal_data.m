function output = ecg_bna_fit_neuronal_data(cfg, phase_bin_centers, spike_phases_histogram2, n_valid_RRinterval, model_type, scaling_factor_sign)
    % Define common inputs
    if numel(spike_phases_histogram2) == length(phase_bin_centers)
        x         = phase_bin_centers;
        currCurve = spike_phases_histogram2';
    else
        x = repmat(phase_bin_centers, 1, n_valid_RRinterval);
        currCurve = nanmean(spike_phases_histogram2, 2);
    end
    
    % Set up initial parameters and perform fitting based on model type
    switch model_type
        case 'linear'
            % Linear model: no specific initial parameters required
            output = perform_fitting('linear', x, spike_phases_histogram2(:), phase_bin_centers, currCurve, cfg);
            
        case 'cosine'
            % Cosine model: calculate initial parameters
            a1 = (nanmax(currCurve) - nanmin(currCurve)) / 2;
            b1 = mod(circ_mean(phase_bin_centers', currCurve), 2 * pi);
            c1 = nanmean(currCurve);
            startPoint = [a1, b1, c1];
            output = perform_fitting('cosine', x', spike_phases_histogram2(:), phase_bin_centers, currCurve, cfg, startPoint);
            
        case 'vonMises'
            % Von Mises model: calculate initial parameters
            a1 = nanmax(currCurve) - nanmin(currCurve);
            
            if scaling_factor_sign > 0
                k1 = circ_kappa(phase_bin_centers, (currCurve - nanmin(currCurve)) / (nanmax(currCurve) - nanmin(currCurve)));
                t1 = mod(circ_mean(phase_bin_centers', currCurve), 2 * pi);
            else
                currCurve_tmp = -1 * currCurve;
                k1 = circ_kappa(phase_bin_centers, (currCurve_tmp - nanmin(currCurve_tmp)) / (nanmax(currCurve_tmp) - nanmin(currCurve_tmp)));
                t1 = mod(circ_mean(phase_bin_centers', currCurve_tmp), 2 * pi);
            end
            
            d1 = nanmean(currCurve);
            
            if isnan(k1)
                k1 = 0;
            end
            
            startPoint = [a1, d1, k1, t1];
            output = perform_fitting('vonMises', x', spike_phases_histogram2(:), phase_bin_centers, currCurve, cfg, startPoint, scaling_factor_sign);
    end
    
    % Output setup
    output.average = currCurve;
    
end

% Define a function to handle fitting and calculating metrics
function [output, fittedmdl, lin_mdl] = perform_fitting(model_type, x, y, phase_bin_centers, currCurve, cfg, startPoint, scaling_factor_sign)

% Initialize output
output = struct();

switch model_type
    case 'linear'
        % Linear fit
        fittedmdl = fitlm(x, y);
        output.yfit        = fittedmdl.Coefficients.Estimate(2) * phase_bin_centers + fittedmdl.Coefficients.Estimate(1);
        output.coefs       = fittedmdl.Coefficients.Estimate;
        output.rsquared    = fittedmdl.Rsquared.Ordinary;
        output.adjrsquared = fittedmdl.Rsquared.Adjusted;
        output.pvalue      = fittedmdl.Coefficients.pValue;
        
    case 'cosine'
        % Cosine fit
        % drop nans
        [y,dropnans] = rmmissing(y);
        x            = x(~dropnans);
        
        % do the non-linear fit
        [fittedmdl,gof] = fit(double(x),y,cfg.fit.cos_mod,'StartPoint', startPoint, 'Lower', cfg.fit.cos_lower, 'Upper', cfg.fit.cos_upper);
        
        coefs = coeffvalues(fittedmdl); % get model coefficients
        coefs(2) = mod(coefs(2),2*pi);
        
        yfit_all = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), x);
        yfit     = cfg.fit.cos_mod(coefs(1), coefs(2), coefs(3), phase_bin_centers);
        
        % employ a linear fit to get a p-value vs. fitting with a
        % constant model
        lin_mdl = fitlm(y(:), yfit_all);
        
        % put data into the 'output' structure
        output.yfit = yfit;
        output.coefs = coefs;
        output.rsquared = gof.rsquare;
        output.adjrsquared = gof.adjrsquare;
        output.sse = gof.sse;
        output.dfe = gof.dfe;
        output.rmse = gof.rmse;
        output.CI = confint(fittedmdl);
        output.startPoint = startPoint;
        output.pvalue = lin_mdl.Coefficients.pValue(2);
    case 'vonMises'
        % Von Mises fit
        
        if scaling_factor_sign > 0
            lower_bounds = cfg.fit.vMpos_lower;
            upper_bounds = cfg.fit.vMpos_upper;
        else
            lower_bounds = cfg.fit.vMneg_lower;
            upper_bounds = cfg.fit.vMneg_upper;
        end
        
        % do the non-linear fit
        [fittedmdl,gof] = fit(x,y,cfg.fit.vonMises_mod,'StartPoint', startPoint, 'Lower', lower_bounds, 'Upper', upper_bounds);
        
        coefs    = coeffvalues(fittedmdl); % get model coefficients
        coefs(4) = mod(coefs(4),2*pi);
        
        yfit_all = cfg.fit.vonMises_mod(coefs(1), coefs(2), coefs(3), coefs(4), x);
        yfit     = cfg.fit.vonMises_mod(coefs(1), coefs(2), coefs(3), coefs(4), phase_bin_centers);
        
        % employ a linear fit to get a p-value vs. fitting with a
        % constant model
        lin_mdl = fitlm(y, yfit_all);
        
        output.yfit = yfit;
        output.coefs = coefs;
        output.rsquared = gof.rsquare;
        output.adjrsquared = gof.adjrsquare;
        output.sse = gof.sse;
        output.dfe = gof.dfe;
        output.rmse = gof.rmse;
        output.CI = confint(fittedmdl);
        output.startPoint = startPoint;
        output.pvalue = lin_mdl.Coefficients.pValue(2);
end

% Calculate AIC and BIC
[output.aic, output.bic] = calculateAICBIC(phase_bin_centers, fittedmdl, currCurve);

end

% Helper functions
function [aic, bic] = calculateAICBIC(phase_bin_centers, mdl, currCurve)
    % Calculate residuals
    ydata = feval(mdl, phase_bin_centers);
    residuals = currCurve(:) - ydata(:);
    
    % Calculate the sum of squared residuals (RSS)
    rss = sum(residuals.^2);
    
    % Number of data points
    n = length(phase_bin_centers);
    
    % Calculate AIC
    try
        aic = n * log(rss / n) + 2 * mdl.NumEstimatedCoefficients;
    end
    try
        aic = n * log(rss / n) + 2 * nargin(mdl);
    end
    
    % Calculate BIC
    try
        bic = n * log(rss / n) + mdl.NumEstimatedCoefficients * log(n);
    end
    try
        bic = n * log(rss / n) + nargin(mdl) * log(n);
    end
end
