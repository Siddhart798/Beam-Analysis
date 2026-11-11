function beam_analysis_smart()
    % BEAM_ANALYSIS_SMART: Analyzes beam parameters with multiple point loads & UDLs, retaining previous entries.
    close all;

    % Persistent variables to store historical user data
    persistent L_old P_old x_p_old w_old x_u1_old x_u2_old E_old I_old

    % Initialize default configurations for the first runtime
    if isempty(L_old), L_old = 10; end
    if isempty(P_old), P_old = [-80]; end
    if isempty(x_p_old), x_p_old = 8; end 
    if isempty(w_old), w_old = []; end      
    if isempty(x_u1_old), x_u1_old = []; end
    if isempty(x_u2_old), x_u2_old = []; end
    if isempty(E_old), E_old = 200; end
    if isempty(I_old), I_old = 8.33e-6; end

    fprintf('\n=== BEAM ANALYSIS INPUTS (Press ENTER to keep Current value) ===\n');

    % 1. Beam Length
    temp = input(sprintf('Enter Beam Length in meters (Current: %.1f m): ', L_old));
    if ~isempty(temp), L_old = temp; end

    % 2. Multiple Point Loads (Interactive Setup)
    fprintf('\n--- Point Loads Configuration ---\n');
    fprintf('Current Point Loads: %s kN at positions: %s m\n', mat2str(P_old), mat2str(x_p_old));
    n_loads = input('Enter number of point loads (Press ENTER to keep current): ');
    
    if isempty(n_loads)
        % Keep existing P_old and x_p_old
    elseif n_loads == 0
        P_old = [];
        x_p_old = [];
    else
        new_P = zeros(1, n_loads);
        new_x = zeros(1, n_loads);
        for i = 1:n_loads
            def_P = 0; if i <= length(P_old), def_P = P_old(i); end
            def_x = 0; if i <= length(x_p_old), def_x = x_p_old(i); end
            
            p_val = input(sprintf('  -> Point Load %d magnitude in kN (Current: %.1f): ', i, def_P));
            if ~isempty(p_val), new_P(i) = p_val; else new_P(i) = def_P; end
            
            x_val = input(sprintf('  -> Point Load %d position from left in meters (Current: %.1f): ', i, def_x));
            if ~isempty(x_val), new_x(i) = x_val; else new_x(i) = def_x; end
        end
        P_old = new_P;
        x_p_old = new_x;
    end

    % 3. UDL Loads
    temp = input(sprintf('\nEnter UDL Magnitudes Array in kN/m (Current: %s): ', mat2str(w_old)));
    if ~isempty(temp), w_old = temp; end

    if ~isempty(w_old)
        temp = input(sprintf('Enter UDL Start Positions in meters (Current: %s): ', mat2str(x_u1_old)));
        if ~isempty(temp), x_u1_old = temp; end
        
        temp = input(sprintf('Enter UDL End Positions in meters (Current: %s): ', mat2str(x_u2_old)));
        if ~isempty(temp), x_u2_old = temp; end
        
        if length(w_old) ~= length(x_u1_old) || length(w_old) ~= length(x_u2_old)
            error('Error: Number of UDL magnitudes, start, and end positions must match!');
        end
    else
        x_u1_old = [];
        x_u2_old = [];
    end

    % 4. Material Properties
    temp = input(sprintf('\nEnter Young''s Modulus E in GPa (Current: %.1f): ', E_old));
    if ~isempty(temp), E_old = temp; end

    temp = input(sprintf('Enter Moment of Inertia I in m^4 (Current: %e): ', I_old));
    if ~isempty(temp), I_old = temp; end

    % Variable Assignment for Calculations
    L = L_old; P = P_old; x_p = x_p_old; w = w_old; x_u1 = x_u1_old; x_u2 = x_u2_old;
    E = E_old * 1e9; I = I_old; EI = E * I;

    %% 2. Calculate Support Reactions
    if isempty(P), M_P = 0; else M_P = sum(P .* x_p); end
    M_U = 0; total_w_load = 0;
    if ~isempty(w)
        for k = 1:length(w)
            L_uk = x_u2(k) - x_u1(k);
            M_U = M_U + w(k) * L_uk * (x_u1(k) + L_uk/2);
            total_w_load = total_w_load + w(k) * L_uk;
        end
    end
    RB = (M_P + M_U) / L;
    if isempty(P), Sum_P = 0; else Sum_P = sum(P); end
    RA = Sum_P + total_w_load - RB;

    fprintf('\n--- RESULTS ---\n');
    fprintf('Reaction RA = %.2f kN, RB = %.2f kN\n', RA, RB);

    %% 3. Discretize and Calculate SF & BM
    n_points = 1000; x = linspace(0, L, n_points);
    SF = zeros(1, n_points); BM = zeros(1, n_points);

    for i = 1:n_points
        xi = x(i); V = RA; M = RA * xi; 
        for j = 1:length(P)
            if xi >= x_p(j), V = V - P(j); M = M - P(j) * (xi - x_p(j)); end
        end
        if ~isempty(w)
            for k = 1:length(w)
                L_uk = x_u2(k) - x_u1(k);
                if xi > x_u1(k)
                    if xi <= x_u2(k)
                        dist = xi - x_u1(k); V = V - w(k) * dist; M = M - w(k) * dist * (dist / 2);
                    else
                        V = V - w(k) * L_uk; M = M - w(k) * L_uk * (xi - (x_u1(k) + L_uk/2));
                    end
                end
            end
        end
        SF(i) = V; BM(i) = M;
    end

    %% 4.5 Slope & Deflection
    curvature = (BM * 1000) / EI; I1 = cumtrapz(x, curvature); I2 = cumtrapz(x, I1);
    C1 = -I2(end) / L; slope = I1 + C1; deflection = (I2 + C1 * x) * 1000; 

    %% 5. Contraflexure
    x_contra = [];
    for i = 1:(n_points - 1)
        if BM(i) * BM(i+1) < 0
            x_contra(end+1) = x(i) - BM(i) * (x(i+1) - x(i)) / (BM(i+1) - BM(i));
        end
    end

    %% 6. Plotting
    figure('Name', 'Beam Analysis: SFD, BMD, Slope & Deflection', 'Position', [50, 20, 1000, 1000]);
    
    % --- Plot 1: Load Diagram (VISUALLY ENHANCED) ---
    subplot(5, 1, 1); hold on; grid on;
    plot([0 L], [0 0], 'k', 'LineWidth', 4); 
    plot(0, 0, '^r', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
    plot(L, 0, '^r', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
    
    % Dynamically scale heights based on loads
    max_load_val = 20;
    if ~isempty(P), max_load_val = max([abs(P), max_load_val]); end
    if ~isempty(w), max_load_val = max([max_load_val, max(w)]); end
    
    % Plot Point Loads (Shorter, cleaner blue arrows)
    for j = 1:length(P)
        arrow_height = max_load_val * 0.45;
        quiver(x_p(j), arrow_height, 0, -arrow_height, 'b', 'LineWidth', 2.5, 'MaxHeadSize', 0.3);
        text(x_p(j), arrow_height * 1.15, sprintf('%d kN', P(j)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % Plot Multiple UDLs
    if ~isempty(w)
        for k = 1:length(w)
            udl_height = max_load_val * 0.35;
            fill([x_u1(k) x_u2(k) x_u2(k) x_u1(k)], [0 0 udl_height udl_height], 'y', 'FaceAlpha', 0.4, 'EdgeColor', 'k', 'LineStyle', ':');
            text((x_u1(k)+x_u2(k))/2, udl_height * 0.5, sprintf('%.1f kN/m', w(k)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
    xlim([-0.5 L+0.5]); 
    ylim([-0.1*max_load_val max_load_val * 1.4]); 
    title('Load Diagram', 'FontWeight', 'bold'); 
    set(gca, 'YTickLabel', []); 
    hold off;

    % --- Plot 2: Shear Force Diagram ---
    subplot(5, 1, 2); hold on; grid on; plot(x, SF, 'b', 'LineWidth', 1.5); fill([x L 0], [SF 0 0], 'b', 'FaceAlpha', 0.1); plot([0 L], [0 0], 'k', 'LineWidth', 1); title('Shear Force Diagram (S.F.D)'); ylabel('SF (kN)'); xlim([0 L]); hold off;

    % --- Plot 3: Bending Moment Diagram ---
    subplot(5, 1, 3); hold on; grid on; plot(x, BM, 'r', 'LineWidth', 1.5); fill([x L 0], [BM 0 0], 'r', 'FaceAlpha', 0.1); plot([0 L], [0 0], 'k', 'LineWidth', 1); title('Bending Moment Diagram (B.M.D)'); ylabel('BM (kNm)'); xlim([0 L]);
    [max_BM, max_idx] = max(BM); plot(x(max_idx), max_BM, 'ko', 'MarkerFaceColor', 'k'); text(x(max_idx), max_BM + (0.15*max_BM), sprintf('Max BM = %.2f kNm (x = %.2fm)', max_BM, x(max_idx)), 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    if ~isempty(x_contra)
        for k = 1:length(x_contra)
            plot(x_contra(k), 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
            text(x_contra(k), -0.2*max(abs(BM)), sprintf('Contraflexure\nx = %.2fm', x_contra(k)), 'HorizontalAlignment', 'center', 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 8);
        end
    end
    hold off;

    % --- Plot 4: Slope Diagram ---
    subplot(5, 1, 4); hold on; grid on; plot(x, slope, 'm', 'LineWidth', 1.5); plot([0 L], [0 0], 'k', 'LineWidth', 1); title('Slope Diagram (\theta)'); ylabel('Slope (rad)'); xlim([0 L]);
    plot(0, slope(1), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 6); text(0.02*L, slope(1), sprintf('\\theta_A = %.6f rad', slope(1)), 'HorizontalAlignment', 'left', 'FontSize', 8, 'BackgroundColor', 'w');
    plot(L, slope(end), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 6); text(0.98*L, slope(end), sprintf('\\theta_B = %.6f rad', slope(end)), 'HorizontalAlignment', 'right', 'FontSize', 8, 'BackgroundColor', 'w');
    [~, max_s_idx] = max(abs(slope)); max_slope_val = slope(max_s_idx); plot(x(max_s_idx), max_slope_val, 'ko', 'MarkerFaceColor', 'k'); text(x(max_s_idx), max_slope_val, sprintf(' Max Slope = %.6f rad (x = %.2fm)', max_slope_val, x(max_s_idx)), 'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'm');
    hold off;

    % --- Plot 5: Deflection Diagram ---
    subplot(5, 1, 5); hold on; grid on; plot(x, deflection, 'g', 'LineWidth', 1.5); fill([x L 0], [deflection 0 0], 'g', 'FaceAlpha', 0.1); plot([0 L], [0 0], 'k', 'LineWidth', 1); title('Deflection Diagram (v)'); xlabel('Beam Length (m)'); ylabel('Deflection (mm)'); xlim([0 L]);
    [max_def_val, def_idx] = min(deflection); plot(x(def_idx), max_def_val, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7); text(x(def_idx), max_def_val - 0.25*abs(max_def_val), sprintf('Max Deflection = %.3f mm\nLocation x = %.2f m', max_def_val, x(def_idx)), 'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', 'w');
    plot(0, 0, 'bs', 'MarkerFaceColor', 'b'); plot(L, 0, 'bs', 'MarkerFaceColor', 'b'); text(0.02*L, 0.05*abs(max_def_val), 'v_A = 0 mm', 'FontSize', 8); text(0.98*L, 0.05*abs(max_def_val), 'v_B = 0 mm', 'HorizontalAlignment', 'right', 'FontSize', 8);
    hold off;
en