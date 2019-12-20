function e_data = rcaEllipses(rc_struct, cond_idx, sub_idx, rc_num, F1, delay_ms, plot_data, plot_offset, cond_labels, cond_colors, error_type, gca_opts)
    
    if nargin < 12 || isempty(error_type)
        error_type = 'SEM';
    else
    end
    
    if nargin < 11 || isempty(gca_opts)
        l_width = 1.5;
        f_size = 12;
        text_opts = {'fontweight','normal','fontname','Helvetica','fontsize', f_size};
        gca_opts = [{'tickdir','out','box','off','fontsize',12,'fontname','Helvetica','linewidth',l_width,'Color','none'}, text_opts{:}];;
    else
        f_size = 12;
        text_opts = {'fontweight','normal','fontname','Helvetica','fontsize', f_size};
    end
    
    if nargin < 10 || isempty(cond_colors)
        c_brewer = load('colorBrewer_new.mat');
        cond_colors = c_brewer.rgb20(1:2:end,:);
        cond_colors([2,3],:) = cond_colors([3,2],:);
    else
    end
   
    if nargin < 9
        cond_labels = [];
    end
    
    if nargin < 8 || isempty(plot_offset)
        plot_offset = 0; %pi/2;
    else
    end
    
    if nargin < 7 || isempty(plot_data)
        plot_data = true;
    else
    end
    
    if nargin < 6 || isempty(delay_ms)
        delay_ms = 66;
    else
    end
    
    if nargin < 5 || isempty(F1)
        F1 = 1; 
    else
    end
    
    if nargin < 4 || isempty(rc_num)
        rc_num = 1;
    else
    end
    
    if nargin < 3 || isempty(sub_idx)
        sub_idx = true(1, size(rc_struct.rca_data,2));
    else
    end
    
    if nargin < 2 || isempty(cond_idx)
        cond_idx = true(1, size(rc_struct.rca_data,1));
    else
    end
    
    rca_real = rc_struct.subjects.real_signal(:,:,:,sub_idx, cond_idx);
    rca_imag = rc_struct.subjects.imag_signal(:,:,:,sub_idx, cond_idx);
    rca_real = squeeze(permute(rca_real,[4,5,1,2,3]));
    rca_imag = squeeze(permute(rca_imag,[4,5,1,2,3]));
    
    if any(~size(rca_real) == size(rca_imag))
        msg = 'real and imaginary matrices are not the same size!';
        error(msg);
    else
    end
    
    n_conds = size(rca_real, 2);
    n_bins = size(rca_real, 3);
    n_comp = size(rca_real, 4);
    
    if isempty(cond_labels)
        cond_labels = arrayfun(@(x) sprintf('condition %01d', x), 1:n_conds, 'uni', false);
    else
    end
    
    harmonic = cell2mat(cellfun(@(x) str2num(x(1)), ... 
        rc_struct.settings.freqLabels(rc_struct.settings.freqIndices), 'uni', false));
    
    if length(unique(harmonic)) > 1
        msg = 'harmonics are different across input';
        error(msg);
    else
        harmonic = harmonic(1);
    end
    
    cycle_len = 1000/(F1*harmonic);
    
    if numel(rc_num) > 1
        msg = 'rcaEllipses only works for one RC at a time';
        error(msg);
    else
    end
    
    % loop over rc
    for r = 1:n_comp   
        for b = 1:n_bins
            for c = 1:n_conds
                xy_vals = cat(2,rca_real(sub_idx, c, b, r), rca_imag(sub_idx, c, b, r));
                nan_vals = sum(isnan(xy_vals),2) > 0;
                [amp_err, phase_err, ~, error_ellipse] = fitErrorEllipse(xy_vals(~nan_vals,:),  error_type);
                e_data(b,c,r).ellipse = error_ellipse;
                % stack mean and ellipse data together
                temp_data = [nanmean(xy_vals); e_data(b,c,r).ellipse];
                % compute phase and amplitude for each data point
                temp_phase = angle(complex(temp_data(:,1), temp_data(:,2)));
                temp_amp = abs(complex(temp_data(:,1),temp_data(:,2)));
                % add delay 
                delay_rad = delay_ms/cycle_len*2*pi;
                temp_phase = temp_phase + delay_rad;
                % unwrap so that there are no negative values
                temp_phase( temp_phase < 0 ) = 2*pi+temp_phase( temp_phase < 0 );
                % convert phase and amplitude back into complex numbers
                temp_complex = temp_amp.*exp(1i*temp_phase);
                temp_data = [real(temp_complex), imag(temp_complex)];
                
                % add values back into the ellipse struct
                e_data(b, c, r).ellipse = temp_data(2:end,:);
                e_data(b, c, r).mean.real = temp_data(1,1);
                e_data(b, c, r).mean.imag = temp_data(1,2);
                
                % vector average amplitude
                e_data(b, c, r).mean.amp = temp_amp(1,:);
                % convert to response latency in ms
                e_data(b, c, r).mean.latency = (temp_phase(1,:)/(2*pi)*cycle_len)';
                
                % now handle the phase error
                phase_err = cat(1, phase_err);
                % convert from degrees to radians
                phase_err = phase_err .* pi / 180; 
                % convert errors to response latency in ms
                e_data(b, c, r).err.latency = (phase_err ./(2*pi)*cycle_len);
                e_data(b, c, r).err.amp = amp_err;
            end
        end
    end
    clear rca_real; clear rca_imag;
    
    if plot_data  
        % PLOT THE DATA

        for r = 1:length(rc_num) 
            % create figure with appropriate subplots
            ellipse_fig = figure;
            set(ellipse_fig,'units','centimeters');
            fig_pos = get(gcf,'pos');
            fig_pos(3) = 50;
            fig_pos(4) = 50;
            set(gcf,'pos', fig_pos);
            amp_plot = subplot(2, 2, 1);
            phase_plot = subplot(2, 2, 2);
            for b = 1:n_bins
                bin_plot(b) = subplot(3,n_bins,n_bins*2+b);
                axis square;
                plot_pos = get(bin_plot(b), 'pos');
                if b == 1
                    y_pos = plot_pos(2)+plot_pos(4)*1.2; 
                    x_pos = plot_pos(1);
                    x_span(1) = x_pos(1);
                else
                    x_pos = x_pos + plot_pos(3);
                    x_span(2) = x_pos(1)+plot_pos(3);
                end
                plot_pos(1) = x_pos;
                plot_pos(2) = y_pos;
                set(bin_plot(b), 'pos', plot_pos);
            end

            amp_pos = get(amp_plot, 'pos');
            amp_pos(1) = x_span(1);
            amp_pos(3) = diff(x_span)*.46;
            set(amp_plot, 'pos', amp_pos);
            phase_pos = get(phase_plot, 'pos');
            phase_pos(3) = diff(x_span)*.46;
            phase_pos(1) = x_span(2)-phase_pos(3);
            set(phase_plot, 'pos', phase_pos);

            % compute bin values and labels
            if ~any(cell2mat(cellfun(@(x) isnumeric(x), rc_struct.settings.binLabels, 'uni', false)))
                log_step = 1;
                bin_vals = (1:length(rc_struct.settings.binLabels))';
                bin_text = rc_struct.settings.binLabels;
                x_min = 0.5;
                if contains('ave', bin_text)
                    bin_vals = bin_vals(1:end-1);
                    extra_bins = bin_vals(end) + 1;
                    x_max = length(bin_vals)+1.5;
                else
                    x_max = length(bin_vals)+.5;
                    extra_bins = [ ];
                end
            else
                bin_vals = cell2mat(cellfun(@(x) str2num(x),rc_struct.settings.binLabels,'uni',false));
                bin_text = arrayfun(@(x) num2str(x,'%.1f'),bin_vals,'uni',false);
                log_step = diff(reallog(bin_vals(1:2))); % step size
                x_min = reallog(bin_vals(1))-log_step*.5;
                if plot_ave
                    bin_text = [bin_text; 'ave'];
                    x_max = reallog(bin_vals(end))+log_step*2.5; % add 2.5 steps
                    extra_bins = arrayfun(@(x) reallog(bin_vals(end))+x, [log_step,log_step*2]);
                else
                    x_max = reallog(bin_vals(end))+log_step*.5; % add 2.5 steps
                    extra_bins = [];
                end
                bin_vals = reallog(bin_vals);
            end

            for b = 1:n_bins
                if contains(bin_text{b}, 'bin')
                    bin_text{b} = bin_text{b}(4:end);
                else
                end
                if b == 1
                    e_max = ceil(max(max(abs(cat(1,e_data(:,:,rc_num(r)).ellipse)))));
                    if e_max == 0
                        e_max = 1;
                    else
                    end 
                else
                end
                if e_max > 5
                    e_units = 4;
                elseif e_max > 1
                    e_units = 1;
                else
                    e_units = .5;
                end
                subplot(bin_plot(b));
                hold on;
                grid_h(1) = plot(zeros(1,2),[-e_max,e_max*.75],'k-','linewidth',l_width);
                grid_h(2) = plot([-e_max,e_max],zeros(1,2),'k-','linewidth',l_width);
                for c = 1:n_conds
                    temp_data = [e_data(b,c,rc_num(r)).mean.real, e_data(b,c,rc_num(r)).mean.imag; ...
                        e_data(b,c,rc_num(r)).ellipse];                

                    if plot_offset > 0
                        % add shift for plotting purposes
                        temp_phase = angle(complex(temp_data(:,1), temp_data(:,2)));
                        temp_amp = abs(complex(temp_data(:,1),temp_data(:,2)));                
                        temp_phase = temp_phase + plot_offset;
                        % unwrap so that there are no negative values
                        temp_phase( temp_phase < 0 ) = 2*pi+temp_phase( temp_phase < 0 );
                        % convert phase and amplitude back into complex numbers
                        temp_complex = temp_amp.*exp(1i*temp_phase);
                        temp_data = [real(temp_complex), imag(temp_complex)];
                    else
                    end

                    p_h(c) = plot([0 temp_data(1,1)],[0 temp_data(1,2)],'-','linewidth',l_width,'Color',cond_colors(c,:));
                    plot(temp_data(2:end,1), temp_data(2:end,2),'-','linewidth',l_width,'Color',cond_colors(c,:));
                end
                text(0,e_max,bin_text(b),'fontsize',12,'fontname','Helvetica','horizontalalignment','center');
                xlim([-e_max,e_max]);
                ylim([-e_max,e_max]);
                axis square;
                if b == 1
                    cycle_steps = linspace(0,2*pi-(pi/2),200)+plot_offset;
                    [cycle_x,cycle_y] = pol2cart(cycle_steps, e_max*.75);
                    cycle_h = plot(cycle_x,cycle_y,'-','color',[200,200,200]./255,'linewidth',l_width);        
                    if abs(cycle_x(end)) < 1.0e-12
                        % horizontal arrow
                        arrow_sign = sign(cycle_x(end-1));
                        for z = 1:2
                            arrow_x(1) = cycle_x(end);
                            arrow_x(2) = cycle_x(end)+e_max*.2*arrow_sign;
                            arrow_y(1) = cycle_y(end);
                            if z == 1
                                arrow_y(2) = cycle_y(end)+e_max*.1;
                            else
                                arrow_y(2) = cycle_y(end)+e_max*.1*-1;
                            end
                            q_h(z) = plot(arrow_x, arrow_y, '-', 'color', [200,200,200]./255, 'linewidth',l_width);
                        end
                    elseif abs(cycle_x(end)) < 1.0e-12
                        % vertical arrow
                        arrow_sign = sign(cycle_y(end-1));
                        for z = 1:2
                            arrow_y(1) = cycle_y(end);
                            arrow_y(2) = cycle_y(end)+e_max*.2*arrow_sign;
                            arrow_x(1) = cycle_x(end);
                            if z == 1
                                arrow_x(2) = cycle_x(end)+e_max*.1;
                            else
                                arrow_x(2) = cycle_x(end)+e_max*.1*-1;
                            end
                            q_h(z) = plot(arrow_x, arrow_y, '-', [200,200,200]./255, 'linewidth',l_width);
                        end
                    else
                        msg = "your cycles does not end at x == 0 or y == 0, something is wrong";
                        error(msg);
                    end
                    uistack(q_h,'bottom');
                    arrayfun(@(x) uistack(x,'bottom'),grid_h);
                    cycle_h = plot(cycle_x,cycle_y,'-','color',[200,200,200]./255,'linewidth',l_width);
                    uistack(cycle_h,'top');
                    set(gca,gca_opts{:},'ticklength',[0.1,0.1],'clipping','off','visible','on','xtick',-e_max:e_units:e_max,'ytick',-e_max:e_units:e_max);
                    text(e_max*1.2,-e_max*1.2,num2str(cycle_len,'cycle length: %.01f ms'), text_opts{:});
                    xlabel('real','fontsize',12,'fontname','Helvetica');
                    ylabel('imag','fontsize',12,'fontname','Helvetica');
                else
                    set(gca,gca_opts{:},'clipping','off','visible','off');
                    if b == 6
                        text(0, e_max*1.5,'displacement (arcmins)','fontsize',12,'fontname','Helvetica','horizontalalignment','center');
                    else
                    end
                end
               hold off
            end
            drawnow;
            for z = 1:2
                for c = 1:n_conds
                    if z == 1
                        % plot amplitudes
                        big_plot(z) = subplot(amp_plot);
                        big_min = 0; big_max = 10; big_unit = 2; %ceil(max(max(squeeze(rc_struct.mean.amp_signal(:,:,rc_num(r))))));
                        plot_vals = ...
                            cell2mat(arrayfun(@(x) e_data(x,c,rc_num(1)).mean.amp, 1:n_bins, 'uni', false));
                        plot_errs = ...
                            cell2mat(arrayfun(@(x) e_data(x,c,rc_num(1)).err.amp', 1:n_bins, 'uni', false));
                        y_text = 'amplitude (\muV)';
                    else
                         % plot phase
                        big_plot(z) = subplot(phase_plot);
                        big_min = 1/F1 * 0; big_max = 1/F1 * 1000; big_unit = 100;
                        plot_vals = ...
                            cell2mat(arrayfun(@(x) e_data(x,c,rc_num(1)).mean.latency, 1:n_bins, 'uni', false));
                        plot_errs = ...
                            cell2mat(arrayfun(@(x) e_data(x,c,rc_num(1)).err.latency', 1:n_bins, 'uni', false));
                        y_text = 'delay (ms)';
                    end
                    hold on

                    if ~isempty(extra_bins)
                        %if c > 2
                        %    cur_x = [bin_vals', extra_bins(2)];
                        %else
                        cur_x = [bin_vals', extra_bins];
                        %end
                        if c == 1
                            plot(ones(1,2)*(bin_vals(end)+log_step*.5) ...
                                ,[big_min,big_max],'k-','linewidth',l_width);
                        else
                        end 
                        p_h(c) = plot(cur_x(1:end-1),plot_vals(1:end-1),'-o','linewidth',l_width,'Color',cond_colors(c,:),'markersize',10,'markerfacecolor','w','markeredgecolor',cond_colors(c,:));
                        plot(cur_x(end),plot_vals(end),'o','linewidth',l_width,'Color',cond_colors(c,:),'markersize',10,'markerfacecolor','w','markeredgecolor',cond_colors(c,:));
                    else
                        cur_x = bin_vals';
                        plot(cur_x,plot_vals,'-o','linewidth',l_width,'Color',cond_colors(c,:),'markersize',10,'markerfacecolor','w','markeredgecolor',cond_colors(c,:));
                    end
                    hold on
                    h_e(c,:) = ErrorBars(cur_x,plot_vals,plot_errs,'type','bar','color',cond_colors(c,:),'cap',false);
                    hold on
                end
                if z == 1
                    legend(p_h,cond_labels,'location','northwest','box','off','fontsize',12,'fontname','Helvetica');
                end

                cellfun(@(x) uistack(x,'bottom'),h_e,'uni',false);
                xlim([x_min,x_max]);
                ylim([big_min,big_max]);

                if isnan(bin_vals(1))
                    xtick_vals = reallog([0.2,0.5,1,2,4,8,16]);
                    xtick_labels = arrayfun(@(x) num2str(x),([0.2,0.5,1,2,4,8,16]),'uni',false);

                    if contains('ave', bin_text)
                        xtick_vals = [xtick_vals,xtick_vals(end) + log_step*1.5];
                        xtick_labels = [xtick_labels,'ave'];
                    else
                    end
                else
                    xtick_vals = 1:length(bin_vals)+1;
                    xtick_labels = bin_text;
                end
                set(gca,gca_opts{:},'ticklength',[0.02, 0.02],'clipping','off','visible','on','xtick',xtick_vals,'xticklabels',xtick_labels, 'ytick', [big_min:big_unit:big_max]);
                ylabel(y_text,'fontsize',12,'fontname','Helvetica');
                hold off
            end
            drawnow;
            tightfig;
            out_path = split(pwd, filesep);
            out_path = cell2mat(join([out_path(1:3); {'Desktop'}], filesep));
            
            if rc_num(r) == n_comp
                export_fig(ellipse_fig, sprintf('%s%sellipseplot_comp.pdf', out_path, filesep) ,'-transparent','-painters');
            else
                export_fig(ellipse_fig, sprintf('%s%sellipseplot_rc%d.pdf', out_path, filesep, rc_num(r)),'-transparent','-painters');
            end
        end
    else
    end
end