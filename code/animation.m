function animation(data, r)
    [num_points, num_times, ~] = size(data);
    
    % 创建图形窗口
    hFig = figure('Name', 'Animation', 'NumberTitle', 'off', 'Position', [100, 100, 600, 500]);
    
    % 创建散点图的初始绘制
    hScatter = scatter(data(:, 1, 1), data(:, 1, 2), 10, 'filled');
    axis([min(data(:)) max(data(:)) min(data(:)) max(data(:))]); % 设置坐标轴范围
    hold on;

    % 如果 r > 0，绘制圆
    if r > 0
        theta = linspace(0, 2*pi, 100);
        x = r * cos(theta);
        y = r * sin(theta);
        plot(x, y, 'r--'); % 绘制圆，使用红色虚线表示
    end
    
    % 用于存储折线图的句柄
    hLines = gobjects(num_points-1, 4);
    
    % 初始化折线图
    for j = 1:num_points-1
        P = generate_dot_(data(j, 1, :), data(j+1, 1, :));
        hLines(j, 1) = plot(P([1, 2], 1), P([1, 2], 2), 'b');
        hLines(j, 2) = plot(P([1, 3], 1), P([1, 3], 2), 'b');
        hLines(j, 3) = plot(P([2, 4], 1), P([2, 4], 2), 'b');
        hLines(j, 4) = plot(P([3, 4], 1), P([3, 4], 2), 'b');
    end
    hold off;

    % 创建进度条
    hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', num_times, 'Value', 1, ...
                        'Position', [100, 20, 400, 20], 'Callback', @slider_callback);
                    
    % 创建+1按钮
    hButtonPlus = uicontrol('Style', 'pushbutton', 'String', '+1', ...
                        'Position', [0, 20, 60, 20], 'Callback', @(~,~) adjust_time(1));


    % 创建-1按钮
    hButtonMinus = uicontrol('Style', 'pushbutton', 'String', '-1', ...
                        'Position', [0, 40, 60, 20], 'Callback', @(~,~) adjust_time(-1));

                    
    % 创建编辑框，用于显示和修改当前时刻
    hEdit = uicontrol('Style', 'edit', 'String', '1', ...
                      'Position', [520, 20, 60, 20], 'Callback', @edit_callback);

    
    % 创建播放按钮
    hButtonPlay = uicontrol('Style', 'pushbutton', 'String', 'Play', ...
                            'Position', [520, 50, 60, 20], 'Callback', @button_callback);
    
    % 创建全局变量控制动画
    playing = false;

    % 进度条回调函数
    function slider_callback(hObject, ~)
        % 获取当前进度条的值
        t = round(get(hObject, 'Value'));
        update_plot(t);
    end

    % 按钮回调函数
    function button_callback(~, ~)
        if ~playing
            % 开始播放
            set(hButtonPlay, 'String', 'Stop');
            playing = true;
            for t = round(get(hSlider, 'Value')):num_times
                if ~playing
                    break;
                end
                update_plot(t);
                set(hSlider, 'Value', t); % 更新进度条位置
                set(hEdit, 'String', num2str(t)); % 更新编辑框显示
                drawnow;
                pause(0.05); % 控制动画播放速度
            end
            % 播放结束
            set(hButtonPlay, 'String', 'Play');
            playing = false;
        else
            % 停止播放
            set(hButtonPlay, 'String', 'Play');
            playing = false;
        end
    end

    % 更新图像显示的函数
    function update_plot(t)
        % 更新散点图
        set(hScatter, 'XData', data(:, t, 1), 'YData', data(:, t, 2));
        % 清除之前的折线图
        delete(hLines);
        % 绘制新的折线图
        hold on;
        for j = 1:num_points-1
            P = generate_dot_(data(j, t, :), data(j+1, t, :));
            hLines(j, 1) = plot(P([1, 2], 1), P([1, 2], 2), 'b');
            hLines(j, 2) = plot(P([1, 3], 1), P([1, 3], 2), 'b');
            hLines(j, 3) = plot(P([2, 4], 1), P([2, 4], 2), 'b');
            hLines(j, 4) = plot(P([3, 4], 1), P([3, 4], 2), 'b');
        end
        hold off;
    end

    % +1和-1按钮的回调函数
    function adjust_time(step)
        current_time = round(get(hSlider, 'Value')) + step;
        current_time = max(1, min(num_times, current_time)); % 确保在合法范围内
        set(hSlider, 'Value', current_time);
        set(hEdit, 'String', num2str(current_time)); % 更新编辑框显示
        update_plot(current_time);
    end

    % 编辑框的回调函数
    function edit_callback(hObject, ~)
        new_time = str2double(get(hObject, 'String'));
        if isnan(new_time) || new_time < 1 || new_time > num_times
            % 如果输入不合法，恢复为当前时间
            new_time = round(get(hSlider, 'Value'));
            set(hObject, 'String', num2str(new_time));
        end
        set(hSlider, 'Value', new_time);
        update_plot(new_time);
    end
end