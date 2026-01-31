% --- 2パラメータモデルによるカーブフィッティング プログラム ---

function main()
    % --- 初期設定 ---
    clear; 
    clc;   
    close all; 

    % --- FEAの解析データ（"正解"データ）を入力 ---
    % 短冊モデルの解析結果
    x_push_data = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    delta_data  = [1.2086, 1.4589, 1.6824, 1.8834, 2.0661, 2.2341, 2.3901, 2.5358, 2.6729, 2.8025];

    % --- 2パラメータモデルの定義 ---
    % x_push = c1 * delta^2 + c2 * delta^4 というモデル式を定義
    model = fittype('c1 * delta^2 + c2 * delta^4', ...
                    'independent', 'delta', 'coefficients', {'c1', 'c2'});
    
    % --- カーブフィッティングの実行 ---
    % 全てのデータに最もよくフィットする c1 と c2 を見つけ出す
    [fit_result, gof] = fit(delta_data', x_push_data', model);
    
    % フィット結果から係数を取得
    C1_fit = fit_result.c1;
    C2_fit = fit_result.c2;

    % --- 結果の表示 ---
    fprintf('=======================================================\n');
    fprintf('<< 2パラメータモデルによるカーブフィット結果 >>\n\n');
    fprintf('フィッティングで得られた係数:\n');
    fprintf('C1 = %.4f\n', C1_fit);
    fprintf('C2 = %.4f\n', C2_fit);
    fprintf('R-squared (決定係数、1に近いほど良い): %.6f\n', gof.rsquare);
    fprintf('=======================================================\n');

    % --- グラフの描画 ---
    figure;
    hold on;
    
    % 1. 元の解析データ（正解）をプロット
    plot(x_push_data, delta_data, 'o-', 'LineWidth', 2, 'DisplayName', 'Analysis (FEA)');
    
    % 2. フィットしたモデルの線をプロット
    %    滑らかな線を描くために細かいdeltaのベクトルを生成
    delta_fine = linspace(min(delta_data), max(delta_data), 100);
    x_push_fit = C1_fit * delta_fine.^2 + C2_fit * delta_fine.^4;
    plot(x_push_fit, delta_fine, '--', 'LineWidth', 2, 'DisplayName', '2-Parameter Fit Model');
    
    hold off;
    grid on;
    legend;
    xlabel('押し込み量 (x_{push}) [mm]');
    ylabel('最大たわみ量 (\delta_0) [mm]');
    title('Curve Fitting with 2-Parameter Model');
end