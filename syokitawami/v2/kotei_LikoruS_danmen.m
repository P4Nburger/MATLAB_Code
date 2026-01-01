% --- 半経験的モデル プログラム (最終版) ---


function main()
    % --- 初期設定 ---
    clear; 
    clc;   
    close all; 

    % --- モデルの選択 ---
    % model_type = 'strip';    % 短冊モデルを使用
     model_type = 'diamond';  % ひし形モデルを使用

    % --- パラメータ設定 ---
    
    % FEAの結果から逆算した有効形状定数C
    %C_strip   = 0.2548; % 短冊モデル用
    C_diamond = 0.2794; % ひし形モデル用
    
    % 選択されたモデルに応じてC_effを決定
    if strcmp(model_type, 'strip')
        C_eff = C_strip;
        model_name = 'Strip Model';
    else
        C_eff = C_diamond;
        model_name = 'Diamond Shape Model';
    end
    
    % --- グラフ描画用の押し込み量ベクトル ---
    x_push_vec = 0.1:0.1:1;

    % --- 計算 ---
    % 補正されたシンプルな式で、たわみ量を直接計算
    delta_0_vec = sqrt(2 * x_push_vec / C_eff);
    
    % --- 結果の表示 ---
    fprintf('=======================================================\n');
    fprintf('<< 半経験的モデルによる予測結果 (%s) >>\n', model_name);
    fprintf('使用した有効形状定数 C_eff = %.4f\n\n', C_eff);
    
    % テーブル形式で表示
    fprintf('押し込み量[mm] | 予測たわみ量[mm]\n');
    fprintf('----------------|------------------\n');
    for i = 1:length(x_push_vec)
        fprintf('     %.2f      |      %.4f\n', x_push_vec(i), delta_0_vec(i));
    end
    fprintf('=======================================================\n');
    
    % --- グラフの描画 ---
    figure;
    plot(x_push_vec, delta_0_vec, '-o', 'LineWidth', 2);
    grid on;
    xlabel('押し込み量 (x_{push}) [mm]');
    ylabel('予測される最大たわみ量 (\delta_0) [mm]');
    title(['Semi-Empirical Model: Push-in vs. Deflection (' model_name ')']);
    xlim([0, 1.1]);
    ylim([0, max(delta_0_vec) * 1.1]);
end