function use_parallel = start_parallel_pool()
%START_PARALLEL_POOL  並列計算の準備をする
%
%   use_parallel = start_parallel_pool()
%
%   Parallel Computing Toolbox が使えれば並列プールを起動して true を返す。
%   使えなければ false を返す（逐次計算になる。結果は同じで、時間が長くかかる）。

use_parallel = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
if use_parallel && isempty(gcp('nocreate'))
    parpool;
end
end
