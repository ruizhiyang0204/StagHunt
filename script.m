%  ==========================================================
%  Stag Hunt Behavioural Modelling  (k-ToM in VBA-toolbox-master)
%  Author : Ruizhi Yang
%  Date   : 2025-11-10
%  Here is a brief explanation of the data:
%  a   : Action (4 directions, 0–3)
%  h   : Hunting outcome (1 = rabbit, 2 = stag)
%  m   : Game properties (col1 = participants’ action order, col2 = maze index)
    % - Stag moves randomly (non-strategic), treated as environmental noise — not modelled as an agent.
    % - Maze index (m(:,2)) indicates environment layout; has no effect on belief modelling.
%  p   : Positions of 3 agents (2 participants + stag), 2D array
%  r   : Points
%  ot  : Trial onset time (absolute)
%  t   : Action onset time (absolute)
    % Each trial t corresponds to one move step (not necessarily the end of a game)
%  rt  : Reaction time (rt = t - ot)
%  ==========================================================

clc;
addpath(genpath('/Users/nicolas/Desktop/StagHunt/VBA-toolbox-master/modules/theory_of_mind'));



%% 2. 构造二值决策序列 (根据p判断朝鹿/朝兔)
T = size(a,1);
SubjectChoices  = zeros(1,T);
OpponentChoices = zeros(1,T);

for t = 1:T
    % ---- 自己 ----
    self_pre  = squeeze(p(t,1,1,1));
    self_post = squeeze(p(t,1,1,2));
    stag_pre  = squeeze(p(t,3,1,1));
    stag_post = squeeze(p(t,3,1,2));

    move_self_to_stag = abs(self_post - stag_post) < abs(self_pre - stag_pre);
    SubjectChoices(t) = move_self_to_stag;

    % ---- 对手 ----
    opp_pre  = squeeze(p(t,2,1,1));
    opp_post = squeeze(p(t,2,1,2));
    move_opp_to_stag = abs(opp_post - stag_post) < abs(opp_pre - stag_pre);
    OpponentChoices(t) = move_opp_to_stag;
end

%% --- Behavioural choices (improved visualization) ---
figure('Color','w','Name','Binarized behavioural choices','Position',[100 100 800 300]);

hold on;
plot(SubjectChoices, '-o', 'Color',[0 0.45 0.74], 'LineWidth',2, 'MarkerSize',4, ...
     'MarkerFaceColor',[0 0.45 0.74], 'DisplayName','Subject');
plot(OpponentChoices, '-s', 'Color',[0.85 0.33 0.10], 'LineWidth',2, 'MarkerSize',4, ...
     'MarkerFaceColor',[0.85 0.33 0.10], 'DisplayName','Opponent');

% 美化坐标轴
xlabel('Trial','FontSize',14,'FontWeight','bold');
ylabel('Choice (1 = Stag, 0 = Rabbit)','FontSize',14,'FontWeight','bold');
title('Binarized behavioural choices','FontSize',16,'FontWeight','bold');

% 调整外观
xlim([1, length(SubjectChoices)]);
ylim([-0.1, 1.1]);
yticks([0 1]);
yticklabels({'Rabbit','Stag'});

legend('Location','northoutside','Orientation','horizontal','FontSize',12,'Box','off');
grid on; box off;

set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out','TickLength',[0.005 0.005]);
axis tight;



%% 3. 定义Stag Hunt的Payoff Table (JN2010b)
% 满足: High > Med >= Low > Fail
PayoffHigh = 10;   % 双方猎鹿成功
PayoffMed  = 6;    % 单方猎兔
PayoffLow  = 6;    % 双方猎兔
PayoffFail = 0;    % 孤身猎鹿失败

% Hunter 1 (Subject) payoff matrix
U_sub = [PayoffHigh, PayoffFail;  % 自己Stag，对手Stag/Rabbit
         PayoffMed,  PayoffLow];  % 自己Rabbit，对手Stag/Rabbit

% Hunter 2 (Opponent) payoff matrix (对称)
U_opp = [PayoffHigh, PayoffMed;
          PayoffFail, PayoffLow];

% 组合成3D数组: payoffTable(:,:,1)=Subject, payoffTable(:,:,2)=Opponent
payoffTable = cat(3, U_sub, U_opp);



%% 4. 设置k-ToM模型结构
K = 4;                  % recursion depth (2-ToM)
role = 1;               % 当前参与者为Hunter 1
[options, dim] = prepare_kToM(K, payoffTable, role, 0);

f_fname = @f_kToM;      % 演化方程
g_fname = @g_kToM;      % 观测方程

% 输入与输出
y = SubjectChoices;     % 被试选择 (要拟合的输出)
u = [zeros(2,1), [OpponentChoices(1:end-1); SubjectChoices(1:end-1)]];

% 设定VBA选项
orderCodes = unique(m(:,1));
if isequal(orderCodes, [0;1;2])
    orderType = m(:,1) + 1;    % map to 1=A→B→Stag, 2=Stag→A→B, 3=A→Stag→B
elseif isequal(orderCodes, [1;2;3])
    orderType = m(:,1);
else
    warning('Unexpected coding in m(:,1): please inspect manually');
    orderType = m(:,1);
end

% ===========================================
% Dynamic gating based on participant order
% m(:,1): 1=A->B->Stag, 2=Stag->A->B, 3=A->Stag->B
% Hunter1 updates after opponent’s move; Hunter2 updates after A’s move
% ===========================================
g1 = false(T,1);   % gating for Hunter1
g2 = false(T,1);   % gating for Hunter2

for t = 1:T
    lastActor = orderType(t);
    % Hunter1 updates only when opponent (2) just acted
    g1(t) = (lastActor == 2);
    % Hunter2 updates only when opponent (1) just acted
    g2(t) = (lastActor == 1);
end

% For Hunter1’s model (subject = Hunter1)
options.skipf = 1 - g1';   % 1=skip, 0=update

% options.skipf = [1, zeros(1, length(y)-1)]; % 第一trial不学习不更新

options.sources.type = 1;  % 1 = Bernoulli output 二项分布输出





%% 5. 运行VBA非线性状态空间模型反演
[posterior, out] = VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options);

% 输出:
% posterior.muX     - trial-by-trial隐藏状态均值
% posterior.SigmaX → 状态方差
% posterior.muTheta - 参数估计 (遗忘/温度等)
% out.F             - 模型自由能 (model evidence)



%% 6. 解包隐藏状态，提取心理变量
hf = unwrapKTOM(posterior.muX, options.inG);


% 原绘图功能不变
% hf 包含:
% data = get(hf,'UserData');
% pK = data.p_sophistication;   % 信念分布 p_t(k')
% Po = data.prediction;         % 对手下步选择概率
% mP = data.mP;                 % 对手参数均值
% vP = data.vP;                 % 对手参数方差


% 对手层级后验概率 p_t(k)
data = get(hf,'UserData');
pK = data.p_sophistication;   % 信念分布 p_t(k)
% 信念熵 (Belief entropy)
H = -sum(pK .* log(pK + eps), 1);
% 主体推理层级 (期望)
% k_sub = sum(pK .* (1:size(pK,1))', 1);  % 平滑更新
% pK: K x T，行对应对手k'=0..K-1
[~, idx_mode] = max(pK, [], 1);     % 1..K （若K=3, 那么对手可能是0,1,2，假设pK(:,t) = [0.1; 0.7; 0.2]，idx_mode=2，所以算k_op_est要索引号去-1）
k_op_est = idx_mode - 1;            % 0..K-1 （对手估计层级）

% ---------- 2) 被试层级（滞后1拍更新）----------
k_sub = zeros(1,T);

% 用首拍的先验/初值来初始化（用当拍的 MAP + 1，再截断到K）
k_sub(1) = min(k_op_est(1) + 1, K);

% 核心：t+1 才用 t 的 \hat{k_op} + 1
for t = 1:T-1
    k_sub(t+1) = min(k_op_est(t) + 1, K);
end



%% ==========================================================
%  Summary plots
%  ==========================================================

figure('Color','w','Name','ToM behavioural dynamics');

% ---------- (A) Opponent sophistication posterior ----------
subplot(3,1,1); hold on;

colors = lines(size(pK,1));  % 生成不同的颜色（每行一个RGB）
for i = 1:size(pK,1)
    plot(pK(i,:), 'Color', colors(i,:), 'LineWidth', 1.5);
end

xlabel('Trial','Interpreter','latex'); 
ylabel('P(oppo = k\text{-}ToM)','Interpreter','latex');  
title('Opponent''s sophistication','Interpreter','latex');

legend(arrayfun(@(k) sprintf('oppo = %d-ToM',k-1), 1:size(pK,1), 'UniformOutput', false), ...
       'Interpreter','latex');

ylim([0 1]);


% ---------- (B) Subject inferred sophistication ----------
subplot(3,1,2); hold on;
plot(k_op_est, 'b','LineWidth',1.8); 
plot(k_sub,    'g','LineWidth',1.8);
xlabel('Trial','Interpreter','latex'); 
ylabel('Level $(k)$','Interpreter','latex');
title('Opponent vs Subject inferred sophistication','Interpreter','latex');
legend({'$\hat{k}_{\mathrm{op}}$','$k_{\mathrm{sub}}$'}, ...
       'Interpreter','latex','Location','best');
ylim([0, K+0.5]); grid on; box off;

% ---------- (C) Belief uncertainty (entropy) ----------
subplot(3,1,3);
plot(H, 'k', 'LineWidth',1.5);
xlabel('Trail','Interpreter','latex'); 
ylabel('Entropy $\,H(p(k_{\mathrm{op}}))$','Interpreter','latex');
title('Uncertainty of belief inference','Interpreter','latex');
ylim([0, log(size(pK,1))]);

sgtitle('Behavioural inference dynamics','Interpreter','latex');


