clc, clear, close all
% 计算三种格式的结果
[W_js, kList] = computeWList(@weno5jse16); % WENO-JS
[W_ze, ~] = computeWList(@weno5ze16); % WENO-ZE
[W_uw5, ~] = computeWList(@uw5); % UW5

% 绘制实部
figure
hold on
plot(kList, real(W_js), 'LineWidth', 1.5, 'DisplayName', 'WENO-JS')
plot(kList, real(W_ze), '--', 'LineWidth', 1.5, 'DisplayName', 'WENO-Z')
plot(kList, real(W_uw5), '-.', 'LineWidth', 1.5, 'DisplayName', 'UW5')

xlabel('$k$', 'Interpreter', 'latex'), ylabel('$\mathrm{Re}(\tilde{k})$', 'Interpreter', 'latex')
title('dispersion'), grid on, legend('Location', 'best')

% 绘制虚部
figure
hold on
plot(kList, imag(W_js), 'LineWidth', 1.5, 'DisplayName', 'WENO-JS')
plot(kList, imag(W_ze), '--', 'LineWidth', 1.5, 'DisplayName', 'WENO-Z')
plot(kList, imag(W_uw5), '-.', 'LineWidth', 1.5, 'DisplayName', 'UW5')
xlabel('$k$', 'Interpreter', 'latex'), ylabel('$\mathrm{Im}(\tilde{k})$', 'Interpreter', 'latex')
title('dissipation'), grid on, legend('Location', 'best')

%%
function [PhiN, phiNList] = computeWList(getDif)
    L = 2;
    tau = 1e-10;
    N = 202;
    h = L / N;

    vj_0 = zeros(1, N);
    vj_tau = zeros(1, N);
    PhiN = zeros(1, N / 2 + 1);

    for n = 0:N / 2
        nIndex = n + 1;
        phiN = 2 * pi * n / N;

        for j = 0:N - 1
            jIndex = j + 1;

            vj_0(jIndex) = exp(1i * j * phiN);

            uavemmm = exp(1i * (j - 3) * phiN);
            uavemm = exp(1i * (j - 2) * phiN);
            uavem = exp(1i * (j - 1) * phiN);
            uave = exp(1i * (j + 0) * phiN);
            uavep = exp(1i * (j + 1) * phiN);
            uavepp = exp(1i * (j + 2) * phiN);

            f_plus = getDif(real([uavemm, uavem, uave, uavep, uavepp])) + ...
                1i * getDif(imag([uavemm, uavem, uave, uavep, uavepp]));
            f_minus = getDif(real([uavemmm, uavemm, uavem, uave, uavep])) + ...
                1i * getDif(imag([uavemmm, uavemm, uavem, uave, uavep]));

            dd =- (f_plus - f_minus) / h;

            vj_tau(jIndex) = vj_0(jIndex) + dd * tau;
        end

        v_0_hat = 1 / N * sum(vj_0 .* exp(-1i * (0:N - 1) * phiN));
        v_tau_hat = 1 / N * sum(vj_tau .* exp(-1i * (0:N - 1) * phiN));

        PhiN(nIndex) = 1i * h / tau * log(v_tau_hat / v_0_hat);

    end

    phiNList = (0:N / 2) * 2 * pi / N;
end

function h = uw5(var)
    uavemm = var(1);
    uavem = var(2);
    uave = var(3);
    uavep = var(4);
    uavepp = var(5);

    h = 1/60 * sum([2, -13, 47, 27, -3] .* [uavemm, uavem, uave, uavep, uavepp]);
end

function h = weno5jse16(var)
    v_mm = var(1);
    v_m = var(2);
    v_0 = var(3);
    v_p = var(4);
    v_pp = var(5);

    epsilon = 1e-16;

    % 计算每个模板上的重构值
    h0 = (2 * v_mm - 7 * v_m + 11 * v_0) / 6;
    h1 = (-v_m + 5 * v_0 + 2 * v_p) / 6;
    h2 = (2 * v_0 + 5 * v_p - v_pp) / 6;

    % 计算光滑指示器
    beta0 = 13/12 * (v_mm - 2 * v_m + v_0) ^ 2 +1/4 * (v_mm - 4 * v_m + 3 * v_0) ^ 2;
    beta1 = 13/12 * (v_m - 2 * v_0 + v_p) ^ 2 +1/4 * (v_m - v_p) ^ 2;
    beta2 = 13/12 * (v_0 - 2 * v_p + v_pp) ^ 2 +1/4 * (3 * v_0 - 4 * v_p + v_pp) ^ 2;

    % 线性权重
    d0 = 1/10;
    d1 = 6/10;
    d2 = 3/10;

    % 计算权重
    w0 = d0 ./ (epsilon + beta0) ^ 2;
    w1 = d1 ./ (epsilon + beta1) ^ 2;
    w2 = d2 ./ (epsilon + beta2) ^ 2;

    S = w0 + w1 + w2;
    w0 = w0 ./ S;
    w1 = w1 ./ S;
    w2 = w2 ./ S;

    % 重构
    h = w0 .* h0 + w1 .* h1 + w2 .* h2;

end

function h = weno5ze16(var)
    v_mm = var(1);
    v_m = var(2);
    v_0 = var(3);
    v_p = var(4);
    v_pp = var(5);

    epsilon = 1e-16;

    % 计算每个模板上的重构值
    h0 = (2 * v_mm - 7 * v_m + 11 * v_0) / 6;
    h1 = (-v_m + 5 * v_0 + 2 * v_p) / 6;
    h2 = (2 * v_0 + 5 * v_p - v_pp) / 6;

    % 计算光滑指示器
    beta0 = 13/12 * (v_mm - 2 * v_m + v_0) ^ 2 +1/4 * (v_mm - 4 * v_m + 3 * v_0) ^ 2;
    beta1 = 13/12 * (v_m - 2 * v_0 + v_p) ^ 2 +1/4 * (v_m - v_p) ^ 2;
    beta2 = 13/12 * (v_0 - 2 * v_p + v_pp) ^ 2 +1/4 * (3 * v_0 - 4 * v_p + v_pp) ^ 2;

    % 线性权重
    d0 = 1/10;
    d1 = 6/10;
    d2 = 3/10;

    % 计算权重
    tau = abs(beta0 - beta2);
    w0 = d0 .* (1 + (tau / (epsilon + beta0)) ^ 2);
    w1 = d1 .* (1 + (tau / (epsilon + beta1)) ^ 2);
    w2 = d2 .* (1 + (tau / (epsilon + beta2)) ^ 2);

    S = w0 + w1 + w2;
    w0 = w0 ./ S;
    w1 = w1 ./ S;
    w2 = w2 ./ S;

    % 重构
    h = w0 .* h0 + w1 .* h1 + w2 .* h2;

end
