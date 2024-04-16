function pdf_values = cauchy_pdf(x, x0, gamma)
    % Cauchy probability density function
    pdf_values = 1./(pi * gamma * (1 + ((x - x0)./gamma).^2));
end